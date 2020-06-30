function [time, freq, prob] = maps_f0(signal, samplerate, blocksize, hopsize, basefrequencies)
    % A base frequency estimation in the magnitude and phase spectrum.
    %
    % MaPS-f0 uses a signal's magnitude STFT and it's IF deviation to
    % estimate the maximum-likelihood base frequency of a tone complex
    % for a given series of base frequencies.
    %
    % blocks: A SignalBlocks instance.
    % basefrequencies: an ordered vector of base frequencies in Hz.
    %
    % Returns a vector of times, a vector of likely base frequencies,
    % and a vector of their likelihood.

    if ~exist('blocksize', 'var') || isempty(blocksize)
        blocksize = 2048;
    end
    if ~exist('hopsize', 'var') || isempty(hopsize)
        hopsize = 1024;
    end
    if ~exist('basefrequencies', 'var') || isempty(basefrequencies)
        basefrequencies = logspace(log10(50), log10(450), 200);
    end
    probabilities = pitch_probability(signal, samplerate, blocksize, hopsize, basefrequencies);
    [time, freq, prob] = max_track_viterbi(probabilities, hopsize/samplerate, basefrequencies);
end

function [prob] = pitch_probability(signal, samplerate, blocksize, hopsize, basefrequencies)
    % A probability-of-pitch estimation in the magnitude and phase spectrum.
    %
    % MaPS-f0 uses a signal's magnitude STFT and it's IF deviation to
    % estimate the likelihood that a given range of base frequencies are
    % base frequencies of a tone complex.
    %
    % blocks: A SignalBlocks instance.
    % basefrequencies: an ordered vector of base frequencies in Hz.
    %
    % Returns a len(blocks) x len(basefrequencies) matrix of likelihood values.

    mc = magnitude_correlation(signal, samplerate, blocksize, hopsize, basefrequencies);
    ifdd = ifd_difference(signal, samplerate, blocksize, hopsize, basefrequencies);
    prob = value2posterior(mc, ifdd);
end

function [correlation] = magnitude_correlation(signal, samplerate, blocksize, hopsize, basefrequencies)
    % Correlate synthetic tone complex spectra with a true spectrum.
    %
    % Generate spectra at a number of given base frequencies, then
    % correlate each of these spectra with the magnitude signal
    % STFT-spectra.
    %
    % Before correlation, each frequency bin is log-weighted to make the
    % correlation perceptually accurate.
    %
    % blocks: A SignalBlocks instance.
    % basefrequencies: An ordered list of base frequencies in Hz.
    %
    % Returns a len(blocks) x len(basefrequencies) matrix of correlation values.

    specsize = floor(blocksize / 2) + 1;

    % weigh differences according to perception:
    f = linspace(0, samplerate/2, specsize);
    log_f_weight =  1 ./ (samplerate/2).^(f / (samplerate/2));

    correlation = zeros(numblocks(signal, blocksize, hopsize), length(basefrequencies));
    synthetic_magnitudes = synthetic_magnitude(samplerate, specsize, basefrequencies);

    for blockidx = 1:numblocks(signal, blocksize, hopsize)
        spectrum = stft(signal, samplerate, blocksize, hopsize, blockidx);
        correlation(blockidx, :) = sum(abs(spectrum)' .* synthetic_magnitudes .* log_f_weight, 2)';
    end
end

function [num] = numblocks(signal, blocksize, hopsize)
    % The number of blocks in a signal

    num = ceil( (size(signal, 1)-blocksize) / hopsize);
end

function [spectrum] = stft(signal, samplerate, blocksize, hopsize, blockidx, nfft, windowfunc)
    % Short-time Fourier Transform of a signal.
    %
    % The signal is cut into short overlapping blocks, and each block is
    % transformed into the frequency domain using the FFT. Before
    % transformation, each block is windowed by windowfunc.
    %
    % blocks: A SignalBlocks instance.
    % nfft: None for blocksize, or a number.
    % windowfunc: A function that returns a window.
    %
    % Returns a complex spectrum.

    if ~exist('nfft', 'var') || isempty(nfft)
        nfft = blocksize;
    end
    if ~exist('windowfunc', 'var') || isempty(windowfunc)
        windowfunc = @hann;
    end

    window = windowfunc(blocksize);
    specsize = floor(nfft/2) + 1;

    block = signal(((blockidx-1)*hopsize+1):((blockidx-1)*hopsize + blocksize));
    spectrum = fft(block.*window, nfft);
    spectrum = spectrum(1:specsize);
end

function [magnitudes] = synthetic_magnitude(samplerate, specsize, basefrequencies)
    % Synthetic magnitude spectra of a range of tone complexes.
    %
    % samplerate: The sampling rate of the tone complexes.
    % specsize: The length of each spectrum.
    % basefrequencies: An ordered vector of tone complex base
    %     frequencies in Hz.
    %
    % Returns a len(basefrequencies) x specsize matrix of tone complex spectra.

    freqs = linspace(0, samplerate/2, specsize);
    magnitudes = zeros(length(basefrequencies), specsize);
    for freqidx = 1:length(basefrequencies)
        magnitudes(freqidx, :) = hannwin_comb(samplerate, basefrequencies(freqidx), specsize);
    end
end

function [comb_spectrum] = hannwin_comb(samplerate, basefrequency, specsize)
    % Approximate a speech-like correlation spectrum of a tone complex.
    %
    % This is an approximation of time_domain_comb that runs much
    % faster.
    %
    % Instead of calculating the FFT of a series of hann-windowed
    % sinuses, this models the spectrum of a tone-complex as a series of
    % hann-window-spectrums.
    %
    % For a perfect reconstruction, this would need to calculate the sum
    % of many hann-window-spectra. Since hann window spectra are very
    % narrow, this assumes that each window spectrum extends from
    % n*basefreq-basefreq/2 to n*basefreq+basefreq/2 and that
    % neighboring spectra do not influence each other.
    %
    % This assumtion holds as long as basefreq >> 1/specsize.
    %
    % Amplitudes are normalized by specsize.
    %
    % To make the spectrum more speech-like, frequencies above 1000 Hz
    % are attenuated by 24 dB/oct.
    %
    % To make the correlation of this spectrum and some other spectrum
    % have a normalized gain, the spectrum is shifted to be zero-mean.
    %
    % samplerate: The sampling rate in Hz of the signal.
    % basefreq: The base frequency in Hz of the tone complex.
    % specsize: The length of the resulting spectrum in bins
    %           (typically 2**N+1 for type(N) == int).
    %
    % Returns a real magnitude spectrum.

    freqs = linspace(0, samplerate/2, specsize);
    % create a local frequency vector around each harmonic, going from
    % -basefreq/2 to basefreq/2 within the area around the nth
    % harmonic n*basefreq-basefreq/2 to n*basefreq+basefreq/2:
    closest_harmonic = floor((freqs + basefrequency/2) / basefrequency);
    % ignore first half-wave:
    closest_harmonic(closest_harmonic==0) = 1;
    local_frequency = closest_harmonic*basefrequency - freqs;
    % convert from absolute frequency to angular frequency:
    local_angular_freq = local_frequency / (samplerate/2) * 2*pi;
    % evaluate hannwin_spectrum at the local frequency vector:
    comb_spectrum = abs(hannwin_spectrum(local_angular_freq, specsize));
    % normalize to zero mean:
    comb_spectrum = comb_spectrum - mean(comb_spectrum);
    % attenuate high frequencies:
    comb_spectrum(freqs>1000) = comb_spectrum(freqs>1000) ./ 10.^(log2(freqs(freqs>1000)/1000)*24/20);
end

function [spectrum] = hannwin_spectrum(angular_freq, specsize)
    % Spectrum of a hann window
    %
    % The hann window is a linear combination of modulated rectangular
    % windows r(n) = 1 for n=[0, N-1]:
    %
    % w(n) = 1/2*(1 - cos((2*pi*n)/(N-1)))
    %      = 1/2*r(n) - 1/4*exp(i*2*pi * n/(N-1))*r(n) - 1/4*exp(-i*2*pi * n/(N-1))*r(n)
    %
    % It's spectrum is then
    %
    % W(omega) = 1/2*R(omega) - 1/4*R(omega + (2*pi)/(N-1)) - 1/4*R(omega - (2*pi/(N-1)))
    %
    % with the spectrum of the rectangular window
    %
    % R(omega) = exp(-i*omega * (N-1)/2) * sin(N*omega/2) / sin(omega/2)
    %
    % (Source: https://en.wikipedia.org/wiki/Hann_function)
    %
    % angular_freq: Angular Frequency omega (0...2*pi), may be a vector.
    % specsize: Length N of the resulting spectrum
    %
    % Returns the spectral magnitude for angular_freq.

    function [spectrum] = rectwin_spectrum(angular_freq)
        % In case of angular_freq == 0, this will calculate NaN. This
        % will be corrected later.
        spectrum = ( exp(-j*angular_freq*(specsize-1)/2) .* ...
                     sin(specsize*angular_freq/2) ./ ...
                     sin(angular_freq/2) );
        spectrum(angular_freq == 0) = specsize;
    end

    delta_f = 2*pi / (specsize-1);
    spectrum = ( 1/2 * rectwin_spectrum(angular_freq) - ...
                 1/4 * rectwin_spectrum(angular_freq + delta_f) - ...
                 1/4 * rectwin_spectrum(angular_freq - delta_f) ) / specsize;
end

function [difference] = ifd_difference(signal, samplerate, blocksize, hopsize, basefrequencies)
    % Compare generated IF deviations with true IF deviation.
    %
    % Generate IF deviations at given base frequencies, then subtract
    % these from the true IF deviation of each block. The minimum
    % difference marks the base frequency of the block.
    %
    % Each difference is log-weighted along the frequency to account for
    % human perception, and bias-corrected along the base frequency to
    % compensate for higher variances at higher base frequencies.
    %
    % blocks: A SignalBlocks instance.
    % basefrequencies: an ordered vector of base frequencies in Hz.
    %
    % Returns a len(blocks) x len(basefrequencies) matrix of difference values.

    specsize = floor(blocksize/2) + 1;

    synthetic_ifds = synthetic_ifd(samplerate, specsize, basefrequencies);

    % weigh differences according to perception:
    freq = linspace(0, samplerate/2, specsize);
    log_f_weight = 1 ./ (samplerate/2).^(freq / (samplerate/2));
    speech_weight = ones(length(freq), 1);

    max_f0 = basefrequencies(end);
    ifds = ifd(signal, samplerate, blocksize, hopsize, max_f0);

    % larger base frequencies lead to larger IFDs:
    difference_bias = nanmean(abs(synthetic_ifds), 2);
    % larger signal variability leads to larger IFDs:
    signal_bias = nanmean(abs(ifds), 2);

    difference = zeros(numblocks(signal, blocksize, hopsize), length(basefrequencies));
    for idx=1:size(ifds, 1)
        this_ifd = ifds(idx, :);
        % the difference between the synthetic baseband instantaneous
        % frequency and the actual baseband instantaneous frequency
        % is minimal at the probable f0.
        difference_matrix = synthetic_ifds - this_ifd;
        bias = sqrt(difference_bias.^2 + signal_bias(idx, :).^2);
        % scale frequencies logarithmically, and correct for bias:
        difference(idx, :) = (nanmean(abs(difference_matrix) .* log_f_weight .* speech_weight', 2) ./ bias)';
    end
end

function [instfreq_deviation] = ifd(signal, samplerate, blocksize, hopsize, max_f0)
    % Instantaneous Frequency Deviation of a signal.
    %
    % Each time-frequency bin in the IFD has as value the difference
    % between the bin's frequency and the most prominent frequency track
    % in the bin's vicinity. As an example, if there is a prominent
    % frequency track 100 Hz above a bin, it's value will be 100. If a
    % bin is situated right on top of a frequency track, it's value will
    % be 0. If it is above a frequency track, it's value is negative.
    %
    % This is equivalent to the frequency differentiation of a
    % baseband-transformed STFT spectrum; aka BPD in [1].
    %
    % [1]: Krawczyk, M.; Gerkmann, T., "STFT Phase Reconstruction in
    %      Voiced Speech for an Improved Single-Channel Speech
    %      Enhancement," in Audio, Speech, and Language Processing,
    %      IEEE/ACM Transactions on , vol.22, no.12, pp.1931-1940, Dec.
    %      2014 doi: 10.1109/TASLP.2014.2354236
    %
    % The signal is cut into overlapping blocks. Each block is
    % differentiated by splitting it in two strongly-overlapping
    % sub-blocks, each sub-block is Fourier-transformed, the phase
    % spectra are calculated, and the difference between the
    % sub-blocks-spectra is calculated, and converted to frequencies.
    %
    % The frequencies thus obtained fall into a narrow range that is
    % limited by the sub-block overlap. The maximum visible difference
    % is given by max_f0, and governs the sub-block overlap. Lower
    % max_f0 decrease the sub-block overlap, and wrap IFD frequencies at
    % max_f0/2.
    %
    % blocks: A SignalBlocks instance.
    % max_f0: The max frequency distance visible in the IFD.
    %
    % Returns a len(blocks) x blocks.samplerate//2+1 matrix of IFD values.

    specsize = floor(blocksize/2) + 1;

    % number of samples that each sub-block-pair overlaps:
    dt = floor(samplerate/max_f0);
    % a time window for each sub-block that maximizes phase accuracy:
    window = hannpoisson(blocksize, [], false);

    longer_blocksize = blocksize+dt;
    instfreq = zeros(numblocks(signal, longer_blocksize, hopsize), specsize);
    for idx = 1:numblocks(signal, longer_blocksize, hopsize)
        block = signal(((idx-1)*hopsize+1):((idx-1)*hopsize+longer_blocksize));
        % differentiate the spectrum phase in the time direction:
        % split the block in two strongly-overlapping blocks, then calculate the
        % angle difference (aka instantaneous frequency).
        spectrum1 = fft(window'.*block(1:end-dt), blocksize);
        spectrum2 = fft(window'.*block(dt+1:end), blocksize);
        instfreq(idx, :) = spectrum2(1:specsize) .* conj(spectrum1(1:specsize));
    end

    % pure sinusoids change phase by this much in the given overlap time:
    baseband_phase_change = exp(j * linspace(0, pi, specsize) * dt);
    % calculate the baseband phase difference / instantaneous frequency deviation:
    instfreq_deviation = angle(instfreq .* conj(baseband_phase_change));
    instfreq_deviation = instfreq_deviation * max_f0/2 / pi; % display as frequencies
end

function [window] = hannpoisson(len, alpha, sym)
    % A window function with no side lobes.
    %
    % The Hann-Poisson window is a Hann Window times a Poisson window.
    % It has the unusual feature of having "no side lobes" in the sense
    % that, for alpha >= 2, the window-transform magnitude has negative
    % slope for all positive frequencies [1][2].
    %
    % [1]:
    % http://www.dsprelated.com/freebooks/sasp/Hann_Poisson_Window.html
    %
    % [2]: eq5.24 on p154 in Window Functions and Their Applications in
    % Signal Processing by K. M. M. Prabhu, CRC Press)
    %
    % length: The length of the window.
    % alpha: A slope constant (smooth slope for alpha >= 2)
    %
    % Returns the window array.

    if ~exist('alpha', 'var') || isempty(alpha)
        alpha = 2;
    end
    if ~exist('sym', 'var') || isempty(sym)
        sym = true;
    end

    if ~sym
        len = len + 1;
    end

    normalization = (len-1) / 2;
    n = -normalization:normalization;
    poisson = exp(-alpha * abs(n) / normalization);
    hann = 1/2 * (1+cos(pi*n / normalization));
    window = poisson .* hann;
    if ~sym
        window = window(1:end-1);
    end
end

function [synthetic_ifds] = synthetic_ifd(samplerate, specsize, basefrequencies)
    % Synthetic instfreq deviations of a range of tone complexes.
    %
    % samplerate: The sampling rate of the tone complexes.
    % specsize: The length of each IFD.
    % basefrequencies: An ordered vector of tone complex base
    %     frequencies in Hz.
    %
    % Returns a len(basefrequencies) x specsize matrix of IFDs.

    max_f0 = basefrequencies(end);
    synthetic_ifds = zeros(length(basefrequencies), specsize);
    freq = linspace(0, samplerate/2, specsize);
    for idx = 1:length(basefrequencies)
        basefrequency = basefrequencies(idx);
        % the baseband_instfreq shows the frequency difference to the closest
        % dominant harmonic:
        closest_harmonic = floor( (freq + basefrequency/2) / basefrequency );
        instfreq = closest_harmonic*basefrequency - freq;
        % since it is derived from the phase, it wraps:
        instfreq = wrap_angles(instfreq, max_f0/2);
        instfreq(freq < basefrequency/2) = 0;
        synthetic_ifds(idx, :) = instfreq;
    end
end

function [angles] = wrap_angles(angles, limit)
    % something like modulo, but wraps n ≷ limit to n ± 2*limit.
    %
    % Useful for confining angular values to a wrapping data range.
    %
    % angles: Some real values.
    % limit: The maximum/minimum valid value.
    % inplace: Whether to overwrite values in angles (faster).
    %
    % Returns wrapped angles.

    too_big = angles > limit;
    angles(too_big) = angles(too_big) - ceil((angles(too_big)-limit)/(2*limit))*2*limit;
    too_small = angles < -limit;
    angles(too_small) = angles(too_small) - floor((angles(too_small)+limit)/(2*limit))*2*limit;
end

function [prob] = value2posterior(mc, ifdd)
    % converts magnitude correlation and instantaneous frequency deviation to posterior probability
    %
    % All values were generated by pitch_probabilities.py
    % Data is encoded as base64 for binary correctness across text files.

    magnitude_correlation = ...
        ['mGrQJey9QMCe/0La69g3wBtUytH+ayzA8lEd3ktMEsBQBFrnZT8UQEqtaNaLZS1ANyySXLJV' ...
         'OEDkAPhmT/xAQKzrpp/FzUVAdNZV2DufSkA9wQQRsnBPQAPW2SQUIVJAZ0sxQc+JVEDMwIhd' ...
         'ivJWQDA24HlFW1lAlKs3lgDEW0D4II+yuyxeQC5Lc2e7SmBA4AWf9Rh/YUCSwMqDdrNiQA=='];
    magnitude_correlation = typecast(base64decode(magnitude_correlation), 'double');

    ifd_difference = ...
        ['QM+ExWYxsj9WRR487eOyP2y7t7JzlrM/gTFRKfpItD+Wp+qfgPu0P6wdhBYHrrU/wZMdjY1g' ...
         'tj/WCbcDFBO3P+x/UHqaxbc/Avbp8CB4uD8XbINnpyq5PyziHN4t3bk/Qli2VLSPuj9Xzk/L' ...
         'OkK7P2xE6UHB9Ls/grqCuEenvD+XMBwvzlm9P62mtaVUDL4/whxPHNu+vj/XkuiSYXG/Pw=='];
    ifd_difference = typecast(base64decode(ifd_difference), 'double');

    posterior = ...
        ['FAvfatwgsj+6t4fUv1awPzCfCvN3Mas/KLtCfM76pD9zS4ObSPmbP5IFkiZpzo4/BomuaYtqfD86' ...
         '0dREyE9mP8tCzWxIulA/zjotfjR2QT89OM4nQVg+P6JaRjc4O0A/Oq/XmvaJQT8b+sytyclCP/z3' ...
         'FqV6+EM/eqPgbJ5hRD8OTYuTNSlEP/LKwhzy90M/TnYNKhRyQz+Pcne44TlCP+yWorIkD7A/cGBM' ...
         'rDFcrT8QjswRE92oPxCxdo0eMaM/KARH6Dh9mT8uDmFlDjmMP32GGH63TXo/pj2yJW/OZD+q/PT8' ...
         'kitPP/iBYIoNOkA/O2hCRdghPD+EGiKIaSg+P5cMtd+tTEA/TFG0lN1HQT+QF4PSOxZCP3bT1HBg' ...
         'dkI/YlvTsad6Qj+cx1tSY1lCP6UrchjeqEE/pL6/hZluQD/MxcByV52tP94d9TS6sqo/W7fTTjlU' ...
         'pj+eVBTCCCShP+ZOfOKV7pY/Ok2O/B/IiT9iKXgBgVx4P/jpwDYAlGM/3KHHceL3TT+DkEXkXbU/' ...
         'P7hCRDQ5Wzs/tmiLVlwBPT9wIYgcPBM/P/jDlc6oPkA/ukLRGpy/QD+dXIP9cAlBP4BPSITsGEE/' ...
         'a8q7o5ftQD/jektu+UFAPzL7pKLQnT4/UMlq2evCuT9qbMhjaQK3P4qHixXe4LI/UgpxGruorD9i' ...
         'mvHRO0KjP5M0f3mu75U/o6uK6kP7hD9F5CyBdh9xP3ia9H+LFVs//COlP3JrTT90jx86NyRJP1CW' ...
         'twsuCEo/JL7i2GZdSz9oQDeDiRhMPwpnqnaobkw/xRv4lOd7TD+QroCy7EJMP2ryWaP5rEs/fP/K' ...
         'wSt/Sj/a97CXYjxJP971ghvX39M/qsNy4wnP0T8sMXxhQzjNPyJhtqNvC8Y/E9ZxtDR3vT9wO6Vi' ...
         'D8iwP95ppFmUC6A/uHGoaEcYij8bQzvzBqt0P1D9qcpTfGY/uqRUZqINYz9WRqZs/IVjP0SdemZw' ...
         'a2Q/OgTNp3r0ZD843ytYES9lP0owKXorMGU/BM/eTyX+ZD/K8uJnloVkP6L0veWBn2M/3hlsRfOr' ...
         'Yj+iYgdleh7qPwJA+O76P+c/tA4Z9Ev+4j8Yp4P6NWDcP5tfidjLxNI/trzXPshRxT+y4KsEI1m0' ...
         'P9Q9MgNyfqA/UvNmR54Iij/yt3oVPCV8PwT9PcqvvHc/CkVZezBCeD/dG7YN1WR5P26cY8ilP3o/' ...
         'IN/qDuTKej/UrtROGhR7P/ah3GIU63o/+jCKZ1w4ej9a35Pd1A15Pywc3Qws5Hc/V6rT2A6P+T+I' ...
         '8E/N0Y/2P9ueKjokcPI/JlxeATKQ6z/I9H8RHSLiP4ZTP7m0itQ/jHRiBMiZwz+sh23mM7yvPx3f' ...
         'Q1rd85g/cEEexpe9ij9GxwE2xoOGP7ahBsXsCYc/b7cSG3IniD/Zujo/BCOJPwfZrvOi5ok/vNFP' ...
         '84Juij+jSjpUmmOKP1FJaXjFrYk/Kv86FmKNiD/WwDNWyH+HP8nmHMCeUARAdJ0dqSWtAUAcK8qO' ...
         'MfT8P5BvmE9eDPY/+v7bA2cn7T/tWMDA12TgP+YPRh7eMM8/VBHFq2I9uT9aNmtK+sKjP6gPVKVK' ...
         'BpU/wNIzBsu8kT8czwUuuzuSP0CCkM61MpM/SC5UjZMVlD8JP84QjMCUPwEmOmjmMJU/hlq0FmQ6' ...
         'lT+g2AmmQL2UP7Kx1L+95ZM/2zJf9IUukz9SRon7gqgNQEfapATtDQlAqZVcC7tHBEAUU4/g7lz/' ...
         'PyLexlZ20fQ/bB3nK70n5z9VkdzsEObVP8QNnJYsvME/YYqoz5mxqz9AD7ZNcE6dPybnhmIuypg/' ...
         '2Gv5Z7KymT+LKRIeSUabP1qpBicNqJw/2jGetjSgnT8IlAgIPhqePwivbuxrFJ4/XBazXNlsnT+f' ...
         '1ek73EWcPygYaNpvWJs/LKSfrNKeE0BKOZw9QXkQQJnl6aCEqgpAzo8l+yPCBEDsE5OHjJr7PzJH' ...
         '79IEkO4/mgH3wJzi3D8BgspRp3PHP3aRKgdLLbI/Sn/hDG4Aoz/xR3pVRQugP0KEiiYHx6A/5aHm' ...
         'MF79oT9fxQ6wdwijP/bh/NOlvaM/AtQ+J3/zoz+oWistMcujPzNayU7lY6M/izTahXivoj/+orur' ...
         'bu6hP6RLib4uDRdAMDOn6QlkFEB+6LFjYgIRQDnd4X7SeApAxAf5NQm+AUCWfFlYyMbzP6Q6nNVX' ...
         'puI/DqcWqpIrzj8YK8hdlwG3PwBhXtuZk6c/XuYYZeHHoz/MSewtfMWkP1+BJWlkiKY/wDwqezAa' ...
         'qD/+T49t0hipP0xZgYgUVak/9gtDKUQSqT9uR8k78qyoP1bY2E1g9qc/x9dbPyjmpj8Y4PJpnBsa' ...
         'QPmBx+QH7hhA1q5zT6WSFUAeFBI1X5YQQGpjITFgdAZA/uu6sOdc+T8QtHTpKXvnP/hdap16hdI/' ...
         'vr6lxNXWuz8YD2iLc0ysP2npu+5Fpqc/4titT7DiqD+iXefg3yqrP0Eiz1QEMq0/ugPQQNeFrj9g' ...
         'a7QisASvP3x7eyow4q4/LCLR29GQrj/+k7NRRvatPyL4vQ3786w/nvlc4IiwG0Bv/nX6Z1kcQCT4' ...
         'OC4O7BlAlfyS8GtKFEBmqOtXLOELQLEC1qe0q/8/pYqQRBSD7D9yccUiyaDVP1Rm0xEHEcA/Tqn3' ...
         'OHeFsD992VDk18mrP57dz2nqT60/7EoCrwHqrz+RRrEMCw2xPx5VfX3N17E/zmvVSdxLsj9QgToX' ...
         'uGOyP/I33nfdS7I/UCVTfxoUsj8ULkl1GcSxPxCogC4xsxxAmjbKIr1cHkA8YTkqNkgdQOWobM2W' ...
         'ChhAuiyIGJ7JEEAEsLl4TfkCQFB/3Kjk2vA/8NU9qxzy2D9cmuz1wSPCP23t1g9EsrI/SP/uG9e+' ...
         'rz8EnEuIf8+wP4AnvKzwPbI/hd/tnAlesz+Whewxn0O0PwQguLUH/rQ/sHlXCe1QtT8NsxZALzu1' ...
         'P7YomyaD6bQ/W19cdryWtD8arxwpwlQfQGZsnVaE1CBAZxLm2+66IEC0xtBweWccQGA1ksFqrBNA' ...
         'qjfBGN+6BUDtH2jsBlbzP8r9Pqqno9w/on5LXtdlxD9zHDz6PqS0P7zo01aJlLE/QuuT9ruysj+t' ...
         'rVbQA0W0P+TYwcRHd7U/osk6yhuGtj9wZsmA25O3PxYDxHsTRbg/xAG10jF3uD+IDP3iOP63P7LM' ...
         '72CLNrc/8bGoqwSmIEBzKmwME+giQJJGBgxjPiNALrpTYnCDIEC8ou1k0FsWQLtX2BIovwdAOMQK' ...
         'Tjj79D9xUwBa+WLfPzQiE2L1YsY/8lqN89uItj94iqg4iDGzP654A3HXdLQ/xXpfd3wztj9ldleI' ...
         'w323P8Q7ykyTvbg/bQtjVgH5uT9GBtZijQq7P340VzKrDrw/juyd9d/2uz/+rqayrOW6PxoqIMh4' ...
         '0iJAiMhHxw45JUAwubZgchMlQG18idYAzSFA7qOP3AEMGEBneUjtfGgJQLjb7zb+IPY/BEw6oOdc' ...
         '4D/YiM77w6PHPwXSUPgSarg/LoDjuW7StD8ErHYgYwy2P406PPr0D7g/kMqCjT2duT9DXk/j8vu6' ...
         'P2B+Pm0kHrw/KrXLmWKXvT/uN8o6+NC/P2B2tpSwb8A/E6hHWe5TwD+ks881EFQmQIwtoexmDihA' ...
         'ekpW8XumJkACH9RGvioiQF5PGPqF1BhAeEPzCOeQC0D0kHNfxiz4P1Qu+z1GXuE/4oztxQzTyD9u' ...
         'D+v7UxC6P+X6hqgyc7Y/rjKiF66qtz/nXAPRRt+5PxzyUekLtLs/tvRW+xsavT+XXjpCqDy+P1Ss' ...
         'RbmAN8A/pPum0f2jwT+Km8SbbkTCP6x1g356gsI//XZXbSPFKkDkfm2dqrkrQByPRPenPihAx/qc' ...
         'op2YIUCox6GHY+8XQMAJzoOCkgxA5i0KkJEj+j/yc03/tNviP/ZpefEJk8o/Kv06pvCauz+yF1ov' ...
         'jge4PyaphgqBuLk/wG3KqGMWvD+8vRKDo6i9P3hgJ0x+xL4/lt/8LNcqwD8WqAUWxdzBP7LYiHxs' ...
         'JcM/VKMAp83Nwj+68dVWpnDCP5RhejaQqDBApFWXOEWHMECsmG48xj4rQHDDR7P45yFABlmUG9aQ' ...
         'F0CpzOYzHl8NQA6RgStXgPs/eMUCr+0B5D8qLgyflzfMP8iIzNv9/7w/OHNobvskuT/yeh/nnYK7' ...
         'P2E310u3N74/HPP3O80pvz/evsf1CJi/P4ml5gnb1cA/oxwI+UtFwz+5hNDy43rEP77MvtLKAsM/' ...
         'RJpAgj24wT8='];
    posterior = reshape(typecast(base64decode(posterior), 'double'), [20,20]);

    prob = interp2(magnitude_correlation, ifd_difference, posterior, mc, ifdd, '*spline');
end

function bytes = base64decode(base64)
    % convert a base64-encoded string to uint8

    % strip line breaks
    base64 = strrep(base64, sprintf('\n'), '');
    % add padding if missing
    if mod(length(base64), 4) ~= 0
        base64 = [base64 repmat('=', [1 mod(length(base64), 4)])];
    end
    % remember padding
    padding = sum(base64 == '=');
    % convert from string representation to base64 bytes
    table = zeros(1, 128, 'uint8');
    table(['A':'Z' 'a':'z' '0':'9' '+' '/']+0) = 0:63;
    base64 = table(base64');

    % convert every four base64 bytes to three uint8 bytes
    bytes = zeros(length(base64)/4*3, 1, 'uint8');
    bytes(1:3:end) = bitor(bitshift(base64(1:4:end), 2), ...
                           bitshift(base64(2:4:end), -4));
    bytes(2:3:end) = bitor(bitshift(bitand(base64(2:4:end), 15), 4), ... % four LSB
                           bitshift(base64(3:4:end), -2));
    bytes(3:3:end) = bitor(bitshift(bitand(base64(3:4:end), 3), 6), ... % two LSB
                           base64(4:4:end));
    % strip padding
    bytes = bytes(1:end-padding);
end

function [time, freq, prob] = max_track_viterbi(probabilities, delta_t, basefrequencies)
    % Extract frequency track from probabilities.
    %
    % The frequency track is the track of maximum probability with a
    % minimum of frequency steps. Frequency steps are penalized
    % propertionally to the multiplicative step size.
    %
    % probabilities: a time x frequency matrix of pitch probabilities.
    % frequency: The second axis of probabilities in Hz.
    %
    % Returns each frame's most likely frequency, and the frequency's
    % probability.

    time = (0:size(probabilities, 1)-1)*delta_t;

    % transition probability between two frequencies is the quotient
    % between those frequencies, normalized to < 1.
    transition = basefrequencies' ./ basefrequencies;
    transition(transition>1) = 1 ./ transition(transition>1);

    % accumulate probabilities for each time step and select the
    % highest cumulative probability path per basefrequencies and time:
    cum_probability = probabilities;
    % save step that lead to this time/basefrequencies:
    idx_probability = zeros(size(probabilities, 1), size(probabilities, 2));
    for idx=2:size(probabilities, 1)
        step_probs = cum_probability(idx-1, :)' .* cum_probability(idx, :) .* transition;
        [~, max_prob_idx] = max(step_probs, [], 1);
        idx_probability(idx, :) = max_prob_idx;
        cum_probability(idx, :) = step_probs(sub2ind(size(step_probs), max_prob_idx, 1:length(basefrequencies)));
        % normalize, so large products of small numbers don't end up zero:
        cum_probability(idx, :) = cum_probability(idx, :) / mean(cum_probability(idx, :));
    end

    % walk backwards and select the path that lead to the maximum
    % cumulative probability. For each step in the path, extract the
    % basefrequencies and the local (non-cumulative) probability:
    freq = zeros(1, size(probabilities, 1));
    prob = zeros(1, size(probabilities, 1));
    [~, f_idx] = max(cum_probability(end, :));
    for t_idx=size(probabilities, 1):-1:1
        freq(t_idx) = basefrequencies(f_idx);
        prob(t_idx) = probabilities(t_idx, f_idx);
        f_idx = idx_probability(t_idx, f_idx);
    end
end
