import scipy.signal
import numpy
from base64 import b64decode
from scipy.interpolate import RectBivariateSpline

_basefrequencies=numpy.logspace(numpy.log10(50), numpy.log10(450), 200)

def estimate_f0(blocks, basefrequencies=_basefrequencies):
    """A base frequency estimation in the magnitude and phase spectrum.

    MaPS-f0 uses a signal's magnitude STFT and it's IF deviation to
    estimate the maximum-likelihood base frequency of a tone complex
    for a given series of base frequencies.

    blocks: A SignalBlocks instance.
    basefrequencies: an ordered vector of base frequencies in Hz.

    Returns a vector of times, a vector of likely base frequencies,
    and a vector of their likelihood.

    """
    probability = pitch_probability(blocks, basefrequencies)
    return max_track_viterbi(probability, blocks.hopsize/blocks.samplerate, basefrequencies)


def pitch_probability(blocks, basefrequencies=_basefrequencies):
    """A probability-of-pitch estimation in the magnitude and phase spectrum.

    MaPS-f0 uses a signal's magnitude STFT and it's IF deviation to
    estimate the likelihood that a given range of base frequencies are
    base frequencies of a tone complex.

    blocks: A SignalBlocks instance.
    basefrequencies: an ordered vector of base frequencies in Hz.

    Returns a len(blocks) x len(basefrequencies) matrix of likelihood values.

    """

    mc = magnitude_correlation(blocks, basefrequencies)
    ifdd = ifd_difference(blocks, basefrequencies)
    return value2posterior(mc, ifdd)


def max_track(probabilities, delta_t, basefrequencies=_basefrequencies):
    """Extract frequency track from probabilities.

    The frequency track is located at the maximum likelihood per time
    frame.

    probabilities: a time x frequency matrix of pitch probabilities.
    delta_t: the time distance between two likelihood frames.
    basefrequencies: The second axis of probabilities in Hz.

    Returns the time of each frame, it's most likely frequency, and
    the frequency's likelihood.

    """
    time = numpy.arange(len(probabilities))*delta_t
    freq = basefrequencies[numpy.argmax(probabilities, axis=1)]
    prob = numpy.max(probabilities, axis=1)
    return time, freq, prob


def max_track_viterbi(probabilities, delta_t, basefrequencies):
    """Extract frequency track from probabilities.

    The frequency track is the track of maximum probability with a
    minimum of frequency steps. Frequency steps are penalized
    propertionally to the multiplicative step size.

    probabilities: a time x frequency matrix of pitch probabilities.
    frequency: The second axis of probabilities in Hz.

    Returns each frame's most likely frequency, and the frequency's
    probability.

    """

    time = numpy.arange(len(probabilities))*delta_t

    # transition probability between two frequencies is the quotient
    # between those frequencies, normalized to < 1.
    transition = basefrequencies[:, None] / basefrequencies[None, :]
    transition[transition>1] = 1/transition[transition>1]

    # accumulate probabilities for each time step and select the
    # highest cumulative probability path per basefrequencies and time:
    cum_probability = probabilities.copy()
    # save step that lead to this time/basefrequencies:
    idx_probability = numpy.empty(probabilities.shape, dtype=int)
    for idx in range(1, len(probabilities)):
        step_probs = cum_probability[idx-1][:,None] * probabilities[idx][None,:] * transition
        max_prob_idx = step_probs.argmax(axis=0)
        idx_probability[idx] = max_prob_idx
        cum_probability[idx] = step_probs[[max_prob_idx, numpy.arange(len(basefrequencies))]]
        # normalize, so large products of small numbers don't end up zero
        cum_probability[idx] /= cum_probability[idx].mean()

    # walk backwards and select the path that lead to the maximum
    # cumulative probability. For each step in the path, extract the
    # basefrequencies and the local (non-cumulative) probability:
    freq = numpy.empty(len(probabilities))
    prob = numpy.empty(len(probabilities))
    f_idx = numpy.argmax(cum_probability[-1])
    for t_idx in reversed(range(len(probabilities))):
        freq[t_idx] = basefrequencies[f_idx]
        prob[t_idx] = probabilities[t_idx, f_idx]
        f_idx = idx_probability[t_idx, f_idx]

    return time, freq, prob


def magnitude_correlation(blocks, basefrequencies=_basefrequencies):
    """Correlate synthetic tone complex spectra with a true spectrum.

    Generate spectra at a number of given base frequencies, then
    correlate each of these spectra with the magnitude signal
    STFT-spectra.

    Before correlation, each frequency bin is log-weighted to make the
    correlation perceptually accurate.

    blocks: A SignalBlocks instance.
    basefrequencies: An ordered list of base frequencies in Hz.

    Returns a len(blocks) x len(basefrequencies) matrix of correlation values.

    """
    specsize = blocks.blocksize//2+1

    # weigh differences according to perception:
    f = numpy.linspace(0, blocks.samplerate/2, specsize)
    log_f_weight = 1 / (blocks.samplerate/2)**(f / (blocks.samplerate/2))

    correlation = numpy.zeros([len(blocks), len(basefrequencies)])
    synthetic_magnitudes = synthetic_magnitude(blocks.samplerate, specsize, basefrequencies)

    for idx, spectrum in enumerate(stft(blocks)):
        # the correlation for real signals does not require the conj():
        correlation[idx] = numpy.sum(numpy.abs(spectrum) *
                                     synthetic_magnitudes *
                                     log_f_weight, axis=1)

    return correlation


def stft(blocks, *, nfft=None, windowfunc=scipy.signal.hann):
    """Short-time Fourier Transform of a signal.

    The signal is cut into short overlapping blocks, and each block is
    transformed into the frequency domain using the FFT. Before
    transformation, each block is windowed by windowfunc.

    blocks: A SignalBlocks instance.
    nfft: None for blocksize, or a number.
    windowfunc: A function that returns a window.

    Returns a complex spectrum.

    """
    window = windowfunc(blocks.blocksize)
    nfft = nfft or blocks.blocksize
    specsize = nfft//2+1
    for idx, block in enumerate(blocks):
        yield numpy.fft.rfft(block*window, nfft)


def synthetic_magnitude(samplerate, specsize, basefrequencies=_basefrequencies):
    """Synthetic magnitude spectra of a range of tone complexes.

    samplerate: The sampling rate of the tone complexes.
    specsize: The length of each spectrum.
    basefrequencies: An ordered vector of tone complex base
        frequencies in Hz.

    Returns a len(basefrequencies) x specsize matrix of tone complex spectra.

    """
    freqs = numpy.linspace(0, samplerate/2, specsize)
    synthetic_spectra = numpy.empty((len(basefrequencies), specsize), numpy.float64)
    for idx, basefrequency in enumerate(basefrequencies):
        synthetic_spectra[idx, :] = hannwin_comb(samplerate, basefrequency, specsize)
    return synthetic_spectra


def hannwin_comb(samplerate, basefreq, specsize):
    """Approximate a speech-like correlation spectrum of a tone complex.

    This is an approximation of time_domain_comb that runs much
    faster.

    Instead of calculating the FFT of a series of hann-windowed
    sinuses, this models the spectrum of a tone-complex as a series of
    hann-window-spectrums.

    For a perfect reconstruction, this would need to calculate the sum
    of many hann-window-spectra. Since hann window spectra are very
    narrow, this assumes that each window spectrum extends from
    n*basefreq-basefreq/2 to n*basefreq+basefreq/2 and that
    neighboring spectra do not influence each other.

    This assumtion holds as long as basefreq >> 1/specsize.

    Amplitudes are normalized by specsize.

    To make the spectrum more speech-like, frequencies above 1000 Hz
    are attenuated by 24 dB/oct.

    To make the correlation of this spectrum and some other spectrum
    have a normalized gain, the spectrum is shifted to be zero-mean.

    samplerate: The sampling rate in Hz of the signal.
    basefreq: The base frequency in Hz of the tone complex.
    specsize: The length of the resulting spectrum in bins
              (typically 2**N+1 for type(N) == int).

    Returns a real magnitude spectrum.

    """

    freqs = numpy.linspace(0, samplerate/2, specsize)
    # create a local frequency vector around each harmonic, going from
    # -basefreq/2 to basefreq/2 within the area around the nth
    # harmonic n*basefreq-basefreq/2 to n*basefreq+basefreq/2:
    closest_harmonic = (freqs + basefreq/2) // basefreq
    # ignore first half-wave:
    closest_harmonic[closest_harmonic==0] = 1
    local_frequency = closest_harmonic*basefreq - freqs
    # convert from absolute frequency to angular frequency:
    local_angular_freq = local_frequency/(samplerate/2)*2*numpy.pi
    # evaluate hannwin_spectrum at the local frequency vector:
    comb_spectrum = numpy.abs(hannwin_spectrum(local_angular_freq, specsize))
    # normalize to zero mean:
    comb_spectrum -= numpy.mean(comb_spectrum)
    # attenuate high frequencies:
    comb_spectrum[freqs>1000] /= 10**(numpy.log2(freqs[freqs>1000]/1000)*24/20)
    return comb_spectrum


def hannwin_spectrum(angular_freq, specsize):
    """Spectrum of a hann window

    The hann window is a linear combination of modulated rectangular
    windows r(n) = 1 for n=[0, N-1]:

    w(n) = 1/2*(1 - cos((2*pi*n)/(N-1)))
         = 1/2*r(n) - 1/4*exp(i*2*pi * n/(N-1))*r(n) - 1/4*exp(-i*2*pi * n/(N-1))*r(n)

    It's spectrum is then

    W(omega) = 1/2*R(omega) - 1/4*R(omega + (2*pi)/(N-1)) - 1/4*R(omega - (2*pi/(N-1)))

    with the spectrum of the rectangular window

    R(omega) = exp(-i*omega * (N-1)/2) * sin(N*omega/2) / sin(omega/2)

    (Source: https://en.wikipedia.org/wiki/Hann_function)

    angular_freq: Angular Frequency omega (0...2*pi), may be a vector.
    specsize: Length N of the resulting spectrum

    Returns the spectral magnitude for angular_freq.

    """

    def rectwin_spectrum(angular_freq):
        # In case of angular_freq == 0, this will calculate NaN. This
        # will be corrected later.
        spectrum = ( numpy.exp(-1j*angular_freq*(specsize-1)/2) *
                     numpy.sin(specsize*angular_freq/2) /
                     numpy.sin(angular_freq/2) )
        # since sin(x) == x for small x, the above expression
        # evaluates to specsize for angular_freq == 0.
        spectrum[angular_freq == 0.0] = specsize
        return spectrum

    angular_freq = numpy.asarray(angular_freq, dtype='float64')
    delta_f = 2*numpy.pi / (specsize-1)
    # don't warn about division by zero, NaNs will be corrected.
    with numpy.errstate(invalid='ignore'):
        return (1/2 * rectwin_spectrum(angular_freq) -
                1/4 * rectwin_spectrum(angular_freq + delta_f) -
                1/4 * rectwin_spectrum(angular_freq - delta_f)) / specsize


def ifd_difference(blocks, basefrequencies=_basefrequencies):
    """Compare generated IF deviations with true IF deviation.

    Generate IF deviations at given base frequencies, then subtract
    these from the true IF deviation of each block. The minimum
    difference marks the base frequency of the block.

    Each difference is log-weighted along the frequency to account for
    human perception, and bias-corrected along the base frequency to
    compensate for higher variances at higher base frequencies.

    blocks: A SignalBlocks instance.
    basefrequencies: an ordered vector of base frequencies in Hz.

    Returns a len(blocks) x len(basefrequencies) matrix of difference values.

    """
    specsize = blocks.blocksize//2+1

    synthetic_ifds = synthetic_ifd(blocks.samplerate, specsize, basefrequencies)

    # weigh differences according to perception:
    f = numpy.linspace(0, blocks.samplerate/2, specsize)
    log_f_weight = 1 / (blocks.samplerate/2)**(f / (blocks.samplerate/2))
    speech_weight = numpy.ones(f.shape)

    max_f0 = basefrequencies[-1]
    ifds = ifd(blocks, max_f0=max_f0)

    # larger base frequencies lead to larger IFDs:
    difference_bias = numpy.nanmean(numpy.abs(synthetic_ifds), axis=1)
    # larger signal variability leads to larger IFDs:
    signal_bias = numpy.nanmean(numpy.abs(ifds), axis=1)

    difference = numpy.zeros([len(blocks), len(basefrequencies)])
    for idx, this_ifd in enumerate(ifds):
        # the difference between the synthetic baseband instantaneous
        # frequency and the actual baseband instantaneous frequency
        # is minimal at the probable f0.
        difference_matrix = synthetic_ifds - this_ifd
        bias = numpy.sqrt(difference_bias**2 + signal_bias[idx]**2)
        # scale frequencies logarithmically, and correct for bias:
        difference[idx] = numpy.nanmean(numpy.abs(difference_matrix) *
                                        log_f_weight *
                                        speech_weight, axis=1) / bias

    return difference


def ifd(blocks, *, max_f0=450):
    """Instantaneous Frequency Deviation of a signal.

    Each time-frequency bin in the IFD has as value the difference
    between the bin's frequency and the most prominent frequency track
    in the bin's vicinity. As an example, if there is a prominent
    frequency track 100 Hz above a bin, it's value will be 100. If a
    bin is situated right on top of a frequency track, it's value will
    be 0. If it is above a frequency track, it's value is negative.

    This is equivalent to the frequency differentiation of a
    baseband-transformed STFT spectrum; aka BPD in [1].

    [1]: Krawczyk, M.; Gerkmann, T., "STFT Phase Reconstruction in
         Voiced Speech for an Improved Single-Channel Speech
         Enhancement," in Audio, Speech, and Language Processing,
         IEEE/ACM Transactions on , vol.22, no.12, pp.1931-1940, Dec.
         2014 doi: 10.1109/TASLP.2014.2354236

    The signal is cut into overlapping blocks. Each block is
    differentiated by splitting it in two strongly-overlapping
    sub-blocks, each sub-block is Fourier-transformed, the phase
    spectra are calculated, and the difference between the
    sub-blocks-spectra is calculated, and converted to frequencies.

    The frequencies thus obtained fall into a narrow range that is
    limited by the sub-block overlap. The maximum visible difference
    is given by max_f0, and governs the sub-block overlap. Lower
    max_f0 decrease the sub-block overlap, and wrap IFD frequencies at
    max_f0/2.

    blocks: A SignalBlocks instance.
    max_f0: The max frequency distance visible in the IFD.

    Returns a len(blocks) x blocks.samplerate//2+1 matrix of IFD values.

    """
    specsize = blocks.blocksize//2+1

    # number of samples that each sub-block-pair overlaps:
    dt = int(blocks.samplerate//max_f0)
    # a time window for each sub-block that maximizes phase accuracy:
    window = hannpoisson(blocks.blocksize, sym=False)

    longer_blocks = SignalBlocks(blocks.data, blocks.samplerate,
                                 blocks.blocksize+dt, blocks.hopsize)
    instfreq = numpy.zeros([len(longer_blocks), specsize], dtype='complex')
    for idx, block in enumerate(longer_blocks):
        # differentiate the spectrum phase in the time direction:
        # split the block in two strongly-overlapping blocks, then calculate the
        # angle difference (aka instantaneous frequency).
        spectrum1 = numpy.fft.rfft(window*block[:-dt], n=blocks.blocksize)
        spectrum2 = numpy.fft.rfft(window*block[dt:], n=blocks.blocksize)
        instfreq[idx] = spectrum2 * spectrum1.conj()

    # pure sinusoids change phase by this much in the given overlap time:
    baseband_phase_change = numpy.exp(1j * numpy.linspace(0, numpy.pi, specsize) * dt)
    # calculate the baseband phase difference / instantaneous frequency deviation:
    instfreq_deviation = numpy.angle(instfreq * baseband_phase_change.conj())
    instfreq_deviation *= max_f0/2 / numpy.pi # display as frequencies

    return instfreq_deviation


def hannpoisson(length, *, alpha=2, sym=True):
    """A window function with no side lobes.

    The Hann-Poisson window is a Hann Window times a Poisson window.
    It has the unusual feature of having "no side lobes" in the sense
    that, for alpha >= 2, the window-transform magnitude has negative
    slope for all positive frequencies [1][2].

    [1]:
    http://www.dsprelated.com/freebooks/sasp/Hann_Poisson_Window.html

    [2]: eq5.24 on p154 in Window Functions and Their Applications in
    Signal Processing by K. M. M. Prabhu, CRC Press)

    length: The length of the window.
    alpha: A slope constant (smooth slope for alpha >= 2)

    Returns the window array.

    """
    if not sym:
        length += 1
    normalization = (length-1) / 2
    n = numpy.arange(-normalization, normalization+1)
    poisson = numpy.exp(-alpha*numpy.abs(n) / normalization)
    hann = 1/2 * (1+numpy.cos(numpy.pi*n / normalization))
    if sym:
        return poisson * hann
    else:
        return (poisson * hann)[:-1]


def synthetic_ifd(samplerate, specsize, basefrequencies=_basefrequencies):
    """Synthetic instfreq deviations of a range of tone complexes.

    samplerate: The sampling rate of the tone complexes.
    specsize: The length of each IFD.
    basefrequencies: An ordered vector of tone complex base
        frequencies in Hz.

    Returns a len(basefrequencies) x specsize matrix of IFDs.

    """
    max_f0 = basefrequencies[-1]
    synthetic_ifds = numpy.zeros((len(basefrequencies), specsize), numpy.float64)
    f = numpy.linspace(0, samplerate/2, specsize)
    for idx, basefrequency in enumerate(basefrequencies):
        # the baseband_instreq shows the frequency difference to the closest
        # dominant harmonic:
        closest_harmonic = (f + basefrequency/2) // basefrequency
        instfreq = closest_harmonic*basefrequency - f
        # since it is derived from the phase, it wraps:
        instfreq = wrap_angles(instfreq, max_f0/2, inplace=True)
        instfreq[f < basefrequency/2] = 0
        synthetic_ifds[idx] = instfreq

    return synthetic_ifds


def wrap_angles(angles, limit, *, inplace=False):
    """something like modulo, but wraps n ≷ limit to n ± 2*limit.

    Useful for confining angular values to a wrapping data range.

    angles: Some real values.
    limit: The maximum/minimum valid value.
    inplace: Whether to overwrite values in angles (faster).

    Returns wrapped angles.

    """
    angles = numpy.array(angles, copy=not inplace)
    too_big = angles > limit
    angles[too_big] -= numpy.ceil((angles[too_big]-limit)/(2*limit))*2*limit
    too_small = angles < -limit
    angles[too_small] -= numpy.floor((angles[too_small]+limit)/(2*limit))*2*limit
    return angles


class SignalBlocks:
    """A generator for short, overlapping signal blocks.

    Each SignalBlocks instance contains the signal `data` and its
    `samplerate`. It generates short, possibly overlapping signal
    blocks of length `blocksize`. Each block starts `hopsize` after
    the previous block.

    The SignalBlocks' `len` is its number of blocks, and its
    `duration` is the signal length in seconds.

    """

    def __init__(self, data, samplerate, blocksize=2048, hopsize=1024):
        self.data = data
        self.samplerate = int(samplerate)
        self.blocksize = int(blocksize)
        self.hopsize = int(hopsize)

    def __iter__(self):
        idx = 0
        while idx+self.blocksize < len(self.data):
            yield(self.data[idx:idx+self.blocksize])
            idx += self.hopsize

    def __len__(self):
        return int(numpy.ceil( (len(self.data)-self.blocksize) / self.hopsize ))

    @property
    def duration(self):
        return len(self.data)/self.samplerate

_magnitude_correlation = numpy.frombuffer(b64decode(
    b'mGrQJey9QMCe/0La69g3wBtUytH+ayzA8lEd3ktMEsBQBFrnZT8UQEqtaNaLZS1ANyySXLJVOEDk'
    b'APhmT/xAQKzrpp/FzUVAdNZV2DufSkA9wQQRsnBPQAPW2SQUIVJAZ0sxQc+JVEDMwIhdivJWQDA2'
    b'4HlFW1lAlKs3lgDEW0D4II+yuyxeQC5Lc2e7SmBA4AWf9Rh/YUCSwMqDdrNiQA=='
    ), dtype='double')

_ifd_difference = numpy.frombuffer(b64decode(
    b'QM+ExWYxsj9WRR487eOyP2y7t7JzlrM/gTFRKfpItD+Wp+qfgPu0P6wdhBYHrrU/wZMdjY1gtj/W'
    b'CbcDFBO3P+x/UHqaxbc/Avbp8CB4uD8XbINnpyq5PyziHN4t3bk/Qli2VLSPuj9Xzk/LOkK7P2xE'
    b'6UHB9Ls/grqCuEenvD+XMBwvzlm9P62mtaVUDL4/whxPHNu+vj/XkuiSYXG/Pw=='
    ), dtype='double')

_posterior = numpy.frombuffer(b64decode(
    b'FAvfatwgsj+6t4fUv1awPzCfCvN3Mas/KLtCfM76pD9zS4ObSPmbP5IFkiZpzo4/BomuaYtqfD86'
    b'0dREyE9mP8tCzWxIulA/zjotfjR2QT89OM4nQVg+P6JaRjc4O0A/Oq/XmvaJQT8b+sytyclCP/z3'
    b'FqV6+EM/eqPgbJ5hRD8OTYuTNSlEP/LKwhzy90M/TnYNKhRyQz+Pcne44TlCP+yWorIkD7A/cGBM'
    b'rDFcrT8QjswRE92oPxCxdo0eMaM/KARH6Dh9mT8uDmFlDjmMP32GGH63TXo/pj2yJW/OZD+q/PT8'
    b'kitPP/iBYIoNOkA/O2hCRdghPD+EGiKIaSg+P5cMtd+tTEA/TFG0lN1HQT+QF4PSOxZCP3bT1HBg'
    b'dkI/YlvTsad6Qj+cx1tSY1lCP6UrchjeqEE/pL6/hZluQD/MxcByV52tP94d9TS6sqo/W7fTTjlU'
    b'pj+eVBTCCCShP+ZOfOKV7pY/Ok2O/B/IiT9iKXgBgVx4P/jpwDYAlGM/3KHHceL3TT+DkEXkXbU/'
    b'P7hCRDQ5Wzs/tmiLVlwBPT9wIYgcPBM/P/jDlc6oPkA/ukLRGpy/QD+dXIP9cAlBP4BPSITsGEE/'
    b'a8q7o5ftQD/jektu+UFAPzL7pKLQnT4/UMlq2evCuT9qbMhjaQK3P4qHixXe4LI/UgpxGruorD9i'
    b'mvHRO0KjP5M0f3mu75U/o6uK6kP7hD9F5CyBdh9xP3ia9H+LFVs//COlP3JrTT90jx86NyRJP1CW'
    b'twsuCEo/JL7i2GZdSz9oQDeDiRhMPwpnqnaobkw/xRv4lOd7TD+QroCy7EJMP2ryWaP5rEs/fP/K'
    b'wSt/Sj/a97CXYjxJP971ghvX39M/qsNy4wnP0T8sMXxhQzjNPyJhtqNvC8Y/E9ZxtDR3vT9wO6Vi'
    b'D8iwP95ppFmUC6A/uHGoaEcYij8bQzvzBqt0P1D9qcpTfGY/uqRUZqINYz9WRqZs/IVjP0SdemZw'
    b'a2Q/OgTNp3r0ZD843ytYES9lP0owKXorMGU/BM/eTyX+ZD/K8uJnloVkP6L0veWBn2M/3hlsRfOr'
    b'Yj+iYgdleh7qPwJA+O76P+c/tA4Z9Ev+4j8Yp4P6NWDcP5tfidjLxNI/trzXPshRxT+y4KsEI1m0'
    b'P9Q9MgNyfqA/UvNmR54Iij/yt3oVPCV8PwT9PcqvvHc/CkVZezBCeD/dG7YN1WR5P26cY8ilP3o/'
    b'IN/qDuTKej/UrtROGhR7P/ah3GIU63o/+jCKZ1w4ej9a35Pd1A15Pywc3Qws5Hc/V6rT2A6P+T+I'
    b'8E/N0Y/2P9ueKjokcPI/JlxeATKQ6z/I9H8RHSLiP4ZTP7m0itQ/jHRiBMiZwz+sh23mM7yvPx3f'
    b'Q1rd85g/cEEexpe9ij9GxwE2xoOGP7ahBsXsCYc/b7cSG3IniD/Zujo/BCOJPwfZrvOi5ok/vNFP'
    b'84Juij+jSjpUmmOKP1FJaXjFrYk/Kv86FmKNiD/WwDNWyH+HP8nmHMCeUARAdJ0dqSWtAUAcK8qO'
    b'MfT8P5BvmE9eDPY/+v7bA2cn7T/tWMDA12TgP+YPRh7eMM8/VBHFq2I9uT9aNmtK+sKjP6gPVKVK'
    b'BpU/wNIzBsu8kT8czwUuuzuSP0CCkM61MpM/SC5UjZMVlD8JP84QjMCUPwEmOmjmMJU/hlq0FmQ6'
    b'lT+g2AmmQL2UP7Kx1L+95ZM/2zJf9IUukz9SRon7gqgNQEfapATtDQlAqZVcC7tHBEAUU4/g7lz/'
    b'PyLexlZ20fQ/bB3nK70n5z9VkdzsEObVP8QNnJYsvME/YYqoz5mxqz9AD7ZNcE6dPybnhmIuypg/'
    b'2Gv5Z7KymT+LKRIeSUabP1qpBicNqJw/2jGetjSgnT8IlAgIPhqePwivbuxrFJ4/XBazXNlsnT+f'
    b'1ek73EWcPygYaNpvWJs/LKSfrNKeE0BKOZw9QXkQQJnl6aCEqgpAzo8l+yPCBEDsE5OHjJr7PzJH'
    b'79IEkO4/mgH3wJzi3D8BgspRp3PHP3aRKgdLLbI/Sn/hDG4Aoz/xR3pVRQugP0KEiiYHx6A/5aHm'
    b'MF79oT9fxQ6wdwijP/bh/NOlvaM/AtQ+J3/zoz+oWistMcujPzNayU7lY6M/izTahXivoj/+orur'
    b'bu6hP6RLib4uDRdAMDOn6QlkFEB+6LFjYgIRQDnd4X7SeApAxAf5NQm+AUCWfFlYyMbzP6Q6nNVX'
    b'puI/DqcWqpIrzj8YK8hdlwG3PwBhXtuZk6c/XuYYZeHHoz/MSewtfMWkP1+BJWlkiKY/wDwqezAa'
    b'qD/+T49t0hipP0xZgYgUVak/9gtDKUQSqT9uR8k78qyoP1bY2E1g9qc/x9dbPyjmpj8Y4PJpnBsa'
    b'QPmBx+QH7hhA1q5zT6WSFUAeFBI1X5YQQGpjITFgdAZA/uu6sOdc+T8QtHTpKXvnP/hdap16hdI/'
    b'vr6lxNXWuz8YD2iLc0ysP2npu+5Fpqc/4titT7DiqD+iXefg3yqrP0Eiz1QEMq0/ugPQQNeFrj9g'
    b'a7QisASvP3x7eyow4q4/LCLR29GQrj/+k7NRRvatPyL4vQ3786w/nvlc4IiwG0Bv/nX6Z1kcQCT4'
    b'OC4O7BlAlfyS8GtKFEBmqOtXLOELQLEC1qe0q/8/pYqQRBSD7D9yccUiyaDVP1Rm0xEHEcA/Tqn3'
    b'OHeFsD992VDk18mrP57dz2nqT60/7EoCrwHqrz+RRrEMCw2xPx5VfX3N17E/zmvVSdxLsj9QgToX'
    b'uGOyP/I33nfdS7I/UCVTfxoUsj8ULkl1GcSxPxCogC4xsxxAmjbKIr1cHkA8YTkqNkgdQOWobM2W'
    b'ChhAuiyIGJ7JEEAEsLl4TfkCQFB/3Kjk2vA/8NU9qxzy2D9cmuz1wSPCP23t1g9EsrI/SP/uG9e+'
    b'rz8EnEuIf8+wP4AnvKzwPbI/hd/tnAlesz+Whewxn0O0PwQguLUH/rQ/sHlXCe1QtT8NsxZALzu1'
    b'P7YomyaD6bQ/W19cdryWtD8arxwpwlQfQGZsnVaE1CBAZxLm2+66IEC0xtBweWccQGA1ksFqrBNA'
    b'qjfBGN+6BUDtH2jsBlbzP8r9Pqqno9w/on5LXtdlxD9zHDz6PqS0P7zo01aJlLE/QuuT9ruysj+t'
    b'rVbQA0W0P+TYwcRHd7U/osk6yhuGtj9wZsmA25O3PxYDxHsTRbg/xAG10jF3uD+IDP3iOP63P7LM'
    b'72CLNrc/8bGoqwSmIEBzKmwME+giQJJGBgxjPiNALrpTYnCDIEC8ou1k0FsWQLtX2BIovwdAOMQK'
    b'Tjj79D9xUwBa+WLfPzQiE2L1YsY/8lqN89uItj94iqg4iDGzP654A3HXdLQ/xXpfd3wztj9ldleI'
    b'w323P8Q7ykyTvbg/bQtjVgH5uT9GBtZijQq7P340VzKrDrw/juyd9d/2uz/+rqayrOW6PxoqIMh4'
    b'0iJAiMhHxw45JUAwubZgchMlQG18idYAzSFA7qOP3AEMGEBneUjtfGgJQLjb7zb+IPY/BEw6oOdc'
    b'4D/YiM77w6PHPwXSUPgSarg/LoDjuW7StD8ErHYgYwy2P406PPr0D7g/kMqCjT2duT9DXk/j8vu6'
    b'P2B+Pm0kHrw/KrXLmWKXvT/uN8o6+NC/P2B2tpSwb8A/E6hHWe5TwD+ks881EFQmQIwtoexmDihA'
    b'ekpW8XumJkACH9RGvioiQF5PGPqF1BhAeEPzCOeQC0D0kHNfxiz4P1Qu+z1GXuE/4oztxQzTyD9u'
    b'D+v7UxC6P+X6hqgyc7Y/rjKiF66qtz/nXAPRRt+5PxzyUekLtLs/tvRW+xsavT+XXjpCqDy+P1Ss'
    b'RbmAN8A/pPum0f2jwT+Km8SbbkTCP6x1g356gsI//XZXbSPFKkDkfm2dqrkrQByPRPenPihAx/qc'
    b'op2YIUCox6GHY+8XQMAJzoOCkgxA5i0KkJEj+j/yc03/tNviP/ZpefEJk8o/Kv06pvCauz+yF1ov'
    b'jge4PyaphgqBuLk/wG3KqGMWvD+8vRKDo6i9P3hgJ0x+xL4/lt/8LNcqwD8WqAUWxdzBP7LYiHxs'
    b'JcM/VKMAp83Nwj+68dVWpnDCP5RhejaQqDBApFWXOEWHMECsmG48xj4rQHDDR7P45yFABlmUG9aQ'
    b'F0CpzOYzHl8NQA6RgStXgPs/eMUCr+0B5D8qLgyflzfMP8iIzNv9/7w/OHNobvskuT/yeh/nnYK7'
    b'P2E310u3N74/HPP3O80pvz/evsf1CJi/P4ml5gnb1cA/oxwI+UtFwz+5hNDy43rEP77MvtLKAsM/'
    b'RJpAgj24wT8='
), dtype='double').reshape((20, 20))

value2posterior = RectBivariateSpline(_magnitude_correlation, _ifd_difference, _posterior).ev
