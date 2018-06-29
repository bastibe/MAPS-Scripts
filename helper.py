import numpy
from scipy.interpolate import interp1d


def create_test_signal(speech, noise, noise_start, snr):
    samplerate = speech.signal.metadata['samplerate']
    speech_signal = speech.signal
    noise_start = int(noise_start*samplerate)
    if isinstance(noise, str) and noise == 'white noise':
        noise_signal = numpy.random.randn(len(speech_signal))
    else:
        noise_signal = noise.signal[noise_start:noise_start+len(speech_signal)]
    if len(noise_signal.shape) > 1:
        noise_signal = noise_signal[:, 0]
    return mix(speech_signal, noise_signal, snr, samplerate)


def mix(signal, noise, snr_db, samplerate):
    signal_power = vad_power(signal, samplerate)
    noise_power = rms(noise)
    mix_signal = signal + noise/noise_power*signal_power / 10**(snr_db/20)
    return mix_signal/numpy.max(mix_signal)


def rms(signal):
    return numpy.mean(signal**2)**0.5


def vad_power(signal, samplerate, blocklength=0.02):
    """The power of a speech signal, only where there is speech.

    Calculate the signal level for short signal blocks, remove
    outliers, then define a threshold at halfway between the min and
    max level. The speech power is calculated from all signal values
    above that threshold.

    """
    # calculate signal level for short signal blocks:
    blocksize = int(samplerate*blocklength)
    levels = [20*numpy.log10(rms(signal[idx:idx+blocksize]))
              for idx in range(0, len(signal)-blocksize, blocksize)]
    # remove outliers:
    low, high = numpy.percentile(levels, [5, 95])
    no_outliers = numpy.array(levels)[(low < levels) & (levels < high)]
    # calculate threshold:
    threshold = (max(no_outliers) + min(no_outliers))/2
    # calculate signal power only where the signal is active:
    powers = [(10**(lvl/20))**2 for lvl in levels
              if lvl > threshold]
    return numpy.mean(powers)**0.5


def generate_htc(duration, basetimes, basefreqs, samplerate, harmonics=10, random_phase=True):
    htc = numpy.zeros(duration*samplerate)
    basefreq = interp1d(basetimes, basefreqs, kind='quadratic', copy=False, assume_sorted=True)
    time = numpy.linspace(0, duration, len(htc))
    for harmonic in range(harmonics):
        freq = basefreq(time)*(harmonic+1)
        if numpy.max(freq) > samplerate/2: break
        phase = numpy.cumsum(2*numpy.pi * freq / samplerate)
        start_phase = 0 if not random_phase else 2*numpy.pi*numpy.random.random()
        htc += numpy.cos(start_phase + phase)
    return htc


def _same_t_in_true_and_est(func):
    def new_func(true_t, true_f, est_t, est_f):
        interpolated_f = interp1d(est_t, est_f, bounds_error=False, kind='nearest', fill_value=0)(true_t)
        return func(true_t, true_f, true_t, interpolated_f)
    return new_func


@_same_t_in_true_and_est
def gross_pitch_error(true_t, true_f, est_t, est_f):
    """The relative frequency in percent of pitch estimates that are
    outside a threshold around the true pitch. Only frames that are
    considered pitched by both the ground truth and the estimator (if
    applicable) are considered.

    """

    correct_frames = _true_voiced_frames(true_t, true_f, est_t, est_f)
    gross_pitch_error_frames = _gross_pitch_error_frames(true_t, true_f, est_t, est_f)
    return numpy.sum(gross_pitch_error_frames)/numpy.sum(correct_frames)


def _gross_pitch_error_frames(true_t, true_f, est_t, est_f):
    voiced_frames = _true_voiced_frames(true_t, true_f, est_t, est_f)
    pitch_error_frames = numpy.abs(est_f/true_f - 1) > 0.2
    return voiced_frames & pitch_error_frames


def _true_voiced_frames(true_t, true_f, est_t, est_f):
    return (est_f != 0) & (true_f != 0)


@_same_t_in_true_and_est
def fine_pitch_error(true_t, true_f, est_t, est_f):
    """The mean error of pitch estimates that are correct according to the
    GPE."""
    voiced_frames = _true_voiced_frames(true_t, true_f, est_t, est_f)
    correct_frames = numpy.abs(est_f/true_f - 1) < 0.2
    correct_voiced_frames = voiced_frames & correct_frames
    return numpy.mean(numpy.abs(est_f[correct_voiced_frames]/true_f[correct_voiced_frames] - 1))


@_same_t_in_true_and_est
def true_positive_rate(true_t, true_f, est_t, est_f):
    true_positives  = numpy.sum((true_f != 0) & (est_f != 0))
    false_negatives = numpy.sum((true_f != 0) & (est_f == 0))
    return true_positives / (true_positives + false_negatives)


@_same_t_in_true_and_est
def false_positive_rate(true_t, true_f, est_t, est_f):
    false_positives = numpy.sum((true_f == 0) & (est_f != 0))
    true_negatives  = numpy.sum((true_f == 0) & (est_f == 0))
    return false_positives / (false_positives + true_negatives)


@_same_t_in_true_and_est
def false_negative_rate(true_t, true_f, est_t, est_f):
    false_negatives = numpy.sum((true_f != 0) & (est_f == 0))
    true_positives  = numpy.sum((true_f != 0) & (est_f != 0))
    return false_negatives / (false_negatives + true_positives)
