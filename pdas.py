import os
import uuid
import pathlib
import numpy
import resampy
import soundfile
from transplant import Matlab
from maps_f0 import SignalBlocks, estimate_f0

# restrict Numpy to a single thread:
os.putenv('OMP_NUM_THREADS', '1')

def maps(signal, samplerate):
    blocks = SignalBlocks(signal, samplerate)
    return estimate_f0(blocks)

def yin(signal, samplerate):
    if samplerate != 16000:
        signal = resampy.resample(signal, samplerate, 16000)
        samplerate = 16000

    tmpfilepath = pathlib.Path.cwd() / (str(uuid.uuid1()) + '.wav')
    try:
        soundfile.write(str(tmpfilepath), signal, samplerate)

        with Matlab(arguments=('-nodesktop', '-nosplash', '-nojvm', '-singleCompThread'), print_to_stdout=False) as matlab:
            algopath = pathlib.Path(__file__).parents[0] / 'YIN'
            matlab.cd(str(algopath))
            struct = matlab.yin(str(tmpfilepath), nargout=1)
    finally:
        if tmpfilepath.exists():
            tmpfilepath.unlink()

    f0 = 440 * 2**struct['f0'].ravel()
    hopsize = struct['hop']
    time = numpy.arange(len(f0))*hopsize/samplerate
    p = numpy.ones(len(f0))

    return time, f0, p


def pefac(signal, samplerate):
    with Matlab(arguments=('-nodesktop', '-nosplash', '-nojvm', '-singleCompThread'), print_to_stdout=False) as matlab:
        algopath = pathlib.Path(__file__).parents[0] / 'PEFAC'
        matlab.cd(str(algopath))
        f0, time, prob = matlab.fxpefac(signal, float(samplerate), nargout=3)

    return time.ravel(), f0.ravel(), prob.ravel()


def rapt(signal, samplerate):
    with Matlab(arguments=('-nodesktop', '-nosplash', '-nojvm', '-singleCompThread'), print_to_stdout=False) as matlab:
        algopath = pathlib.Path(__file__).parents[0] / 'RAPT'
        matlab.cd(str(algopath))
        f0, times = matlab.fxrapt(signal, float(samplerate), 'u', nargout=2)
        f0 = f0.ravel()
        time = numpy.empty(len(f0)*2)
        time[::2] = (times[:,0]-1)/samplerate
        time[1::2] = times[:,1]/samplerate
        f0s = numpy.empty(len(f0)*2)
        f0s[::2] = f0s[1::2] = f0
        p = numpy.ones(len(f0s))
        p[numpy.isnan(f0s)] = 0

    return time, f0s, p
