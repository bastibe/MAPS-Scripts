This is part of the dissertation [Pitch of Voiced Speech in the Short-Time Fourier Transform: Algorithms, Ground Truths, and Evaluation Methods](https://bastibe.github.io/Dissertation-Website/), on the topic of [Fundamental Frequency Estimation](https://bastibe.github.io/Dissertation-Website/maps/index.html)  
(Preprint Manuscript)  
© 2020, Bastian Bechtold, Jade Hochschule & Carl von Ossietzky Universität Oldenburg, Germany.

# MAPS: A Fundamental Frequency Estimation Algorithm for Speech in the Magnitude and Phase Spectrogram

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3923585.svg)](https://doi.org/10.5281/zenodo.3923585) ![GitHub](https://img.shields.io/github/license/bastibe/MAPS-Scripts)

This repository contains source code for MAPS, the *Magnitude and Phase Spectrogram* fundamental frequency estimation algorithm.

Implementations of the algorithm are provided in Python (*maps_f0.py*), Matlab (*maps_f0.m*), and Julia (*maps_f0.jl*). MAPS is provided under the terms of the GPLv3 license. [PEFAC](http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html)&nbsp;[1], [RAPT](http://www.speech.kth.se/wavesurfer/links.html)&nbsp;[2], and [YIN](http://audition.ens.fr/adc/)&nbsp;[3] are covered by their respective licenses.

Additionally, *MAPS Evaluation.ipynb* contains a reproducible research notebook for comparing MAPS to the well-known fundamental frequency estimation algorithms [PEFAC](http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html)&nbsp;[1], [RAPT](http://www.speech.kth.se/wavesurfer/links.html)&nbsp;[2], and [YIN](http://audition.ens.fr/adc/)&nbsp;[3] on the [PTDB-TUG](https://www.spsc.tugraz.at/databases-and-tools/ptdb-tug-pitch-tracking-database-from-graz-university-of-technology.html)&nbsp;[4] speech corpus and the [QUT-NOISE](https://research.qut.edu.au/saivt/databases/qut-noise-databases-and-protocols/)&nbsp;[5] noise corpus.

## References:

1. Sira Gonzalez and Mike Brookes. PEFAC - A Pitch Estimation Algorithm Robust to High Levels of Noise. IEEE/ACM Transactions on Audio, Speech, and Language Processing, 22(2):518—530, February 2014.
2. David Talkin. A robust algorithm for pitch tracking (RAPT). Speech coding and synthesis, 495:518, 1995.
3. Alain de Cheveigné and Hideki Kawahara. YIN, a fundamental frequency estimator for speech and music. The Journal of the Acoustical Society of America, 111(4):1917, 2002.
4. Gregor Pirker, Michael Wohlmayr, Stefan Petrik, and Franz Pernkopf. A Pitch Tracking Corpus with Evaluation on Multipitch Tracking Scenario. page 4, 2011.
5. David B. Dean, Sridha Sridharan, Robert J. Vogt, and Michael W. Mason. The QUT-NOISE-TIMIT corpus for the evaluation of voice activity detection algorithms. Proceedings of Interspeech 2010, 2010.
