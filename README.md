# Audio mosaicing based on nonnegative matrix factorization (NMF)

## About
Sample MATLAB script for audio mosaicing based on NMF coded by Daichi Kitamura

## Reference Paper
J. Driedger, T. Pratzlich, and M. Muller, "Let it bee — towards NMF-inspired audio mosaicing," in Proc. ISMIR, pp. 350–356, 2015.
J. Driedger and M. Muller, "TSM toolbox: MATLAB implementations of time-scale modification algorithms," in Proc. DAFx, pp. 249–256, 2014.

## Contents
- input [dir]:		        includes test audio signals
- audioMosaicNmf.m:         NMF-based audio mosaicing
- ISTFT.m:			        inverse short-time Fourier transform
- main.m:			        main script with parameter settings
- STFT.m:			        short-time Fourier transform

## Usage Note
When the source signal does not include enough pitches for mosaicing, you must poduce various pitches by applying a pitch shifting algorithm. 
According to the reference paper, "pitchShiftViaTSM.m" can be used, which is a part of Toolbox called "TSM toolbox" and can be downloaded at https://www.audiolabs-erlangen.de/resources/MIR/TSMtoolbox
Useful information can be found at https://www.audiolabs-erlangen.de/resources/MIR/2015-ISMIR-LetItBee