# Audio mosaicing based on nonnegative matrix factorization (NMF)

## About
Sample MATLAB script for audio mosaicing based on NMF coded by Daichi Kitamura

## Reference Paper
J. Driedger, T. Pratzlich, and M. Muller, "Let it bee — towards NMF-inspired audio mosaicing," in Proc. ISMIR, pp. 350–356, 2015.
J. Driedger and M. Muller, "TSM toolbox: MATLAB implementations of time-scale modification algorithms," in Proc. DAFx, pp. 249–256, 2014.

## Contents
- input [dir]:		                includes test audio signals
- input_withPitchShift [dir]:	    includes test audio signals that require pitch shifting before applying mosaicing
- MATLAB_TSM-Toolbox_2.02 [dir]:    inludes TSM toolbox that is utilized for pitch shifting (redistributed under GNU license without any modifications)
- audioMosaicNmf.m:                 NMF-based audio mosaicing
- ISTFT.m:			                inverse short-time Fourier transform
- main.m:			                sample main script
- main_withPitchShift.m:		    sample main script with pitch shifter
- STFT.m:			                short-time Fourier transform

## Usage Note
When the source signal does not include enough pitches for mosaicing, you must poduce various pitches by applying a pitch shifting algorithm. This is implemented in "main_withPitchShift.m" using TSM toolbox, which can be downloaded at https://www.audiolabs-erlangen.de/resources/MIR/TSMtoolbox and redistributed under GNU license.  
For NMF-based audio mosaicing, please visit the authors' website: https://www.audiolabs-erlangen.de/resources/MIR/2015-ISMIR-LetItBee