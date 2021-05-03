%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample program for audio mosaicing based on NMF                         %
%                                                                         %
% Coded by D. Kitamura (d-kitamura@ieee.org)                              %
%                                                                         %
% # Original paper                                                        %
% J. Driedger, T. Pratzlich, and M. Muller,                               %
% "Let it bee — towards NMF-inspired audio mosaicing,"                    %
% in Proc. ISMIR, pp. 350–356, 2015.                                      %
%                                                                         %
% # Use pitch shifting algorithm for source signal if you need            %
% J. Driedger and M. Muller,                                              %
% "TSM toolbox: MATLAB implementations of time-scale modification         %
% algorithms," in Proc. DAFx, pp. 249–256, 2014.                          %
%                                                                         %
% See also:                                                               %
% https://www.audiolabs-erlangen.de/resources/MIR/2015-ISMIR-LetItBee     %
% https://www.audiolabs-erlangen.de/resources/MIR/TSMtoolbox              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
addpath('./MATLAB_TSM-Toolbox_2.02');

% Select audio data
dirName = './input/ClToVn'; % directory name of input audio files
% dirName = './input/ObToHr'; % directory name of input audio files
% dirName = './input/VcToTb'; % directory name of input audio files

% Set parameters
seed = 1; % pseudo random seed
fftSize = 4096; % window length in STFT [points]
shiftSize = 2048; % shift length in STFT [points]
windowType = "hamming"; % window function used in STFT
repeatLim = 3; % length of horizontal neighborhood elements to restrict repetition [frames] (denoted as "r" in the paper)
polyphLim = 10; % number of simultaneous activations to restrict polyphony [bases] (denoted as "p" in the paper)
diagKer = 3; % length of diagonal neighborhood elements to enhance continuity (denoted as "c" in the paper)
nIter = 100; % number of iterations
drawConv = true; % true or false (true: show convergence behavior and estimated activation matrix, false: faster and do not plot any figures)

% Fix random seed
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',seed))

% Input data and convert to monaural
[tarSig, sampFreq] = audioread(sprintf('%s/target.mp3', dirName)); % signal x channel x source (source image)
[srcSig, sampFreq] = audioread(sprintf('%s/source.mp3', dirName)); % signal x channel x source (source image)
tarSig = tarSig(:,1); % convert to monaural
srcSig = srcSig(:,1); % convert to monaural

% Audio Mosaicing based on NMF
[estSig, cost] = audioMosaicNmf(tarSig, srcSig, fftSize, shiftSize, windowType, repeatLim, polyphLim, diagKer, nIter, drawConv);

% Output separated signals
outputDir = sprintf('./output');
if ~isdir( outputDir )
    mkdir( outputDir );
end
audiowrite(sprintf('%s/estSig.wav', outputDir), estSig, sampFreq); % estimated signal

fprintf('The files are saved in "./output".\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%