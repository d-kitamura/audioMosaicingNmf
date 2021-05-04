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
addpath('./MATLAB_TSM-Toolbox_2.02'); % toolbox for pitch shifter

% Select audio data
dirName = './input_withPitchShift/ClToVn'; % directory name of input audio files

% Set parameters
seed = 1; % pseudo random seed
fftSize = 4096; % window length in STFT [points]
shiftSize = 2048; % shift length in STFT [points]
windowType = "hamming"; % window function used in STFT
sigDomain = 1; % signal domain in NMF decomposition (1: amplitude spectrogram, 2: power spectrogram)
nmfCost = "KL"; % similarity measure of NMF (IS: Itakura-Saito divergence, KL: Kullback-Leibler divergence, EU: squared Euclidean distance)
repeatLim = 3; % length of horizontal neighborhood elements to restrict repetition [frames] (denoted as "r" in the paper)
polyphLim = 10; % number of simultaneous activations to restrict polyphony [bases] (denoted as "p" in the paper)
diagKer = 3; % length of diagonal neighborhood elements to enhance continuity (denoted as "c" in the paper)
nIter = 100; % number of iterations
drawConv = true; % true or false (true: show convergence behavior and estimated activation matrix, false: faster and do not plot any figures)
minPitch = 900; % minimum pitch produced in pitch shifter [cent] (relative to source signal, -1200 indicates 1 octave lower)
maxPitch = 1600; % maxmum pitch produced in pitch shifter [cent] (relative to source signal, 1200 indicates 1 octave upper)
intPitch = 50; % pitch interval produced in pitch shifter [cent] (100 corresponds to a semitone in the chromatic scale)
breakLen = 0.1; % break length [s] between each of pitch-shifted tones

% Fix random seed
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',seed))

% Input data and convert to monaural
[tarSig, sampFreq] = audioread(sprintf('%s/target.wav', dirName)); % signal x channel x source (source image)
[srcSig, sampFreq] = audioread(sprintf('%s/source.wav', dirName)); % signal x channel x source (source image)
tarSig = tarSig(:,1); % convert to monaural
srcSig = srcSig(:,1); % convert to monaural

% Produce various pitch sounds of source signal using pitch shifter
param.fsAudio = sampFreq; param.algTSM = @hpTSM; % parameters for "pitchShiftViaTSM.m"
srcLen = size(srcSig, 1); % length of srcSig
tmpSig = zeros(srcLen+breakLen*sampFreq, ((maxPitch-minPitch)/intPitch)); % temporal array for pitch-shifted signals
pitches = minPitch:intPitch:maxPitch; % pitches produced by pitch shifter
for np = 1:size(pitches, 2)
    tmpSig(1:srcLen,np) = pitchShiftViaTSM(srcSig, pitches(np), param);
end
srcSig = tmpSig(:); % concatenate all the pitch-shifted signals with breaks

% Audio Mosaicing based on NMF
[estSig, cost] = audioMosaicNmf(tarSig, srcSig, fftSize, shiftSize, windowType, sigDomain, nmfCost, repeatLim, polyphLim, diagKer, nIter, drawConv);

% Output separated signals
outputDir = sprintf('./output');
if ~isdir( outputDir )
    mkdir( outputDir );
end
audiowrite(sprintf('%s/estSig.wav', outputDir), estSig, sampFreq); % estimated signal

fprintf('The files are saved in "./output".\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%