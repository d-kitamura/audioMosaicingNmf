function [est, cost] = audioMosaicNmf(tarSig, srcSig, fftLen, shiftLen, winType, dom, criterion, repLim, polLim, diaKer, nIter, isDraw)
% Audio Mosaicing based on NMF
%
% Coded by D. Kitamura (d-kitamura@ieee.org)
%
% # Original paper
% J. Driedger, T. Pratzlich, and M. Muller,
% "Let it bee — towards NMF-inspired audio mosaicing,"
% in Proc. ISMIR, pp. 350–356, 2015.
%
% See also:
% https://www.audiolabs-erlangen.de/resources/MIR/2015-ISMIR-LetItBee
% https://www.audiolabs-erlangen.de/resources/MIR/TSMtoolbox
%
% [syntax]
%   [est,cost] = audioMosaicNmf(tar, src, fftLen, shiftLen, winType, dom, criterion, repLim, polLim, diaKer, nIter, isDraw)
%
% [inputs]
%      tarSig: target signal (sigLen x 1)
%      srcSig: source signal (sigLen x 1)
%      fftLen: window length [points] in STFT (scalar, default: next higher power of 2 that exceeds 0.256*sampFreq)
%    shiftLen: shift length [points] in STFT (scalar, default: fftSize/2)
%     winType: window function used in STFT (name of window function, default: 'hamming')
%         dom: signal domain in NMF decomposition (scalar, default: 1 (amplitude spectrogram))
%   criterion: similarity measure of NMF (name of similarity function, default: KL (Kullback-Leibler divergence))
%      repLim: length of horizontal neighborhood elements to restrict repetition (scalar, default: 3, denoted as "r" in the paper)
%      polLim: number of simultaneous activations to restrict polyphony (scalar, default:10, denoted as "p" in the paper)
%      diaKer: length of diagonal neighborhood elements to enhance continuity (scalar, default: 3, denoted as "c" in the paper)
%       nIter: number of iterations in the parameter update in NMF (scalar, default: 100)
%      isDraw: plot cost function values in each iteration or not (true or false, default: false)
%
% [outputs]
%         est: estimated signals (sigLen x 1)
%        cost: convergence behavior of cost function in NMF (nIter+1 x 1)
%

% Arguments check and set default values
arguments
    tarSig (:,1) double
    srcSig (:,1) double
    fftLen (1,1) double {mustBeInteger(fftLen)} = 2^nextpow2(0.256*sampFreq)
    shiftLen (1,1) double {mustBeInteger(shiftLen)} = fftLen/2
    winType char {mustBeMember(winType,{'hamming','hann','rectangular','blackman','sine'})} = 'hamming'
    dom (1,1) double = 1
    criterion char {mustBeMember(criterion,{'IS','KL','EU'})} = 'KL'
    repLim (1,1) double {mustBeInteger(repLim)} = 3
    polLim (1,1) double {mustBeInteger(polLim)} = 10
    diaKer (1,1) double {mustBeInteger(diaKer)} = 3
    nIter (1,1) double {mustBeInteger(nIter)} = 100
    isDraw (1,1) logical = false
end

% Error check
[tarLen, tarCh] = size(tarSig); % tarLen: length of target signal, tarCh: number of channels
[srcLen, srcCh] = size(srcSig); % srcLen: length of source signal, srcCh: number of channels
if tarCh ~= 1 || srcCh ~= 1; error("The input audio signals must be monaural.\n"); end
if fftLen < 1; error("The FFT length in STFT (fftSize) must be a positive integer value.\n"); end
if shiftLen < 1; error("The shift length in STFT (shiftSize) must be a positive integer value.\n"); end
if repLim < 1; error("The parameter (repLim) must be a positive integer value.\n"); end
if polLim < 1; error("The parameter (polLim) must be a positive integer value.\n"); end
if diaKer < 1; error("The parameter (diaKer) must be a positive integer value.\n"); end
if nIter < 1; error("The number of iterations (nIter) must be a positive integer value.\n"); end

% Apply short-time Fourier transform (STFT)
[tarSpecgram, winInStft] = STFT(tarSig, fftLen, shiftLen, winType);
[srcSpecgram, winInStft] = STFT(srcSig, fftLen, shiftLen, winType);

% Apply NMF
[estSpecgram, cost] = local_NMF(tarSpecgram, srcSpecgram, dom, criterion, repLim, polLim, diaKer, nIter, isDraw);

% Apply inverse STFT
est = ISTFT(estSpecgram, shiftLen, winInStft, tarLen);
end

%% Local function for NMF
function [cY, cost] = local_NMF(cV, cW, dom, sim, r, p, c, nIter, isDraw)
% [inputs]
%       cV: target complex-valued spectrogram (I x J1)
%       cW: source complex-valued spectrogram (I x J2)
%      dom: signal domain in NMF decomposition (scalar)
%      sim: similarity measure in NMF decomposition
%        r: length of horizontal neighborhood elements to restrict repetition
%        p: number of simultaneous activations to restrict polyphony
%        c: length of diagonal neighborhood elements to enhance continuity
%    nIter: number of iterations of the parameter updates
%   isDraw: plot cost function values in each iteration and estimated activation matrix or not (true or false)
%
% [outputs]
%       cY: mosaiced complex-valued spectrogram (I x J1)
%     cost: convergence behavior of cost function in NMF (nIter+1 x 1)
%

% Initialization
[I,J1] = size(cV); % I:frequency bins, J1: time frames
[I,J2] = size(cW); % I:frequency bins, J2: time frames
V = max(abs(cV).^dom, eps); % amplitude spectrogram of V
W = max(abs(cW).^dom, eps); % amplitude spectrogram of W
H = max(rand(J2,J1), eps); % activation matrix
E = ones(I,J1); % ones matrix
ker = eye(2*c+1); % kernel matrix for enhancing continuity by 2D convolution
cost = zeros(nIter+1, 1);

% Calculate initial cost function value
if isDraw; cost(1,1) = local_calcCost(V, W*H, sim); end

% Optimize activation matrix H
fprintf('Iteration:    ');
for iIter = 1:nIter
    fprintf('\b\b\b\b%4d', iIter);
    
    %%%%% Update parameters %%%%%
    mu = movmax(H, 2*r+1, 2); % moving max using (2r+1)-length filter for dim=2 (Eq. 4)
    idxR = (H==mu); % indices that satisfy H = mu (true or false)
    R = H.*idxR + H.*(~idxR)*(1-iIter/nIter); % repetition restricted activation matrix (Eq. 3)
    topVal = maxk(R, p, 1); % find p-highest values for dim=1
    idxP = (R >= min(topVal,[],1)); % indices that satisfy R >= minumum of top-p value (true or false)
    P = R.*idxP + R.*(~idxP)*(1-iIter/nIter); % polyphony restricted activation matrix (Eq. 5)
    C = conv2(P, ker); % 2-D convolution with diagonal (unit) matrix kernel
    C = C(1+c:end-c, 1+c:end-c); % continuity enhancing activation matrix (Eq. 6)
    switch sim
        case "IS"
            H = C .* sqrt( (W.'*(V./(W*C)./(W*C))) ./ (W.'*(E./(W*C))) ); % update activation matrix (IS div.)
        case "KL"
            H = C .* ( (W.'*(V./(W*C))) ./ (W.'*E) ); % update activation matrix (KL div., Eq. 7)
        case "EU"
            H = C .* ( (W.'*V) ./ (W.'*(W*C)) ); % update activation matrix (EU dist.)
    end
    
    %%%%% Calculate cost function value %%%%%
    if isDraw; cost(iIter+1,1) = local_calcCost(V, W*H, sim); end
end

% Reconstruct mosaiced complex-valued spectrogram
cY = cW*H;

% Draw convergence behavior and estimated activation matrix
if isDraw
    figure; semilogy((0:nIter), cost);
    set(gca, 'FontName', 'Times', 'FontSize', 16);
    xlabel('Number of iterations', 'FontName', 'Arial', 'FontSize', 16);
    ylabel('Value of cost function', 'FontName', 'Arial', 'FontSize', 16);
    
    figure; imagesc(H); caxis([0, 1]);
    set(gca, 'FontName', 'Times', 'FontSize', 16);
    title('Estimated activation matrix');
    xlabel('Target time frame', 'FontName', 'Arial', 'FontSize', 16);
    ylabel('Source time frame', 'FontName', 'Arial', 'FontSize', 16);
end
fprintf(' NMF done.\n');
end

%% Local function for calculating cost function value in NMF
function [ cost ] = local_calcCost(V, WH, sim)
switch sim
    case "IS"
        cost = sum((V./WH)-log(V./WH)-1, 'all'); % Itakura-Saito divergence
    case "KL"
        cost = sum(V.*log(V./WH)-V+WH, 'all'); % Generalized Kullback-Leibler divergence
    case "EU"
        cost = sum((V - WH).^2, 'all'); % Squared Euclidean distance
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%