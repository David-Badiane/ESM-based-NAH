function [out] = whiteNoise(in, SNR)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% WHITENOISE this function adds white gaussian noise to the input %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   INPUTS
%   in   (array)   = input array;
%   SNR  (double)  = signal to noise ratio;

%   OUTPUT
%   out  (array)   = output array;

snr = SNR; % [dB] Signal to Noise ratio

out = awgn( in , snr, 'measured');

end

