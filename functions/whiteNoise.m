function [out] = whiteNoise(in, SNR)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% WHITENOISE this function adds white gaussian noise to the input %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   INPUTS
%   in   (1Darray)   = input array;
%   SNR  (double)  = signal to noise ratio;

%   OUTPUT
%   out  (1Darray)   = output array;

snr = SNR; % [dB] Signal to Noise ratio
out = awgn( in , snr, 'measured');

end

