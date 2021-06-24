function frequencySignal = FFT(timeSignal)
% FFT computes the single sided spectrum of the time signal
%   INPUT
%   timeSignal      (1DArray) = input time signal to transform
%   OUTPUT
%   frequencySignal (1DArray) = signal transformed
    L = length(timeSignal);
    frequencySignal = fft(timeSignal);
    frequencySignal = frequencySignal(1:L/2);
    frequencySignal(2:end-1) = 2*frequencySignal(2:end-1);
end


