function frequencySignal = FFT(timeSignal)
% computes the single sided magnitude of the spectrum of the time signal
    L = length(timeSignal);
    frequencySignal = fft(timeSignal);
    frequencySignal = frequencySignal(1:L/2);
    frequencySignal(2:end-1) = 2*frequencySignal(2:end-1);
end


