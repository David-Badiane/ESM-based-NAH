function [Hv,f0, fLocs] = EMASimple(HvSVD, fAxis,minPeakVal, minPeakWidth, plotData)
%EMASimple
%
%   Simplified modal analysis algorithm retrieving eigenfrequencies, modal
%   damping ratios and modeshape value in the point
%
%   HvSVD   (array)       = spectrum to analyse;
%   fAxis   (array)       = frequency axis of the spectrum;
%   minPeakVal (double)   = minimum value of the peaks for peak analysis;
%   minPeakWidth (double) = minimum value of the width of the maximum;
    
    Hv = HvSVD(1:length(fAxis));
    fHigh = fAxis(end);
    [fAmps, fLocs] = findpeaks(abs(Hv),'MinPeakProminence',minPeakVal,'MinPeakWidth', minPeakWidth);
    fLocs = fLocs(:);
    Hv = Hv(:);
    f0 = fAxis(fLocs);
    
    if plotData
         figure()
         semilogy(fAxis, abs(Hv));
         hold on;
         for ii = 1:length(fLocs)
             stem(fAxis(fLocs(ii)),fAmps(ii));
         end
    end  
end