function [Hv,f0, fLocs, csis, Q] = EMAPoly(HvSVD, fAxis,minPeakVal, minPeakWidth)
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
    w_nat = 2*pi*f0; 
    deltafLocs = fLocs - [0;fLocs(1:end-1)];
    deltafLocs = [deltafLocs;length(fAxis)-fLocs(end)];

%      figure()
%      semilogy(fAxis, abs(Hv));
%      hold on;
%      for ii = 1:length(fLocs)
%          stem(fAxis(fLocs(ii)),fAmps(ii));
%      end
    [freqIndexes, coarseIndexes] = findSubBands(Hv, fAxis, fAmps, fLocs, deltafLocs); 
    bandsLength = freqIndexes(:,2)-freqIndexes(:,1);
    
    csis   = 100*ones(length(fLocs),1);
    shapes = zeros(length(fLocs),1);
    figure()
    for ii = 1:length(freqIndexes(:,1))
        subBand = Hv(freqIndexes(ii,1):freqIndexes(ii,2));
        subBand = abs(subBand)/max(abs(subBand));
        subBandFreq = fAxis(freqIndexes(ii,1):freqIndexes(ii,2)); 
        
        ind = find(subBand>=0.65);
        subBand = subBand(ind);
        subBandFreq = subBandFreq(ind);
        p = polyfit(subBandFreq,subBand,2);
        fAxisFit = linspace(subBandFreq(1)-5,subBandFreq(end)+5, 2000);
        y1 = polyval(p,fAxisFit);
        
        ind = find(y1>=0.6);
        fAxisFit = fAxisFit(ind);
        y1 = y1(ind);
        
        lossFunction = 1./(abs(y1-0.707*ones(size(y1))));
        [minvals,minlocs] = findpeaks(lossFunction, 'MinPeakProminence', 10);
        
        if length(minlocs) ==2 
            w1 = 2*pi*fAxisFit(minlocs(1));
            w2 = 2*pi*fAxisFit(minlocs(2));
            csis(ii) =(w2^2 - w1^2)/(4*(w_nat(ii))^2);
            c=2*w_nat(ii)*csis(ii);
            shapes(ii) = -imag(Hv(fLocs(ii))*c*w_nat(ii));
        end
        
        subplot(4,3,ii)
          
        plot(subBandFreq, subBand,'LineWidth',1.8)
        hold on 
        plot(fAxisFit ,y1, fAxisFit, 0.707*ones(size(fAxisFit)));
        
        ylim([0,1]);
        legend('band', 'polyfit', '0.707 line')
        xlabel('f  [Hz]')
        ylim([0.2,1]);
        yticks([0,0.707, 1]);
    end
    Q = 1./(2*csis);
end

