function [Hv,f0, fLocs, csis, Q, modeShapes] = EMASimple(HvSVD, fAxis,minPeakVal, minPeakWidth, plotData)
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
    [freqIndexes, coarseIndexes] = findSubBands(Hv, fAxis, fAmps, fLocs, deltafLocs, plotData); 
    bandsLength = freqIndexes(:,2)-freqIndexes(:,1);
    
    
    csis   = 100*ones(length(fLocs),1);
    modeShapes = zeros(length(fLocs),1);
    % Fit parabola into bands to retrieve csi with half power point method
    
    if plotData
     figure(33)
    end
    for ii = 1:length(freqIndexes(:,1))
        % get sub band and its frequency axis
        subBand = Hv(freqIndexes(ii,1):freqIndexes(ii,2));
        subBand = abs(subBand)/max(abs(subBand));
        subBandAxis = fAxis(freqIndexes(ii,1):freqIndexes(ii,2)); 
        
        % fit parabola
        ind = find(subBand>=0.5);
        subBand = subBand(ind);
        subBandAxis = subBandAxis(ind);
        p = polyfit(subBandAxis,subBand,2);
        fAxisFit = linspace(subBandAxis(1)-3,subBandAxis(end)+3, 2000);
        y1 = polyval(p,fAxisFit);
        
        % cut parabola
        ind = find(y1>=0.5);
        fAxisFit = fAxisFit(ind);
        y1 = y1(ind);
        
        % search 0.707 point with a proper loss function
        lossFunction = 1./(abs(y1-0.707*ones(size(y1))));
        [minvals,minlocs] = findpeaks(lossFunction, 'MinPeakProminence', 10);
        
        % calculate modal parameters of subBands with intercept at 0.707
        if length(minlocs) ==2 
            w1 = 2*pi*fAxisFit(minlocs(1));
            w2 = 2*pi*fAxisFit(minlocs(2));
            csis(ii) =(w2^2 - w1^2)/(4*(w_nat(ii))^2);
            c=2*w_nat(ii)*csis(ii);
            modeShapes(ii) = -imag(Hv(fLocs(ii))*c*w_nat(ii));
        end
        
       % plot images
       if plotData 
        subplot(5,6,ii)
          
        plot(subBandAxis, subBand,'LineWidth',1.8)
        hold on 
        plot(fAxisFit ,y1, fAxisFit, 0.707*ones(size(fAxisFit)));
        
        ylim([0,1]);
        legend('band', 'polyfit', '0.707 line')
        xlabel('f  [Hz]')
        ylim([0.2,1]);
        yticks([0,0.707, 1]);
       end
    end
    Q = 1./(2*csis);
    
end

