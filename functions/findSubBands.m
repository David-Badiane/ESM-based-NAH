 function [freqIndexes, coarseIndexes] = findSubBands(Hv, fAxis, fAmps, fLocs, deltafLocs, plotData)
%FIND707BAND Summary of this function goes here
%   Detailed explanation goes here
    if plotData
     figure()
     semilogy(fAxis, abs(Hv), 'lineWidth', 1.2);
     xlabel('frequency  [Hz]')
     ylabel('|H_v| / max(|H_v|)')
     hold on;
     stem(fAxis(fLocs),fAmps);
     hold on;
    end
     freqIndexes = zeros(length(fLocs),2);
     coarseIndexes = zeros(length(fLocs),2);
     
     grad = sign(gradient(abs(Hv)));
     bwConst = 0.45;
     for ii = 1:length(fLocs)
         token = true;
         nearestNeighbour = min([deltafLocs(ii),deltafLocs(ii+1)]);
         bw =floor(bwConst * nearestNeighbour );
         counter = bw;
        while token
            leftCheck = false;
            rightCheck = false;
            if mean(grad(fLocs(ii) - counter :fLocs(ii)-1)  >= 0 ) ==1
                bwLeft = fLocs(ii)-counter;
                leftCheck = true;
            end
            if  mean(grad(fLocs(ii)+1:fLocs(ii)+ counter) <= 0) ==1
                bwRight = fLocs(ii)+counter;
                rightCheck = true;
            end
            if(leftCheck && rightCheck)
                token = false;
            end
            counter = counter -1;
        end
        freqIndexes(ii,:) = [bwLeft,bwRight];
        coarseIndexes(ii,:) = [max(fLocs(ii)-bw,1), min(fLocs(ii)+bw,fAxis(end))];
     end
     
     %stem(fAxis(freqIndexes(:)),ones(size(freqIndexes(:))));
     
end

