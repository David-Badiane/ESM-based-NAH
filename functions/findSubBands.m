 function [freqIndexes, coarseIndexes] = findSubBands(Hv, fAxis, fAmps, fLocs, deltafLocs, plotData)
%findSubBands find the bands around the peaks of an FRF, the bands are
%             individuated through the computation of the gradient and 
%             its change of sign
%   INPUTS
%   Hv          (1DArray) = FRF to analyse;
%   fAxis       (1DArray) = frequency axis of the FRF;
%   fAmps       (1DArray) = amplitudes of the peaks;
%   fLocs       (1DArray) = indexes of the peaks;
%   deltafLocs  (1DArray) = (fLocs - circshift(fLocs,1)) to coarsely get
%                           the subbands boundaries
%   plotData    (boolean) = see images
%   OUTPUTS
%   freqIndexes   (1DArray) = frequency indexes of the bands after gradient
%       computation - the limit bwLeft is independent of bwRight - subbands 
%       may be asymmetric;
%   coarseIndexes (1DArray) = coarse indexes obtained by evaluating deltafLocs -
%                             subbands are symmetric;


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
end

