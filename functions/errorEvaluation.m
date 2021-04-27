function [NCC,NMSE] = errorEvaluation(x,xr)

x = abs(x);
x(isnan(x)) = [];
xr = abs(xr);

% normalize the input vector
x = x./max(x);
xr = xr./max(xr);

% compute the normalized cross correlation
NCC = (xr'*x)/(norm(x)*norm(xr));

% compute the Normalized Mean Square Error 
NMSE = 10*log10(norm(xr - x)^2/(norm(x)^2));
end

