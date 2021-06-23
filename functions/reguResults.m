function [Qs, v_TSVD, v_TIK] = reguResults( k, lambda, pressure, omega, rho, G_p, G_v, virtualPoints )
% reguResults - Carries on ESM calculations - 
% (regularization inverse problem) + (solution direct problem) 
%   

%   INPUTS
%   k             (double)  = regularization parameter TSVD;
%   lambda        (double)  = regularization parameter TIK;
%   pressure      (1Darray) = hologram vector;
%   omega         (double)  = radians frequency;
%   rho           (double)  = medium density;
%   G_p           (2Darray) = Green matrix - pressure;
%   G_v           (2Darray) = normal gradient Green matrix - velocity;
%   virtualPoints (2Darray) = virtual points matrix (nPts x 3);

%   OUTPUTS
%   Qs            (struct) =  eq sources weight struct - members = qTIK  qTSVD;
%   v_TSVD        (1Darray) = estimated velocity TSVD;
%   v_TIK         (1Darray) = estimated velocity TIK;

    [U,s,V] = csvd (G_p);
    q_TSVD = (1/(1i*omega*rho)).*tsvd (U,s,V, pressure , k); % perform the TSVD -> estimate the source strength
    q_TIK= (1/(1i*omega*rho)).*tikhonov(U,s,V, pressure,lambda);  % perform the Tikhonov SVD -> estimate the source strength
    v_TSVD = - G_v*q_TSVD; % reconstructed velocity with truncated SVD
    v_TIK = - G_v*q_TIK; % reconstructed velocity with Tikhonov
    Qs = struct('qTSVD', q_TSVD, 'qTIK', q_TIK);
    if nargin >= 8
        figure(10)
        virtualPoints(~isnan(virtualPoints(:,3)),3) = abs(q_TIK);
        plot3(virtualPoints(:,1), virtualPoints(:,2), virtualPoints(:,3), '.', 'markerSize', 8 );
    end
end

