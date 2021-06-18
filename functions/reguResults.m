function [Qs, v_TSVD, v_TIK] = reguResults( k, lambda, pressure, omega, rho, G_p, G_v, virtualPoints )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    [U,s,V] = csvd (G_p);
    q_TSVD = (1/(1i*omega*rho)).*tsvd (U,s,V, pressure , k); % perform the TSVD -> estimate the source strength
    q_TIK= (1/(1i*omega*rho)).*tikhonov(U,s,V, pressure,lambda);  % perform the Tikhonov SVD -> estimate the source strength
    v_TSVD = G_v*q_TSVD; % reconstructed velocity with truncated SVD
    v_TIK = G_v*q_TIK; % reconstructed velocity with Tikhonov
    Qs = struct('qTSVD', q_TSVD, 'qTIK', q_TIK);
    figure(10)
    virtualPoints(~isnan(virtualPoints(:,3)),3) = abs(q_TIK);
    plot3(virtualPoints(:,1), virtualPoints(:,2), virtualPoints(:,3), '.', 'markerSize', 8 );
end

