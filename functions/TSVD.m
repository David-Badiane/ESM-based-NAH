function [x_rec, svals] = TSVD(A,b,k)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THIS FUNCTION PERFORMS THE TRUNCATED SINGULAR VALUE %%%%%%%%%%%%
%%% DECOMPOSITION AT THE K-TH SINGULAR VALUE            %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT: A is the matrix of the system Ax = b         %%%%%%%%%%%%
%%%        b is the known term of the system            %%%%%%%%%%%%
%%%        k is the order of the truncation             %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OUTPUT: x_rec is the reconstructed solution         %%%%%%%%%%%%
%%%         svals is a vector containing the singular values %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[U,D,V] = svd(A); % SVD decomposition

svals = diag(D); % singular values of the matrix A

D_k = zeros(size(D.')); % truncated D initialization

for j = 1:k % truncate D (TSVD)
    D_k(j,j) = 1/svals(j);
end

x_rec = V*D_k*U.'*b; 

end

