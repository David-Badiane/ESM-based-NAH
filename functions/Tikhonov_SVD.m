function [x_rec] = Tikhonov_SVD(A , b, alpha)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THIS FUNCTION PERFORMS THE SINGULAR VALUE DECOMPOSITION  %%%%%%%
%%% WITH THE TIKHONOV REGULARIZATION                    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT: A is the matrix of the system Ax = b         %%%%%%%%%%%%
%%%        b is the known term of the system            %%%%%%%%%%%%
%%%        alpha is the regularization parameter        %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OUTPUT: x_rec is the reconstructed solution         %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[U,D,V] = svd(A); % SVD decomposition

svals = diag(D); % singular values of the matrix A

D_alpha = zeros(size(D.')); % regularized D initialization

for j = 1:k % regularized D
    D_alpha(j,j) = svals(j)/(svals(j)^2 + alpha);
end

x_rec = V*D_alpha*U.'*b; 

end

