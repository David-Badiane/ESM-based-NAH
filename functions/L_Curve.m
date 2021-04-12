function [outputArg1,outputArg2] = L_Curve(A, b, range, N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THIS FUNCTION PLOTS THE L-CURVE. USEFUL TO CHOOSE THE    %%%%%%%
%%% REGULARIZATION PARAMETER.                                %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT: A is the matrix of the system Ax = b         %%%%%%%%%%%%
%%%        b is the known term of the system            %%%%%%%%%%%%
%%%        range is the interval of the regularization  %%%%%%%%%%%%
%%%        parameter                                    %%%%%%%%%%%%
%%%        N is the number of parameters                %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OUTPUT: plot of the L curve                         %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

min = range(1); % the lowest regularization parameter
max = range(2); % the highest regularization parameter

alpha_array = linspace(min, max, N);

figure %L curve
hold on
for j = 1:N                 %norms are all norm-2
    alpha = alpha_array(j); %current regularization parameter
    [x_rec] = Tikhonov_SVD(A , b, alpha) 
    seminorm = norm(x_rec); % seminorm of the regularized solution
    res_norm = norm(A*x_rec - b); % norm the residual
    loglog(res_norm , seminorm) 
end

end
