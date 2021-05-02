function [points] = L_CurveV(A, b, G_v_omega, v_ex_vector, range, N, rho, omega)

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

alpha_array = logspace(min, max, N);

figure(112) %L curve
hold on
points = zeros(N,3);

for j = 1:N                 %norms are all norm-2
    alpha = alpha_array(j); %current regularization parameter
    [x_rec] = (1/(1i*omega*rho))*Tikhonov_SVD(A , b, alpha);
    seminorm = norm(alpha*x_rec); % seminorm of the regularized solution
    res_norm = norm(G_v_omega*x_rec - v_ex_vector); % norm the residual
    s = scatter(res_norm , seminorm); 
    points(j,1:3) = [res_norm, seminorm, alpha];
end

set(gca,'xscale','log','yscale','log')
xlabel('||A xr - b||_2')
ylabel('||xr||_2')

grid on
hold off
end

