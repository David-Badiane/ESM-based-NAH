function [outputArg1,outputArg2] = plotErrorPressure(measuredPressure, G_p_omega , measuredPressureN , omega , rho , range , numberParameters  )

min = range(1); % the lowest regularization parameter
max = range(2); % the highest regularization parameter

alpha_array = linspace(min, max, numberParameters);
reconstructionError = zeros(1, length(alpha_array));


for j = 1:numberParameters                 %norms are all norm-2
    alpha = alpha_array(j); %current regularization parameter
    q_TIK = (1/(1i*omega*rho)).*Tikhonov_SVD(G_p_omega , measuredPressureN  , alpha); % reconstructed source streght
    p_TIK = 1i*omega*rho*G_p_omega*q_TIK; % reconstructed velocity with Tikhonov
    reconstructionError(j) = (p_TIK'*measuredPressure)/(norm(measuredPressure)*norm(p_TIK)) ;
end

figure(112) %Reconstruction error for pressure
plot(alpha_array, reconstructionError)
xlabel('alpha')
ylabel('NCC')
grid on
%ylim([-1 1])

end

