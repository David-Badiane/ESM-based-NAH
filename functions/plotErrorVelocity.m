function [outputArg1,outputArg2] = plotErrorVelocity(v_ex_vector, measuredPressureN, G_p_omega, G_v_omega, range, numberParameters, omega, rho)

min = range(1); % the lowest regularization parameter
max = range(2); % the highest regularization parameter

alpha_array = linspace(min, max, numberParameters);
NCC_array = zeros(1, length(alpha_array));
NMSE_array = zeros(1, length(alpha_array));

for j = 1:numberParameters                 %norms are all norm-2
    alpha = alpha_array(j); %current regularization parameter
    q_TIK = (1/(1i*omega*rho)).*Tikhonov_SVD(G_p_omega , measuredPressureN  , alpha); % reconstructed source streght
    v_TIK = G_v_omega*q_TIK; % reconstructed velocity with Tikhonov
%     v_ex_vector(v_ex_vector==0) = [];
%     v_TIK(v_TIK==0)=[];
%     v_TIK(1) = [];
    [NCC , NMSE] = errorEvaluation(v_ex_vector, v_TIK);
    NCC_array(j) = NCC;
    NMSE_array(j) = NMSE;
end

figure(111) %NCC 
plot(alpha_array, NCC_array)
hold on
plot(alpha_array, NMSE_array)
hold off
legend('NCC' , 'NMSE')
xlabel('alpha')
ylabel('NCC')

end

