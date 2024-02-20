Pcyl_mmH20= [165, 164, 164, 161, 152, 142, 133, 120, 108, 97, 85, 72, 64, 56, 55, 55, 55, 58, 59, 58, 55, 55, 55, 56, 64, 72, 85, 97, 108, 120, 133, 142,152, 161, 164, 164, 165];
Pinf_mmH20= 115 ;

Pcyl_pa= Pcyl_mmH20*9.80665 ;
Pinf_pa= Pinf_mmH20*9.80665;
Pstag= Pcyl_pa(1);

theta= [0:10:360].';
cos_theta= cosd(theta);

Cp_m= ((Pcyl_pa-Pinf_pa)/(Pstag-Pinf_pa)).';
Cp_m_cos_theta= (Cp_m.*cos_theta);

Cd_m= 0.5*trapz(theta, Cp_m_cos_theta)


x= [0:360];
Cp_t= 1-4*sind(x).^2;
Cp_t_cos_theta = (1-4*sind(x).^2).*cosd(x);
Cp_t_cos_theta_func= @(x) (1-4*sind(x).^2).*cosd(x);
Cd_t=0.5*integral(Cp_t_cos_theta_func,0,360)



hold on
plot(theta,Cp_m_cos_theta)
plot(x,Cp_t_cos_theta)
xlabel('θ(°)','FontSize', 10)
ylabel('Cpcos(θ)','FontSize', 10)
title('Cpcos(θ) vs θ', 'FontSize', 13)  
legend('Measured', 'Theoritical')
