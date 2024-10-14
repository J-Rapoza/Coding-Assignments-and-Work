clc
clear all

%Known Temperatures/pressures for heat exchangers (temperature in celsius,
%pressure in bar)
T4=679.6943-273.15;
T5=500-273.15;
T6=40.0228;
T7=400;
P67=10;
P45=1;
m_dot_vapor=9.69;
m_dot_air=30E3/(T4-T5);


%%%%Calculations for overall Q
%superheater
Qs=m_dot_vapor*(XSteam("h_pT",P67,T7)-XSteam("hV_p",P67));
%preheater
Qp=m_dot_vapor*XSteam("CpL_p", P67)*(XSteam("Tsat_p",P67)-T6);
%evaporator
Qe=m_dot_vapor*(XSteam("hV_p",P67)-XSteam("hL_p",P67));
%total
Q_tot=Qe+Qp+Qs

%%%%Calculations for lengths
%%% superheater
%known/easily found stuff
pr_air=.690;
pr_vapor=1.05; %take average between the two 
Dh=(80-50)*10^(-3);
k_air=26.3e-3;
k_vapor=(40.5e-3);
Cp_air=1.004; %assuming air is an ideal gas
Cp_vapor=(XSteam("CpV_p",P67)+XSteam("Cp_pT",P67,180))/2;
density_air=1.169;
density_vapor=(XSteam("rho_pT",P67,T7)+XSteam("rhoV_p",P67))/2;
N=200;
A_air_s=pi*(25e-3)^2;
A_vapor_s=(pi*((40e-3)^2-(25e-3)^2));
dynamic_air=2.573e-5;
dynamic_vapor=(XSteam("my_pT",P67, T7)+XSteam("my_pT", P67, 180))/2;

T4b=T4-(Qs/(m_dot_air*Cp_air));
U_m_air_s=m_dot_air/(density_air*(A_air_s)*N);
U_m_vapor_s=m_dot_vapor/(density_vapor*(A_vapor_s)*N);
Re_air_s=((density_air*U_m_air_s*50e-3)/dynamic_air);
Re_vapor_s=((density_vapor*U_m_vapor_s*(30e-3))/dynamic_vapor);
Nu_air_s=(.023)*(Re_air_s)^(4/5)*(pr_air)^(.4);
Nu_vapor_s=(.023)*(Re_vapor_s)^(4/5)*(pr_vapor)^(.3);
h_vapor_s=(Nu_vapor_s*k_vapor)/(30*10^(-3));
h_air_s=(Nu_air_s*k_air)/(50e-3);
U_s=((1/h_air_s)+(1/h_vapor_s))^-1;
Delta_T1_s=T4b-XSteam("Tsat_p",P67);
Delta_T2_s=T4-T7;
Delta_Lm_s=(Delta_T1_s-Delta_T2_s)/log(Delta_T1_s/Delta_T2_s);

Ls=(Qs*10^3)/(U_s*Delta_Lm_s*N*(pi*(50e-3)));


%%%Evaporator
%Pool boiling so only air needs to be accounted for in h for U
%same process as superheater
A_vapor_e=A_vapor_s;
A_air_e=A_air_s;
T4c=T4b-(Qe/m_dot_air*Cp_air);
U_m_air_e=U_m_air_s;
Re_air_e=Re_air_s;
Nu_air_e=Nu_air_s;
h_air_e=h_air_s;
U_e=(1/h_air_e)^-1;
Delta_T1_e=T4c-XSteam("Tsat_p",P67);
Delta_T2_e=T4b-XSteam("Tsat_p",P67);
Delta_Lm_e=(Delta_T1_e-Delta_T2_e)/log(Delta_T1_e/Delta_T2_e);

Le=(Qe*10^3)/(U_e*Delta_Lm_e*N*(pi*(50e-3)));

%Still need to check that critical heat flux is below the threshold
g=9.81;
Cz=.131;
hfg=2015.29;
sigma_lv=42.9e-3;
rho_v=XSteam('rhoV_p',P67);
rho_l=XSteam('rhoL_p',P67);
crit_heat=(Cz*hfg*rho_v)*(((sigma_lv*g*(rho_l-rho_v))/(rho_v^2))^(1/4));

%actual heat flux
h_flux=Qe/(N*pi*(50e-3)*Le);

%%%Preheater
%Need to use epsilon-NTU
%know it is a shell and tube with one shell and 2 passes
Cp_water_p=(XSteam('CpL_p',P67)+XSteam('Cp_pT',P67,T6))/2;
density_water_p=(XSteam('rho_pT',P67,T6)+XSteam('rhoL_p',P67))/2;
dynamic_water_p=(XSteam('my_pT',P67, T6)+XSteam('my_pT', P67, 179.5))/2;
k_water_p=(.673+(.645))/2;
pr_water_p=2.384;
%

Qp_max=Cp_water_p*m_dot_vapor*(T4c-T6);
C_r=(Cp_water_p*m_dot_vapor)/(Cp_air*m_dot_air);
%Q_p2=m_dot_air*Cp_air*(T4c-T5);
epsilon=Qp/Qp_max;
E=((2/epsilon)-(1+C_r))/((1+C_r^2)^.5);
NTU=(-(1+C_r^2)^(-.5))*log((E-1)/(E+1));
Dp=20e-3;
SD=1.5*Dp;
D_square=(19*SD)+Dp;
D_shell= sqrt(D_square^2+D_square^2);
bs=.35*D_shell;
Sm=(bs)*(SD-Dp);
Gs=(m_dot_air)/(Sm);
De=(4*(1*SD^2-((pi*Dp^2)/4)))/(pi*Dp); %assuming square pitch so Cpd=1
Res_air_p=(De*Gs)/(dynamic_air);
Nu_air_p=.36*(Res_air_p^.55)*(pr_air^(1/3));
h_air_p=(k_air*Nu_air_p)/Dp;
Re_water_p=(4*m_dot_vapor)/(pi*Dp*dynamic_water_p*N);
Nu_water_p=.023*(Re_water_p^(4/5))*(pr_water_p)^(.4);
h_water_p=(Nu_water_p*k_water_p)/(Dp);
U_p=((1/h_air_p)+(1/h_water_p))^-1;

Lp=(NTU*Cp_water_p*m_dot_vapor*1000)/(U_p*pi*Dp*N)


