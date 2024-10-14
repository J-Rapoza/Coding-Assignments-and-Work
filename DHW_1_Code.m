clc
clear all

% known variables
P23=10;
P14=.07384;
x4=.88;
T3U=400;
T1=40;
h1=XSteam('hL_T',T1);
nt=.9;
Qb=30*10^6;
S1=.5724;
D=25*10^-3;
L=10;
U=1000;
Tci=18.3;
Tco=21.1;
epsilon=.0015*10^-3;
Cp=4.18*10^3;
rho=999;
mu=1.12*10^-3;

%solving for lower bound of turbine inlet temp
s4f=XSteam('sL_p',P14);
s4g=XSteam('sV_p',P14);
s4=(1-x4)*s4f+x4*s4g;
T3L=XSteam('T_ps',P23,s4);

%initialize output array
temp=[];
eff=[];
var_min=T3L;
var_max=T3U;
step=.1;

for T3= var_min:step:var_max
    %solving for the boiler
    T2=XSteam('T_ps',P23,S1);
    h2=XSteam('h_ps',P23,S1);
    h3=XSteam('h_pT',P23, T3);
    m_dot=Qb/(h3-h2);
    
    % solving the turbine portion of the problem
    S3=XSteam('s_pT',P23,T3);
    h4s=XSteam('h_ps',P14,S3);
    h4=nt*(h3-h4s)+h3;
    Wt=m_dot*(h4-h3);
    
    %solving for Rankine pump
    Wp1=m_dot*(h2-h1);
    
    %energy balance to find Qc
    Qc=Qb+Wp1-Wt;
    
    %solving for condenser
    delta_T1=abs(T1-Tci);
    delta_T2=abs(T1-Tco);
    delta_T=(delta_T1-delta_T2)/(log(delta_T1/delta_T2));
    N=Qc/(U*delta_T*pi*D*L);
    
    
    %finding major loss and friction coefficient
    m_dot_cooling_total=Qc/(Cp*(Tco-Tci));
    m_dot_cooling=m_dot_cooling_total/N;
    V=(m_dot_cooling)/(rho*pi*(D/2)^2);
    Reynolds= (rho*V*D)/mu; 
    f=((1)/(-1.8*log10(((epsilon/D)/3.7)^1.11*(6.9/Reynolds))))^2;
    Maj_loss=((f*L*rho*V^2)*N)/(2*D);
    
    %cooling pump
    flow_rate=V*pi*(D/2)^2;
    Wp2=flow_rate*Maj_loss;
    
    %overall rankine efficiency
    nth= (Wt-Wp1-Wp2)/Qb;
    
    %output array building
    temp=[temp,T3];
    eff=[eff,nth];
end

%efficiency plot
plot(temp,eff, "blue");
hold on
plot(T3,nth,"red*");
legend("Efficiency as a Function of Temperature","Maximum Efficiency")
xlabel("Turbine Inlet Temperature (°C)");
ylabel("Rankine Efficiency");
title("Rankine Cycle Efficiency as a Function of Turbine Inlet Temperature");
grid on;
