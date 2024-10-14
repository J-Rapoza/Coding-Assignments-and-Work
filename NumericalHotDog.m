clear, close all
clc

%% Given Parameters
% Constants
rho = 880 ;        % kg/m3
k = 0.52 ;       % W/mK
c = 3350 ;       % J/kg-K
D = .0254  ;      % m
r = D/2;
alpha = k/c/rho; %thermal diffusivity [m2/s]

% Initial Constraints
Ti = 10; % C
Tsmax = 100 ; % C
Tinf = 250;  % C
Tcoal = 450; % C
ecoal = .8;
ehd = .45;

% Finite elements and iterations
m = 100; % such that finite elements number from 0 to m
dr = r/m;
dt = 0.1;
Tvec = ones(m+1, 1)*Ti; %initial temperature distribution
n = 4720; %number of iteration
time = (1:n)*dt;

%% Initialization
Center_T = zeros(1,n); %Temperature of centerline
Surface_T = zeros(1,n); % Temperature of surface
T30 = zeros(1,m+1); %temperature distribution at 30s
T120 = zeros(1,m+1); %temperature distribution at 120s
TEnd = zeros(1,m+1); %final temperature
hconvvec = zeros(1,n);
hradvec = zeros(1,n);

% Properties of air around the coal -- used to calculate hconv
g = 9.81;
Tfcoal = (Tinf+Tcoal)/2;
beta = 1/Tfcoal;
nu_coal = findNu(Tfcoal); %interpolate the viscosity of air
D_coal = 0.0254; %diameter of coal -- same as hotdog
Gr = g*beta.*(Tcoal-Tinf)*D_coal^3./nu_coal^2;
Re_coal = Gr^0.5;
vaircoal = Re_coal*nu_coal/D_coal;


for nn = 1:n
    Center_T(nn) = Tvec(1); 
    Surface_T(nn) = Tvec(end);%Record temp at this time step
    
    %% Calculate convection heat transfer coefficient
    % Things happened to the hotdog
    vairdog = vaircoal/2;
    Tfhd = (Tvec(end)+Tinf)/2; %Film temperature around hotdog
    nu_HD = findNu(Tfhd); %vsicosity around HD
    Pr = findPr(Tfhd); %Prandtl number of air
    k_air = findK(Tfhd);
    ReHD = vairdog*D/nu_HD;
    Nu = 0.3+0.62*ReHD^(1/2)*Pr^(1/3)/(1+(0.4/Pr)^(2/3))^(1/4)...
        *(1+(ReHD/282000)^(5/8))^(4/5);
    hconvvec(nn) = Nu*k_air/D;
    
    %% Calculate Radiation Heat Transfer Coefficients
    T1 = 450+273;
    T2 = Tvec(end)+273;
    A1 = 0.6; % width of grill
    A2 = pi*D; % perimeter of cylinder
    s1 = A1/2;
    s2 = -A1/2;
    y = 0.098; % distance from hot dog to coal
    F12 = r/(s1-s2)*(2*atan(s1/y));
    F21 = A1*F12/A2;
    T2_star = T2./(ecoal*F21)^(1/4);
    sigma = 5.67e-8;
    hradvec(nn) = ecoal*ehd*F21*sigma*(T1+T2_star).*(T1^2+T2_star.^2);
    
    h_tot = hconvvec(nn)+hradvec(nn);
    
    %% Calculate Temperature Change
    Fo = alpha*dt/dr^2;
    Foh = h_tot*dt/rho/c/dr;
    heatM = zeros(m+1,m+1);
    heatM(1, [1 2]) = heatM(1, [1 2])+[4*Fo+1 -4*Fo];
    heatM(m+1, [m m+1]) = [-4*Fo*(2*m-1)/(4*m-1)...
        (4*(2*m-1)*Fo+8*m*Foh+4*m-1)/(4*m-1)];
    for ii = 1:m-1
        heatM(ii+1, ii:ii+2) = [-(1-1/(2*ii))*Fo (2*Fo+1) -(1+1/(2*ii))*Fo];
    end

    Tvec(end) = Tvec(end)+8*m/(4*m-1)*Foh*Tinf;
    Tvec = heatM\Tvec;
    if nn*dt==30
        T30 = Tvec;
    elseif nn*dt == 120
        T120 = Tvec;
    end
end
TEnd = Tvec;
Tvec(1)
Tvec(end)

% %% Temp vs Time Graph
% temp_t_figure = figure;
% centerP = plot(time, Center_T, "LineWidth", 1.5);
% hold on
% surfaceP = plot(time, Surface_T, "LineWidth", 1.5);
% title('Temperature Over Time')
% xlabel('Time (s)')
% xlim([0 time(end)])
% ylabel('Centerline Temperature (C)')
% legend("Centerline Temp", "Surface Temp", "location", "best")
% hold off
% saveas(gcf, 'TempOverTime_num.png')
% 
 %% Graph Temp as a function of r
% rvec = linspace(0,r,m+1);
% temp_r_figure = figure;
% T30p = plot(rvec, T30, "LineWidth", 1.5);
% hold on
% T120p = plot(rvec, T120, "LineWidth", 1.5);
% TEndp = plot(rvec, TEnd, "LineWidth", 1.5);
% xlabel('Distance (m)')
% ylabel('Temperature (C)')
% legend("T = 30s", "T = 2min", "Final Time", "location", "best")
% xlim([0 r])
% ylim([0 110])
% saveas(gcf, 'TempDistribution_num.png')

function nu = findNu(Temp)
    % input Temp has a unit of C and need to be converted to K
    knownTemp = 200:50:700;
    knownNu = [7.59 11.44 15.89 20.92 26.41 32.39 38.79 45.57 52.69...
        60.21 68.1]/10^6;
    nu = interp1(knownTemp, knownNu, Temp+273.15);
end

function k = findK(Temp)
    % input Temp has a unit of C and need to be converted to K
    knownTemp = 200:50:700;
    knownK = [18.1 22.3 26.3 30.0 33.8 37.3 40.7 43.9 46.9 49.7 52.4]/10^3;
    k = interp1(knownTemp, knownK, Temp+273.15);
end

function Pr = findPr(Temp)
    % input Temp has a unit of C and need to be converted to K
    knownTemp = 200:50:700;
    knownPr = [0.737 0.720 0.707 0.700 0.690 0.686 0.684 0.683 0.685...
        0.690 0.695];
    Pr = interp1(knownTemp, knownPr, Temp+273.15);
end
