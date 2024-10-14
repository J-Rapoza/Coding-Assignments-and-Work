clear, close all
clc

% Heat Transfer Final Project
%% Given Parameters
% Constants
    rho = 880 ;        % kg/m3 
    k = 0.52 ;       % W/mK
    c = 3350 ;       % J/kg-K 
    D = .0254  ;      % m 
% Iinitial Constraints
    Ti = 10; % C
    Tsmax = 100 ; % C
    Tinf = 250;  % C
    Tcoal = 450; % C
    ecoal = .8;
    ehd = .45;
    time = 0:455;
% Heat Transfer values
    hconv = 9.06; % W/m^2 * K
    hrad = 7.64;  % W/m^2 * K 
    h_tot = hconv+hrad; % total effective heat transfer rate
    
%% Calculating Relevant Coefficients
% Biot and Fourier Number 
    Bio = (h_tot*(D/2))/k ;
    Fo = ((k/(rho*c))*(time))/((D/2)^2) ;

% Find eta's
    etaLength = 100; % self-defined: multiple eta to derive correct
                      % initial temperature
    etaVec = zeros(1, etaLength);
    for n = 1:etaLength
        syms eta_n
        eqn = eta_n*besselj(1, eta_n)/besselj(0, eta_n) == Bio;
        guess = 0.1+(n-1)*pi;
        etaVec(n) = vpasolve(eqn, eta_n, guess);
    end

% Find C's
    CVec = zeros(1, etaLength); %vector of coefficients, same length
    for n = 1:etaLength
        CVec(n) = 2/etaVec(n)*besselj(1, etaVec(n))...
                    /(besselj(0, etaVec(n))^2+besselj(1, etaVec(n))^2);
    end

%% Infinite Cylinder Approximate Solution Equations
    Surface_Theta = zeros(1, length(time));
    Center_Theta = zeros(1, length(time));
    for n = 1:length(CVec)
        Surface_Theta = Surface_Theta+CVec(n)*exp(-etaVec(n)^2*Fo)...
                                      *besselj(0, etaVec(n));
        Center_Theta = Center_Theta+CVec(n)*exp((-(etaVec(n)^2)*(Fo)));
    end

% Caluclating Temperatures from thetas
    Surface_T = (Surface_Theta)*(Ti-Tinf) + (Tinf);
    Center_T = (Center_Theta)*(Ti-Tinf) + (Tinf);

%% Temp vs Time Graph
    temp_t_figure = figure;
    centerP = plot(time, Center_T, "LineWidth", 1.5);
    hold on
    surfaceP = plot(time, Surface_T, "LineWidth", 1.5);
    title('Temperature Over Time -- Analytical')
    xlabel('Time (s)')
    xlim([0 time(end)])
    ylabel('Centerline Temperature (C)')
    legend("Centerline Temp", "Surface Temp", "location", "best")
    hold off
    saveas(gcf, 'TempOverTime.png')
    
%% Find Relationship between Temp and r
    rStar = 0:0.01:1;
% When t = 30s
    Fo30 = Fo(time==30);
    theta30 = zeros(1, length(rStar)); %find shape distribution
    for n = 1:length(rStar)
        for m = 1:length(CVec)
            theta30(n) = theta30(n)+CVec(m)*exp(-etaVec(m)^2*Fo30)...
                                           *besselj(0, etaVec(m)*rStar(n));
        end
    end
    T30 = theta30*(Ti-Tinf)+Tinf;
    
% When t = 2min = 120s
    Fo120 = Fo(time==120);
    theta120 = zeros(1, length(rStar)); %find shape distribution
    for n = 1:length(rStar)
        for m = 1:length(etaLength)
            theta120(n) = theta120(n)+CVec(m)*exp(-etaVec(m)^2*Fo120)...
                                           *besselj(0, etaVec(m)*rStar(n));
        end
    end
    T120 = theta120*(Ti-Tinf)+Tinf;

% When t = final time
    FoEnd = Fo(end);
    thetaEnd = zeros(1, length(rStar)); %find shape distribution
    for n = 1:length(rStar)
        for m = 1:length(etaLength)
            thetaEnd(n) = thetaEnd(n)+CVec(m)*exp(-etaVec(m)^2*FoEnd)...
                                           *besselj(0, etaVec(m)*rStar(n));
        end
    end
    TEnd = thetaEnd*(Ti-Tinf)+Tinf;

%% Graph Temp as a function of r
    r = rStar*D/2;
    temp_r_figure = figure;
    T30p = plot(r, T30, "LineWidth", 1.5);
    hold on
    T120p = plot(r, T120, "LineWidth", 1.5);
    TEndp = plot(r, TEnd, "LineWidth", 1.5); 
    title("Temperature distribution along radius - Analytical")
    xlabel('Distance (m)')
    ylabel('Temperature (C)')
    legend("T = 30s", "T = 2min", "Final Time", "location", "best")
    xlim([0 r(end)])
    ylim([0 110])
    saveas(gcf, 'TempDistribution120Seconds.png')









