% This file is used to get the value of specific run off and surface
% temperature of glacier for a specific value of albedo
% This also gives the variation of snow and ice temperatures with depth
% into the glacier and during the duration of observation

global bi bs C ca ci le lf mui mus B di Ta rh Rs U Pp e Rl da q G_q t ds xsnow xice A

% CONSTANTS

bi  = 0.05;            % thickness interval of ice for heat conduction and irradiance absorption (m)
bs  = 0.1;             % thickness interval of snow for heat conduction and irradiance absorption (m)
C   = 0.002;           % bulk coefficient for sensible and latent heat
ca  = 1006;            % specific heat of air (J/kg-K)
ci  = 2100;            % specific heat of ice (J/kg-K)
le  = 2.50*10^6;       % latent heat of evaporation of ice (J/kg)
lf  = 0.334*10^5;      % latent heat of fusion of ice (J/kg)
mui = 10;              % extinction coefficient of ice (/m) 
mus = 0;               % extinction coefficient of snow (/m)
B   = 5.67*10^-8;      % stefan-boltzmann constant (W/m^2-K^4)
di  = 900;             % density of ice (kg/m^3)

% Observed data from Various Citations (mainly from Fujita and Ageta Paper) 

Ta    = -5;        % Air temperature (degree Celcius)
rh    = 0.779;         % Relative Humidity
Rs    = 280;           % Observed downward shortwave radiation (W/m^2)
U     = 4.1;           % Wind speed (m/s)
Pp    = 670;           % Daily precipitaion (mm.w.e.d^-1)
e     = 0.91;          % Emissivity of snow/ice surface
Rl    = 307;           % Downward longwave radiation
da    = 1.35;          % Density of air (kg/m^3) from Fujita and Ageta Paper
q     = 0.014;         % Saturated specific humidity- from kondo paper
G_q   = 1;             % assumption- specific humidity linearly varies with temperature 
ds    = 400;           % density of snow (kg/m^3)from Fujita and Ageta Paper

t     = linspace(0, 86400, 3600);     % duration of the observation, taken as 1 day = 86400s, divided into 3600 divisions
xsnow = linspace(0, 0.0075, 50);        % snow layer thickness (m)
xice  = linspace(0, 0.1925, 1300);     % ice layer thickness (m)

% main input- Albedo:-

prompt1 = 'Albedo = ';
A = input(prompt1);
%A = 0.85;

%-----------------------------Literature----------------------------------%

% mi, mi1 = Mass of water due to melting of ice (kg) (both mi and mi1)
% mi0 = Initial Mass of water due to melting of ice (kg)

% Hg  = Heat of into glacier (W/m^2)
% Hg0 = Initial Heat of into glacier (W/m^2)

% Ts, Ts1 = Surface temperature of glacier (degree celcius)
% Ts0 = Initial Surface temperature of glacier (degree celcius)

% Fs = Mass of superimposed ice on ice surface below snow (mm.w.e)
% Fc = Mass of water that is refrozen within snow layer (mm.w.e)
% Fc_capillary = Mass of capillary water that is refrozen within wet snow(mm.w.e)

% Pr = Rainfall (mm.w.e.d^-1)
% Ps = Snowfall (mm.w.e.d^-1)
% DS = Runoff   (mm.w.e.d^-1)

%------------------Initialisation------------------------------------------%

mi0 = 0; %initially there won't be any melt water ie, in the initial stage, ice won't melt

Hg0 = 0; % Initial Heat transfer into the glacier (W/m^2)

% Calculation of Ts = Surface Temperature (degree Celcius)

Ts0 = Surface_Temp(A, Rs, Rl, e, Ta, da, U, rh, q, G_q, Hg0) %with Hg0 = 0

mi1 = Melted_ice(Ts0, A, Rs, Rl, e, da, U, Ta, rh, q, Hg0)

Hg = Heat_Into_Glacier(mi1, Ts0, t, ds, xsnow, xice) %with initial conditions

Ts1 = Surface_Temp(A, Rs, Rl, e, Ta, da, U, rh, q, G_q, Hg) %with obtained Hg


%-----------------------------Code By Aryan-------------------------------%


while(abs(Ts1 - Ts0) >= 0.1)
    mi1 = Melted_ice(Ts1, A, Rs, Rl, e, da, U, Ta, rh, q, Hg)
    Hg = Heat_Into_Glacier(mi1, Ts1, t, ds, xsnow, xice)
    Ts0 = Ts1 %1st loop Ts
    Ts1 = Surface_Temp(A, Rs, Rl, e, Ta, da, U, rh, q, G_q, Hg) %2nd loop Ts
end
mi = mi1
Ts = Ts1


%-----------------------------Mass Balance--------------------------------%
%----------------------------Code By Rikitha------------------------------%
Fc_capillary=0;
Fc=0;
Fs=0;
if Ts>=0 %ice melts but some water might be refrozen in 3 forms
    %mi
    Fs = Superimposed_ice(ds, xice);
    Fc_capillary = Mass_of_Refrozen_Capillary_Water(ds, Ts);
    Fc = Mass_of_Snow_Percolated_Water(xsnow, t, ds);
%else % Ts<0 => no meltwater   
end

%Precipitations
[Pr,Ps] = Precipitations(Pp,Ta); %rainfall and snowfall

%Final Output of the code:-
DS = Pr - (Fc_capillary + Fs + Fc) + mi; % Runoff from the glacier
disp(strcat('The final runoff  obtained is  ',num2str(DS),' mm.w.e'))

%--------------------------------Functions--------------------------------%

%1. Surface_Temp(A, Rs, Rl, e, Ta, da, U, rh, q, G_q, Hg)
% Returns the value of surface temperature in degree Celcius
function Ts = Surface_Temp(A, Rs, Rl, e, Ta, da, U, rh, q, G_q, Hg)
    global C ca le B
    Ts = Ta + ((1-A)*Rs + e*Rl - e* B *(Ta + 273.2)^4 - le*da*C*U*(1-rh)*q*Ta + Hg)/(4*e*B*(Ta + 273.2)^3 + (G_q*le + ca)*da*C*U);
    if Ts>0
        Ts = 0;
    end
end

%2. Melted_ice(Ts, A, Rs, Rl, e, da, U, Ta, rh, q, Hg)
% Return the value of ice that is melted in the glacier in mm.w.e.
function mi = Melted_ice(Ts, A, Rs, Rl, e, da, U, Ta, rh, q, Hg)
   Rn = Net_Radiation(A, Rs, Rl, e, Ts);
   Hs = Sensible_heat(da, U, Ta, Ts);% see the function comments for actual parameters
   Hl = Latent_Heat_Flux(da, U, rh, q, Ta, Ts);
   Hm = Rn+Hs+Hl+Hg; 
   L = 3.36e+5;  %latent heat of fusion of ice (J/kg)
   mi = Hm/(L*1000);
   if mi<0
       mi = 0;
   end
end

%3. Heat_Into_Glacier(mi, miprev, Ts, t, ds, xsnow, xice)
% Returns the value of Hg = heat into the glacier in W/m^2
function Hg = Heat_Into_Glacier(mi, Ts, t, ds, xsnow, xice)
    global bs ci
    if mi == 0 || Ts<0 % if no ice is melted ie, when surface temperature is negative
        [temp_s, t] = Snow_Temperature();
        [temp_ice, t] = Ice_Temperature();
        snow_trap = trapz(temp_s(end,:)); %integration of the temperature variation of snow with limits of depth of snowlayer- from snow-ice interface to glacier surface
        ice_trap = trapz(temp_ice(end,:)); %integration of the temperature variation of ice with limits of depth of snowlayer- from critical depth - zc(from Fujita paper) to snow-ice interface
        Hg = -1*ci*((ice_trap + snow_trap))/t(end) ;
    elseif mi == 0 || Ts == 0 % if no ice is melted but surface temperature = 0 degree celcius
        Hg = 0;
    elseif Ts<0 && mi>0
        Ks = Snow_Thermal_Conductivity(ds);
        Hg = -Ks*Ts/bs;
    %elseif no water in snow
    % then integration formula for Hg
    end
end

%4. Net_Radiation(A, Rs, Rl, e, Ts)
% Returns net radiation (Rn) in W/m^2 to the function - Melted_ice
function Rn = Net_Radiation(A, Rs, Rl, e, Ts)
    global B
    Rn = ((1-A)*Rs)+e*Rl - e*B*(Ts + 273.2).^4;
end

%5. Sensible_heat(da, U, Ta, Ts)
% Returns Sensible heat (Hs) in W/m^2 to the function - Melted_ice
function Hs = Sensible_heat(da, U, Ta, Ts)
    global ca C
    Hs = ca*da*C*U*(Ta-Ts);
end

%6.Latent_Heat_Flux(da, U, rh, q, Ta, Ts)
% Returns Latent heat (Hl) in W/m^2 to the function - Melted_ice
function Hl = Latent_Heat_Flux(da, U, rh, q, Ta, Ts)
    global le C
    Hl = le*da*C*U*((rh*q*Ta)-(q*Ts));
end

%7. Superimposed_ice(ds, xice)
% Returns Fs = amount of superimposed ice on ice layer below snow layer in mm.w.e
function Fs = Superimposed_ice(ds, xice)
    global di ci lf
    [temp_ice, t] = Ice_Temperature();
    ice_trap = trapz(temp_ice(end,:)); %integration of the temperature variation of ice with limits of depth of icelayer- from critical depth - zc(from Fujita paper) to snow-ice interface
    Fs = di * ci *(ice_trap)/lf ; 
end

%8. Mass_of_Refrozen_Capillary_Water(ds, Ts)
% Returns Fc_capillary =  Mass of capillary water that is refrozen within wet snow(mm.w.e)
function Fc_capillary = Mass_of_Refrozen_Capillary_Water(ds, Ts)
global lf bs    
    Ks = Snow_Thermal_Conductivity(ds); % Themal conductivity of snow in W/mK
    Fc_capillary = (-1*Ks*Ts)/(lf*bs);
end

%9. Mass_of_Snow_Percolated_Water(xsnow, t, ds)
% Returns Fc = Mass of water that is refrozen within snow layer (mm.w.e)
function Fc = Mass_of_Snow_Percolated_Water(xsnow, t, ds)
    global ci lf    
    [temp_s, ~] = Snow_Temperature();%should use for this function- xice and t
    snow_trap = trapz(temp_s(end,:));%integration of the temperature variation of snow with limits of depth of snowlayer- from snow-ice interface to glacier surface
    Fc = ds*ci*(snow_trap)/lf; 
end


%10. Snow_Temperature(ds, xsnow, t, Ts)
% Returns the variation of temperature of snow layer with time and depth
function [temp_s,t] = Snow_Temperature()
    m = 0;    
    x = linspace(0, 0.75, 50); %snow thickness
    t = linspace(0, 86400, 3600);

    sol_snow = pdepe(m, @temp_var_snow, @snow_ic, @snow_bc, x, t); % partial differential equation

    temp_s = sol_snow(:,:,1);
    
    % plots of snow temperature with time and depth
    
    figure, plot(x,temp_s(end,:))
    title(strcat('Temperature vs Depth for Snow Layer at t =  ', num2str(86400),'s'))
    xlabel('Depth (cm)')
    ylabel('Temperature (C)')

    figure, plot(t,temp_s(:,1))
    title('Surface Temperature variation with Time')
    xlabel('Time (s)')
    ylabel('Temperature (C)')
end

%11. temp_var_snow(x,t,u,dudx) - heat variation in snow
function [c,f,s] = temp_var_snow(x,t,u,dudx)
    global ds ci A Rs mus
    c = ds*ci;
    f = 0.493*dudx;
    s = (1-A)*Rs*exp(-1*mus*x);
end

%12. snow_ic() - initial conditions of snow
% Returns initial conditions of the differential question to the function
% Snow_Temperature(ds, xsnow, t, Ts) for pde
function u0 = snow_ic(x)
    u0 = 0;
end

%13. snow_bc(xl,ul,xr,ur,t,Ts) - boundary conditions of snow
% Returns left and right boundary conditions of the differential question
% to the function Snow_Temperature(ds, xsnow, t, Ts) for pde
function [pl,ql,pr,qr] = snow_bc(xl,ul,xr,ur,t)
    pl = ul; % write (ul - Ts0) here
    ql = 0;
    pr = ur + 5; %-5 C from Fujita Paper
    qr = 0;
end

%14. Ice_Temperature(xice, t)
% Returns the variation of temperature of ice layer with time and depth
% upto critical depth (zc)
function [temp_ice, t] = Ice_Temperature()
    m = 0;
    x = linspace(0, 19.25, 1300); %ice thickness
    t = linspace(0, 86400, 3600);

    sol = pdepe(m, @temp_var_ice, @ice_ic, @ice_bc, x, t);% partial differential equation

    temp_ice = sol(:,:,1);

    % plots of snow temperature with time and depth
    
    figure, plot(x,temp_ice(end,:))
    title(strcat('Temperature vs Depth for Ice Layer at t = ', num2str(86400),'s'))
    xlabel('Depth (cm)')
    ylabel('Temperature (C)')

    figure, plot(t,temp_ice(:,1))
    title('Snow-Ice Interface Temperature variation with Time')
    xlabel('Time (s)')
    ylabel('Temperature (C)')
end


%15. temp_var_ice(x,t,u,dudx) - heat variation in ice
function [c,f,s] = temp_var_ice(x,t,u,dudx)
    global di ci
    c = di*ci*(273.2 + u)/(616.6 + 0.47*u);
    f = dudx;
    s = 0;
end

%16. ice_ic() - initial conditions of ice
% Returns initial conditions of the differential question to the function
% Ice_Temperature(ds, xice, t, Ts) for pde
function u0 = ice_ic(x)
    u0 = -5;
end

%17. ice_bc(xl,ul,xr,ur,t) - boundary conditions of ice
% Returns left and right boundary conditions of the differential question
% to the function Ice_Temperature(ds, xsnow, t, Ts) for pde
function [pl,ql,pr,qr] = ice_bc(xl,ul,xr,ur,t)
    pl = ul + 5;
    ql = 0;
    pr = ur + 7;
    qr = 0;
end

%18.Snow_Thermal_Conductivity(ds)
% Returns the thermal conductivity of snow in (W/mK)
function Ks = Snow_Thermal_Conductivity(ds)
    Ks = 0.029*(1+(0.001*ds*ds));
end

%19. Precipitations(Pp, Ta)
% Returns the rainfall and snowfall in mm.w.e/d
% Pr = Rainfall, Ps = Snowfall
function [Pr, Ps] = Precipitations(Pp, Ta) %units- mm.w.e/d
    if Ta<=0
        Ps = Pp;
    elseif Ta<6
        Ps = Pp * (1-(Ta/6));
    else
        Ps = 0;
    end
    Pr = Pp - Ps;
end
