%function Rl = Longwave_Rad(Rs,dew_temp)
% CONSTANTS
global bi bs C ca ci le lf mui mus B di

bi = 0.5;           % thickness interval of ice for heat conduction and irradiance absorption (m)
bs = 0.1;           % thickness interval of snow for heat conduction and irradiance absorption (m)
C = 0.002;          % bulk coefficient for sensible and latent heat
ca = 1006;          % specific heat of air (J/kg-K)
ci = 2100;          % specific heat of ice (J/kg-K)
le = 2.50*10^6;     %latent heat of evaporation of ice (J/kg)
lf = 3.34*10^5;     %latent heat of fusion of ice (J/kg)
mui = 10;           %extinction coefficient of ice (/m) 
mus = 40;           %extinction coefficient of snow (/m)
B = 5.67*10^-8;     %stefan-boltzmann constant (W/m^2-K^4)
di = 900;           %density of ice (kg/m^3)

function Ks = Snow_Thermal_Conductivity(ds)
    Ks = 0.029( 1 + 0.001*(ds^2))
end

function Ts = Surface_Temp(Ta,A,Rs,e,Rl,da,U,rh,q,Hg,G_q)
    Ts = Ta + ((1-A)*Rs + e*Rl - e*B*(Ta + 273.2)^4 - le*da*C*U*(1-rh)*q*Ta + Hg)/(4*e*B*(Ta + 273.2)^3 + (G_q*le + ca)*da*C*U);
end

function Hg = Heat_Into_Glacier(ds,A,Rs,Ts)
    if Ts = 0
        Hg = 0
    elseif Ts<0
        Ks = Snow_Thermal_Conductivity(ds)
        Hg = -Ks*Ts/bs
    %elseif no water in snow
    % then integration formula for Hg
    end
