% CONSTANTS

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


% Observed data from ECMWF website

prompt1 = "Air Temperature = ";                          % deg celcius 
Ta = input(prompt1);
prompt2 = 'Relative humidity = ';                        
rh = input(prompt2);
prompt3 = 'Observed downward shortwave radiation = ';    % W/m^2
Rs = input(prompt3);
prompt4 = 'Wind speed = ';                               % m/s
U = input(prompt4);
prompt5 = 'Daily precipitaion = ';                        % mm w.e per day
Pp = input(prompt5);

prompt7 = 'Emissivity of snow/ice surface = ';
e = input(prompt7);
prompt8 = 'Downward longwave radiation = ';
Rl = input(prompt8);
prompt9 = 'Density of air = ';
da = input(prompt9);
prompt10 = 'Saturated specific humidity = ';
q = input(prompt10);
prompt11 = 'Gradient of saturated specific humidity = ';
G_q = input(prompt11);
% q and G_q are functions of air temperature


% initial conditions

% snow density
% Thickness of snow layer
% Ice temperature profile


% main input

prompt6 = 'Albedo = ';
A = input(prompt6);


Hg = 0; % heat transfer into the glacier

Ts1 = Ta + ((1-A)*Rs + e*Rl - e*B*(Ta + 273.2)^4 - le*da*C*U*(1-rh)*q*Ta + Hg)/(4*e*B*(Ta + 273.2)^3 + (G_q*le + ca)*da*C*U);

% using the above value of Ts, we calculate the new value of Hg


Rn = (1-A)*Rs + e*Rl - e*B*(Ts + 273.2)^4;





