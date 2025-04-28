%% Finding T-v & T-ω for DUSC2025

CdA = 0.05836913;
rho_air = 1.292; %kg/m^3
tyre_rad = 0.588./2;
m_car = 270;
g = 9.81;

top_speed = 162;

v_car_kph = linspace(0,top_speed,top_speed.*10 + 1);

v_car = v_car_kph ./ 3.6;

w = v_car./tyre_rad;

mu_RR = 0.0085;


F_Drag = 0.5 .* rho_air .* v_car.^2 .* CdA;

F_RR = mu_RR .* m_car .* g;

F_total = F_Drag + F_RR;

T_car = F_total .* tyre_rad;

%plot(v_car_kph, T_car);



