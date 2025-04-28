clc; clear; close all;

run('DUSC2025_TorqueSpeedCharacteristics.m');
% Define parameters
ri_values = [95]; % mm
n_p_values = (16);
gamma_values = (2);
mag_values = (12);
 % core gap in mm 
Maximum_Efficiencies = zeros(length(n_p_values), length(gamma_values), length(ri_values), length(mag_values), 3);
for ci = (1)
for m = 1:size(mag_values,2)
mag = mag_values(m);
for c = 1:size(ri_values,2)
ri = ri_values(c);
ro = 183; % mm
n_phi = 3;


%for z = 1:1:size(gamma_values)
for z = 1:1:size(n_p_values,2)
n_p = n_p_values(z);
theta_p = 2 * pi ./ n_p;
for gamma_index = 1:size(gamma_values,2)
gamma = gamma_values(gamma_index);
wd_max = theta_p * ri / (n_phi * gamma);
N = 2 * n_p .* gamma; % Number of turns per Phase
tyre_dia = 0.558; % in m
tyre_rad = tyre_dia./2; % in m

velocity_index = 751;

v_req_kph = v_car_kph(velocity_index);

v_req = v_req_kph / 3.6; % m/s

w_req = v_req ./ (0.279); % rad/s

n_s = 60 .* w_req ./ 2 ./ pi;

f = n_p .* n_s ./ 120; % Electrical Frequency (Hz)

theta_a_max = theta_p/n_phi/gamma  - (0.1*pi/180);
 % theta_a = pi ./ 180 .* 2.1222;
 % theta_c = pi ./ 180 .* 59;
theta_a = linspace(0,theta_a_max,20); % End winding
theta_c = pi ./ 180 .* (3:1:89)'; % Column vector

% Create a meshgrid for 2D data
[Theta_A, Theta_C] = meshgrid(theta_a, theta_c);



% Compute ri_prime
ri_prime = ri ./ ((1 + 0.5 .* (theta_p - Theta_A) .* tan(Theta_C)));
% ri_prime = ri.*(1 - tan(Theta_C).*(theta_p - Theta_A)./2);


% Compute wd numerator and denominator
% wd_num = 1 .* (ri_prime .* (theta_p - n_phi .* gamma .* Theta_A) .* sin(2 .* Theta_C));
% wd_den = n_phi .* gamma .* cos(Theta_C); %.* (1 + tan(Theta_C ./ 2) .* sin(2 .* Theta_C) + cos(2 .* Theta_C));
%wd_den = n_phi .* gamma .* cos(Theta_C) .* (1 + tan(Theta_C ./ 2) .* sin(2 .* Theta_C) + cos(2 .* Theta_C));
 
mu_1 = ri_prime .* (theta_p - n_phi .* gamma .* Theta_A) ./ (2 .* cos(Theta_C) .* n_phi .* gamma);
mu_2 = ri_prime .* ( 2.* theta_p - n_phi .* gamma .* Theta_A) ./ (2 .* cos(Theta_C) .* n_phi .* gamma);

A =  mu_1 .* sin(2.*Theta_C);
B =  ci .* Theta_A .* ri_prime .* sin(Theta_C);
T = tan(3 .* Theta_C ./ 2);
K = (1./2) + (1./(2.*T));
% wd_checker = -((ci - A) - sqrt((ci - A).^2 + 4.*(A.*ci+B)))./2;
% wd = min(-((ci - A) - sqrt((ci - A).^2 + 4.*(A.*ci+B)))./2, (mu_2 - mu_1).*sin(2.*Theta_C));
a = K;
b = (K .* ci - A);
ce = -(A .* ci + B);

% Solve for w_d
wd = min((-b + sqrt(b.^2 - 4.*a.*ce)) ./ (2.*a), (mu_2-mu_1).*sin(2.*Theta_C) );
wd_A = (mu_2-mu_1).*sin(2.*Theta_C);
wd_B = (-b + sqrt(b.^2 - 4.*a.*ce)) ./ (2.*a);
% Adjust for Theta_C >= 45 degrees
% wd_den(Theta_C >= (pi/4)) = n_phi .* gamma .* cos(Theta_C(Theta_C >= (pi/4))) .* (1 + tan(Theta_C(Theta_C >= (pi/4)) ./ 2) .* sin(2 .* Theta_C(Theta_C >= (pi/4))) - cos(2 .* Theta_C(Theta_C >= (pi/4))));

% Compute wd
% wd = wd_num ./ wd_den;
% 
% invalid_indices = (Theta_A .* ri_prime) < wd;
% ri_prime(invalid_indices) = NaN;
% Theta_A(invalid_indices) = NaN;
% wd(invalid_indices) = inf;
% theta_a(invalid_indices) = NaN;
% Finding Correct AWG with Insulation
%AWG = [1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22; 7.348, 6.544, 5.189, 4.115, 3.264, 2.588, 2.053, 1.628, 1.291, 1.024, 0.812, 0.644; 14.4 , 11.6, 10.6, 9.2, 6.8, 5.1, 4.3, 3.7, 3.3, 3, 2.8, 2.6];
%AWG = [22, 20, 18, 16, 14, 12; 0.644, 0.812, 1.024, 1.291, 1.628, 2.053; 2.6, 2.8, 3, 3.3, 3.7, 4.3]; % Shitty RS stuff 
%AWG = [22, 20, 18, 17, 14, 13, 12, 11, 10, 9 , 8; 0.644, 0.812, 1.024, 1.291, 1.628, 1.828, 2.053, 2.305, 2.588, 2.906, 3.264; 0.9652, 1.0668, 1.3716,  1.6764, 2.4384, 2.5908, 3.0734, 3.2258, 3.9878, 4.699, 5.3086]; % Litz Wire % AWG;CD;WD
AWG = [22, 20, 18, 17, 14, 13, 12, 11, 10, 9 , 8; 0.644, 0.812, 1.024, 1.291, 1.628, 1.828, 2.053, 2.305, 2.588, 2.906, 3.264; NaN, NaN, NaN,  NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]; % Litz Wire % AWG;CD;WD

AWG(3,:) = AWG(2,:).*1.2;

AWG_OD = AWG(3,:);

wd(wd < AWG_OD(:,1)) = inf;

%wd_rounded = arrayfun(@(x) max(AWG_OD(AWG_OD <= x)), wd);

% wd_rounded = wd;
% cd_rounded = wd;

wd_rounded = zeros(size(wd));
cd_rounded = zeros(size(wd));
indices = zeros(size(wd));  % Initialize with NaN to ensure proper size



% Loop through each value in wd
for i = 1:numel(wd)
    if isnan(wd(i))
        continue;  % Skip invalid wd values
    end

    % Find AWG sizes that are <= wd(i)
    valid_indices = find(AWG_OD <= wd(i));
    if ~isempty(valid_indices)
        index = valid_indices(end);  % Use largest valid one
        indices(i) = index;
        
        % Assign rounded values
        wd_rounded(i) = AWG(3,index);
        cd_rounded(i) = AWG(2,index);
    end
end

% Set wd_rounded values greater than wd_max to NaN
wd_rounded(wd_rounded > wd_max) = NaN;
cd_rounded(isnan(wd_rounded)) = NaN;  % Ensure cd_rounded matches


% Stator thickness
t_min = 2 .* wd_rounded;

% Compute L (total conductor length)
L = 2 .* (ro - ri) + theta_p .* (ro + ri);

cd_rounded = cd_rounded ./ 1000; % in m

% Compute del_L (Inner length variation)
del_L1 = (ri_prime .* (theta_p - Theta_A)) ./ cos(theta_c) + ri_prime .* Theta_A - ri .* theta_p;

wd_rounded_matrix = wd_rounded;
cd_rounded_matrix = cd_rounded;



% Outer Chamfer
theta_co_values = (3:1:44);
theta_co_size = size(theta_co_values);

imax = theta_co_size(2);
jmax = numel(wd_rounded_matrix);

% Allocating Matrices

theta_ao = zeros(imax,jmax);
ro_prime = zeros(imax,jmax);
L_prime_tot = zeros(imax,jmax);
R_tot_phase = zeros(imax,jmax);
By_avg = zeros(imax,jmax);
k_T = zeros(imax,jmax);
vol = zeros(imax,jmax);
P_Loss_Eddy = zeros(imax,jmax);
I_req = zeros(imax,jmax);
P_Loss_R = zeros(imax,jmax);
P_Loss = zeros(imax,jmax);
Efficiency = zeros(imax,jmax);
P_in = zeros(imax,jmax);
V_req = zeros(imax,jmax);
lambda = zeros(imax,jmax);
winding_eff = zeros(imax,jmax);
L_prime = zeros(imax,jmax);
del_L2 = zeros(imax,jmax);
g_record = zeros(imax,jmax);
phi_f = zeros(imax,jmax);
wd_check = zeros(imax,jmax);

for j = 1:jmax
for i = 1:imax
theta_co = theta_co_values(i) .* pi ./ 180; % 
wd_rounded = wd_rounded_matrix(j);
cd_rounded =  cd_rounded_matrix(j);



taoA = 1 ./ (n_phi .* gamma);
taoB = wd_rounded .* cos(theta_co) .* n_phi .* gamma ./ ro;

%if theta_co <= 45
taoC = (1 + cos(2 .* theta_co) + tan(theta_co ./ 2) .* sin(2 .* theta_co));

%else
%    taoC = (1 - cos(2 .* theta_co) + tan(theta_co ./ 2) .* sin(2 .* theta_co));
%end
% 
% taoA = 1 - (wd_rounded * n_phi * gamma * cos(theta_co))./(2 * ro * theta_p * sin(2*theta_co));
% taoB = 2 / tan(theta_co);

% theta_ao(i,j) = theta_p - taoB*taoA;

% theta_ao(i,j) = taoA .* (theta_p - wd_rounded .* cos(theta_co) .* n_phi .* gamma ./ ro ./ sin(2 .* theta_co)); % good one
% theta_ao(i,j) = taoA .* (theta_p - taoB .* (taoC));

% Ao = ro;
% Bo = theta_p;
% Co = n_phi * gamma;
% Do = cos(theta_co);
% Eo = ci;
% To = tan(theta_co / 2);
% Wo = wd_rounded;

%Fo = Wo + (Wo * To / 2);
% 
% numerator = Fo - (Ao * Bo) / (2 * Do * Co);
% denominator = -Ao / (2 * Do) + (Ao * Eo) / (Wo + Eo);
A = ro;
B = theta_p;
C = n_phi * gamma;
D = cos(theta_co);
S = sin(2 * theta_co);
E = ci;
W = wd_rounded;

% ==== CLOSED-FORM CALCULATION ====
numerator = W - (A * B * S) / (D * C);
denominator = - (A * S) / D + ((A - W) * E) / (W + E);

theta_ao(i,j) = numerator/ denominator;


% theta_ao_comp = 180./pi.*(taoA .* (theta_p - wd_rounded .* cos(theta_co) .* n_phi .* gamma ./ ro ./ sin(2 .* theta_co)));



ro_prime(i,j) = min((ro .* (1 - 0.5 .* (theta_p - theta_ao(i,j)) .* tan(theta_co))) - (wd_rounded./2 .* tan(theta_co)) , ro - wd_rounded);
%ro_prime(i,j) = (ro .* (1 - 0.5 .* (theta_p - theta_ao(i,j))));

%ro_prime(i,j) = wd * n_phi * gamma * cos(theta_co)./ (theta_p * sin(2*theta_co));

if ro * theta_ao(i,j) < wd_rounded
    wd_rounded = NaN;
    cd_rounded = NaN;
end

wd_check(i,j) = (ro_prime(i,j) * theta_p * sin(2*theta_co))./(2 * n_phi * gamma * cos(theta_co));

if wd_check(i,j) < wd_rounded
theta_co = NaN;
end

del2A = (theta_p - theta_ao(i,j)) ./ cos(theta_co);
del2B = (theta_p - theta_ao(i,j)) .* tan(theta_co);

del_L2(i,j) = ro .* (del2A + theta_ao(i,j) - del2B - theta_p);

L_prime(i,j) = L + del_L1(j) + del_L2(i,j);

L_tot = L .* gamma .* n_p ./ 1000; % in m


L_prime_tot(i,j) = L_prime(i,j) .* gamma .* n_p ./ 1000; % in m

% Compute istive Losses

rho_copper = 1.724e-8; % in Ohm*m
dens_copper = 8960;

R_tot_phase(i,j) = rho_copper .* 4 .* L_prime_tot(i,j) ./ (pi .* cd_rounded.^2); 

% Find the minimum resistance and corresponding wd
[R_min, idx] = min(R_tot_phase(:));

% Find the Average Magnetic Flux Density as a Function of wd

% Winding factor

k_dist = sin(pi/2/n_phi) ./ gamma ./ sin(pi/2/n_phi/gamma);

Qc = gamma .* n_phi .* n_p;

theta_m = pi .* n_p ./ Qc;

k_pitch = 1;

k_w = k_dist .* k_pitch;

% Variables
n_values = (1:2:15); % harmonic number (i acc dk what to put this up to)
r_avg = (0.5 .* (ri + ro_prime(i,j))) ./ 1000; % average radius
k_d = ri ./ ro_prime(i,j); % radius ratio
tm = mag./1000; % magnet thickness
tau_p = theta_p; % pole pitch in radians
d_p = r_avg .* tau_p; % pole pitch in m
tau_m = tau_p - (1 * pi/180); % magnet pitch
mu_o = 1.25663e-6; % permenability of air
mu_rec = 1.05; % recoil multiplier
wd_rounded_in_m = wd_rounded./1000; % changing wire diameter to m
core_ironless = ci./1000; % ironless core to wrap wires around
t = 2 .* wd_rounded_in_m + core_ironless; % stator thickness
g = (t + 4./1000 + 2 .* tm);
d_m = r_avg .* tau_m;
B_rem = 1.3; % remenant magnetic field of NdFeB magnets (probably?)
air_gap = t + 4./1000;
d = 3000;

B_n = zeros(numel(n_values,1));
By_n = zeros(d,numel(n_values));

g_record(i,j) = g;
    for k = 1:numel(n_values)
        
        n = n_values(k);
        u_n = pi .* n ./ d_p;

        theta = linspace(-tau_p/2, tau_p/2, d);

        x = (r_avg .* theta)';

        J_n = 4 ./ d_p .* B_rem ./ mu_o ./ mu_rec .* sin(u_n .* d_m ./ 2);

        B_n(k) = J_n .* mu_o ./ u_n .* sinh(u_n .* tm) ./ sinh( u_n .* g ./ 2);
        
        B_yn(:,k) = B_n(k) .* cos(u_n .* x);

    end
% Find Average Flux Density (Root mean square)
    
    By = sum(B_yn, 2);

   % plot(By);
    
    dx = r_avg .* tau_p ./ d;



    By_avg(i,j) = sum(By .* dx)./ (r_avg .* tau_p);


% Find Torque Constant
    
    
   % phi_f(i,j) = By_avg(i,j) .* pi ./ 4 ./ n_p .* (2 .* ro_prime(i,j)./1000).^2 .* (1 - k_d.^2);
    phi_f(i,j) = By_avg(i,j) .* pi / n_p .* ((ro_prime(i,j)./1000).^2 - (ri./1000).^2);



    %k_T(i,j) = 1 .* n_p ./ pi .* N .* k_w .* phi_f(i,j); % Nm/A
   % k_T(i,j) = n_p ./ pi .* N .* k_w .* phi_f(i,j);
   k_T(i,j) = 1 ./ sqrt(2) .* n_p ./ 2 .* N .* k_w .* phi_f(i,j);


% Eddy Current Losses
% for v = 1:2:7
%     f_v = f .* v ;
%     w_v = 2 .* pi .* f_v;
%     k_v = sqrt(w_v .* mu_o .* mu_rec ./ 2 ./ rho_copper);
%     beta_v = v .* pi ./ d_p;
% 
%     B_z1 = max(By(:));
% 
% 
% end
    
    %vol(i,j) = N .* pi .* cd_rounded.^2 ./ 4 .* (ro_prime(i,j) - ri) ./1000;
   
    
    B_z1 = max(By(:));
    nu_d = 1;
    num_strands = 340;
    strand_diameter = cd_rounded./sqrt(num_strands);
    strand_radius = strand_diameter ./ 2;
    w_v = 2 .* pi .* f; % electrical frequency rad/s
    lambda(i,j) = (ro_prime(i,j) - ri)./1000;

    
   % P_Loss_Eddy(i,j) = pi.^2 ./ 4 ./ rho_copper .* f^2 .* strand_diameter.^2 .* N .* n_phi .* vol(i,j) .* B_z1.^2 .* nu_d.^2 ;
    P_Loss_Eddy(i,j) = pi .* B_z1.^2 .* w_v.^2 .* strand_radius.^4 .* num_strands.* N .* n_phi .* lambda(i,j) ./ 8 ./ rho_copper;

 winding_eff(i,j) = L_prime_tot(i,j) ./ lambda(i,j) ./ N;

    %P_Loss_Eddy(i,j) = 2.1;
% Find Power Loss at a Given Torque 


v_req_kph = v_car_kph(velocity_index);

v_req = v_req_kph / 3.6; % m/s

w_req = v_req ./ (0.279); % rad/s

T_req = T_car(velocity_index); % Nm

I_req(i,j) = T_req ./ k_T(i,j);

P_Loss_R(i,j) = 3 .* I_req(i,j).^2 .* R_tot_phase(i,j);

P_Loss_FW = 2.1;

P_Loss_Whl = 18.8/111 .* v_req_kph;

P_Loss(i,j) = P_Loss_R(i,j) + P_Loss_FW + P_Loss_Eddy(i,j) + P_Loss_Whl;

P_out = w_req .* T_req;


Efficiency(i,j) = (P_out ./ (P_out + P_Loss(i,j)));


P_in(i,j) = P_out ./ Efficiency(i,j);


V_req(i,j) = P_in(i,j) ./ 3 ./ I_req(i,j);

 % if V_req(i,j) > 200
 %        Efficiency(i, j) = NaN;
 % end
end
end
Efficiency(Efficiency > 1) = 0; 
% Find the maximum efficiency and its index
[Maximium_Efficiency, max_eff_idx] = max(Efficiency(:), [], 'omitnan');
[max_eff_i, max_eff_j] = ind2sub(size(Efficiency), max_eff_idx);


% Get the corresponding input voltage at maximum efficiency
V_in_max_eff = V_req(max_eff_idx);

% Display results
fprintf('Maximum Efficiency: %.2f%%\n', Maximium_Efficiency * 100);
fprintf('Input Voltage at Maximum Efficiency: %.2f V\n', V_in_max_eff);

Maximum_Efficiencies(z, gamma_index, c,m,ci) = Maximium_Efficiency;

%max(V_req(:))

P_Loss_at_max_eff = P_Loss(max_eff_idx)

P_in_at_max_eff = P_in(max_eff_idx)

theta_co_at_max_eff = theta_co_values(max_eff_i)

theta_ao_at_max_eff = 180 ./ pi .* theta_ao(max_eff_i, max_eff_j)

wd_at_max_eff = wd_rounded_matrix(max_eff_j)

ro_prime_at_max_eff = ro_prime(max_eff_idx)

theta_A = Theta_A;

theta_a_at_max_eff =  180 ./ pi .* theta_A(max_eff_j)

theta_c_at_max_eff =  180 ./ pi .* Theta_C(max_eff_j)

I_req_at_max_eff = I_req(max_eff_idx)

R_tot_phase_at_max_eff = R_tot_phase(max_eff_idx)

k_T_at_max_eff = k_T(max_eff_idx)

ri_prime_at_max_eff = ri_prime(max_eff_j)

P_Loss_Eddy_at_max_eff = P_Loss_Eddy(max_eff_idx)

P_Loss_R_at_max_eff = P_Loss_R(max_eff_idx)

g_at_max_eff = 1000 .* g_record(max_eff_idx)

L_prime_tot_max_eff = L_prime_tot(max_eff_idx)
end
end
end
end
end