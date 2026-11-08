%% ========================================================================
%  FILE: main script to run the simulink model and plots
%% Description of the script %%
% 1. The first part defines:
%                       - geometry and inertia of the spacecraft
%                       - orbital parameters and constants
%                       - initial conditions of angular velocity and DCM
%                       - sensors and actuators characteristics
%                       - solution of Riccati equation for control
% 2. Second part:
%                - defines simulation options
%                - runs the model
%                - extracts the results
% 3. Third part computes and plots the eclipse
% 4. Fourth part shows the plots
%
% AUTHORS:
%   Simone Frescura
%   Nicola Pifano
%   Luca Puzzi
%   Lorenzo Veronese
clc; clear; close all
rng(1); % fixed random seed: sensor errors (and therefore results) are reproducible
%%
%% Principal moments of inertia of the spacecraft
% geometry of the spacecraft
a_SC = 0.300;  % [ m ]
b_SC = 0.200;  % [ m ]
c_SC = 0.100;  % [ m ]
mass = 7.23; % [kg]
% Inertia moments
Ix= (1/12) * mass * (b_SC^2 + c_SC^2);   % [ kg m^2 ]
Iy = (1/12) * mass * (a_SC^2 + c_SC^2);   % [ kg m^2 ]
Iz = (1/12) * mass * (a_SC^2 + b_SC^2);   % [ kg m^2 ]
% Inertia matrix
J_B=diag([Ix,Iy,Iz]); %[kg*m^2]
%% Orbital parameter of satellite's orbit in the ECI frame
R_earth = 6371; % km
mu_earth = 398600.433; %Gravity constant of the Earth [km^3/s^2]
a=6978; %Semimajor axis of the orbit [km]
e=0; %Eccentricity of the orbit
i=97.8; %Inclination of the orbit [deg]
Omega=0; %Right ascension of the ascending node [deg]
omega=0; %Argument of pericenter [deg]
theta=0; %Initial true anomaly [deg]
% Conversion from degrees to radians
i_rad=deg2rad(i);
RAAN_rad=deg2rad(Omega);
omega_rad=deg2rad(omega);
[r0_ECI,v0_ECI] = COEtor0v0(a,e,i,Omega,omega,theta); %Initial position and velocity of the spacecraft expressed in ECI [km] [km/s]
T = 2 * pi * sqrt(a^3 / mu_earth) ; %Orbital period [ s ]
%% Position of the Earth wrt the Sun in the ECI frame
mu_sun = 132712440017.987; % [km^3/s^2]	
AU = 149597870.691;     % astronomic unit [km]
rES0_ECI = [0 AU 0];
vES0_ECI = [(-sqrt(mu_sun/norm(rES0_ECI))) 0 0 ] ;
% Initial Sun relative to Earth
rSun0_ECI = -rES0_ECI;
vSun0_ECI = -vES0_ECI;
%% Initial DCM (direction cosine matrix) from ECI to Body
A_BN_0=eye(3);
%% Initial angular velocity
wx=0.5; %[rad/s]
wy=0.3; %[rad/s]
wz=0.4; %[rad/s]
omega0=[wx;wy;wz]; %Initial angular velocity vector [rad/s]
%% Local magnetic dipole of the spacecraft expressed in body-frame [A*m^2]
m = [0.01; 0.01; 0.01];
%% Sensors specification
seed=0;
delay=1;
D=500*10^-12; % [T/sqrt(Hz)] Noise density
f_m=10; % [Hz] max sampling frequency
Ts_m=1/f_m; % [s] magnetometer sampling time
non_ort=deg2rad(-1+2*rand); % [rad] orthogonality error
%non_ort=deg2rad(-1+2);
Res=6.5*10^-9;  %[T] resolution
fc_pink = 0.01;
k_pink=1;
a_pink=exp(-2*pi*fc_pink*Ts_m*k_pink);
%Sun sensor
acc=deg2rad(0.5*randn); %[rad] accuracy (gaussian)
f_ss=10; %[hz]
Ts_ss=1/f_ss; % [s]
%% Reaction wheels specification
Irx = 0.2*10^-3; %[kg*m^2]
Iry = 0.2*10^-3; %[kg*m^2]
Irz = 0.2*10^-3; %[kg*m^2]
Ir = diag([Irx,Iry,Irz]); %Relative inertia matrix of the rotors [kg*m^2]
inv_Ir = inv(Ir);
A = eye(3); %Matrix whose columns indicate the direction of the axis of rotation of the wheels
Ap = eye(3); %Pseudoinverse of A
sigma = (2*pi/60)/3; %Standard deviation of the rotors noise
w_max = [250; 250; 250]; %Saturation speed of the rotors [rad/s]
wr_0 = [0; 0; 0]; %Initial speed of the rotors [rad/s]
Mr_max = [0.01; 0.01; 0.01]; %Maximum relative torque applied to the wheels [N m]
%%Magnetorque
m_max=[0.65; 0.65; 0.5]; %Maximum local magnetic dipole that the magnetorquer can generate [A*m^2]
%% Riccati equation solver
Ar = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
Br = [0 0 0; 0 0 0; 0 0 0; 1/Ix 0 0; 0 1/Iy 0; 0 0 1/Iz];
S = zeros(6,3);
ximax = [deg2rad(5) deg2rad(5) deg2rad(5) deg2rad(2) deg2rad(2) deg2rad(2)];
uimax = [10^-4 10^-4 10^-4]; %Maximum desired torque
ar = 1./(ximax.^2);
br = 1./(uimax.^2);
Q = diag(ar);
R = diag(br);
[~,K,~] = icare(Ar,Br,Q,R);
%% Simulation Options
sim_options.SolverType = 'Fixed-step';      % Set the solver type to Fixed-step
sim_options.Solver = 'ode4';                % Select ode4 as solver
sim_options.FixedStep = '0.01';              % Select a time step of 0.01 s
sim_options.StartTime = '0';                % Start from 0 seconds [default]
sim_options.StopTime = 'T';               % End the simulation
%% Simulate the model
results = sim('Complete_SAD_model', sim_options);
%% Results extraction
time = results.tout;
N=length(time);
% DKE
DCM= results.DCM;
w_real=results.omega;
M_gg=results.gravity_torque;
M_B= results.Magnetic_Torque;
M_tot=results.Tot_disturbances; %add to workspace
% sensors
b_B=results.b_B;
b_B_measured=results.b_B_measured;
S_B=results.S_B;
S_B_measured=results.S_B_measured;
% attitude determination
DCM_mis= results.A_mis;
w_mis=results.w_mis;
% Actuators
Mc_ideal_mgt = results.Mc_ideal_mgt;
M_mgt=results.M_mgt;
M_ideal_rw=results.M_ideal_rw;
M_rw=results.M_rw;
Mc_tot=results.Mc_tot;
pointing_error = results.pointing_error;
%% ==============================================================
% SUN-EARTH-SC Simulation and Data Extraction
% ==============================================================
time_s    = time;      % [s]
time_days = time_s / 86400;           % [giorni]
r_eci     = results.r_ECI;        % SC relative to Earth [km]
r_ES_eci  = results.r_ES_ECI;     % Earth relative to the SUN [km]
r_sun_eci = results.r_SE_ECI;     % Sun relative to Earth [km]
r_sun_sc  = results.r_SATSUN_ECI; % Sun relative to SC [km]
r_ES_norm    = sqrt(sum(r_sun_eci.^2, 2));    % (vector EARTH to SUN  norm)[km]
r_SCsun_norm = sqrt(sum(r_sun_sc.^2, 2));     % (vector SC to SUN norm)[km]
T_orbit    = T;                               % [s]
mask_orbit1 = (time_s >= time_s(1)) & (time_s <= time_s(1) + T_orbit);
%% ==============================================================
%  Eclipse (Cilindric Shadow) + period over a single orbit
%  Use r_eci (Earth->SC) and r_sun_eci (Earth->Sun) from Simulink
% ==============================================================
% Sun direction from the Earth (versor Earth->Sun)
r_ES_norm2 = sqrt(sum(r_sun_eci.^2, 2));                 % [km]
u_sun_E    = r_sun_eci ./ r_ES_norm2;                    % [N x 3] unit vector
% Condition 1: satellite "behind" the Earth relative to the Sun
% (r_sc · u_sun < 0)
proj_sc = sum(r_eci .* u_sun_E, 2);                      % [km]
behind  = (proj_sc < 0);
% Condition 2: distance of the satellite from SUN-EARTH vector <= R_earth
% d_perp = || r_sc - (r_sc·u) u ||
r_perp  = r_eci - u_sun_E .* proj_sc;                    % [N x 3] km
d_perp  = sqrt(sum(r_perp.^2, 2));                       % [km]
inCyl   = (d_perp <= R_earth);
% Eclipse (cilindric shadow)
inEclipse = behind & inCyl;                              % [N x 1] logical
%% Eclipse duration over ONE single orbit (mask_orbit1)
t1   = time_s(mask_orbit1);
ecl1 = inEclipse(mask_orbit1);
Tecl_orbit1 = trapz(t1, double(ecl1));                   % [s]
fprintf('\n--- ECLIPSE (orbit 1) ---\n');
fprintf('Orbital Period T        = %.3f min\n', T/60);
fprintf('Time in Eclipse (Orbit) = %.3f min (%.1f %% of the orbit)\n', ...
   Tecl_orbit1/60, 100*Tecl_orbit1/T);
%% Initial/final eclipse time guess in the first orbit
e = ecl1(:);
tt = t1(:);
de = diff([0; e; 0]);         % padding to capture complete segments
iEntry = find(de == 1);       % initial segment indices (in e)
iExit  = find(de == -1) - 1;  % final segment indices
if isempty(iEntry)
fprintf('None eclipses in the first orbit.\n');
else
% If there is more than one eclipse segment (rare), take the longest one
segDur = tt(iExit) - tt(iEntry);
[~,k] = max(segDur);
tEntry = tt(iEntry(k));
tExit  = tt(iExit(k));
fprintf('Eclipse start ~ t = %.3f min after orbit start\n', (tEntry-tt(1))/60);
fprintf('Eclipse finish~ t = %.3f min after orbit start\n', (tExit -tt(1))/60);
end
%% Plot: eclipse over a single orbit
figure;
plot(time_days(mask_orbit1), double(ecl1), 'LineWidth', 1.6);
grid on;
xlabel('Time [days]');
ylabel('In Eclipse (1=yes)');
title('Eclipse over a single orbit ');
ylim([-0.1 1.1]);
%% DKE plots
% Angular velocity evolution
figure
plot(time,w_real); grid on;
xlabel('time [s]')
ylabel(' w [rad/s]')
title('Angular velocity evolution')
legend('w_x','w_y','w_z')
% Body axis evolution in inertial frame
x_body=zeros(N,3);
y_body=zeros(N,3);
z_body=zeros(N,3);
for i=1:N
  x_body(i, :) = DCM(1,:,i); % Body x-axis in inertial frame
  y_body(i, :) = DCM(2,:,i); % Body y-axis in inertial frame
  z_body(i, :) = DCM(3,:,i); % Body z-axis in inertial frame
end
figure
subplot(3,1,1)
plot(time,x_body); grid on;
xlabel('time [s]')
ylabel('x body [-]')
legend('x component','y component','z component')
title('time evolution of x axis body frame')
subplot(3,1,2)
plot(time,y_body); grid on;
xlabel('time [s]')
ylabel('y body [-]')
legend('x component','y component','z component')
title('time evolution of y axis body frame')
subplot(3,1,3)
plot(time,z_body); grid on;
xlabel('time [s]')
ylabel('z body [-]')
legend('x component','y component','z component')
title('time evolution of z axis body frame')
% DCM orthornormalization check
DCM_orth_err = zeros(N,1);
for i=1:N
  DCM_orth_err(i)=norm(DCM(:,:,i)'*DCM(:,:,i)-eye(3));
end
figure
plot(time,DCM_orth_err); grid on;
xlabel('time [s]')
ylabel('Orthonormality error [-]')
title('DCM orthonormality error')
% verify angular momentum and kinetic energy
h_norm=zeros(N,1);
Ec=zeros(N,1);
for k=1:N
   h=J_B*w_real(k,:)';
   h_norm(k)=norm(h);
   Ec(k)=0.5*w_real(k,:)*(J_B*w_real(k,:)');
end
figure
subplot(2,1,1)
plot(time, h_norm); grid on;
xlabel('time [s]')
ylabel ('h [kg*m^2/s]')
title('Angular momentum norm')
subplot(2,1,2)
plot(time,Ec); grid on;
xlabel('time [s]')
ylabel ('Ec [J]')
title("Kinetic energy")
%% torques plots
Mgg_norm=zeros(N,1);
MB_norm=zeros(N,1);
Mtot_norm=zeros(N,1);
for k=1:N
  Mgg_norm(k)=norm(M_gg(k,:));
  MB_norm(k)=norm(M_B(k,:));
  Mtot_norm(k)=norm(M_tot(k,:));
end
figure
plot(time,Mgg_norm,time,MB_norm,time,Mtot_norm); grid on;
xlabel('time [s]')
ylabel('M [Nm]')
title('Perturbation torques')
legend('Gravity gradient', 'Magnetic field','Total disturbance torque');
%% Sensors plot
%  Measured magnetic field
figure
plot(time,b_B_measured); grid on;
xlabel('time [s]')
ylabel('b measured [T]')
title('Measured magnetic field (body)')
legend('x component','y component','z component')
% Real magnetic field
figure
plot(time,b_B); grid on;
title('Real magnetic field')
xlabel('time [s]')
ylabel('b_B [T]')
legend('x component','y component','z component')
% Calculate the error between measured and actual magnetic field
err_b=zeros(N,1);
for k=1:N
 err_b(k)=norm(b_B(k,:)-b_B_measured(k,:));
end
figure
plot(time,err_b); grid on
xlabel('time [s]')
ylabel('error [T]')
title('Magnetic field measurement error')
%  Measured Sun direction
figure
plot(time,S_B_measured); grid on;
xlabel('time [s]')
ylabel('S_B measured [-]')
title('Measured Sun direction (body)')
legend('x component','y component','z component')
% Real Sun direction
figure
plot(time,S_B); grid on;
title('Real Sun direction')
xlabel('time [s]')
ylabel('S_B [-]')
legend('x component','y component','z component')
% Calculate the error between measured and actual sun direction
err_s=zeros(N,1);
for k=1:N
  err_s(k)=rad2deg(acos(dot(S_B(k,:),S_B_measured(k,:))/(norm(S_B(k,:)*norm(S_B_measured(k,:))))));
end
figure
plot(time, err_s); grid on;
xlabel('time [s]')
ylabel('Error [deg]')
title('Sun direction error')
figure
plot(time, Mc_ideal_mgt); grid on;
xlabel('time [s]')
ylabel('Torque [N m]')
legend('x component','y component','z component')
title('Torque sent from attitude controller to magnetorquer')
figure
plot(time, M_mgt); grid on;
xlabel('time [s]')
ylabel('Torque [N m]')
legend('x component','y component','z component')
title('Real torque generated by the magnetic torquer')
figure
plot(time, M_ideal_rw); grid on;
xlabel('time [s]')
ylabel('Torque [N m]')
legend('x component','y component','z component')
title('Torque sent from attitude controller to reaction wheels')
figure
plot(time, M_rw); grid on;
xlabel('time [s]')
ylabel('Torque [N m]')
legend('x component','y component','z component')
title('Torque generated by the reaction wheels')
figure
plot(time, Mc_tot); grid on;
xlabel('time [s]')
ylabel('Torque [N m]')
legend('x component','y component','z component')
title('Total torque generated by the actuators')
%% Attitude determination plots
% algebraic method
err_DCM=zeros(N,1);
for k=1:N
 err_DCM(k) = rad2deg(acos((trace(DCM(:,:,k)'*DCM_mis(:,:,k))-1)/2));
end
figure
plot(time,err_DCM); grid on;
xlabel('time [s]')
ylabel('A error [deg]')
title('Error between A measured and real DCM')
% angular velocity error
err_w=zeros(N,1);
for k=1:N
  err_w(k)=rad2deg(norm(w_real(k,:)-w_mis(k,:)));
end
figure
plot(time,err_w); grid on;
xlabel('time [s]')
ylabel('w err [deg/s] ')
title('Angular velocity error between real w and measured w')
figure
plot(time,pointing_error); grid on;
xlabel('time [s]')
ylabel('Pointing error [deg] ')
title('Pointing error between x axis of the BF and the satellite-Sun direction')


%% Local COEtor0v0 function
function [r0_COE,v0_COE] = COEtor0v0(a,ecc,i,raan,w,v)
mu=398600.4355067;
i=deg2rad(i);
raan=deg2rad(raan);
w=deg2rad(w);
v=deg2rad(v);
R3_RAAN = [cos(raan), sin(-raan), 0; sin(raan), cos(raan), 0; 0, 0, 1];
R1_inc = [1, 0, 0; 0, cos(i), -sin(i); 0, sin(i), cos(i)];
R3_omega = [cos(w), sin(-w), 0; sin(w), cos(w), 0; 0, 0, 1];
R = R3_RAAN*R1_inc*R3_omega;
r = a*(1-ecc^2)/(1+ecc*cos(v));
h = sqrt(a*mu*(1-ecc^2));
r0_COE = R*r*[cos(v); sin(v); 0];
v0_COE = R*[-mu*sin(v)/h; mu*(ecc+cos(v))/h; 0];
end
