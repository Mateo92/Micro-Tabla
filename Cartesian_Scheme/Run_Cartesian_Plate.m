% Some material properties %%%%%%%
rho_au = 19.3 * 1e3;  % gold density kg/m3
E_au = 79 * 1e9;  % gold Young modulus
nu_au = 0.4;  % gold Poisson's ratio

rho_al = 2.7 * 1e3;  % Aluminium density kg/m3
E_al = 70 * 1e9;  % Aluminium Young modulus 
nu_al = 0.35;  % Aluminium Poisson's ratio

rho_cytop = 2.03 * 1e3;  % CYTOP density Kg/m3
E_cytop = 7.9 * 1e9;  % CYTOP Young modulus Pa=kg/(ms^2)
nu_cytop = 0.42;  % CYTOP Poisson's ratio

%%% External radius features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Di_ext = 150 * 1e-6;
H_ext = 15 * 1e-6;
E_ext = E_cytop;
rho_ext = rho_cytop;
nu_ext = nu_cytop;

%%% Internal radious features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Di_int = Di_ext * 0.3;
H_int = H_ext * 0.1;
E_int = E_au;
rho_int = rho_au;
nu_int = nu_au;

%%% Homogeneous cytop plate first eigenfrequencies %%%%%%%%%%%%%%%%%%%%%%%%
D = E_cytop * H_ext^3 / (12 * (1 - nu_cytop^2));
kappa = sqrt(D / (rho_cytop * H_ext)) / (Di_ext / 2)^2;
f_fund = kappa * 10.22 * 1e-6 / (2 * pi);
f_11 = kappa * 21.26 * 1e-6 / (2 * pi);
f_12 = kappa * 34.88 * 1e-6 / (2 * pi);
f_20 = kappa * 39.77 * 1e-6 / (2 * pi);
f_13 = kappa * 51.04 * 1e-6 / (2 * pi);
f_21 = kappa * 60.82 * 1e-6 / (2 * pi);

%%% Channels measures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ch_len = 14 * 1e-5;
ch_wid = 5e-5;

ch_wid = 0;
ch_len = 0;

%%% Initial conditions parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Initial position %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PE = 0.001;  % Initial potential energy [Joules]

% u0 = Heigth_PE(PE,Di_int,nu_int,E_int,H_int,Di_ext,nu_ext,E_ext,H_ext); % %initial position height [m]

u0 = 0;  
% u0 = 4e-6;

%%%% Initial velocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ctr = [0.5 0.75];  % center location 
wid = 0.25;  % width of excitation
KE = 1e-6;  % Maximum initial kinetic energy [Joules]

v0 = init_vel_KE(KE, H_ext, rho_ext, Di_int, H_int, rho_int, wid * Di_ext);  % maximum initial kinetic energy

%%% Logistic loading function parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_sigmoid = 2.0;  % Logistic growth parameter

%%% Spatial mesh parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = 128;  % number of spatial subdivisions in the x and y axes

%%% Time domain and loss parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma_time = 1e-6;  % Time scaling factor

v0 = v0 * sigma_time;  % Scaled Velocity

TF = 1;  % simulation time (10^-4s)
T60 = TF * 2;  % 60 dB decay time (10^-4s)
sig0 = 6 * log(10) / T60;  % loss parameter
% sig0 = 0;  % Undamped

%%% Readout position %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rp_mat = [[0.5 0.5]; [0.9, 0.5]; [0.85, 0.5]];  % position of readout([0-1,0-1])

%%% Plot and animation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ploting = 1;
anim = 1;
n_anim = 500;
plot_fft = 1;
limz = 1e-6;

%%%% Write metadata %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% T = table(Di_ext,H_ext,E_ext,rho_ext,nu_ext,Di_int,H_int,E_int,rho_int,nu_int,ch_len,ch_wid,...
%     ctr,wid,KE,PE,u0,v0,sig0,k_sigmoid,Nx,TF);
% folder=strcat('Experiments/Single_experiment/',datestr(datetime('now')),'/');
% mkdir(folder)
% metadatafolder=strcat(folder,'/metadata/');
% mkdir(metadatafolder)
% writetable(T,strcat(metadatafolder,'metadata.csv'))
% writematrix(rp_mat,strcat(metadatafolder,'rp_mat.csv'))

%%% Run simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[out_cart, SR, h, ss, k_cart] = ...
    Plate_cartesian(Di_ext, H_ext, E_ext, rho_ext, nu_ext, Di_int, H_int, E_int, rho_int, nu_int, ...
    ch_len, ch_wid, ctr, wid, u0, v0, sig0, k_sigmoid, Nx, TF, sigma_time, rp_mat, ploting, anim, n_anim, limz, plot_fft);

%%% Write outputs:

% writematrix(out_cart,strcat(folder,'out.csv'));
% writematrix(SR,strcat(folder,'SR.csv'));
% writematrix(h,strcat(folder,'h.csv'));
% writematrix(k_cart,strcat(folder,'k.csv'));