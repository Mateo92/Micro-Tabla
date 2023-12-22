% Some material properties %%%%%%%

rho_au=19.3*1e3 ; %gold density kg/m3
E_au=79*1e9 ; %gold Young modulus Pa=kg/(ms^2)
nu_au=0.4 ; %gold Poisson's ratio

rho_al=2.7*1e3 ;  %Aluminium density kg/m3
E_al=70*1e9 ; %Aluminium Young modulus Pa=kg/(ms^2)
nu_al=0.35 ; %Aluminium Poisson's ratio

rho_cytop=2.03*1e3 ; %CYTOP density Kg/m3
E_cytop=7.9*1e9 ; %CYTOP Young modulus Pa=kg/(ms^2)
nu_cytop=0.42 ; %CYTOP Poisson's ratio

%%% External radius features
Di_ext=150*1e-6;   %External Diameter in [m]
H_ext=15*1e-6;     %External thickness in [m]
E_ext=E_cytop;     %External Young modulus
rho_ext=rho_cytop; %External density
nu_ext=nu_cytop;  %External Poisson's ratio

% %Homogeneous cytop plate first eigenfrequency
D=E_cytop*H_ext^3/(12*(1-nu_cytop^2));
kappa=sqrt(D/(rho_cytop*H_ext))/(Di_ext/2)^2
f_res=kappa*10.22*1e-6/(2*pi)
f_2=kappa*21.26*1e-6/(2*pi)
f_3=kappa*34.88*1e-6/(2*pi)
f_11=kappa*39.77*1e-6/(2*pi)
f_40=kappa*51.04*1e-6/(2*pi)
f_22=kappa*60.82*1e-6/(2*pi)


% f_2/f_res
%%
%%% Internal radious features
Di_int=Di_ext*0.5;  %Internal Diameter in [m]
H_int=H_ext*1.54;      %Internal thickness in [m]

% Di_int=0;         %Homogeneous case
% H_int=0;          %Homogeneous case

E_int=E_au;         %Internal Young modulus
rho_int=rho_au;     %Internal density
nu_int=nu_au;       %Internal Poisson's ratio

%% Initial conditions parameters

%%%  Initial position

% PE=0; %Initial potencial enegy [Joules]    
% 
% u0=Heigth_PE(PE,Di_int,nu_int,E_int,H_int,Di_ext,nu_ext,E_ext,H_ext); %initial position height [mm]

 u0=0.0;       %initial position height [m]

%%% Initial velocity

ctr = [0 0]*Di_ext/2; % center location in cartecian coordinates [x,y]
wid = 0.25*Di_ext;  % width of excitation for initial velocity
KE=1e-6;         %   Maximum initial kinetic energy [Joules]

v0=init_vel_KE(KE,H_ext,rho_ext,Di_int,H_int,rho_int,wid); %maximum initial kinetic energy
v0=130;
% v0=5;         %maximum initial velocity [m/s]

%%% Read out points

 rp_mat = [[0 0];[0.7 0.5*pi];[0.9 1*pi];];   % position of readout([0-1,0-2*pi])
%rp_mat = [[0 0];];   % position of readout([0-1,0-2*pi])

% [x,y]=pol2cart(0.5, 0.5*pi);
% rp_mat_cart = [[0.3 0.2];[0.4 0.3];[0.45 0.25];...
%               [0.5 0.5]; [0.6 0.15]; [0.7 0.8]]-0.5; % position of readout([x,y])
% 
% rp_mat=zeros(size(rp_mat_cart));
% for i=1:size(rp_mat,1)
%     rp=rp_mat_cart(i,:);
%     rp_x=rp(1);
%     rp_y=rp(2);
%     [rp_th,rp_r]=cart2pol(rp_x,rp_y);
%     while rp_th<0
%     rp_th=rp_th+2*pi;
%     end
%     rp_mat(i,:)=[rp_r rp_th];
% end
% rp_mat(7,:)=[0.5, 0.5*pi]; % position of readout([0-1,0-2pi])
% rp_mat
%%

%%% Stifnes parameter's parameters

logistic_fun=0;  %Decide is using logistic growth 1=true, 0= false
k_sigmoid=55;   %Logistic growth parameter (only affects the simulation if logistic_fun==0)

%%% Spacial mesh parameters

Nr=30; %Number of radial points -1 := Nr-1 (because 0 is added later)

Nt=30; %Number of angular points := Nt

%%% Time domain and loss parameter

sigma_time=1e-6;  %Time scaling factor

% E_scale_param=1e-12;

v0=v0*sigma_time; %Scaled Velocity

E_scale_param=sigma_time^2; % Scale parameter for the Young's modulus

kstability_constant=1;

TF=10;             % simulation time [sqrt(E_scale_param)s]
T60 = TF;       % 60 dB decay time [sqrt(E_scale_param)s]
sig0 = 6*log(10)/T60;           % loss parameter

% sig0=0;            % For free vibrations

%%% Plot and animation

ploting=1;  %Plots are shown if ploting==1
anim=1;     %Animation of the simulation is shown if anim==1
n_anim=5000;   %Number of time steps between each frame in the animation
plot_fft=1;  %Plot Fourier transform of the output signals  
limz=1e-6;
time_series=0;
%%%Write metadata

% T = table(PE,Di_ext,H_ext,E_ext,rho_ext,nu_ext,Di_int,H_int,E_int,rho_int,nu_int,...
%     ctr,wid,u0,v0,sig0,k_sigmoid,logistic_fun,Nr,Nt,TF);
% folder=strcat('CSV/Single_simulations/',datestr(datetime('now')),'/');
% mkdir(folder)
% metadatafolder=strcat(folder,'/metadata/');
% mkdir(metadatafolder)
% writetable(T,strcat(metadatafolder,'metadata.csv'))
% writematrix(rp_mat,strcat(metadatafolder,'rp_mat.csv'))

% %%%%% Run simulation

 [out_polar,SR,hr,ht,k_polar]=...
      Polar_plate_2lapl(Di_ext,H_ext,E_ext,rho_ext,nu_ext,Di_int,H_int,E_int,rho_int,nu_int,E_scale_param,...
kstability_constant,ctr,wid,u0,v0,sig0,k_sigmoid,logistic_fun,Nr,Nt,TF,rp_mat,ploting,anim,n_anim,limz,plot_fft,time_series);

%Write outputs
% writematrix(out,strcat(folder,'out.csv'));
% writematrix(SR,strcat(folder,'SR.csv'));
% writematrix(hr,strcat(folder,'hr.csv'));
% writematrix(ht,strcat(folder,'ht.csv'));
% writematrix(u0,strcat(folder,'u0.csv'));


