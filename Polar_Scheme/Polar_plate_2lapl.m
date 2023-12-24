function [out,SR,hr,ht,k]=...
Polar_plate_2lapl(Di_ext,H_ext,E_ext,rho_ext,nu_ext,Di_int,H_int,E_int,...
rho_int,nu_int,sigma_time,k_stability_constant,ctr,wid,u0,v0,sig0,...
k_sigmoid,logistic_fun,Nr,Nt,TF,rp_mat,ploting,anim,n_anim,limz,plot_fft)
%%%%%%%%%%%%%%%
% Simulate a non-homogeneous plate using finite difference method
%
% Inputs:
% Di_ext, Di_int: External and internal diameters of the plate.
% H_ext, H_int: External and internal heights of the plate.
% E_ext, E_int: Young's moduli of the external and internal materials.
% rho_ext, rho_int: Densities of the external and internal materials.
% nu_ext, nu_int: Poisson's ratios of the external and internal materials.
% sigma_time: %Time scaling factor.
% kstability_constant: Constant used for stability calculation.
% ctr: Center coordinates for a velocity distribution calculation.
% wid: Width parameter for the velocity distribution.
% u0: Initial displacement.
% v0: Initial velocity height.
% sig0: Parameter for the velocity distribution function.
% k_sigmoid: Parameter for a sigmoid function.
% logistic_fun: Flag for selecting a logistic function.(1=true, 0= false)
% Nr: Number of radial steps for discretization.
% Nt: Number of angular steps for discretization.
% TF: Final simulation time.
% rp_mat: Matrix specifying readout positions.([0-1,0-2*pi])
% ploting: Flag for enabling plot generation. (1=true, 0= false)
% anim: Flag for enabling animation generation. (1=true, 0= false)
% n_anim: Parameter for animation frequency.
% limz: Parameter for plot z-axis limits.
% plot_fft: Flag for enabling Fast Fourier Transform (FFT) plots. (1=true, 0= false)
%
% Outputs:
% out: Simulation output for readout positions over time.
% SR: Sample Rate calculated based on stability condition.
% hr: Radial step size.
% ht: Angular step size.
% k: Time step for the simulation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic % Start measuring execution time

E_scale_param=sigma_time^2; % Scale parameter for the Young's modulus

% Calculate time scale for display
time_scale = strcat('10^{', num2str(log10(sqrt(E_scale_param))), '}s');

% Scale Young's modulus
E_ext_scale = E_ext * E_scale_param;
E_int_scale = E_int * E_scale_param;

% Compute Stiffness parameters for exterior and interior
D_ext = (E_ext_scale) * H_ext^3 / (12 * (1 - nu_ext^2));
K_ext = sqrt(D_ext / (rho_ext * H_ext * (Di_ext / 2)^4));
HT = H_int + H_ext;
rho_l = HT / (H_ext / rho_ext + H_int / rho_int);
E_l = HT / (H_ext / E_ext_scale + H_int / E_int_scale);
nu_l = HT / (H_ext / nu_ext + H_int / nu_int);
D_int = E_l * HT^3 / (12 * (1 - nu_l^2));
K_int = sqrt(D_int / (rho_l * HT * (Di_ext / 2)^4));

K_max = max([K_int, K_ext]);
K_min = min([K_int, K_ext]);
R_int = Di_int / Di_ext;

% Radial coordinates step size and count
hr = 1 / Nr;
Nr = floor(1 / hr);

% Angular coordinates step size and count
ht = 2 * pi / Nt;
Nt = round(2 * pi / ht);

% Stability condition - calculate time step 'k'
k = k_stability_constant * (hr^2 / (2 * K_max)) * (1 / (1 + 1 / (ht^2)));

% Calculate Sample Rate
SR =floor(1 / k);

% Number of Time Steps
NF = floor(SR * TF);

% Polar and Cartesian coordinates axes
[R, T] = meshgrid([0:Nr] * hr * Di_ext / 2, [0:Nt - 1] * ht);
[X, Y] = pol2cart(T, R);

% Compute readout parameters index for convenience
rp_mat_index = zeros(size(rp_mat));
for i = 1:size(rp_mat, 1)
    rp = rp_mat(i, :);
    rp_r = rp(1);
    rp_t = rp(2);
    rpr_index = floor(rp_r / hr) + 1;
    rpt_index = floor(rp_t / ht) + 1;
    rp_mat_index(i, :) = [rpr_index rpt_index];
end

% Compute rising cosine velocity distribution
%%%% Cartesian coordinates

dist = sqrt((X-ctr(1)).^2 +(Y-ctr(2)).^2);
ind = sign(max(-dist+wid/2,0)); 
rc = 0.5*ind'.*(1+cos(2*pi*dist'/wid));

%%%% Polar coordinates
[t,r,rc_pol]=cart2pol(X,Y,rc);

% Stifness parameters matrix

K_mat=zeros(Nr+1,Nt); 
for i=1:Nr+1
    l=i-1;
    r=hr*l;
    step_func=K_int;
    if r>R_int
        step_func=K_ext;
    end
    for j=1:Nt
        if logistic_fun==1
        logis=K_int-(K_int-K_ext)*(1/(1+exp(-k_sigmoid*(r-R_int))));
        K_mat(i,j)=logis;
        else
        K_mat(i,j)=step_func;
        end

    end
end
[x,y,K_mat_cart]=pol2cart(T,R,K_mat);


% Fourth order polinomial Initial position
qu=zeros(Nr+1,Nt);
for i=1:Nr+1
    l=i-1;
    r=hr*l;
    for j=1:Nt
       qu(i,j)=(1-(r)^2)^2;
    end
end
[x,y,qu_cart]=pol2cart(T,R,qu);

% Plot initial conditions
if ploting==1
figure(1)
tiledlayout(1,5)
nexttile
mesh(X',Y',qu_cart*u0)
title('Initial position','Fontsize',20,'Interpreter','latex')
xlabel('$x[m]$','Fontsize',20,'Interpreter','latex')
ylabel('$y[m]$','Fontsize',20,'Interpreter','latex')
zlabel('$u(x,y,t=0)[m]$','Fontsize',20,'Interpreter','latex')
nexttile
mesh(X',Y',rc*v0)
title('Initial velocity','Fontsize',20,'Interpreter','latex')
xlabel('$x[m]$','Fontsize',20,'Interpreter','latex')
ylabel('$y[m]$','Fontsize',20,'Interpreter','latex')
zlabel('$u_t(x,y,t=0)[m\sigma^{-1}]$','Fontsize',20,'Interpreter','latex')
nexttile
mesh(R',T',rc_pol*v0)
tit=['Initial velocity' newline '(Polar coord.)'];
title(tit,'Fontsize',20,'Interpreter','latex')
xlabel('$x[m]$','Fontsize',20,'Interpreter','latex')
ylabel('$y[m]$','Fontsize',20,'Interpreter','latex')
zlabel('$u_t(x,y,t=0)[m\sigma^{-1}]$','Fontsize',20,'Interpreter','latex')
nexttile
mesh(X',Y',K_mat_cart)
tit=['Stiffness distribution' newline '(Cartecian coord.)'];
title(tit,'Fontsize',20,'Interpreter','latex')
xlabel('$x[m]$','Fontsize',20,'Interpreter','latex')
ylabel('$y[m]$','Fontsize',20,'Interpreter','latex')
zlabel('$\kappa(x,y)[\sigma^{-1}]$','Fontsize',20,'Interpreter','latex')
nexttile
mesh(R',T',K_mat)
tit=['Stiffness distribution' newline '(Polar coord.)'];
title(tit,'Fontsize',20,'Interpreter','latex')
xlabel('$r[m]$','Fontsize',20,'Interpreter','latex')
ylabel('$\theta$','Fontsize',20,'Interpreter','latex')
zlabel('$\kappa(r,\theta)[\sigma^{-1}]$','Fontsize',20,'Interpreter','latex')

x0=10;
y0=500;
width=1400;
height=350;
set(gcf,'position',[x0,y0,width,height])
end

% Initialize grid functions
u_nm1 = u0*qu;
u_n = u0*qu+k*v0*rc_pol;

for m=0:Nt-1
    M=mod(m,Nt)+1;
    u_nm1(1,M)=u_nm1(1,1);
    u_n(1,M)=u_n(1,1);
end

u_np1 = zeros(Nr+1,Nt); 

% Initialise output

out = zeros(NF,length(rp_mat));
for i=1:size(rp_mat,1)
   out(1,i) = u_nm1(rp_mat_index(i,1),rp_mat_index(i,2));
   out(2,i) = u_n(rp_mat_index(i,1),rp_mat_index(i,2));
end

% Initialise test function v

v=zeros(Nr+1,Nt);
%% Main loop
for n=3:NF
if mod(n,50000)==0
NF-n
end
    %%%% v=∆u, compute Laplacian of u
   for l=1:Nr-1
       L=l+1;
     for m=0:Nt-1
         M=mod(m,Nt)+1;

         v(L,M)=((- 2* l* u_n(L,mod(m,Nt)+1) + ...
    (l - 0.5)* u_n(L-1,mod(m,Nt)+1) + (l + 0.5)* u_n(L+1,mod(m,Nt)+1))/(hr^2* l) + ...
    (- 2* u_n(L,mod(m,Nt)+1) + u_n(L,mod(m-1,Nt)+1) + u_n(L,mod(m+1,Nt)+1))/(hr^2* ht^2* l^2));

    end
    
   end

    %%% Boundary condition at l=Nr
    %v=∂rr(u)+∂thth(u)/r^2

   for l=Nr
       L=l+1;
     for m=0:Nt-1
         M=mod(m,Nt)+1;

 v(L,M)=(2*u_n(L-1,mod(m,Nt)+1))/(hr^2); %u_n(L+1,mod(m,Nt)+1)->u_n(L-1,mod(m,Nt)+1) as ∂r(u)=0


    end
    
   end

%.........Center point v..............

lapl_00=0;
  for m=0:Nt-1
      M=mod(m,Nt)+1;
 lapl_00=lapl_00+(4/((Nt)*hr^2))*(u_n(2,M)-u_n(1,M)); %forward first order aproximation

  end

for m=1:Nt
  v(1,m)=lapl_00;
end

%%% ∆∆u=∆v, compute biharmonic of u as the Laplacian of v
    for l=1:Nr-1
       L=l+1;
     for m=0:Nt-1
         M=mod(m,Nt)+1;

u_np1(L,M)=(-k^2*K_mat(L,M)^2)/(1+k*sig0)*((- 2* l* v(L,mod(m,Nt)+1) + ...
    (l - 0.5)* v(L-1,mod(m,Nt)+1) + (l + 0.5)* v(L+1,mod(m,Nt)+1))/(hr^2* l) + ...
    (- 2* v(L,mod(m,Nt)+1) + v(L,mod(m-1,Nt)+1) + v(L,mod(m+1,Nt)+1))/(hr^2* ht^2* l^2))...
    +2/(k*sig0+1)*u_n(L,M)+(k*sig0-1)/(k*sig0+1)*u_nm1(L,M);

    end
    
    end

%.........Center point u..............

%%%%%%% l=0, L=1
  %u00
%l=1,L=2
lapl_00=0;
  for m=0:Nt-1
      M=mod(m,Nt)+1;
 lapl_00=lapl_00+(4/((Nt)*hr^2))*(v(2,M)-v(1,M));

  end

 for m=1:Nt
  u_np1(1,m)=(-(k^2*K_mat(1,m)^2)/(1+k*sig0))*lapl_00+...
      2/(k*sig0+1)*u_n(1,m)+(k*sig0-1)/(k*sig0+1)*u_nm1(1,m);
 end


     if size(unique(u_n(1,:)))~=size(1)

         stop
     end

     if any(abs(u_n)>1e3)
         u_n(1,1)
         'The scheme diverges'
         stop
     end

% Store information in the output

for i=1:size(rp_mat,1)
   out(n,i) = u_np1(rp_mat_index(i,1),rp_mat_index(i,2));
end

% Update grid funtion values

u_nm1 = u_n; u_n = u_np1;

% Animation

if anim==1
   if mod(n,n_anim)==0
         NF-n
    [X,Y,u_n_cart]=pol2cart(T,R,u_n);
      figure(2)
      tiledlayout(1,3)
        nexttile
   mesh(X',Y',u_n_cart)
    xlabel('$x[m]$','Fontsize',17,'Interpreter','latex')
   ylabel('$y[m]$','Fontsize',17,'Interpreter','latex')
   zlabel('$u(x,y,t)[m]$','Fontsize',17,'Interpreter','latex')
   zlim([-limz,limz])
nexttile
   mesh(R',T',u_n)
   xlabel('$r[m]$','Fontsize',17,'Interpreter','latex')
   ylabel('$\theta$','Fontsize',17,'Interpreter','latex')
   zlabel('$u(r,\theta,t)$','Fontsize',17,'Interpreter','latex')
   zlim([-limz,limz])
   title(strcat('Time $t=',num2str(n*k,'%.4f'),'\times',time_scale,'$'),'Fontsize',20,'Interpreter','latex')
nexttile
hold on
   plot(hr*[0:Nr]*Di_ext/2,u_n(1:Nr+1,Nt),'o-','DisplayName','$u(r,\theta=0,t)$')
   plot(hr*[0:Nr]*Di_ext/2,K_mat(1:Nr+1,1)/K_max*limz,'*-','DisplayName','$\kappa(r,\theta=0)$ (scale)')
   legend('Fontsize',17,'Interpreter','latex')
   ylim([-limz,limz])
   xlabel('$r[m]$','Fontsize',17,'Interpreter','latex')
   ylabel('$u(r,\theta=0,t)$','Fontsize',17,'Interpreter','latex')
   hold off
x0=10;
y0=50;
width=1400;
height=350;
set(gcf,'position',[x0,y0,width,height])

   drawnow limitrate
   end
end

end
%%%%%% end main loop

%% Plot waveforms
if ploting==1
out_sum=zeros(size(out(:,i)));
t=[1:NF]*k;
figure(3)
tiledlayout(2,1)
nexttile
hold on
for i=1:size(rp_mat,1)
    leyenda=strcat('$(r,\theta)=$','(',num2str(rp_mat(i,1),'%.4f'),',',num2str(rp_mat(i,2),'%.4f'),')');
    plot(t,out(:,i),'DisplayName',leyenda)
    out_sum=out_sum+out(:,i);
end
lgd=legend('Fontsize',17,'Interpreter','latex');
title(lgd,'Read out position');
ylabel('Amplitude $[m]$','Fontsize',20,'Interpreter','latex')
xlabel(strcat('Time $[\sigma]=[',time_scale,']$'),'Fontsize',20,'Interpreter','latex')
nexttile
plot(t,out_sum,'DisplayName','Sum of all readout points')
ylabel('Amplitude $[m]$','Fontsize',20,'Interpreter','latex')
xlabel(strcat('Time $[\sigma]=[',time_scale,']$'),'Fontsize',20,'Interpreter','latex')
legend('Fontsize',17,'Interpreter','latex')
x0=10;
y0=500;
width=1300;
height=8500;
set(gcf,'position',[x0,y0,width,height])
end
%% Plot FFT
if plot_fft==1
figure(4)
tiledlayout(size(rp_mat,1),1)
out_sum=zeros(size(out(:,i)));
for i=1:size(rp_mat,1)
out_sum=out_sum+out(:,i);
Y = fft(out(:,i));       % Fourier transform
P2 = abs(Y/NF);      
P1 = P2(1:NF/2+1);
P1(2:end-1) = 2*P1(2:end-1); 

f = 1e-6/(sqrt(E_scale_param))*SR*(0:(NF/2))/NF;  % if [E]˜=E_scale_param


[M,In]=max(P1);

nexttile

findpeaks(P1(1:floor(end/5)),f(1:floor(end/5)),'MinPeakDistance',0.5,...
            'MinPeakHeight',M/70)
set(gca,'FontSize',80/size(rp_mat,1));
xlabel('Frequency (MHz)','FontSize',88/size(rp_mat,1),'Interpreter','latex'); 
ylabel('$|\hat{u}|$','FontSize',88/size(rp_mat,1),'Interpreter','latex');
xlim([0,f(In)+50])
read_out_point=strcat('$(r,\theta)=$','(',num2str(rp_mat(i,1),'%.4f'),',',num2str(rp_mat(i,2),'%.4f'),')');
freq_res=strcat('Fourier Transform $f_r$=',num2str(f(In),'%.4f'),'MHz');
tit=[freq_res newline read_out_point];
title(tit,'FontSize',80/size(rp_mat,1),'Interpreter','latex');
end
x0=500;
y0=500;
width=800;
height=1300;
set(gcf,'position',[x0,y0,width,height])
Y = fft(out_sum);       %  Fourier transform
P2 = abs(Y/NF);      
P1 = P2(1:NF/2+1);
P1(2:end-1) = 2*P1(2:end-1); 

f = 1e-6/(sqrt(E_scale_param))*SR*(0:(NF/2))/NF;  % if [E]˜=E_scale_param

[M,In]=max(P1);

figure(5)
findpeaks(P1(1:floor(end/5)),f(1:floor(end/5)),'MinPeakDistance',0.5,...
            'MinPeakHeight',M/50,'Annotate','extents')

[pks,locs]=findpeaks(P1(1:floor(end/5)),f(1:floor(end/5)),'MinPeakDistance',2,...
            'MinPeakHeight',M/10,'Annotate','extents');
set(gca,'FontSize',20);
xlabel('Frequency (MHz)','FontSize',22,'Interpreter','latex'); 
ylabel('$|\hat{u}|$','FontSize',22,'Interpreter','latex');
xlim([0,max(locs)+5])
freq_res=strcat('$f_r$=',num2str(f(In),'%.4f'),'MHz');
tit=[freq_res newline 'Sum of all readout poins'];
title(tit,'FontSize',24,'Interpreter','latex');
x0=0;
y0=500;
width=800;
height=300;
set(gcf,'position',[x0,y0,width,height])
end
toc % Stop measuring execution time
