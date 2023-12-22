function [out,SR,hr,ht,k]=...
Polar_plate_2lapl(Di_ext,H_ext,E_ext,rho_ext,nu_ext,Di_int,H_int,E_int,rho_int,nu_int,E_scale_param,...
kstability_constant,ctr,wid,u0,v0,sig0,k_sigmoid,logistic_fun,Nr,Nt,TF,rp_mat,ploting,anim,n_anim,limz,plot_fft,time_series)

tic % For measuring the executing time of the program
time_scale=strcat('10^{',num2str(log10(sqrt(E_scale_param))),'}s');

E_ext_scale=E_ext*E_scale_param; %Scale the Young modulus
E_int_scale=E_int*E_scale_param;%Scale the Young modulus

%%% Compute Stiffness parameters

D_ext=(E_ext_scale)*H_ext^3/(12*(1-nu_ext^2));
K_ext=sqrt(D_ext/(rho_ext*H_ext*(Di_ext/2)^4));

HT=H_int+H_ext;

rho_l=HT/(H_ext/rho_ext+H_int/rho_int);
E_l=HT/(H_ext/E_ext_scale+H_int/E_int_scale);
nu_l=HT/(H_ext/nu_ext+H_int/nu_int);

D_int=E_l*HT^3/(12*(1-nu_l^2));
K_int=sqrt(D_int/(rho_l*HT*(Di_ext/2)^4));


K_max=max([K_int,K_ext]);

K_min=min([K_int,K_ext]);

R_int=Di_int/Di_ext;

%%% Radial coordinates

hr = 1/Nr;  % radial step size

Nr=floor(1/hr);

ht=2*pi/Nt; % angular step size

Nt=round(2*pi/ht);

%%% Stability condition 

 %plate free center
 k=kstability_constant*(hr^2/(2*K_max))*(1/(1+1/(ht^2)));
%Sample Rate
SR=1/k;

%Number of Time Steps
NF = floor(SR*TF);  

%Polar coordinates axes
[R, T] = meshgrid([0:Nr]*hr*Di_ext/2, [0:Nt-1]*ht);

%Cartesian coordinates axes 

[X,Y]=pol2cart(T,R);

% Compute readout parameters index

rp_mat_index=zeros(size(rp_mat));
for i=1:size(rp_mat,1)
    rp=rp_mat(i,:);
    rp_r=rp(1);
    rp_t=rp(2);
    rpr_index=floor(rp_r/hr)+1;
    rpt_index=floor(rp_t/ht)+1;
    rp_mat_index(i,:)=[rpr_index rpt_index];
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


%Fourth order polinomial Initial position
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

%%

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

if time_series==1
   figure(8)
       tiledlayout(2,4)
end
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
%  lapl_00=lapl_00+(4/((Nt)*hr^2))*(-(11/6)*u_n(1,M)+3*u_n(2,M)-(3/2)*u_n(3,M)+(1/3)*u_n(4,M)); %forward third order aproximation
 
%  lapl_00=lapl_00+(4/((Nt)*hr^2))*(-(25/12)*u_n(1,M)+4*u_n(2,M)-(3)*u_n(3,M)+(4/3)*u_n(4,M)-(1/4)*u_n(5,M)); %forward fourth order aproximation

%  lapl_00=lapl_00+(2/((Nt)*hr^2))*(u_n(2,M)-u_n(2,mod(m,floor(Nt/2))+1)); %centered first order aproximation

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
%  lapl_00=lapl_00+(4/((Nt)*hr^2))*(-(11/6)*v(1,M)+3*v(2,M)-(3/2)*v(3,M)+(1/3)*v(4,M)); %forward second order aproximation

%  lapl_00=lapl_00+(4/((Nt)*hr^2))*(-(25/12)*v(1,M)+4*v(2,M)-(3)*v(3,M)+(4/3)*v(4,M)-(1/4)*v(5,M)); %forward fourth order aproximation


%  lapl_00=lapl_00+(2/((Nt)*hr^2))*(v(2,M)-v(2,mod(m,floor(Nt/2))+1)); %centered first order aproximation

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
         'El esquema diverge'
         stop
     end

% Store information in the output

for i=1:size(rp_mat,1)
   out(n,i) = u_np1(rp_mat_index(i,1),rp_mat_index(i,2));
end



u_nm1 = u_n; u_n = u_np1;


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
%    time_scale=strcat('10^{',num2str(log10(sqrt(E_scale_param))),'}s$');
   title(strcat('Time $t=',num2str(n*k,'%.4f'),'\times',time_scale,'$'),'Fontsize',20,'Interpreter','latex')
nexttile
hold on
   plot(hr*[0:Nr]*Di_ext/2,u_n(1:Nr+1,Nt),'o-','DisplayName','$u(r,\theta=0,t)$')
   plot(hr*[0:Nr]*Di_ext/2,K_mat(1:Nr+1,1)/K_max*limz,'*-','DisplayName','$\kappa(r,\theta=0)$ (scaled)')
   legend('Fontsize',17,'Interpreter','latex')
   ylim([-limz,limz])
   xlabel('$r[m]$','Fontsize',17,'Interpreter','latex')
   ylabel('$u(r,\theta=0,t)$','Fontsize',17,'Interpreter','latex')
   hold off
x0=10;
y0=75;
width=1400;
height=300;
set(gcf,'position',[x0,y0,width,height])

   drawnow limitrate
   end
end

if time_series==1
%        a=[50,100,200,300,400,500,600,700];
%         if ismember(n,a)
[X,Y,u_n_cart]=pol2cart(T,R,u_n);
for j=0:7
    if n*k>=0.002+j*0.03 && n*k<0.002+j*0.03+k
       %if (n*k>=0.002 && n*k<0.002+k) || (n*k>=0.004 && n*k<0.04+k) || (n*k>=0.006 && n*k<0.006+k) || (n*k>=0.003 && n*k<0.003+k) || (n*k>=0.004 && n*k<0.004+k) || (n*k>=0.005 && n*k<0.005+k) || (n*k>=0.006 && n*k<0.006+k) || (n*k>=0.007 && n*k<0.007+k)
       nexttile
       
       mesh(X',Y',u_n_cart)
    xlabel('$x[m]$','Fontsize',17,'Interpreter','latex')
   ylabel('$y[m]$','Fontsize',17,'Interpreter','latex')
   zlabel('$u(x,y,t)[m]$','Fontsize',17,'Interpreter','latex')
   title(strcat('$t=',num2str(n*k,'%.3f'),'\times',time_scale,'$'),'Fontsize',20,'Interpreter','latex')
   zlim([-limz,limz])

    end
end
end

end
%%%%%% end main loop
if time_series==1
x0=0;
        y0=50;
        width=1500;
        height=600;
        set(gcf,'position',[x0,y0,width,height])
end
%% Plot waveforms
if ploting==1
out_sum=zeros(size(out(:,i)));
t=[1:NF]*k;
figure(3)
% tiledlayout(2,1)
% nexttile
hold on
for i=1:size(rp_mat,1)
    leyenda=strcat('$(r,\theta)=$','(',num2str(rp_mat(i,1),'%.4f'),',',num2str(rp_mat(i,2),'%.4f'),')');
    plot(t,out(:,i),'DisplayName',leyenda)
    out_sum=out_sum+out(:,i);
end
lgd=legend('Fontsize',20,'Interpreter','latex');
title(lgd,'Read out position');
ylabel('Amplitude $(\approx u(r,\theta))$ $[m]$','Fontsize',30,'Interpreter','latex')
xlabel(strcat('Time $[\sigma]=[',time_scale,']$'),'Fontsize',30,'Interpreter','latex')
% nexttile
% plot(t,out_sum,'DisplayName','Sum of all readout points')
% ylabel('Amplitude $[m]$','Fontsize',20,'Interpreter','latex')
% xlabel(strcat('Time $[\sigma]=[',time_scale,']$'),'Fontsize',20,'Interpreter','latex')
% legend('Fontsize',17,'Interpreter','latex')
x0=10;
y0=500;
width=1300;
height=500;
set(gcf,'position',[x0,y0,width,height])
end
%% Plot FFT
if plot_fft==1
figure(4)
tiledlayout(size(rp_mat,1),1)
out_sum=zeros(size(out(:,i)));
for i=1:size(rp_mat,1)
out_sum=out_sum+out(:,i);
Y = fft(out(:,i));       %Transformada de Fourier
P2 = abs(Y/NF);      
P1 = P2(1:NF/2+1);
P1(2:end-1) = 2*P1(2:end-1); 

f = 1e-6/(sqrt(E_scale_param))*SR*(0:(NF/2))/NF;  % if [E]˜=E_scale_param


[M,In]=max(P1);

nexttile
hold on
read_out_point=strcat('$(r,\theta)=$','(',num2str(rp_mat(i,1),'%.4f'),',',num2str(rp_mat(i,2),'%.4f'),')');
p=plot(f(1:floor(end/5)),P1(1:floor(end/5)));
p.Annotation.LegendInformation.IconDisplayStyle = 'off';
[pks,locs] =findpeaks(P1(1:floor(end/5)),f(1:floor(end/5)),'MinPeakDistance',0.5,...
            'MinPeakHeight',M/70);

[Mp,Inp]=max(pks);
% Index numbers of elements to remove
indicesToRemove = [1, Inp];
% Create a logical index to keep elements not in the indicesToRemove list
keepIndices = true(1, numel(pks));
keepIndices(indicesToRemove) = false;
pks_other=pks(keepIndices);
locs_other=locs(keepIndices);
% set(gca,'FontSize',80/size(rp_mat,1));
xlabel('Frequency (MHz)','FontSize',88/size(rp_mat,1),'Interpreter','latex'); 
ylabel('$|\hat{u}|$','FontSize',88/size(rp_mat,1),'Interpreter','latex');
xlim([0,f(In)+50])
freq_res=strcat('$f_{max}$=',num2str(f(In),'%.4f'),'MHz');
f_1st=strcat('$f_{1^{st}}$=',num2str(locs(1),'%.4f'),'MHz');
scatter(locs_other,pks_other,'^','DisplayName','Other modes')
scatter(locs(1),pks(1),'o','DisplayName',f_1st)
scatter(locs(Inp),pks(Inp),'*','DisplayName',freq_res)
lgd=legend('Fontsize',17,'Interpreter','latex');
title(lgd,read_out_point);
hold off
% tit=[freq_res newline read_out_point];
% title(tit,'FontSize',80/size(rp_mat,1),'Interpreter','latex');
end
x0=500;
y0=500;
width=800;
height=1300;
set(gcf,'position',[x0,y0,width,height])
% Y = fft(out_sum);       %Transformada de Fourier
% P2 = abs(Y/NF);      
% P1 = P2(1:NF/2+1);
% P1(2:end-1) = 2*P1(2:end-1); 
% 
% f = 1e-6/(sqrt(E_scale_param))*SR*(0:(NF/2))/NF;  % if [E]˜=E_scale_param
% 
% [M,In]=max(P1);
% 
% figure(5)
% findpeaks(P1(1:floor(end/5)),f(1:floor(end/5)),'MinPeakDistance',0.5,...
%             'MinPeakHeight',M/50,'Annotate','extents')
% 
% [pks,locs]=findpeaks(P1(1:floor(end/5)),f(1:floor(end/5)),'MinPeakDistance',2,...
%             'MinPeakHeight',M/10,'Annotate','extents');
% set(gca,'FontSize',20);
% xlabel('Frequency (MHz)','FontSize',22,'Interpreter','latex'); 
% ylabel('$|\hat{u}|$','FontSize',22,'Interpreter','latex');
% xlim([0,max(locs)+5])
% freq_res=strcat('$f_r$=',num2str(f(In),'%.4f'),'MHz');
% tit=[freq_res newline 'Sum of all readout poins'];
% title(tit,'FontSize',24,'Interpreter','latex');
% x0=0;
% y0=500;
% width=800;
% height=300;
% set(gcf,'position',[x0,y0,width,height])
end
toc
