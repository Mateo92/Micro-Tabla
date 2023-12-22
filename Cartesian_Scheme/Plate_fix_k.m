function [out,SR,h,ss,k]=...
     Plate_fix_k(Di_ext,H_ext,E_ext,rho_ext,nu_ext,Di_int,H_int,E_int,rho_int,nu_int,...
    ch_len,ch_wid,ctr,wid,u0,v0,sig0,k_sigmoid,k,TF,E_scale_param,rp_mat,iv,ploting,anim,n_anim,plot_fft)

tic
time_scale=strcat('10^{',num2str(log10(sqrt(E_scale_param))),'}s');

E_ext_scale=E_ext*E_scale_param;
E_int_scale=E_int*E_scale_param;

D_ext=(E_ext_scale)*H_ext^3/(12*(1-nu_ext^2));
K_ext=sqrt(D_ext/(rho_ext*H_ext*(Di_ext)^4));

HT=H_int+H_ext;

rho_l=HT/(H_ext/rho_ext+H_int/rho_int);
E_l=HT/(H_ext/E_ext_scale+H_int/E_int_scale);
nu_l=HT/(H_ext/nu_ext+H_int/nu_int);

D_int=E_l*HT^3/(12*(1-nu_l^2));
K_int=sqrt(D_int/(rho_l*HT*(Di_ext)^4));


%%begin global parameters
K1=K_ext;  % unloaded plate stiffness parameter (1/sig)
K2 =K_int; % loaded plate stiffness parameter (1/sig)
% stability condition/scheme parameters
K=max([K1,K2]);    
K_min=min([K1,K2]);
mu = 0.25;          % scheme free parameter
% begin derived parameters
h = sqrt(K*k/mu); % find grid spacing
Nx = floor(1/h)+4; % number of x-subdivisions of spatial domain
SR=floor(1/k); % sample rate

%h = sqrt(K*k/mu); % find grid spacing
%Nx = floor(1/h)+4; % number of x-subdivisions of spatial domain 
%Ny = floor(1/h)+4; % number of y-subdivisions of spatial domain
%Ny=Nx;
y_extra=floor(ch_len/(h*Di_ext));
Ny = Nx+y_extra; 

%T60 = 3;       % 60 dB decay time (s)
%epsilon = 1; % domain aspect ratio 
       % duration of simulation(sig)
NF = floor(SR*TF);              % duration of simulation (samples)
[X, Y] = meshgrid([1:Nx-1]*h, [1:Ny-1]*h);
Xaxis=[1:Nx-1]*h;
Yaxis=[1:Ny-1]*h;
ss = (Nx-1)*(Ny-1); % total grid size

%%%Boundary condition
Bo=zeros((Ny-1),(Nx-1));
for i=1:Nx-1
    for j=1:Ny-1
        if sqrt((i-round((Nx-1)/2))^2+(j-round((Ny-1)/2))^2)<=(Nx-5)/2
            Bo(j,i)=1;
        end
    end
end

Bo_rect=zeros((Ny-1),(Nx-1));
for i=1:Nx-1
%     for j=1:Ny-1
    for j=2:Ny-2
        if abs(h*Di_ext*(i-(Nx-1)/2))<ch_wid/2
            Bo_rect(j,i)=1;
        end
    end
end
 Bo=Bo+Bo_rect;
 for i=1:Nx-1
    for j=1:Ny-1
        if Bo(j,i)>0
            Bo(j,i)=1;
        end
    end
 end
%Stiffness parameter distribution
K_Mat=zeros((Ny-1),(Nx-1));
R2_mat=(Di_int/Di_ext)*(Nx-1)/2;
%ka=2;
for i=1:Nx-1
    for j=1:Ny-1

        r=sqrt((i-(Nx-1)/2)^2+(j-(Ny-1)/2)^2);

        logis=K-(K-K_min)*(1/(1+exp(-k_sigmoid*(r-R2_mat))));

        K_Mat(j,i)=logis;

        % Step Function Distribution
%         if sqrt((i-(Nx-1)/2)^2+(j-(Ny-1)/2)^2)<=Nx*(D2/2)
%             K_Mat(i,j)=K;
%         end

    end
end
Xax=Di_ext*X;
Yax=Di_ext*Y;
K_vect=reshape(K_Mat, ss,1);
K_vect2=K_vect.^2;
%% generate difference matrix/scheme matrices
Dxx = sparse(toeplitz([-2/h^2;1/h^2;zeros(Nx-3,1)]));
Dyy = sparse(toeplitz([-2/h^2;1/h^2;zeros(Ny-3,1)]));
D = kron(eye(Nx-1), Dyy)+kron(Dxx, eye(Ny-1)); DD = D*D;
B = (1/(1+sig0*k))*(2*eye(ss)-k^2*K^2*DD);
I2=sparse((1/(1+sig0*k))*2*eye(ss));
DD2=sparse((1/(1+sig0*k))*k^2*DD);
C = ((1-sig0*k)/(1+sig0*k))*sparse(eye(ss));
%No Damping:
%C = 1*sparse(eye(ss));
% readout interpolation parameters
rp_index_vect=zeros(size(rp_mat,1));
for i=1:size(rp_mat,1)
    rp=rp_mat(i,1:end);
rp_index= (Ny-1)*floor(rp(1)*Nx)+floor(rp(2)*Ny); 
rp_index_vect(i)=rp_index;
end
% create 2D raised cosine
%dist = sqrt((X-ctr(1)).^2+(Y-ctr(2)).^2); 
dist = sqrt((X-ctr(1)*max(max(X))).^2+(Y-ctr(2)*max(max(Y))).^2); 
ind = sign(max(-dist+wid/2,0)); 
rc = 0.5*ind.*(1+cos(2*pi*dist/wid)); 
rc = reshape(rc, ss,1);
RC=reshape(rc,[Ny-1,Nx-1]);
%Cuadratic I.C.
cu=zeros((Ny-1),(Nx-1));
for i=1:Nx-1
    for j=1:Ny-1
        r=sqrt((i-(Nx-1)/2)^2+(j-(Ny-1)/2)^2);
        if r*h*Di_ext<=Di_ext/2
        cu(j,i)=(Di_ext/2)^2-(r*h*Di_ext)^2;
        end
    end
end
cu = reshape(cu, ss,1);
%Fourth I.C.
qu=zeros((Ny-1),(Nx-1));
for i=1:Nx-1
    for j=1:Ny-1
        r=sqrt((i-(Nx-1)/2)^2+(j-(Ny-1)/2)^2);
        if r*h*Di_ext<=Di_ext/2
        qu(j,i)=((Di_ext/2)^2-(r*h*Di_ext)^2)^2;
        end
    end
end
%Draw IC:
% figure(1)
% tiledlayout(2,1)
% nexttile
% mesh(Xax,Yax,qu)
% nexttile
% mesh(Xax,Yax,u0*qu)
% 
% qu = reshape(qu,ss,1);
% 
% max(u0*qu)
% stop

qu=reshape(qu,ss,1);
%Poit Initial velocity
point_iv=zeros((Ny-1),(Nx-1));
point_iv(floor((Nx-1)*ctr(1)),floor((Ny-1)*ctr(2)))=1;
point_iv=reshape(point_iv,[Ny-1,Nx-1]);
% set initial conditions/initialize output
u2 = u0*rc; 
if iv==1
    u1 = (u0+k*v0)*point_iv;
else
u1 = (u0+k*v0)*rc;
end
if iv==2
    u1 = 0*cu;
    %u1 = (u0+k*v0)*cu;
    u2 = u0*cu; 
else
u1 = (u0+k*v0)*rc;
end
if iv==3
    u1 = u0*qu+k*v0*rc;
    %u1=0*qu;
    u2 = u0*qu; 
else
u1 = (u0+k*v0)*rc;
end
% limitey=max([max(u2),max(u1),1e-5]);
limitey=1e-4;
u = zeros(ss,1);
out = zeros(NF,length(rp_mat));

U1=reshape(u1,[Ny-1,Nx-1]);
BU1=Bo.*U1;
u1=reshape(BU1, ss,1);

if ploting==1
    figure(1)
    tiledlayout(1,4)
    nexttile
%     mesh(Xax,Yax,Bo)
    imagesc(Xaxis,Yaxis,Bo)
    colorbar
    title('Boundary Matrix','Fontsize',20,'Interpreter','latex')
    xlabel('$x[m]$','Fontsize',20,'Interpreter','latex')
    ylabel('$y[m]$','Fontsize',20,'Interpreter','latex')
    zlabel('Boundary','Fontsize',20,'Interpreter','latex')
    nexttile
    mesh(Xax,Yax,K_Mat)
    title('Stiffness parameter distribution','Fontsize',20,'Interpreter','latex')
    xlabel('$x[m]$','Fontsize',20,'Interpreter','latex')
    ylabel('$y[m]$','Fontsize',20,'Interpreter','latex')
    zlabel('$\kappa(x,y)[\sigma^{-1}]$','Fontsize',20,'Interpreter','latex')
    nexttile
    U=reshape(u2,[Ny-1,Nx-1]);
    mesh(Xax,Yax,U,'FaceColor','interp')
    title('Initial Position','Fontsize',20,'Interpreter','latex')
    xlabel('$x[m]$','Fontsize',20,'Interpreter','latex')
    ylabel('$y[m]$','Fontsize',20,'Interpreter','latex')
    zlabel('$u(x,y,t=0)[m]$','Fontsize',20,'Interpreter','latex')
    nexttile
    mesh(Xax,Yax,RC*v0)%,'FaceColor','interp')
    title('Initial velocity distribution','Fontsize',20,'Interpreter','latex')
    xlabel('$x[m]$','Fontsize',20,'Interpreter','latex')
    ylabel('$y[m]$','Fontsize',20,'Interpreter','latex')
    zlabel('$u_t(x,y,t=0)[m\sigma^{-1}]$','Fontsize',20,'Interpreter','latex')
    x0=10;
    y0=500;
    width=1400;
    height=350;
    set(gcf,'position',[x0,y0,width,height])
end
%% %%%% start main loop
for n=3:NF
if mod(n,50000)==0
NF-n
end   

   %Finnite difference scheme:
   u = I2*u1 - DD2*(u1.*K_vect2) - C*u2;% -D*u1;
   
   %Clamped Boundary:
   U=reshape(u,[Ny-1,Nx-1]);
   BU=Bo.*U;
   u=reshape(BU,ss,1);
   
  %Ploting
   if anim==1
   if mod(n,n_anim)==0
        figure(2)
       tiledlayout(1,2)
       nexttile
       U=reshape(u,[Ny-1,Nx-1]);
%        mesh(Xax,Yax,U,'FaceColor','interp')
       mesh(Xax,Yax,U)

   title(strcat('Time $t=',num2str(n*k,'%.4f'),'\times',time_scale,'$'),'Fontsize',20,'Interpreter','latex')
       xlabel('$x[m]$','Fontsize',20,'Interpreter','latex')
       ylabel('$y[m]$','Fontsize',20,'Interpreter','latex')
       zlabel('$u(x,y,t)$','Fontsize',20,'Interpreter','latex')
       %zlim([-1e-5,1e-5])
       zlim([-limitey*1.5,limitey*1.5])
       nexttile
       hold on
       plot([1:Nx-1]*Di_ext*h,U(round(Ny/2),1:Nx-1),'o-','DisplayName',strcat('$u(x,y=',num2str(0.5*Di_ext,'%.2e'),',t)$'))
       plot([1:Nx-1]*Di_ext*h,K_Mat(floor(Ny/2),1:Nx-1)/K*limitey*1.5,'*-','DisplayName',strcat('$\kappa(x,y=',num2str(0.5*Di_ext,'%.2e'),',t)$(scaled)'))
       xlabel('$x[m]$','Fontsize',20,'Interpreter','latex')
       legend('Fontsize',17,'Interpreter','latex','Location','southeast')
       %ylim([-1e-5,1e-5])
       ylim([-limitey*1.5,limitey*1.5])
       title('Centred side view','Fontsize',20,'Interpreter','latex')
       hold off
       x0=250;
        y0=50;
        width=1000;
        height=500;
        set(gcf,'position',[x0,y0,width,height])
       drawnow limitrate
    end
   end


   u2 = u1; u1 = u;
    for i=1:size(rp_mat,1)
   out(n,i) = u(rp_index_vect(i));
    end
  
end

%%%%%% end main loop
%%
if ploting==1

t=[1:NF]*k;
figure(3)
hold on
for i=1:size(rp_mat,1)
    leyenda=strcat('$(x,y)=$','(',num2str(rp_mat(i,1)*Di_ext,'%.2e'),',',num2str(rp_mat(i,2)*Di_ext,'%.2e'),')');
    plot(t,out(:,i),'DisplayName',leyenda)

end
lgd=legend('Fontsize',17,'Interpreter','latex');
title(lgd,'Read out position');
ylabel('Amplitude $[m]$','Fontsize',20,'Interpreter','latex')
% xlabel('Time $[\sigma]=[10^{-4}s]$','Fontsize',20,'Interpreter','latex')
xlabel(strcat('Time $[\sigma]=[',time_scale,']$'),'Fontsize',20,'Interpreter','latex')

ylim([-limitey,limitey])
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
for i=1:size(rp_mat,1)
Y = fft(out(:,i));       %Transformada de Fourier
P2 = abs(Y/NF);      
P1 = P2(1:NF/2+1);
P1(2:end-1) = 2*P1(2:end-1); 
% f = 1e-2*SR*(0:(NF/2))/NF;
f = 1e-6/(sqrt(E_scale_param))*SR*(0:(NF/2))/NF;  % if [E]˜=E_scale_param
[M,In]=max(P1);

nexttile

findpeaks(P1(1:floor(end/5)),f(1:floor(end/5)),'MinPeakDistance',2,...
            'MinPeakHeight',M/10,'Annotate','extents')
set(gca,'FontSize',20);
xlabel('Frequency (MHz)','FontSize',20,'Interpreter','latex'); 
ylabel('$|\hat{u}|$','FontSize',20,'Interpreter','latex');
xlim([0,f(In)+100])
read_out_point=strcat('$(x,y)=$','(',num2str(rp_mat(i,1)*Di_ext,'%.2e'),',',num2str(rp_mat(i,2)*Di_ext,'%.2e'),')');
freq_res=strcat('Fourier Transform $f_r$=',num2str(f(In)),'MHz');
tit=[freq_res newline read_out_point];
title(tit,'FontSize',20,'Interpreter','latex');
end
x0=500;
y0=500;
width=800;
height=1300;
set(gcf,'position',[x0,y0,width,height])
end
toc





