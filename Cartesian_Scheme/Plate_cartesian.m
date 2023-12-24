function [out, SR, h, ss, k] = Plate_cartesian(Di_ext, H_ext, E_ext, rho_ext, nu_ext, ...
    Di_int, H_int, E_int, rho_int, nu_int, ch_len, ch_wid, ctr, wid, u0, v0, ...
    sig0, k_sigmoid, Nx, TF, sigma_time, rp_mat, ploting, anim, n_anim, limz, plot_fft)
%
% Simulate a plate using finite difference method in Cartesian coordinates
%
% Inputs:
% Di_ext, Di_int: External and internal diameters of the plate [m].
% H_ext, H_int: External and internal heights of the plate [m].
% E_ext, E_int: Young's moduli of the external and internal materials [Pa=kg/(ms^2)].
% rho_ext, rho_int: Densities of the external and internal materials [kg/m^3].
% nu_ext, nu_int: Poisson's ratios of the external and internal materials.
% ch_len, ch_wid: Length and width of the chanels in the y direction
% ctr: Center coordinates for boundary conditions [0-1,0-1] as a fraction 
% of the plate's diameter.
% wid: Width parameter for the boundary conditions [0-1,0-1] as a fraction
% of the plate's diameter.
% u0: Initial displacement in [m].
% v0: Initial velocity distribution in height in [m/sigma^2].
% sig0: Parameter for the finite difference scheme.
% k_sigmoid: Parameter for a sigmoid function.
% Nx: Number of spatial subdivisions in the x-direction.
% TF: Final simulation time.
% sigma_time: Scale parameter for Young's modulus.
% rp_mat: Matrix specifying readout positions ([0-1,0-1]) as a fraction 
% of the plate's diameter.
% ploting: Flag for enabling plot generation.
% anim: Flag for enabling animation generation.
% n_anim: Parameter for animation frequency.
% limitey: Limit for y-axis in plots.
% plot_fft: Flag for enabling Fast Fourier Transform (FFT) plots.
%
% Outputs:
% out: Simulation output for readout positions over time.
% SR: Sample Rate calculated based on stability condition.
% h: Grid spacing.
% ss: Total grid size.
% k: Time step for the simulation.
%

% Start measuring execution time
tic;

% Calculate scale parameter for Young's modulus
E_scale_param = sigma_time^2;

% Calculate time scale for display
time_scale = strcat('10^{', num2str(log10(sqrt(E_scale_param))), '}s');

% Calculate scaled Young's moduli
E_ext_scale = E_ext * E_scale_param;
E_int_scale = E_int * E_scale_param;

% Calculate stiffness parameters for external and internal materials
D_ext = (E_ext_scale) * H_ext^3 / (12 * (1 - nu_ext^2));
K_ext = sqrt(D_ext / (rho_ext * H_ext * (Di_ext)^4));

HT = H_int + H_ext;

rho_l = HT / (H_ext / rho_ext + H_int / rho_int);
E_l = HT / (H_ext / E_ext_scale + H_int / E_int_scale);
nu_l = HT / (H_ext / nu_ext + H_int / nu_int);

D_int = E_l * HT^3 / (12 * (1 - nu_l^2));
K_int = sqrt(D_int / (rho_l * HT * (Di_ext)^4));

% Set global parameters
K1 = K_ext;  % unloaded plate stiffness parameter (1/sig)
K2 = K_int;  % loaded plate stiffness parameter (1/sig)

% Stability condition/scheme parameters
K = max([K1, K2]);
K_min = min([K1, K2]);
mu = 0.25;  % scheme free parameter

% Calculate derived parameters
h = 1 / (Nx - 4);
k = mu * h^2 / K;
SR = floor(1 / k); % sample rate

% Calculate number of subdivisions in the y-direction
y_extra = floor(ch_len / (h * Di_ext));
Ny = Nx + y_extra;

% Calculate simulation duration
NF = floor(SR * TF);

% Create meshgrid for spatial coordinates
[X, Y] = meshgrid([1:Nx-1] * h, [1:Ny-1] * h);
Xaxis = [1:Nx-1] * h;
Yaxis = [1:Ny-1] * h;
ss = (Nx - 1) * (Ny - 1); % total grid size

% Boundary condition matrix
Bo = zeros((Ny - 1), (Nx - 1));
for i = 1:Nx-1
    for j = 1:Ny-1
        if sqrt((i - (Nx - 1) / 2)^2 + (j - (Ny - 1) / 2)^2) <= (Nx - 5) / 2
            Bo(j, i) = 1;
        end
    end
end

% Additional rectangular boundary condition
Bo_rect = zeros((Ny - 1), (Nx - 1));
for i = 1:Nx-1
    for j = 2:Ny-2
        if abs(h * Di_ext * (i - (Nx - 1) / 2)) < ch_wid / 2
            Bo_rect(j, i) = 1;
        end
    end
end
Bo = Bo + Bo_rect;
Bo(Bo > 0) = 1;

% Stiffness parameter distribution
K_Mat = zeros((Ny - 1), (Nx - 1));
R2_mat = (Di_int / Di_ext) * (Nx - 1) / 2;

for i = 1:Nx-1
    for j = 1:Ny-1
        r = sqrt((i - (Nx - 1) / 2)^2 + (j - (Ny - 1) / 2)^2);
        logis = K - (K - K_min) * (1 / (1 + exp(-k_sigmoid * (r - R2_mat))));
        K_Mat(j, i) = logis;
    end
end

Xax = Di_ext * X;
Yax = Di_ext * Y;
K_vect = reshape(K_Mat, ss, 1);
K_vect2 = K_vect.^2;

% Finite difference scheme matrices
Dxx = sparse(toeplitz([-2 / h^2; 1 / h^2; zeros(Nx - 3, 1)]));
Dyy = sparse(toeplitz([-2 / h^2; 1 / h^2; zeros(Ny - 3, 1)]));
D = kron(eye(Nx - 1), Dyy) + kron(Dxx, eye(Ny - 1));
DD = D * D; % Laplacian matrix operator
I2 = sparse((1 / (1 + sig0 * k)) * 2 * eye(ss));
DD2 = sparse((1 / (1 + sig0 * k)) * k^2 * DD); %Biharmonic matrix operator
C = ((1 - sig0 * k) / (1 + sig0 * k)) * sparse(eye(ss));

% Readout interpolation parameters
rp_index_vect = zeros(size(rp_mat, 1));
for i = 1:size(rp_mat, 1)
    rp = rp_mat(i, 1:end);
    rp_index = (Ny - 1) * floor(rp(1) * Nx) + floor(rp(2) * Ny);
    rp_index_vect(i) = rp_index;
end

% Create 2D raised cosine distribution
dist = sqrt((X - ctr(1) * max(max(X))).^2 + (Y - ctr(2) * max(max(Y))).^2);
ind = sign(max(-dist + wid / 2, 0));
rc = 0.5 * ind .* (1 + cos(2 * pi * dist / wid));
rc = reshape(rc, ss, 1);
RC = reshape(rc, [Ny - 1, Nx - 1]);

% Initial conditions for different velocity distributions
qu = zeros((Ny - 1), (Nx - 1));

for i = 1:Nx - 1
    for j = 1:Ny - 1
        r = sqrt((i - (Nx - 1) / 2)^2 + (j - (Ny - 1) / 2)^2);
        
        % Fourth Power Initial Position
        if r * h * Di_ext <= Di_ext / 2
            qu(j, i) = ((1 / 2)^2 - (r * h * 1)^2)^2;
        end
        
    end
end

% Reshape initial conditions
qu = reshape(qu, ss, 1);

% Initialize grid functions

u_nm1 = u0 * qu;
u_n = u_nm1 + k * v0 * rc;

% Initialise output signals

out = zeros(NF, length(rp_mat));

for i=1:size(rp_mat,1)
   out(1,i) = u_nm1(rp_index_vect(i));
   out(2,i) = u_n(rp_index_vect(i));
end

U1 = reshape(u_n, [Ny - 1, Nx - 1]);
BU1 = Bo .* U1;
u_n = reshape(BU1, ss, 1);

% Plotting if enabled
if ploting == 1
    figure(1)
    tiledlayout(1, 4)
    nexttile
    imagesc(Xaxis, Yaxis, Bo)
    colorbar
    title('Boundary Matrix', 'Fontsize', 20, 'Interpreter', 'latex')
    xlabel('$x[m]$', 'Fontsize', 20, 'Interpreter', 'latex')
    ylabel('$y[m]$', 'Fontsize', 20, 'Interpreter', 'latex')
    zlabel('Boundary', 'Fontsize', 20, 'Interpreter', 'latex')
    
    nexttile
    mesh(Xax, Yax, K_Mat)
    title('Stiffness parameter distribution', 'Fontsize', 20, 'Interpreter', 'latex')
    xlabel('$x[m]$', 'Fontsize', 20, 'Interpreter', 'latex')
    ylabel('$y[m]$', 'Fontsize', 20, 'Interpreter', 'latex')
    zlabel('$\kappa(x,y)[\sigma^{-1}]$', 'Fontsize', 20, 'Interpreter', 'latex')
    
    nexttile
    U = reshape(u_nm1, [Ny - 1, Nx - 1]);
    mesh(Xax, Yax, U, 'FaceColor', 'interp')
    title('Initial Position', 'Fontsize', 20, 'Interpreter', 'latex')
    xlabel('$x[m]$', 'Fontsize', 20, 'Interpreter', 'latex')
    ylabel('$y[m]$', 'Fontsize', 20, 'Interpreter', 'latex')
    zlabel('$u(x,y,t=0)[m]$', 'Fontsize', 20, 'Interpreter', 'latex')
    
    nexttile
    mesh(Xax, Yax, RC * v0)
    title('Initial velocity distribution', 'Fontsize', 20, 'Interpreter', 'latex')
    xlabel('$x[m]$', 'Fontsize', 20, 'Interpreter', 'latex')
    ylabel('$y[m]$', 'Fontsize', 20, 'Interpreter', 'latex')
    zlabel('$u_t(x,y,t=0)[m\sigma^{-1}]$', 'Fontsize', 20, 'Interpreter', 'latex')
    
    x0 = 10;
    y0 = 500;
    width = 1400;
    height = 350;
    set(gcf, 'position', [x0, y0, width, height])
end

% Main loop
for n = 3:NF
    if mod(n, 50000) == 0
        NF - n
    end

    % Finite difference scheme
    u_np1 = I2 * u_n - DD2 * (u_n .* K_vect2) - C * u_nm1;

    % Clamped Boundary
    U = reshape(u_np1, [Ny - 1, Nx - 1]);
    BU = Bo .* U;
    u_np1 = reshape(BU, ss, 1);

    % Plotting for animation
    if anim == 1
        if mod(n, n_anim) == 0
            figure(2)
            tiledlayout(1, 2)
            nexttile
            U = reshape(u_np1, [Ny - 1, Nx - 1]);
            mesh(Xax, Yax, U)

            title(strcat('Time $t=', num2str(n * k, '%.4f'), '\times', time_scale, '$'), 'Fontsize', 20, 'Interpreter', 'latex')
            xlabel('$x[m]$', 'Fontsize', 20, 'Interpreter', 'latex')
            ylabel('$y[m]$', 'Fontsize', 20, 'Interpreter', 'latex')
            zlabel('$u(x,y,t)$', 'Fontsize', 20, 'Interpreter', 'latex')
            zlim([-limz * 1.5, limz * 1.5])
            
            nexttile
            hold on
            plot([1:Nx - 1] * Di_ext * h, U(round(Ny / 2), 1:Nx - 1), 'o-', 'DisplayName', strcat('$u(x,y=', num2str(0.5 * Di_ext, '%.2e'), ',t)$'))
            plot([1:Nx - 1] * Di_ext * h, K_Mat(floor(Ny / 2), 1:Nx - 1) / K * limz * 1.5, '*-', 'DisplayName', strcat('$\kappa(x,y=', num2str(0.5 * Di_ext, '%.2e'), ',t)$(scaled)'))
            xlabel('$x[m]$', 'Fontsize', 20, 'Interpreter', 'latex')
            legend('Fontsize', 17, 'Interpreter', 'latex', 'Location', 'southeast')
            ylim([-limz * 1.5, limz * 1.5])
            title('Centred side view', 'Fontsize', 20, 'Interpreter', 'latex')
            hold off
            
            x0 = 250;
            y0 = 50;
            width = 1000;
            height = 500;
            set(gcf, 'position', [x0, y0, width, height])
            drawnow limitrate
        end
    end

    u_nm1 = u_n;
    u_n = u_np1;

    % Record output for readout positions
    for i = 1:size(rp_mat, 1)
        out(n, i) = u_np1(rp_index_vect(i));
    end
end

% End of the main loop

% Plotting simulation output if enabled
if ploting == 1
    t = [1:NF] * k;
    figure(3)
    hold on
    for i = 1:size(rp_mat, 1)
        leyenda = strcat('$(x,y)=(', num2str(rp_mat(i, 1) * Di_ext, '%.2e'), ', ', num2str(rp_mat(i, 2) * Di_ext, '%.2e'), ')$');
        plot(t, out(:, i), 'DisplayName', leyenda)
    end
    lgd = legend('Fontsize', 17, 'Interpreter', 'latex');
    title(lgd, 'Readout position');
    ylabel('Amplitude $[m]$', 'Fontsize', 20, 'Interpreter', 'latex')
    xlabel(strcat('Time $[\sigma]=[', time_scale, ']$'), 'Fontsize', 20, 'Interpreter', 'latex')
    
    x0 = 10;
    y0 = 500;
    width = 1300;
    height = 500;
    set(gcf, 'position', [x0, y0, width, height])
end

% Plot Fast Fourier Transform (FFT) if enabled
if plot_fft == 1
    figure(4)
    tiledlayout(size(rp_mat, 1), 1)
    
    for i = 1:size(rp_mat, 1)
        Y = fft(out(:, i));
        P2 = abs(Y / NF);
        P1 = P2(1:NF / 2 + 1);
        P1(2:end - 1) = 2 * P1(2:end - 1);
        f = 1e-6 / (sqrt(E_scale_param)) * SR * (0:(NF / 2)) / NF;  % if [E]˜=E_scale_param
        [M, In] = max(P1);
        
        nexttile
        
        findpeaks(P1(1:floor(end / 5)), f(1:floor(end / 5)), 'MinPeakDistance', 2, ...
            'MinPeakHeight', M / 10, 'Annotate', 'extents')
        set(gca, 'FontSize', 20);
        xlabel('Frequency (MHz)', 'FontSize', 20, 'Interpreter', 'latex');
        ylabel('$|\hat{u}|$', 'FontSize', 20, 'Interpreter', 'latex');
        xlim([0, f(In) + 100])
        read_out_point = strcat('$(x,y)=(', num2str(rp_mat(i, 1) * Di_ext, '%.2e'), ', ', num2str(rp_mat(i, 2) * Di_ext, '%.2e'), ')$');
        freq_res = strcat('Fourier Transform $f_r$=', num2str(f(In)), 'MHz');
        tit = [freq_res newline read_out_point];
        title(tit, 'FontSize', 20, 'Interpreter', 'latex');
    end
    
    x0 = 500;
    y0 = 500;
    width = 800;
    height = 1300;
    set(gcf, 'position', [x0, y0, width, height])
end

% End measuring execution time
toc;
end
