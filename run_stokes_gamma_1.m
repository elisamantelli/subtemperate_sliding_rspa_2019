N = 8;
maxNumCompThreads(N)

T_bed = -0.5;
choice = 'mono_small';
gamma = 1;
index_wavelength = 1;

% grid parameters
bd_x = 4*pi;
bd_z.ice = 4;            %top of the boundary layer
bd_z.bed = -4;           %bottom of the bed boundary layer

n_z.psi = 200;
n_z.Tice = 600;
n_z.Tbed = 600;
n_x = 300;

%timestepping parameters
dt = 0.005;
t_plot = 0.1;
stepmax = 1000;
t_init = 0;
step_init = 1;

scale_pert = 1e-03; %amplitude of initial perturbation

network_script_cluster(T_bed, choice, gamma,bd_x,n_x, bd_z, n_z, dt,t_plot, stepmax,t_init, step_init, index_wavelength, scale_pert)

quit
