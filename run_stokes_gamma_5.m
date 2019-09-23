N = 8;
maxNumCompThreads(N)
T_bed = -0.5;
choice = 'wn';
gamma = 5;
index_wavelength = 2;



% grid parameters
bd_x = 4*pi;
bd_z.ice = 4;           %top of the boundary layer
bd_z.bed = -4;          %bottom of the bed boundary layer

n_z.psi = 100;
n_z.Tice = 100;
n_z.Tbed = 100;
n_x = 1050;

upwind = 1/2; %set to zero for upwinding hor fluxes in heat equation and Q equation

%timestepping parameters
dt = 7e-03;
t_plot = 10*dt;
stepmax = 2000;
t_init = 0;
step_init = 1;

scale_pert = 1e-03; %amplitude of initial perturbation

network_script_cluster(T_bed, choice, gamma,bd_x,n_x, bd_z, n_z, dt,t_plot, stepmax,t_init, step_init, index_wavelength, scale_pert, upwind)
quit
