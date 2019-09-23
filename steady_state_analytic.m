function fout = steady_state_analytic(parameters)
%STEP 1: CONSTRUCT ANALYTICAL STEADY STATE .
z_ice = parameters.grid.T.ice.coor_nodes.z(parameters.grid.T.ice.bdy_nodes.inflow);
z_bed = parameters.grid.T.bed.coor_nodes.z(parameters.grid.T.bed.bdy_nodes.inflow);
[f_slide_Tbed, ~] = regularization_old(parameters.T_bed, parameters);
gamma = parameters.gamma;
F0 = f_slide_Tbed;
u0 = 3*F0./(gamma + 3*F0);
tau0 = ( F0 *gamma^-1 ).^(-1) .*u0;
flux_bed = parameters.nu + parameters. alpha.*u0.*tau0;

%construct linear temperature profiles
T_nodes_ver = parameters.grid.T.ice.n_nodes.vert;
T_nodes_hor = parameters.grid.T.ice.n_nodes.hor;

Z_ice = repmat(z_ice,[1,T_nodes_hor]);
Z_bed = repmat(z_bed,[1,T_nodes_hor]);
T_bed_matrix = repmat(parameters.T_bed.',[T_nodes_ver, 1]);
fluxbed_matrix = repmat(flux_bed.',[T_nodes_ver, 1]);
Tb = T_bed_matrix - parameters.nu*Z_bed;
Tice = T_bed_matrix - fluxbed_matrix.*Z_ice;

%construct analytical steady state for psi and omega
z_psi = parameters.grid.psi.coor_nodes.z(parameters.grid.psi.bdy_nodes.inflow);
psi_nodes_hor = parameters.grid.psi.n_nodes.hor;
Z_psi = repmat(z_psi,[1,psi_nodes_hor]);

omega_matrix = tau0(1) - tau0(1)*Z_psi;
psi_matrix = Z_psi + 1/6* (-2 *Z_psi + 3*Z_psi.^2 - Z_psi.^3)* tau0(1);
fout =[reshape(psi_matrix.',[],1); reshape(omega_matrix.',[],1);reshape(Tice.',[],1); reshape(Tb.',[],1); flux_bed; parameters.T_bed];
end