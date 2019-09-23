function [fout, faux] = network_subtemp_slab_timedep(v_in,parameters)
%finite-volume, steady solver of the Stokes problem in the stream function
%(psi)-vorticity(omega) formulation over a rectangular domain with T-dep
%sliding. Heat conservation is an advection-diffusion problem solved in a
%boundary layer near the bed. Same boundary layer treatment for the heat
%equation (diffusion only) in the bed. Boundary conditions are periodic in
%the along flow direction

%The equations being solved are:

%div(grad psi) = omega, div(grad omega) = 0,

%with stress free ice surface (omega = 0 on z = 1), velocity strengthening friction law (omega = f(psi_z) on z=0) with two different friction coefficients, impermeability of the
%bed (psi = 0 on z=0), and constant flux (psi = 1 on z=1)

%dT/dt +u dT/dx + w dT/dz -d^2T/dz^2 = 0 on z>0
%dT/dt -d^2T/dz^2 = 0 on z<0
% 
% with far field flux in the boundary layer (z->inf) prescribed by the Q-equation
% dQ/dt +u dQ/dx -Q du/dx = 0,
% 
% [dT/dz]^+_- + alpha taub u_b = 0 on z= 0
% 
% dT/dz -> - nu on z-> -inf

%Jacobian is in network_subtemp_slab_jacobian.m 
%
% Input variables are:
% v_in:          concatenated vector of size (????,1) of the form v_in =
% [psi;omega; T_ice;T_bed; Q;T_bdy_bed];
% parameters:    Parameter structure with the following fields

% grid:
%               substructure with fields
%               n_nodes: substructure with fields 'tot' (number of nodes in network), 'ver' (number of nodes for one ice column). 
%               n_edges: substructure with fields 'hor' and 'vert'. Number of network edges
%               up_node: substructure with fields 'hor' and 'vert'. List of upstream nodes for each edge
%               down_node: substructure with fields 'hor' and 'vert'. List of downstream nodes for each edge
%                   (up_node and down_node together give the data expected
%                   in a connectivity array)
%               bdy_nodes: substructure with fields 'flux', 'dir'. Indexes
%               of nodes where bdy conditions apply
%               bed_nodes
%               coor_nodes: substructure with fields 'sigma' and 'eta'. List of
%               coordinates of nodes
%               id_node: list of flags (1-2-3) identifying the subdomain
%               Delta_z: scalar, vertical spacing of nodes
%               Delta_x_ scalar, hor spacing of nodes

% n_x:          half the number of nodes in the hor direction

% alpha:        strength of strain heating
% nu:           geothermal heating
% Pe:           Peclet number

%Tested against numerical jacobian. Elisa Mantelli, 13 Nov 2018

%unpack parameters
%ice thickness scale grid
psi_nodes = parameters.grid.psi.n_nodes.tot;                                   %number of nodes
psi_up_node_ver = parameters.grid.psi.up_node.vert;                            %list (n_edges-by-1 vector) of 'upstream' node for each vertical edge
psi_down_node_ver = parameters.grid.psi.down_node.vert;                        %list (n_edges-by-1 vector) of 'downstream' nodes for each vertical edge
psi_up_node_hor = parameters.grid.psi.up_node.hor;                             %list (n_edges-by-1 vector) of 'upstream' node for each horizontal edge
psi_down_node_hor = parameters.grid.psi.down_node.hor;                         %list (n_edges-by-1 vector) of 'downstream' nodes for each horizontal edge
psi_bdy_nodes_top = parameters.grid.psi.bdy_nodes.top;                         %list of nodes adjacent to the top of the box;
psi_bdy_nodes_bed = parameters.grid.psi.bdy_nodes.bed;
psi_bdy_nodes_inflow = parameters.grid.psi.bdy_nodes.inflow;
psi_bdy_nodes_outflow = parameters.grid.psi.bdy_nodes.outflow;

psi_Delta_z_cell = parameters.grid.psi.Delta_z_cell;                           %length of cells, ver (list)
psi_Delta_x_cell = parameters.grid.psi.Delta_x_cell;                           %length of cells, hor (list)

%boundary layer grid
T_nodes = parameters.grid.T.ice.n_nodes.tot;                                   %number of nodes
T_up_node_ver = parameters.grid.T.ice.up_node.vert;                            %list (n_edges-by-1 vector) of 'upstream' node for each vertical edge
T_down_node_ver = parameters.grid.T.ice.down_node.vert;                        %list (n_edges-by-1 vector) of 'downstream' nodes for each vertical edge
T_up_node_hor = parameters.grid.T.ice.up_node.hor;                             %list (n_edges-by-1 vector) of 'upstream' node for each horizontal edge
T_down_node_hor = parameters.grid.T.ice.down_node.hor;                         %list (n_edges-by-1 vector) of 'downstream' nodes for each horizontal edge
T_bdy_nodes_top = parameters.grid.T.ice.bdy_nodes.top;                         %list of nodes adjacent to the top of the box;
T_bdy_nodes_bed = parameters.grid.T.ice.bdy_nodes.bed;

T_Delta_z_edge = parameters.grid.T.ice.Delta_z_edge;                             %length of cells, ver (list)
T_Delta_z_cell = parameters.grid.T.ice.Delta_z_cell;                             %length of cells, ver (list)
T_Delta_x_cell = parameters.grid.T.ice.Delta_x_cell;                             %length of cells, hor (list)

index_bed_to_T_hor_edges = parameters.grid.T.ice.index_bed_to_T_hor_edges;
index_bed_to_T_ver_edges = parameters.grid.T.ice.index_bed_to_T_ver_nodes;

T_coor_ver_edges_z = parameters.grid.T.ice.coor_veredges.z;
T_coor_ver_edges_z_top = parameters.grid.T.ice.coord_veredges_top.z;

%T bed grid
Tb_up_node_ver = parameters.grid.T.bed.up_node.vert;                            %list (n_edges-by-1 vector) of 'upstream' node for each vertical edge
Tb_down_node_ver = parameters.grid.T.bed.down_node.vert;                        %list (n_edges-by-1 vector) of 'downstream' nodes for each vertical edge
Tb_bdy_nodes_top = parameters.grid.T.bed.bdy_nodes.top;                         %list of nodes adjacent to the top of the box;
Tb_bdy_nodes_bed = parameters.grid.T.bed.bdy_nodes.bed;
Tb_Delta_z_cell = parameters.grid.T.bed.Delta_z_cell; 

%1D grid
Q_nodes = parameters.grid.Q.n_nodes.tot;                                   %number of nodes
Q_up_node = parameters.grid.Q.up_node.hor;                                 %list (n_edges-by-1 vector) of 'upstream' node for each horizontal edge
Q_down_node = parameters.grid.Q.down_node.hor;                             %list (n_edges-by-1 vector) of 'downstream' nodes for each horizontal edge
Q_Delta_x_cell = parameters.grid.Q.Delta_x_cell;                           %length of cells, hor (list)

%physical parameters
gamma = parameters.gamma;
alpha = parameters.alpha;
nu = parameters.nu;
dt = parameters.dt;

%unpack input variable v_in 
psi = v_in(1:psi_nodes);
omega = v_in(psi_nodes+1:2*psi_nodes);
T = v_in(2*psi_nodes+1:2*psi_nodes+T_nodes);
Tb = v_in(2*psi_nodes+T_nodes+1:2*psi_nodes+2*T_nodes);
Q = v_in(2*psi_nodes+2*T_nodes+1:2*psi_nodes+2*T_nodes+Q_nodes);
T_bed = v_in(2*psi_nodes+2*T_nodes+Q_nodes+1:end);

%unpack input variable at previous time step
v_in_prev = parameters.v_in_prev;
T_prev = v_in_prev(2*psi_nodes+1:2*psi_nodes+T_nodes);
Tb_prev = v_in_prev(2*psi_nodes+T_nodes+1:2*psi_nodes+2*T_nodes);
Q_prev = v_in_prev(2*psi_nodes+2*T_nodes+1:2*psi_nodes+2*T_nodes+Q_nodes);
%initialize output
fout = zeros(length(v_in),1);
%%  PRELIMINARIES 

%DISCRETE DEP. VARIABLES
[discvar, ~] = discretisation_v3(v_in, parameters);

dpsi_dx = discvar.dpsi_dx;
dpsi_dz = discvar.dpsi_dz;
% flux_psi_in = v_in(2*psinodes+1:2*nodes+length(bdy_nodes_inflow));
% flux_omega_in = v_in(2*nodes+length(bdy_nodes_inflow)+1:2*nodes+2*length(bdy_nodes_inflow));

domega_dx = discvar.domega_dx;
domega_dz = discvar.domega_dz;

T_hor = discvar.T.T_hor;
T_vert = discvar.T.T_vert;
T_top = discvar.T.T_top;

dTb_dx = discvar.Tb.dT_dx;
dTb_dz = discvar.Tb.dT_dz;

Q_hor = discvar.Q.Q_hor;

T_bed_psigrid = discvar.T_bed.Tpsigrid;

%construct regularized bedwater content
[f_slide_Tbed, ~] = regularization_old(T_bed, parameters);
[f_slide_Tbedpsigrid, ~] = regularization_old(T_bed_psigrid, parameters);
%% STREAM FUNCTION
net_psi_horflux = (accumarray(psi_up_node_hor, dpsi_dx, [psi_nodes,1])-accumarray(psi_down_node_hor, dpsi_dx, [psi_nodes,1]));
net_psi_verflux = (accumarray(psi_up_node_ver, dpsi_dz, [psi_nodes,1])-accumarray(psi_down_node_ver, dpsi_dz, [psi_nodes,1]));

%enforce Dirichlet conditions on top and bottom boundaries (psi_x=0 on
%inflow and outflow bdies)
%ice surface: psi = 1
%net_psi_verflux(psi_bdy_nodes_top) = net_psi_verflux(psi_bdy_nodes_top) + (1-psi(psi_bdy_nodes_top))./(psi_Delta_z_cell(psi_bdy_nodes_top)/2);
psitop = 1;
dpsidztop = (-9*psi(psi_bdy_nodes_top) + psi(psi_bdy_nodes_top + length(psi_bdy_nodes_top)) +8*psitop)./(3*psi_Delta_z_cell(psi_bdy_nodes_top));
net_psi_verflux(psi_bdy_nodes_top) = net_psi_verflux(psi_bdy_nodes_top) + dpsidztop;

%bed: psi = 0
%net_psi_verflux(psi_bdy_nodes_bed) = net_psi_verflux(psi_bdy_nodes_bed) - (psi(psi_bdy_nodes_bed))./(psi_Delta_z_cell(psi_bdy_nodes_bed)/2);
net_psi_verflux(psi_bdy_nodes_bed) = net_psi_verflux(psi_bdy_nodes_bed) -(9*psi(psi_bdy_nodes_bed) - psi(psi_bdy_nodes_bed - length(psi_bdy_nodes_bed)))./(3*psi_Delta_z_cell(psi_bdy_nodes_bed));
%periodic boundary conditions
flux_psi_in = (psi(psi_bdy_nodes_inflow) -psi(psi_bdy_nodes_outflow))./(psi_Delta_x_cell(psi_bdy_nodes_inflow)/2 + psi_Delta_x_cell(psi_bdy_nodes_outflow)/2);
%inflow boundary
net_psi_horflux(psi_bdy_nodes_inflow) = net_psi_horflux(psi_bdy_nodes_inflow) - flux_psi_in;
%outflow boundary
net_psi_horflux(psi_bdy_nodes_outflow) = net_psi_horflux(psi_bdy_nodes_outflow) + flux_psi_in;

%conservation law
div_psifluxes = net_psi_horflux./psi_Delta_x_cell  +  net_psi_verflux./psi_Delta_z_cell;
fout(1:psi_nodes) = div_psifluxes - omega;

%% SLIDING LAW
%at psi cell centres
u_bed = (9*psi(psi_bdy_nodes_bed) - psi(psi_bdy_nodes_bed - length(psi_bdy_nodes_bed)))./(3*psi_Delta_z_cell(psi_bdy_nodes_bed)); 
tau_bed = gamma*(u_bed)./f_slide_Tbedpsigrid; 

%at T cell centres
%derivative of sliding velocity
u_bed_full = [u_bed; u_bed(1)];
dubed_dx = (u_bed_full(2:end) - u_bed_full(1:end-1))./Q_Delta_x_cell;

%sliding velocity and stress
u_bed_Tcentre = (u_bed_full(1:end-1)+u_bed_full(2:end))./2;
tau_bed_Tcentre = gamma*u_bed_Tcentre./f_slide_Tbed;
%% VORTICITY
net_omega_horflux = (accumarray(psi_up_node_hor, domega_dx, [psi_nodes,1])-accumarray(psi_down_node_hor, domega_dx, [psi_nodes,1]));
net_omega_verflux = (accumarray(psi_up_node_ver, domega_dz, [psi_nodes,1])-accumarray(psi_down_node_ver, domega_dz, [psi_nodes,1]));

%enforce boundary conditions
%ice surface: omega = 0
net_omega_verflux(psi_bdy_nodes_top) = net_omega_verflux(psi_bdy_nodes_top) + (-omega(psi_bdy_nodes_top))./(psi_Delta_z_cell(psi_bdy_nodes_top)/2);
%bed: omega = omega_bed
%domegadz_bed = (omega(psi_bdy_nodes_bed)-tau_bed)./(psi_Delta_z_cell(psi_bdy_nodes_bed)/2);
domegadz_bed =  (9*omega(psi_bdy_nodes_bed) - omega(psi_bdy_nodes_bed - length(psi_bdy_nodes_bed))-8*tau_bed)./(3*psi_Delta_z_cell(psi_bdy_nodes_bed));

net_omega_verflux(psi_bdy_nodes_bed) = net_omega_verflux(psi_bdy_nodes_bed) - domegadz_bed;

%periodic boundary conditions
flux_omega_in = (omega(psi_bdy_nodes_inflow) -omega(psi_bdy_nodes_outflow))./(psi_Delta_x_cell(psi_bdy_nodes_inflow)/2 + psi_Delta_x_cell(psi_bdy_nodes_outflow)/2);
%inflow boundary
net_omega_horflux(psi_bdy_nodes_inflow) = net_omega_horflux(psi_bdy_nodes_inflow) - flux_omega_in;
%outflow boundary
net_omega_horflux(psi_bdy_nodes_outflow) = net_omega_horflux(psi_bdy_nodes_outflow) + flux_omega_in;

%conservation law
fout(psi_nodes+1:2*psi_nodes) = net_omega_horflux./psi_Delta_x_cell  +  net_omega_verflux./psi_Delta_z_cell;
%inferred tau_b 
% taub_balance = -(fout(psi_nodes + psi_bdy_nodes_bed) .* psi_Delta_z_cell(psi_bdy_nodes_bed)/2 - omega(psi_bdy_nodes_bed));

%% VELOCITY FIELD IN THE BOUNDARY LAYER
%horizontal velocity at T hor edges
u_bed_Tedges = u_bed;
u_bed_edges_Tgrid = u_bed_Tedges(index_bed_to_T_hor_edges);

%vertical velocity at T cell centres w = - z du_bed/dx
W = -dubed_dx(index_bed_to_T_ver_edges).*T_coor_ver_edges_z;
W_top = -dubed_dx.*T_coor_ver_edges_z_top;
%% Q EQUATION
% solve dQ/dt + d/dx(Q *u) -2Qdu/dx = 0 on 1D grid with periodic boundary
% conditions
Qbed_flux = Q_hor.*(u_bed);
net_Qbedflux = (accumarray(Q_up_node, Qbed_flux, [Q_nodes,1])-accumarray(Q_down_node, Qbed_flux, [Q_nodes,1]));

%source term
source_Q = 2*Q.*dubed_dx;
fout(2*psi_nodes+2*T_nodes+1:2*psi_nodes+2*T_nodes+Q_nodes) =  (Q-Q_prev)/dt +net_Qbedflux./Q_Delta_x_cell-source_Q;%
%% HEAT EQUATION ICE
%solve dT/dt + d/dX((u_bed)T) + d/dZ[-Z du_b/dX T - dT/dZ] = 0
%nb: T grid is shifted with respect to psi grid, so T edges correspond to
%psi centres

%HORIZONTAL HEAT FLUX
QHa = (u_bed_edges_Tgrid).*T_hor; 
net_horflux = (accumarray(T_up_node_hor, QHa, [T_nodes,1])-accumarray(T_down_node_hor, QHa, [T_nodes,1]));

%VERTICAL HEAT FLUX 
QVa = W.*T_vert; 
net_a_verflux = (accumarray(T_up_node_ver, QVa, [T_nodes,1])-accumarray(T_down_node_ver, QVa, [T_nodes,1]));
%correct for advective mass flux at the top boundary (assume a row of ghost
%cells at prescribed temperature T_top
net_a_verflux(T_bdy_nodes_top) =  net_a_verflux(T_bdy_nodes_top) + (W_top.*(T_top+ T(T_bdy_nodes_top))./2);

%diffusive flux
Qvd = - (T(T_down_node_ver)- T(T_up_node_ver))./T_Delta_z_edge;
net_d_verflux = (accumarray(T_up_node_ver, Qvd, [T_nodes,1])-accumarray(T_down_node_ver, Qvd, [T_nodes,1]));

%basal energy budget
bedflux_bl_ice = - (T(T_bdy_nodes_bed)-T_bed)./(T_Delta_z_cell(T_bdy_nodes_bed)/2);   
net_d_verflux(T_bdy_nodes_bed) = net_d_verflux(T_bdy_nodes_bed) - bedflux_bl_ice;
% %Neumann condition at the top of the box
net_d_verflux(T_bdy_nodes_top) =  net_d_verflux(T_bdy_nodes_top) + Q;

%sum advective and diffusive fluxes
net_verflux = net_d_verflux + net_a_verflux;

%enforce conservation law 
fout(2*psi_nodes+1:2*psi_nodes+T_nodes) =  (T-T_prev)/dt  + net_verflux./T_Delta_z_cell + net_horflux./T_Delta_x_cell;
%% HEAT EQUATION BED

%VERTICAL HEAT FLUX
QVdbed = -dTb_dz;
net_bedverflux = (accumarray(Tb_up_node_ver, QVdbed, [T_nodes,1])-accumarray(Tb_down_node_ver, QVdbed, [T_nodes,1]));

%correct for geothermal heat flux at the bottom
net_bedverflux(Tb_bdy_nodes_bed) = net_bedverflux(Tb_bdy_nodes_bed) - nu;
%correct for bed flux
bedflux_bl_bed = - (T_bed-Tb(Tb_bdy_nodes_top))./(Tb_Delta_z_cell(Tb_bdy_nodes_top)/2);   
net_bedverflux(Tb_bdy_nodes_top) = net_bedverflux(Tb_bdy_nodes_top) + bedflux_bl_bed;

%enforce conservation law
fout(2*psi_nodes+T_nodes+1:2*psi_nodes+2*T_nodes) = (Tb-Tb_prev)/dt + net_bedverflux./Tb_Delta_z_cell ;%+net_bedhorflux./Tb_Delta_x_cell; 
%% BASAL ENERGY BUDGET
%this is the boundary condition for basal temperature
m = -bedflux_bl_ice + bedflux_bl_bed +alpha*tau_bed_Tcentre.*u_bed_Tcentre;

%conservation law
fout(2*psi_nodes+2*T_nodes+Q_nodes+1:2*psi_nodes+2*T_nodes+2*Q_nodes) = m;
%% AUXILIARY VARIABLES
faux.bedflux = bedflux_bl_ice;
faux.heating = alpha*tau_bed_Tcentre.*u_bed_Tcentre;
faux.u_bed = u_bed;
faux.tau_bed = tau_bed;
faux.W = W;
faux.Q =Q;
faux.vorticityflux = domegadz_bed;

































