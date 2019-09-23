function [fout, faux] = network_sstate_v4(v_in,parameters)
%Computes functions whose zeros define the solution to the steady state problem for a flowline of a marine ice sheet
%with basal thermal transitions.

%Input variables are:
%v_in:          Concatenated column vectors for H, T, xc, xs, xt

%parameters:    Parameter structure with the following fields
% grid_h:       substructure with fields
%               n_nodes: number of nodes in network
%               n_edges: number of network edges
%               up_node: list of upstream nodes for each edge
%               down_node: list of downstream nodes for each edge
%                   (up_node and down_node together give the data expected
%                   in a connectivity array)
%               coor_nodes: List of coordinates of nodes
%               coor_edges: List of coordinates of edges
%               id_node: list of flags (1-2-3) identifying the subdomain
%               subtemp_slid: list of subtemperate sliding nodes

% grid_T, grid_u:
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
%               Delta_eta: scalar, vertical spacing of nodes

% n_x:          substructure with fields 'c','s','t'. Number of h,T nodes along
%               one row in each subdomain

%flag:          substructure with fields 'discr'(steady/unsteady- the latter not yet implemented), setting 
%               the upwinding scheme for the horizontal heat flux. 'u_sl' (flag setting the
%               averaging scheme for the vertical heat flux in the computation of the
%               sliding velocity. 'upwind'-> heat flux is upwinded like T_hor; 'flux_av'-> flux along hor 
%               edge is averaged between adjacent cells; 'T_av'-> T is averaged between adjacent cells )

% alpha:        strength of strain heating
% gamma:        friction parameter
% nu:           geothermal heating
% Pe:           Peclet number
% a:            accumulation rat

%bed:           substructure with fields
%               b0: depth to bed at the divide (<0)
%               b1: slope (>0)

%gr:            substructure with parameters for flux bc at the grounding
%               line. subfields k and m for the relationship between h and
%               q at the grounding line, Q_g=k*H_g^m

%
%The vectors H and T in v_in have dimensions h_nodes-by-one,T_nodes-by-one
%respectively, xc,xs and xt are one-by-one, and v_in = [H; T; xc; xs; xc].
% H and T are defined at nodes.
%
%Physically, H is ice thickness, T is temperature, and x_i are the
%boundaries between cold-subtemperate, sub-temperate and temperate
%subdomains, and the grounding line.
%
%Conservation of mass at each node gives (constraint #1)
% d(storage at node) / dt + 1/2* gamma * sum (dS_R/dt+(n_c-1)*dS_K/dt)-  ( q_in +
% sum Q_in - sum Q_out ) = 0
%where Q_in is  discharges along an edge into the node, and Q_out is
%discharge along an edge out of the node, and the sum over dS_R/dt nd
%S_K/dt is taken over all edges into or out of the node. Storage at a node
%has an 'elastic' componennt gamma_store*(N-N_0), where N_0 is a local
%constant.
%

%Conservation of heat at each node gives (constraint #2)
% d(storage at node) / dt + 1/2* gamma * sum (dS_R/dt+(n_c-1)*dS_K/dt)-  ( q_in +
% sum Q_in - sum Q_out ) = 0
%where Q_in is  discharges along an edge into the node, and Q_out is
%discharge along an edge out of the node, and the sum over dS_R/dt nd
%S_K/dt is taken over all edges into or out of the node. Storage at a node
%has an 'elastic' componennt gamma_store*(N-N_0), where N_0 is a local
%constant.
%

%Boundary conditions for the free boundaries are
% (constraint #3)

%The output of this function is a concatenated list of the left-hand
%sides of the constraints given
%above evaluated at each node.  The order is: mass conservation; heat
%conservation/Dirichlet conditions for temperature; 3 b.c. for free
%boundaries. 
%
%Elisa Mantelli, Sept. 2019


%unpack parameters
h_nodes = parameters.grid_h.n_nodes;                  %number of nodes,h
h_up_node = parameters.grid_h.up_node;                %list (h_n_edge-by-1 vector) of 'upstream' node, h
h_down_node = parameters.grid_h.down_node;            %list (h_n_edges-by-1 vector) of 'downstream' node,h
h_n_edges = parameters.grid_h.n_edges;                %number of edges in the h network 


T_nodes = parameters.grid_T.n_nodes.tot;                                   %number of temperature nodes
nodes_ver = parameters.grid_T.n_nodes.vert;                                %number of nodes in the vertical direction
T_n_edges_hor = parameters.grid_T.n_edges.hor;        %number of horizontal edges, T
T_up_node_ver = parameters.grid_T.up_node.vert;                            %list (n_edges-by-1 vector) of 'upstream' node for each vertical edge,T
T_down_node_ver = parameters.grid_T.down_node.vert;                        %list (n_edges-by-1 vector) of 'downstream' nodes for each vertical edge,T
T_up_node_hor = parameters.grid_T.up_node.hor;                             %list (n_edges-by-1 vector) of 'upstream' node for each horizontal edge,T
T_down_node_hor = parameters.grid_T.down_node.hor;                         %list (n_edges-by-1 vector) of 'downstream' nodes for each horizontal edge,T
T_Delta_eta = parameters.grid_T.Delta_eta;                                 %scalar, vertical spacing between cell centres
T_bdy_nodes_flux = parameters.grid_T.bdy_nodes.flux;                       %list (n_bdy_nodes-by-1 vector) of nodes at domain boundaries where Neumann conditions
T_bdy_nodes_dir_surf = parameters.grid_T.bdy_nodes.dir.surf;               %list  of nodes at the upper surface where Dirichlet conditions apply
T_bdy_nodes_dir_bed = parameters.grid_T.bdy_nodes.dir.bed;                 %list  of nodes at the bed where Dirichlet conditions apply
T_bdy_nodes_inflow = parameters.grid_T.bdy_nodes.inflow;                     
T_bdy_nodes_subtemp = parameters.grid_T.bdy_nodes.subtemp_slid;            %list of basal T nodes in the subtemperate subdomain
T_coor_nodes_eta = parameters.grid_T.coor_nodes.eta;

alpha = parameters.alpha;   % strain heating parameter
gamma = parameters.gamma;   % friction coefficient
Pe = parameters.Pe;         % Peclet number
nu = parameters.nu;         % basal heating parameter
a = parameters.a;           %accumulation rate
r_rho = parameters.r_rho;   %water to ice density ratio

T_surf = parameters.T_surf; %surface temperature

flag_plot = parameters.flag_plot;

%unpack input variable v_in
H = v_in(1:h_nodes);
T = v_in(h_nodes+1:h_nodes+T_nodes);
xc = v_in(h_nodes+T_nodes+1);
xs = v_in(h_nodes+T_nodes+2);
xt = v_in(h_nodes+T_nodes+3);

%initialize fout
fout = zeros(h_nodes+T_nodes+3,1);
%%  PRELIMINARIES 

%horizontal length of subdomains, cells and edges 
length_h = fvlength([xc;xs;xt], parameters,'h');
length_T = fvlength([xc;xs;xt], parameters,'T');
%unpack
L_h = length_h.L;
Delta_sigma_h = length_h. Delta_sigma;
L_T = length_T.L;
Delta_sigma_T = length_T. Delta_sigma;

%BED
st_bed = bed([xc;xs;xt], parameters);
%unpack: depth to bed at grounding line
B_g = st_bed.B.g;

%DISCRETE DEP. VARIABLES
[discvar, ~] = discretisation_v2(v_in, parameters);

%ICE THICKNESS at the grounding line
H_g = discvar.H.H_g;

%TEMPERATURE 
%along T edges 
T_hor = discvar.T.T_hor;
T_vert = discvar.T.T_vert;
%along hor edges at the grounding line
T_hor_g = discvar.T.T_hor_g;

%FLUX
QH_v_edge = discvar.flux.QH_v_edge;

%SHEARING VELOCITY FIELD
[hor_velocity, ~] = hor_velocity_SIA_v2(v_in, parameters);
QMd = hor_velocity.QMd;
QHa_shear = hor_velocity.QHa_shear;
QHa_shear_g = hor_velocity.QHa_shear_g;
S = hor_velocity.S;
net_flux_shear_eta_sum = hor_velocity.net_flux_shear_eta_sum;
%grounding line
QMd_g = hor_velocity.QMd_g;

%SLIDING VELOCITY
% compute sliding velocity at the grounding line
flux = fluxg(H_g,parameters);
QMg = flux.flux;
U_SL_g = (QMg-QMd_g)/H(end);
QHha_slide_g = U_SL_g*H_g.*T_hor_g;
%along h edges
[usl, ~] = sliding_v3(v_in, parameters);
QMa = usl.QM_a;
QHha_slide = usl.QHh_a;
U_SL_temp_st = usl.U_SL_temp_st ; %temperate sliding velocity at x = xs
%% MASS CONSERVATION
%sum advective and diffusive mass flux along each edge
QM = QMd + QMa;
%net mass flux out of each node
netmassflux = (accumarray(h_up_node, QM, [h_nodes,1])-accumarray(h_down_node, QM, [h_nodes,1]));
% correct for grounding line outflow
netmassflux(end) = netmassflux(end) + QMd_g + U_SL_g*H(end);
%impose mass conservation at each node
fout(1:h_nodes) = netmassflux./Delta_sigma_h - a*L_h;
%% HEAT CONSERVATION

%HORIZONTAL HEAT FLUX (only advection)
QHh = Pe*(QHha_slide+QHa_shear); 
net_H_horflux = (accumarray(T_up_node_hor, QHh, [T_nodes,1])-accumarray(T_down_node_hor, QHh, [T_nodes,1]));
%correct for advective mass flux at the grounding line
net_H_horflux(end-nodes_ver+1:end) =  net_H_horflux(end-nodes_ver+1:end)+ Pe*(QHha_slide_g + QHa_shear_g);
%VERTICAL HEAT FLUX along vertical edges + flux bc at the bed
%diffusive component
[QHvd,~] = diffver_heatflux(v_in, parameters);
QH_v_d = QHvd.QH_v_d;
%advective component
netmassflux_shear = (accumarray(h_up_node, QMd, [h_nodes,1])-accumarray(h_down_node, QMd, [h_nodes,1]));
netmassflux_shear(end) = netmassflux_shear(end)+ QMd_g;
W_eff =T_Delta_eta.*reshape((1:nodes_ver)'*(netmassflux_shear./(L_h.*Delta_sigma_h) - a).',[T_nodes,1])-(net_flux_shear_eta_sum./L_T);
index_u_to_T = repmat([ones(nodes_ver-1,1); 0], [h_nodes,1]);%reindex from u to T ver network
QH_v_a = L_T(index_u_to_T==1).*W_eff(index_u_to_T==1).*T_vert;
QH_v_a_surf = L_T(index_u_to_T==0).*W_eff(index_u_to_T==0).*T_surf;

%ver flux
QH_v = (Pe)* QH_v_a +QH_v_d;
net_H_verflux = (accumarray(T_up_node_ver, QH_v, [T_nodes,1])-accumarray(T_down_node_ver, QH_v, [T_nodes,1]));

%correct for geothermal heat flux at cold bed cells
net_H_verflux(T_bdy_nodes_flux) = net_H_verflux(T_bdy_nodes_flux)- nu*L_T(T_bdy_nodes_flux);

%correct ver flux in order to account for basal temperature at the melting
%point in subtemperate and temperate subdomain
net_H_verflux(T_bdy_nodes_dir_bed) = net_H_verflux(T_bdy_nodes_dir_bed)- ( - L_T(T_bdy_nodes_dir_bed).*T(T_bdy_nodes_dir_bed)./(H((T_bdy_nodes_dir_bed-1)/nodes_ver+1).*T_Delta_eta/2))  ;

%correct for ver flux in surface cells to account for T_surf = -1. n.b.:
%both advective and conductive contributions matter
net_H_verflux(T_bdy_nodes_dir_surf) = net_H_verflux(T_bdy_nodes_dir_surf) + ( - L_T(T_bdy_nodes_dir_surf).*(-1 - T(T_bdy_nodes_dir_surf))./(H((T_bdy_nodes_dir_surf)/nodes_ver).*T_Delta_eta/2))...
    + Pe* QH_v_a_surf;

%ENFORCE HEAT CONS at each node
fout(h_nodes+1:h_nodes+T_nodes) = net_H_horflux./Delta_sigma_T  +  net_H_verflux./T_Delta_eta - alpha*L_T.*S;

%% FREE BOUNDARIES 
% flux as computed at subtemperate nodes equals geothermal heat fluxfout(h_nodes+T_nodes+1) = QH_v_edge(parameters.n_x.c) - nu;  

% basal flux such that first temperate edge satisfies basal energy budget=0;
fout(h_nodes+T_nodes+2) =  -(QH_v_edge(parameters.n_x.c+parameters.n_x.s)) + alpha*gamma*(U_SL_temp_st).^2 + nu;

% flotation at the grounding line
fout(h_nodes+T_nodes+3) = H_g - r_rho*B_g;

%computation of the vertical velocity

index_h_to_T = reshape(repmat((1:h_n_edges),[nodes_ver,1]),[T_n_edges_hor,1]);
%the quantities below live on T hor edges, whereas w_eff lives on T ver
%edges
ub = usl.U_sl_Thoredges;
ushear = hor_velocity.U_Thoredges;
Uedges = ub+ushear;
dbdx = st_bed.dB.edges(index_h_to_T);
ds =  discvar.dH_dx;
eta_Tedges_inflow = T_coor_nodes_eta(T_bdy_nodes_inflow);
eta_ds_matrix = eta_Tedges_inflow*ds.';
eta_ds = reshape(eta_ds_matrix,[T_n_edges_hor,1]);
u_stretch = Uedges.*(-dbdx +eta_ds);
%construct W_eff at T nodes
W_eff_surf = W_eff(index_u_to_T==0);
W_eff_bed = zeros(length(W_eff_surf),1);
W_eff_body = W_eff(index_u_to_T==1);

W_eff_nodes =  1/2*(accumarray(T_up_node_ver, W_eff_body, [T_nodes,1])+accumarray(T_down_node_ver, W_eff_body, [T_nodes,1]));
W_eff_nodes(T_bdy_nodes_dir_surf) =  W_eff_nodes(T_bdy_nodes_dir_surf)+1/2*W_eff_surf;
W_eff_nodes([T_bdy_nodes_flux; T_bdy_nodes_dir_bed]) =  W_eff_nodes([T_bdy_nodes_flux; T_bdy_nodes_dir_bed])+1/2*W_eff_bed;

%costruct u at T nodes
u_stretch_nodes =  1/2*(accumarray(T_up_node_hor, u_stretch, [T_nodes,1])+accumarray(T_down_node_hor, u_stretch, [T_nodes,1]));

faux.w_nodes = (W_eff_nodes + u_stretch_nodes);
faux.u_nodes = 1/2*(accumarray(T_up_node_hor, Uedges, [T_nodes,1])+accumarray(T_down_node_hor, Uedges, [T_nodes,1]));
faux.W_eff = W_eff_nodes.*L_T;

faux.u_edges = Uedges.*reshape(repmat(discvar.H.H_edge.',[nodes_ver,1]),[T_n_edges_hor,1]);
%compute w at T ver edges
w_matrix = reshape(faux.w_nodes, [nodes_ver,h_nodes]);
%w_matrix = [w_matrix(1,:);  w_matrix; w_matrix(end,:)];
w_Tvedges_matrix = 1/2*(w_matrix(1:end-1,:)+w_matrix(2:end,:) );
faux.w.body = reshape(w_Tvedges_matrix,[(nodes_ver-1)*h_nodes,1]);
faux.w.surf = w_matrix(end,:).';
faux.w.bed = 1/2*(accumarray(h_up_node, -usl.U_SL.*st_bed.dB.edges, [h_nodes,1])+accumarray(h_down_node, -usl.U_SL.*st_bed.dB.edges, [h_nodes,1]));
faux.u_bed = 1/2*(accumarray(h_up_node, usl.U_SL, [h_nodes,1])+accumarray(h_down_node, usl.U_SL, [h_nodes,1]));
faux.w_eff.body = W_eff_body.*L_T(index_u_to_T==1);
faux.w_eff.bed = W_eff_bed;
faux.w_eff.surf = W_eff_surf.*L_T(index_u_to_T==0);
faux.u_g = H_g.*U_SL_g*ones(nodes_ver,1) + H_g*hor_velocity.U_g;
faux.Q=discvar.flux.QH_v;








