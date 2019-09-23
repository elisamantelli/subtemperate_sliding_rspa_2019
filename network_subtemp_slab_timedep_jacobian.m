function [fout, faux] = network_subtemp_slab_timedep_jacobian(v_in,parameters)
%Jacobian of network_subtemp_slab_timedep.m
%Tested against numerical jacobian. Elisa Mantelli, 13 Nov 2018

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
Tb_Delta_z_cell = parameters.grid.T.bed.Delta_z_cell;  

%1D grid
Q_nodes = parameters.grid.Q.n_nodes.tot;                                   %number of nodes
Q_up_node = parameters.grid.Q.up_node.hor;                                 %list (n_edges-by-1 vector) of 'upstream' node for each horizontal edge
Q_down_node = parameters.grid.Q.down_node.hor;                             %list (n_edges-by-1 vector) of 'downstream' nodes for each horizontal edge

Q_Delta_x_cell = parameters.grid.Q.Delta_x_cell;                           %length of cells, hor (list)

%physical parameters
gamma = parameters.gamma;
alpha = parameters.alpha;

%timestepping
dt = parameters.dt; 

%unpack input variable v_in 
psi = v_in(1:psi_nodes);
T = v_in(2*psi_nodes+1:2*psi_nodes+T_nodes);
Q = v_in(2*psi_nodes+2*T_nodes+1:2*psi_nodes+2*T_nodes+Q_nodes);
T_bed = v_in(2*psi_nodes+2*T_nodes+Q_nodes+1:end);

%initialize output
fout = sparse(length(v_in),length(v_in));
faux = fout;
%%  PRELIMINARIES 

%DISCRETE DEP. VARIABLES
[discvar, Ddiscvar] = discretisation_v3(v_in, parameters);

ddpsidx_dpsiu = Ddiscvar.ddpsidx_dpsiup;
ddpsidx_dpsid = Ddiscvar.ddpsidx_dpsidown;
ddpsidz_dpsiu = Ddiscvar.ddpsidz_dpsiup;
ddpsidz_dpsid = Ddiscvar.ddpsidz_dpsidown;

ddomegadx_domegau = Ddiscvar.ddomegadx_domegaup;
ddomegadx_domegad = Ddiscvar.ddomegadx_domegadown;
ddomegadz_domegau = Ddiscvar.ddomegadz_domegaup;
ddomegadz_domegad = Ddiscvar.ddomegadz_domegadown;

T_hor = discvar.T.T_hor;
T_vert = discvar.T.T_vert;

dThor_dTup = Ddiscvar.T.dThor_dTup;
dThor_dTdown = Ddiscvar.T.dThor_dTdown;

dTver_dTup = Ddiscvar.T.dTver_dTup;
dTver_dTdown = Ddiscvar.T.dTver_dTdown;

T_top = discvar.T.T_top;
dTtop_dTtop = Ddiscvar.T.dTtop_dTtop;
dTtop_dQ = Ddiscvar.T.dTtop_dQ;

ddTbdz_dTup = Ddiscvar.Tb.ddTdz_dTup;
ddTbdz_dTdown = Ddiscvar.Tb.ddTdz_dTdown;

Q_hor = discvar.Q.Q_hor;
dQhor_dQup = Ddiscvar.Q.dQhor_dQup;
dQhor_dQdown = Ddiscvar.Q.dQhor_dQdown;

T_bed_psigrid = discvar.T_bed.Tpsigrid;%where up and down are identified by v_up_node and v_down_node applied to (T_bdy_nodes_bed) 
dTbedpsi_dTbed_up = Ddiscvar.T_bed.dTpsigrid_dTbed_up;
dTbedpsi_dTbed_down = Ddiscvar.T_bed.dTpsigrid_dTbed_down;

%construct regularized bedwater content
[f_slide_Tbed, df_slide_Tbed] = regularization_old(T_bed, parameters);
[f_slide_Tbedpsigrid, df_slide_Tbedpsigrid] = regularization_old(T_bed_psigrid, parameters);
%% STREAM FUNCTION
dpsihorflux_dpsi = sparse(psi_up_node_hor,psi_up_node_hor,ddpsidx_dpsiu,psi_nodes, psi_nodes)+...
    sparse(psi_up_node_hor,psi_down_node_hor,ddpsidx_dpsid,psi_nodes, psi_nodes)-...
    (sparse(psi_down_node_hor,psi_up_node_hor,ddpsidx_dpsiu,psi_nodes, psi_nodes)+...
    sparse(psi_down_node_hor, psi_down_node_hor,ddpsidx_dpsid ,psi_nodes, psi_nodes));

dpsiverflux_dpsi = sparse(psi_up_node_ver,psi_up_node_ver,ddpsidz_dpsiu,psi_nodes, psi_nodes)+...
    sparse(psi_up_node_ver,psi_down_node_ver,ddpsidz_dpsid,psi_nodes, psi_nodes)-...
    (sparse(psi_down_node_ver,psi_up_node_ver,ddpsidz_dpsiu,psi_nodes, psi_nodes)+...
    sparse(psi_down_node_ver, psi_down_node_ver,ddpsidz_dpsid ,psi_nodes, psi_nodes));

dpsidztop_dpsitop =(-9)./(3*psi_Delta_z_cell(psi_bdy_nodes_top));
dpsidztop_dpsibtop = 1./(3*psi_Delta_z_cell(psi_bdy_nodes_top));

dpsiverflux_dpsi =  dpsiverflux_dpsi + (sparse(psi_bdy_nodes_top,psi_bdy_nodes_top,dpsidztop_dpsitop,psi_nodes,psi_nodes)+...
    sparse(psi_bdy_nodes_top ,psi_bdy_nodes_top + length(psi_bdy_nodes_top),dpsidztop_dpsibtop,psi_nodes,psi_nodes));

dpsiverflux_dpsi =  dpsiverflux_dpsi - (sparse(psi_bdy_nodes_bed,psi_bdy_nodes_bed,9./(3*psi_Delta_z_cell(psi_bdy_nodes_bed)),psi_nodes,psi_nodes))...
        - (sparse(psi_bdy_nodes_bed,psi_bdy_nodes_bed -length(psi_bdy_nodes_bed),-1./(3*psi_Delta_z_cell(psi_bdy_nodes_bed)),psi_nodes,psi_nodes));

%periodic boundary conditions
%flux_psi_in = (psi(psi_bdy_nodes_inflow) -psi(psi_bdy_nodes_outflow))./(psi_Delta_x_cell(psi_bdy_nodes_inflow)/2 + psi_Delta_x_cell(psi_bdy_nodes_outflow)/2);
dfluxpsiin_dpsiin = 1./(psi_Delta_x_cell(psi_bdy_nodes_inflow)/2 + psi_Delta_x_cell(psi_bdy_nodes_outflow)/2);
dfluxpsiin_dpsiout = -1./(psi_Delta_x_cell(psi_bdy_nodes_inflow)/2 + psi_Delta_x_cell(psi_bdy_nodes_outflow)/2);
%inflow boundary
%net_psi_horflux(psi_bdy_nodes_inflow) = net_psi_horflux(psi_bdy_nodes_inflow) - flux_psi_in;
dpsihorflux_dpsi =  dpsihorflux_dpsi - (sparse(psi_bdy_nodes_inflow,psi_bdy_nodes_inflow,dfluxpsiin_dpsiin,psi_nodes,psi_nodes) +...
    sparse(psi_bdy_nodes_inflow,psi_bdy_nodes_outflow,dfluxpsiin_dpsiout,psi_nodes,psi_nodes));
%outflow boundary
%net_psi_horflux(psi_bdy_nodes_outflow) = net_psi_horflux(psi_bdy_nodes_outflow) + flux_psi_in;
dpsihorflux_dpsi =  dpsihorflux_dpsi + (sparse(psi_bdy_nodes_outflow,psi_bdy_nodes_outflow,dfluxpsiin_dpsiout,psi_nodes,psi_nodes)+...
   sparse(psi_bdy_nodes_outflow,psi_bdy_nodes_inflow,dfluxpsiin_dpsiin,psi_nodes,psi_nodes));

%conservation law
Dfoutpsi_dpsi = spdiags(1./psi_Delta_x_cell,0,psi_nodes,psi_nodes)*dpsihorflux_dpsi + spdiags(1./psi_Delta_z_cell,0,psi_nodes,psi_nodes)*dpsiverflux_dpsi;
Dfoutpsi_domega = -speye(psi_nodes,psi_nodes);
fout(1:psi_nodes,:) = [Dfoutpsi_dpsi Dfoutpsi_domega sparse(psi_nodes,2*T_nodes+2*Q_nodes)];
faux(1:psi_nodes,:) =fout(1:psi_nodes,:);

%% SLIDING LAW
u_bed = (9*psi(psi_bdy_nodes_bed) - psi(psi_bdy_nodes_bed - length(psi_bdy_nodes_bed)))./(3*psi_Delta_z_cell(psi_bdy_nodes_bed));
tau_bed = gamma*(u_bed)./f_slide_Tbedpsigrid; 
dubed_dpsibed = (9)./(3*psi_Delta_z_cell(psi_bdy_nodes_bed));
dubed_dpsiabed = (-1)./(3*psi_Delta_z_cell(psi_bdy_nodes_bed));

dtaubed_dpsibed = gamma*dubed_dpsibed./f_slide_Tbedpsigrid;
dtaubed_dpsiabed = gamma*dubed_dpsiabed./f_slide_Tbedpsigrid;
dtaubed_dTbed = -(tau_bed./f_slide_Tbedpsigrid).*df_slide_Tbedpsigrid;

%at T cell centres
%derivative of sliding velocity
u_bed_full = [u_bed; u_bed(1)];
dubedfull_dpsibed = [dubed_dpsibed;dubed_dpsibed(1)];
dubedfull_dpsiabed = [dubed_dpsiabed;dubed_dpsiabed(1)];
dubed_dx = (u_bed_full(2:end) - u_bed_full(1:end-1))./Q_Delta_x_cell;

dubeddx_dpsiup_bed = -dubedfull_dpsibed(1:end-1)./Q_Delta_x_cell;
dubeddx_dpsiup_abed = -dubedfull_dpsiabed(1:end-1)./Q_Delta_x_cell;

dubeddx_dpsidown_bed = dubedfull_dpsibed(2:end)./Q_Delta_x_cell;
dubeddx_dpsidown_abed = dubedfull_dpsiabed(2:end)./Q_Delta_x_cell;

%sliding velocity and stress
u_bed_Tcentre = (u_bed_full(1:end-1)+u_bed_full(2:end))./2;
dubedTcentre_dpsiup_bed = dubedfull_dpsibed(1:end-1)./2;
dubedTcentre_dpsiup_abed = dubedfull_dpsiabed(1:end-1)./2;

dubedTcentre_dpsidown_bed = dubedfull_dpsibed(2:end)./2;
dubedTcentre_dpsidown_abed = dubedfull_dpsiabed(2:end)./2;

tau_bed_Tcentre = gamma*u_bed_Tcentre./f_slide_Tbed;  

dtaubedTcentre_dpsiup_bed = gamma*dubedTcentre_dpsiup_bed./f_slide_Tbed ;
dtaubedTcentre_dpsiup_abed = gamma*dubedTcentre_dpsiup_abed./f_slide_Tbed ;
dtaubedTcentre_dpsidown_bed = gamma*dubedTcentre_dpsidown_bed./f_slide_Tbed;
dtaubedTcentre_dpsidown_abed = gamma*dubedTcentre_dpsidown_abed./f_slide_Tbed;

dtaubedTcentre_dTbed = -(tau_bed_Tcentre./f_slide_Tbed).*df_slide_Tbed;
%% VORTICITY
domegahorflux_domega = sparse(psi_up_node_hor,psi_up_node_hor,ddomegadx_domegau,psi_nodes, psi_nodes)+...
    sparse(psi_up_node_hor,psi_down_node_hor,ddomegadx_domegad,psi_nodes, psi_nodes)-...
    (sparse(psi_down_node_hor,psi_up_node_hor,ddomegadx_domegau,psi_nodes, psi_nodes)+...
    sparse(psi_down_node_hor, psi_down_node_hor,ddomegadx_domegad ,psi_nodes, psi_nodes));

domegaverflux_domega = sparse(psi_up_node_ver,psi_up_node_ver,ddomegadz_domegau,psi_nodes, psi_nodes)+...
    sparse(psi_up_node_ver,psi_down_node_ver,ddomegadz_domegad,psi_nodes, psi_nodes)-...
    (sparse(psi_down_node_ver,psi_up_node_ver,ddomegadz_domegau,psi_nodes, psi_nodes)+...
    sparse(psi_down_node_ver, psi_down_node_ver,ddomegadz_domegad ,psi_nodes, psi_nodes));

%enforce boundary conditions
%ice surface: omega = 0
domegaverflux_domega =  domegaverflux_domega + (sparse(psi_bdy_nodes_top,psi_bdy_nodes_top,-1./(psi_Delta_z_cell(psi_bdy_nodes_top)/2),psi_nodes,psi_nodes));
    domegadzbed_dpsibed = (-8*dtaubed_dpsibed)./(3*psi_Delta_z_cell(psi_bdy_nodes_bed));
    domegadzbed_dpsiabed = (-8*dtaubed_dpsiabed)./(3*psi_Delta_z_cell(psi_bdy_nodes_bed));
    domegadzbed_dobed = 9./(3*psi_Delta_z_cell(psi_bdy_nodes_bed));
    domegadzbed_doabed = -1./(3*psi_Delta_z_cell(psi_bdy_nodes_bed));
    
    domegaverflux_dpsi = - (sparse(psi_bdy_nodes_bed,psi_bdy_nodes_bed,domegadzbed_dpsibed,psi_nodes,psi_nodes))...
        - (sparse(psi_bdy_nodes_bed,psi_bdy_nodes_bed- length(psi_bdy_nodes_bed),domegadzbed_dpsiabed,psi_nodes,psi_nodes));
    domegaverflux_domega = domegaverflux_domega - (sparse(psi_bdy_nodes_bed,psi_bdy_nodes_bed,domegadzbed_dobed,psi_nodes,psi_nodes)...
        + sparse(psi_bdy_nodes_bed,psi_bdy_nodes_bed - length(psi_bdy_nodes_bed),domegadzbed_doabed,psi_nodes,psi_nodes));
    domegaverflux_dT_bed = -(sparse(psi_bdy_nodes_bed,Q_up_node,-8*dtaubed_dTbed.*dTbedpsi_dTbed_up./(3*psi_Delta_z_cell(psi_bdy_nodes_bed)), psi_nodes,Q_nodes)+...
    sparse(psi_bdy_nodes_bed,Q_down_node,-8*dtaubed_dTbed.*dTbedpsi_dTbed_down./(3*psi_Delta_z_cell(psi_bdy_nodes_bed)), psi_nodes,Q_nodes));

%periodic boundary conditions
%flux_omega_in = (omega(psi_bdy_nodes_inflow) -omega(psi_bdy_nodes_outflow))./(psi_Delta_x_cell(psi_bdy_nodes_inflow)/2 + psi_Delta_x_cell(psi_bdy_nodes_outflow)/2);
dfluxomegain_domegain = 1./(psi_Delta_x_cell(psi_bdy_nodes_inflow)/2 + psi_Delta_x_cell(psi_bdy_nodes_outflow)/2);
dfluxomegain_domegaout = -1./(psi_Delta_x_cell(psi_bdy_nodes_inflow)/2 + psi_Delta_x_cell(psi_bdy_nodes_outflow)/2);

%inflow boundary
%net_omega_horflux(psi_bdy_nodes_inflow) = net_psi_horflux(psi_bdy_nodes_inflow) - flux_omega_in;
domegahorflux_domega =  domegahorflux_domega - (sparse(psi_bdy_nodes_inflow,psi_bdy_nodes_inflow,dfluxomegain_domegain,psi_nodes,psi_nodes) +...
    sparse(psi_bdy_nodes_inflow,psi_bdy_nodes_outflow,dfluxomegain_domegaout,psi_nodes,psi_nodes));
%outflow boundary
%net_omega_horflux(psi_bdy_nodes_outflow) = net_psi_horflux(psi_bdy_nodes_outflow) + flux_omega_in;
domegahorflux_domega =  domegahorflux_domega + (sparse(psi_bdy_nodes_outflow,psi_bdy_nodes_inflow,dfluxomegain_domegain,psi_nodes,psi_nodes) +...
    sparse(psi_bdy_nodes_outflow,psi_bdy_nodes_outflow,dfluxomegain_domegaout,psi_nodes,psi_nodes));

%conservation law
Dfoutomega_domega = spdiags(1./psi_Delta_x_cell,0,psi_nodes,psi_nodes)*domegahorflux_domega + spdiags(1./psi_Delta_z_cell,0,psi_nodes,psi_nodes)*domegaverflux_domega;
Dfoutomega_dpsi = spdiags(1./psi_Delta_z_cell,0,psi_nodes,psi_nodes)*domegaverflux_dpsi;
Dfoutomega_dT_bed = spdiags(1./psi_Delta_z_cell,0,psi_nodes,psi_nodes)*domegaverflux_dT_bed;
fout(psi_nodes+1:2*psi_nodes,:) = [ Dfoutomega_dpsi Dfoutomega_domega sparse(psi_nodes,2*T_nodes+Q_nodes) Dfoutomega_dT_bed];
faux(psi_nodes+1:2*psi_nodes,:) = fout(psi_nodes+1:2*psi_nodes,:);

%% VELOCITY FIELD IN THE BOUNDARY LAYER
%horizontal velocity at T hor edges
u_bed_Tedges = u_bed;
u_bed_edges_Tgrid = u_bed_Tedges(index_bed_to_T_hor_edges);

%vertical velocity at T cell centres w = - z du_bed/dx
W = -dubed_dx(index_bed_to_T_ver_edges).*T_coor_ver_edges_z;
dW_dpsi_up_bed = -dubeddx_dpsiup_bed(index_bed_to_T_ver_edges).*T_coor_ver_edges_z;
dW_dpsi_up_abed = -dubeddx_dpsiup_abed(index_bed_to_T_ver_edges).*T_coor_ver_edges_z;
dW_dpsi_down_bed = -dubeddx_dpsidown_bed(index_bed_to_T_ver_edges).*T_coor_ver_edges_z;
dW_dpsi_down_abed = -dubeddx_dpsidown_abed(index_bed_to_T_ver_edges).*T_coor_ver_edges_z;
%construct psi index
index_psibed_full = [psi_bdy_nodes_bed;psi_bdy_nodes_bed(1)];
index_psiabed_full = [psi_bdy_nodes_bed - length(psi_bdy_nodes_bed);psi_bdy_nodes_bed(1)- length(psi_bdy_nodes_bed)];

index_psibed_up = index_psibed_full(1:end-1);
index_psibed_down = index_psibed_full(2:end);
index_psiabed_up = index_psiabed_full(1:end-1);
index_psiabed_down = index_psiabed_full(2:end);

W_top = -dubed_dx.*T_coor_ver_edges_z_top;
dWtop_dpsi_up_bed = -dubeddx_dpsiup_bed.*T_coor_ver_edges_z_top;
dWtop_dpsi_up_abed = -dubeddx_dpsiup_abed.*T_coor_ver_edges_z_top;
dWtop_dpsi_down_bed = -dubeddx_dpsidown_bed.*T_coor_ver_edges_z_top;
dWtop_dpsi_down_abed = -dubeddx_dpsidown_abed.*T_coor_ver_edges_z_top;
%% Q EQUATION 

% Qbed_flux = Q_hor.*(u_bed);
dQbedflux_dpsi_bed = Q_hor.*dubed_dpsibed;
dQbedflux_dpsi_abed = Q_hor.*dubed_dpsiabed;
%jacobian
dQQ_dQ = sparse(Q_up_node,Q_up_node,(u_bed).*dQhor_dQup,Q_nodes, Q_nodes)+...
        sparse(Q_up_node,Q_down_node,(u_bed).*dQhor_dQdown ,Q_nodes, Q_nodes)-...
        (sparse(Q_down_node,Q_up_node,(u_bed).*dQhor_dQup,Q_nodes, Q_nodes)+...
        sparse(Q_down_node, Q_down_node,(u_bed).*dQhor_dQdown ,Q_nodes, Q_nodes));

dQQ_dpsi = sparse(Q_up_node,psi_bdy_nodes_bed,dQbedflux_dpsi_bed,Q_nodes, psi_nodes) +...
    sparse(Q_up_node,psi_bdy_nodes_bed -length(psi_bdy_nodes_bed),dQbedflux_dpsi_abed,Q_nodes, psi_nodes)-...
    sparse(Q_down_node,psi_bdy_nodes_bed,dQbedflux_dpsi_bed,Q_nodes, psi_nodes)-...
    sparse(Q_down_node,psi_bdy_nodes_bed-length(psi_bdy_nodes_bed),dQbedflux_dpsi_abed,Q_nodes, psi_nodes);

%source term
% source_Q = 2*Q.*dubed_dx;
dsourceterm_dQ = spdiags(2*dubed_dx, 0, Q_nodes,Q_nodes);
dsourceterm_dpsi = (sparse(1:Q_nodes, index_psibed_up, 2*Q.*dubeddx_dpsiup_bed,Q_nodes,psi_nodes) + ...
    sparse(1:Q_nodes, index_psiabed_up, 2*Q.*dubeddx_dpsiup_abed,Q_nodes,psi_nodes)+...
    sparse(1:Q_nodes, index_psibed_down, 2*Q.*dubeddx_dpsidown_bed,Q_nodes,psi_nodes)+...
    sparse(1:Q_nodes, index_psiabed_down, 2*Q.*dubeddx_dpsidown_abed,Q_nodes,psi_nodes));

%conservation law 
%fout(...) = net_Qbedflux./v_Delta_x_cell-source_Q;
DfoutQ_dQ =  spdiags(ones(Q_nodes,1)./dt, 0, Q_nodes,Q_nodes) + spdiags(1./Q_Delta_x_cell,0,Q_nodes,Q_nodes)*dQQ_dQ - dsourceterm_dQ;%
DfoutQ_dpsi =  spdiags(1./Q_Delta_x_cell,0,Q_nodes,Q_nodes)*dQQ_dpsi - dsourceterm_dpsi;
fout(2*psi_nodes+2*T_nodes+1:2*psi_nodes+2*T_nodes+Q_nodes, :) = [DfoutQ_dpsi sparse(Q_nodes, psi_nodes+2*T_nodes) DfoutQ_dQ sparse(Q_nodes,Q_nodes)];

DfoutQ_dQ_steady =   spdiags(1./Q_Delta_x_cell,0,Q_nodes,Q_nodes)*dQQ_dQ - dsourceterm_dQ;
faux(2*psi_nodes+2*T_nodes+1:2*psi_nodes+2*T_nodes+Q_nodes, :) = [DfoutQ_dpsi sparse(Q_nodes, psi_nodes+2*T_nodes) DfoutQ_dQ_steady sparse(Q_nodes,Q_nodes)];
%% HEAT EQUATION 
%solve dT/dt + d/dX((u_bed-V)T) + d/dZ[-Z du_b/dX T - dT/dZ] = 0
%nb: T grid is shifted with respect to psi grid, so T edges correspond to
%psi centres

%HORIZONTAL HEAT FLUX
dQHa_dT_u = (u_bed_edges_Tgrid).*dThor_dTup;
dQHa_dT_d = (u_bed_edges_Tgrid).*dThor_dTdown;
dQHa_dpsi_edge_bed = dubed_dpsibed(index_bed_to_T_hor_edges).*T_hor;
dQHa_dpsi_edge_abed = dubed_dpsiabed(index_bed_to_T_hor_edges).*T_hor;
index_abed = psi_bdy_nodes_bed - length(psi_bdy_nodes_bed);

dQHa_dT = sparse(T_up_node_hor,T_up_node_hor,dQHa_dT_u,T_nodes, T_nodes)+...
        sparse(T_up_node_hor,T_down_node_hor,dQHa_dT_d ,T_nodes, T_nodes)-...
        (sparse(T_down_node_hor,T_up_node_hor,dQHa_dT_u,T_nodes, T_nodes)+...
        sparse(T_down_node_hor, T_down_node_hor,dQHa_dT_d ,T_nodes, T_nodes));
dQHa_dpsi = sparse(T_up_node_hor,psi_bdy_nodes_bed(index_bed_to_T_hor_edges),dQHa_dpsi_edge_bed,T_nodes,psi_nodes)+...
    sparse(T_up_node_hor,index_abed(index_bed_to_T_hor_edges),dQHa_dpsi_edge_abed,T_nodes,psi_nodes)-...
    sparse(T_down_node_hor,psi_bdy_nodes_bed(index_bed_to_T_hor_edges),dQHa_dpsi_edge_bed,T_nodes,psi_nodes)-...
    sparse(T_down_node_hor,index_abed(index_bed_to_T_hor_edges),dQHa_dpsi_edge_abed,T_nodes,psi_nodes);
%VERTICAL HEAT FLUX 
%advective flux
dQVa_dT_u = W.*dTver_dTup;
dQVa_dT_d = W.*dTver_dTdown;
dQVa_dpsi_up_bed = dW_dpsi_up_bed.*T_vert;
dQVa_dpsi_up_abed = dW_dpsi_up_abed.*T_vert;
dQVa_dpsi_down_bed = dW_dpsi_down_bed.*T_vert;
dQVa_dpsi_down_abed = dW_dpsi_down_abed.*T_vert;
dQVa_dT = sparse(T_up_node_ver,T_up_node_ver,dQVa_dT_u,T_nodes, T_nodes)+...
        sparse(T_up_node_ver,T_down_node_ver,dQVa_dT_d ,T_nodes, T_nodes)-...
        (sparse(T_down_node_ver,T_up_node_ver,dQVa_dT_u,T_nodes, T_nodes)+...
        sparse(T_down_node_ver, T_down_node_ver,dQVa_dT_d ,T_nodes, T_nodes));
dQVa_dpsi = sparse(T_up_node_ver,index_psibed_up(index_bed_to_T_ver_edges),dQVa_dpsi_up_bed,T_nodes,psi_nodes)+...
    sparse(T_up_node_ver,index_psiabed_up(index_bed_to_T_ver_edges),dQVa_dpsi_up_abed,T_nodes,psi_nodes)+...
    sparse(T_up_node_ver,index_psibed_down(index_bed_to_T_ver_edges),dQVa_dpsi_down_bed,T_nodes,psi_nodes)+...
    sparse(T_up_node_ver,index_psiabed_down(index_bed_to_T_ver_edges),dQVa_dpsi_down_abed,T_nodes,psi_nodes)-...
    (sparse(T_down_node_ver,index_psibed_up(index_bed_to_T_ver_edges),dQVa_dpsi_up_bed,T_nodes,psi_nodes)+...
    sparse(T_down_node_ver,index_psiabed_up(index_bed_to_T_ver_edges),dQVa_dpsi_up_abed,T_nodes,psi_nodes)+...
    sparse(T_down_node_ver,index_psibed_down(index_bed_to_T_ver_edges),dQVa_dpsi_down_bed,T_nodes,psi_nodes)+...
     sparse(T_down_node_ver,index_psiabed_down(index_bed_to_T_ver_edges),dQVa_dpsi_down_abed,T_nodes,psi_nodes));
%correct for advective mass flux at the top boundary  
%net_a_verflux(T_bdy_nodes_top) =  net_a_verflux(T_bdy_nodes_top) + (W_top.*(T_top+ T(T_bdy_nodes_top))./2);
dQVa_dT =  dQVa_dT + (sparse(T_bdy_nodes_top,T_bdy_nodes_top,W_top./2,T_nodes,T_nodes))+...
    (sparse(T_bdy_nodes_top,T_bdy_nodes_top,W_top.*dTtop_dTtop./2,T_nodes,T_nodes));
dQVa_dpsi =  dQVa_dpsi + sparse(T_bdy_nodes_top,index_psibed_up,dWtop_dpsi_up_bed.*(T_top+ T(T_bdy_nodes_top))./2,T_nodes,psi_nodes)...
+    sparse(T_bdy_nodes_top,index_psiabed_up,dWtop_dpsi_up_abed.*(T_top+ T(T_bdy_nodes_top))./2,T_nodes,psi_nodes)... 
+ sparse(T_bdy_nodes_top,index_psibed_down,dWtop_dpsi_down_bed.*(T_top+ T(T_bdy_nodes_top))./2,T_nodes,psi_nodes)...
+ sparse(T_bdy_nodes_top,index_psiabed_down,dWtop_dpsi_down_abed.*(T_top+ T(T_bdy_nodes_top))./2,T_nodes,psi_nodes);
dQVa_dQ = sparse(T_bdy_nodes_top,1:Q_nodes,W_top.*dTtop_dQ./2,T_nodes,Q_nodes);
%diffusive flux
dQVd_dT_u = 1./T_Delta_z_edge;
dQVd_dT_d = -1./T_Delta_z_edge;
dQVd_dT = sparse(T_up_node_ver,T_up_node_ver,dQVd_dT_u,T_nodes, T_nodes)+...
        sparse(T_up_node_ver,T_down_node_ver,dQVd_dT_d ,T_nodes, T_nodes)-...
        (sparse(T_down_node_ver,T_up_node_ver,dQVd_dT_u,T_nodes, T_nodes)+...
        sparse(T_down_node_ver, T_down_node_ver,dQVd_dT_d ,T_nodes, T_nodes));
%correct for top flux
%net_d_verflux(T_bdy_nodes_top) =  net_d_verflux(T_bdy_nodes_top) + Q;    
dQVd_dQ = sparse(T_bdy_nodes_top,1:Q_nodes,ones( Q_nodes,1),T_nodes, Q_nodes);


% %Dirichlet condition at the bed 
% %basal energy budget
% bedflux_bl_ice = - (T(T_bdy_nodes_bed)-T_bed)./(T_Delta_z_cell(T_bdy_nodes_bed)/2);   
dbedfluxblice_dT = -1./(T_Delta_z_cell(T_bdy_nodes_bed)/2);
dbedfluxice_dTbed = 1./(T_Delta_z_cell(T_bdy_nodes_bed)/2);
% net_d_verflux(T_bdy_nodes_bed) = net_d_verflux(T_bdy_nodes_bed) - bedflux_bl_ice;
dQVd_dT =  dQVd_dT - sparse(T_bdy_nodes_bed,T_bdy_nodes_bed,dbedfluxblice_dT,T_nodes,T_nodes) ;
dQVd_dT_bed = -sparse(T_bdy_nodes_bed,1:Q_nodes,dbedfluxice_dTbed,T_nodes,Q_nodes) ;

%sum advective and diffusive fluxes
dQV_dT = dQVd_dT + dQVa_dT;
dQV_dQ = dQVd_dQ + dQVa_dQ;
dQV_dTbed = dQVd_dT_bed;
dQV_dpsi = dQVa_dpsi;
dQH_dT = dQHa_dT;%sparse(size(dQHa_dT,1), size(dQHa_dT,2));%
dQH_dpsi = dQHa_dpsi;%sparse(size(dQHa_dpsi,1), size(dQHa_dpsi,2));%
%conservation law
DfoutT_dT = spdiags(ones(T_nodes,1)./dt, 0, T_nodes,T_nodes)  + spdiags(1./T_Delta_z_cell,0,T_nodes,T_nodes)*dQV_dT+ spdiags(1./T_Delta_x_cell,0,T_nodes,T_nodes)*dQH_dT;
DfoutT_dT_steady =  spdiags(1./T_Delta_z_cell,0,T_nodes,T_nodes)*dQV_dT+ spdiags(1./T_Delta_x_cell,0,T_nodes,T_nodes)*dQH_dT;
DfoutT_dpsi =  spdiags(1./T_Delta_z_cell,0,T_nodes,T_nodes)*dQV_dpsi + spdiags(1./T_Delta_x_cell,0,T_nodes,T_nodes)*dQH_dpsi;
DfoutT_dQ = spdiags(1./T_Delta_z_cell,0,T_nodes,T_nodes)*dQV_dQ;
DfoutT_dT_bed = spdiags(1./T_Delta_z_cell,0,T_nodes,T_nodes)*dQV_dTbed;
fout(2*psi_nodes+1:2*psi_nodes+T_nodes,:) = [DfoutT_dpsi sparse(T_nodes,psi_nodes) DfoutT_dT sparse(T_nodes,T_nodes) DfoutT_dQ DfoutT_dT_bed];
faux(2*psi_nodes+1:2*psi_nodes+T_nodes,:) = [DfoutT_dpsi sparse(T_nodes,psi_nodes) DfoutT_dT_steady sparse(T_nodes,T_nodes) DfoutT_dQ DfoutT_dT_bed];
%% BED HEAT EQUATION

% %VERTICAL HEAT FLUX
% QVdbed = -dTb_dz;
dQVbed_dT_u = -ddTbdz_dTup;
dQVbed_dT_d = -ddTbdz_dTdown;
dQVbed_dTb = sparse(Tb_up_node_ver,Tb_up_node_ver,dQVbed_dT_u,T_nodes, T_nodes)+...
        sparse(Tb_up_node_ver,Tb_down_node_ver,dQVbed_dT_d ,T_nodes, T_nodes)-...
        (sparse(Tb_down_node_ver,Tb_up_node_ver,dQVbed_dT_u,T_nodes, T_nodes)+...
        sparse(Tb_down_node_ver, Tb_down_node_ver,dQVbed_dT_d,T_nodes, T_nodes));
% %correct for bed flux
% bedflux_bl_bed = - (T_bed-Tb(Tb_bdy_nodes_top))./(Tb_Delta_z_cell(Tb_bdy_nodes_top)/2);   
dbedfluxbed_dTb = 1./(Tb_Delta_z_cell(Tb_bdy_nodes_top)/2);
dbedfluxbed_dTbed = -1./(Tb_Delta_z_cell(Tb_bdy_nodes_top)/2);
% net_bedverflux(Tb_bdy_nodes_top) = net_bedverflux(Tb_bdy_nodes_top) + bedflux_bl_bed;
dQVbed_dTb = dQVbed_dTb + sparse(Tb_bdy_nodes_top,Tb_bdy_nodes_top,dbedfluxbed_dTb, T_nodes,T_nodes);
dQVbed_dTbed = sparse(Tb_bdy_nodes_top,1:length(Tb_bdy_nodes_top),dbedfluxbed_dTbed, T_nodes,Q_nodes);

% %enforce conservation law
% fout(2*psi_nodes+T_nodes+1:2*psi_nodes+2*T_nodes) = net_bedhorflux./Tb_Delta_x_cell + net_bedverflux./Tb_Delta_z_cell;
DfoutTb_dT = spdiags(ones(T_nodes,1)./dt, 0, T_nodes,T_nodes)  + spdiags(1./Tb_Delta_z_cell,0,T_nodes,T_nodes)*dQVbed_dTb;
DfoutTb_dT_steady =  spdiags(1./Tb_Delta_z_cell,0,T_nodes,T_nodes)*dQVbed_dTb;
DfoutTb_dT_bed = spdiags(1./Tb_Delta_z_cell,0,T_nodes,T_nodes)*dQVbed_dTbed;
fout(2*psi_nodes+T_nodes+1:2*psi_nodes+2*T_nodes,:) = [sparse(T_nodes,2*psi_nodes+T_nodes) DfoutTb_dT sparse(T_nodes,Q_nodes) DfoutTb_dT_bed];
faux(2*psi_nodes+T_nodes+1:2*psi_nodes+2*T_nodes,:) = [sparse(T_nodes,2*psi_nodes+T_nodes) DfoutTb_dT_steady sparse(T_nodes,Q_nodes) DfoutTb_dT_bed];
%% BASAL ENERGY BUDGET
% solve m=0 on z =0, 1D version of the T grid
% m = -bedflux_bl_ice + bedflux_bl_bed +alpha*tau_bed_Tcentre.*u_bed_Tcentre;

dsourceterm_dT = sparse(Q_nodes,T_nodes) - (sparse(1:Q_nodes,T_bdy_nodes_bed, dbedfluxblice_dT,Q_nodes,T_nodes) );
dsourceterm_dTb = sparse(Q_nodes,T_nodes) + (sparse(1:Q_nodes,Tb_bdy_nodes_top, dbedfluxbed_dTb,Q_nodes,T_nodes) );
dsourceterm_dTbed = -sparse(1:Q_nodes,1:Q_nodes,dbedfluxice_dTbed, Q_nodes,Q_nodes)+sparse(1:Q_nodes,1:Q_nodes,dbedfluxbed_dTbed, Q_nodes,Q_nodes);

dfriction_dpsiup_bed = alpha*(tau_bed_Tcentre.*dubedTcentre_dpsiup_bed+ u_bed_Tcentre.*dtaubedTcentre_dpsiup_bed);
dfriction_dpsiup_abed = alpha*(tau_bed_Tcentre.*dubedTcentre_dpsiup_abed + u_bed_Tcentre.*dtaubedTcentre_dpsiup_abed);
dfriction_dpsidown_bed = alpha*(tau_bed_Tcentre.*dubedTcentre_dpsidown_bed + u_bed_Tcentre.*dtaubedTcentre_dpsidown_bed);
dfriction_dpsidown_abed = alpha*(tau_bed_Tcentre.*dubedTcentre_dpsidown_abed + u_bed_Tcentre.*dtaubedTcentre_dpsidown_abed);

dfriction_dTbed = alpha*dtaubedTcentre_dTbed.*u_bed_Tcentre;

dsourceterm_dTbed = dsourceterm_dTbed + sparse(1:Q_nodes,1:Q_nodes,dfriction_dTbed,Q_nodes,Q_nodes);
dsourceterm_dpsi = sparse(1:Q_nodes,index_psibed_up,dfriction_dpsiup_bed,Q_nodes,psi_nodes)+...
    sparse(1:Q_nodes,index_psiabed_up,dfriction_dpsiup_abed,Q_nodes,psi_nodes)+...
    sparse(1:Q_nodes,index_psibed_down,dfriction_dpsidown_bed,Q_nodes,psi_nodes)+...
    sparse(1:Q_nodes,index_psiabed_down,dfriction_dpsidown_abed,Q_nodes,psi_nodes);
% %conservation law
% fout(2*psi_nodes+2*T_nodes+Q_nodes+1:2*psi_nodes+2*T_nodes+2*Q_nodes) = m;
Dfoutv_dT = dsourceterm_dT;
Dfoutv_dTb = dsourceterm_dTb;
Dfoutv_dpsi = dsourceterm_dpsi;
Dfoutv_dTbed = dsourceterm_dTbed;
fout(2*psi_nodes+2*T_nodes+Q_nodes+1:2*psi_nodes+2*T_nodes+2*Q_nodes,:) =[ Dfoutv_dpsi sparse(Q_nodes,psi_nodes) Dfoutv_dT Dfoutv_dTb sparse(Q_nodes, Q_nodes) Dfoutv_dTbed];
faux(2*psi_nodes+2*T_nodes+Q_nodes+1:2*psi_nodes+2*T_nodes+2*Q_nodes,:) = fout(2*psi_nodes+2*T_nodes+Q_nodes+1:2*psi_nodes+2*T_nodes+2*Q_nodes,:);
