function [fout, Dfout] = discretisation_v3(v_in, parameters)
%account for variable cell edge and network edge lengths in both
%directions; implements centered scheme for both vorticity and stream
%function

%psi grid
psi_nodes = parameters.grid.psi.n_nodes.tot;                                   %number of nodes
psi_up_node_ver = parameters.grid.psi.up_node.vert;                            %list (n_edges-by-1 vector) of 'upstream' node for each vertical edge
psi_down_node_ver = parameters.grid.psi.down_node.vert;                        %list (n_edges-by-1 vector) of 'downstream' nodes for each vertical edge
psi_up_node_hor = parameters.grid.psi.up_node.hor;                             %list (n_edges-by-1 vector) of 'upstream' node for each horizontal edge
psi_down_node_hor = parameters.grid.psi.down_node.hor;                         %list (n_edges-by-1 vector) of 'downstream' nodes for each horizontal edge

psi_Delta_z_cell = parameters.grid.psi.Delta_z_cell;                         %list (n_edges-by-1 vector) of lengths of vertical cell edges;
psi_Delta_x_cell = parameters.grid.psi.Delta_x_cell;                         %length of cells, hor (list)
psi_Delta_z_edge = parameters.grid.psi.Delta_z_edge;                         %list (n_edges-by-1 vector) of lengths of vertical cell edges;
psi_Delta_x_edge = parameters.grid.psi.Delta_x_edge;                         %length of cells, hor (list)

%T ice grid
T_nodes = parameters.grid.T.ice.n_nodes.tot;                                   %number of nodes
T_n_edges_hor = parameters.grid.T.ice.n_edges.hor;
T_up_node_ver = parameters.grid.T.ice.up_node.vert;                            %list (n_edges-by-1 vector) of 'upstream' node for each vertical edge
T_down_node_ver = parameters.grid.T.ice.down_node.vert;                        %list (n_edges-by-1 vector) of 'downstream' nodes for each vertical edge
T_up_node_hor = parameters.grid.T.ice.up_node.hor;                             %list (n_edges-by-1 vector) of 'upstream' node for each horizontal edge
T_down_node_hor = parameters.grid.T.ice.down_node.hor;                         %list (n_edges-by-1 vector) of 'downstream' nodes for each horizontal edge

T_bdy_nodes_top = parameters.grid.T.ice.bdy_nodes.top;        

T_Delta_z_edge = parameters.grid.T.ice.Delta_z_edge;                             %length of cells, ver (list)
T_Delta_z_cell = parameters.grid.T.ice.Delta_z_cell; 
T_Delta_x_cell = parameters.grid.T.ice.Delta_x_cell; 
T_Delta_x_edge = parameters.grid.T.ice.Delta_x_edge;                             %length of cells, hor (list)

%T bed grid
Tb_n_edges_hor = parameters.grid.T.bed.n_edges.hor;
Tb_up_node_ver = parameters.grid.T.bed.up_node.vert;                            %list (n_edges-by-1 vector) of 'upstream' node for each vertical edge
Tb_down_node_ver = parameters.grid.T.bed.down_node.vert;                        %list (n_edges-by-1 vector) of 'downstream' nodes for each vertical edge
Tb_up_node_hor = parameters.grid.T.bed.up_node.hor;                             %list (n_edges-by-1 vector) of 'upstream' node for each horizontal edge
Tb_down_node_hor = parameters.grid.T.bed.down_node.hor;                         %list (n_edges-by-1 vector) of 'downstream' nodes for each horizontal edge

Tb_Delta_z_edge = parameters.grid.T.bed.Delta_z_edge;                             %length of cells, ver (list)
Tb_Delta_z_cell = parameters.grid.T.bed.Delta_z_cell; 
Tb_Delta_x_edge = parameters.grid.T.bed.Delta_x_edge;                             %length of cells, hor (list)

%FLUX GRID
Q_nodes = parameters.grid.Q.n_nodes.tot;  
Q_up_node = parameters.grid.Q.up_node.hor;                                 %list (n_edges-by-1 vector) of 'upstream' node for each horizontal edge
Q_down_node = parameters.grid.Q.down_node.hor;  
Q_n_edges_hor = parameters.grid.Q.n_edges.hor;

%unpack input variable v_in
psi = v_in(1:psi_nodes);
omega = v_in(psi_nodes+1:2*psi_nodes);
T = v_in(2*psi_nodes+1:2*psi_nodes+T_nodes);
Tb = v_in(2*psi_nodes+T_nodes+1:2*psi_nodes+2*T_nodes);
Q = v_in(2*psi_nodes+2*T_nodes+1:2*psi_nodes+2*T_nodes+Q_nodes);
T_bed = v_in(2*psi_nodes+2*T_nodes+Q_nodes+1:end);

%STREAM FUNCTION
%along edges
fout.psi_hor = (psi(psi_up_node_hor).*psi_Delta_x_cell(psi_down_node_hor)+psi(psi_down_node_hor).*psi_Delta_x_cell(psi_up_node_hor))./(psi_Delta_x_cell(psi_up_node_hor)+psi_Delta_x_cell(psi_down_node_hor));
fout.psi_vert = (psi(psi_up_node_ver).*psi_Delta_z_cell(psi_down_node_ver)+psi(psi_down_node_ver).*psi_Delta_z_cell(psi_up_node_ver))./(psi_Delta_z_cell(psi_up_node_ver)+psi_Delta_z_cell(psi_down_node_ver));
fout.dpsi_dx = (psi(psi_down_node_hor)-psi(psi_up_node_hor))./psi_Delta_x_edge;
fout.dpsi_dz = (psi(psi_down_node_ver)-psi(psi_up_node_ver))./psi_Delta_z_edge;

%derivatives
Dfout.dpsihor_dpsiup = psi_Delta_x_cell(psi_down_node_hor)./(psi_Delta_x_cell(psi_up_node_hor)+psi_Delta_x_cell(psi_down_node_hor));
Dfout.dpsihor_dpsidown = psi_Delta_x_cell(psi_up_node_hor)./(psi_Delta_x_cell(psi_up_node_hor)+psi_Delta_x_cell(psi_down_node_hor));
Dfout.dpsiver_dpsiup = psi_Delta_z_cell(psi_down_node_ver)./(psi_Delta_z_cell(psi_up_node_ver)+psi_Delta_z_cell(psi_down_node_ver));
Dfout.dpsiver_dpsidown = psi_Delta_z_cell(psi_up_node_ver)./(psi_Delta_z_cell(psi_up_node_ver)+psi_Delta_z_cell(psi_down_node_ver));

Dfout.ddpsidx_dpsiup = -1./psi_Delta_x_edge;
Dfout.ddpsidx_dpsidown = 1./psi_Delta_x_edge;
Dfout.ddpsidz_dpsiup = -1./psi_Delta_z_edge;
Dfout.ddpsidz_dpsidown = 1./psi_Delta_z_edge;

%VORTICITY
%along edges
fout.omega_hor = (omega(psi_up_node_hor).*psi_Delta_x_cell(psi_down_node_hor)+omega(psi_down_node_hor).*psi_Delta_x_cell(psi_up_node_hor))./(psi_Delta_x_cell(psi_up_node_hor)+psi_Delta_x_cell(psi_down_node_hor));
fout.omega_vert = (omega(psi_up_node_ver).*psi_Delta_z_cell(psi_down_node_ver)+omega(psi_down_node_ver).*psi_Delta_z_cell(psi_up_node_ver))./(psi_Delta_z_cell(psi_up_node_ver)+psi_Delta_z_cell(psi_down_node_ver));

fout.domega_dx = (omega(psi_down_node_hor)-omega(psi_up_node_hor))./psi_Delta_x_edge;
fout.domega_dz = (omega(psi_down_node_ver)-omega(psi_up_node_ver))./psi_Delta_z_edge;

%derivatives
Dfout.domegahor_dpsiup = psi_Delta_x_cell(psi_down_node_hor)./(psi_Delta_x_cell(psi_up_node_hor)+psi_Delta_x_cell(psi_down_node_hor));
Dfout.domegahor_dpsidown = psi_Delta_x_cell(psi_up_node_hor)./(psi_Delta_x_cell(psi_up_node_hor)+psi_Delta_x_cell(psi_down_node_hor));
Dfout.domegaver_dpsiup = psi_Delta_z_cell(psi_down_node_ver)./(psi_Delta_z_cell(psi_up_node_ver)+psi_Delta_z_cell(psi_down_node_ver));
Dfout.domegaver_dpsidown = psi_Delta_z_cell(psi_up_node_ver)./(psi_Delta_z_cell(psi_up_node_ver)+psi_Delta_z_cell(psi_down_node_ver));

Dfout.ddomegadx_domegaup = -1./psi_Delta_x_edge;
Dfout.ddomegadx_domegadown = 1./psi_Delta_x_edge;
Dfout.ddomegadz_domegaup = -1./psi_Delta_z_edge;
Dfout.ddomegadz_domegadown = 1./psi_Delta_z_edge;

%TEMPERATURE FIELD IN THE ICE
%along T edges 
theta_T = 0; %upwinding
fout.T.T_hor = T(T_up_node_hor).*(1-theta_T)+T(T_down_node_hor).*theta_T;
fout.T.T_vert = (T(T_up_node_ver).*T_Delta_z_cell(T_down_node_ver)+T(T_down_node_ver).*T_Delta_z_cell(T_up_node_ver))./(T_Delta_z_cell(T_up_node_ver)+T_Delta_z_cell(T_down_node_ver));
fout.T.dT_dx = (T(T_down_node_hor)-T(T_up_node_hor))./T_Delta_x_edge;
fout.T.dT_dz = (T(T_down_node_ver)-T(T_up_node_ver))./T_Delta_z_edge;


Dfout.T.dThor_dTup = (1-theta_T)*ones(T_n_edges_hor,1);
Dfout.T.dThor_dTdown = theta_T*ones(T_n_edges_hor,1);
Dfout.T.dTver_dTup = T_Delta_z_cell(T_down_node_ver)./(T_Delta_z_cell(T_up_node_ver)+T_Delta_z_cell(T_down_node_ver));
Dfout.T.dTver_dTdown = T_Delta_z_cell(T_up_node_ver)./(T_Delta_z_cell(T_up_node_ver)+T_Delta_z_cell(T_down_node_ver));
Dfout.T.ddTdx_dTup = -1./T_Delta_x_edge;
Dfout.T.ddTdx_dTdown = 1./T_Delta_x_edge;
Dfout.T.ddTdz_dTup = -1./T_Delta_z_edge;
Dfout.T.ddTdz_dTdown = 1./T_Delta_z_edge;

fout.T.T_top = T(T_bdy_nodes_top) - Q.*T_Delta_z_cell(T_bdy_nodes_top);
Dfout.T.dTtop_dTtop = ones(length(T_bdy_nodes_top),1);
Dfout.T.dTtop_dQ = -T_Delta_z_cell(T_bdy_nodes_top);

%TEMPERATURE FIELD IN THE BED
fout.Tb.dT_dx = (Tb(Tb_down_node_hor)-Tb(Tb_up_node_hor))./Tb_Delta_x_edge;
fout.Tb.dT_dz = (Tb(Tb_down_node_ver)-Tb(Tb_up_node_ver))./Tb_Delta_z_edge;

Dfout.Tb.ddTdx_dTup = -1./Tb_Delta_x_edge;
Dfout.Tb.ddTdx_dTdown = 1./Tb_Delta_x_edge;
Dfout.Tb.ddTdz_dTup = -1./Tb_Delta_z_edge;
Dfout.Tb.ddTdz_dTdown = 1./Tb_Delta_z_edge;

%FLUX IN THE BED
fout.Q.Q_hor = Q(Q_up_node);
Dfout.Q.dQhor_dQup = ones(Q_n_edges_hor,1);
Dfout.Q.dQhor_dQdown = zeros(Q_n_edges_hor,1);

%compute T bed at psi grid centres
fout.T_bed.Tpsigrid = (T_bed(Q_up_node) +T_bed(Q_down_node))./2;
Dfout.T_bed.dTpsigrid_dTbed_up = ones(length(Q_up_node),1)./2;
Dfout.T_bed.dTpsigrid_dTbed_down = ones(length(Q_up_node),1)./2;

end


