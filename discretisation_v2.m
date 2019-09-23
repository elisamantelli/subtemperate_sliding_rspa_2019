function [fout, Dfout] = discretisation_v2(v_in, parameters)

%unpack parameters
h_nodes = parameters.grid_h.n_nodes;                  %number of nodes,h
h_edges = parameters.grid_h.n_edges;                  %number of nodes,h
h_up_node = parameters.grid_h.up_node;                %list (h_n_edge-by-1 vector) of 'upstream' node, h
h_down_node = parameters.grid_h.down_node;            %list (h_n_edges-by-1 vector) of 'downstream' node,h

T_nodes = parameters.grid_T.n_nodes.tot;                                   %number of temperature nodes
nodes_ver = parameters.grid_T.n_nodes.vert;                                %number of nodes in the vertical direction
T_n_edges_hor = parameters.grid_T.n_edges.hor;                             %number of horizontal edges, T
T_n_edges_ver = parameters.grid_T.n_edges.vert;                            %number of vertical edges,T
T_up_node_ver = parameters.grid_T.up_node.vert;                            %list (n_edges-by-1 vector) of 'upstream' node for each vertical edge,T
T_down_node_ver = parameters.grid_T.down_node.vert;                        %list (n_edges-by-1 vector) of 'downstream' nodes for each vertical edge,T
T_up_node_hor = parameters.grid_T.up_node.hor;                             %list (n_edges-by-1 vector) of 'upstream' node for each horizontal edge,T
T_down_node_hor = parameters.grid_T.down_node.hor;                         %list (n_edges-by-1 vector) of 'downstream' nodes for each horizontal edge,T
T_Delta_eta = parameters.grid_T.Delta_eta;                                 %scalar, vertical spacing between cell centres
T_bed_nodes = parameters.grid_T.bdy_nodes.bed_nodes;                       %list of bed nodes


n_x.c = parameters.n_x.c;    %number of h nodes, cold subdomain
n_x.s = parameters.n_x.s;    %number of h nodes, subtemperate subdomain

%unpack input variable v_in
H = v_in(1:h_nodes);
T = v_in(h_nodes+1:h_nodes+T_nodes);
xc = v_in(h_nodes+T_nodes+1);
xs = v_in(h_nodes+T_nodes+2);
xt = v_in(h_nodes+T_nodes+3);

%horizontal length of subdomains, cells and edges 
[length_h, Dlength_h] = fvlength([xc;xs;xt], parameters,'h');
%unpack
Delta_x_h = length_h.Delta_x;
dDeltax_dxc = Dlength_h. dDeltax_dxc;
dDeltax_dxs = Dlength_h. dDeltax_dxs;
dDeltax_dxt = Dlength_h. dDeltax_dxt;

Delta_sigma_h = length_h.Delta_sigma;

L_h = length_h.L;
dL_dxc = Dlength_h.dL_dxc;
dL_dxs = Dlength_h.dL_dxs;
dL_dxt = Dlength_h.dL_dxt;

%BED
[st_bed,Dst_bed] = bed([xc;xs;xt], parameters);
%unpack: bed slope along hor edges (list h_n_edges-by-one), at grounding
%line, and depth to bed at grounding line
dBx_edges = st_bed.dB.edges;
dBx_g = st_bed.dB.g;
ddBxedges_dxc = Dst_bed.ddBx_dxc;
ddBxedges_dxs = Dst_bed.ddBx_dxs;
ddBxedges_dxt = Dst_bed.ddBx_dxt;
ddBg_dxs = Dst_bed.ddBg_dxs;
ddBg_dxt = Dst_bed.ddBg_dxt;

%UPWINDING
theta = upwinding(parameters);
theta_T = theta.theta_T;
theta_h = theta.theta_h;
dthetah_dHup = theta.dtheta_h_dhu;
dthetah_dHdown = theta.dtheta_h_dhd;

%ICE THICKNESS AND DRIVING STRESS 
%implement linear interpolation

%averaged scheme
den = (L_h(h_up_node).*Delta_sigma_h(h_up_node)+L_h(h_down_node).*Delta_sigma_h(h_down_node));
dden_dxc = (dL_dxc(h_up_node).*Delta_sigma_h(h_up_node)+dL_dxc(h_down_node).*Delta_sigma_h(h_down_node));
dden_dxs = (dL_dxs(h_up_node).*Delta_sigma_h(h_up_node)+dL_dxs(h_down_node).*Delta_sigma_h(h_down_node));
dden_dxt = (dL_dxt(h_up_node).*Delta_sigma_h(h_up_node)+dL_dxt(h_down_node).*Delta_sigma_h(h_down_node));

num = H(h_up_node).*L_h(h_down_node).*Delta_sigma_h(h_down_node)+H(h_down_node).*L_h(h_up_node).*Delta_sigma_h(h_up_node);
dnum_dxc = H(h_up_node).*dL_dxc(h_down_node).*Delta_sigma_h(h_down_node)+H(h_down_node).*dL_dxc(h_up_node).*Delta_sigma_h(h_up_node);
dnum_dxs = H(h_up_node).*dL_dxs(h_down_node).*Delta_sigma_h(h_down_node)+H(h_down_node).*dL_dxs(h_up_node).*Delta_sigma_h(h_up_node);
dnum_dxt = H(h_up_node).*dL_dxt(h_down_node).*Delta_sigma_h(h_down_node)+H(h_down_node).*dL_dxt(h_up_node).*Delta_sigma_h(h_up_node);

H_edge = num./den;
dHedge_dHup = L_h(h_down_node).*Delta_sigma_h(h_down_node)./den;
dHedge_dHdown = L_h(h_up_node).*Delta_sigma_h(h_up_node)./den;
dHedge_dxc = (dnum_dxc.*den-num.*dden_dxc)./den.^2;
dHedge_dxs = (dnum_dxs.*den-num.*dden_dxs)./den.^2;
dHedge_dxt = (dnum_dxt.*den-num.*dden_dxt)./den.^2;

%assemble output
fout.H.H_edge = H_edge;
Dfout.H.dHedge_dHup = dHedge_dHup;
Dfout.H.dHedge_dHdown = dHedge_dHdown;
Dfout.H.dHedge_dxc = dHedge_dxc;
Dfout.H.dHedge_dxs = dHedge_dxs;
Dfout.H.dHedge_dxt = dHedge_dxt;

%upwinding scheme
fout.H.H_edge_adv = H(h_up_node).*(1-theta_h)+H(h_down_node).*theta_T;
Dfout.H.dHedgea_dHup = (1-theta_h)+ H(h_up_node)*(-dthetah_dHup)+H(h_down_node).*dthetah_dHup;
Dfout.H.dHedgea_dHdown = H(h_up_node)*(-dthetah_dHdown)+theta_h+H(h_down_node).*dthetah_dHdown;

fout.dS.dS_edge = (H(h_down_node)-H(h_up_node))./Delta_x_h - dBx_edges;
fout.dH_dx =  (H(h_down_node)-H(h_up_node))./Delta_x_h;
Dfout.dS.ddSedge_dHup = -1./Delta_x_h;
Dfout.dS.ddSedge_dHdown = 1./Delta_x_h;
Dfout.dS.ddSedge_dxc = -(H(h_down_node)-H(h_up_node))./Delta_x_h.^2.*dDeltax_dxc - ddBxedges_dxc; 
Dfout.dS.ddSedge_dxs = -(H(h_down_node)-H(h_up_node))./Delta_x_h.^2.*dDeltax_dxs - ddBxedges_dxs;
Dfout.dS.ddSedge_dxt = -(H(h_down_node)-H(h_up_node))./Delta_x_h.^2.*dDeltax_dxt - ddBxedges_dxt;

%at the grounding line
fout.H.H_g = (3*H(end)-H(end-1))/2;
Dfout.H.dHg_dH = [sparse(h_nodes-2,1); -1/2 ; 3/2];

fout.dS.dS_g = (H(end)-H(end-1))/Delta_x_h(end) - dBx_g;
Dfout.dS.ddSg_dH = [sparse(h_nodes-2,1); -1/Delta_x_h(end) ; 1/Delta_x_h(end)];
Dfout.dS.ddSg_dxc = 0;
Dfout.dS.ddSg_dxs = -(H(end)-H(end-1))/Delta_x_h(end).^2*dDeltax_dxs(end) - ddBg_dxs;
Dfout.dS.ddSg_dxt = -(H(end)-H(end-1))/Delta_x_h(end).^2*dDeltax_dxt(end) - ddBg_dxt;

%TEMPERATURE 
%along T edges 
fout.T.T_hor = T(T_up_node_hor).*(1-theta_T)+T(T_down_node_hor).*theta_T;
% fout.T.T_edge = 1/2*(T(T_up_node_hor)+T(T_down_node_hor));
fout.T.T_vert = 1/2*(T(T_up_node_ver)+T(T_down_node_ver));

Dfout.T.dThor_dTup = (1-theta_T)*ones(T_n_edges_hor,1);
Dfout.T.dThor_dTdown = theta_T*ones(T_n_edges_hor,1);
% Dfout.T.dTedge_dTup = 1/2*ones(T_n_edges_hor,1);
% Dfout.T.dTedge_dTdown = 1/2*ones(T_n_edges_hor,1);
Dfout.T.dTver_dTup = 1/2*ones(T_n_edges_ver,1);
Dfout.T.dTver_dTdown = 1/2*ones(T_n_edges_ver,1);

%at the grounding line
fout.T.T_hor_g = T(end-nodes_ver+1:end);
T_up_node_hor_g = T_down_node_hor(end-nodes_ver+1:end);
dTupnode_dT_g = ones(nodes_ver,1);
Dfout.T.dThorg_dT_up = sparse(1:nodes_ver, T_up_node_hor_g, dTupnode_dT_g, nodes_ver, T_nodes);

%BED FLUX
T_bed = T(T_bed_nodes);
QH_v = -1/(T_Delta_eta/2)*(T_bed)./H;
dQHv_dTbed = -1./(T_Delta_eta/2*H);
dQHv_dH = 1/(T_Delta_eta/2)*(T_bed)./H.^2;

fout.flux.QH_v = QH_v;

if strcmp(parameters.flag.u_sl, 'flux_upwind')==1
    
    %assemble output
    fout.flux.QH_v_edge = fout.flux.QH_v(h_up_node);
    Dfout.flux.dQHv_dQHvup = ones(h_edges,1);
    Dfout.flux.dQHv_dQHvdown = zeros(h_edges,1);
    
    Dfout.flux.dQHv_dTbedup = Dfout.flux.dQHv_dQHvup.*dQHv_dTbed(h_up_node);
    Dfout.flux.dQHv_dTbeddown = zeros(h_edges,1);
    Dfout.flux.dQHv_dHup = Dfout.flux.dQHv_dQHvup.*dQHv_dH(h_up_node);
    Dfout.flux.dQHv_dHdown = zeros(h_edges,1);
    
    Dfout.flux.dQHv_dxc = zeros(h_edges,1);
    Dfout.flux.dQHv_dxs = zeros(h_edges,1);
    Dfout.flux.dQHv_dxt = zeros(h_edges,1);
    
elseif strcmp(parameters.flag.u_sl, 'flux_av')==1
    den = (L_h(h_up_node).*Delta_sigma_h(h_up_node)+L_h(h_down_node).*Delta_sigma_h(h_down_node));
    dden_dxc = (dL_dxc(h_up_node).*Delta_sigma_h(h_up_node)+dL_dxc(h_down_node).*Delta_sigma_h(h_down_node));
    dden_dxs = (dL_dxs(h_up_node).*Delta_sigma_h(h_up_node)+dL_dxs(h_down_node).*Delta_sigma_h(h_down_node));
    dden_dxt = (dL_dxt(h_up_node).*Delta_sigma_h(h_up_node)+dL_dxt(h_down_node).*Delta_sigma_h(h_down_node));
    
    num = QH_v(h_up_node).*L_h(h_down_node).*Delta_sigma_h(h_down_node)+QH_v(h_down_node).*L_h(h_up_node).*Delta_sigma_h(h_up_node);
    dnum_dxc = QH_v(h_up_node).*dL_dxc(h_down_node).*Delta_sigma_h(h_down_node)+QH_v(h_down_node).*dL_dxc(h_up_node).*Delta_sigma_h(h_up_node);
    dnum_dxs = QH_v(h_up_node).*dL_dxs(h_down_node).*Delta_sigma_h(h_down_node)+QH_v(h_down_node).*dL_dxs(h_up_node).*Delta_sigma_h(h_up_node);
    dnum_dxt = QH_v(h_up_node).*dL_dxt(h_down_node).*Delta_sigma_h(h_down_node)+QH_v(h_down_node).*dL_dxt(h_up_node).*Delta_sigma_h(h_up_node);
    
    %assemble output
    fout.flux.QH_v_edge = num./den;
    Dfout.flux.dQHv_dQHvup = L_h(h_down_node).*Delta_sigma_h(h_down_node)./den;
    Dfout.flux.dQHv_dQHvdown = L_h(h_up_node).*Delta_sigma_h(h_up_node)./den;
    
    Dfout.flux.dQHv_dTbedup = Dfout.flux.dQHv_dQHvup.*dQHv_dTbed(h_up_node);
    Dfout.flux.dQHv_dTbeddown = Dfout.flux.dQHv_dQHvdown.*dQHv_dTbed(h_down_node);
    Dfout.flux.dQHv_dHup = Dfout.flux.dQHv_dQHvup.*dQHv_dH(h_up_node);
    Dfout.flux.dQHv_dHdown = Dfout.flux.dQHv_dQHvdown.*dQHv_dH(h_down_node);
 
    Dfout.flux.dQHv_dxc = (dnum_dxc.*den-num.*dden_dxc)./den.^2;
    Dfout.flux.dQHv_dxs = (dnum_dxs.*den-num.*dden_dxs)./den.^2;
    Dfout.flux.dQHv_dxt = (dnum_dxt.*den-num.*dden_dxt)./den.^2;
       
end


end


