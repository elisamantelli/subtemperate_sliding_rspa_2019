function [fout, Dfout] = hor_velocity_SIA_v2(v_in, parameters)
% compute shallow ice velocity at velocity nodes, uses velocity network to
% compute mass flux by shearing, horizontal advective heat flux through
% shearing velocity, strain heating, shearing velocity at the grounding
% line (along hot T edges), and the shearing component of int_0^(eta') uh d eta'.
%If nargout == 2 produces also jacobian on net fluxes, already corrected for groiunding line outflow, strain heating
% and grounding line velocity with respect to the concatenated vector v_in. 

%input: concatenated vector with H, T at each network node, and xc, xs, xt

%unpack parameters
u_n_nodes = parameters.grid_u.n_nodes.tot;            %number of u nodes
u_nodes_ver = parameters.grid_u.n_nodes.ver;          %number of nodes in the vertical direction
u_up_node_vert = parameters.grid_u.up_node.vert;      %list (n_edges-by-1 vector) of 'upstream' node for each vertical edge,u
u_down_node_vert = parameters.grid_u.down_node.vert;  %list (n_edges-by-1 vector) of 'downstream' nodes for each vertical edge,u
u_eta_nodes = parameters.grid_u.coor_nodes.eta;       %vertical coordinates of nodes (n_nodes-by-one vector),u
Delta_eta = parameters.grid_u.Delta_eta;              %vertical spacing of u nodes

h_nodes = parameters.grid_h.n_nodes;                  %number of nodes,h
h_n_edges = parameters.grid_h.n_edges;                %number of edges in the h network 
h_up_node = parameters.grid_h.up_node;                %list (h_n_edge-by-1 vector) of 'upstream' node, h
h_down_node = parameters.grid_h.down_node;            %list (h_n_edges-by-1 vector) of 'downstream' node,h

T_nodes = parameters.grid_T.n_nodes.tot;              %number of temperature nodes
nodes_ver = parameters.grid_T.n_nodes.vert;           %number of nodes in the vertical direction
T_n_edges_hor = parameters.grid_T.n_edges.hor;        %number of horizontal edges, T
T_up_node_hor = parameters.grid_T.up_node.hor;        %list (n_edges-by-1 vector) of 'upstream' node for each horizontal edge,T
T_down_node_hor = parameters.grid_T.down_node.hor;    %list (n_edges-by-1 vector) of 'downstream' nodes for each horizontal edge,T

n_x_c = parameters.n_x.c;                             %number of horizontal nodes, cold subd
n_x_s = parameters.n_x.s;                             %number of horizontal nodes, subt subd
n_x_t = parameters.n_x.t;                             %number of horizontal nodes, temp subd

%unpack discretized variables
[discvar, Ddiscvar] = discretisation(v_in, parameters);

H_edge = discvar.H.H_edge;
dHedgedHup = Ddiscvar.H.dHedge_dHup;
dHedgedHdown = Ddiscvar.H.dHedge_dHdown;
dHedgedxc = Ddiscvar.H.dHedge_dxc;
dHedgedxs = Ddiscvar.H.dHedge_dxs;
dHedgedxt = Ddiscvar.H.dHedge_dxt;

dS_edge = discvar.dS.dS_edge;
ddS_edgedHup = Ddiscvar.dS.ddSedge_dHup;
ddS_edgedHdown = Ddiscvar.dS.ddSedge_dHdown;
ddSedge_dxc = Ddiscvar.dS.ddSedge_dxc;
ddSedge_dxs = Ddiscvar.dS.ddSedge_dxs;
ddSedge_dxt = Ddiscvar.dS.ddSedge_dxt;

Hg = discvar.H.H_g;
dHg_dH = Ddiscvar.H.dHg_dH;
dSg = discvar.dS.dS_g;
ddSg_dH = Ddiscvar.dS.ddSg_dH;
ddSg_dxc = Ddiscvar.dS.ddSg_dxc;
ddSg_dxs = Ddiscvar.dS.ddSg_dxs;
ddSg_dxt = Ddiscvar.dS.ddSg_dxt;

T_hor = discvar.T.T_hor;
dThor_dTup = Ddiscvar.T.dThor_dTup;
dThor_dTdown = Ddiscvar.T.dThor_dTdown;

T_hor_g = discvar.T.T_hor_g ;
dThorg_dTup = Ddiscvar.T.dThorg_dT_up;

%reindexings
index_h_to_u = reshape(repmat((1:length(H_edge)),[u_nodes_ver,1]),[u_n_nodes,1]);
index_h_to_T = reshape(repmat((1:length(H_edge)),[nodes_ver,1]),[T_n_edges_hor,1]);
index_h_to_u_Thor = reshape(repmat((1:length(H_edge)),[nodes_ver,1]),[T_n_edges_hor,1]);

index_h_to_u_g = ones(u_nodes_ver,1);
u_up_ver_g = 1:u_nodes_ver-1;
u_down_ver_g = 2:u_nodes_ver;

%horizontal length of cells
length_T = fvlength(v_in(end-2:end), parameters,'T');
length_h = fvlength(v_in(end-2:end), parameters,'h');
Delta_sigma_T = length_T. Delta_sigma;
Delta_sigma_h = length_h. Delta_sigma;
%SHEARING VELOCITY AT T HOR EDGES compute u_{j+1/2,k}=1/2 (u_{j+1/2,k+1/2}+u_{j+1/2,k-1/2}), with jk indexes
%in T network

%compute shear velocity at u nodes;
H_edge_u = H_edge(index_h_to_u);
dS_edge_u = dS_edge(index_h_to_u);
U_unodes =  -1/2*(1-(1-u_eta_nodes).^2).*(H_edge_u.^2.*dS_edge_u);
UH_unodes = U_unodes.*H_edge_u;
tau_unodes = - (1-u_eta_nodes).*dS_edge_u;% this is really tau/h. Make it clear in the documentation

%grounding line
Hg_unodes = Hg(index_h_to_u_g);
dSg_unodes = dSg(index_h_to_u_g);
U_unodes_g = -1/2*(1-(1-u_eta_nodes(1:u_nodes_ver)).^2).*(Hg_unodes.^2.*dSg_unodes);

%reindexing to T hor edges and split grounding line edges
U_unode_up_ver = U_unodes(u_up_node_vert);
U_unode_down_ver = U_unodes(u_down_node_vert);
UH_unodes_up_ver = UH_unodes(u_up_node_vert);
UH_unodes_down_ver = UH_unodes(u_down_node_vert);
tau_unodes_up = tau_unodes(u_up_node_vert);
tau_unodes_down = tau_unodes(u_down_node_vert);

U_unodes_g_up = U_unodes_g(u_up_ver_g);
U_unodes_g_down = U_unodes_g(u_down_ver_g);

if nargout == 2
    dUn_dHe = -1/2*(1-(1-u_eta_nodes).^2).*(2*H_edge_u.*dS_edge_u);
    dUn_ddSe = -1/2*(1-(1-u_eta_nodes).^2).*(H_edge_u.^2);
    dUn_dxc = dUn_ddSe.*ddSedge_dxc(index_h_to_u)+dUn_dHe.*dHedgedxc(index_h_to_u);
    dUn_dxs = dUn_ddSe.*ddSedge_dxs(index_h_to_u)+dUn_dHe.*dHedgedxs(index_h_to_u);
    dUn_dxt = dUn_ddSe.*ddSedge_dxt(index_h_to_u)+dUn_dHe.*dHedgedxt(index_h_to_u);
    
    indexcomposite_hu = h_up_node(index_h_to_u);
    indexcomposite_hd = h_down_node(index_h_to_u);    

    dtaun_dHup = - (1-u_eta_nodes).*ddS_edgedHup(index_h_to_u);
    dtaun_dHdown = - (1-u_eta_nodes).*ddS_edgedHdown(index_h_to_u);
    dtaun_dxc = - (1-u_eta_nodes).*ddSedge_dxc(index_h_to_u);
    dtaun_dxs = - (1-u_eta_nodes).*ddSedge_dxs(index_h_to_u);
    dtaun_dxt = - (1-u_eta_nodes).*ddSedge_dxt(index_h_to_u);
    dUHn_dHu = U_unodes.*dHedgedHup(index_h_to_u) + H_edge(index_h_to_u).* (dUn_dHe.*dHedgedHup(index_h_to_u) + dUn_ddSe.*ddS_edgedHup(index_h_to_u));
    dUHn_dHd = U_unodes.*dHedgedHdown(index_h_to_u) + H_edge(index_h_to_u).* (dUn_dHe.*dHedgedHdown(index_h_to_u)+ dUn_ddSe.*ddS_edgedHdown(index_h_to_u));
    dUHn_dxc = H_edge_u.*dUn_dxc + U_unodes.*dHedgedxc(index_h_to_u);
    dUHn_dxs = H_edge_u.*dUn_dxs + U_unodes.*dHedgedxs(index_h_to_u);
    dUHn_dxt = H_edge_u.*dUn_dxt + U_unodes.*dHedgedxt(index_h_to_u);

    dUng_dHg = -1/2*(1-(1-u_eta_nodes(1:u_nodes_ver)).^2).*(2*Hg_unodes.*dSg_unodes);
    dUng_ddSg = -1/2*(1-(1-u_eta_nodes(1:u_nodes_ver)).^2).*(Hg_unodes.^2);
    dUng_dxc = dUng_ddSg.*ddSg_dxc(index_h_to_u_g);
    dUng_dxs = dUng_ddSg.*ddSg_dxs(index_h_to_u_g);
    dUng_dxt = dUng_ddSg.*ddSg_dxt(index_h_to_u_g);
end
%n.b. what here is called tau is actually tau/h
U_Thoredges = 1/2*(U_unode_up_ver + U_unode_down_ver);
UH_Thoredges = U_Thoredges.*H_edge(index_h_to_T);
dUH_deta = 1/Delta_eta*(UH_unodes_down_ver-UH_unodes_up_ver); % dUH/deta_{j+1/2, k} = (UH_{j+1/2,k+1/2}-UH_{j+1/2,k-1/2})/Delta_eta
tau_Thoredges = 1/2*(tau_unodes_up+tau_unodes_down); %Tau_{j+1/2,k}=1/2 (Tau_{j+1/2,k+1/2}+Tau_{j+1/2,k-1/2})
S_Thoredges = dUH_deta.*tau_Thoredges;

fout.U_Thoredges = U_Thoredges;
U_Thoredges_g = 1/2*(U_unodes_g_up+U_unodes_g_down);
fout.U_g = U_Thoredges_g;


if nargout == 2
    dUThoredges_dHedge = 1/2*dUn_dHe(u_up_node_vert)+1/2*dUn_dHe(u_down_node_vert);
    dUThoredges_ddSedge = 1/2*dUn_ddSe(u_up_node_vert)+1/2*dUn_ddSe(u_down_node_vert);
    dUThoredges_dxc = 1/2*dUn_dxc(u_up_node_vert)+1/2*dUn_dxc(u_down_node_vert);
    dUThoredges_dxs = 1/2*dUn_dxs(u_up_node_vert)+1/2*dUn_dxs(u_down_node_vert);
    dUThoredges_dxt = 1/2*dUn_dxt(u_up_node_vert)+1/2*dUn_dxt(u_down_node_vert);
    
    dUHThoredges_dHedge = dUThoredges_dHedge.*H_edge(index_h_to_T) + U_Thoredges;
    dUHThoredges_ddSedge = dUThoredges_ddSedge.*H_edge(index_h_to_T);
    dUHThoredge_dHup = dUHThoredges_dHedge.*dHedgedHup(index_h_to_T) + dUHThoredges_ddSedge.*ddS_edgedHup(index_h_to_T) ;
    dUHThoredge_dHdown = dUHThoredges_dHedge.*dHedgedHdown(index_h_to_T) + dUHThoredges_ddSedge.*ddS_edgedHdown(index_h_to_T);
    dUHThoredge_dxc = dUThoredges_dxc.*H_edge(index_h_to_T) + U_Thoredges.*dHedgedxc(index_h_to_T);
    dUHThoredge_dxs = dUThoredges_dxs.*H_edge(index_h_to_T)+ U_Thoredges.*dHedgedxs(index_h_to_T);
    dUHThoredge_dxt = dUThoredges_dxt.*H_edge(index_h_to_T)+ U_Thoredges.*dHedgedxt(index_h_to_T);
    
    indexcomposite_2_u = indexcomposite_hu(u_up_node_vert);
    indexcomposite_2_d = indexcomposite_hd(u_down_node_vert);
    ddUHdeta_dHup = 1/Delta_eta*(dUHn_dHu(u_down_node_vert)-dUHn_dHu(u_up_node_vert));
    ddUHdeta_dHdown = 1/Delta_eta*(dUHn_dHd(u_down_node_vert)-dUHn_dHd(u_up_node_vert));
    ddUHdeta_dxc = 1/Delta_eta*(dUHn_dxc(u_down_node_vert)-dUHn_dxc(u_up_node_vert));
    ddUHdeta_dxs = 1/Delta_eta*(dUHn_dxs(u_down_node_vert)-dUHn_dxs(u_up_node_vert));
    ddUHdeta_dxt = 1/Delta_eta*(dUHn_dxt(u_down_node_vert)-dUHn_dxt(u_up_node_vert));
    
    dtauThoredges_dHup = 1/2*dtaun_dHup(u_up_node_vert)+1/2*dtaun_dHup(u_down_node_vert)  ;
    dtauThoredges_dHdown = 1/2*dtaun_dHdown(u_up_node_vert)+1/2*dtaun_dHdown(u_down_node_vert);
    dtauThoredges_dxc = 1/2*dtaun_dxc(u_up_node_vert)+1/2*dtaun_dxc(u_down_node_vert);
    dtauThoredges_dxs = 1/2*dtaun_dxs(u_up_node_vert)+1/2*dtaun_dxs(u_down_node_vert);
    dtauThoredges_dxt = 1/2*dtaun_dxt(u_up_node_vert)+1/2*dtaun_dxt(u_down_node_vert);
    
    dSThoredges_dHu = dUH_deta.*dtauThoredges_dHup + tau_Thoredges.* ddUHdeta_dHup;
    dSThoredges_dHd = dUH_deta.*dtauThoredges_dHdown+ tau_Thoredges.* ddUHdeta_dHdown;
    dSThoredges_dxc = dUH_deta.*dtauThoredges_dxc + tau_Thoredges.*ddUHdeta_dxc;
    dSThoredges_dxs = dUH_deta.*dtauThoredges_dxs + tau_Thoredges.*ddUHdeta_dxs;
    dSThoredges_dxt = dUH_deta.*dtauThoredges_dxt + tau_Thoredges.*ddUHdeta_dxt;
    
    %grounding line  
    dUThoredgesg_dHg = 1/2*dUng_dHg(u_up_ver_g)+1/2*dUng_dHg(u_down_ver_g);
    dUThoredgesg_ddSg =  1/2*dUng_ddSg(u_up_ver_g) +1/2*dUng_ddSg(u_down_ver_g);
    dUThoredgesg_dxc = 1/2*dUng_dxc(u_up_ver_g)+1/2*dUng_dxc(u_down_ver_g);
    dUThoredgesg_dxs = 1/2*dUng_dxs(u_up_ver_g)+1/2*dUng_dxs(u_down_ver_g);
    dUThoredgesg_dxt = 1/2*dUng_dxt(u_up_ver_g)+1/2*dUng_dxt(u_down_ver_g);

    Dfout.dUg_dHg = dUThoredgesg_dHg;
    Dfout.dUg_ddSg = dUThoredgesg_ddSg;
    Dfout.dUg_dxc = dUThoredgesg_dxc;
    Dfout.dUg_dxs = dUThoredgesg_dxs;
    Dfout.dUg_dxt = dUThoredgesg_dxt;
end

%ADVECTIVE HEAT FLUX AND ITS DERIVATIVES
fout.QHa_shear = U_Thoredges.* H_edge(index_h_to_T).*T_hor;
fout.QHa_shear_g = U_Thoredges_g* Hg.*T_hor_g;
if nargout == 2
    %temperature
    dQHa_shear_dT_u = U_Thoredges.*H_edge(index_h_to_T).*dThor_dTup;
    dQHa_shear_dT_d = U_Thoredges.*H_edge(index_h_to_T).*dThor_dTdown;
    dQHa_shear_dT = sparse(T_up_node_hor,T_up_node_hor,dQHa_shear_dT_u,T_nodes, T_nodes)+...
        sparse(T_up_node_hor,T_down_node_hor,dQHa_shear_dT_d ,T_nodes, T_nodes)-...
        (sparse(T_down_node_hor,T_up_node_hor,dQHa_shear_dT_u,T_nodes, T_nodes)+...
        sparse(T_down_node_hor, T_down_node_hor,dQHa_shear_dT_d ,T_nodes, T_nodes));

    dQHa_shearg_dT = Hg* sparse(1:nodes_ver,1:nodes_ver, U_Thoredges_g, nodes_ver, nodes_ver)*dThorg_dTup;
    %ice thickness
    indexcomposite_u = h_up_node(index_h_to_T);
    indexcomposite_d = h_down_node(index_h_to_T);
    
    dQHa_shear_dH_u = U_Thoredges.*T_hor.*dHedgedHup(index_h_to_T) + ...
        T_hor.*H_edge(index_h_to_T).*(dUThoredges_dHedge.*dHedgedHup(index_h_to_T)+dUThoredges_ddSedge.*ddS_edgedHup(index_h_to_T)) ;
    dQHa_shear_dH_d = U_Thoredges.*T_hor.*dHedgedHdown(index_h_to_T)+...
        T_hor.*H_edge(index_h_to_T).*(dUThoredges_dHedge.*dHedgedHdown(index_h_to_T)+dUThoredges_ddSedge.*ddS_edgedHdown(index_h_to_T)) ;
    dQHa_shear_dH = sparse(T_up_node_hor,indexcomposite_u,dQHa_shear_dH_u,T_nodes, h_nodes)+...
        sparse(T_up_node_hor,indexcomposite_d,dQHa_shear_dH_d ,T_nodes, h_nodes)-...
        (sparse(T_down_node_hor,indexcomposite_u,dQHa_shear_dH_u,T_nodes, h_nodes)+...
        sparse(T_down_node_hor, indexcomposite_d,dQHa_shear_dH_d ,T_nodes, h_nodes));
    dQHa_shearg_dH = diag(U_Thoredges_g.*T_hor_g)*repmat(dHg_dH', [nodes_ver,1])+diag(Hg*T_hor_g.*dUThoredgesg_dHg)*repmat(dHg_dH', [nodes_ver,1])+...
        diag(Hg*T_hor_g.*dUThoredgesg_ddSg)*repmat(ddSg_dH', [nodes_ver,1]);
    %length of subdomains
    dQHa_shear_dxc = (accumarray(T_up_node_hor, (dUThoredges_dxc.*H_edge(index_h_to_T)+U_Thoredges.*dHedgedxc(index_h_to_T)).*T_hor, [T_nodes,1])...
        -accumarray(T_down_node_hor, (dUThoredges_dxc.*H_edge(index_h_to_T)+U_Thoredges.*dHedgedxc(index_h_to_T)).*T_hor, [T_nodes,1]));
    dQHa_shear_dxs = (accumarray(T_up_node_hor, (dUThoredges_dxs.*H_edge(index_h_to_T)+U_Thoredges.*dHedgedxs(index_h_to_T)).*T_hor, [T_nodes,1])...
        -accumarray(T_down_node_hor, (dUThoredges_dxs.*H_edge(index_h_to_T)+U_Thoredges.*dHedgedxs(index_h_to_T)).*T_hor, [T_nodes,1]));
    dQHa_shear_dxt = (accumarray(T_up_node_hor, (dUThoredges_dxt.*H_edge(index_h_to_T)+U_Thoredges.*dHedgedxt(index_h_to_T)).*T_hor, [T_nodes,1])...
        -accumarray(T_down_node_hor, (dUThoredges_dxt.*H_edge(index_h_to_T)+U_Thoredges.*dHedgedxt(index_h_to_T)).*T_hor, [T_nodes,1]));
    dQHa_shearg_dxc = Hg* T_hor_g.*dUThoredgesg_dxc;
    dQHa_shearg_dxs = Hg* T_hor_g.*dUThoredgesg_dxs;
    dQHa_shearg_dxt = Hg* T_hor_g.*dUThoredgesg_dxt;
    %correct for grounding line
    dQHa_shear_dT(end-nodes_ver+1:end, : ) =  dQHa_shear_dT(end-nodes_ver+1:end, : ) + dQHa_shearg_dT;
    dQHa_shear_dH(end-nodes_ver+1:end, : ) =  dQHa_shear_dH(end-nodes_ver+1:end, : ) + dQHa_shearg_dH;
    dQHa_shear_dxc(end-nodes_ver+1:end ) =  dQHa_shear_dxc(end-nodes_ver+1:end ) + dQHa_shearg_dxc;
    dQHa_shear_dxs(end-nodes_ver+1:end ) =  dQHa_shear_dxs(end-nodes_ver+1:end ) + dQHa_shearg_dxs;
    dQHa_shear_dxt(end-nodes_ver+1:end ) =  dQHa_shear_dxt(end-nodes_ver+1:end ) + dQHa_shearg_dxt;
    %assemble output
    Dfout.dQHa_shear_dT_net = dQHa_shear_dT ;
    Dfout.dQHa_shear_dH_net = dQHa_shear_dH;
    Dfout.dQHa_shear_dxc_net = dQHa_shear_dxc;
    Dfout.dQHa_shear_dxs_net = dQHa_shear_dxs;
    Dfout.dQHa_shear_dxt_net = dQHa_shear_dxt;
end

%STRAIN HEATING AND DERIVATIVES
T_up_edge_hor = 1:length(S_Thoredges)-nodes_ver;
T_down_edge_hor = nodes_ver+1:length(S_Thoredges);
S_up = S_Thoredges(T_up_edge_hor);
S_down = S_Thoredges(T_down_edge_hor);
S = 1/2*(S_up+S_down);
%correct for grounding line (upwinding from upstream edge) and divide(
%averaging with divide edge, where no flux c. applies)
fout.S = [ 1/2*S_up(1:nodes_ver); S; S_down(end-nodes_ver+1:end)];

if nargout == 2
    dSu_dHu = dSThoredges_dHu(T_up_edge_hor);
    dSu_dHd = dSThoredges_dHd(T_up_edge_hor);
    dSd_dHu = dSThoredges_dHu(T_down_edge_hor);
    dSd_dHd = dSThoredges_dHd(T_down_edge_hor);
    
    indexcomposite_3_uu = indexcomposite_2_u(T_up_edge_hor);
    indexcomposite_3_ud = indexcomposite_2_d(T_up_edge_hor);
    indexcomposite_3_du = indexcomposite_2_u(T_down_edge_hor);
    indexcomposite_3_dd = indexcomposite_2_d(T_down_edge_hor);
    indexrow = T_up_edge_hor+nodes_ver;
    
    dS_dH = sparse(indexrow,indexcomposite_3_uu,1/2*dSu_dHu,T_nodes,h_nodes )+...
        sparse(indexrow,indexcomposite_3_ud,1/2*dSu_dHd,T_nodes,h_nodes )+...
        sparse(indexrow,indexcomposite_3_du,1/2*dSd_dHu,T_nodes,h_nodes )+...
        sparse(indexrow,indexcomposite_3_dd,1/2*dSd_dHd,T_nodes,h_nodes );
    dS_dxc = 1/2*(dSThoredges_dxc(T_up_edge_hor)+ dSThoredges_dxc(T_down_edge_hor));
    dS_dxs = 1/2*(dSThoredges_dxs(T_up_edge_hor)+ dSThoredges_dxs(T_down_edge_hor));
    dS_dxt = 1/2*(dSThoredges_dxt(T_up_edge_hor)+ dSThoredges_dxt(T_down_edge_hor));
    %divide
    indexrow_d = 1:nodes_ver;
    indexcomposite_3_du_div = indexcomposite_3_uu(indexrow_d);
    indexcomposite_3_dd_div = indexcomposite_3_ud(indexrow_d);
    dS_dHu_div = dSu_dHu(indexrow_d);
    dS_dHd_div = dSu_dHd(indexrow_d);
    dS_dH_div = sparse(indexrow_d,indexcomposite_3_du_div,1/2*dS_dHu_div,T_nodes,h_nodes )+...
        sparse(indexrow_d,indexcomposite_3_dd_div,1/2*dS_dHd_div,T_nodes,h_nodes );
    dS_dxc_div = 1/2*dSThoredges_dxc(T_up_edge_hor(indexrow_d));
    dS_dxs_div = 1/2*dSThoredges_dxs(T_up_edge_hor(indexrow_d));
    dS_dxt_div = 1/2*dSThoredges_dxt(T_up_edge_hor(indexrow_d));
    %grounding line: upwinding from upstream EDGE
    indexrow_g = T_nodes-nodes_ver+1:T_nodes;
    indexcomposite_3_uu_g = indexcomposite_3_du(end-nodes_ver+1:end);
    indexcomposite_3_ud_g = indexcomposite_3_dd(end-nodes_ver+1:end);
    dSu_dHu_g = dSd_dHu(end-nodes_ver+1:end);
    dSu_dHd_g = dSd_dHd(end-nodes_ver+1:end);
    dS_dH_g = sparse(indexrow_g,indexcomposite_3_uu_g,dSu_dHu_g,T_nodes,h_nodes )+...
        sparse(indexrow_g,indexcomposite_3_ud_g,dSu_dHd_g,T_nodes,h_nodes );
    dS_dxc_g = dSThoredges_dxc(T_down_edge_hor(end-nodes_ver+1:end));
    dS_dxs_g = dSThoredges_dxs(T_down_edge_hor(end-nodes_ver+1:end));
    dS_dxt_g = dSThoredges_dxt(T_down_edge_hor(end-nodes_ver+1:end));
    %assemble final jacobian
    Dfout.dS_dH = dS_dH + dS_dH_div + dS_dH_g;
    Dfout.dS_dxc = [dS_dxc_div; dS_dxc; dS_dxc_g];
    Dfout.dS_dxs = [dS_dxs_div; dS_dxs; dS_dxs_g];
    Dfout.dS_dxt = [dS_dxt_div; dS_dxt; dS_dxt_g];
end

%DIFFUSIVE MASS FLUX AT H HOR EDGES AND ITS DERIVATIVES
fout.QMd = Delta_eta*accumarray (index_h_to_u_Thor, U_Thoredges, [h_n_edges,1]).* H_edge;
fout.QMd_g = Delta_eta*sum(U_Thoredges_g* Hg);

if nargout == 2
    dQMd_dHedge = Delta_eta*accumarray (index_h_to_T, U_Thoredges, [h_n_edges,1]) + Delta_eta*accumarray (index_h_to_T, dUThoredges_dHedge, [h_n_edges,1]).* H_edge;
    dQMd_ddSedge = Delta_eta*accumarray (index_h_to_T, dUThoredges_ddSedge, [h_n_edges,1]).* H_edge;
    dQMd_dxc = Delta_eta*(accumarray (index_h_to_T, dUThoredges_dxc, [h_n_edges,1]).* H_edge + accumarray (index_h_to_u_Thor, U_Thoredges, [h_n_edges,1]).* dHedgedxc);
    dQMd_dxs = Delta_eta*(accumarray (index_h_to_T, dUThoredges_dxs, [h_n_edges,1]).* H_edge + accumarray (index_h_to_u_Thor, U_Thoredges, [h_n_edges,1]).* dHedgedxs);
    dQMd_dxt = Delta_eta*(accumarray (index_h_to_T, dUThoredges_dxt, [h_n_edges,1]).* H_edge + accumarray (index_h_to_u_Thor, U_Thoredges, [h_n_edges,1]).* dHedgedxt);
    
    dQMd_dH_u = dQMd_dHedge.*dHedgedHup + dQMd_ddSedge.*ddS_edgedHup;
    dQMd_dH_d = dQMd_dHedge.*dHedgedHdown + dQMd_ddSedge.*ddS_edgedHdown;
    
    dQMd_dH = sparse(h_up_node,h_up_node,dQMd_dH_u,h_nodes, h_nodes)+...
        sparse(h_up_node,h_down_node,dQMd_dH_d ,h_nodes, h_nodes)-...
        (sparse(h_down_node,h_up_node,dQMd_dH_u,h_nodes, h_nodes)+...
        sparse(h_down_node, h_down_node,dQMd_dH_d ,h_nodes, h_nodes));
    dQMd_dxc = (accumarray(h_up_node, dQMd_dxc, [h_nodes,1])-accumarray(h_down_node, dQMd_dxc, [h_nodes,1]));
    dQMd_dxs = (accumarray(h_up_node, dQMd_dxs, [h_nodes,1])-accumarray(h_down_node, dQMd_dxs, [h_nodes,1]));
    dQMd_dxt = (accumarray(h_up_node, dQMd_dxt, [h_nodes,1])-accumarray(h_down_node, dQMd_dxt, [h_nodes,1]));
    
    dQMddown_dHdown = sparse(h_up_node,h_down_node,dQMd_dH_d ,h_nodes, h_nodes);
    dQMddown_dHup = sparse(h_up_node,h_up_node,dQMd_dH_u ,h_nodes, h_nodes);
    %grounding line
    dQMdg_dH = Delta_eta*sum(U_Thoredges_g)*dHg_dH + Delta_eta*sum(dUThoredgesg_ddSg* Hg)*ddSg_dH + Delta_eta*sum(dUThoredgesg_dHg* Hg)*dHg_dH;
    dQMdg_dxc = Delta_eta*sum(dUThoredgesg_ddSg* Hg)*ddSg_dxc;
    dQMdg_dxs = Delta_eta*sum(dUThoredgesg_ddSg* Hg)*ddSg_dxs;
    dQMdg_dxt = Delta_eta*sum(dUThoredgesg_ddSg* Hg)*ddSg_dxt;
    %correct for grounding line
    dQMd_dH(end,:) = dQMd_dH(end,:) + dQMdg_dH.';
    dQMd_dxc(end) = dQMd_dxc(end) + dQMdg_dxc;
    dQMd_dxs(end) = dQMd_dxs(end) + dQMdg_dxs;
    dQMd_dxt(end) = dQMd_dxt(end) + dQMdg_dxt;
    
    dQMddown_dHdown(end,:) = dQMddown_dHdown(end,:) + dQMdg_dH.';
    %assemble output
    Dfout.dQMd_dH_net = dQMd_dH;
    Dfout.dQMd_dxc_net = dQMd_dxc;
    Dfout.dQMd_dxs_net = dQMd_dxs;
    Dfout.dQMd_dxt_net = dQMd_dxt;
    Dfout.dQMdg_dH = dQMdg_dH;
    Dfout.dQMdg_dxc = dQMdg_dxc;
    Dfout.dQMdg_dxs = dQMdg_dxs;
    Dfout.dQMdg_dxt = dQMdg_dxt;
    Dfout.dQMddown_dHdown = dQMddown_dHdown;
    Dfout.dQMddown_dHup = dQMddown_dHup;
end

% FLUX (eta) = INT_0^eta u_shear(eta') deta' FOR VERTICAL VELOCITY SOLVER
UHThoredges_up = UH_Thoredges(T_up_edge_hor);
UHThoredges_down = UH_Thoredges(T_down_edge_hor);
net_flux_shear_eta = Delta_eta./Delta_sigma_T.*[UH_Thoredges(1:nodes_ver) ;...
    (UHThoredges_down-UHThoredges_up);U_Thoredges_g*Hg - UH_Thoredges(end-nodes_ver+1:end)];
%sum along the vertical 
ind_sum = kron(speye(n_x_c+n_x_s+n_x_t,n_x_c+n_x_s+n_x_t), sparse(tril(ones(nodes_ver))));
fout.net_flux_shear_eta_sum = ind_sum*net_flux_shear_eta;
if nargout == 2
    dfseu_dHu = dUHThoredge_dHup(T_up_edge_hor);
    dfseu_dHd = dUHThoredge_dHdown(T_up_edge_hor);
    dfsed_dHu = dUHThoredge_dHup(T_down_edge_hor);
    dfsed_dHd = dUHThoredge_dHdown(T_down_edge_hor);
    
    dnfse_dH = sparse(indexrow,indexcomposite_3_du,dfsed_dHu,T_nodes,h_nodes )+...
        sparse(indexrow,indexcomposite_3_dd,dfsed_dHd,T_nodes,h_nodes )-...
        (sparse(indexrow,indexcomposite_3_uu,dfseu_dHu,T_nodes,h_nodes )+...
        sparse(indexrow,indexcomposite_3_ud,dfseu_dHd,T_nodes,h_nodes ));
    dnfse_dxc = dUHThoredge_dxc(T_down_edge_hor)-dUHThoredge_dxc(T_up_edge_hor);
    dnfse_dxs = dUHThoredge_dxs(T_down_edge_hor)-dUHThoredge_dxs(T_up_edge_hor);
    dnfse_dxt = dUHThoredge_dxt(T_down_edge_hor)-dUHThoredge_dxt(T_up_edge_hor);
    %divide
    dfsed_dHu_div = dfseu_dHu(indexrow_d);
    dfsed_dHd_div = dfseu_dHd(indexrow_d);
    dnfse_dH_div = sparse(indexrow_d,indexcomposite_3_du_div,dfsed_dHu_div,T_nodes,h_nodes )+...
        sparse(indexrow_d,indexcomposite_3_dd_div,dfsed_dHd_div,T_nodes,h_nodes );
    dnfse_dxc_div = dUHThoredge_dxc(T_up_edge_hor(indexrow_d));
    dnfse_dxs_div = dUHThoredge_dxs(T_up_edge_hor(indexrow_d));
    dnfse_dxt_div = dUHThoredge_dxt(T_up_edge_hor(indexrow_d));
    %grounding line
    dfseu_dHu_g = dfsed_dHu(end-nodes_ver+1:end);
    dfseu_dHd_g = dfsed_dHd(end-nodes_ver+1:end);
    dfsed_dH = sparse([sparse((T_nodes-nodes_ver),h_nodes);
        diag(U_Thoredges_g+Hg*dUThoredgesg_dHg )*repmat(dHg_dH', [nodes_ver,1])]) +...
    sparse([sparse((T_nodes-nodes_ver),h_nodes);
        diag(Hg*dUThoredgesg_ddSg )*repmat(ddSg_dH', [nodes_ver,1])]);
    dnfse_dH_g = dfsed_dH-(sparse(indexrow_g,indexcomposite_3_uu_g,dfseu_dHu_g,T_nodes,h_nodes )+...
        sparse(indexrow_g,indexcomposite_3_ud_g,dfseu_dHd_g,T_nodes,h_nodes ));
    dnfse_dxc_g = dUThoredgesg_dxc*Hg - dUHThoredge_dxc(T_down_edge_hor(end-nodes_ver+1:end));
    dnfse_dxs_g = dUThoredgesg_dxs*Hg - dUHThoredge_dxs(T_down_edge_hor(end-nodes_ver+1:end));
    dnfse_dxt_g = dUThoredgesg_dxt*Hg - dUHThoredge_dxt(T_down_edge_hor(end-nodes_ver+1:end));
    %assemble full jacobian accounting for grounding line
    dnfse_dH = sparse(1:T_nodes, 1:T_nodes,Delta_eta./Delta_sigma_T,T_nodes,T_nodes)*(dnfse_dH + dnfse_dH_div + dnfse_dH_g);
    dnfse_dxc =  Delta_eta./Delta_sigma_T.*[dnfse_dxc_div; dnfse_dxc; dnfse_dxc_g];
    dnfse_dxs =  Delta_eta./Delta_sigma_T.*[dnfse_dxs_div; dnfse_dxs; dnfse_dxs_g];
    dnfse_dxt =  Delta_eta./Delta_sigma_T.*[dnfse_dxt_div; dnfse_dxt; dnfse_dxt_g];
    %sum along the vertical
    Dfout.dnfse_sum_dH = ind_sum*dnfse_dH;
    Dfout.dnfse_sum_dxc = ind_sum*dnfse_dxc;
    Dfout.dnfse_sum_dxs = ind_sum*dnfse_dxs;
    Dfout.dnfse_sum_dxt = ind_sum*dnfse_dxt;
end
end










