function [fout, Dfout] = sliding_v3(v_in,parameters)
%returns advective mass flux and horizontal advective heat flux not
%corrected for grounding line outflow, plus subtemperate sliding velocity
%and temperate sliding velocity
%at the subtemperate-temperate boundary. If nargout == 2, function also returns jacobian
%of net fluxes (not corrected for grounding line outflow) and slide velocity at xs with respect to v_in

%unpack parameters
h_nodes = parameters.grid_h.n_nodes;                  %number of nodes,h
h_n_edges = parameters.grid_h.n_edges;                %number of horizontal edges,h
h_up_node = parameters.grid_h.up_node;                %list (h_n_edge-by-1 vector) of 'upstream' node, h
h_down_node = parameters.grid_h.down_node;            %list (h_n_edges-by-1 vector) of 'downstream' node,h
h_id_node = parameters.grid_h.id_node;                %list (n_nodes_1D-by-one vector) of identifiers of subdomain for the 1D grid

T_nodes = parameters.grid_T.n_nodes.tot;                                   %number of temperature nodes
T_id_node = parameters.grid_T.id_node;                                     %list (n_T_nodes-by-one vector) of identifiers of subdomain for the T grid
nodes_ver = parameters.grid_T.n_nodes.vert;                                %number of nodes in the vertical direction
T_n_edges_hor = parameters.grid_T.n_edges.hor;                             %number of horizontal edges, T
T_up_node_hor = parameters.grid_T.up_node.hor;                             %list (n_edges-by-1 vector) of 'upstream' node for each horizontal edge,T
T_down_node_hor = parameters.grid_T.down_node.hor;                        %list (n_edges-by-1 vector) of 'downstream' nodes for each horizontal edge,T

gamma = parameters.gamma;   %friction coefficient

flag = parameters.flag.u_sl;

%interpolation
[discvar, Ddiscvar] = discretisation_v2(v_in, parameters);

%ICE THICKNESS AND DRIVING STRESS
%along hor edges
H_edge = discvar.H.H_edge;
dHedgedHup = Ddiscvar.H.dHedge_dHup;
dHedgedHdown = Ddiscvar.H.dHedge_dHdown;
dHedgedxc = Ddiscvar.H.dHedge_dxc;
dHedgedxs = Ddiscvar.H.dHedge_dxs;
dHedgedxt = Ddiscvar.H.dHedge_dxt;

H_edge_adv = discvar.H.H_edge_adv;
dHedgeadv_dHup = Ddiscvar.H.dHedgea_dHup;
dHedgeadv_dHdown = Ddiscvar.H.dHedgea_dHdown;
dS_edge = discvar.dS.dS_edge;
ddS_edgedHup = Ddiscvar.dS.ddSedge_dHup;
ddS_edgedHdown = Ddiscvar.dS.ddSedge_dHdown;

%TEMPERATURE
%along T edges
T_hor = discvar.T.T_hor;
dThor_dTup = Ddiscvar.T.dThor_dTup;
dThor_dTdown = Ddiscvar.T.dThor_dTdown;
T_edge = discvar.T.T_hor;

%SLIDING VELOCITY
U_SL = sparse(h_n_edges,1);
dUSL_dxc = sparse(h_n_edges,1);
dUSL_dxs = sparse(h_n_edges,1);
dUSL_dxt = sparse(h_n_edges,1);

%COLD SUB_DOMAIN
U_SL(h_id_node(h_up_node) == 1) = 0;

%SUBTEMP SUB-DOMAIN
%distinguish between different schemes for the basal heat flux
if strcmp(flag,'flux_upwind')==1 || strcmp(flag,'flux_av')==1 %flux at hor edge is either averaged or upwinded between adjacent cells
    
    %SUBTEMPERATE SUBDOMAIN
    %pack parameters
    %nodes and edges needed for the computation of sliding velocity. averaging of flux between adjacent cells requires the first temperate node as well
    h_subt_edge_index = find(logical((h_id_node(h_up_node) == 2))==1);
    h_subt_edge_index = [h_subt_edge_index(1)-1; h_subt_edge_index];
    h_subt_node_index = [h_up_node(h_subt_edge_index(1)); find(h_id_node == 2); h_down_node(h_subt_edge_index(end))];
    T_subt_edge_index = find(logical((T_id_node(T_up_node_hor) == 2))==1);
    T_subt_edge_index = [T_subt_edge_index(1:nodes_ver)-nodes_ver; T_subt_edge_index];
    T_subt_node_index = [T_up_node_hor(T_subt_edge_index(1:nodes_ver));find(T_id_node == 2); T_down_node_hor(T_subt_edge_index(end-nodes_ver+1:end))];
    parameters.slide.h_node_index = h_subt_node_index;
    parameters.slide.h_edge_index = h_subt_edge_index;
    parameters.slide.T_edge_index = T_subt_edge_index;
    parameters.slide.T_node_index = T_subt_node_index;
    
    if nargout == 1
        subt_subd = subtemp_slid_flux_v4(v_in, parameters);
        U_SL(h_subt_edge_index) = subt_subd.U_SL;
        
    elseif nargout == 2
        %compute jacobian
        flag = 1; %edges inside the subtemperate subdomain
        [subt_subd, Dsubt_subd] = subtemp_slid_flux_v4(v_in, parameters, flag);
        U_SL(h_subt_edge_index) = subt_subd.U_SL;
        dUSL_dxc(h_subt_edge_index) = Dsubt_subd.dUSL_dxc;
        dUSL_dxs(h_subt_edge_index) = Dsubt_subd.dUSL_dxs;
        dUSL_dxt(h_subt_edge_index) = Dsubt_subd.dUSL_dxt;
        dUSL_dH_sub = Dsubt_subd.dUSL_dH_sub;
        dUSL_dT_sub = Dsubt_subd.dUSL_dT_sub;
        dQMa_dT_sub = Dsubt_subd.dQMa_dT_sub;
        dQMa_dH_sub = Dsubt_subd.dQMa_dH_sub;
        dQHa_dT_sub_basalvelocityjac = Dsubt_subd.dQHa_dT_sub_basalvelocityjac;
        dQHa_dH_sub_basalvelocityjac = Dsubt_subd.dQHa_dH_sub_basalvelocityjac;
    end
    
    %SUBTEMPERATE-TEMPERATE BOUNDARY
    %pack parameters
    h_subt_edge_index = find(logical((h_id_node(h_up_node) == 2).*(h_id_node(h_down_node) == 3))==1);
    h_subt_node_index = [h_up_node(h_subt_edge_index); h_down_node(h_subt_edge_index)];
    T_subt_edge_index = find(logical((T_id_node(T_up_node_hor) == 2).*(T_id_node(T_down_node_hor) == 3))==1);
    T_subt_node_index = [T_up_node_hor(T_subt_edge_index); T_down_node_hor(T_subt_edge_index)];
    parameters.slide.h_node_index = h_subt_node_index;
    parameters.slide.h_edge_index = h_subt_edge_index;
    parameters.slide.T_edge_index = T_subt_edge_index;
    parameters.slide.T_node_index = T_subt_node_index;
    
    if nargout == 1
        st_bdy = subtemp_slid_flux_v4(v_in, parameters);
        U_SL_s_sub = st_bdy.U_SL;
        
    elseif nargout == 2
        flag = 0; %compute U_SL and its jacobian at the sub-temp/temp boundary, for the relevant b.c.
        [st_bdy, Dst_bdy] = subtemp_slid_flux_v4(v_in, parameters,flag);
        U_SL_s_sub = st_bdy.U_SL;
        dUSL_dxc_st = Dst_bdy.dUSL_dxc;
        dUSL_dxs_st = Dst_bdy.dUSL_dxs;
        dUSL_dxt_st = Dst_bdy.dUSL_dxt;
        dUSLssub_dT = Dst_bdy.dUSL_dT_sub;
        dUSLssub_dH = Dst_bdy.dUSL_dH_sub;
        
        %add derivatives with respect to lengths of subdomains
    end
    
elseif strcmp(flag,'T_av')==1 %T is averaged

end

%TEMPERATE SUB-DOMAIN
index_temp_edge = find(h_id_node(h_up_node) == 3);
index_temp_edge_T = find(T_id_node(T_up_node_hor) == 3);

index_stb_edge = find(logical((h_id_node(h_up_node) == 2).*(h_id_node(h_down_node) == 3))==1);

H_edge_temp = H_edge(index_temp_edge);
dS_edge_temp = dS_edge(index_temp_edge);
T_hor_temp = T_hor(index_temp_edge_T);

H_edge_st = H_edge(index_stb_edge);
dS_edge_st = dS_edge(index_stb_edge);

U_SL(index_temp_edge) = -1/gamma*H_edge_temp.*dS_edge_temp;
%h_subt_edge_index
U_SL_s_temp = -1/gamma*H_edge_st.*dS_edge_st;

%jacobian
if nargout == 2
    %mass conservation
    indexcomposite_up_temp = h_up_node(index_temp_edge);
    indexcomposite_down_temp = h_down_node(index_temp_edge);
    
    indexcomposite_up_st = h_up_node(h_subt_edge_index);
    indexcomposite_down_st = h_down_node(h_subt_edge_index);
    
    dU_SL_dHedge_temp = -1/gamma.*dS_edge_temp;
    dU_SL_ddSedge_temp = -1/gamma.*H_edge_temp;
    dU_SL_dHedge_st = -1/gamma.*dS_edge_st;
    dU_SL_ddSedge_st = -1/gamma.*H_edge_st;
    
    dHedgedHup_temp = dHedgedHup(index_temp_edge);
    dHedgedHdown_temp = dHedgedHdown(index_temp_edge);
    ddS_edgedHup_temp = ddS_edgedHup(index_temp_edge);
    ddS_edgedHdown_temp = ddS_edgedHdown(index_temp_edge);
    dHedgeadv_dHup_temp = dHedgeadv_dHup(index_temp_edge);
    dHedgeadv_dHdown_temp = dHedgeadv_dHdown(index_temp_edge);
    
    dHedgedHup_st = dHedgedHup(index_stb_edge);
    dHedgedHdown_st = dHedgedHdown(index_stb_edge);
    ddS_edgedHup_st = ddS_edgedHup(index_stb_edge);
    ddS_edgedHdown_st = ddS_edgedHdown(index_stb_edge);
    
    dU_SL_dHup_temp = dU_SL_dHedge_temp.*dHedgedHup_temp + dU_SL_ddSedge_temp.*ddS_edgedHup_temp;
    dU_SL_dHdown_temp = dU_SL_dHedge_temp.*dHedgedHdown_temp + dU_SL_ddSedge_temp.*ddS_edgedHdown_temp;
    
    dU_SL_dHup_st = dU_SL_dHedge_st.*dHedgedHup_st + dU_SL_ddSedge_st.*ddS_edgedHup_st;
    dU_SL_dHdown_st = dU_SL_dHedge_st.*dHedgedHdown_st + dU_SL_ddSedge_st.*ddS_edgedHdown_st;
    
    dQMa_dHu_temp = H_edge_adv(index_temp_edge).*dU_SL_dHup_temp + U_SL(index_temp_edge).*dHedgeadv_dHup_temp;
    dQMa_dHd_temp = H_edge_adv(index_temp_edge).*dU_SL_dHdown_temp + U_SL(index_temp_edge).*dHedgeadv_dHdown_temp;
    
    dQMa_dH_temp = sparse(indexcomposite_up_temp,indexcomposite_up_temp,dQMa_dHu_temp,h_nodes, h_nodes)+...
        sparse(indexcomposite_up_temp,indexcomposite_down_temp,dQMa_dHd_temp ,h_nodes, h_nodes)-...
        (sparse(indexcomposite_down_temp,indexcomposite_up_temp,dQMa_dHu_temp,h_nodes, h_nodes)+...
        sparse(indexcomposite_down_temp,indexcomposite_down_temp,dQMa_dHd_temp ,h_nodes, h_nodes));
    dUSL_dH_temp = sparse(indexcomposite_up_temp,indexcomposite_up_temp,dU_SL_dHup_temp,h_nodes, h_nodes)+...
        sparse(indexcomposite_up_temp,indexcomposite_down_temp,dU_SL_dHdown_temp ,h_nodes, h_nodes);
    
    %heat conservation
    T_up_node_hor_temp = T_up_node_hor(index_temp_edge_T);
    T_down_node_hor_temp = T_down_node_hor(index_temp_edge_T);

    index_h_to_T_temp = reshape(repmat((1:length(index_temp_edge)),[nodes_ver,1]),[length(index_temp_edge_T),1]);
    T_indexcomposite_up_temp = indexcomposite_up_temp(index_h_to_T_temp);
    T_indexcomposite_down_temp = indexcomposite_down_temp(index_h_to_T_temp);
    
    H_edge_T = H_edge_temp(index_h_to_T_temp);
    dUSLdHu_T = dU_SL_dHup_temp(index_h_to_T_temp);
    dUSLdHd_T = dU_SL_dHdown_temp(index_h_to_T_temp);
    dQHa_dH_temp_basalvelocity_u = H_edge_T.*T_hor_temp.* dUSLdHu_T;
    dQHa_dH_temp_basalvelocity_d = H_edge_T.*T_hor_temp.* dUSLdHd_T;
    
    dQHa_dH_temp_basalvelocityjac = sparse(T_up_node_hor_temp,T_indexcomposite_up_temp,dQHa_dH_temp_basalvelocity_u,T_nodes, h_nodes)+...
                sparse(T_up_node_hor_temp,T_indexcomposite_down_temp,dQHa_dH_temp_basalvelocity_d ,T_nodes, h_nodes)-...
                (sparse(T_down_node_hor_temp,T_indexcomposite_up_temp,dQHa_dH_temp_basalvelocity_u,T_nodes, h_nodes)+...
                sparse(T_down_node_hor_temp,T_indexcomposite_down_temp,dQHa_dH_temp_basalvelocity_d ,T_nodes, h_nodes));
            
    %length of subdomains
    dUSL_dxc(index_temp_edge) = dU_SL_ddSedge_temp.*Ddiscvar.dS.ddSedge_dxc(index_temp_edge) + dU_SL_dHedge_temp.*dHedgedxc(index_temp_edge);
    dUSL_dxs(index_temp_edge) = dU_SL_ddSedge_temp.*Ddiscvar.dS.ddSedge_dxs(index_temp_edge) + dU_SL_dHedge_temp.*dHedgedxs(index_temp_edge);
    dUSL_dxt(index_temp_edge) = dU_SL_ddSedge_temp.*Ddiscvar.dS.ddSedge_dxt(index_temp_edge) + dU_SL_dHedge_temp.*dHedgedxt(index_temp_edge);
    %sliding velocity at sub-temp/temp boundary  
    dUSLstemp_dH = sparse(1, indexcomposite_up_st,dU_SL_dHup_st,1,h_nodes)+...
        sparse(1,indexcomposite_down_st,dU_SL_dHdown_st,1,h_nodes) ;
    dUSLtemp_dxc_st = dU_SL_ddSedge_st.*Ddiscvar.dS.ddSedge_dxc(index_stb_edge) + dU_SL_dHedge_st.*dHedgedxc(index_stb_edge);
    dUSLtemp_dxs_st = dU_SL_ddSedge_st.*Ddiscvar.dS.ddSedge_dxs(index_stb_edge) + dU_SL_dHedge_st.*dHedgedxs(index_stb_edge);
    dUSLtemp_dxt_st = dU_SL_ddSedge_st.*Ddiscvar.dS.ddSedge_dxt(index_stb_edge) + dU_SL_dHedge_st.*dHedgedxt(index_stb_edge);
end

%jacobian of the non u_sl-dependent parts of the horizontal advective heat
%flux by sliding, for all subdomains
index_h_to_T = reshape(repmat((1:h_n_edges),[nodes_ver,1]),[T_n_edges_hor,1]);
if nargout == 2
    %temperature
    dQHa_sl2_dT_u = U_SL(index_h_to_T).*H_edge(index_h_to_T).*dThor_dTup;
    dQHa_sl2_dT_d = U_SL(index_h_to_T).*H_edge(index_h_to_T).*dThor_dTdown;
    
    dQHa_sl2_dT = sparse(T_up_node_hor,T_up_node_hor,dQHa_sl2_dT_u,T_nodes, T_nodes)+...
        sparse(T_up_node_hor,T_down_node_hor,dQHa_sl2_dT_d ,T_nodes, T_nodes)-...
        (sparse(T_down_node_hor,T_up_node_hor,dQHa_sl2_dT_u,T_nodes, T_nodes)+...
        sparse(T_down_node_hor, T_down_node_hor,dQHa_sl2_dT_d ,T_nodes, T_nodes));
    %ice thickness
    indexcomposite_u = h_up_node(index_h_to_T);
    indexcomposite_d = h_down_node(index_h_to_T);
    
    dQHa_sl2_dH_u = U_SL(index_h_to_T).*T_hor.*dHedgedHup(index_h_to_T);
    dQHa_sl2_dH_d = U_SL(index_h_to_T).*T_hor.*dHedgedHdown(index_h_to_T);
    
    dQHa_sl2_dh = sparse(T_up_node_hor,indexcomposite_u,dQHa_sl2_dH_u,T_nodes, h_nodes)+...
        sparse(T_up_node_hor,indexcomposite_d,dQHa_sl2_dH_d ,T_nodes, h_nodes)-...
        (sparse(T_down_node_hor,indexcomposite_u,dQHa_sl2_dH_u,T_nodes, h_nodes)+...
        sparse(T_down_node_hor, indexcomposite_d,dQHa_sl2_dH_d ,T_nodes, h_nodes));
    %length of subdomains
    dQMa_dxc = (accumarray(h_up_node, dUSL_dxc.*H_edge_adv, [h_nodes,1])-accumarray(h_down_node, dUSL_dxc.*H_edge_adv, [h_nodes,1]));
    dQMa_dxs = (accumarray(h_up_node, dUSL_dxs.*H_edge_adv, [h_nodes,1])-accumarray(h_down_node, dUSL_dxs.*H_edge_adv, [h_nodes,1]));
    dQMa_dxt = (accumarray(h_up_node, dUSL_dxt.*H_edge_adv, [h_nodes,1])-accumarray(h_down_node, dUSL_dxt.*H_edge_adv, [h_nodes,1]));
    
    dQHa_dxc_edge = T_hor.*(dUSL_dxc(index_h_to_T).*H_edge(index_h_to_T) + U_SL(index_h_to_T).*dHedgedxc(index_h_to_T));
    dQHa_dxs_edge = T_hor.*(dUSL_dxs(index_h_to_T).*H_edge(index_h_to_T) + U_SL(index_h_to_T).*dHedgedxs(index_h_to_T));
    dQHa_dxt_edge = T_hor.*(dUSL_dxt(index_h_to_T).*H_edge(index_h_to_T) + U_SL(index_h_to_T).*dHedgedxt(index_h_to_T));
    
    dQHa_dxc = (accumarray(T_up_node_hor, dQHa_dxc_edge, [T_nodes,1])-accumarray(T_down_node_hor, dQHa_dxc_edge, [T_nodes,1]));
    dQHa_dxs = (accumarray(T_up_node_hor, dQHa_dxs_edge, [T_nodes,1])-accumarray(T_down_node_hor, dQHa_dxs_edge, [T_nodes,1]));
    dQHa_dxt = (accumarray(T_up_node_hor, dQHa_dxt_edge, [T_nodes,1])-accumarray(T_down_node_hor, dQHa_dxt_edge, [T_nodes,1]));
    
end
%assemble output
fout.U_SL = U_SL;
fout.QM_a = U_SL.*H_edge_adv;
fout.QHh_a = U_SL(index_h_to_T).*H_edge(index_h_to_T).*T_hor;
fout.U_SL_sub_st = U_SL_s_sub;
fout.U_SL_temp_st = U_SL_s_temp;

%assemble output for calculation of vertical velocity
fout.U_sl_Thoredges = U_SL(index_h_to_T);


if nargout == 2
    Dfout.dQMa_dT = dQMa_dT_sub;
    Dfout.dQMa_dH = dQMa_dH_sub + dQMa_dH_temp;
    Dfout.dQMa_dxc = dQMa_dxc;
    Dfout.dQMa_dxs = dQMa_dxs;
    Dfout.dQMa_dxt = dQMa_dxt;
    
    Dfout.dQHa_dT = dQHa_dT_sub_basalvelocityjac + dQHa_sl2_dT;
    Dfout.dQHa_dH = dQHa_dH_temp_basalvelocityjac + dQHa_dH_sub_basalvelocityjac+ dQHa_sl2_dh;
    Dfout.dQHa_dxc = dQHa_dxc;
    Dfout.dQHa_dxs = dQHa_dxs;
    Dfout.dQHa_dxt = dQHa_dxt;
    
    Dfout.dUSLssub_dT = dUSLssub_dT;
    Dfout.dUSLssub_dH = dUSLssub_dH;
    Dfout.dUSLssub_dxc = dUSL_dxc_st;
    Dfout.dUSLssub_dxs = dUSL_dxs_st;
    Dfout.dUSLssub_dxt = dUSL_dxt_st;
    
    Dfout.dUSLstemp_dH = dUSLstemp_dH;
    Dfout.dUSLstemp_dxc = dUSLtemp_dxc_st;
    Dfout.dUSLstemp_dxs = dUSLtemp_dxs_st;
    Dfout.dUSLstemp_dxt = dUSLtemp_dxt_st; 
    
    dUSL_dH = dUSL_dH_sub + dUSL_dH_temp;
    Dfout.dUSL = [dUSL_dH dUSL_dT_sub [dUSL_dxc dUSL_dxs dUSL_dxt; sparse(1,3)]];
end
end






















