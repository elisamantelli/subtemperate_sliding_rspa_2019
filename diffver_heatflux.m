function [fout, Dfout] = diffver_heatflux(v_in, parameters)
% function computes vertical diffusive heat flux as well as its jacobian
% Not yet tested or benchmarked. Elisa Mantelli, Feb 9th, 2017

%unpack parameters
h_nodes = parameters.grid_h.n_nodes;                                       %number of nodes,h
T_nodes = parameters.grid_T.n_nodes.tot;                                   %number of temperature nodes
nodes_ver = parameters.grid_T.n_nodes.vert;                                %number of nodes in the vertical direction
T_n_edges_ver = parameters.grid_T.n_edges.vert;                            %number of vertical edges,T
T_up_node_ver = parameters.grid_T.up_node.vert;                            %list (n_edges-by-1 vector) of 'upstream' node for each vertical edge,T
T_down_node_ver = parameters.grid_T.down_node.vert;                        %list (n_edges-by-1 vector) of 'downstream' nodes for each vertical edge,T
T_Delta_eta = parameters.grid_T.Delta_eta;                                 %scalar, vertical spacing between cell centres

%unpack input variable v_in
H = v_in(1:h_nodes);
T = v_in(h_nodes+1:h_nodes+T_nodes);
xc = v_in(end-2);
xs = v_in(end-1);
xt = v_in(end);

[length_T,Dlength_T] = fvlength([xc;xs;xt], parameters,'T');
L_T = length_T.L;
dLT_dxc = Dlength_T.dL_dxc;
dLT_dxs = Dlength_T.dL_dxs;
dLT_dxt = Dlength_T.dL_dxt;

% reindexing
index_h_to_T = reshape(repmat(1:h_nodes,[nodes_ver,1]),[T_nodes,1]);
index_h_to_T_veredges = reshape(repmat(1:h_nodes,[nodes_ver-1,1]),[T_n_edges_ver,1]);
index_Tnodes_to_Tveredges = repmat([ones(nodes_ver-1,1);0],[h_nodes,1]);

T_down = T(T_down_node_ver);
T_up = T(T_up_node_ver);
H_T = H(index_h_to_T);
H_Tup = H_T(T_up_node_ver);
H_Tdown = H_T(T_down_node_ver);
L_T_veredges = L_T(logical(index_Tnodes_to_Tveredges));

%diffusive heat flux
fout.QH_v_d = - L_T_veredges.*(T_down-T_up)./(T_Delta_eta*H(index_h_to_T_veredges));
%jacobian
if nargout == 2
    h_node_list = (1:h_nodes)';
    h_node_list_Tveredges = h_node_list(index_h_to_T_veredges);
    
    dQHvd_dTup = - L_T_veredges.*(-1)./(T_Delta_eta*H_Tup);
    dQHvd_dTdown = - L_T_veredges.*(1)./(T_Delta_eta*H_Tdown);
    dQHvd_dH = - fout.QH_v_d./H(index_h_to_T_veredges);
    dQHvd_dxc = - dLT_dxc(logical(index_Tnodes_to_Tveredges)).*(T_down-T_up)./(T_Delta_eta*H(index_h_to_T_veredges));
    dQHvd_dxs = - dLT_dxs(logical(index_Tnodes_to_Tveredges)).*(T_down-T_up)./(T_Delta_eta*H(index_h_to_T_veredges));
    dQHvd_dxt = - dLT_dxt(logical(index_Tnodes_to_Tveredges)).*(T_down-T_up)./(T_Delta_eta*H(index_h_to_T_veredges));
    
    Dfout.dQHvd_dT = sparse(T_up_node_ver,T_up_node_ver,dQHvd_dTup,T_nodes, T_nodes)+...
        sparse(T_up_node_ver,T_down_node_ver,dQHvd_dTdown ,T_nodes, T_nodes)-...
        (sparse(T_down_node_ver,T_up_node_ver,dQHvd_dTup,T_nodes, T_nodes)+...
        sparse(T_down_node_ver,T_down_node_ver,dQHvd_dTdown,T_nodes, T_nodes));
    Dfout.dQHvd_dH = sparse(T_up_node_ver,h_node_list_Tveredges,dQHvd_dH,T_nodes, h_nodes)-...
    sparse(T_down_node_ver,h_node_list_Tveredges,dQHvd_dH,T_nodes, h_nodes);
    Dfout.dQHvd_dxc = sparse(T_up_node_ver,ones(length(T_up_node_ver),1),dQHvd_dxc,T_nodes,1)-...
        sparse(T_down_node_ver,ones(length(T_up_node_ver),1),dQHvd_dxc,T_nodes,1);
    Dfout.dQHvd_dxs = sparse(T_up_node_ver,ones(length(T_up_node_ver),1),dQHvd_dxs,T_nodes,1)-...
        sparse(T_down_node_ver,ones(length(T_up_node_ver),1),dQHvd_dxs,T_nodes,1);
    Dfout.dQHvd_dxt = sparse(T_up_node_ver,ones(length(T_up_node_ver),1),dQHvd_dxt,T_nodes,1)-...
        sparse(T_down_node_ver,ones(length(T_up_node_ver),1),dQHvd_dxt,T_nodes,1);
end
end