function fout = network_sstate_jacobian_v5(v_in,parameters)
%returns jacobian of network_sstate_v4
%Tested against numerical jacobian, Jun 8th 2017. 

%unpack parameters
h_nodes = parameters.grid_h.n_nodes;                  %number of nodes,h
h_up_node = parameters.grid_h.up_node;                %list (h_n_edge-by-1 vector) of 'upstream' node, h
h_down_node = parameters.grid_h.down_node;            %list (h_n_edges-by-1 vector) of 'downstream' node,h

T_nodes = parameters.grid_T.n_nodes.tot;                                   %number of temperature nodes
nodes_ver = parameters.grid_T.n_nodes.vert;                                %number of nodes in the vertical direction
T_n_edges_ver = parameters.grid_T.n_edges.vert;                            %number of vertical edges,T
T_up_node_ver = parameters.grid_T.up_node.vert;                            %list (n_edges-by-1 vector) of 'upstream' node for each vertical edge,T
T_down_node_ver = parameters.grid_T.down_node.vert;                        %list (n_edges-by-1 vector) of 'downstream' nodes for each vertical edge,T
T_Delta_eta = parameters.grid_T.Delta_eta;                                 %scalar, vertical spacing between cell centres
T_bdy_nodes_flux = parameters.grid_T.bdy_nodes.flux;                       %list (n_bdy_nodes-by-1 vector) of nodes at domain boundaries where Neumann conditions
T_bdy_nodes_dir_surf = parameters.grid_T.bdy_nodes.dir.surf;               %list  of nodes at the upper surface where Dirichlet conditions apply
T_bdy_nodes_dir_bed = parameters.grid_T.bdy_nodes.dir.bed;                 %list  of nodes at the bed where Dirichlet conditions apply
T_bdy_nodes_melt = parameters.grid_T.bdy_nodes.T_melt;                     %list of edges needed for T=0 at the bed at the cold/subtemp boundary
T_bdy_nodes_subtemp = parameters.grid_T.bdy_nodes.subtemp_slid;            %list of basal T nodes in the subtemperate subdomain

alpha = parameters.alpha;   % strain heating parameter
gamma = parameters.gamma;   %friction coefficient
Pe = parameters.Pe;         % Peclet number
nu = parameters.nu;         % basal heating parameter
a = parameters.a;           %accumulation rate
r_rho = parameters.r_rho;   %water to ice density ratio

T_surf = parameters.T_surf; %surface temperature

%unpack input variable v_in
H = v_in(1:h_nodes);
T = v_in(h_nodes+1:h_nodes+T_nodes);
xc = v_in(h_nodes+T_nodes+1);
xs = v_in(h_nodes+T_nodes+2);
xt = v_in(h_nodes+T_nodes+3);

%initialize fout
fout = sparse(h_nodes+T_nodes+3,h_nodes+T_nodes+3);
%%  PRELIMINARIES 

%horizontal length of subdomains, cells and edges 
[length_h,Dlength_h] = fvlength([xc;xs;xt], parameters,'h');
[length_T,Dlength_T] = fvlength([xc;xs;xt], parameters,'T');
%unpack
L_h = length_h.L;
dLh_dxc = Dlength_h.dL_dxc;
dLh_dxs = Dlength_h.dL_dxs;
dLh_dxt = Dlength_h.dL_dxt;
Delta_sigma_h = length_h. Delta_sigma;
L_T = length_T.L;
Delta_sigma_T = length_T. Delta_sigma;
dLT_dxc = Dlength_T.dL_dxc;
dLT_dxs = Dlength_T.dL_dxs;
dLT_dxt = Dlength_T.dL_dxt;
%BED
[st_bed,Dst_bed] = bed([xc;xs;xt], parameters);
%unpack: depth to bed at grounding line
dBg_dxc = Dst_bed.dBg_dxc;
dBg_dxs = Dst_bed.dBg_dxs;
dBg_dxt = Dst_bed.dBg_dxt;
%DISCRETE DEP. VARIABLES
[discvar, Ddiscvar] = discretisation_v2(v_in, parameters);
%upwinding for bed flux
upwind = upwinding(parameters);
theta_QHv = upwind.theta_QHv;

%FLUX 
dQHvbed_dTbedup = Ddiscvar.flux.dQHv_dTbedup;
dQHvbed_dTbeddown = Ddiscvar.flux.dQHv_dTbeddown;
dQHvbed_dHup = Ddiscvar.flux.dQHv_dHup;
dQHvbed_dHdown = Ddiscvar.flux.dQHv_dHdown;
dQHvbed_dxc = Ddiscvar.flux.dQHv_dxc;
dQHvbed_dxs = Ddiscvar.flux.dQHv_dxs;
dQHvbed_dxt = Ddiscvar.flux.dQHv_dxt;

%ICE THICKNESS at the grounding line
H_g = discvar.H.H_g;
dHg_dH = Ddiscvar.H.dHg_dH;

%TEMPERATURE 
%along T edges 
T_vert = discvar.T.T_vert;
dTver_dTup = Ddiscvar.T.dTver_dTup;
dTver_dTdown = Ddiscvar.T.dTver_dTdown;

%along hor edges at the grounding line
T_hor_g = discvar.T.T_hor_g;
dThorg_dT = Ddiscvar.T.dThorg_dT_up;
%SHEARING VELOCITY FIELD
[hor_velocity, Dhor_velocity] = hor_velocity_SIA_v2(v_in, parameters);
QMd = hor_velocity.QMd;
dQMd_dH = Dhor_velocity.dQMd_dH_net;
dQMd_dxc = Dhor_velocity.dQMd_dxc_net;
dQMd_dxs = Dhor_velocity.dQMd_dxs_net;
dQMd_dxt = Dhor_velocity.dQMd_dxt_net;

S = hor_velocity.S;
dS_dH = Dhor_velocity.dS_dH;
dS_dxc = Dhor_velocity.dS_dxc;
dS_dxs = Dhor_velocity.dS_dxs;
dS_dxt = Dhor_velocity.dS_dxt;

dQHa_shear_dH = Dhor_velocity.dQHa_shear_dH_net;
dQHa_shear_dT = Dhor_velocity.dQHa_shear_dT_net;
dQHa_shear_dxc = Dhor_velocity.dQHa_shear_dxc_net;
dQHa_shear_dxs = Dhor_velocity.dQHa_shear_dxs_net;
dQHa_shear_dxt = Dhor_velocity.dQHa_shear_dxt_net;

net_flux_shear_eta_sum = hor_velocity.net_flux_shear_eta_sum;
Dnet_flux_shear_eta_sum_dH = Dhor_velocity.dnfse_sum_dH;
Dnet_flux_shear_eta_sum_dxc = Dhor_velocity.dnfse_sum_dxc;
Dnet_flux_shear_eta_sum_dxs = Dhor_velocity.dnfse_sum_dxs;
Dnet_flux_shear_eta_sum_dxt = Dhor_velocity.dnfse_sum_dxt;
%grounding line
QMd_g = hor_velocity.QMd_g;
dQMdg_dH = Dhor_velocity.dQMdg_dH;
dQMdg_dxc = Dhor_velocity.dQMdg_dxc;
dQMdg_dxs = Dhor_velocity.dQMdg_dxs;
dQMdg_dxt = Dhor_velocity.dQMdg_dxt;

%SLIDING VELOCITY
% compute sliding velocity at the grounding line
[flux,Dflux] = fluxg(H_g,parameters);
QMg = flux.flux;
dQMg_dHg = Dflux.dflux_dHg;
U_SL_g = (QMg-QMd_g)/H(end);
dUSLg_dH = -(QMg-QMd_g)/H(end)^2*[sparse(h_nodes-1,1); 1] + (dQMg_dHg/H(end))*dHg_dH - dQMdg_dH./H(end);
dUSLg_dxc = (-dQMdg_dxc)/H(end);
dUSLg_dxs = (-dQMdg_dxs)/H(end);
dUSLg_dxt = (-dQMdg_dxt)/H(end);

%along h edges
[usl, Dusl] = sliding_v3(v_in, parameters);
dQMa_dH = Dusl.dQMa_dH;
dQMa_dT = Dusl.dQMa_dT;
dQMa_dxc = Dusl.dQMa_dxc;
dQMa_dxs = Dusl.dQMa_dxs;
dQMa_dxt = Dusl.dQMa_dxt;

dQHaslide_dH = Dusl.dQHa_dH;
dQHaslide_dT = Dusl.dQHa_dT;
dQHaslide_dxc = Dusl.dQHa_dxc;
dQHaslide_dxs = Dusl.dQHa_dxs;
dQHaslide_dxt = Dusl.dQHa_dxt;

%temperate sliding velocity at x = xs
U_SL_temp_st = usl.U_SL_temp_st ; 
dUSLstemp_dH = Dusl.dUSLstemp_dH ;
dUSLstemp_dxc = Dusl.dUSLstemp_dxc ;
dUSLstemp_dxs = Dusl.dUSLstemp_dxs ;
dUSLstemp_dxt = Dusl.dUSLstemp_dxt ;
%% MASS CONSERVATION
%jacobian of net mass flux
dQM_dH = dQMd_dH + dQMa_dH;
dQM_dT = dQMa_dT;
dQM_dxc = dQMd_dxc + dQMa_dxc;
dQM_dxs = dQMd_dxs+ dQMa_dxs;
dQM_dxt = dQMd_dxt + dQMa_dxt;
% correct for grounding line outflow, only sliding component U_SL_g*H(end) needs correction
dQM_dH(end,:) = dQM_dH(end,:) + dUSLg_dH.'*H(end) + U_SL_g*[sparse(h_nodes-1,1); 1].';
dQM_dxc(end,:) = dQM_dxc(end,:) + dUSLg_dxc.'*H(end);
dQM_dxs(end,:) = dQM_dxs(end,:) + dUSLg_dxs.'*H(end);
dQM_dxt(end,:) = dQM_dxt(end,:) + dUSLg_dxt.'*H(end);
%jacobian of mass conservation fout(1:h_nodes) = netmassflux - a*L_h.*Delta_sigma_h;
fout(1:h_nodes,:) = [ spdiags(1./(Delta_sigma_h),0,h_nodes,h_nodes)*dQM_dH,  spdiags(1./(Delta_sigma_h),0,h_nodes,h_nodes)*dQM_dT, dQM_dxc./Delta_sigma_h-a*(dLh_dxc), dQM_dxs./Delta_sigma_h-a*(dLh_dxs), dQM_dxt./Delta_sigma_h-a*(dLh_dxt) ];
%% HEAT CONSERVATION
%HORIZONTAL HEAT FLUX (only advection) 
%jacobian
dQHh_dH = dQHaslide_dH+dQHa_shear_dH;
dQHh_dT = dQHaslide_dT+dQHa_shear_dT;
dQHh_dxc = dQHaslide_dxc+dQHa_shear_dxc;
dQHh_dxs = dQHaslide_dxs+dQHa_shear_dxs;
dQHh_dxt = dQHaslide_dxt+dQHa_shear_dxt;
%correct for outflow at the grounding line (flux by shearing already
%corrected)
dQHh_dH(end-nodes_ver+1:end, :) = dQHh_dH(end-nodes_ver+1:end, :)+ spdiags(T_hor_g,0,nodes_ver,nodes_ver)*repmat((dUSLg_dH*H_g).',[nodes_ver,1])+...
   spdiags(T_hor_g,0,nodes_ver,nodes_ver)*repmat((U_SL_g*dHg_dH).',[nodes_ver,1]);
dQHh_dT(end-nodes_ver+1:end, :) = dQHh_dT(end-nodes_ver+1:end, :)+ dThorg_dT*U_SL_g*H_g;
dQHh_dxc(end-nodes_ver+1:end) = dQHh_dxc(end-nodes_ver+1:end)+ T_hor_g.*(H_g*dUSLg_dxc);
dQHh_dxs(end-nodes_ver+1:end) = dQHh_dxs(end-nodes_ver+1:end)+ T_hor_g.*(H_g*dUSLg_dxs);
dQHh_dxt(end-nodes_ver+1:end) = dQHh_dxt(end-nodes_ver+1:end)+ T_hor_g.*(H_g*dUSLg_dxt);
%VERTICAL HEAT FLUX along vertical edges + flux bc at the bed
%diffusive component
[QHvd,dQHvd] = diffver_heatflux(v_in, parameters);
dQHvd_dT = dQHvd.dQHvd_dT;
dQHvd_dH = dQHvd.dQHvd_dH;
dQHvd_dxc = dQHvd.dQHvd_dxc;
dQHvd_dxs = dQHvd.dQHvd_dxs;
dQHvd_dxt = dQHvd.dQHvd_dxt;

%correct ver flux in order to account for basal temperature at the melting
%point in subtemperate and temperate subdomain
%net_H_verflux(T_bdy_nodes_dir_bed) = net_H_verflux(T_bdy_nodes_dir_bed)- ( - L_T(T_bdy_nodes_dir_bed).*T(T_bdy_nodes_dir_bed)./(H((T_bdy_nodes_dir_bed-1)/nodes_ver+1).*T_Delta_eta/2));

index_T_to_h = (T_bdy_nodes_dir_bed-1)/nodes_ver+1;
flux_bed =  - L_T(T_bdy_nodes_dir_bed).*T(T_bdy_nodes_dir_bed)./(H(index_T_to_h).*T_Delta_eta/2);
dfluxbed_dT =  - L_T(T_bdy_nodes_dir_bed)./(H(index_T_to_h).*T_Delta_eta/2);
dfluxbed_dH = - flux_bed./H(index_T_to_h);
dfluxbed_dxc = - T(T_bdy_nodes_dir_bed)./(H(index_T_to_h).*T_Delta_eta/2).*dLT_dxc(T_bdy_nodes_dir_bed);
dfluxbed_dxs = - T(T_bdy_nodes_dir_bed)./(H(index_T_to_h).*T_Delta_eta/2).*dLT_dxs(T_bdy_nodes_dir_bed);
dfluxbed_dxt = - T(T_bdy_nodes_dir_bed)./(H(index_T_to_h).*T_Delta_eta/2).*dLT_dxt(T_bdy_nodes_dir_bed);

dQHvd_dT(T_bdy_nodes_dir_bed,:) = dQHvd_dT(T_bdy_nodes_dir_bed,:) - sparse(1:length(T_bdy_nodes_dir_bed), T_bdy_nodes_dir_bed, dfluxbed_dT,length(T_bdy_nodes_dir_bed),T_nodes);
dQHvd_dH(T_bdy_nodes_dir_bed,:) = dQHvd_dH(T_bdy_nodes_dir_bed,:) - sparse(1:length(T_bdy_nodes_dir_bed), index_T_to_h, dfluxbed_dH ,length(T_bdy_nodes_dir_bed),h_nodes);
dQHvd_dxc(T_bdy_nodes_dir_bed) = dQHvd_dxc(T_bdy_nodes_dir_bed) - dfluxbed_dxc;
dQHvd_dxs(T_bdy_nodes_dir_bed) = dQHvd_dxs(T_bdy_nodes_dir_bed) - dfluxbed_dxs;
dQHvd_dxt(T_bdy_nodes_dir_bed) = dQHvd_dxt(T_bdy_nodes_dir_bed) - dfluxbed_dxt;

%%correct ver flux at the surface!!!
%net_H_verflux(T_bdy_nodes_dir_surf) = net_H_verflux(T_bdy_nodes_dir_surf) + ( - L_T(T_bdy_nodes_dir_surf).*(-1 - T(T_bdy_nodes_dir_surf))./(H((T_bdy_nodes_dir_surf)/nodes_ver).*T_Delta_eta/2))
index_T_to_h_surf = (T_bdy_nodes_dir_surf)/nodes_ver;
flux_surf =  - L_T(T_bdy_nodes_dir_surf).*(-1 - T(T_bdy_nodes_dir_surf))./(H(index_T_to_h_surf).*T_Delta_eta/2);
dfluxsurf_dT =   L_T(T_bdy_nodes_dir_surf)./(H(index_T_to_h_surf).*T_Delta_eta/2);
dfluxsurf_dH = - flux_surf./H(index_T_to_h_surf);
dfluxsurf_dxc =  -(-1 - T(T_bdy_nodes_dir_surf))./(H(index_T_to_h_surf).*T_Delta_eta/2).*dLT_dxc(T_bdy_nodes_dir_surf);
dfluxsurf_dxs =  -(-1 - T(T_bdy_nodes_dir_surf))./(H(index_T_to_h_surf).*T_Delta_eta/2).*dLT_dxs(T_bdy_nodes_dir_surf);
dfluxsurf_dxt =  -(-1 - T(T_bdy_nodes_dir_surf))./(H(index_T_to_h_surf).*T_Delta_eta/2).*dLT_dxt(T_bdy_nodes_dir_surf);

dQHvd_dT(T_bdy_nodes_dir_surf,:) = dQHvd_dT(T_bdy_nodes_dir_surf,:) + sparse(1:length(T_bdy_nodes_dir_surf), T_bdy_nodes_dir_surf, dfluxsurf_dT,length(T_bdy_nodes_dir_surf),T_nodes);
dQHvd_dH(T_bdy_nodes_dir_surf,:) = dQHvd_dH(T_bdy_nodes_dir_surf,:) + sparse(1:length(T_bdy_nodes_dir_surf), index_T_to_h_surf, dfluxsurf_dH ,length(T_bdy_nodes_dir_surf),h_nodes);
dQHvd_dxc(T_bdy_nodes_dir_surf) = dQHvd_dxc(T_bdy_nodes_dir_surf) + dfluxsurf_dxc;
dQHvd_dxs(T_bdy_nodes_dir_surf) = dQHvd_dxs(T_bdy_nodes_dir_surf) + dfluxsurf_dxs;
dQHvd_dxt(T_bdy_nodes_dir_surf) = dQHvd_dxt(T_bdy_nodes_dir_surf) + dfluxsurf_dxt;

%advective component
netmassflux_shear = (accumarray(h_up_node, QMd, [h_nodes,1])-accumarray(h_down_node, QMd, [h_nodes,1]));
netmassflux_shear(end) = netmassflux_shear(end)+ QMd_g;
%computation of vertical velocity and its jac
W_eff =T_Delta_eta.*reshape((1:nodes_ver)'*(netmassflux_shear./(L_h.*Delta_sigma_h) - a).',[T_nodes,1])-(net_flux_shear_eta_sum./L_T);
index_u_to_T = repmat([ones(nodes_ver-1,1); 0], [h_nodes,1]);%reindex from u to T ver network
L_T_veredges = L_T(index_u_to_T==1);

%jacobian
index_red_veredges = (1:T_nodes)';
index_red_veredges(~mod(index_red_veredges,nodes_ver))=[]; %exclude ice surface

dWeff_dH = T_Delta_eta*kron(speye(h_nodes),(1:nodes_ver)')*(spdiags(1./(L_h.*Delta_sigma_h),0,h_nodes,h_nodes)*dQMd_dH)...
    - spdiags(1./L_T,0,T_nodes,T_nodes)*Dnet_flux_shear_eta_sum_dH;
dWeff_dxc = T_Delta_eta.*reshape((1:nodes_ver)'*(dQMd_dxc./(L_h.*Delta_sigma_h) - netmassflux_shear./(L_h.^2.*Delta_sigma_h).*dLh_dxc).',[T_nodes,1])...
    - (Dnet_flux_shear_eta_sum_dxc./L_T - net_flux_shear_eta_sum./L_T.^2.*dLT_dxc);
dWeff_dxs = T_Delta_eta.*reshape((1:nodes_ver)'*(dQMd_dxs./(L_h.*Delta_sigma_h) - netmassflux_shear./(L_h.^2.*Delta_sigma_h).*dLh_dxs).',[T_nodes,1])...
    - (Dnet_flux_shear_eta_sum_dxs./L_T - net_flux_shear_eta_sum./L_T.^2.*dLT_dxs);
dWeff_dxt = T_Delta_eta.*reshape((1:nodes_ver)'*(dQMd_dxt./(L_h.*Delta_sigma_h) - netmassflux_shear./(L_h.^2.*Delta_sigma_h).*dLh_dxt).',[T_nodes,1])...
    - (Dnet_flux_shear_eta_sum_dxt ./L_T - net_flux_shear_eta_sum./L_T.^2.*dLT_dxt);
%ice thickness QH_v_a = L_T_veredges.*W_eff(index_u_to_T==1).*T_vert;
%separate surface contribution
dWeff_dH_surf = dWeff_dH(index_u_to_T==0,:);
dWeff_dH = dWeff_dH(index_u_to_T==1,:);

dQHva_dH_matrix = spdiags(L_T_veredges.*T_vert,0,T_n_edges_ver,T_n_edges_ver)*dWeff_dH;
[i,j,s_dQHva_dH] = find(dQHva_dH_matrix);
down_node_ver =  i+ceil(i./(nodes_ver-1));
up_node_ver = i+floor((i-1)./(nodes_ver-1));
dQHvadH_up = sparse(up_node_ver,j,s_dQHva_dH,T_nodes,h_nodes);
dQHvadH_down = sparse(down_node_ver,j,s_dQHva_dH,T_nodes,h_nodes);
dQHva_dH = dQHvadH_up - dQHvadH_down;
%correct for surface temperature
dQHva_dH(T_bdy_nodes_dir_surf,:) = dQHva_dH(T_bdy_nodes_dir_surf,:) + dWeff_dH_surf*T_surf.*spdiags(L_T(T_bdy_nodes_dir_surf),0,length(T_bdy_nodes_dir_surf), length(T_bdy_nodes_dir_surf));

%temperature
W_eff_nosurf = W_eff(index_u_to_T==1);
W_eff_surf = W_eff(index_u_to_T==0);

dQHva_dT = (sparse(T_up_node_ver,T_up_node_ver,L_T_veredges.*W_eff_nosurf.*dTver_dTup,T_nodes, T_nodes)+...
        sparse(T_up_node_ver,T_down_node_ver,L_T_veredges.*W_eff_nosurf.*dTver_dTdown ,T_nodes, T_nodes)-...
        (sparse(T_down_node_ver,T_up_node_ver,L_T_veredges.*W_eff_nosurf.*dTver_dTup,T_nodes, T_nodes)+...
        sparse(T_down_node_ver, T_down_node_ver,L_T_veredges.*W_eff_nosurf.*dTver_dTdown ,T_nodes, T_nodes)));
%length of subdomains
%separate surface contribution
dWeff_dxc_surf = dWeff_dxc(T_bdy_nodes_dir_surf);
dWeff_dxc = dWeff_dxc(index_u_to_T==1);

dWeff_dxs_surf = dWeff_dxs(T_bdy_nodes_dir_surf);
dWeff_dxs = dWeff_dxs(index_u_to_T==1);

dWeff_dxt_surf = dWeff_dxt(T_bdy_nodes_dir_surf);
dWeff_dxt = dWeff_dxt(index_u_to_T==1);

dQHva_dxc = L_T_veredges.*T_vert.*dWeff_dxc + T_vert.*W_eff_nosurf.*dLT_dxc(index_red_veredges);
dQHva_dxs = L_T_veredges.*T_vert.*dWeff_dxs + T_vert.*W_eff_nosurf.*dLT_dxs(index_red_veredges);
dQHva_dxt = L_T_veredges.*T_vert.*dWeff_dxt + T_vert.*W_eff_nosurf.*dLT_dxt(index_red_veredges);

dQHva_dxc_net = (accumarray(T_up_node_ver,dQHva_dxc , [T_nodes,1])-accumarray(T_down_node_ver,dQHva_dxc, [T_nodes,1]));
dQHva_dxs_net = (accumarray(T_up_node_ver, dQHva_dxs, [T_nodes,1])-accumarray(T_down_node_ver, dQHva_dxs, [T_nodes,1]));
dQHva_dxt_net = (accumarray(T_up_node_ver, dQHva_dxt, [T_nodes,1])-accumarray(T_down_node_ver, dQHva_dxt, [T_nodes,1]));

%correct for surface nodes,QH_v_a_surf = L_T(index_u_to_T==0).*W_eff(index_u_to_T==0).*T_surf;
dQHva_dxc_net(T_bdy_nodes_dir_surf) = dQHva_dxc_net(T_bdy_nodes_dir_surf) + L_T(T_bdy_nodes_dir_surf)*T_surf.*dWeff_dxc_surf + T_surf.*W_eff_surf.*dLT_dxc(T_bdy_nodes_dir_surf);
dQHva_dxs_net(T_bdy_nodes_dir_surf) = dQHva_dxs_net(T_bdy_nodes_dir_surf) + L_T(T_bdy_nodes_dir_surf)*T_surf.*dWeff_dxs_surf + T_surf.*W_eff_surf.*dLT_dxs(T_bdy_nodes_dir_surf);
dQHva_dxt_net(T_bdy_nodes_dir_surf) = dQHva_dxt_net(T_bdy_nodes_dir_surf) + L_T(T_bdy_nodes_dir_surf)*T_surf.*dWeff_dxt_surf + T_surf.*W_eff_surf.*dLT_dxt(T_bdy_nodes_dir_surf);

%sum derivatives of ver net flux QH_v = (Pe)*QH_v_a+QH_v_d;
dQHv_dH = Pe*dQHva_dH + dQHvd_dH;
dQHv_dT = Pe*dQHva_dT + dQHvd_dT;
dQHv_dxc = Pe*dQHva_dxc_net+ dQHvd_dxc;
dQHv_dxs = Pe*dQHva_dxs_net + dQHvd_dxs;
dQHv_dxt = Pe*dQHva_dxt_net + dQHvd_dxt;

%correct for flux condition at the bed in cold subdomain: net_H_verflux(T_bdy_nodes_flux) = net_H_verflux(T_bdy_nodes_flux)- nu*L_T(T_bdy_nodes_flux);
dQHv_dxc(T_bdy_nodes_flux) = dQHv_dxc(T_bdy_nodes_flux) - nu*dLT_dxc(T_bdy_nodes_flux);
dQHv_dxs(T_bdy_nodes_flux) = dQHv_dxs(T_bdy_nodes_flux) - nu*dLT_dxs(T_bdy_nodes_flux);
dQHv_dxt(T_bdy_nodes_flux) = dQHv_dxt(T_bdy_nodes_flux) - nu*dLT_dxt(T_bdy_nodes_flux);

%JACOBIAN OF HEAT CONS. fout(h_nodes+1:h_nodes+T_nodes) = Pe*net_H_horflux./Delta_sigma_T + net_H_verflux./T_Delta_eta - alpha*L_T.*S;
Dheatcons_dH = Pe*spdiags(1./Delta_sigma_T,0,T_nodes,T_nodes)*dQHh_dH + spdiags(L_T,0,T_nodes,T_nodes)* (- alpha*dS_dH) + dQHv_dH/T_Delta_eta;
Dheatcons_dT = Pe*spdiags(1./Delta_sigma_T,0,T_nodes,T_nodes)*dQHh_dT + dQHv_dT/T_Delta_eta;
Dheatcons_dxc =Pe./Delta_sigma_T.*dQHh_dxc + dQHv_dxc/T_Delta_eta + L_T.* ( - alpha*dS_dxc)+...
    dLT_dxc.* (- alpha*S);
Dheatcons_dxs = Pe./Delta_sigma_T.*dQHh_dxs + dQHv_dxs/T_Delta_eta + L_T.* ( - alpha*dS_dxs)+...
    dLT_dxs.* (- alpha*S);
Dheatcons_dxt = Pe./Delta_sigma_T.*dQHh_dxt + dQHv_dxt/T_Delta_eta + L_T.*(- alpha*dS_dxt)+...
    dLT_dxt.* (- alpha*S);
fout(h_nodes+1:h_nodes+T_nodes,:) = [Dheatcons_dH, Dheatcons_dT, Dheatcons_dxc, Dheatcons_dxs, Dheatcons_dxt];

%% FREE BOUNDARIES
% flux as computed at subtemperate nodes equals geothermal heat flux
index_up_T = T_bdy_nodes_melt(1);
index_down_T = T_bdy_nodes_melt(1)+nodes_ver;
index_up_h = parameters.n_x.c;
index_down_h = parameters.n_x.c+1;

index_edge = parameters.n_x.c;

dbedflux_up_dH = dQHvbed_dHup(index_edge);
dbedflux_up_dT = dQHvbed_dTbedup(index_edge);

dbedflux_down_dH = dQHvbed_dHdown(index_edge);
dbedflux_down_dT = dQHvbed_dTbeddown(index_edge);

fout(h_nodes+T_nodes+1,:) = [sparse(1,index_up_h, dbedflux_up_dH,1,h_nodes) sparse(1,index_up_T,dbedflux_up_dT,1,T_nodes)  sparse(1,3)]+...
    [ sparse(1,index_down_h, dbedflux_down_dH,1,h_nodes), sparse(1,index_down_T,dbedflux_down_dT,1,T_nodes) sparse(1,3)]+...
    [sparse(1,h_nodes+T_nodes) dQHvbed_dxc(index_edge) dQHvbed_dxs(index_edge) dQHvbed_dxt(index_edge)];

% basal flux such that first temperate edge satisfies basal energy budget=0;
%fout(h_nodes+T_nodes+2) =  -QH_v_edge(index_edge) + alpha*gamma*(U_SL_temp_st).^2 + nu;
index_edge = parameters.n_x.c+parameters.n_x.s;

index_h_up = parameters.n_x.c+parameters.n_x.s;
index_h_down = parameters.n_x.c+parameters.n_x.s+1;
index_T_up = T_bdy_nodes_subtemp(end);
index_T_down = T_bdy_nodes_subtemp(end)+nodes_ver;

%ice thickness
dbctemp_dU = alpha*gamma*2*(U_SL_temp_st);

dbctemp_dH_up = dQHvbed_dHup(index_edge);
dbctemp_dH_down = dQHvbed_dHdown(index_edge);

dbctemp_dH = -(sparse(1,index_h_up, dbctemp_dH_up, 1, h_nodes) + sparse(1,index_h_down, dbctemp_dH_down, 1, h_nodes)) ...
             + dbctemp_dU*dUSLstemp_dH;
         
%temperature
dbctemp_dT_up = dQHvbed_dTbedup(index_edge);
dbctemp_dT_down = dQHvbed_dTbeddown(index_edge);
dbctemp_dT = -(sparse(1,index_T_up, dbctemp_dT_up, 1, T_nodes) + sparse(1,index_T_down, dbctemp_dT_down, 1, T_nodes));

%domain bpoundaries
dbctemp_dxc = -dQHvbed_dxc(index_edge) + dbctemp_dU*dUSLstemp_dxc;
dbctemp_dxs = -dQHvbed_dxs(index_edge) + dbctemp_dU*dUSLstemp_dxs;
dbctemp_dxt = -dQHvbed_dxt(index_edge) + dbctemp_dU*dUSLstemp_dxt;

fout(h_nodes+T_nodes+2,:) = [dbctemp_dH, dbctemp_dT, dbctemp_dxc, dbctemp_dxs, dbctemp_dxt];
% flotation at the grounding line fout = H_g - r_rho*B_g
fout(h_nodes+T_nodes+3,:) = [dHg_dH.', sparse(1,T_nodes),  - r_rho*dBg_dxc, - r_rho*dBg_dxs, - r_rho*dBg_dxt];
end



