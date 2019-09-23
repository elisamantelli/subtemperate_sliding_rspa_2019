function fout = interp_sstate(grid_h,grid_T,filename)

st = load([filename '.mat']);
parameters = st.parameters;
parameters.gamma = st.gamma_list(6);
v_in = st.v_in_matrix(:,5);

h_nodes = parameters.grid_h.n_nodes;                                       %number of nodes,h
T_nodes = parameters.grid_T.n_nodes.tot;                                   %number of temperature nodes
nodes_ver = parameters.grid_T.n_nodes.vert;

%unpack input variable v_in
H = v_in(1:h_nodes);
T = v_in(h_nodes+1:h_nodes+T_nodes);
T_matrix = reshape(T,[nodes_ver,h_nodes]);

xc = v_in(h_nodes+T_nodes+1);
xs = v_in(h_nodes+T_nodes+2);
xt = v_in(h_nodes+T_nodes+3);

%coordinate of nodes in the solution to be interpolated
grid_h_90 = parameters.grid_h;
coord_sigma_h_90 = grid_h_90.coor_nodes;

grid_T_90 = parameters.grid_T;
coord_sigma_T_90 = reshape(grid_T_90.coor_nodes.sigma,[nodes_ver,h_nodes]);
coord_eta_T_90 = reshape(grid_T_90.coor_nodes.eta,[nodes_ver,h_nodes]);

%coordinates of nodes where initial guess is computed
coord_sigma_h = grid_h.coor_nodes;
coord_sigma_T = reshape(grid_T.coor_nodes.sigma,[grid_T.n_nodes.vert,length(coord_sigma_h) ]);
coord_eta_T = reshape(grid_T.coor_nodes.eta,[grid_T.n_nodes.vert,length(coord_sigma_h) ]);

%interpolation
h_guess = interp1(coord_sigma_h_90,H,coord_sigma_h,'spline');
T_guess = interp2(coord_sigma_T_90,coord_eta_T_90, T_matrix ,coord_sigma_T,coord_eta_T,'spline');

fout =[h_guess; reshape(T_guess,[size(T_guess,1)*size(T_guess,2),1]); xc; xs; xt]; 

end
