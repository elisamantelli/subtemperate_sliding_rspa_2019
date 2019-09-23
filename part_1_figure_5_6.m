
st = load(['steady_state_Pe_4_b1_01_k_new_lowPe.mat']);

parameters = st.parameters;
v_in = real(st.v_in_matrix(:,1));
parameters.Pe = st.Pe_list(1);

%recompute grid
[grid_h, grid_T, grid_u] = fv_grid(parameters.n_x,parameters.n_x.c);

%parameter structures
parameters.grid_h = grid_h;
parameters.grid_T = grid_T;
parameters.grid_u = grid_u;

h_nodes = parameters.grid_h.n_nodes;                                       %number of nodes,h
h_n_edges = parameters.grid_h.n_edges; 
T_nodes = parameters.grid_T.n_nodes.tot;                                   %number of temperature nodes
nodes_ver = parameters.grid_T.n_nodes.vert;
T_n_edges_hor = parameters.grid_T.n_edges.hor;                             %number of horizontal edges, T

%destination grid parameters
T_bed_nodes = parameters.grid_T.bdy_nodes.bed_nodes;                       %list of bed node
T_Delta_eta = parameters.grid_T.Delta_eta;                                 %scalar, vertical spacing between cell centres

%unpack input variable v_in
H = v_in(1:h_nodes);
T = v_in(h_nodes+1:h_nodes+T_nodes);
xc = v_in(end-2);
xs = v_in(end-1);
xt = v_in(end);

length_h = fvlength([xc;xs;xt], parameters,'h');
Delta_sigma_h = length_h. Delta_sigma;

%velocity field
[~, faux] = network_sstate_v4(v_in,parameters);
w = faux.w_nodes;
u = faux.u_nodes;
w_matrix = reshape(w, [nodes_ver, h_nodes]);


parameters.w.body = faux.w_eff.body;
parameters.w.surf = faux.w_eff.surf;
parameters.w.bed = faux.w_eff.bed;
parameters.u.body = faux.u_edges;
parameters.u.g = faux.u_g;

%compute stream lines for given geometry
v_in_stream = T;
parameters.v_in = v_in;

%coordinate of nodes in the solution to be interpolated
grid_h = parameters.grid_h;
coord_sigma_h = grid_h.coor_nodes;

x_coord_h = [xc*(Delta_sigma_h(1)/2:Delta_sigma_h(1): 1-Delta_sigma_h(1)/2)  xc+(xs-xc)*Delta_sigma_h(parameters.n_x.c+1)/2 ...
        xc+(xs-xc)*(3*Delta_sigma_h(parameters.n_x.c+1)/2:Delta_sigma_h(parameters.n_x.c+1): 1-Delta_sigma_h(parameters.n_x.c+1)/2) ...
        xs+(xt-xs)*Delta_sigma_h(parameters.n_x.c+parameters.n_x.s+1)/2 ...
        xs+(xt-xs)*(3*Delta_sigma_h(parameters.n_x.c+parameters.n_x.s+1)/2:Delta_sigma_h(parameters.n_x.c+parameters.n_x.s+1): 1-Delta_sigma_h(parameters.n_x.c+parameters.n_x.s+1)/2)];    

grid_T = parameters.grid_T;
coord_sigma_T = reshape(grid_T.coor_nodes.sigma,[nodes_ver,h_nodes]);
coord_eta_T = reshape(grid_T.coor_nodes.eta,[nodes_ver,h_nodes]);

srch.verbose=1;
parameters.T_surf_adv = -x_coord_h;
stream = Newton_v2(@network_sstate_streamlines_v2,@network_sstate_streamline_jacobian_v2,v_in_stream,parameters,srch);



figure;
B=@(x) parameters. bed. b0+parameters. bed. b1.*x;

parplot.n_rows = 2;
parplot.n_cols = 2; 
parplot.box = 'true';
parplot.hold = 'true';
parplot.col_tickoffset = 0.05;
parplot.row_tickoffset = 0.025;
fontsize = 20;
%VELOCITY FIELD AND STREAMLINES

%PANEL 1: VELOCITY
labels.xlab = '$x$';
labels.ylab = '$z$';
labels.panel_lab = '(a)';
lims.x = [0,8.5];
lims.y = [-1,3];
plotpanel(1,parplot,labels,lims, fontsize);

set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(gca,'Layer','top')

H_matrix = repmat(H.', [nodes_ver,1]);
X = repmat(x_coord_h,[nodes_ver, 1]);
Z = reshape(coord_eta_T, [nodes_ver,h_nodes]).*(H_matrix) - B(X) ;
U_matrix = reshape(u, [nodes_ver, h_nodes]);

map = brewermap([],'YlGn');
colormap(gca, map)
contourf(real(X(:,1:end-2)), real(Z(:,1:end-2)),U_matrix(:,1:end-2),100, 'LineColor', 'none')
hold on

contour(real(X(:,1:end-2)), real(Z(:,1:end-2)),real(U_matrix(:,1:end-2)),0.1:0.5:max(max(U_matrix(:,1:end-2))), '-', 'LineWidth', 1, 'LineColor', [0.83 0.83 0.83]);

hold on
plot([x_coord_h(320) x_coord_h(320)], [H(320)-B(x_coord_h(320)) -B(x_coord_h(320)) ], ':', [x_coord_h(640) x_coord_h(640)], [H(640)-B(x_coord_h(640)) -B(x_coord_h(640)) ], ':','LineWidth', 1, 'Color', [0.83 0.83 0.83]  )

c = colorbar('eastoutside');
c.Label.String = '$u$';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize', 20)

%PANEL 2: STREAMLINES
labels.xlab = '$x$';
labels.ylab = '$z$';
lims.x = [0,8.5];
lims.y = [-1,3];
labels.panel_lab = '(c)';
plotpanel(2,parplot,labels,lims, fontsize);
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

map = brewermap([],'YlGn');
colormap(gca, map)
S = reshape(stream, [nodes_ver,h_nodes]);
contour(real(X(:,1:end-2)), real(Z(:,1:end-2)),real(S(:,1:end-2)),[ -linspace(0.01, max(x_coord_h),15)], '-', 'LineWidth', 1);
hold on
plot(x_coord_h, H-B(x_coord_h.'),'-k', x_coord_h, -B(x_coord_h),'-k', 'LineWidth', 1 )

plot([x_coord_h(320) x_coord_h(320)], [H(320)-B(x_coord_h(320)) -B(x_coord_h(320)) ], ':', [x_coord_h(640) x_coord_h(640)], [H(640)-B(x_coord_h(640)) -B(x_coord_h(640)) ], ':','LineWidth', 1, 'Color', 'k'  )
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize', 20)

%PANEL 3: TEMPERATURE FIELD
labels.xlab = '$x$';
labels.ylab = '$z$';
lims.x = [0,8.5];
lims.y = [-1,3];
labels.panel_lab = '(b)';
plotpanel(3,parplot,labels,lims, fontsize);
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(gca,'Layer','top')

map2 = brewermap([],'*YlGnBu');
colormap(gca, map2)

T_matrix = reshape(T, [nodes_ver,h_nodes]);
contourf(real(X(:,1:end-2)), real(Z(:,1:end-2)),T_matrix(:,1:end-2),100, 'LineColor', 'none')
set(gca,'FontSize', 20)
hold on

contour(real(X(:,1:end-2)), real(Z(:,1:end-2)),real(T_matrix(:,1:end-2)),-0.2:-0.2:min(min(T_matrix(:,1:end-2))), '-', 'LineWidth', 1, 'LineColor', [0.83 0.83 0.83]);
hold on
contour(real(X(:,1:end-2)), real(Z(:,1:end-2)),real(T_matrix(:,1:end-2)),[-2e-03 -3e-03], '-', 'LineWidth', 3, 'LineColor', [0.83 0.83 0.83]);

plot([x_coord_h(320) x_coord_h(320)], [H(320)-B(x_coord_h(320)) -B(x_coord_h(320)) ], ':', [x_coord_h(640) x_coord_h(640)], [H(640)-B(x_coord_h(640)) -B(x_coord_h(640)) ], ':','LineWidth', 1, 'Color', [0.83 0.83 0.83]  )

c = colorbar('eastoutside');
c.Label.String = '$T$';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
set(gca,'TickLabelInterpreter','latex')

%PANEL 4: FLUX AND SLIDING VELOCITY

labels.xlab = 'x';
labels.ylab = '$u_b$';
lims.y = [-.05 3.6];
labels.panel_lab = '(d)';
hout = plotpanel(4,parplot,labels,lims, fontsize);
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

slid = sliding_v3(v_in,parameters);
U_SL = slid.U_SL;
QH_v = - 1/(T_Delta_eta/2)*(T(T_bed_nodes))./H;
QH_v(1:parameters.n_x.c) = parameters.nu;

c = brewermap(10,'Blues');

plot(x_coord_h(2:end),U_SL, '-','Color',c(5,:), 'LineWidth', 1)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize', 20)
%set(gca,'ycolor',c(5,:));
hold on
plot(x_coord_h,QH_v, '-','Color',c(9,:), 'LineWidth', 1)
ll = legend('sliding velocity', 'bed flux');
ll.Interpreter = 'latex';
ll.AutoUpdate = 'off';
plot([x_coord_h(320) x_coord_h(320)], [-.1 3.5], ':', [x_coord_h(640) x_coord_h(640)], [-.1 3.5], ':','LineWidth', 1, 'Color', 'k'  )
set(gca,'FontSize', 20)

hold on
plot(x_coord_h(2:end),-QH_v(2:end)+parameters.nu + parameters.alpha *parameters.gamma*U_SL.^2, '-','Color',c(9,:), 'LineWidth', 1)
ll = legend('sliding velocity', 'bed flux');
ll.Interpreter = 'latex';
ll.AutoUpdate = 'off';
plot([x_coord_h(320) x_coord_h(320)], [-.1 3.5], ':', [x_coord_h(640) x_coord_h(640)], [-.1 3.5], ':','LineWidth', 1, 'Color', 'k'  )
set(gca,'FontSize', 20)



%PANEL 5: LENGTH AS A FUNCTION OF FRICTION COEFFICIENT

figure;
parplot.n_rows = 1;
parplot.n_cols = 2; 
parplot.box = 'true';
parplot.hold = 'true';
parplot.col_tickoffset = 0.05;
parplot.row_tickoffset = 0.025;
fontsize = 20;


labels.xlab = '$\gamma$';
labels.ylab = ['L/x_g'];
labels.panel_lab = '(a)';
lims.x = [0,20];
lims.y = [0 1];
hout = plotpanel(1,parplot,labels,lims, fontsize);

st = load(['steady_state_Pe_4_b1_01_k_new_lowgamma.mat']);
gamma_list_small = st.gamma_list(2:end-2);
v_in_matrix = st.v_in_matrix(:,1:end-3); 
L_c = (v_in_matrix(end-2,1:end-1))./v_in_matrix(end,1:end-1);
L_s = (v_in_matrix(end-1,1:end-1)-v_in_matrix(end-2,1:end-1))./v_in_matrix(end,1:end-1);
L_t = (v_in_matrix(end,1:end-1)-v_in_matrix(end-1,1:end-1))./v_in_matrix(end,1:end-1);

L_small = [L_c;L_s;L_t;];

st = load(['steady_state_Pe_4_b1_01_k_new_higamma.mat']);
gamma_list_large = st.gamma_list;
v_in_matrix = st.v_in_matrix(:,2:end-1); 
L_c = (v_in_matrix(end-2,:))./v_in_matrix(end,:);
L_s = (v_in_matrix(end-1,:)-v_in_matrix(end-2,:))./v_in_matrix(end,:);
L_t = (v_in_matrix(end,:)-v_in_matrix(end-1,:))./v_in_matrix(end,:);

L_large = [L_c;L_s;L_t;];

st = load(['steady_state_Pe_4_b1_01_k_new_hihigamma.mat']);
gamma_list_llarge = st.gamma_list(1:51);
v_in_matrix = st.v_in_matrix(:,1:52); 
L_c = (v_in_matrix(end-2,1:end-1))./v_in_matrix(end,1:end-1);
L_s = (v_in_matrix(end-1,1:end-1)-v_in_matrix(end-2,1:end-1))./v_in_matrix(end,1:end-1);
L_t = (v_in_matrix(end,1:end-1)-v_in_matrix(end-1,1:end-1))./v_in_matrix(end,1:end-1);

L_llarge = [L_c;L_s;L_t];

L= [flip(L_small(:,2:end),2) L_large(:,3:end-1) L_llarge];
gamma_list = [flip(gamma_list_small(2:end)) gamma_list_large(2:end-2) gamma_list_llarge];
set(0,'DefaultAxesColorOrder',brewermap(3,'Paired'))
for n = 1:3
    gamma_int = linspace(min(gamma_list), max(gamma_list),30);
    L_int = interp1(gamma_list,L(n,:),gamma_int);
    c = plot(gamma_int,L_int, '.-','LineWidth',1, 'MarkerSize',15);
    hold all
end
set(gca,'TickLabelInterpreter','latex')

%PANEL 6: LENGTH AS A FUNCTION OF Pe
labels.xlab = '$Pe$';
labels.ylab = ['L/x_g'];
labels.panel_lab = '(b)';
lims.x = [1,10.5];
lims.y = [0 1];
hout = plotpanel(2,parplot,labels,lims, fontsize);

st = load(['steady_state_Pe_4_b1_01_k_new_lowPe.mat']);
gamma_list_small = st.Pe_list(1:end-1);
v_in_matrix = st.v_in_matrix(:,2:end-1); 
L_c = (v_in_matrix(end-2,1:end-1))./v_in_matrix(end,1:end-1);
L_s = (v_in_matrix(end-1,1:end-1)-v_in_matrix(end-2,1:end-1))./v_in_matrix(end,1:end-1);
L_t = (v_in_matrix(end,1:end-1)-v_in_matrix(end-1,1:end-1))./v_in_matrix(end,1:end-1);

L_small = [L_c;L_s;L_t;];

st = load(['steady_state_Pe_4_b1_01_k_new_hiPe.mat']);
gamma_list_large = st.Pe_list(1:end-1);
v_in_matrix = st.v_in_matrix(:,2:end-1); 
L_c = (v_in_matrix(end-2,1:end-1))./v_in_matrix(end,1:end-1);
L_s = (v_in_matrix(end-1,1:end-1)-v_in_matrix(end-2,1:end-1))./v_in_matrix(end,1:end-1);
L_t = (v_in_matrix(end,1:end-1)-v_in_matrix(end-1,1:end-1))./v_in_matrix(end,1:end-1);

L_large = [L_c;L_s;L_t;];

st = load(['steady_state_Pe_4_b1_01_k_new_hihiPe.mat']);
gamma_list_llarge = st.Pe_list(2:10);
v_in_matrix = st.v_in_matrix(:,2:end-1); 
L_c = (v_in_matrix(end-2,1:end-1))./v_in_matrix(end,1:end-1);
L_s = (v_in_matrix(end-1,1:end-1)-v_in_matrix(end-2,1:end-1))./v_in_matrix(end,1:end-1);
L_t = (v_in_matrix(end,1:end-1)-v_in_matrix(end-1,1:end-1))./v_in_matrix(end,1:end-1);

L_llarge = [L_c(:,1:9);L_s(:,1:9);L_t(:,1:9)];


L= [flip(L_small,2) L_large(:,2:end) L_llarge(:,2:end)];
gamma_list = [flip(gamma_list_small(2:end)) gamma_list_large(2:end) gamma_list_llarge(2:end)];
set(0,'DefaultAxesColorOrder',brewermap(3,'Paired'))
for n = 1:3
    gamma_int = linspace(min(gamma_list), max(gamma_list),30);
    L_int = interp1(gamma_list,L(n,:),gamma_int);
    c = plot(gamma_int(1:end),L_int, '.-','LineWidth',1, 'MarkerSize',15);
    hold all
end
set(gca,'TickLabelInterpreter','latex')


