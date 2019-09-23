%setup parplot
parplot.n_rows = 3;
parplot.n_cols = 2; 
parplot.box = 'true';
parplot.hold = 'true';
parplot.col_tickoffset = 0.05;
parplot.row_tickoffset = 0.025;
fontsize = 20;

st_s = load('stability_p_150_steady_state_Pe_4_b1_01_k_new_lowgamma_jac_v5.mat');
parameters_s = st_s.parameters;
st_f = load('stability_p_150_steady_state_Pe_4_b1_01_k_new_higamma_jac_v5.mat');
parameters_f = st_f.parameters;

st_steady_s = load('steady_state_Pe_4_b1_01_k_new_lowgamma.mat');
st_steady_f = load('steady_state_Pe_4_b1_01_k_new_higamma.mat');
vin_ss = st_steady_s.v_in_matrix(:,2);
vin_sf = st_steady_f.v_in_matrix(:,5);

%coordinate of nodes
grid_h = st_s.parameters.grid_h;
coord_sigma_h = grid_h.coor_nodes;
length_h_s = fvlength(vin_ss(end-2:end), parameters_s,'h');
Delta_sigma_h = length_h_s.Delta_sigma;
parameters = parameters_s;
epsilon = 10;
xc = vin_ss(end-2) + epsilon* st_s.xc;
xs = vin_ss(end-1) + epsilon* st_s.xs;
xt = vin_ss(end) + epsilon* st_s.xt;
x_coord_h_s = [xc*(Delta_sigma_h(1)/2:Delta_sigma_h(1): 1-Delta_sigma_h(1)/2)  xc+(xs-xc)*Delta_sigma_h(parameters.n_x.c+1)/2 ...
        xc+(xs-xc)*(3*Delta_sigma_h(parameters.n_x.c+1)/2:Delta_sigma_h(parameters.n_x.c+1): 1-Delta_sigma_h(parameters.n_x.c+1)/2) ...
        xs+(xt-xs)*Delta_sigma_h(parameters.n_x.c+parameters.n_x.s+1)/2 ...
        xs+(xt-xs)*(3*Delta_sigma_h(parameters.n_x.c+parameters.n_x.s+1)/2:Delta_sigma_h(parameters.n_x.c+parameters.n_x.s+1): 1-Delta_sigma_h(parameters.n_x.c+parameters.n_x.s+1)/2)];    

length_h_s = fvlength(vin_sf(end-2:end), parameters_f,'h');
Delta_sigma_h = length_h_s.Delta_sigma;
parameters = parameters_f;
xc = vin_sf(end-2) + epsilon* st_f.xc;
xs = vin_sf(end-1) + epsilon* st_f.xs;
xt = vin_sf(end) + epsilon* st_f.xt;
x_coord_h_f = [xc*(Delta_sigma_h(1)/2:Delta_sigma_h(1): 1-Delta_sigma_h(1)/2)  xc+(xs-xc)*Delta_sigma_h(parameters.n_x.c+1)/2 ...
        xc+(xs-xc)*(3*Delta_sigma_h(parameters.n_x.c+1)/2:Delta_sigma_h(parameters.n_x.c+1): 1-Delta_sigma_h(parameters.n_x.c+1)/2) ...
        xs+(xt-xs)*Delta_sigma_h(parameters.n_x.c+parameters.n_x.s+1)/2 ...
        xs+(xt-xs)*(3*Delta_sigma_h(parameters.n_x.c+parameters.n_x.s+1)/2:Delta_sigma_h(parameters.n_x.c+parameters.n_x.s+1): 1-Delta_sigma_h(parameters.n_x.c+parameters.n_x.s+1)/2)];
clear parameters;     

%second row: eigenfunctions
labels.ylab = '$h$';
labels.xlab = '$x$';
labels.panel_lab = '(a)';
lims.y = [min(st_s.H) max(st_s.H)];
lims.x = [0, max(x_coord_h_s)];
plotpanel(1,parplot,labels, lims, fontsize);

set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(gca,'Layer','top')

plot(x_coord_h_s, st_s.H)

labels.panel_lab = '(b)';
lims.y = [min(st_f.H) max(st_f.H)];
lims.x = [0, max(x_coord_h_f)];
plotpanel(2,parplot,labels, lims, fontsize);

set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(gca,'Layer','top')

plot(x_coord_h_f, st_f.H)

%third row: sliding velocity
labels.ylab = '$u_b$';
labels.xlab = '$x$';
labels.panel_lab = '(c)';

[~, Dslid] = sliding_v3(vin_ss,parameters_s);
dUSL_s = Dslid.dUSL;
us_s = dUSL_s*[st_s.H;st_s.T; st_s.xc;st_s.xs;st_s.xt];

lims.x = [0, max(x_coord_h_s)];
lims.y = [min(us_s), max(us_s)];
plotpanel(3,parplot,labels, lims, fontsize);

set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(gca,'Layer','top')

plot(x_coord_h_s, us_s)

[~, Dslid] = sliding_v3(vin_sf,parameters_f);
dUSL_f = Dslid.dUSL;
us_f = dUSL_f*[st_f.H;st_f.T; st_f.xc;st_f.xs;st_f.xt];

labels.panel_lab = '(d)';
lims.y = [min(us_f), max(us_f)];
lims.x = [0, max(x_coord_h_f)];
plotpanel(4,parplot,labels, lims, fontsize);

set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(gca,'Layer','top')

plot(x_coord_h_f, us_f)

%fourth row: flux
labels.ylab = '$Q$';
labels.xlab = '$x$';
labels.panel_lab = '(e)';

T_bed_nodes = st_s.parameters.grid_T.bdy_nodes.bed_nodes;  
T_bed = st_s.T(T_bed_nodes);
H = st_s.H;
T_Delta_eta = st_s.parameters.grid_T.Delta_eta; 

T_bed_0 = vin_ss(st_s.parameters.n_x.c+st_s.parameters.n_x.s+st_s.parameters.n_x.t+T_bed_nodes);
T_abed_0 = vin_ss(st_s.parameters.n_x.c+st_s.parameters.n_x.s+st_s.parameters.n_x.t+T_bed_nodes+1);
H_0 = vin_ss(st_s.parameters.n_x.c+st_s.parameters.n_x.s+st_s.parameters.n_x.t);

QH_v0 = -2/T_Delta_eta*(T_bed_0)./H_0;
Qs =-2/T_Delta_eta*(T_bed)./H_0 - QH_v0./H_0.*H ;
Qs(1:st_s.parameters.n_x.c)=0;

lims.y = [min(Qs), max(Qs)];
lims.x = [0, max(x_coord_h_s)];
plotpanel(5,parplot,labels, lims, fontsize);

set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(gca,'Layer','top')

plot(x_coord_h_s, Qs)


labels.ylab = '$Q$';
labels.xlab = '$x$';
labels.panel_lab = '(f)';

T_bed_nodes = st_f.parameters.grid_T.bdy_nodes.bed_nodes;  
T_bed = st_f.T(T_bed_nodes);
H = st_f.H;
T_Delta_eta = st_f.parameters.grid_T.Delta_eta; 

T_bed_0 = vin_sf(st_f.parameters.n_x.c+st_f.parameters.n_x.s+st_f.parameters.n_x.t+T_bed_nodes);
T_abed_0 = vin_sf(st_f.parameters.n_x.c+st_f.parameters.n_x.s+st_f.parameters.n_x.t+T_bed_nodes+1);
H_0 = vin_sf(st_f.parameters.n_x.c+st_f.parameters.n_x.s+st_f.parameters.n_x.t);

QH_v0 = -2/T_Delta_eta*(T_bed_0)./H_0;
Qf =-2/T_Delta_eta*(T_bed)./H_0 - QH_v0./H_0.*H ;
Qf(1:st_s.parameters.n_x.c)=0;

lims.y = [min(Qf), max(Qf)];
lims.x = [0, max(x_coord_h_f)];
plotpanel(6,parplot,labels, lims, fontsize);

set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(gca,'Layer','top')

plot(x_coord_h_f, Qf)




