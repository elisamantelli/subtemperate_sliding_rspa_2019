%setup parplot
parplot.n_rows = 1;
parplot.n_cols = 2; 
parplot.box = 'true';
parplot.hold = 'true';
parplot.col_tickoffset = 0.05;
parplot.row_tickoffset = 0.025;
fontsize = 20;

st_s = load('stability_gamma_17.mat');
st_f = load('stability_gamma_07.mat');

st_steady_s = load('steady_state_Pe_4_b1_01_k_new_lowgamma.mat');
st_steady_f = load('steady_state_Pe_4_b1_01_k_new_higamma.mat');
vin_ss = st_steady_s.v_in_matrix(:,1);
vin_sf = st_steady_f.v_in_matrix(:,1);
    
%first panel: spectrum 
parplot.semilogy = 1;
labels.ylab = '$Re[\sigma]$';
labels.xlab = '$Im[\sigma]$';
lims.y = [0 300];
lims.x = [-5,5];
labels.panel_lab = '(a)';
plotpanel(1,parplot,labels,lims,fontsize);

% set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
% set(gca,'Layer','top')

plot(imag(st_s.fout_eig_160.eig),real(st_s.fout_eig_160.eig), 'o','MarkerSize', 10)

lims.y = [10000 50000];
plotpanel(2,parplot,labels,lims, fontsize);

set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(gca,'Layer','top')

plot(imag(st_f.fout_eig_160.eig),real(st_f.fout_eig_160.eig),  'o','MarkerSize', 10)

parplot.semilogy = 1;


