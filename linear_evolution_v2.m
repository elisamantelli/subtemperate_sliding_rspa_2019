
function [fout,faux] = linear_evolution_v2(lin,parameters,t)

fout_disp = disp_analytic_v2(lin, parameters);
amplitude_Q = repmat(lin.fouriercoeff_Q.',[size(fout_disp.Atheta,1) 1])./fout_disp.Atheta;
amplitude_psi = repmat(lin.fouriercoeff_psi.',[size(fout_disp.Atheta,1) 1])./fout_disp.Atheta;

AQ = amplitude_Q.*fout_disp.AQ;
Atheta = amplitude_Q.*fout_disp.Atheta;
Atheta_psi = amplitude_psi.*fout_disp.Atheta;
AQ(isinf(AQ)) = 0;
Atheta(isinf(Atheta)) = 0;
    
sqrt_ice = fout_disp.sqrt_ice;
sqrt_bed = fout_disp. sqrt_bed;
sigma = fout_disp.sigma;

%compute initial condition for amplitudes independent of z
%1) amplitudes independent of ver coordinate: Q and Theta_0)
sum_AQ = nansum(AQ.*exp(sigma*t),1);
sum_Atheta = nansum(Atheta.*exp(sigma*t),1);

T_bed_init = ifft(ifftshift(sum_Atheta));
Q_init = ifft(ifftshift(sum_AQ));

%2) amplitudes dependent of ver coordinate

%temperature field
zlist_ice= parameters.grid.T.ice.coor_nodes.z(1:parameters.n_x:end);
zlist_bed = parameters.grid.T.bed.coor_nodes.z(1:parameters.n_x:end);
T_init = zeros(length(zlist_ice), length(parameters.grid.Q.coor_nodes.x));
Tbed_init = zeros(length(zlist_ice), length(parameters.grid.Q.coor_nodes.x));
for jj = 1:length(zlist_ice)
    z_ice = zlist_ice(jj);
    z_bed = zlist_bed(jj);
    
    sum_Vice = nansum((Atheta.* exp(-z_ice.*sqrt_ice)-AQ.*z_ice).*exp(sigma*t), 1);
    
    sum_Vbed = nansum(Atheta.* exp(z_bed.*sqrt_bed).*exp(sigma*t), 1);
    
    T_init(jj,:) = ifft(ifftshift(sum_Vice));
    Tbed_init(jj,:) = ifft(ifftshift(sum_Vbed));
end
%stream function and vorticity
zlist= parameters.grid.psi.coor_nodes.z(1:parameters.n_x:end);
psi_init = zeros(length(zlist), length(parameters.grid.psi.coor_nodes.x(1:parameters.grid.psi.n_nodes.hor)));
omega_init = zeros(length(zlist), length(parameters.grid.psi.coor_nodes.x(1:parameters.grid.psi.n_nodes.hor)));
%compute A in such a way to satisfy the slding law
[f_slide_Tbed, dfslide_Tbed] = regularization_old(parameters.T_bed(1), parameters);
gamma = parameters.gamma;
F0 = f_slide_Tbed; 
DF = dfslide_Tbed; 

u0 = 3*F0/(gamma + 3*F0);
tau0 = ( F0 *gamma^-1 )^(-1) *u0;
kmatrix = repmat(lin.klist, [size(Atheta_psi,1) 1]);
A =  gamma^(-1)*tau0*DF .*Atheta_psi./(sinh(2*kmatrix) -2*kmatrix - gamma^(-1).* F0 .*(-4*kmatrix.* (sinh(kmatrix)).^2));
B = (-A/2.* (1-exp(2*kmatrix)));
C = (-A);
D = (A/2.* (1-exp(-2*kmatrix)));
%A*k^2 *(-2 *k + 3 *sinh(2 k))

for jj = 1:length(zlist)
    z = zlist(jj);
   
    Vpsi = ((A+B*z).*exp(-kmatrix*z) + (C+D*z).*exp(kmatrix*z));
    
    sum_Vpsi = nansum(Vpsi.*exp(sigma*t), 1);
    
    sum_Vpsi(isinf(sum_Vpsi)) = 0;
    psi_init(jj,:) =ifft(ifftshift(sum_Vpsi));
    
    Vomega = -2*B.*exp(-kmatrix*z).*kmatrix +2*D.*exp(kmatrix*z).*kmatrix;
    sum_Vomega = nansum(Vomega.*exp(sigma*t), 1);
    
    sum_Vomega(isinf(sum_Vomega)) = 0;
    omega_init(jj,:) = ifft(ifftshift(sum_Vomega));
end

%compute analytical basal shear stress and sliding velocity
V_ub = A.*sinh(2*kmatrix) -2*kmatrix.*A;
V_taub = -4*kmatrix.* A.*(sinh(kmatrix)).^2;
V_vortictyflux = 2*kmatrix.^2.*(B + D);

sum_Vub = nansum(V_ub.*exp(sigma*t), 1);
sum_Vtaub = nansum(V_taub.*exp(sigma*t), 1);
sum_Vvortictyflux = nansum(V_vortictyflux.*exp(sigma*t), 1);

sum_Vub(isinf(sum_Vub)) = 0;
ub_init =ifft(ifftshift(sum_Vub));

sum_Vtaub(isinf(sum_Vtaub)) = 0;
taub_init =ifft(ifftshift(sum_Vtaub));

sum_Vvortictyflux(isinf(sum_Vvortictyflux)) = 0;
vorticityflux_init =ifft(ifftshift(sum_Vvortictyflux));

fout = [reshape(psi_init.',[],1); reshape(omega_init.',[],1); reshape(T_init.',[],1);reshape(Tbed_init.',[],1); Q_init.'; T_bed_init.'];
% faux.ub = ub_init;
% faux.taub = taub_init;

%check sliding law
[f_slide_Tbed, dfslide_Tbed] = regularization_old(parameters.T_bed(1), parameters);
gamma = parameters.gamma;
F0 = f_slide_Tbed; 
DF = dfslide_Tbed; 

u0 = 3*F0/(gamma + 3*F0);
tau0 = ( F0 *gamma^-1 )^(-1) *u0;

residual = ub_init - gamma^(-1)*(F0 * taub_init + tau0*DF .*T_bed_init);%this is not going to be zero because T_bed_init is computed at T_nodes.
faux.ub = u0 + ub_init;% 
faux.taub = tau0+ taub_init;
faux.vorticityfluxbed = -3*(1-u0) + vorticityflux_init;%
%3*(1-u0)*(1-zlist(end)) + omega_init(end,:) - faux.taub;%
faux.u0 = u0;
faux.vor0 = -3*(1-u0);
