function network_script_cluster (T_bed, choice, gamma,bd_x,n_x, bd_z, n_z, dt,t_plot, stepmax,t_init, step_init, index_wavelength, scale)
                              
%Updated Nov 16th. Uses second order accurate biharmonic solver for both
%steady state and non linear problem. 

parameters.choice = choice; 
parameters.index_wavelength = index_wavelength;
%consruct grids
%boundaries of the box
parameters.bd_x = bd_x;                %horizontal width
parameters.bd_z = bd_z;               
parameters.n_x = n_x;         %must be even
parameters.n_z = n_z;
parameters.n = 1;             %vertical grid refinement(must be >=1, 1 is a even grid)

grid = fv_grid_periodic(parameters.n_x,parameters.n_z,parameters.bd_x, parameters.bd_z, parameters);

%parameter structures
parameters.grid = grid;

parameters.alpha = 1;                               % strain heating parameter
parameters.nu = 1;                                  % basal heating parameter
parameters.gamma = gamma;                               % >3 to ensure that dhearing is larger than sliding
parameters.u_temp = 1;                            %3/(parameters.gamma + 3);
parameters.delta = 0.1;
parameters.T_bed = T_bed*ones(length(parameters.grid.Q.coor_nodes.x),1);

parameters.reg.epsilon_f = parameters.delta;
parameters.reg.epsilon =  parameters.reg.epsilon_f/5;

%INITIAL CONDITION

%STEP 1: CONSTRUCT STEADY STATE . 
%VELOCITY FIELD AT PRESCRIBED BED TEMPERATURE (good only for steady
%state)
srchparams.verbose = 1;
srchparams.toldelta = 1e-07;
srchparams.tolF = 1e-05;

steady_state_a = steady_state_analytic(parameters); %for more complicated initial conditions make sure that the numeric steady state is used.
steady_state_numeric =  Newton_v2(@network_subtemp_slab_steady,@network_subtemp_slab_steady_jacobian,  steady_state_a(1:end-n_x),parameters,srchparams);
steady_state_numeric = [steady_state_numeric; parameters.T_bed];


%STEP 3: CONSTRUCT PERTURBATION AS SUPERPOSITION OF NORMAL MODES
fshift = 2*pi/bd_x*(-parameters.n_x/2:parameters.n_x/2-1);
lin.klist = fshift;
k = fshift(length(fshift)/2+1+index_wavelength); 
parameters.scale = scale;
parameters.forcing_psinodes = scale*sin(k*parameters.grid.psi.coor_nodes.x(1:n_x));
parameters.forcing_Qnodes = scale*sin(k*parameters.grid.Q.coor_nodes.x);
%fourier transform the forcing
ypsi = fft(parameters.forcing_psinodes);
yQ = fft(parameters.forcing_Qnodes);
%compute wavenumbers and amplitudes centered on k = 0;
ypsishift = fftshift(ypsi);
yQshift = fftshift(yQ);
lin.fouriercoeff_psi = ypsishift;%yshift(index:end).';%
lin.fouriercoeff_Q = yQshift;
v_in_pert = linear_evolution_v2(lin,parameters,t_init);
%4) SUM STEADY STATE AND PERTURBATION
%scale_pert = scale/max(v_in_pert);
v_in = steady_state_numeric + real(v_in_pert);

parameters.lin = lin;
%% TIME STEPPING

parameters.dt = dt;
parameters.tplot = t_plot;
parameters.stepmax = stepmax;
parameters.v_in_prev = v_in;

parameters.filename = ['slab_Tbed_' num2str(parameters.T_bed(1)) '_choice_' choice '_delta_' num2str(parameters.delta) '_alpha_' num2str(parameters.alpha)...
    '_gamma_' num2str(parameters.gamma) '_nu_' num2str(parameters.nu) 'bd_x' num2str(bd_x) '_indexk_' num2str(index_wavelength) '_n_x_' num2str(parameters.n_x) '_dt_' num2str(parameters.dt)];

%initialization
step = 1;
step_in = step_init;
parameters.t_init = t_init;
t = t_init;

%save initial condition
psi_nodes = parameters.grid.psi.n_nodes.tot;
T_nodes = parameters.grid.T.ice.n_nodes.tot;
Q_nodes = parameters.grid.Q.n_nodes.tot;

fout.psi(:,step) = v_in(1:psi_nodes);
fout.omega(:,step)  = v_in(psi_nodes+1:2*psi_nodes);
fout.T(:,step) = v_in(2*psi_nodes+1:2*psi_nodes+T_nodes);
fout.Tb(:,step) = v_in(2*psi_nodes+T_nodes+1:2*psi_nodes+2*T_nodes);
fout.Q(:,step) = v_in(2*psi_nodes+2*T_nodes+1:2*psi_nodes+2*T_nodes+Q_nodes);
fout.T_bed(:,step) = v_in(2*psi_nodes+2*T_nodes+Q_nodes+1:end);
fout.error.tot(step) = 0;


while step <= parameters.stepmax-1
    t = t+parameters.dt;
    display(['t = ' num2str(t)])
    tic
    parameters.v_in_prev = v_in;
    [v_in,error_flag] = Newton_v2(@network_subtemp_slab_timedep,@network_subtemp_slab_timedep_jacobian,v_in,parameters,srchparams);
    toc
    
    %check deviation from linear solution
    if strcmp(choice, 'mono_small') == 1 || strcmp(choice, 'random') == 1|| strcmp(choice, 'wn') == 1
        v_in_pert = linear_evolution_v2(lin,parameters,t);
        v_in_analytic = steady_state_a + real(v_in_pert);
        
        analitic.error_psi = norm(v_in_analytic(1:psi_nodes))*sqrt(parameters.bd_x/parameters.n_x)*sqrt(1/parameters.n_z.psi);
        analitic.error_omega = norm(v_in_analytic(psi_nodes+1:2*psi_nodes))*sqrt(parameters.bd_x/parameters.n_x)*sqrt(1/parameters.n_z.psi);
        analitic.error_T = norm(v_in_analytic(2*psi_nodes+1:2*psi_nodes+T_nodes))*sqrt(parameters.bd_x/parameters.n_x)*sqrt(parameters.bd_z.ice/parameters.n_z.Tice);
        analitic.error_Tb = norm(v_in_analytic(2*psi_nodes+T_nodes+1:2*psi_nodes+2*T_nodes))*sqrt(parameters.bd_x/parameters.n_x)*sqrt(-parameters.bd_z.bed/parameters.n_z.Tbed);
        analitic.error_Q = norm(v_in_analytic(2*psi_nodes+2*T_nodes+1:2*psi_nodes+2*T_nodes+Q_nodes))*sqrt(parameters.bd_x/parameters.n_x);
        analitic.error_Tbed = norm(v_in_analytic(2*psi_nodes+2*T_nodes+Q_nodes+1:2*psi_nodes+2*T_nodes+2*Q_nodes))*sqrt(parameters.bd_x/parameters.n_x);
        analitic.error_tot = analitic.error_psi + analitic.error_omega + analitic.error_T + analitic.error_Tb + analitic.error_Q + analitic.error_Tbed;
        
        numeric.error_psi = norm(v_in(1:psi_nodes))*sqrt(parameters.bd_x/parameters.n_x)*sqrt(1/parameters.n_z.psi);
        numeric.error_omega = norm(v_in(psi_nodes+1:2*psi_nodes))*sqrt(parameters.bd_x/parameters.n_x)*sqrt(1/parameters.n_z.psi);
        numeric.error_T = norm(v_in(2*psi_nodes+1:2*psi_nodes+T_nodes))*sqrt(parameters.bd_x/parameters.n_x)*sqrt(parameters.bd_z.ice/parameters.n_z.Tice);
        numeric.error_Tb = norm(v_in(2*psi_nodes+T_nodes+1:2*psi_nodes+2*T_nodes))*sqrt(parameters.bd_x/parameters.n_x)*sqrt(-parameters.bd_z.bed/parameters.n_z.Tbed);
        numeric.error_Q = norm(v_in(2*psi_nodes+2*T_nodes+1:2*psi_nodes+2*T_nodes+Q_nodes))*sqrt(parameters.bd_x/parameters.n_x);
        numeric.error_Tbed = norm(v_in(2*psi_nodes+2*T_nodes+Q_nodes+1:2*psi_nodes+2*T_nodes+2*Q_nodes))*sqrt(parameters.bd_x/parameters.n_x);
        numeric.error_tot = numeric.error_psi + numeric.error_omega + numeric.error_T + numeric.error_Tb + numeric.error_Q + numeric.error_Tbed;
        
        diff.error_psi = norm(v_in(1:psi_nodes)-v_in_analytic(1:psi_nodes))*sqrt(parameters.bd_x/parameters.n_x)*sqrt(1/parameters.n_z.psi);
        diff.error_omega = norm(v_in(psi_nodes+1:2*psi_nodes)-v_in_analytic(psi_nodes+1:2*psi_nodes))*sqrt(parameters.bd_x/parameters.n_x)*sqrt(1/parameters.n_z.psi);
        diff.error_T = norm(v_in(2*psi_nodes+1:2*psi_nodes+T_nodes)-v_in_analytic(2*psi_nodes+1:2*psi_nodes+T_nodes))*sqrt(parameters.bd_x/parameters.n_x)*sqrt(parameters.bd_z.ice/parameters.n_z.Tice);
        diff.error_Tb = norm(v_in(2*psi_nodes+T_nodes+1:2*psi_nodes+2*T_nodes)-v_in_analytic(2*psi_nodes+T_nodes+1:2*psi_nodes+2*T_nodes))*sqrt(parameters.bd_x/parameters.n_x)*sqrt(-parameters.bd_z.bed/parameters.n_z.Tbed);
        diff.error_Q = norm(v_in(2*psi_nodes+2*T_nodes+1:2*psi_nodes+2*T_nodes+Q_nodes)-v_in_analytic(2*psi_nodes+2*T_nodes+1:2*psi_nodes+2*T_nodes+Q_nodes))*sqrt(parameters.bd_x/parameters.n_x);
        diff.error_Tbed = norm(v_in(2*psi_nodes+2*T_nodes+Q_nodes+1:2*psi_nodes+2*T_nodes+2*Q_nodes)-v_in_analytic(2*psi_nodes+2*T_nodes+Q_nodes+1:2*psi_nodes+2*T_nodes+2*Q_nodes))*sqrt(parameters.bd_x/parameters.n_x);
        diff.error_tot = diff.error_psi + diff.error_omega + diff.error_T + diff.error_Tb + diff.error_Q + diff.error_Tbed;
        
        
        error_rel = diff.error_tot/analitic.error_tot;
        error_analytic = analitic.error_tot;
        error_numeric = numeric.error_tot;
        disp(['error = ' num2str(error_rel)])
    else
        error_rel = NaN;
    end
    
    if error_flag
        warning('Iteration failure: convergence not achieved')
        t = t-parameters.dt;
        parameters.dt = parameters.dt/2;
        disp(['dt = ' num2str(parameters.dt)])
    end
    
    %saving 
    if t-parameters.t_init >= parameters.tplot*(step-step_in+1)
        step = step+1;
        display(['step # ' num2str(step) ' of ' num2str(parameters.stepmax)])
        
        fout.t(step) = t;
        
        fout.psi(:,step) = v_in(1:psi_nodes);
        fout.omega(:,step)  = v_in(psi_nodes+1:2*psi_nodes);
        fout.T(:,step) = v_in(2*psi_nodes+1:2*psi_nodes+T_nodes);
        fout.Tb(:,step) = v_in(2*psi_nodes+T_nodes+1:2*psi_nodes+2*T_nodes);
        fout.Q(:,step) = v_in(2*psi_nodes+2*T_nodes+1:2*psi_nodes+2*T_nodes+Q_nodes);
        fout.T_bed(:,step) = v_in(2*psi_nodes+2*T_nodes+Q_nodes+1:end);
        fout.error.rel(step) = error_rel;
        fout.error.analytic(step) = error_analytic;
        fout.error.numeric(step) = error_numeric;
        if floor(step/10) == step/10
            save([parameters.filename '.mat'],'fout','parameters','-v7.3')
        end
        
        
    end
    
    %stop when freezing happens
    [~, faux] = network_subtemp_slab_timedep(v_in,parameters);
    umin = min(faux.u_bed);
    if umin < 0.005 
        fout.t(step) = t;
        
        psi_nodes = parameters.grid.psi.n_nodes.tot;
        T_nodes = parameters.grid.T.ice.n_nodes.tot;
        Q_nodes = parameters.grid.Q.n_nodes.tot;
        
        fout.psi(:,step) = v_in(1:psi_nodes);
        fout.omega(:,step)  = v_in(psi_nodes+1:2*psi_nodes);
        fout.T(:,step) = v_in(2*psi_nodes+1:2*psi_nodes+T_nodes);
        fout.Tb(:,step) = v_in(2*psi_nodes+T_nodes+1:2*psi_nodes+2*T_nodes);
        fout.Q(:,step) = v_in(2*psi_nodes+2*T_nodes+1:2*psi_nodes+2*T_nodes+Q_nodes);
        fout.T_bed(:,step) = v_in(2*psi_nodes+2*T_nodes+Q_nodes+1:end);
        save([parameters.filename '.mat'],'fout','parameters','-v7.3')
        
        break
    end
        
end










