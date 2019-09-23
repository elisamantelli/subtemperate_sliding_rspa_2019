%steady state solver and and numerical stability analysis for the flow line model
%with subtemperate sliding in "Ice sheet flow with thermally activated sliding. Part 1&2"

%Elisa Mantelli, Sept 2019


st = load(['steady_state_Pe_4_b1_01_k_new_higamma.mat']);
parameters = st.parameters;

%construct computational grid
parameters.gamma = st.gamma_list(6);
parameters.n_x.c = 160;               % # nodes cold subdomain
parameters.n_x.s = 160;               % # nodes subtemp subdomain
parameters.n_x.t = 160;               % # nodes temp subdomain

[grid_h, grid_T, grid_u] = fv_grid(parameters.n_x,160);

%parameter structures
parameters.grid_h = grid_h;
parameters.grid_T = grid_T;
parameters.grid_u = grid_u;

%construct steady state flow line with subtemperate sliding starting from
%initial guess
v_in_interp = interp_sstate_v2(grid_h,grid_T,'steady_state_guess');
srch.verbose = 1;
srch.tolF = 5.10^(-8);
[v_in_new,error_flag] = Newton_v2(@network_sstate_v4,@network_sstate_jacobian_v5,v_in_interp,parameters,srch);

% Solve eigenvalue problem arising from numerical stability analysis

%options for the eigenvalue solver
parameters_eig_solver.opts.p = 150;     %dimension of the search subspace
parameters_eig_solver.opts.disp = 0;    %set to 1 for verbose eig solver
parameters_eig_solver.n_eig = 10;       %number of computed eigenvalues
parameters_eig_solver.srch_type = 'lr'; %search eigv with largest real part
parameters.k=0;

%eigenvalue solver
fout_eig = eig_solve(v_in_new, parameters, parameters_eig_solver);

%plot spectrum
figure; hold on
plot(imag(fout_eig.eig), real(fout_eig.eig),'o')

