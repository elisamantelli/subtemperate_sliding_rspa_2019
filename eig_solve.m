function fout = eig_solve(v_in, parameters, parameters_eig_solver,varargin)
%compute eigenvalue and eigenfunctions for a prescribed steady state

%unpack parameters
h_nodes = parameters.grid_h.n_nodes;                  %number of nodes,h
T_nodes = parameters.grid_T.n_nodes.tot;              %number of temperature nodes

matrix = network_eigproblem_v3(v_in,parameters);
S = matrix.S;
M = matrix.M;
P = matrix.P;
list_el = matrix.list_el;

opts = parameters_eig_solver.opts;

[V,d] = eigs(S,M, parameters_eig_solver.n_eig,parameters_eig_solver.srch_type ,opts);

%extract most unstable eigenvalue
I = find(max(diag(d)));
fout.eig = diag(d);
fout.eigmax = d(I,I);

%construct eigenfunction
v_pert = V(:,I);
fout.v_pert = v_pert;

%reconstruct eliminated components and re assemble the eigenvector
v_el = P*V(:,I);
[list_el, I_sort] = sort(list_el,'ascend');
v_el = v_el(I_sort);

for j = 1:length(list_el)
    v_new = [v_pert(1:list_el(j)-1); v_el(j); v_pert(list_el(j):end)];
    v_pert = v_new;
end

%eigenfunctions
fout.H = v_pert(1:h_nodes);
fout.T = v_pert(h_nodes+1:h_nodes+T_nodes);
fout.xc = v_pert(end-2);
fout.xs = v_pert(end-1);
fout.xt = v_pert(end);

end