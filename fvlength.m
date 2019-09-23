function [fout, Dfout] = fvlength(v_in, parameters, flag)
%length of cells/subdomains/horizontal edges

% input: flag identifies the grid ('h' or 'T') that is being considered. 

%output: 
%        L: list (n_nodes-by-one) of x-length of the subdomain
%        Delta_sigma: sigma-length of cells
%        Delta_x: x-length of horizontal edges
%        Delta_x_g : x-length of horizontal edge for the cell adjacent to
%        grounding line in u_grid

%unpack parameters
n_x_c = parameters.n_x.c;                             %number of horizontal nodes, cold subd
n_x_s = parameters.n_x.s;                             %number of horizontal nodes, subt subd
n_x_t = parameters.n_x.t;                             %number of horizontal nodes, temp subd

xc = v_in(1);
xs = v_in(2);
xt = v_in(3);

%set number of nodes in the vertical
if strcmp(flag,'h') == 1
    n_z = 1;
elseif strcmp(flag,'T') == 1
    n_z = parameters.grid_T.n_nodes.vert;             %number of nodes in the vertical direction
else
    error('Flag incompatible. Flag must be set to h or T');
end

%compute lengths
fout.L = [ xc*ones(n_z*(n_x_c),1); (xs-xc)*ones(n_z*(n_x_s),1) ;(xt-xs)*ones(n_z*(n_x_t),1)];%cells x
fout.Delta_sigma = [1/n_x_c*ones(n_z*(n_x_c),1); 1/n_x_s*ones(n_z*(n_x_s),1) ; 1/n_x_t*ones(n_z*(n_x_t),1) ];%cells sigma
if strcmp(flag,'h') == 1 %edges x
    fout.Delta_x = [ xc/n_x_c*ones(n_z*(n_x_c-1),1) ; 1/2*(xc/n_x_c+(xs-xc)/n_x_s)*ones(n_z,1) ; (xs-xc)/n_x_s*ones(n_z*(n_x_s-1),1) ; 1/2*((xs-xc)/n_x_s+(xt-xs)/n_x_t)*ones(n_z,1) ;...
        (xt-xs)/n_x_t*ones(n_z*(n_x_t-1),1) ];
end

%jacobian

if nargout == 2
    Dfout.dL_dxc = [ ones(n_z*(n_x_c),1); -ones(n_z*(n_x_s),1) ; zeros(n_z*(n_x_t),1)];
    Dfout.dL_dxs = [ zeros(n_z*(n_x_c),1); ones(n_z*(n_x_s),1) ; -ones(n_z*(n_x_t),1)];
    Dfout.dL_dxt = [ zeros(n_z*(n_x_c),1); zeros(n_z*(n_x_s),1) ; ones(n_z*(n_x_t),1)];
    if strcmp(flag,'h') == 1
        Dfout.dDeltax_dxc = [ 1/n_x_c*ones(n_z*(n_x_c-1),1) ; 1/2*(1/n_x_c-1/n_x_s)*ones(n_z,1) ; -1/n_x_s*ones(n_z*(n_x_s-1),1) ; 1/2*(-1/n_x_s)*ones(n_z,1) ;...
            zeros(n_z*(n_x_t-1),1)];
        Dfout.dDeltax_dxs = [ zeros(n_z*(n_x_c-1),1) ; 1/2*(1/n_x_s)*ones(n_z,1) ; 1/n_x_s*ones(n_z*(n_x_s-1),1) ; 1/2*(1/n_x_s-1/n_x_t)*ones(n_z,1) ;...
            -1/n_x_t*ones(n_z*(n_x_t-1),1) ];
        Dfout.dDeltax_dxt = [ zeros(n_z*(n_x_c-1),1) ; zeros(n_z,1) ; zeros(n_z*(n_x_s-1),1) ; 1/2*(1/n_x_t)*ones(n_z,1) ;...
            1/n_x_t*ones(n_z*(n_x_t-1),1) ];
    end
end





end