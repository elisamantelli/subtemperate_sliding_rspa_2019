
function [grid_h, grid_T, grid_u] = fv_grid(n_x,n_z)
grid_h = h_grid(n_x);
grid_T = T_grid(n_x,n_z);
grid_u = u_grid(n_x,n_z);
end

function fout  = h_grid(n_x)
%Generate line of n.c+n.s+n.t nodes, and return spacing, up_node,
%down_node for each edge. Describes 3 subdomains
%Modified Nov. 25th, 2016. Elisa Mantelli

%input: structure containing the fields c,s,t. Number of nodes for each
%subdomain

%define topology and boundary nodes. Include edges between subdomains, does
%NOT include edges onto lateral boundaries
fout.n_edges = sum(cell2mat(struct2cell(n_x)))-3+2;
fout.n_nodes = sum(cell2mat(struct2cell(n_x)));
fout.up_node = (1:fout.n_nodes-1)';
fout.down_node = (2:fout.n_nodes)';
fout.subtemp_slid = (n_x.c+1:(n_x.c+n_x.s))';

%number of nodes in the vertical direction
n_z=1;

%define length of edges (list n_edges-by-one)
L.sigma = [ 1/n_x.c*ones(n_z*(n_x.c-1),1) ; 1/2*(1/n_x.c+1/n_x.s)*ones(n_z,1) ; 1/n_x.s*ones(n_z*(n_x.s-1),1) ; 1/2*(1/n_x.s+1/n_x.t)*ones(n_z,1) ;...
    1/n_x.t*ones(n_z*(n_x.t-1),1)];

%define locations of nodes and cell edges
fout.coor_nodes = reshape(repmat([((1:n_x.c)'-1/2)*1/n_x.c; 1+((1:n_x.s)'-1/2)*1/n_x.s; 2+ ((1:n_x.t)'-1/2)*1/n_x.t]',n_z,1),[fout.n_nodes 1]);
fout.coor_edges = reshape(repmat([((1:n_x.c)')*1/n_x.c; ((1:n_x.s)')*1/n_x.s;  ((1:n_x.t-1)')*1/n_x.t]',n_z,1),[fout.n_edges 1]);

%identify subdomain for each node
fout.id_node = [ones(n_x.c*n_z,1) ; 2*ones(n_x.s*n_z,1) ; 3*ones(n_x.t*n_z,1)];

end

function fout  = T_grid(n_x,n_z)
%Generates rectangular grid for T describing 3 adjacent subdomains, with domain boundaries corresponding to cell boundaries.
%Vertical edges are oriented towards the ice surface, horizontal edges
%towards the grounding line. Includes list of boundary nodes and boundary edges(needed for T=0 st the cold-subtemp boundary).

% input:     n_x: structure (3 fields) including the number of horizontal nodes
%                 for each subdomain. Must be the same as n_x in grid_h
%
%            n_z: number of nodes in the vertical direction


%define topology and boundary nodes (does NOT include boundary edges) 
fout.n_edges.hor = (sum(cell2mat(struct2cell(n_x)))-3+2)*n_z;
fout.n_edges.vert = (n_z-1)*(sum(cell2mat(struct2cell(n_x))));
fout.n_nodes.tot = sum(cell2mat(struct2cell(n_x)))*n_z;
fout.n_nodes.vert = n_z;

index_nodes = (1:fout.n_nodes.tot)';
up = index_nodes; up(n_z:n_z:fout.n_nodes.tot)=[];
down = index_nodes; down(1:n_z:fout.n_nodes.tot)=[];
%define up and down nodes for each edge
fout.up_node.vert = index_nodes(up); 
fout.down_node.vert = index_nodes(down);

fout.up_node.hor = (1:fout.n_nodes.tot-n_z)';
fout.down_node.hor = (n_z+1:fout.n_nodes.tot)';

%list of nodes where Dirichlet or Neumann conditions apply
fout.bdy_nodes.flux = (1:n_z:n_x.c*n_z-n_z+1)' ;
fout.bdy_nodes.dir.surf = (n_z:n_z:fout.n_nodes.tot)';
fout.bdy_nodes.dir.bed = (n_x.c*n_z+1:n_z:fout.n_nodes.tot-n_z+1)';
fout.bdy_nodes.subtemp_slid = (n_x.c*n_z+1:n_z:(n_x.c+n_x.s)*n_z)';
fout.bdy_nodes.bed_nodes = (1:n_z:(n_x.c+n_x.s+n_x.t)).';
fout.bdy_nodes.T_melt = [(n_x.c-1)*n_z+1;(n_x.c-1)*n_z+2]; %nodes needed for cold-subtemp free boundary condition
fout.bdy_nodes.bed_nodes = (1:n_z:fout.n_nodes.tot).'; %all bed nodes
fout.bdy_nodes.inflow = (1:n_z).'; 

fout.Delta_eta = 1/n_z;

%define locations of nodes 
fout.coor_nodes.sigma = reshape(repmat([((1:n_x.c)'-1/2)*1/n_x.c; 1+((1:n_x.s)'-1/2)*1/n_x.s;  2+((1:n_x.t)'-1/2)*1/n_x.t]',n_z,1),[fout.n_nodes.tot 1]);
fout.coor_nodes.eta = reshape(repmat(((1:n_z)'-1/2)*1/n_z,1,sum(cell2mat(struct2cell(n_x)))),[fout.n_nodes.tot 1]);

%identify subdomain for each node
fout.id_node = [ones(n_x.c*n_z,1) ; 2*ones(n_x.s*n_z,1) ; 3*ones(n_x.t*n_z,1)];
end 

function fout  = u_grid(n_x,n_z)
%Generates rectangular grid for u describing 3 adjacent subdomains, with nodes on horizontal domain boundaries, but NO DIVIDE and NO GROUNDING LINE.
%Vertical edges are oriented towards the ice surface, horizontal edges
%towards the grounding line. 

% input:     n_x: structure (3 fields) including the number of horizontal nodes
%                 for each subdomain in the h grid. Must be the same as n_x
%                 in grid_h
%
%            n_z: number of nodes along the vertical in T grid. Must be the
%            same as n_z in grid_T


%define topology and boundary nodes 
fout.n_edges.hor = (sum(cell2mat(struct2cell(n_x)))-1)*(n_z+1);
fout.n_edges.vert = ((n_z)*(sum(cell2mat(struct2cell(n_x)))-1));
fout.n_nodes.tot = (sum(cell2mat(struct2cell(n_x)))-1)*(n_z+1);
fout.n_nodes.ver = n_z+1;

index_nodes = (1:fout.n_nodes.tot)';
up = index_nodes; up(n_z+1:n_z+1:fout.n_nodes.tot)=[];
down = index_nodes; down(1:n_z+1:fout.n_nodes.tot)=[];
%define up and down nodes for each edge
fout.up_node.vert = index_nodes(up); 
fout.down_node.vert = index_nodes(down);

fout.up_node.hor = (1:fout.n_nodes.tot-(n_z+1))';
fout.down_node.hor = (n_z+2:fout.n_nodes.tot)';

fout.Delta_eta = 1/n_z;

%define locations of nodes 
fout.coor_nodes.sigma = reshape(repmat([ ((1:n_x.c)')*1/n_x.c; 1+((1:n_x.s)')*1/n_x.s;  2+((1:n_x.t-1)')*1/n_x.t]',n_z+1,1),[fout.n_nodes.tot 1]);
fout.coor_nodes.eta = reshape(repmat(((0:n_z)')*1/n_z,1,(sum(cell2mat(struct2cell(n_x)))-1)),[fout.n_nodes.tot 1]);

end 


