
function grid = fv_grid_periodic(n_x,n_z,bd_x, bd_z, parameters)
grid.Q  = grid_1D(n_x,bd_x);
n = parameters.n;
grid.T.ice  = grid_T_bl(n_x,n_z.Tice,bd_x, bd_z.ice, n);
grid.T.bed  = grid_T_bl(n_x,n_z.Tice,bd_x, bd_z.bed, n);
grid.psi = grid_psi(n_x,n_z.psi,bd_x,1, n);

end


function fout  = grid_psi(n_x,n_z,bd_x, bd_z, n)
%Generates rectangular grid for T in the domain (0<z<z_p, 0<x<xp), with domain boundaries corresponding to cell edges.
%Vertical edges are oriented upwards (numbering along the vertical, horizontal edges
%downstream from the inflow boundary. Includes list of boundary nodes and
%boundary edges. Only cells completely inside the domain are considered.

% input:     n_x: half the number of nodes in the horizontal
%            n_z: scalar, number of nodes in the vertical
%            bd_x: x domain boundary 
%            bd_z: upper boundary of the domain. lower bdy is set to zero
%            by default

% Last modified: Nov 20th, 2017


%define topology and boundary nodes (does NOT include boundary edges) 
fout.n_nodes.vert = n_z;
fout.n_nodes.hor = n_x;
fout.n_nodes.tot = fout.n_nodes.hor*fout.n_nodes.vert;
fout.n_edges.hor = (fout.n_nodes.hor -1)*fout.n_nodes.vert;
fout.n_edges.vert = (fout.n_nodes.vert -1)*fout.n_nodes.hor;

index_nodes = (1:fout.n_nodes.tot)';
up = index_nodes; up(fout.n_nodes.hor:fout.n_nodes.hor:fout.n_nodes.tot)=[];
down = index_nodes; down(1:fout.n_nodes.hor:fout.n_nodes.tot)=[];

%list of nodes where Dirichlet or Neumann conditions apply
fout.bdy_nodes.top = (1:fout.n_nodes.hor)';
fout.bdy_nodes.inflow = (1:fout.n_nodes.hor:fout.n_nodes.tot)';
fout.bdy_nodes.outflow = (fout.n_nodes.hor:fout.n_nodes.hor:fout.n_nodes.tot)';
fout.bdy_nodes.bed = (fout.n_nodes.hor*(fout.n_nodes.vert-1)+1:fout.n_nodes.tot)';

%define up and down nodes for each edge 
fout.up_node.hor = index_nodes(up); 
fout.down_node.hor = index_nodes(down);

fout.down_node.vert = reshape(reshape( (1:fout.n_nodes.tot-fout.n_nodes.hor),[ fout.n_nodes.hor, fout.n_nodes.vert-1]).',[fout.n_edges.vert,1]);
fout.up_node.vert = reshape(reshape( (fout.n_nodes.hor+1:fout.n_nodes.tot),[fout.n_nodes.hor,fout.n_nodes.vert-1]).',[fout.n_edges.vert,1]);

fout.bdy_edges_ver.outflow = fout.n_edges.vert -(fout.n_nodes.vert-1) +1:fout.n_edges.vert;
fout.bdy_edges_ver.inflow = 1:(fout.n_nodes.vert-1);
fout.bdy_edges_hor.bed = fout.n_edges.hor -(fout.n_nodes.hor-1) +1:fout.n_edges.hor;
%define vectors of length of ver edges and length of cells in the ver
%direction

%start from position of edges
f_z = @(z) z.^n;
%vertical
z_even = flip(linspace(0,bd_z,n_z+1).');

z_uneven = f_z(z_even);
z_ver_edges = z_uneven;

%horizontal
x_even = linspace(0,bd_x,n_x+1);
x_hor_edges =  x_even;

%compute coordinate of nodes, length of cells and edges
z_nodes = (z_ver_edges(1:end-1)+z_ver_edges(2:end))./2;
ver_cell_length = (z_ver_edges(1:end-1)-z_ver_edges(2:end));
ver_edge_length = z_nodes(1:end-1)-z_nodes(2:end);

x_nodes = (x_hor_edges(1:end-1)+x_hor_edges(2:end))./2;
hor_cell_length = (x_hor_edges(2:end)-x_hor_edges(1:end-1));
hor_edge_length = x_nodes(2:end)-x_nodes(1:end-1);

%construct list of cell lengths, edge lengths, coordinates
fout.Delta_z_edge = reshape(repmat(ver_edge_length, [1,fout.n_nodes.hor]),[fout.n_edges.vert,1]);
fout.Delta_z_cell = reshape(repmat(ver_cell_length, [1,fout.n_nodes.hor]).',[fout.n_nodes.tot,1]);

fout.Delta_x_edge = reshape(repmat(hor_edge_length, [fout.n_nodes.vert,1]).',[fout.n_edges.hor,1]);
fout.Delta_x_cell = reshape(repmat(hor_cell_length, [fout.n_nodes.vert,1]).',[fout.n_nodes.tot,1]);

%define locations of nodes (nodes are indexed starting from the top left of
%the box to the bottom right, horizontally)
xline =  x_nodes;
zline =  z_nodes;
fout.coor_nodes.x = reshape(repmat(xline,[fout.n_nodes.vert 1]).', [ fout.n_nodes.tot,1]);
fout.coor_nodes.z = reshape( repmat([ zline.'],[ fout.n_nodes.hor,1]), [ fout.n_nodes.tot,1]);

%define location of edges, keeping the boundaries of the box separate
%hor edges correspond to u nodes
xline_edge = x_hor_edges(2:end-1);
fout.coor_horedges.x = reshape(repmat(xline_edge,[fout.n_nodes.vert 1]).', [ fout.n_edges.hor,1]);
fout.coor_horedges.z = reshape( repmat([ zline.'],[ (fout.n_nodes.hor -1),1]), [ fout.n_edges.hor,1]);

fout.coord_horedges_inflow.x = ones(fout.n_nodes.vert,1)*x_hor_edges(1);
fout.coord_horedges_inflow.z = zline;

fout.coord_horedges_outflow.x = ones(fout.n_nodes.vert,1)*x_hor_edges(end);
fout.coord_horedges_outflow.z = zline;

%ver edges correspond to w nodes
z_line_edge =  z_ver_edges(2:end-1);
fout.coor_veredges.x = reshape(repmat( xline,[length(z_line_edge) 1]), [ fout.n_edges.vert,1]);
fout.coor_veredges.z = reshape( repmat([ z_line_edge.'],[ fout.n_nodes.hor,1]).', [ fout.n_edges.vert,1]);

fout.coord_veredges_top.x =  xline.';
fout.coord_veredges_top.z = ones(fout.n_nodes.hor,1)*z_ver_edges(1);

fout.coord_veredges_inflow.z = z_line_edge;
fout.index_bed_to_psi_hor_edges = reshape(repmat((1:fout.n_nodes.hor-1).', [ 1 fout.n_nodes.vert]), [fout.n_edges.hor 1]);
end 

function fout  = grid_1D(n_x,bd_x)
%Generates 1D grid in the domain  0<x<xp, with domain boundaries corresponding to cell centres.
% horizontal edges
%downstream from the inflow boundary. Includes list of boundary nodes and
%boundary edges. Only cells completely inside the domain are considered.

% input:     n_x: half the number of nodes in the horizontal
%            n_z: scalar, number of nodes in the vertical
%            bd_x: x domain boundary 
%            bd_z: upper boundary of the domain. lower bdy is set to zero
%            by default

% Last modified: Nov 20th, 2017


%define topology and boundary nodes 
fout.n_nodes.hor = n_x;
fout.n_nodes.tot = fout.n_nodes.hor;
fout.n_edges.hor = (fout.n_nodes.hor);

index_nodes = (1:fout.n_nodes.tot)';
up = circshift(index_nodes,1); %up(fout.n_nodes.hor:fout.n_nodes.hor:fout.n_nodes.tot)=[];
down = index_nodes; %down(1:fout.n_nodes.hor:fout.n_nodes.tot)=[];

%list of nodes where Dirichlet or Neumann conditions apply
fout.bdy_nodes.inflow = 1';
fout.bdy_nodes.outflow = fout.n_nodes.tot';

%define up and down nodes for each edge 
fout.up_node.hor = index_nodes(up); 
fout.down_node.hor = index_nodes(down);

%define vectors of length of ver edges and length of cells in the ver
%direction

%start from position of edges
f_z = @(z) z;

%horizontal
x_even =linspace(bd_x/(2*(n_x)),bd_x+bd_x/(2*(n_x)),n_x+1);
x_uneven = f_z(x_even);
x_hor_edges =  x_uneven;

%compute coordinate of nodes, length of cells and edges
x_nodes = (x_hor_edges(1:end-1)+x_hor_edges(2:end))./2;
hor_cell_length = (x_hor_edges(2:end)-x_hor_edges(1:end-1));
hor_edge_length = x_nodes(1:end)-[0 x_nodes(1:end-1)];

%construct list of cell lengths, edge lengths, coordinates
fout.Delta_x_edge = reshape(repmat(hor_edge_length, [1,1]).',[fout.n_edges.hor,1]);
fout.Delta_x_cell = reshape(repmat(hor_cell_length, [1,1]).',[fout.n_nodes.tot,1]);

%define locations of nodes (nodes are indexed starting from the top left of
%the box to the bottom right, horizontally)
xline =  x_nodes;
fout.coor_nodes.x = reshape(repmat(xline,[1 1]).', [ fout.n_nodes.tot,1]);

%define location of edges, keeping the boundaries of the box separate
%hor edges correspond to u nodes
xline_edge = x_hor_edges(1:end-1);
fout.coor_horedges.x = reshape(repmat(xline_edge,[1 1]).', [ fout.n_edges.hor,1]);

fout.coord_horedges_inflow.x = x_hor_edges(1);
fout.coord_horedges_outflow.x = x_hor_edges(end);

end 


function fout  = grid_T_bl(n_x,n_z,bd_x, bd_z, n)
%Generates rectangular grid for T in the domain (0<z<z_p, 0<x<xp), with T nodes lying on domain boundaries.
%Vertical edges are oriented upwards (numbering along the vertical, horizontal edges
%downstream from the inflow boundary. Includes list of boundary nodes and
%boundary edges. Only cells completely inside the domain are considered.

% input:     n_x: half the number of nodes in the horizontal
%            n_z: scalar, number of nodes in the vertical
%            bd_x: x domain boundary 
%            bd_z: upper boundary of the domain. lower bdy is set to zero
%            by default

% Last modified: Dec 8th, 2017


%define topology and boundary nodes (does NOT include boundary edges) 
fout.n_nodes.vert = n_z;
fout.n_nodes.hor = n_x;
fout.n_nodes.tot = fout.n_nodes.hor*fout.n_nodes.vert;
fout.n_edges.hor = (fout.n_nodes.hor )*fout.n_nodes.vert;
fout.n_edges.vert = (fout.n_nodes.vert -1)*fout.n_nodes.hor;

index_nodes = (1:fout.n_nodes.tot)';
up = index_nodes; %up(fout.n_nodes.hor:fout.n_nodes.hor:fout.n_nodes.tot)=[];
down = index_nodes; %down(1:fout.n_nodes.hor:fout.n_nodes.tot)=[];

%list of nodes where Dirichlet or Neumann conditions apply
fout.bdy_nodes.top = (1:fout.n_nodes.hor)';
fout.bdy_nodes.inflow = (1:fout.n_nodes.hor:fout.n_nodes.tot)';
fout.bdy_nodes.outflow = (fout.n_nodes.hor:fout.n_nodes.hor:fout.n_nodes.tot)';
fout.bdy_nodes.bed = (fout.n_nodes.hor*(fout.n_nodes.vert-1)+1:fout.n_nodes.tot)';
fout.bdy_edges_hor.outflow = fout.n_edges.hor - (fout.n_nodes.vert)+1 :fout.n_edges.hor;
fout.bdy_edges_ver.bed = fout.n_edges.vert - (fout.n_nodes.hor)+1:fout.n_edges.vert;
fout.bdy_edges_hor.bed = fout.n_nodes.vert: fout.n_nodes.vert : fout.n_edges.hor;
%define up and down nodes for each edge 
% fout.up_node.hor = index_nodes(up); 
% fout.down_node.hor = index_nodes(down);
% 
% fout.down_node.vert = reshape(reshape( (1:fout.n_nodes.tot-fout.n_nodes.hor),[ fout.n_nodes.hor, fout.n_nodes.vert-1]).',[fout.n_edges.vert,1]);
% fout.up_node.vert = reshape(reshape( (fout.n_nodes.hor+1:fout.n_nodes.tot),[fout.n_nodes.hor,fout.n_nodes.vert-1]).',[fout.n_edges.vert,1]);

%numbering of edges so that T hor edges are numbered as psi ver edges, and
%T ver edges are numbered as psi hor edges

fout.down_node.vert = (1:fout.n_nodes.tot-fout.n_nodes.hor).';
fout.up_node.vert = (fout.n_nodes.hor+1:fout.n_nodes.tot).' ;

fout.down_node.hor = reshape(reshape(index_nodes(down), [fout.n_nodes.hor,fout.n_nodes.vert]).', [fout.n_edges.hor,1]);
fout.up_node.hor = reshape(reshape(index_nodes(up), [fout.n_nodes.hor,fout.n_nodes.vert]).', [fout.n_edges.hor,1]);
fout.up_node.hor = circshift(fout.down_node.hor, fout.n_nodes.vert);

index_edge_hor = 1:fout.n_edges.hor;
up_edge = index_edge_hor(1:end-fout.n_nodes.vert);
down_edge = index_edge_hor(fout.n_nodes.vert+1:end);

fout.up_edge.hor = reshape(reshape(up_edge,[fout.n_nodes.vert,fout.n_nodes.hor-1]).', [length(up_edge),1]);
fout.down_edge.hor = reshape(reshape(down_edge,[fout.n_nodes.vert,fout.n_nodes.hor-1]).', [length(down_edge),1]);


%define vectors of length of ver edges and length of cells in the ver
%direction

%start from position of psi_edges
f_z = @(z) z.^n;
%vertical
z_even = flip(linspace(0,bd_z,fout.n_nodes.vert+1).');

if bd_z>0
    z_uneven = f_z(z_even);
elseif bd_z<0
    z_uneven = -flip(f_z(abs(z_even)));
end

z_psi_ver_edges = z_uneven;

%horizontal
x_even =linspace(bd_x/(2*(n_x)),bd_x+bd_x/(2*(n_x)),n_x+1);
%x_uneven = f_z(x_even);
x_hor_edges =  x_even;

%compute coordinate of nodes, length of cells and edges
z_psi_nodes = (z_psi_ver_edges(1:end-1)+z_psi_ver_edges(2:end))./2;
psi_ver_cell_length = (z_psi_ver_edges(1:end-1)-z_psi_ver_edges(2:end));
psi_ver_edge_length = z_psi_nodes(1:end-1)-z_psi_nodes(2:end);

% if bd_z >0
%     ver_cell_length = [2*(z_psi_nodes(1)-z_psi_nodes(2));psi_ver_edge_length;2*(z_psi_nodes(end))];
% elseif bd_z <0
%     ver_cell_length = [2*(-z_psi_nodes(1));psi_ver_edge_length;2*(z_psi_nodes(end-1)- z_psi_nodes(end))];
% end
% 
% ver_edge_length = psi_ver_cell_length;

x_nodes = (x_hor_edges(1:end-1)+x_hor_edges(2:end))./2;
hor_cell_length = (x_hor_edges(2:end)-x_hor_edges(1:end-1));
hor_edge_length = x_nodes(1:end)-[0 x_nodes(1:end-1)];

%define locations of nodes (nodes are indexed starting from the top left of
%the box to the bottom right, horizontally)
xline =  x_nodes;
zline =  z_psi_nodes;
fout.coor_nodes.x = reshape(repmat(xline,[fout.n_nodes.vert 1]).', [ fout.n_nodes.tot,1]);
fout.coor_nodes.z = reshape( repmat([ zline.'],[ fout.n_nodes.hor,1]), [ fout.n_nodes.tot,1]);


%construct list of cell lengths, edge lengths, coordinates  CHECK!!!!
fout.Delta_z_edge = reshape(repmat(psi_ver_edge_length, [1,fout.n_nodes.hor]).',[fout.n_edges.vert,1]);
fout.Delta_z_cell = reshape(repmat(psi_ver_cell_length, [1,fout.n_nodes.hor]).',[fout.n_nodes.tot,1]);

fout.Delta_x_edge = reshape(repmat(hor_edge_length, [fout.n_nodes.vert,1]),[fout.n_edges.hor,1]);
fout.Delta_x_cell = reshape(repmat(hor_cell_length, [fout.n_nodes.vert,1]).',[fout.n_nodes.tot,1]);


%define location of edges, keeping the boundaries of the box separate
%hor edges correspond to u nodes
xline_edge = x_hor_edges(1:end-1);
fout.coor_horedges.x = reshape(repmat(xline_edge,[fout.n_nodes.vert 1]), [ fout.n_edges.hor,1]);
fout.coor_horedges.z = reshape( repmat([ zline],[ 1 fout.n_nodes.hor]), [ fout.n_edges.hor,1]);

fout.coord_horedges_inflow.x = ones(fout.n_nodes.vert,1)*x_hor_edges(1);
fout.coord_horedges_inflow.z = zline;

fout.coord_horedges_outflow.x = ones(fout.n_nodes.vert,1)*x_hor_edges(end);
fout.coord_horedges_outflow.z = zline;

fout.coor_veredges.x = reshape(repmat(xline,[fout.n_nodes.vert-1 1]).', [ fout.n_edges.vert,1]);
fout.coor_veredges.z = reshape( repmat([ z_psi_ver_edges(2:end-1)],[ 1, fout.n_nodes.hor]).', [ fout.n_edges.vert,1]);

fout.coord_veredges_top.x =  xline.';
fout.coord_veredges_top.z = ones(fout.n_nodes.hor,1)*z_psi_ver_edges(1);

%mapping for boundary layer sliding velocity
fout.index_bed_to_T_hor_edges = reshape(repmat((1:fout.n_nodes.hor), [ fout.n_nodes.vert 1]), [fout.n_edges.hor 1]);
fout.index_bed_to_T_ver_nodes = reshape(repmat((1:fout.n_nodes.hor), [fout.n_nodes.vert-1,1]).', [fout.n_edges.vert 1]);

fout.index_1D_to_2D = ones(fout.n_nodes.vert,1);
end 

