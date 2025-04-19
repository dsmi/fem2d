function edges = make_edges(tri)
% edges = make_edges(tri)
%
% Given a triangular mesh, makes a set of unique edges of size nedges-by-2
%

% Number of triangles
ntri = size(tri,1);

% All edges - ones different by direction only not merged yet.
all_edge_tris = repmat(1:ntri,3,1);
all_edge_tris = all_edge_tris(:);
all_free_vert_loc = repmat([1; 2; 3],ntri,1);
edgev = [ 2 3 ; 3 1 ; 1 2 ];
idx = sub2ind(size(tri), repmat(all_edge_tris,1,2), repmat(edgev,ntri,1) );
all_edges = tri(idx);

% Unify the edges direction - vertex with the lesser index comes first.
to_flip = find(all_edges(:,1) > all_edges(:,2));
v1 = all_edges(to_flip,1);
all_edges(to_flip,1) = all_edges(to_flip,2);
all_edges(to_flip,2) = v1;

edges = unique( all_edges, 'rows' );
