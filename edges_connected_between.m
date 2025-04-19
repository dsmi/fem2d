function cei = edges_connected_between(edges, verts)
% cei = edges_connected_between(edges, verts)
%
% Indices of the edges connected between the specified vertices.
%

sortedv = sort(verts);
    
cei = find( ismember( edges(:,1), sortedv ) & ismember( edges(:,2), sortedv ) );
