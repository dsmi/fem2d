function v = integ_edge_uu( r, edges, ni, nj )
% function v = integ_edge_uu( r, edges, ni, nj )
%
%  Evaluates integral of shape times test functions associated with the
%  specified edges.
%  ni and nj for each edge
%   v = \integ N_i N_j dl
%

% edge endpoints
r1 = r( edges(:,1),:);
r2 = r( edges(:,2),:);

l = sqrt( sum( ( r1 - r2 ).^2, 2 ) );

k = 1/6 + ( ni == nj ) * 1/6;

v = k*l;
