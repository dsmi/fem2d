function v = integ_gg( r, t, ni, nj )
% function v = integ_gg( r, t, ni, nj )
%
%  Evaluates integral of product of gradients of the basis functions
%  Ni and nj for each triangle:
%   v = \integ \delta N_i \delta N_j dx dy
%

% Vertices of the triangle, N-by-2
r1 = r(t(:,1),:);
r2 = r(t(:,2),:);
r3 = r(t(:,3),:);

% Triangle jacobians, 2-by-2-by-N, last index is triangle
J = cat( 1, permute( r1-r3, [ 3 2 1 ] ), ...
            permute( r2-r3, [ 3 2 1 ] ) );

% det(J), 1-by-1-by-N
detj = J(1,1,:).*J(2,2,:) - J(1,2,:).*J(2,1,:);

% inv(J)
invj = repmat( 1.0./detj, 2, 2 ) .* [  J(2,2,:) , -J(1,2,:) ; ...
                                      -J(2,1,:) ,  J(1,1,:) ];
    
% Gradients of the basis functions with respect to
% u and v, 3-by-2 matrix, first index is the basis function
gradn = [ 1 0 ; 0 1 ; -1 -1 ];

% Number of triangles
nt = size( t, 1 );

% Gradients of i'th and j'th basis functions repeated for each triangle
gradni = repmat( gradn(ni,:), [ 2 1 nt ] );
gradnj = repmat( gradn(nj,:), [ 2 1 nt ] );

v = permute( detj*(1/2).*dot( dot( invj, gradni, 2 ), ...
                              dot( invj, gradnj, 2 ), 1 ), [ 3 1 2 ] );
