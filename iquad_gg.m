function v = iquad_gg( r, ni, nj )
% function v = iquad_gg( r, ni, nj )
%
%  Evaluates the same integral as integ_qq for one triangle using a fine
%  quadrature -- to be used by the tests.
%   r is 3-by-2 tree vertices of a triangle
%  

% We are integrating a constant -- we don't need the quadrature here. 
%   qN = 2;
%   [ qr qw ] = simplexquad(qN, r)

% Basis functions in triangle, using local coordinates where side 3-1 lies
% on u, side 3-2 lies on v and vertex 3 is at zero:
%  N1 = u
%  N2 = v
%  N3 = 1 - u - v
% Conversion is
%  x = N1*x1 + N2*x2 + N3*x3
%  y = N1*y1 + N2*y2 + N3*y3

% Global to local jacobian
% J = [ dx/du dy/du ;
%       dx/dv dy/dv ]
% J = [ x(1) - x(3) , y(1) - y(3) ;
%       x(2) - x(3) , y(2) - y(3) ];
J = [ r(1,:) - r(3,:) ;...
      r(2,:) - r(3,:) ];

% Gradients of the basis functions with respect to
% u and v, 2-by-3 matrix, second index is the basis function
gradn = transpose( [ 1 0 ; 0 1 ; -1 -1 ] );

% Both gradients are pulled outside of the integral, det( J ) * (1/2) gives
% the integral over the triangle.
v = dot( J \ gradn(:,ni), J \ gradn(:,nj) ) * det( J ) * (1/2);
