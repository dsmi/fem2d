function v = iquad_uu( r, ni, nj )
% function v = iquad_gg( r, ni, nj )
%
%  Evaluates the same integral as integ_uu for one triangle using a fine
%  quadrature -- to be used by the tests.
%   r is 3-by-2 tree vertices of a triangle
%  

qN = 50;
[ qr qw ] = simplexquad( qN, 2 );

% Basis functions in triangle, using local coordinates where side 3-1 lies
% on u, side 3-2 lies on v and vertex 3 is at zero:
%   N1 = u
%   N2 = v
%   N3 = 1 - u - v
% Conversion is
%    x = N1*x1 + N2*x2 + N3*x3
%    y = N1*y1 + N2*y2 + N3*y3
% Vertices of the triangle in this new coordunate system are:
%    p1=(1,0), p2=(0,1), p3=(0,0).
n = zeros( size( qr, 1 ), 3 );
n(:,1) = qr(:,1);
n(:,2) = qr(:,2);
n(:,3) = 1 - n(:,1) - n(:,2);

% Global to local jacobian
% J = [ dx/du dy/du ;
%       dx/dv dy/dv ]
% J = [ x(1) - x(3) , y(1) - y(3) ;
%       x(2) - x(3) , y(2) - y(3) ];
J = [ r(1,:) - r(3,:) ;...
      r(2,:) - r(3,:) ];

v = det( J ) * dot( n(:,ni).*n(:,nj), qw, 1 );
