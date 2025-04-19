function v = iquad_edge_uu( r1, r2, ni, nj )
% function v = iquad_edge_uu( r1, r2, ni, nj )
%
%  Same as integ_edge_uu but using a quadrature, used for testing.
%
    
l = sqrt( sum( ( r2 - r1 ).^2, 2 ) );

% 2-point quadrature as we are integrating the second order polynomial
x = w = ones( 1, 2 );
x(1) = -1/sqrt(3);
x(2) = -x(1);

v  = r2 - r1;
t  = ( x + 1 ) ./ 2;
v1 = ( 2 - ni ) * t + ( ni - 1 ) * ( 1 - t );
v2 = ( 2 - nj ) * t + ( nj - 1 ) * ( 1 - t );

v = sum( v1 .* v2 .* w ) * ( l/2 );

