function vin = verts_by( r, fd, tol )
% vin = verts_by( r, fd, tol )
%
% Indices of the vertices for which the specified distance function is zero
% within tol.
%

vin = find( abs( feval( fd, r ) ) < tol );
