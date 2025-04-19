function vin = verts_in_box( r, box )
% vin = verts_in_box( r, box )
%
% Indices of the vertices which are in the specified box.
%   box = [ xmin,ymin ; xmax,ymax ]
%

vin = find( r(:,1) >= box(1,1) & r(:,1) < box(2,1) ...
             & r(:,2) >= box(1,2) & r(:,2) < box(2,2) );
