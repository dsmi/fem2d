function gr = trigrad( r, tri, u )
% function gr = trigrad( r, tri, u )
%
%  Given a scalar function defined at the vertices of the unstructured
%  mesh made of triangles, calculate gradient in each trianble.
%

% Vertices of the triangle
v = cell(3,1);
for i = 1:3    
    v{i} = r(tri(:,i),:);
end

% Edges, edge i is opposite to vertex i
e = cell(3,1);
for i = 1:3    
    e{i} = v{ rem(i+1,3)+1 } - v{ rem(i,3)+1 };
end

% triangle areas, signed
A = e{1}(:,1).*e{2}(:,2) - e{2}(:,1).*e{1}(:,2);

% Edges rotated 90 degrees. Inside for CCW triangle.
j = cell(3,1);
for i = 1:3    
    j{i} = [ -e{i}(:,2) e{i}(:,1) ];
end

% We can now calculate the gradients
gr = zeros( size(tri, 1), 2 );
for i = 1:3
    gr = gr + repmat( 1./(2*A).*u(tri(:,i)), 1, 2 ).*j{i};
end    
