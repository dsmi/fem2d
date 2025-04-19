function v = integ_uu( r, t, ni, nj )
% function v = integ_uu( r, t, ni, nj )
%
%  Evaluates integral of the product of the basis functions
%  ni and nj for each triangle:
%   v = \integ N_i N_j dx dy
%

% Vertices of the triangle, N-by-2
r1 = r(t(:,1),:);
r2 = r(t(:,2),:);
r3 = r(t(:,3),:);

% Edges of the triangle, N-by-2
r12 = r2 - r1;
r23 = r3 - r2;
r31 = r1 - r3;

% Triangle areas
A = abs( r12(:,1).*r23(:,2) - r12(:,2).*r23(:,1) ) ./ 2;

% Integrals
v = A * ( 1 + (ni == nj) ) ./ 12;
