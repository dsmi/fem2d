
% microstrip parameters
w = 1e-4; % width
t = 2e-5; % thickness
h = 1e-4; % height above ground
eps1 = eps0 * 4.2; % dielectric

% meshed area dimensions
mw = max( w, h ) * 6;
mh =   ( h + t ) * 4;

% signed distance function for a rectangle
drectangle = @(p,x1,x2,y1,y2) -min(min(min(-y1+p(:,2),y2-p(:,2)),-x1+p(:,1)),x2-p(:,1));

% Outer boundary
outb = @(p) drectangle(p,-mw/2,mw/2,0,mh);

% Microstrip boundary
mstripb = @(p) drectangle(p,-w/2,w/2,h,h+t);

% Boundary signed distance function.
fd = @(p) max( outb(p), -mstripb(p) );

% Signed distance function for (relative) mesh distribution.
% We want it closer to the microstrip.
fh = @(p) 0.0001+1*mstripb(p);

h0 = min( w, t ) / 3;

nbv = round((mw-w)/2/h0);
bv = linspace( -mw/2, -w/2, nbv );

% Constraint edges not supported on Octave -- add multiple fixed
% vertices along the dielectric interface.
fv = [-mw/2,0;-mw/2,mh;mw/2,0;mw/2,mh;...
     -w/2,h;-w/2,h+t;w/2,h;w/2,h+t;[bv',bv'*0+h];[-bv',bv'*0+h]];

[r,tri] = distmesh( fd, fh, h0, [-mw/2,0;mw/2,mh], fv, [] );

nt = size( tri, 1 );

% Dielectric permeability
eps_tri = eps0 * ones( nt, 1 );

% triangle centers.
tc = ( r(tri(:,1),:) + r(tri(:,2),:) + r(tri(:,3),:) ) / 3;

% Dielectric between microstrip and ground         
eps_tri( tc(:,2) < h ) = eps1;

% total number of vertices, known and unknown
nv = size( r, 1 );

% Vertices where we know the potential
%% vground = verts_by( r, @(p) drectangle(p,-mw/2,mw/2,0,0), 1.0e-8 );
vground = verts_by( r, outb, 1.0e-8 ); % the outer boundary
vmstrip = verts_by( r, mstripb, 1.0e-8 );

% The known potentials
uk = zeros( nv, 1 );
uk( vground ) = 0.0; % 0 at the outer boundary
uk( vmstrip ) = 1.0; % 1 at the microstrip
                       
%% patch( 'vertices', r, 'faces', tri, 'facevertexcdata', uk,
%%        'FaceColor','interp', 'cdatamapping', 'scaled' )

edges = make_edges( tri );

nv = size( r, 1 );

% build full stifness matrix
S = sparse( nv, nv );
for ni=1:3
    for nj=1:3
        gg = eps_tri.*integ_gg( r, tri, ni, nj );
        S = S + sparse( tri(:,ni), tri(:,nj), gg, nv, nv );
    end
end

% Indices of the boundary vertices (both outer and microstrip)
vbound = [ vground ; vmstrip ];

% Indices of the vertices with u unknown
vunk = find( !ismember( 1:nv, sort( vbound ) ) );

nvunk = length(vunk);

% Matrix to get rows/columns for the unknown vertices from S
K = sparse( vunk, transpose( 1:nvunk ), ones( nvunk, 1 ), nv, nvunk );

% This gives us weights of the shape functions for the unknown vertices
x = (K'*S*K)\(-K'*S*uk);

% Potentials of all vertices!
u = uk + K*x;

% all boundary edges
bedges = edges( edges_connected_between( edges, vbound ), : );

% Matrix of the edge basis functions products for all boundary edges
% connected to the corresponding vertex.
U = sparse( nv, nv );
for ni=1:2
    for nj=1:2
        uu = integ_edge_uu( r, bedges, ni, nj );
        U = U + sparse( bedges(:,ni), bedges(:,nj), uu, nv, nv );
    end
end

nvbound = length(vbound);

% Matrix to get rows/columns for the known/boundary vertices from S
T = sparse( vbound, transpose( 1:nvbound ), ones( nvbound, 1 ), nv, nvbound );

% This gives us weights of the flux shape functions for the boundary vertices
x = (T'*U*T)\(-T'*S*u);

% Flux weights for all vertices, zero for non-boundary ones
vf = T*x;

% edge lengths
l = sqrt( sum( ( r( edges(:,1),:) - r( edges(:,2),:) ).^2, 2 ) );

% fluxes for all edges. only the boundary ones make sense.
ef = l./2.*( vf( edges(:,1) ) + vf( edges(:,2) ) );

% Edges forming the outer boundary and boundary of the microstrip.
eground = edges_connected_between( edges, vground );
emstrip = edges_connected_between( edges, vmstrip );

flux_ground = sum(ef(eground));
flux_mstrip = sum(ef(emstrip));

C = -flux_mstrip

gradu = trigrad( r, tri, u );

colormap('jet')
patch( 'vertices', r, 'faces', tri, 'facevertexcdata', u,
       'FaceColor','interp', 'cdatamapping', 'scaled' )
%% hold on
%% quiver( tc(:,1), tc(:,2), gradu(:,1), gradu(:,2), 'm' );
%% hold off

