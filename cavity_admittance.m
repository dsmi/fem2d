%
% Rectangular cavity with two ports on the short sides.
%

% The geometry
inch2meter = 2.54e-2;
mil2meter = 1.0e-3*inch2meter;
d = (10.0 + 1.35 + 10.0)*mil2meter; % plane-to-plane separation
w = 1e-2;  % cavity width  (x size)
l = 2e-2;  % cavity length (y size)
x1 = -w/6;        % first port
y1 = -l/2 + l/6;
x2 =  w/4;        % second port
y2 =  l/2 - l/6;
rp =  w/80;

% Dielectric params
lt = 0.02;
er0 = 4.3;
fr = 1e9;

% metal conductivity
sigma = 5.8e7;

% signed distance function for a rectangle
drectangle = @(p,x1,x2,y1,y2) ...
             -min(min(min(-y1+p(:,2),y2-p(:,2)),-x1+p(:,1)),x2-p(:,1));

% signed distance function for a circle
dcirc = @(p,xc,yc,r) sqrt((p(:,1)-xc).^2+(p(:,2)-yc).^2)-r;

% Outer boundary
foutb = @(p) drectangle(p,-w/2,w/2,-l/2,l/2);

% ports
fp1 = @(p) dcirc( p, x1, y1, rp );
fp2 = @(p) dcirc( p, x2, y2, rp );

% Boundary signed distance function.
fd = @(p) max( foutb(p), max( -fp1(p), -fp2(p) ) );

% Signed distance function for (relative) mesh distribution.
% We want it denser closer to the ports.
fh = @(p) 0.0002 + 0.3*min(fp1(p), fp2(p));

h0 = rp/2;

% Fixed vertices
fv = [-w/2,-l/2;-w/2,l/2;w/2,-l/2;w/2,l/2 ];

[ r, tri ] = distmesh( fd, fh, h0, [-w/2,-l/2;w/2,l/2], fv, [] );

nt = size( tri, 1 );

%% patch( 'vertices', r, 'faces', tri, 'facecolor', [.9, .9, .9] )

% total number of vertices, known and unknown
nv = size( r, 1 );

% Port vertices
vp1 = verts_by( r, fp1, 1.0e-8 );
vp2 = verts_by( r, fp2, 1.0e-8 );

% The known potentials. Each column applies excitation to
% the corresponding port.
uk = zeros( nv, 2 );
uk( vp1, 1 ) = 1.0;
uk( vp2, 2 ) = 1.0;
                       
%% patch( 'vertices', r, 'faces', tri, 'facevertexcdata', uk(:,1),
%%        'FaceColor','interp', 'cdatamapping', 'scaled' )

edges = make_edges( tri );

nv = size( r, 1 );

% Indices of the port vertices -- nonzero flux and fixed potential
vports = [ vp1 ; vp2 ];

% Indices of the vertices with u unknown
vunk = find( !ismember( 1:nv, sort( vports ) ) );

nvunk = length(vunk);

% Matrix to get rows/columns for the unknown vertices from S
K = sparse( vunk, transpose( 1:nvunk ), ones( nvunk, 1 ), nv, nvunk );
% boundary where we want to know the flux
fedges = edges( edges_connected_between( edges, vports ), : );

nvports = length(vports);

% Matrix to get rows/columns for the known/port vertices from S
T = sparse( vports, transpose( 1:nvports ), ones( nvports, 1 ), nv, nvports );

% Edges forming the boundary of each port
ep1 = edges_connected_between( edges, vp1 );
ep2 = edges_connected_between( edges, vp2 );

nports = 2;
nedges = size( edges, 1 )

% Matrix to sum the fluxes of the port edges to get the port currents
P = sparse( [ ep1*0+1 ; ep2*0+2 ], [ ep1 ; ep2 ], 1.0, nports, nedges );

% angular frequencies
% freqs = 1e9*2*pi;
freqs = linspace(1e5, 1e11, 201)*2*pi;

Yf = [] ; % Simulated admittance for all frequency points

for freq = freqs,
    freq

    % Calculate cavity admittance, impedance and wavenumber
    er = debye(er0, lt, fr, freq/(2*pi));
    Yplane = j*freq*eps0*er/d;
    Zs = 2*sqrt(j*freq*mu0/sigma); % surface impedance
    Zplane = Zs + j*freq*mu0*d;
    k = sqrt(-Yplane*Zplane);

    % build full stifness matrix
    S = sparse( nv, nv );
    for ni=1:3
        for nj=1:3
            gg = integ_gg( r, tri, ni, nj );
            S = S + sparse( tri(:,ni), tri(:,nj), gg, nv, nv );
            uu = integ_uu( r, tri, ni, nj );
            S = S - k*k*sparse( tri(:,ni), tri(:,nj), uu, nv, nv );
        end
    end

    A = (K'*S*K);
    x = (-K'*S*uk);

    % This gives us weights of the shape functions for the unknown vertices
    x = A\x;

    % Potentials of all vertices!
    u = uk + K*x;


    % Matrix of the edge basis functions products for all flux edges
    % connected to the corresponding vertex.
    U = sparse( nv, nv );
    for ni=1:2
        for nj=1:2
            uu = integ_edge_uu( r, fedges, ni, nj );
            U = U + sparse( fedges(:,ni), fedges(:,nj), uu, nv, nv );
        end
    end

    % This gives us weights of the flux shape functions for the boundary vertices
    x = (T'*U*T)\(-T'*S*u);

    % Flux weights for all vertices, zero for non-boundary ones
    vf = T*x;

    % edge lengths
    l = sqrt( sum( ( r( edges(:,1),:) - r( edges(:,2),:) ).^2, 2 ) );

    % edge fluxes, nedges-by-nports. only the boundary fluxes make sense.
    ef = repmat(l, 1, nports)/2.*( vf( edges(:,1), : ) + vf( edges(:,2), : ) );

    % Port currents, which is also the admittance matrix since we applied
    % unit voltages to the ports.
    Y = -P*ef/Zplane;

    Yf = cat(3, Yf, Y);
    
end

tswrite( 'cavity_admittance.y2p', freqs/(2*pi), Yf, 'Y', 50 );

%% colormap('jet')
%% patch( 'vertices', r, 'faces', tri, 'facevertexcdata', u( :, 1 ),
%%        'FaceColor','interp', 'cdatamapping', 'scaled' )


