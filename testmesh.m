
dcirc = @(p,xc,yc,r) sqrt((p(:,1)-xc).^2+(p(:,2)-yc).^2)-r;
drectangle = @(p,x1,x2,y1,y2) -min(min(min(-y1+p(:,2),y2-p(:,2)),-x1+p(:,1)),x2-p(:,1));
fd = @(p) max( drectangle(p,-0.01,0.01,-0.01,0.01), max( -dcirc(p,0.004,0,0.001), -dcirc(p,-0.004,0,0.001) )  );
fh = @(p) 0.0005 + 0.3*min(dcirc(p,0.004,0,0.001), dcirc(p,-0.004,0,0.001));
[p,t] = distmesh( fd, fh, 0.0005, [-1,-1;1,1], [-0.01,-0.01;-0.01,0.01;0.01,-0.01;0.01,0.01] );
patch( 'vertices', p, 'faces', t, 'facecolor', [.9, .9, .9] )
