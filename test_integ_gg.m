
fd = @(p) -min(min(min(1+p(:,2),1-p(:,2)),1+p(:,1)),1-p(:,1));
fh = @(p) ones(size(p,1),1);
[ r, tri ] = ...
distmesh( fd, fh, 0.5, [-1,-1;1,1], [-1,-1;-1,1;1,-1;1,1] );

% patch( 'vertices', r, 'faces', tri, 'facecolor', [.9, .9, .9] )

for ni=1:3
    for nj=1:3
        v  = integ_gg( r, tri, ni, nj );
        vq = v*0;
        for t=1:size(tri,1);
            rt = r( tri(t,:), : );
            vq(t) = iquad_gg( rt, ni, nj );
        end
        assert( norm(v-vq)<1.0e-10, "test_integ_gg failed" )
    end
end
