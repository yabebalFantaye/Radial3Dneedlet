psi_jqk,outbeta,nside=nside,lmax=lmax_in,nmax=nmax_in,nj=nj,nv=nv,jint=jint,bb=bb,pow=pow


;;fix radial point
ir=5
map_shell=reform(outbeta[*,ir,*])

;;for the north pole pixel, get all pixels along theta [0,pi] 
ip=500

nring=4l*long(nside)
thvec=findgen(nring)*!pi/float(nring-1)
phvec=0*thvec

ang2pix_ring,nside, thvec,phvec,listpix

;pix2vec_ring, nside, ip, vec
;query_disc,nside,vec,20d0,listpix,/deg


map_thvec=map_shell[listpix,*] 


end
