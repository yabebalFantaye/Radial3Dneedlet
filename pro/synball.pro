pro synball,an,cl,ball,almn=almn,almr=almr,beam=beam,lmax=lmax_in,nmax=nmax_in,nshell=nshell,nside=nside

usage='usage: synball,an,cl,ball,almn=almn,almr=almr,beam=beam,lmax=lmax,nmax=nmax,nshell=nshell,nside=nside'
if n_params() lt 3 then message,usage
if not keyword_set(nside) then message,usage
if not keyword_set(nshell) then message,usage

nmax=n_elements(an)-1
lmax=n_elements(cl)-1

if keyword_set(lmax_in) then lmax=lmax_in
if keyword_set(nmax_in) then nmax=nmax_in

;if (lmax gt lmax) then message,'given lmax is larger than n_elements(cl)'
;if (nmax gt nmax_an) then message,'given nmax is larger than n_elements(an)'


create_alm,cl,alm,lmax=lmax,12312

almn=dcomplexarr(1,lmax+1,lmax+1,nmax+1)

for nn=0,nmax do begin
   for ll=0,lmax do begin
      almn[0,ll,0:lmax,nn]=alm[0,ll,0:lmax]*an[nn]
   endfor
endfor

;;almn2ball_sincos,almn,ball,nside=nside,nshell=nshell,lmax=lmax,nmax=nmax,nproc=nproc,almr=almr,rmin=0d0,rmax=1d0

almn2ball,almn,ball,nside=nside,nshell=nshell,lmax=lmax,nmax=nmax,nproc=nproc,almr=almr,rmin=0d0,rmax=1d0


end
