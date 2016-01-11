pro simball,ball,beam=beam,lmax=lmax_in,nmax=nmax_in,nshell=nshell,nside=nside

usage='usage: synball,an,cl,ball,almn=almn,almr=almr,beam=beam,lmax=lmax,nmax=nmax,nshell=nshell,nside=nside'
if n_params() lt 1 then message,usage
if not keyword_set(nside) then message,usage
if not keyword_set(nshell) then message,usage

;; nmax=n_elements(an)-1
;; lmax=n_elements(cl)-1

;; if keyword_set(lmax_in) then lmax=lmax_in
;; if keyword_set(nmax_in) then nmax=nmax_in

;if (lmax gt lmax) then message,'given lmax is larger than n_elements(cl)'
;if (nmax gt nmax_an) then message,'given nmax is larger than n_elements(an)'

npix=nside2npix(nside)
lpix=lindgen(npix)

pix2vec_ring, nside, lpix, vec

ball=dblarr(npix,nshell+1)

r=dindgen(nshell+1)/double(nshell)
for ii=0,nshell-1 do begin
   ball[*,ii]=(1d0-r[ii])*r[ii]*exp(-3d0*r[ii]^2)*(-vec[*,0]^2+0.001*vec[*,1]+r[ii]*vec[*,2]+2)  ;
   ;;ball[*,ii]=cos(2d0*!pi*r[ii])
   ;;+cos(4*2d0*!pi*r[ii])+cos(6*2d0*!pi*r[ii])+cos(8*2d0*!pi*r[ii])+sin(4*2d0*!pi*r[ii])

   ;;ball[*,ii]=(r[ii]-1d0)^2d0;;cos(2d0*!pi*r[ii])*sin(4*2d0*!pi*r[ii])
endfor



end
