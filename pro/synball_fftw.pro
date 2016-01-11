pro synball_fftw,an,cl,ball,almn=almn,almr=almr,beam=beam,lmax=lmax_in,grcn=grcn,nmax_sim=nmax_sim,$
                 nshell=nshell,nside=nside,napo=napo,fr=fr,cn=cn,rseed=rseed,iseed=iseed,rfr=rfr,ran=ran,rcn=rcn

  usage='usage: synball,an,cl,ball,almn=almn,almr=almr,beam=beam,lmax=lmax,nmax=nmax,nshell=nshell,nside=nside'
  if n_params() lt 3 then message,usage
  if not keyword_set(nside) then message,usage
  if not keyword_set(nshell) then message,usage
  if not keyword_set(lmax_in) then message,usage
  if not keyword_set(grcn) then grcn=1 ;; generate g(r) from power law C_n

  lmax=lmax_in
  nmax=nshell/2
  if keyword_set(nmax_sim) then nmax=nmax_sim

  create_alm,cl,alm,lmax=lmax,12312
  almn=dcomplexarr(nshell,lmax+1,lmax+1)

;; Nf = double(nshell)
;; r = dindgen(nshell)/Nf
;; fr = [randomu(123,nmax/2), 0*indgen(nmax/2)]*cos(2d0*!pi*r) 


  if grcn eq 1 then begin
     print, 'g(r) will be generated from a power law C_n with no dipole and up to nmax=',nshell/2-1
     r = (dindgen(nshell/2+1)+1)
     g=10*(r)^(-1./2.)/sqrt(2d0)
     g[0]=0d0
     g[nshell/2]=0d0

     ran = complex(g,g) ;; an such that c_n=|a_n|^2=n^-alpha
     cn_sim = 100*r^(-1)
     cn_sim[0]=0d0
     ;;print, 'n^-2 - abs(ran)^2',cn-abs(ran)^2

     ;;convert from a_nmax to IDL a_n format
     a=findgen(nshell/2+1)
     b=-reverse(findgen(nshell/2-1)+1)
     c=[a,b]                    ;[0,1,...,N/2,-N/2-1,..-1]


     an_sim = [ran[abs(c[0:nshell/2])],conj(ran[abs(c[(nshell/2+1):nshell-1])])]
     fr_sim = double(fft(an_sim,1))

     ;;to check my an order is similar to idl one. If so an-an2=0, which it is
     an=fft(fr_sim,-1)
     fr=double(fft(an,1))
     cn = abs(an[0:nshell/2])^2

     rfr=double(fft(an,1))
     ran=fft(rfr,-1)
     rcn = abs(ran[0:nshell/2])^2
     ;;print, 'n^-2 - abs(ran)^2',cn-abs(an2[0:nshell/2])^2

  endif else begin
     print, 'g(r) will be generated as (r-0.5)^n a power law C_n with no dipole and up to nmax=',nmax
     TT=1
     Npts=nshell
     t=dindgen(Npts)/(Npts-1)*TT
     ;;t[0]=1d0 ;;0 padding
     ;;t[Npts-1]=1d0 ;;0 padding
     
     fr = 0d0*t
     for ii=1,nmax do begin
        fr = fr+ gaussian(t,[1.,0.5,1./double(ii)]);;(1d0-1.5*t)^double(ii)
     endfor
     ;;fr[0]=0d0
     ;;fr[nshell-1]=0d0
     an=fft(fr,-1)              ; complex fourier coefficients
     cn = abs(an[0:nshell/2])^2

     rfr=double(fft(an,1))
     ran=fft(rfr,-1)
     rcn = abs(ran[0:nshell/2])^2
  endelse
  
;;Gaussian realization
;;cn = dindgen(nmax)^2d0 
;;ran = sqrt(cn)*dcomplex(randomn(rseed,nmax),randomn(rseed,nmax))

  fr = abs(fft(an,1))

  for ll=0,lmax do begin
     for mm=0,lmax do begin
        almn[0:nshell-1,ll,mm]=alm[0,ll,mm]*an[0:nshell-1]
     endfor
  endfor

;;almn2ball_sincos,almn,ball,nside=nside,nshell=nshell,lmax=lmax,nmax=nmax,nproc=nproc,almr=almr,rmin=0d0,rmax=1d0

  almn2ball,almn,ball,nside=nside,nshell=nshell,lmax=lmax,nmax=nmax,nproc=nproc,almr=almr

;;appodize the edge in the radial direction

  if keyword_set(napo) then begin
     ball[*,0]=0.
     ball[*,nshell-1]=0.
     
     if napo gt 1 then begin
        factor = abs(dindgen(napo)/double(napo-1) - 1d0)*!pi/2d0
        
        for i=0,napo-1 do begin
           
           apo_factor = cos(factor[i])
           
           print, i,apo_factor
           
           ball[*,i]=apo_factor*ball[*,napo-1]
           ball[*,nmax-1-i]=apo_factor*ball[*,nmax-1-napo]
        endfor
     endif
  endif


end
