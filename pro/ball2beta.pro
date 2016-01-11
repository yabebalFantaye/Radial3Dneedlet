

;;;IDL extension to HEALPix [TESTVERSION]: Written by Frode K. Hansen, Nov. 2006
;;yabebal fantaye 12 Jan 2013
pro ball2beta,map,betamap,lmax=lmax_in,nmax=nmax_in,cl=cl_in,nj=nj,jint=jint,bb=bb,pow=pow,$
              mask=mask_in,pol=pol,all=all,w8=w8,feedback=feedback,gln=gln,$
              nproc=nproc_in,almr=almr,outbeta=outbeta,$
              almn=almn,norun=norun,float=float,itype=itype,norm=norm

                                ;dir='$HEALPIX/src/idl/almball/'
  dir='./temp/'
  direxe='./'
  spawn,'mkdir -p '+dir


;direxe='$HEALPIX/bin/'
  unit=99
  status = FSTAT(unit)
  if (status.open ne 0) then begin
     WHILE(status.open ne 0) do begin
        unit=unit-1
        status = FSTAT(unit)
     ENDWHILE
  endif

  nproc = 1
  if keyword_set(nproc_in) then nproc=nproc_in

  routine = 'ball2beta'
  uroutine = strupcase(routine)
  print, '--------- '+uroutine+' ----------'


  if (n_params() ne 2) then begin
     PRINT, 'Wrong number of arguments in ' + uroutine
     print,'Syntax : '
     print,uroutine+', map, alm, '
     print,'              [LMAX=, NMAX=, CL=, MASK= ,POL= ,ALL=,W8= ] '
;    print,' Type '+uroutine+', /help '
;    print,'   for an extended help'
     return
  endif

  polar=0l
  if (KEYWORD_SET(pol)) then polar=1l

  if polar ne 0l then begin
     print, 'ball2almn can not yet handle polarized data.'
     return
  endif


  if size(map,/type) eq 7 then begin
     ;;read pixelized ball map
     print, 'reading file: '+map
     openr,unit,map,/f77

     npix = lonarr(1)
     nshell=lonarr(1)
     maxcount=lonarr(1)
     
     readu,unit,npix,nshell,prec,maxcount
     npix = npix[0] & nshell=nshell[0] & maxcount=maxcount[0]

     if not keyword_set(float) then map=dblarr(npix,nshell)
     if keyword_set(float) then  map=fltarr(npix,nshell)
     
     readu,unit,map
     close,unit
     map = double(map)
  endif 

  if keyword_set(norm) then map=norm*map
  ;;ckeck for map array validity  
  if (polar eq 0) then map=reform(map)

  smap=size(map)
  smap=LONG(smap[0])
  
  npix=n_elements(map[*,0])
  nshell=n_elements(map[0,*])
  nmax=nshell

  if (KEYWORD_SET(lmax_in)) then lmax=LONG(lmax_in)
  if (KEYWORD_SET(nmax_in)) then nmax=LONG(nmax_in)

  ;;if lmax and nmax not passed, extract it from almn array or issue
  ;;an error
  if  not KEYWORD_SET(lmax_in)  then begin
     print, 'lmax must be set when passing filename to almn.'
     return
  endif

  if (lmax lt 1) then begin
     print,'lmax MUST BE GREATER THAN 0!'
     return
  endif

  if (nmax lt 1) then begin
     print,'nmax MUST BE GREATER THAN 0!'
     return
  endif

  ;;map[*,1:nshell-1] = map[*,1:nshell-1]-mean(map[*,1:nshell-1])

  nside=sqrt(DOUBLE(npix)/12d0)
  cns=alog(DOUBLE(nside))/alog(2d0)
  if ((abs(cns-NINT(cns)) gt 1d-6) or (nside lt 1)) then begin
     print,'nside MUST BE of the form 2^N'
     return
  endif

  if (KEYWORD_SET(mask_in)) then begin
     mask=fltarr(npix,1+2*polar)
     mask[*,*]=1.
     mask=mask_in

     if (n_elements(mask(*,0)) ne npix) then begin
        print,'map and mask MUST HAVE THE SAME NSIDE!'
        return
     endif
     smask=size(mask)
     smask=LONG(smask[0])
     if ((smask lt 1) or (smask gt 2)) then begin
        print,'MASK must be of the form mask=fltarr(npix,1+pol*2)'
        return
     endif
     if (polar eq 0) then mask=reform(mask(*,0))
     if ((polar eq 1) and (smask eq 1)) then begin
        mask=fltarr(npix,3)
        mask(*,0)=mask_in
        mask(*,1)=mask_in
        mask(*,2)=mask_in
     endif
     if ((polar eq 1) and (smask eq 2)) then begin
        mask=fltarr(npix,3)
        mask(*,0)=mask_in(*,0)
        mask(*,1)=mask_in(*,1)
        mask(*,2)=mask_in(*,1)
     endif
     
     for ii=0,nshell-1 do begin
        map[*,ii]=map[*,ii]*mask[*,0]
     endfor

  endif

  if (KEYWORD_SET(w8)) then begin
     wfile='weight_ring_n'
     wnumfile=trim('00000'+trim(string(LONG(nside)),2),2)
     wsize=strlen(wnumfile)
     wfile='$HEALPIX/data/'+wfile+strmid(wnumfile,wsize-5,wsize)+'.fits'
     read_fits_map,wfile,wl
     wl=wl+1d0
  endif else begin
     wl=dblarr(2*nside,3)
     wl(*,*)=1d0
  endelse
  wunf,wl,dir+'w8.unf'

  glnpow=1
  if keyword_set(pow) then glnpow=pow
  if not keyword_set(jint) then jint=0
  if not keyword_set(bb) then bb=2d0
  if not keyword_set(nj) then nj=6
  if not keyword_set(feedback) then feedback=3


  fnamepar = dir+'params_ball2beta.unf'
  openw,unit,fnamepar,/f77
  writeu,unit,LONG(glnpow)
  writeu,unit,double(bb)
  writeu,unit,LONG(jint)
  writeu,unit,LONG(nj)
  writeu,unit,LONG(nside)
  writeu,unit,LONG(polar)
  writeu,unit,LONG(nmax)
  writeu,unit,LONG(lmax)
  writeu,unit,LONG(nshell)
  writeu,unit,LONG(feedback)
  close,unit


  openw,unit,dir+'ballmap.unf',/f77
  writeu,unit,double(map)
  close,unit

  f90code='./idl_ball2beta_par '
  print, 'parameter file to '+f90code


  str='cd '+direxe+'; time mpirun -np '+strtrim(string(nproc),2)+' -hostfile $HOME/mpihosts --prefix $MPI_HOME  '+f90code+' '+fnamepar
  print, '************************* spawning **********************'
  print, '** ', str, ' **'
  print, '*********************************************************'
  if not keyword_set(norun) then  SPAWN,str ;,stderr
                                ;print, stderr

  ;;    if (keyword_set(nsim_in)) then begin
  ;; endif else begin
  ;;    SPAWN,'cd '+direxe+' ; ./idl_ball2almn'
  ;; endelse



  print, 'npix, nshell: ',npix,nshell

  outbeta = dblarr(npix,nshell,nj) ;;2-for sin and cos basis
  single_betamap = dblarr(npix,nshell)


  for jj=0,nj-1 do begin

     filename=dir+'beta_ball_j'+strn(jj)+'.unf'
     if glnpow gt 1 then filename=dir+'beta_ball_j'+strn(jj)+'_glnpow'+strn(glnpow)+'.unf'
     
     openr,unit,filename,/f77
     readu,unit,single_betamap 
     close,unit
     print, 'jj, sum(beta(0:npix-1,0:nshell) ',jj,total(single_betamap)
     outbeta[0:npix-1,0:nshell-1,jj] = single_betamap[0:npix-1,0:nshell-1]
  endfor

;;endfor

  if size(betamap,/type) eq 7 then begin
     wunf,outbeta,betamap
  endif 


  almr=dcomplexarr(nshell,lmax+1,lmax+1)
  openr,unit,dir+'almr_from_ball.unf',/f77
  readu,unit,almr
  close,unit

;;if keyword_set(gln) then begin

  gln = dblarr(lmax+1,nshell,nj)
  fname=dir+'gln.unf'
  if (glnpow gt 1) then fname=dir+'gln2.unf'
  openr,unit,fname,/f77
  readu,unit,gln
  close,unit


  almn=dcomplexarr(nshell,lmax+1,lmax+1)
  openr,unit,dir+'almn_from_ball.unf',/f77
  readu,unit,almn
  close,unit
  

;;endif

end


