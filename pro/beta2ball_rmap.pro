

;;;IDL extension to HEALPix [TESTVERSION]: Written by Frode K. Hansen, Nov. 2006
;;yabebal fantaye 12 Jan 2013
pro beta2ball_rmap,betamap,map,lmax=lmax_in,nmax=nmax_in,cl=cl_in,nj=nj,nv=nv,jint=jint,bb=bb,$
              mask=mask_in,pol=pol,all=all,w8=w8,feedback=feedback,nside=nside,nrpix=nrpix,nout=nout,$
              nproc=nproc_in,rmin=rmin,rmax=rmax,nrsample=nrsample_in,almr=almr,outmap=outmap,almn=almn,$
              outbeta=outbeta,norun=norun

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

  nproc = 4
  if keyword_set(nproc_in) then nproc=nproc_in

  routine = 'beta2ball'
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


  if (KEYWORD_SET(lmax_in)) then lmax=LONG(lmax_in)
  if (KEYWORD_SET(nmax_in)) then nmax=LONG(nmax_in)

  ;;if lmax and nmax not passed, extract it from almn array or issue
  ;;an error
  if  not KEYWORD_SET(lmax_in) and not KEYWORD_SET(nmax_in) then begin
     print, 'lmax and nmax must be set when passing filename to almn.'
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

  if not keyword_set(rmin) then rmin=0d0
  if not keyword_set(rmax) then rmax=1d0

  if not keyword_set(nrpix) then nrpix=16  
  nshell=nrpix-1

  if not keyword_set(nside) then nside=16
  npix=nside2npix(nside)

  nside_out=nside
  if keyword_set(nout) then nside_out=nout



  ;; nside=sqrt(DOUBLE(npix)/12d0)
  ;; cns=alog(DOUBLE(nside))/alog(2d0)
  ;; if ((abs(cns-NINT(cns)) gt 1d-6) or (nside lt 1)) then begin
  ;;    print,'nside MUST BE of the form 2^N'
  ;;    return
  ;; endif

  ;; if (KEYWORD_SET(mask_in)) then begin
  ;;    mask=fltarr(npix,1+2*polar)
  ;;    mask[*,*]=1.
  ;;    mask=mask_in

  ;;    if (n_elements(mask(*,0)) ne npix) then begin
  ;;       print,'map and mask MUST HAVE THE SAME NSIDE!'
  ;;       return
  ;;    endif
  ;;    smask=size(mask)
  ;;    smask=LONG(smask[0])
  ;;    if ((smask lt 1) or (smask gt 2)) then begin
  ;;       print,'MASK must be of the form mask=fltarr(npix,1+pol*2)'
  ;;       return
  ;;    endif
  ;;    if (polar eq 0) then mask=reform(mask(*,0))
  ;;    if ((polar eq 1) and (smask eq 1)) then begin
  ;;       mask=fltarr(npix,3)
  ;;       mask(*,0)=mask_in
  ;;       mask(*,1)=mask_in
  ;;       mask(*,2)=mask_in
  ;;    endif
  ;;    if ((polar eq 1) and (smask eq 2)) then begin
  ;;       mask=fltarr(npix,3)
  ;;       mask(*,0)=mask_in(*,0)
  ;;       mask(*,1)=mask_in(*,1)
  ;;       mask(*,2)=mask_in(*,1)
  ;;    endif
     
  ;;    for ii=0,nrpix-1 do begin
  ;;       map[*,ii]=map[*,ii]*mask[*,0]
  ;;    endfor

  ;; endif

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


if not keyword_set(jint) then jint=1
if not keyword_set(bb) then bb=1.1
if not keyword_set(nj) then nj=5
if not keyword_set(feedback) then feedback=3
if not keyword_set(nv) then nv=nshell
;;ceil(double(bb)^(jint+nj-1))



fnamepar = dir+'params_beta2ball.unf'
  openw,unit,fnamepar,/f77
  writeu,unit,double(bb)
  writeu,unit,LONG(jint)
  writeu,unit,LONG(nj)
  writeu,unit,LONG(nv)
  writeu,unit,LONG(nside)
  writeu,unit,LONG(nside_out)
  writeu,unit,LONG(polar)
  writeu,unit,LONG(lmax)
  writeu,unit,LONG(nmax)
  writeu,unit,LONG(nrpix)
  writeu,unit,double(rmin),double(rmax)
  writeu,unit,LONG(feedback)
  close,unit


  str='cd '+direxe+'; mpirun -np '+strtrim(string(nproc),2)+' -hostfile $HOME/mpihosts --prefix $MPI_HOME  ./idl_beta2ball_rmap '+fnamepar
  print, '************************* spawning **********************'
  print, '** ', str, ' **'
  print, '*********************************************************'
  if not keyword_set(norun) then SPAWN,str ;,stderr


  if size(map, /type) eq 7 then begin
     print, 'beta2ball.pro: copying  '+dir+'map_ball_beta.unf '+map
     spawn, 'cp '+dir+'map_ball_beta.unf '+map
  endif else map=outmap


;;  if keyword_set(outmap) then begin
     print,'beta2ball.pro: reading: ',dir+'map_ball_beta.unf'
     openr,unit,dir+'map_ball_beta.unf',/f77
     
     npix = lonarr(1)
     nrpix=npix
     maxcount=npix
     
     readu,unit,npix,nrpix,prec,maxcount
     npix = npix[0] & nrpix=nrpix[0] & maxcount=maxcount[0]
     outmap=fltarr(npix,nrpix)
     readu,unit,outmap
     close,unit
;;  endif


;;  outmap=fltarr(npix,nrpix)
 ;; runf, outmap,dir+'map_ball_beta.unf'
  ;;outmap=outmap-mean(outmap) ;;remove monopole


;;--------------------
;; outbeta = fltarr(npix,nv+1,nj)
;; single_betamap = fltarr(npix,1,nv+1)
;; for jj=0,nj-1 do begin
;;    openr,unit,'temp_t3/beta_ball_j'+strn(jj)+'.unf',/f77
;;    readu,unit,single_betamap 
;;    close,unit
;;    outbeta[0:npix-1,0:nv,jj] = single_betamap[0:npix-1,0,0:nv]
;; endfor


  
  
if keyword_set(almr) then begin
     almr=complexarr(nrpix, 1+2*polar,lmax+1,lmax+1)
     openr,unit,dir+'arlm_from_beta.unf',/f77
     readu,unit,almr
     close,unit
;;(0:nshell,1:1+2*polar,0:lmax,0:lmax)
endif

;; if keyword_set(almn) then begin
;;      almn=complexarr(1+2*polar,lmax+1,lmax+1,nmax+1)
;;      openr,unit,dir+'almn_from_beta.unf',/f77
;;      readu,unit,almn
;;      close,unit
;; ;;(0:nshell,1:1+2*polar,0:lmax,0:lmax)
;; endif


end


