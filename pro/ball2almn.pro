

;;;IDL extension to HEALPix [TESTVERSION]: Written by Frode K. Hansen, Nov. 2006
;;yabebal fantaye 12 Jan 2013
pro ball2almn,map,almn,lmax=lmax_in,nmax=nmax_in,cl=cl_in,$
              mask=mask_in,pol=pol,all=all,w8=w8,large=large_in,$
              nproc=nproc_in,almr=almr,outalmn=outalmn

  ;dir='$HEALPIX/src/idl/almball/'
  dir='./'
  large_in= 1;'large/'
  spawn,'mkdir -p large'

  dir0=dir
  ldir=''
  if (keyword_set(large_in)) then begin
     if (size(large_in,/tn) eq 'INT') then begin
        if (large_in eq 1) then ldir='large/'
     endif
     if (ldir eq '') then ldir=trim(string(large_in),2)+'/'
     if (ldir ne 'large/') then dir=dir+ldir
  endif
  direxe=dir
  if (ldir eq 'large/') then dir=dir+ldir

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


  routine = 'ball2almn'
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
  if not (KEYWORD_SET(lmax_in)) then lmax=60l

  if (lmax lt 1) then begin
     print,'lmax MUST BE GREATER THAN 0!'
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
     map=dblarr(npix,nshell)
     readu,unit,map
     close,unit
     map = map
  endif

  npix=n_elements(map[*,0])
  nshell=n_elements(map[0,*])
  nmax=nshell
  
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


  openw,unit,dir+'params_ball2almn.unf',/f77
  writeu,unit,LONG((ldir eq 'large/'))
  writeu,unit,LONG(nside)
  writeu,unit,LONG(polar)
  writeu,unit,LONG(lmax)
  writeu,unit,LONG(nmax)
  writeu,unit,LONG(nshell)
  close,unit

  openw,unit,dir+'ballmap.unf',/f77
  writeu,unit,double(map)
  close,unit


  if ((ldir ne '') and (ldir ne 'large/')) then SPAWN,'cp '+dir0+'idl_ball2almn '+direxe

  str='cd '+direxe+'; time mpirun -np '+strtrim(string(nproc),2)+' -hostfile $HOME/mpihosts --prefix $MPI_HOME  ./idl_ball2almn'
  print, '************************* spawning **********************'
  print, '** ', str, ' **'
  print, '*********************************************************'
  SPAWN,str ;,stderr
  ;print, stderr

  ;;    if (keyword_set(nsim_in)) then begin
  ;; endif else begin
  ;;    SPAWN,'cd '+direxe+' ; ./idl_ball2almn'
  ;; endelse


  if size(almn,/type) eq 7 then begin
     ;;spawn, 'ls -l '+dir+'almn.unf '
     print, 'spawning: ','cp '+dir+'almn.unf '+ almn 
     spawn, 'cp '+dir+'almn.unf '+ almn 
  endif else begin
     almn=dcomplexarr(nshell,lmax+1,lmax+1)
     openr,unit,dir+'almn.unf',/f77
     readu,unit,almn
     close,unit
     ;;cl_in=almn2cl(almn,all=all)
  endelse

;if keyword_set(outalmn) then begin
     outalmn=dcomplexarr(nshell,lmax+1,lmax+1)
     openr,unit,dir+'almn.unf',/f77
     readu,unit,outalmn
     close,unit
;endif

;if keyword_set(almr) then begin
     almr=dcomplexarr(nshell,lmax+1,lmax+1)
     openr,unit,dir+'almr.unf',/f77
     readu,unit,almr
     close,unit
;;(0:nshell,1:1+2*polar,0:lmax,0:lmax)
;endif

end


