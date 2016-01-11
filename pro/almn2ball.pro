

;;;IDL extension to HEALPix [TESTVERSION]: Written by Frode K. Hansen, Nov. 2006
;;yabebal fantaye 12 Jan 2013

pro almn2ball,almn,map,nside=nside_in,beam=beam_in,fwhm=fwhm_in,nopixwin=nopixwin,pol=pol,$
              silent=silent,d1=d1,d2=d2,der1=der1_in,der2=der2_in,large=large_in,$
              lmax=lmax_in,nmax=nmax_in,nshell=nshell_in,nproc=nproc_in,almr=almr,monopole=monopole

                                ;dir='$HEALPIX/src/idl/almball/'
  dir='./'
  large_in=1
  spawn,'mkdir -p large'

  if not keyword_set(monopole) then monopole=0.

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

  nproc = 1
  if keyword_set(nproc_in) then nproc=nproc_in


  routine = 'almn2ball'
  uroutine = strupcase(routine)
  print, '--------- '+uroutine+' ----------'

  if (n_params() ne 2) then begin
     PRINT, 'Wrong number of arguments in ' + uroutine
     print,'Syntax : '
     print,uroutine+', almn, map,'
     print,'              [NSIDE=, BEAM=, FWHM=, NOPIXWIN=, POL= , SILENT] '
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

  ;; load almn if a file name is passed
  if size(almn,/type) eq 7 then begin

     if not KEYWORD_SET(lmax_in) or not KEYWORD_SET(nmax_in) then begin
        print, 'lmax and nmax must be set when passing filename to almn.'
        return
     endif
     
     fname = almn
     print, 'lmax, nmax,polar = ',lmax,nmax,polar

     almn=dcomplexarr(nmax,lmax+1,lmax+1)

     openr,unit,fname,/f77
     readu,unit, almn
     close,unit

  endif else begin

     ;;check for lmax, nmax consitency given alm
     if (n_elements(almn) lt 2) then begin
        print,'almn MUST BE of the form almn=complexarr(nmax,lmax,lmax) OR LMAX must be specified!'
        return
     endif

     salmn=size(almn)
     salmn=LONG(salmn[0])
     if ((salmn lt 3) and (salmn gt 4)) then begin
        print,'almn MUST BE of the form almn=complexarr(nmax+1, 1+2*pol,lmax+1,lmax+1) OR NMAX, LMAX must be specified'
        return
     endif

     if (salmn eq 3) then begin
        nmax=n_elements(reform(almn(*,0,0)))
        lmax=n_elements(reform(almn(0,*,0)))-1
     endif
     print, 'almn2beta.pro lmax, nmax',lmax, nmax
     ;; lmax_in=lmax
     ;; nmax_in=nmax

     if (salmn eq 4) then begin
        polar=n_elements(reform(almn(*,0,0,0)))
        if polar eq 1 then polar=0l
        if polar gt 1 then polar=1l
     endif
     
  endelse

  if (lmax lt 1) then begin
     print,'lmax MUST BE GREATER THAN 0!'
     return
  endif

  if (nmax lt 1) then begin
     print,'nmax MUST BE GREATER THAN 0!'
     return
  endif

  nshell=nmax
  if keyword_set(nshell_in) then begin
     if nshell ne nshell_in then begin
        print, 'nshell_in and nmax from almn are not the same'
        return
     endif
  endif

  if not (KEYWORD_SET(nside_in)) then begin
     print,'nside must be passed when passing filename to map.'
     return
  endif
  nside=nside_in



  cns=alog(DOUBLE(nside))/alog(2d0)
  if ((abs(cns-NINT(cns)) gt 1d-6) or (nside lt 1)) then begin
     print,'nside MUST BE of the form 2^N'
     return
  endif
  if (lmax gt 3*nside) then begin
     if (NOT (KEYWORD_SET(silent))) then print,'WARNING:lmax>3*nside will use lmax=3*nside'
     lmax=3*nside
  endif

  fwhm0=0.
  if (KEYWORD_SET(fwhm_in)) then fwhm0=double(fwhm_in)
  beam0=dblarr(lmax+1,1l+polar*2l)
  beam0(*,*)=1d0
  if (KEYWORD_SET(beam_in)) then begin
     beam=beam_in
     sbeam=SIZE(beam0)
     sbeam=LONG(sbeam[0])
     if ((sbeam lt 1) or (sbeam gt 2)) then begin
        print,'BEAM must be of the form beam=dblarr(lmax+1,1+pol)'
        return
     endif
     if (sbeam eq 1) then beam0(*,0)=beam(0:lmax)
     if (sbeam eq 2) then beam0(*,0)=beam(0:lmax,0)
     if (polar eq 1l) then begin
        if (sbeam eq 1) then beam0(*,1)=beam(0:lmax)
        if (sbeam eq 1) then beam0(*,2)=beam(0:lmax)
        if (sbeam eq 2) then beam0(*,1)=beam(0:lmax,1)
        if (sbeam eq 2) then beam0(*,2)=beam(0:lmax,1)
     endif
  endif
  if (NOT (KEYWORD_SET(nopixwin))) then begin
     wl=pixwin(nside)
     beam0(0:lmax,0)=beam0(0:lmax,0)*wl(0:lmax,0)
     if (polar eq 1l) then begin
        beam0(0:lmax,1)=beam0(0:lmax,1)*wl(0:lmax,0)
        beam0(0:lmax,2)=beam0(0:lmax,2)*wl(0:lmax,0)
     endif
     if (NOT (KEYWORD_SET(silent))) then print,'WARNING! Using pixel window function'
  endif

  openw,unit,dir+'beam.unf',/f77
  writeu,unit,beam0
  close,unit
  openw,unit,dir+'fwhm.unf',/f77
  writeu,unit,fwhm0
  close,unit

  ;;write parameter file 
  openw,unit,dir+'params_almn2ball.unf',/f77
  writeu,unit,LONG((ldir eq 'large/'))
  writeu,unit,LONG(nside)
  writeu,unit,LONG(polar)
  writeu,unit,LONG(lmax)
  writeu,unit,LONG(nmax)
  writeu,unit,LONG(nshell)
  close,unit

  wunf, almn,dir+'almn.unf'

  ;; der1=0l
  ;; der2=0l
  ;; if (KEYWORD_SET(d1)) then der1=1l
  ;; if (KEYWORD_SET(d2)) then der2=1l
  ;; openw,unit,dir+'der1.unf',/f77
  ;; writeu,unit,der1
  ;; close,unit
  ;; openw,unit,dir+'der2.unf',/f77
  ;; writeu,unit,der2
  ;; close,unit

  if ((ldir ne '') and (ldir ne 'large/')) then SPAWN,'cp '+dir0+'idl_almn2ball '+direxe

  str='cd '+direxe+'; time mpirun -np '+strtrim(string(nproc),2)+' -hostfile $HOME/mpihosts --prefix $MPI_HOME  ./idl_almn2ball'
  print, '************************* spawning **********************'
  print, '** ', str, ' **'
  print, '*********************************************************'
  SPAWN,str


  if size(map, /type) eq 7 then begin
     spawn, 'cp '+dir+'ballmap.unf '+map
  endif 

  openr,unit,dir+'ballmap.unf',/f77
  
  npix = lonarr(1)
  nshell=npix
  maxcount=npix
  
  readu,unit,npix,nshell,prec,maxcount
  npix = npix[0] & nshell=nshell[0] & maxcount=maxcount[0]
  map=dblarr(npix,nshell)
  readu,unit,map
  close,unit

  map = map+monopole  

  ;; if (der1 eq 1l) then begin
  ;;    der1_in=dblarr(nside^2*12l,2l*(1l+polar*2l))
  ;;    openr,unit,dir+'map_der1.unf',/f77
  ;;    readu,unit,der1_in
  ;;    close,unit
  ;; endif
  ;; if (der2 eq 1l) then begin
  ;;    der2_in=dblarr(nside^2*12l,3l*(1l+polar*2l))
  ;;    openr,unit,dir+'map_der2.unf',/f77
  ;;    readu,unit,der2_in
  ;;    close,unit
  ;; endif


;if keyword_set(almr) then begin
     almr=dcomplexarr(nshell,lmax+1,lmax+1)
     alm=dcomplexarr(1,lmax+1,lmax+1)
     for ii=0,nshell-1 do begin

        openr,unit,dir+'almr.unf_'+strn(ii),/f77
        readu,unit,alm
        close,unit
        almr[ii,*,*] = alm[0,*,*]
     endfor

;  endif

end





