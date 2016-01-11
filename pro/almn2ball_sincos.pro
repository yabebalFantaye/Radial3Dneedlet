

;;;IDL extension to HEALPix [TESTVERSION]: Written by Frode K. Hansen, Nov. 2006
;;yabebal fantaye 12 Jan 2013

pro almn2ball_sincos,almn,map,nside=nside_in,beam=beam_in,fwhm=fwhm_in,nopixwin=nopixwin,pol=pol,$
                     silent=silent,d1=d1,d2=d2,der1=der1_in,der2=der2_in,large=large_in,rmin=rmin,rmax=rmax,$
                     lmax=lmax_in,nmax=nmax_in,nshell=nshell_in,nproc=nproc_in,almr=almr,$
                     sinoutmap=sinoutmap,cosoutmap=cosoutmap,monopole=monopole

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

  nproc = 4
  if keyword_set(nproc_in) then nproc=nproc_in


  routine = 'almn2ball_sincos'
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
                                ;print, nmax*(1l+polar*2l)*(lmax+1l)^2
     almn=dcomplexarr(1l+polar*2l,lmax+1,lmax+1,nmax+1,2)

     openr,unit,fname,/f77
     readu,unit, almn
     close,unit
  endif 


  lmax=n_elements(reform(almn(0,*,0,0,0)))-1
  nmax=n_elements(reform(almn(0,0,0,*,0)))-1


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



  ;;check for nside consistency of the output map
  if keyword_set(nside_in) then nside=LONG(nside_in)
  if keyword_set(nshell_in) then nshell=LONG(nshell_in)
  if keyword_set(nshell) then nrpix=nshell+1

  if not (KEYWORD_SET(nside_in) and KEYWORD_SET(nshell_in)) then begin

     if size(map,/type) eq 7 then begin
        print,'nside and nrpix (odd integer) must be passed when passing filename to map.'
        return
     endif

     smap=size(map)
     smap=LONG(smap[0])

     if (polar eq 0) and (smap ne 2) then begin
        print,'MAP must be of the form map=fltarr(npix,nrpix,1+pol*2)'
        return
     endif

     if ((polar eq 1) and (smap ne 3)) then begin
        print,'MAP must be of the form map=fltarr(npix,nrpix,1+pol*2)'
        return
     endif
     ;; Check polarized map is TQU
     if (polar eq 1) then begin
        if (n_elements(map(0,0,*)) ne 3) then begin
           print,'MAP must be of the form map=fltarr(npix,nrpix,1+pol*2)'
           return
        endif
     endif

     nside=sqrt(DOUBLE(n_elements(map[*,0]))/12d0)
     nrpix=n_elements(map[0,*])
     nshell=nrpix-1
  endif


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
  writeu,unit,LONG(nrpix)
  writeu,unit,double(rmin),double(rmax)
  close,unit

  for ii=0,1 do begin
     itype=2+ii
     wunf, almn[*,*,*,*,ii],dir+'almn_stype'+strn(itype)+'.unf'
  endfor

  der1=0l
  der2=0l
  if (KEYWORD_SET(d1)) then der1=1l
  if (KEYWORD_SET(d2)) then der2=1l
  openw,unit,dir+'der1.unf',/f77
  writeu,unit,der1
  close,unit
  openw,unit,dir+'der2.unf',/f77
  writeu,unit,der2
  close,unit

  if ((ldir ne '') and (ldir ne 'large/')) then SPAWN,'cp '+dir0+'idl_almn2ball_sincos '+direxe

  str='cd '+direxe+'; time mpirun -np '+strtrim(string(nproc),2)+' -hostfile $HOME/mpihosts --prefix $MPI_HOME  ./idl_almn2ball_sincos'
  print, '************************* spawning **********************'
  print, '** ', str, ' **'
  print, '*********************************************************'
  SPAWN,str


  for ii=0,1 do begin
     itype=2+ii
     openr,unit,dir+'ballmap_stype'+strn(itype)+'.unf',/f77
     
     npix = lonarr(1)
     nrpix=npix
     maxcount=npix
     
     readu,unit,npix,nrpix,prec,maxcount
     npix = npix[0] & nrpix=nrpix[0] & maxcount=maxcount[0]
     map=dblarr(npix,nrpix)
     readu,unit,map
     if ii eq 0 then sinoutmap=map+monopole
     if ii eq 1 then cosoutmap=map+monopole
     close,unit
  endfor
  map = sinoutmap + cosoutmap


;if keyword_set(almr) then begin
  almr=dcomplexarr(nrpix, 1+2*polar,lmax+1,lmax+1)
  alm=dcomplexarr(1+2*polar,lmax+1,lmax+1)
  itype=2



  for ii=0,nshell-1 do begin
     openr,unit,dir+'almr_stype'+strn(itype)+'.unf_'+strn(ii),/f77
     readu,unit,alm
     close,unit
     almr[ii,*,*,*] = alm
  endfor

;  endif

end





