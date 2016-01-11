

;;;IDL extension to HEALPix [TESTVERSION]: Written by Frode K. Hansen, Nov. 2006
;;yabebal fantaye 12 Jan 2013
pro psi_jqk,outbeta,nside=nside,lmax=lmax_in,nmax=nmax_in,nj=nj,nv=nv,jint=jint,bb=bb,pow=pow,$
              feedback=feedback,gln=gln,$
              nproc=nproc_in,rmin=rmin,rmax=rmax,norun=norun,foutbeta=foutbeta

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

  routine = 'psi_jqk'
  uroutine = strupcase(routine)
  print, '--------- '+uroutine+' ----------'
  if (n_params() ne 1) then begin
     PRINT, 'Wrong number of arguments in ' + uroutine
     print,'Syntax : '
     print,uroutine+', outbeta, '
     print,'              [LMAX=, NMAX=, CL=, MASK= ,POL= ,ALL=,W8= ] '
;    print,' Type '+uroutine+', /help '
;    print,'   for an extended help'
     return
  endif


  ;;r sampling per shel
  nrsample=0l
  polar=0l

  lmax=60
  nmax=20
  if (KEYWORD_SET(lmax_in)) then lmax=LONG(lmax_in)
  if (KEYWORD_SET(nmax_in)) then nmax=LONG(nmax_in)


  glnpow=1
  if keyword_set(pow) then glnpow=pow
  if not keyword_set(jint) then jint=0
  if not keyword_set(bb) then bb=2d0
  if not keyword_set(nj) then nj=6
  if not keyword_set(feedback) then feedback=3
  
  if not keyword_set(rmin) then rmin=0d0
  if not keyword_set(rmax) then rmax=1d0
  if not keyword_set(nshell) then nshell=64
  if not keyword_set(nside) then nside=16
  if not keyword_set(nv) then nv=nshell 
  nrpix=nshell+1


  ;;map[*,1:nshell-1] = map[*,1:nshell-1]-mean(map[*,1:nshell-1])







;;ceil(double(bb)^(jint+nj-1))

fnamepar = dir+'params_psi_jqk.unf'
  openw,unit,fnamepar,/f77
  writeu,unit,LONG(glnpow)
  writeu,unit,double(bb)
  writeu,unit,LONG(jint)
  writeu,unit,LONG(nj)
  writeu,unit,LONG(nv)
  writeu,unit,LONG(nside)
  writeu,unit,LONG(polar)
  writeu,unit,LONG(lmax)
  writeu,unit,LONG(nmax)
  writeu,unit,LONG(nrpix)
  writeu,unit,double(rmin),double(rmax)
  writeu,unit,LONG(feedback)
  close,unit

  f90code='./idl_psi_jqk '
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


npix=nside2npix(nside)

print, 'npix, nv: ',npix,nv

outbeta = dblarr(npix,nv+1,nj)
single_betamap = dblarr(npix,1,nv+1)

for jj=0,nj-1 do begin
   filename=dir+'beta_ball_j'+strn(jj)+'.unf'
   if glnpow gt 1 then filename=dir+'beta_ball_j'+strn(jj)+'_glnpow'+strn(glnpow)+'.unf'

   openr,unit,filename,/f77
   readu,unit,single_betamap 
   close,unit
   outbeta[0:npix-1,0:nv,jj] = single_betamap[0:npix-1,0,0:nv]
endfor

if keyword_set(foutbeta) then begin
   wunf,outbeta,foutbeta
endif 

end


