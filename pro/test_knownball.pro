
cl_file='/mn/owl1/d3/yabebalf/NG-27Oct2011/output/Dolney2004/model/Cl_model_bin5.dat'
readcol,cl_file, ll,clmm,clmm_err


maxcount_mocbox=1
bb=2

nrint=2l^7
nside=64
nshell=2*nrint ;;number of radial shells 

nmax=25 ;;maximum radial multipole
lmax= 65 ;;angular multipole
nrsample=1 ;;r sampling factor for ball2almn

;;input angular power spectra and radial coefficients
clmm=1e5*clmm
anmm=[0,1d0/(indgen(nmax)+1)]


j0=-1
nj = ceil(alog(sqrt(1.*(lmax*(lmax+1.)+(1.*nmax)^2)))/alog(1.*bb)+1)
print, 'nj = ',nj
nproc=1

fkey = '_lmax'+strn(lmax)+'_nmax'+strn(nmax)+'_nr'+strn(nshell)

figdir='../figures_clmm'+fkey+'/'
spawn, 'mkdir -p '+figdir

fball = '../output/ball_clmm'+fkey+'.unf'
fbeta = '../output/beta_clmm'+fkey+'.unf'
fbeta2 = '../output/beta_gln2_clmm'+fkey+'.unf'
falmn = '../output/almn_clmm'+fkey+'.unf'


lmax_sim=60
nmax_sim=20

npix=nside2npix(nside)
fball_orig=dblarr(npix,nshell+1)
r=dindgen(nshell+1)/double(nshell)


;;get an example cl and an to construct input almn
cl_sim=dblarr(lmax+1)
cl_sim[2:lmax_sim]=clmm[1:lmax_sim-1]

an_sin =  anmm[0:nmax_sim] 
an_cos =  max(an_sin)*randomn(123,nmax_sim+1)
an_sim = an_sin+an_cos

;;create angular map
create_alm,cl_sim,alm_sim,lmax=lmax,12312
alm2map,alm_sim,y20,/nopix,nside=nside 


;;simulate band limited ball
synball_sincos,an_sin,an_cos,cl_sim,ball_sim,almn=almn_sim,lmax=lmax_sim,nmax=nmax_sim,nshell=nshell,nside=nside,almr=almr_sim


argpi=!pi
.r
snrvec=dblarr(nshell+1)
mono=snrvec
for ii=0,nshell-1 do begin

   snr_last=0d0
   snr=0d0
   for jj=0,nmax_sim do begin
      ntilda=jj
      snr_last = snr_last + an_sim[jj]*sinc_x(argpi*ntilda*r[nshell])*argpi*ntilda 
      snr = snr + an_sim[jj]*sinc_x(argpi*ntilda*r[ii])*argpi*ntilda   
   endfor

   snrvec[ii]=snr  ;(snr-snr_last)
   fball_orig[*,ii]=y20[*,0]*(snr-snr_last)
   mono[ii] = mean(fball_orig[*,ii])
endfor
end

;;==========================================
;;         Radial line reconstruction
;;==========================================
;;check radial part only
ball2an,snrvec,fan,nmax=nmax,nshell=nshell
an2ball,fan,fout,nshell=nshell,nmax=nmax

tek_color &  plot, indgen(nmax+1),an_sim-fan,xra=[1,nmax_sim],xtitle='n',ytitle=textoidl('a_n') 
tek_color &  plot, indgen(nmax+1),an_sim,xtitle='n',ytitle=textoidl('a_n'),xra=[1,nmax_sim] & oplot, indgen(nmax+1),fan, color=2
tek_color &  plot, r,snrvec,xtitle='r',ytitle='f(r)' & oplot, r, fout, color=2
tek_color &  plot, r,snrvec-fout,xtitle='r',ytitle=textoidl('\Delta f(r)')
tek_color &  plot, r, (snrvec-fout),yra=[-1e-3,1e-3],xtitle='r',ytitle='\Delta f(r)'

planck_plot, indgen(nmax+1),fan-an_sim,xra=[1,nmax_sim],yra=[-1,10]*1e-7,ax=['xttl','yttl'],ay=['$n$','$\Delta a_n$'],fname=figdir+'an_diff_orig_recon.ps',/psend
planck_plot, indgen(nmax+1),an_sim,xra=[1,nmax_sim],ax=['xttl','yttl'],ay=['$n$','$a_n$'],yover=fan,fname=figdir+'an_orig_recon.ps',/psend
planck_plot, r,snrvec,ax=['xttl','yttl'],ay=['$r$','$f(r)$'],yover=fout,fname=figdir+'fr_orig_recon.ps',/psend
planck_plot, r,snrvec-fout,ax=['xttl','yttl'],ay=['$r$','$\Delta f(r)$'],fname=figdir+'fr_diff_orig_recon.ps',/psend


;;check angular part only

;;==========================================
;;         Harmonic space  reconstruction
;;==========================================
;;checl radial and angular part without needlet
ball2almn,ball_sim,falmn,lmax=lmax,nmax=nmax,/large,nproc=nproc,almr=almr,outalmn=outalmn,rmin=0d0,rmax=1d0,/w8
almn2ball,falmn,fballx,nside=nside,nshell=nshell,lmax=lmax,nmax=nmax,nproc=nproc,almr=ralmr,rmin=0d0,rmax=1d0

fballxx=fballx

ratio_map=abs(fball_orig[*,1:nshell-1]-fballx[*,1:nshell-1]) ;*100/fball_orig[*,1:nshell-1]
rmax=5
ir=50
ipix=200
fname=figdir+'fr_ball_almn2ball_ipix'+strn(ipix+1)+'.ps'
planck_plot, r,reform(fball_orig[ipix,*]),ax=['xttl','yttl'],ay=['$r$','$f(r)$'],yover=reform(fballx[ipix,*]),fname=fname,/psend
fname=figdir+'fr_ball_almn2ball_diff_ipix'+strn(ipix+1)+'.ps'
planck_plot, r,reform(fball_orig[ipix,*]-fballx[ipix,*]),ax=['xttl','yttl'],ay=['$r$','$\Delta f(r)$'],fname=fname,/psend

fname=figdir+'ball_r'+strn(ir+1)+'.ps'
plot_mollview,fball_orig[*,ir],file_ps=fname ,title=' ',ctitle='zzttl',ax=['zzttl'],ay=['Original: r index='+strn(ir)]
fname=figdir+'almn2ball_r'+strn(ir+1)+'.ps'
plot_mollview,fballx[*,ir],file_ps=fname ,title=' ',ctitle='zzttl',ax=['zzttl'],ay=['almn2ball: r index='+strn(ir)]
fname=figdir+'ball_almn2ball_diff_r'+strn(ir+1)+'.ps'
plot_mollview,ratio_map[*,ir],min=0,max=0.1,file_ps=fname ,title=' ',ctitle='zzttl',ax=['zzttl'],ay=['Diff: r index='+strn(ir)]



map2alm,fball_orig[*,20],alm,lmax=lmax,/w8
alm2map,alm,falm,nside=nside
cl1=alm2cl(alm)



alm_map2=reform(almr[20,*,*,*])
cl2=alm2cl(alm_map2)
cl2_ratio=(cl1(2:lmax_sim)-cl2(2:lmax_sim))/cl1(2:lmax_sim)
print,'ball2almn to orig cl(r=20) check',cl2_ratio


alm_map22=reform(ralmr[20,*,*,*])
cl22=alm2cl(alm_map22)
cl22_ratio=(cl1(2:lmax_sim)-cl22(2:lmax_sim))/cl1(2:lmax_sim)
print,'almn2ball reverse cl(r=20) check',cl22_ratio

inn=2
outalm_n=reform(outalmn[0,*,*,inn])
cln=alm2cl(outalm_n)
print,'almn2ball reverse cln(n=2) check',cln


tek_color &  plot, r,reform(fball_orig[20,*]) & oplot, r,reform(fballxx[20,*]), color=2
tek_color &  plot, (reform(fball_orig[20,*])-reform(fballxx[20,*]))*100/reform(fballxx[20,*]),yra=[-0.5,0.5]


;;==========================================
;;         Needlet reconstruction
;;==========================================

;;check needlet reconstruction
ball2beta,fball_orig,fbeta,lmax=lmax,nmax=nmax,nproc=nproc,feedback=1,almn=almn_orig,outbeta=outbeta,almr=almr_ball,$
          rmin=0d0,rmax=1d0,nv=nshell,jint=j0,nj=nj,bb=bb,gln=gln,/w8, /norun

beta2ball,fbeta,fball,nrpix=nshell+1,nside=nside,lmax=lmax,nmax=nmax,nproc=nproc,feedback=3,almn=almn_beta,almr=almr_beta,$
               outmap=outmap,rmin=0d0,rmax=1d0,jint=j0,nv=nshell,nj=nj,bb=bb,/w8   ;,outbeta=betamap2 ;,/norun



alm_map3=reform(almr_ball[20,*,*,*])
cl3=alm2cl(alm_map3)
cl3_ratio=(cl1(2:lmax_sim)-cl3(2:lmax_sim))/cl1(2:lmax_sim)
print,'ball2beta to orig cl(r=20) check',cl3_ratio

alm_map4=reform(almr_beta[20,*,*,*])
cl4=alm2cl(alm_map4)
cl4_ratio=(cl1(2:lmax_sim)-cl4(2:lmax_sim))/cl1(2:lmax_sim)
print,'ball2beta to reverse cl(r=20) check',cl4_ratio

inn=2
origalm_n=reform(almn_orig[0,*,*,inn])
cln2=alm2cl(origalm_n)
print,'orig almn to cln(n=2) check',cln2

inn=2
betaalm_n=reform(almn_beta[0,*,*,inn])
cln3=alm2cl(betaalm_n)
print,'beta almn to cln(n=2) check',cln3



ratio_map=abs(fball_orig[*,1:nshell-1]-outmap[*,1:nshell-1]) ;*100/fball_orig[*,1:nshell-1]
rmax=5
ir=50
ipix=200
fname=figdir+'fr_ball_beta2ball_ipix'+strn(ipix+1)+'.ps'
planck_plot, r,reform(fball_orig[ipix,*]),ax=['xttl','yttl'],ay=['$r$','$f(r)$'],yover=reform(outmap[ipix,*]),fname=fname,/psend
fname=figdir+'fr_ball_beta2ball_diff_ipix'+strn(ipix+1)+'.ps'
planck_plot, r,reform(fball_orig[ipix,*]-outmap[ipix,*]),ax=['xttl','yttl'],ay=['$r$','$\Delta f(r)$'],fname=fname,/psend
fname=figdir+'beta2ball_r'+strn(ir+1)+'.ps'
plot_mollview,outmap[*,ir],file_ps=fname ,title=' ',ctitle='zzttl',ax=['zzttl'],ay=['beta2ball: r index='+strn(ir)]
fname=figdir+'ball_beta2ball_diff_r'+strn(ir+1)+'.ps'
plot_mollview,ratio_map[*,ir],max=0.1,file_ps=fname ,title=' ',ctitle='zzttl',ax=['zzttl'],ay=['Diff: r index='+strn(ir)]



tek_color &  plot, reform(fball_orig[20,*]) & oplot, reform(outmap[20,*]), color=2
tek_color &  plot, (reform(fball_orig[20,*])-reform(outmap[20,*]))/reform(outmap[20,*]),yra=[-1,1]


;;=============================================
;;    needlet gln^2 maps
;;===============================================
ball2beta,ball_sim,fbeta2,lmax=lmax,nmax=nmax,nproc=nproc,outbeta=outbeta3,feedback=1,almn=almn2_sim,almr=almr2_sim,$
          rmin=0d0,rmax=1d0,nv=nshell,jint=j0,nj=nj,bb=bb,gln=gln2,/w8,pow=2, /norun

;; get_3dneedlet,bb,j0,nj,lmax,nmax,gln
;; gln2=gln^2

beta2_sim=outbeta2
.r
almn_sim_jj=almn_sim
for jmap=0,nj-1 do begin
   for nn=0,nmax_sim do begin
      for ll=0,lmax_sim do begin
         almn_sim_jj[0,ll,0:lmax_sim,nn]=almn2_sim[0,ll,0:lmax_sim,nn]*gln2[ll,nn,jmap]
      endfor
   endfor
   
   almn2ball,almn_sim_jj,beta2,nside=nside,nshell=nshell,lmax=lmax,nmax=nmax,nproc=nproc,almr=almr_beta2,rmin=0d0,rmax=1d0
   beta2_sim[*,*,jmap]=beta2
endfor
wunf, beta2_sim,'temp/beta_gln2_sim.unf'
end


rmax=5
ir=50
ipix=200
.r
for jmap=0,nj-1 do begin
fname=figdir+'fr_ball2beta_ipix'+strn(ipix+1)+'_jmap'+strn(jmap)+'.ps'
planck_plot, r,reform(outbeta2[ipix,*,jmap]),ax=['xttl','yttl'],ay=['$r$','$f(r)$'],fname=fname,/psend

fname=figdir+'ball2beta_r'+strn(ir+1)+'_jmap'+strn(jmap)+'.ps'
plot_mollview,outbeta2[*,ir,jmap],file_ps=fname ,title=' ',ctitle='zzttl',ax=['zzttl'],ay=['r index='+strn(ir)+', j='+strn(jmap)]
endfor
end








for ii=0,nshell-1 do begin

   snr_last=0d0
   snr=0d0
   for jj=0,nmax_sim do begin
      ntilda=jj
      snr_last = snr_last + an_sim[jj]*sinc_x(argpi*ntilda*r[nshell])*argpi*ntilda 
      snr = snr + an_sim[jj]*sinc_x(argpi*ntilda*r[ii])*argpi*ntilda   
   endfor

   snrvec[ii]=snr  ;(snr-snr_last)
   fball_orig[*,ii]=y20[*,0]*(snr-snr_last)
   mono[ii] = mean(fball_orig[*,ii])
endfor
end

;; browse_ballpix,fball
;; browse_ballpix,fball_orig

;; fball_diff = (abs(fball-fball_orig))
;; print, 'mean diff = ',mean(fball_diff)

;; browse_ballpix,fball_diff

end




;;=======================

;; density, z, x, y, zmin=zmin, zmax=zmax, zlog=zlog, colorbar=colorbar, $
;;              nlevels=nlevels, fill=fill, bar_charsize=bar_charsize, $
;;              bar_pad=bar_pad, bar_title=bar_title, bar_format=bar_format, scientific=scientific, $
;;              _extra=extra
