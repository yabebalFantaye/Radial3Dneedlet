
cl_file='../cls/Cl_model_bin5.dat'
readcol,cl_file, ll,clmm,clmm_err

grcn=2  ;grcn=1 - g(r) from power law C_n; grcn=2 - g(r) = (r-0.5)^n
nmax_sim=20

nside=64
nshell=128 ;;number of radial shells 

nmax=nshell ;;maximum radial multipole
lmax= 120 ;;angular multipole
lmax_sim=lmax

;;input angular power spectra and radial coefficients
clmm=1e5*clmm


bb=2
j0=-1
nj = ceil(alog(sqrt(1.*(lmax*(lmax+1.)+(1.*nmax_sim)^2)))/alog(1.*bb)+1)
print, 'nj = ',nj
nproc=1


fkey = '_grcn'+strn(grcn)+'_lmax'+strn(lmax)+'_nmax'+strn(nshell/2)+'_nr'+strn(nshell)

figdir='../figsOct14/fftw'+fkey+'/'
spawn, 'mkdir -p '+figdir

fball = '../output/ball_clmm'+fkey+'.unf'
fbeta = '../output/beta_clmm'+fkey+'.unf'
fbeta2 = '../output/beta_gln2_clmm'+fkey+'.unf'
falmn = '../output/almn_clmm'+fkey+'.unf'


npix=nside2npix(nside)
r=dindgen(nshell)/double(nshell)
rmax=5
ir=50
ipix=200

;;get an example cl and an to construct input almn
cl_sim=dblarr(lmax+1)
cl_sim[2:lmax_sim]=clmm[1:lmax_sim-1]

;;simul_fftwate band limited ball
synball_fftw,an,cl_sim,ball_sim,almn=almn_sim,lmax=lmax_sim,nshell=nshell,nside=nside,almr=almr_sim,$
             fr=fr,cn=cn,rfr=rfr,ran=ran,rcn=rcn,grcn=grcn,nmax_sim=nmax_sim


;;create angular map
create_alm,cl_sim,alm_sim,lmax=lmax,12312




;;==========================================
;;         Radial space  reconstruction
;;==========================================

nvec = indgen(nshell/2+1)
lvec = indgen(lmax+1)

;;Radial plots: C_n and g(r)
ind_an=where(abs(an-ran) gt 1e-20)
planck_plot, ind_an,abs(an(ind_an)-ran(ind_an))*1e10/abs(an(ind_an)),y=[0,0.2],xra=[1,nshell/2],ax=['xttl','yttl'],ay=['$n$','$10^{10}\Delta |a_n|/|a_n|$'],fname=figdir+'an_diff_orig_recon.ps',/psend
planck_plot, nvec,cn,yra=[0,10],xra=[1,nshell/2],ax=['xttl','yttl'],ay=['$n$','$C_n$'],yover=rcn,fname=figdir+'Cn_orig_recon.ps',/psend,/xlog,/xstyle
planck_plot, r,fr,ax=['xttl','yttl'],ay=['$r$','$g(r)$'],yover=rfr,fname=figdir+'fr_orig_recon.ps',/psend
;;planck_plot, r,(fr-rfr)*1e3/fr,ax=['xttl','yttl'],ay=['$r$','$10^{3}\Delta g(r)/g(r)$'],fname=figdir+'fr_diff_orig_recon.ps',/psend


;;==========================================
;;         Harmonic space  reconstruction
;;==========================================
;;checl radial and angular part without needlet
ball2almn,ball_sim,falmn,lmax=lmax,nmax=nshell,/large,nproc=nproc,almr=almr,outalmn=outalmn,/w8
monopole=0d0;mean(ball_sim)
almn2ball,falmn,fballx,nside=nside,nshell=nshell,lmax=lmax,nmax=nshell,nproc=nproc,almr=ralmr,monopole=monopole

tek_color & plot, r,reform(ball_sim[ipix,*]),color=1 & oplot, r,reform(fballx[ipix,*]),color=2


;;max_diff_almn2ball = max(abs(almn_sim-falmn),dimension=1) 

cln = avg(abs(almn_sim)^2,2)
rcln = avg(abs(falmn)^2,2)

ratio_cln = ((cln)-rcln)/(cln)
ratio_cln[where(finite(ratio_cln) eq 0)]=0d0

ratio_cn = max(ratio_cln,dimension=2) 
ratio_cl = max(ratio_cln,dimension=1) 

nvec = indgen(nshell/2+1)
lvec = indgen(lmax+1)
planck_plot, nvec,ratio_cn[nvec]*1e10,yra=[1.0,1.5],xra=[0,nshell/2],ax=['xttl','yttl'],ay=['$n$','$10^{10}(\Delta \hat{C}_n/\hat{C}_n)^{max}$'],fname=figdir+'Xn_diff_ball2alm.ps',/psend
planck_plot, lvec,ratio_cl[lvec]*1e10,yra=[0,2],xra=[0,lmax],ax=['xttl','yttl'],ay=['$\ell$','$10^{10}(\Delta \hat{C}_{\ell}/\hat{C}_{\ell})^{max}$'],fname=figdir+'Xl_diff_ball2alm.ps',/psend




;;ball plots
ratio_map=(ball_sim[*,0:nshell-1]-fballx[*,0:nshell-1])/ball_sim[*,0:nshell-1]

;;fr plots
fname=figdir+'fr_ball_almn2ball_ipix'+strn(ipix+1)+'.ps'
planck_plot, r,reform(ball_sim[ipix,*]),ax=['xttl','yttl'],ay=['$r$','$f(r)$'],yover=reform(fballx[ipix,*]),fname=fname,/psend

fname=figdir+'fr_ball_almn2ball_diff_ipix'+strn(ipix+1)+'.ps'
planck_plot, r,(reform(ratio_map[ipix,*]))*1e10,ax=['xttl','yttl'],ay=['$r$','$10^{12}\Delta g(r)/g(r)$'],fname=fname,/psend

;;mollview plots
fname=figdir+'ball_r'+strn(ir+1)+'.ps'
plot_mollview,ball_sim[*,ir],file_ps=fname ,title=' ',ctitle='zzttl',ax=['zzttl'],ay=['Original: r index='+strn(ir)]
fname=figdir+'almn2ball_r'+strn(ir+1)+'.ps'
plot_mollview,fballx[*,ir],file_ps=fname ,title=' ',ctitle='zzttl',ax=['zzttl'],ay=['almn2ball: r index='+strn(ir)]
fname=figdir+'ball_almn2ball_diff_r'+strn(ir+1)+'.ps'
plot_mollview,ratio_map[*,ir]*1e9,min=-1,max=1,file_ps=fname ,title=' ',ctitle='zzttl',ax=['zzttl'],ay=['Frac. Diff [$\times 10^{9}$]: r index='+strn(ir)]


;;angular plots: C_l
ellvec = indgen(lmax-1)+2

map2alm,ball_sim[*,20],alm,lmax=lmax,/w8
alm2map,alm,falm,nside=nside
cl0=alm2cl(alm)

alm_map2=reform(almr[20,*,*])
cl1=alm2cl(alm_map2)
cl10_ratio=(cl0(2:lmax_sim)-cl1(2:lmax_sim))/cl0(2:lmax_sim)
print,'ball2almn to orig cl(r=20) check',cl10_ratio


alm_map22=reform(ralmr[20,*,*])
cl2=alm2cl(alm_map22)
cl20_ratio=(cl0(2:lmax_sim)-cl2(2:lmax_sim))/cl0(2:lmax_sim)
print,'almn2ball reverse cl(r=20) check',cl20_ratio

planck_plot, ellvec,cl10_ratio,xra=[1,lmax],ax=['xttl','yttl'],ay=['$\ell$','$\Delta C_{\ell}/C_{\ell}$'],fname=figdir+'Cl_diff_orig_vs_ball2almr.ps',/psend
planck_plot, ellvec,cl20_ratio,xra=[1,lmax],ax=['xttl','yttl'],ay=['$\ell$','$\Delta C_{\ell}/C_{\ell}$'],fname=figdir+'Cl_diff_orig_vs_almn2almr.ps',/psend





tek_color &  plot, r,reform(ball_sim[20,*]) & oplot, r,reform(fballxx[20,*]), color=2
tek_color &  plot, (reform(ball_sim[20,*])-reform(fballxx[20,*]))*100/reform(fballxx[20,*]),yra=[-0.5,0.5]


;;==========================================
;;         Needlet reconstruction
;;==========================================


ball2beta,ball_sim,fbeta,lmax=lmax,nmax=nmax,nproc=nproc,feedback=1,almn=almn_orig,outbeta=outbeta,almr=almr_ball,$
          jint=j0,nj=nj,bb=bb,gln=gln,/w8 ;,pow=2

;monopole=mean(ball_sim)
beta2ball,fbeta,fball,nshell=nshell,nside=nside,lmax=lmax,nmax=nmax,nproc=nproc,feedback=3,almn=almn_beta,almr=almr_beta,$
               outmap=outmap,jint=j0,nj=nj,bb=bb,/w8  ;,monopole=monopole

tek_color &  plot, reform(ball_sim[20,*]),color=1 & oplot, reform(outmap[20,*]), color=2
tek_color &  plot, reform(ball_sim[20,*]-1.3591841*outmap[20,*]),color=1 

;;--
almn2ball,almn_orig,fballx1,nside=nside,nshell=nshell,lmax=lmax,nmax=nmax,nproc=nproc,almr=ralmr,monopole=monopole
almn2ball,almn_beta,fballx2,nside=nside,nshell=nshell,lmax=lmax,nmax=nmax,nproc=nproc,almr=ralmr,monopole=monopole
plot, fballx2[ipix,*],color=1 & oplot, fballx1[ipix,*],color=2 & oplot, fballx[ipix,*],color=3


;;----
 
alm_map3=reform(almr_ball[20,*,*])
cl3=alm2cl(alm_map3)

alm_map4=reform(almr_beta[20,*,*])
cl4=alm2cl(alm_map4)


cl3_ratio=abs(cl3(2:lmax_sim)-cl4(2:lmax_sim))
print,'ball2beta to orig cl(r=20) check',cl3_ratio



inn=2
origalm_n=reform(almn_sim[inn,*,*])
cln1=alm2cl(origalm_n)
print,'sim almn to cln(n=2) check',cln1

origalm_n=reform(almn_orig[inn,*,*])
cln2=alm2cl(origalm_n)
print,'ball2almn almn to cln(n=2) check',cln2

inn=2
betaalm_n=reform(almn_beta[inn,*,*])
cln3=alm2cl(betaalm_n)
print,'beta almn to cln(n=2) check',cln3



ratio_map=(ball_sim[*,1:nshell-1]-outmap[*,1:nshell-1])/ball_sim[*,1:nshell-1]

fname=figdir+'fr_ball_beta2ball_ipix'+strn(ipix+1)+'.ps'
planck_plot, r,reform(ball_sim[ipix,*]),ax=['xttl','yttl'],ay=['$r$','$g(r)$'],yover=reform(outmap[ipix,*]),fname=fname,/psend

fname=figdir+'fr_ball_beta2ball_diff_ipix'+strn(ipix+1)+'.ps'
planck_plot, r,reform((ratio_map[*,ir]*1e3),ax=['xttl','yttl'],ay=['$r$','$10^{3}\Delta g(r)/g(r)$'],fname=fname,/psend

fname=figdir+'beta2ball_r'+strn(ir+1)+'.ps'
plot_mollview,outmap[*,ir],file_ps=fname ,title=' ',ctitle='zzttl',ax=['zzttl'],ay=['beta2ball: r index='+strn(ir)]

fname=figdir+'ball_beta2ball_diff_r'+strn(ir+1)+'.ps'
plot_mollview,ratio_map[*,ir]*1e2,min=-0.1,max=0.1,file_ps=fname ,title=' ',ctitle='zzttl',ax=['zzttl'],ay=['Frac. Diff [\%]: r index='+strn(ir)]




tek_color &  plot, (reform(ball_sim[20,*])-reform(outmap[20,*]))/reform(outmap[20,*]),yra=[-1,1]


;;=============================================
;;    needlet gln^2 maps
;;===============================================
;; ball2beta,ball_sim,fbeta2,lmax=lmax,nmax=nmax,nproc=nproc,outbeta=outbeta3,feedback=1,almn=almn2_sim,almr=almr2_sim,$
;;           rmin=0d0,rmax=1d0,nv=nshell,jint=j0,nj=nj,bb=bb,gln=gln2,/w8,pow=2, /norun



ball2beta,ball_sim,fbeta,lmax=lmax,nmax=nmax,nproc=nproc,feedback=1,outbeta=outbeta2,$
          jint=j0,nj=nj,bb=bb,gln=gln2,/w8,pow=2

;; get_3dneedlet,bb,j0,nj,lmax,nmax,gln
;; gln2=gln^2

;; monopole=mean(ball_sim)

;; beta2_sim=sc_outbeta
;; .r
;; almn_sim_jj=outalmn

;; for jmap=0,nj-1 do begin

;;    for nn=0,nmax_sim do begin
;;       for ll=0,lmax_sim do begin
;;          almn_sim_jj[0,ll,0:lmax_sim,nn,*]=outalmn[0,ll,0:lmax_sim,nn,*]*gln2[ll,nn,jmap]
;;       endfor
;;    endfor


;;    almn2ball_sincos,almn_sim_jj,beta2,nside=nside,nshell=nshell,lmax=lmax,nmax=nmax,nproc=nproc,almr=almr_b2,rmin=0d0,rmax=1d0,$
;;                     sinoutmap=sinoutmap,cosoutmap=cosoutmap,monopole=monopole
   
;;    beta2_sim[*,*,jmap]=beta2
;; endfor

;; wunf, beta2_sim,'temp/beta_gln2_sim.unf'
;; end



.r
for jmap=0,nj-1 do begin
fname=figdir+'fr_ball2beta_ipix'+strn(ipix+1)+'_jmap'+strn(jmap)+'.ps'
planck_plot, r,reform(outbeta2[ipix,*,jmap]),ax=['xttl','yttl'],ay=['$r$','$f(r)$'],fname=fname,/psend

fname=figdir+'ball2beta_r'+strn(ir+1)+'_jmap'+strn(jmap)+'.ps'
plot_mollview,outbeta2[*,ir,jmap],file_ps=fname ,title=' ',ctitle='zzttl',ax=['zzttl'],ay=['r index='+strn(ir)+', j='+strn(jmap)]
endfor
end







end




;;=======================

;; density, z, x, y, zmin=zmin, zmax=zmax, zlog=zlog, colorbar=colorbar, $
;;              nlevels=nlevels, fill=fill, bar_charsize=bar_charsize, $
;;              bar_pad=bar_pad, bar_title=bar_title, bar_format=bar_format, scientific=scientific, $
;;              _extra=extra
