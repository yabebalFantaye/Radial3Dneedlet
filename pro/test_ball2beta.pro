maxcount_mocbox=519
bb=2

nrint=64
nside=16
nshell=2*nrint ;;number of radial shells 

nmax=50 ;;maximum radial multipole
lmax=60 ;;angular multipole
nproc=1 
nrsample=1 ;;r sampling factor for ball2almn
figdir='../vls_figures_sincos2_nshell'+strn(nshell)+'_nmax'+strn(nmax)+'_nside'+strn(nside)+'_lmax'+strn(lmax)+'/'
spawn, 'mkdir -p '+figdir

j0=0
nj = ceil(alog(sqrt(1.*(lmax*(lmax+1.)+(1.*nmax)^2)))/alog(1.*bb)+1)
print, 'nj = ',nj


r=dindgen(nshell+1)/double(nshell)

use_orig=2
.r
if use_orig eq 1 then begin
   fball_orig = '../output/ballmap_mockbox_nrint'+strn(nrint)+'.unf'
   fball = '../output/ballmap_mockbox_almn2ball_nrint'+strn(nrint)+'_1.unf'
   fbeta = '../output/beta_mockbox_1.unf'
   float=1
endif else begin
   ;;fball_orig =
   ;;'../output/ballmap_mockbox_almn2ball_nrint'+strn(nrint)+'_2.unf'
   fball_orig = '../output/ballmap_mockbox_almn2ball_nrint'+strn(nrint)+'_1.unf'
   fball = '../output/ballmap_mockbox_almn2ball_nrint'+strn(nrint)+'_2.unf'
   fbeta = '../output/beta_mockbox_3.unf'
   float=0
endelse
end
fkey = '_nshell'+strn(nshell)+'_nmax'+strn(nmax)




;;ball2ball,fball_orig,fbeta,lmax=lmax,nmax=nnmax,nproc=nproc,feedback=1,almn=almn_orig,outbeta=outbeta,$
;;          rmin=0d0,rmax=1d0,jint=j0,nj=nj,bb=bb,gln=gln,mapjkn=mapjkn,dmapjkn=dmapjkn,snsnp=snsnp ;,/norun


ball2beta,fball_orig,fbeta,lmax=lmax,nmax=nmax,nproc=nproc,feedback=1,almn=almn_orig,outbeta=outbeta,almr=almr_ball,$
          rmin=0d0,rmax=1d0,nv=nshell,jint=j0,nj=nj,bb=bb,gln=gln,/w8, float=float, /norun

beta2ball,fbeta,fball,nrpix=nshell+1,nside=nside,lmax=lmax,nmax=nmax,nproc=nproc,feedback=3,almn=almn_beta,almr=almr_beta,$
               outmap=outmap,rmin=0d0,rmax=1d0,jint=j0,nv=nshell,nj=nj,bb=bb,/w8   ;,outbeta=betamap2 ;,/norun



sincos_ball2beta,fball_orig,fbeta,lmax=lmax,nmax=nmax,nproc=nproc,feedback=1,almn=sc_almn_orig,outbeta=sc_outbeta,almr=sc_almr_ball,$
          rmin=0d0,rmax=1d0,nv=nshell,jint=j0,nj=nj,bb=bb,gln=gln,float=float,/w8

monopole=mean(fball_orig)
sincos_beta2ball,fbeta,fball,nrpix=nshell+1,nside=nside,lmax=lmax,nmax=nmax,nproc=nproc,feedback=3,almn=sc_almn_beta,almr=sc_almr_beta,$
               outmap=outmap,rmin=0d0,rmax=1d0,jint=j0,nv=nshell,nj=nj,bb=bb,/w8,sinoutmap=sinoutmap,cosoutmap=cosoutmap,monopole=monopole




ratio_map=abs(fball_orig[*,1:nshell-1]-outmap[*,1:nshell-1]) ;*100/fball_orig[*,1:nshell-1]
rmax=5
ir=40
ipix=800



map2alm,fball_orig[*,ir],alm,lmax=lmax,/w8
alm2map,alm,map_ir,nside=nside

fname=figdir+'fr_ball_beta2ball_ipix'+strn(ipix+1)+'.ps'
planck_plot, r,reform(fball_orig[ipix,*]),ax=['xttl','yttl'],ay=['$r$','$f(r)$'],yover=reform(outmap[ipix,*]),fname=fname,/psend
fname=figdir+'fr_ball_beta2ball_diff_ipix'+strn(ipix+1)+'.ps'
planck_plot, r,reform(fball_orig[ipix,*]-outmap[ipix,*]),ax=['xttl','yttl'],ay=['$r$','$\Delta f(r)$'],fname=fname,/psend


fname=figdir+'ball_r'+strn(ir+1)+'.ps'
plot_mollview,fball_orig[*,ir],file_ps=fname ,title=' ',ctitle='zzttl',ax=['zzttl'],ay=['Original: r index='+strn(ir)]
fname=figdir+'ball_healpix_r'+strn(ir+1)+'.ps'
plot_mollview,map_ir,file_ps=fname ,title=' ',ctitle='zzttl',ax=['zzttl'],ay=['alm2map: r index='+strn(ir)]
fname=figdir+'beta2ball_r'+strn(ir+1)+'.ps'
plot_mollview,outmap[*,ir],file_ps=fname ,title=' ',ctitle='zzttl',ax=['zzttl'],ay=['beta2ball: r index='+strn(ir)]
fname=figdir+'ball_beta2ball_diff_r'+strn(ir+1)+'.ps'
plot_mollview,ratio_map[*,ir],max=0.1,file_ps=fname ,title=' ',ctitle='zzttl',ax=['zzttl'],ay=['Diff: r index='+strn(ir)]



tek_color &  plot, reform(fball_orig[20,*]) & oplot, reform(outmap[20,*]), color=2
tek_color &  plot, (reform(fball_orig[20,*])-reform(outmap[20,*]))/reform(outmap[20,*]),yra=[-1,1]




;;;=======================
;;save 2D projected density maps

axisx = ['xlbl','ylbl','zlbl']
axisy = ['Healpix pixel number', 'Radial pixel number','$\Delta N/N [\%]$']
CTFILE = 'Planck_CT.tbl'

;;===================================
rmax=5
ir=10
fname=figdir+'seedlet_orig_ball_r'+strn(ir+1)+'.ps'
plot_mollview,fball_orig[*,ir],file_ps=fname ,title=' ',ctitle='zzttl',ax=['zzttl'],ay=['Original ball: r index='+strn(ir)]

ir=10
fname=figdir+'seedlet_reconstruct_ball_r'+strn(ir+1)+'.ps'
plot_mollview,outmap[*,ir],file_ps=fname ,title=' ',ctitle='zzttl',ax=['zzttl'],ay=['(Reconstructed ball: r index='+strn(ir)]

ir=100
fname=figdir+'seedlet_orig_ball_r'+strn(ir+1)+'.ps'
plot_mollview,fball_orig[*,ir],file_ps=fname ,title=' ',ctitle='zzttl',ax=['zzttl'],ay=['Original ball: r index='+strn(ir)]

ir=100
fname=figdir+'seedlet_reconstruct_ball_r'+strn(ir+1)+'.ps'
plot_mollview,outmap[*,ir],file_ps=fname ,title=' ',ctitle='zzttl',ax=['zzttl'],ay=['(Reconstructed ball: r index='+strn(ir)]


;;==========================

rmax=5
ir=0
fname=figdir+'seedlet_orig_vs_reconstruct_ball_r'+strn(ir+1)+'.ps'
plot_mollview,ratio_map[*,ir],0,rmax,file_ps=fname ,title=' ',ctitle='zzttl',ax=['zzttl'],ay=['(In-Rec)/In [\%]']

ir=20
fname=figdir+'seedlet_orig_vs_reconstruct_ball_r'+strn(ir+1)+'.ps'
plot_mollview,ratio_map[*,ir],0,rmax,file_ps=fname ,title=' ',ctitle='zzttl',ax=['zzttl'],ay=['(In-Rec)/In [\%]']

ir=60
fname=figdir+'seedlet_orig_vs_reconstruct_ball_r'+strn(ir+1)+'.ps'
plot_mollview,ratio_map[*,ir],0,rmax,file_ps=fname ,title=' ',ctitle='zzttl',ax=['zzttl'],ay=['(In-Rec)/In [\%]']

ir=100
fname=figdir+'seedlet_orig_vs_reconstruct_ball_r'+strn(ir+1)+'.ps'
plot_mollview,ratio_map[*,ir],0,rmax,file_ps=fname ,title=' ',ctitle='zzttl',ax=['zzttl'],ay=['(In-Rec)/In [\%]']

ir=nshell-2
fname=figdir+'seedlet_orig_vs_reconstruct_ball_r'+strn(ir+1)+'.ps'
plot_mollview,ratio_map[*,ir],0,rmax,file_ps=fname ,title=' ',ctitle='zzttl',ax=['zzttl'],ay=['(In-Rec)/In [\%]']




;; browse_ballpix,fball
;; browse_ballpix,fball_orig

;; fball_diff = (abs(fball-fball_orig))
;; print, 'mean diff = ',mean(fball_diff)

;; browse_ballpix,fball_diff

end




;;----------------------------------

;; ls_mollview, outbeta[*,1,1],ps=figdir+'beta_j1_v1.ps',title='B_jkv map: j=1, v=1'
;; ls_mollview, outbeta[*,2,2],ps=figdir+'beta_j2_v2.ps',title='B_jkv map: j=2, v=2'
;; ls_mollview, outbeta[*,2,3],ps=figdir+'beta_j2_v3.ps',title='B_jkv map: j=2, v=3'
;; ls_mollview, outbeta[*,1,4],ps=figdir+'beta_j1_v4.ps',title='B_jkv map: j=1, v=4'
;; ls_mollview, outbeta[*,4,1],ps=figdir+'beta_j4_v1.ps',title='B_jkv map: j=4, v=1'





;fball_orig[*,0]=0.01
;fball[*,0]=0.01

;; ;;--------------------------
;; fname=figdir+'test_original_ball'+fkey+'.ps'
;; ps_start_planck,file=fname,/large,xmar=xmar,ymar=ymar,ytdx=ytdx,xtdy=xtdy,xsize=xsize,ysize=ysize

;; hfi_ct, CTDIR='./', CTFILE=CTFILE, /LOAD

;; density, fball_orig, x, y, /colorbar, title='Original 3D ball: from VLS simulation Ngal = 1841599',$
;;              bar_charsize=2, xtitle='xlbl',ytitle='ylbl',xra=[1,nside2npix(nside)],yra=[1,nshell],/xs,/ys

;; ps_end_planck, /png,feps=fname,/latex,xsize=xsize,ysize=ysize,axisx=axisx,axisy=axisy

;; ;;--------------------------
;; fname=figdir+'test_reconstructed_ball'+fkey+'.ps'
;; ps_start_planck,file=fname,/large,xmar=xmar,ymar=ymar,ytdx=ytdx,xtdy=xtdy,xsize=xsize,ysize=ysize

;; hfi_ct, CTDIR='./', CTFILE=CTFILE, /LOAD

;; density, fball, x, y, /colorbar, title='Reconstructed 3D ball',$
;;              bar_charsize=2, xtitle='xlbl',ytitle='ylbl',xra=[1,nside2npix(nside)],yra=[1,nshell],/xs,/ys

;; ps_end_planck, /png,feps=fname,/latex,xsize=xsize,ysize=ysize,axisx=axisx,axisy=axisy

;; ;;--------------------------
;; fname=figdir+'test_difference_ball'+fkey+'.ps'
;; ps_start_planck,file=fname,/large,xmar=xmar,ymar=ymar,ytdx=ytdx,xtdy=xtdy,xsize=xsize,ysize=ysize

;; hfi_ct, CTDIR='./', CTFILE=CTFILE, /LOAD

;; density, (fball_orig - fball), x, y, zmin=-2, zmax=2, /colorbar, title='Difference  3D ball',$
;;              bar_charsize=2, xtitle='xlbl',ytitle='ylbl',xra=[1,nside2npix(nside)],yra=[1,nshell],/xs,/ys

;; ;*100./(fball_orig)

;; ps_end_planck, /png,feps=fname,/latex,xsize=xsize,ysize=ysize,axisx=axisx,axisy=axisy
;; ;;---------------------------


;;=======================

;; density, z, x, y, zmin=zmin, zmax=zmax, zlog=zlog, colorbar=colorbar, $
;;              nlevels=nlevels, fill=fill, bar_charsize=bar_charsize, $
;;              bar_pad=bar_pad, bar_title=bar_title, bar_format=bar_format, scientific=scientific, $
;;              _extra=extra
