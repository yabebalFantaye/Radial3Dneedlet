
ipix=20
fname=figdir+'knownball_orig_vs_recons_pix'+strn(ipix)+'.ps'
ps_start,file=fname
plot,outmap(ipix,*)*500. & oplot,fball_orig(ipix,*),li=1
ps_end

ipix=50
fname=figdir+'knownball_orig_vs_recons_pix'+strn(ipix)+'.ps'
ps_start,file=fname
plot,outmap(ipix,*)*500. & oplot,fball_orig(ipix,*),li=1
ps_end


;mollview, outmap[*,50]/fball_orig[*,50] 

;outmap2=total(betamap2,3)
ratio_map=abs(fball_orig[*,1:nshell-1]-outmap[*,1:nshell-1])*100/fball_orig[*,1:nshell-1]
;almn2ball,falmn,fball,nside=nside,nshell=nshell,lmax=lmax,nmax=nmax,/large,nproc=nproc,almr=ralmr,rmin=0.2,rmax=0.9d0

;return

x = indgen(nside2npix(nside))+1
y = indgen(nshell+1)

;;save 2D projected density maps

axisx = ['xlbl','ylbl','zlbl']
axisy = ['Healpix pixel number', 'Radial pixel number','$\Delta N/N [\%]$']
CTFILE = 'Planck_CT.tbl'

;;===================================
figdir='../knownball_nshell'+strn(nshell)+'_nmax'+strn(nmax)+'_nside'+strn(nside)+'_lmax'+strn(lmax)+'/'
spawn, 'mkdir -p '+figdir

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


;;=============================
fname=figdir+'seedlet_frac_diff_raidal_orig_vs_reconstruct_ball'+fkey+'.ps'
ps_start_planck,file=fname,/large,xmar=xmar,ymar=ymar,ytdx=ytdx,xtdy=xtdy,xsize=xsize,ysize=ysize

tek_color

fball=outmap
plot, y,(fball_orig[10,*]-fball[10,*])*100./fball_orig[10,*], title=' ',color=0,$
             xtitle='ylbl',ytitle='zlbl',xra=[0,nshell],yra=[-100,100],/xs,/ys

oplot, y,(fball_orig[100,*]-fball[100,*])*100./fball_orig[100,*],color=2
oplot, y,(fball_orig[500,*]-fball[500,*])*100./fball_orig[500,*],color=3
oplot, y,(fball_orig[1000,*]-fball[1000,*])*100./fball_orig[1000,*],color=4
oplot, y,(fball_orig[2000,*]-fball[2000,*])*100./fball_orig[2000,*],color=5


legend,['pixel=10','pixel=100','pixel=500','pixel=1000','pixel=2000'],color=[0,2,3,4,5],textcolor=[0,2,3,4,5],/top,/left,box=0,$ 
       line=[0,0,0,0,0],pspacing=1,thick=4,charsize=lcharsize

ps_end_planck, /png,feps=fname,/latex,xsize=xsize,ysize=ysize,axisx=axisx,axisy=axisy
;;--------------------------




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


end
