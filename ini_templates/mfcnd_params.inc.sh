cat > $parsname <<EOF

feedback=1

#------------
#the following parameters are for post-processing
figdir=$figdir
fnameplot=$fnameplot
exename=$exename

#------------

nside =$nside
lmax = $lmax
lmaxgen=$lmax
lmaxin=$lmax

tmaskext=$tmaskext
pmaskext=$pmaskext
tempext=$tempext
polext=$polext
dot=$dot
dop=$dop
doqu=$doqu
mpolar=$mpolar

mfsinglemode=0
mfneedlet=1
nj=$nj
j0=$j0
bb=$bb

glpow=$glpow
mpow=$mpow

nsteps=26
nu_min=-3d0
nu_max=3d0

cmbnorm = $simcmbnorm
noisenorm = $noisenorm
wmapmode = 0

modulate = $modulate
#in degree
modulfwhm = $modulfwhm
modulth = $modulth
modulph = $modulph
modulamp = $modulamp

boost = $boost
boost_factor = $bv 

apodisation = $apodisation
#nj_sigma used is read (apod1) /written(apod.ne.2) on dirvecfile
dirvecfile=$odir/nj_sigma_${data}_${mask}_nmaps${nnsim}.unf_norm_njsigma_used

polar = $polar
clfile = $clfile

mask_map=${mask_map}
nmask = $nmask
maskfile = $maskfile


fwhm_arcmin = $fwhm
nwindow = 1
windowfile =  $windowfile

nsim=$nnsim
nsim_start=${nsim_start}
iseed= -1282287

mapfile =  $mapfile

noisemapfile=$noisemapfile

datafile = $datafile
noisedatafile=$noisedatafile

datafileout = $datafileout

tempdir=$tempdir

EOF

#=========================================
