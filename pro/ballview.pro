;pro ballview,filename



maxcount_mocbox=519

filename='../output/ballmap_mockbox.unf'


browse_ballpix,filename,norm=1./float(maxcount_mocbox)



;window, /free
;gLoadCT, 0
;gSurf, arr, /Shaded, Shades=BytScl(arr)
;gLoadCT, 25, /Brewer, /Reverse
;gSurf, arr, Shades=BytScl(arr), /NoErase, Title='Surface Title'

end
