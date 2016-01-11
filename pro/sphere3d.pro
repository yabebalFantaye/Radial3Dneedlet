PRO sphere3d

;; Create random data in spherical coordinates.

omega = !RaDeg * ASIN(-1.0 + 2.0*RANDOMU(seed, 100))
rho = RANDOMU(seed, 100) * 360.0
radius = RANDOMU(seed, 100)
value = COS(omega)^2/SIN(rho)^2

;; Convert to retangular coordinates.

sphericalCoords = FLTARR(3,100)
sphericalCoords(0,*) = rho
sphericalCoords(1,*) = omega
sphericalCoords(2,*) = radius
rectCoords = CV_COORD(From_Sphere=sphericalCoords, /To_Rect, /Degrees)

;; Grid the rectangular coordinates.

gridded = GRID3(rectCoords(0,*), rectCoords(1,*), rectCoords(2,*), value)
COMMON Volume_Data, xx
xx = gridded
SLICER
END
