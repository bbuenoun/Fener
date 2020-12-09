#!/usr/bin/env python
# -*- coding: utf-8 -*-


from numpy import ones, zeros
import Ztran


def ztran(ncap, esp, con, volCap, dt):

    den = ones(ncap) * 1000.0
    res = ones(ncap) * 0.0
    ces = zeros(ncap)
    for i in range(ncap):
        ces[i] = volCap[i] / 1000.0

    umuro, a, c, b, d, ns = Ztran.ztran(esp, con, den, ces, res, dt)

    #  Truncamiento
    for i in range(ns):
        if abs(d[i]) > 1e-5:
            nd = i + 1
    for i in range(ns):
        if abs(a[i]) > 1e-5:
            nc = i + 1

    # Correccion para el permanente
    sumaa = 0.0
    sumab = 0.0
    sumac = 0.0
    for i in range(nc - 1):
        sumaa = sumaa + a[i]
        sumab = sumab + b[i]
        sumac = sumac + c[i]
    sumad = 0.0
    for i in range(nd):
        sumad = sumad + d[i]
    a[nc - 1] = umuro * sumad - sumaa
    b[nc - 1] = umuro * sumad - sumab
    c[nc - 1] = umuro * sumad - sumac

    return umuro, a, c, b, d, nc, nd, ns


# ~ dt = 900.

# ~ esp = 0.1 #[.100/ncap,.100/ncap,.100/ncap,.100/ncap,.100/ncap,.100/ncap,.100/ncap,.100/ncap,.100/ncap,.100/ncap]
# ~ con =  1. #[1.000/ncap,1.000/ncap,1.000/ncap,1.000/ncap,1.000/ncap,1.000/ncap,1.000/ncap,1.000/ncap,1.000/ncap,1.000/ncap]
# ~ den = 1000.0 #[1000.0/ncap,1000.0/ncap,1000.0/ncap,1000.0/ncap,1000.0/ncap,1000.0/ncap,1000.0/ncap,1000.0/ncap,1000.0/ncap,1000.0/ncap]
# ~ ces = 1000. #[1000.0/ncap,1000.0/ncap,1000.0/ncap,1000.0/ncap,1000.0/ncap,1000.0/ncap,1000.0/ncap,1000.0/ncap,1000.0/ncap,1000.0/ncap]
# ~ res = 0.0 #[.100/ncap,.100/ncap,.100/ncap,.100/ncap,.100/ncap,.100/ncap,.100/ncap,.100/ncap,.100/ncap,.100/ncap]

# esp = [ 0.01,0.0475,0.095,0.0475,0.01]
# con = [ 0.25,2.3,2.3,2.3 ,0.25]
# volCap = [  945000.,2300000.,2300000.,2300000.,945000.]

# ~ esp = [ 0.01  , 0.0475 ,0.095 , 0.0475,0.035  ,0.07, 0.035 ,0.01  ]
# ~ con = [ 0.25, 2.3,  2.3,  2.3,  0.03, 0.03,0.03, 0.25]
# ~ volCap = [  945000., 2300000., 2300000., 2300000.  , 29400. ,  29400. , 29400. , 945000.]

# ~ [ 0.01  0.05  0.01] [ 0.25  2.3   0.25] [  945000.  2300000.   945000.]
# ~ [ 0.01  0.05  0.01] [ 0.25  2.3   0.25] [  945000.  2300000.   945000.]
# ~ [ 0.01  0.05  0.01] [ 0.25  2.3   0.25] [  945000.  2300000.   945000.]
# ~ [ 0.01    0.0475  0.095   0.0475  0.01  ] [ 0.25  2.3   2.3   2.3   0.25] [  945000.  2300000.  2300000.  2300000.   945000.]

# ~ ncap = 5
# ~ esp = [ 0.01 ,  0.0475,  0.095,  0.0475, 0.01  ]
# ~ con = [ 0.25, 2.3,  2.3,  2.3,  0.25]
# ~ volCap = [  945000., 2300000., 2300000., 2300000.,  945000.]

# ~ u = 0
# ~ for i in range(ncap):
# ~ u = u +esp[i]/con[i]
# ~ print(1/u)
# ~ ncap = 3
# ~ esp = [0.0475,  0.095,  0.0475]
# ~ con = [2.3,  2.3,  2.3]
# ~ volCap = [2300000., 2300000., 2300000.]

# ~ umuro,a,c,b,d,nc,nd,ns = ztran(ncap,esp[::-1],con[::-1],volCap[::-1],dt)
# ~ print(umuro)
# ~ print(nd,d)
# ~ print(nc,a)
# ~ print(nc,b)
# ~ print(nc,c)
# ~ print(ns)

# ~ print((298-293)/umuro)

# ~ qse = zeros(145) #ones(145)*(293-287)/umuro
# ~ qsi = zeros(145) #ones(145)*(287-293)/umuro
# ~ tse = ones(145)*281
# ~ tsi = ones(145)*281
# ~ tse[1] = 281.
# ~ tsi[1] = 284.

# ~ i=1
# ~ for j in range(nc):
# ~ qse[i] = qse[i]-a[j]*tse[i-j]+b[j]*tsi[i-j]
# ~ qsi[i] = qsi[i]+b[j]*tse[i-j]-c[j]*tsi[i-j]
# ~ for j in range(1,nd):
# ~ qse[i] = qse[i]-d[j]*qse[i-j]
# ~ qsi[i] = qsi[i]-d[j]*qsi[i-j]
# ~ print(qse[1],qsi[1])

# ---------------------------------------------
# ~ qse = zeros(145) #*(293-287)/umuro
# ~ qsi = zeros(145) #*(287-293)/umuro
# ~ tse = ones(145)*298
# ~ tsi = ones(145)*293
# ~ tse[0] = 293

# ---------------------------------------------
# ~ for i in range(145):
# ~ for j in range(nc):
# ~ qse[i] = qse[i]-a[j]*tse[i-j]+b[j]*tsi[i-j]
# ~ qsi[i] = qsi[i]+b[j]*tse[i-j]-c[j]*tsi[i-j]
# ~ for j in range(1,nd):
# ~ qse[i] = qse[i]-d[j]*qse[i-j]
# ~ qsi[i] = qsi[i]-d[j]*qsi[i-j]

# ~ if i < 5:
# ~ for j in range(5):
# ~ print(qsi[j])
# ~ print('---')
# ~ print(qsi)
# ---------------------------------------------

# ~ tempExt = zeros(ns)
# ~ tempInt = zeros(ns)

# ~ for i in range(145):
# ~ conHeatFluxExt = qse[i]
# ~ conHeatFluxInt = qsi[i]
# ~ for j in range(ns-1,-1,-1):
# ~ tempExt[j] = tse[i-j]
# ~ tempInt[j] = tsi[i-j]

# ~ for j in range(nc):
# ~ conHeatFluxExt = conHeatFluxExt-a[j]*tempExt[j]+b[j]*tempInt[j]
# ~ conHeatFluxInt = conHeatFluxInt+b[j]*tempExt[j]-c[j]*tempInt[j]
# ~ for j in range(1,nd):
# ~ conHeatFluxExt = conHeatFluxExt-d[j]*qse[j-1]
# ~ conHeatFluxInt = conHeatFluxInt-d[j]*qsi[j-1]
# ~ for j in range(145-1,0,-1):
# ~ qse[j] = qse[j-1]
# ~ qsi[j] = qsi[j-1]
# ~ qse[0] = conHeatFluxExt
# ~ qsi[0] = conHeatFluxInt
# ~ if i < 10:
# ~ for j in range(10):
# ~ print(qsi[j])
# ~ print('---')
# ~ print(qsi)
