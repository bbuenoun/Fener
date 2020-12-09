#! /usr/bin/python
# -*- coding: utf-8 -*-

from numpy import zeros, exp, log, atleast_2d, sqrt
import psychutils as psych

# ====================================================
#
# Fraunhofer-Institut für Solare Energiesysteme
# Thermal Systems and Buildings
# Heidenhofstr. 2, 79110 Freiburg, Germany
# http://www.ise.fraunhofer.de
#
# Author: Luis Camilo Londono
# Team Passive Systems and Daylighting
# luislondono@ise.fraunhofer.de
#
# ====================================================


# -----Extracts the area of the windows constructions for further calculations in Radiant temperature


def winAreas(winVect, frames, surfaces, tWindows, tSurfaces):

    WinAreas = zeros((winVect.shape[0], 4))
    idFrames = zeros([frames.shape[0]])
    surfAreas = zeros((surfaces.shape[0], 3))
    tPerWin = zeros((1, tWindows.shape[0]))
    tPerWall = zeros((1, tSurfaces.shape[0]))

    for i in range(frames.shape[0]):
        idFrames[i] = frames[i, 3]

    # -----Finds out to which frame  and wall belongs each window

    for i in range(winVect.shape[0]):
        WinAreas[i, 0] = winVect[i, 0] * winVect[i, 1]  # Area
        WinAreas[i, 1] = winVect[i, 2]  # frame id
        WinAreas[i, 2] = idFrames[int(WinAreas[i, 1])]  # surface id

    # -----Creates a vector with the area of each surface.

    for i in range(surfaces.shape[0]):
        surfAreas[i, 0] = surfaces[i, 0] * surfaces[i, 1]
        surfAreas[i, 1] = surfaces[i, 0] * surfaces[i, 1]

    # -----Calculates the area percentage of each window (WinAreas[i,3]), and substracts the windows areas from the wall it belongs (surfAreas[i,1])

    for i in range(winVect.shape[0]):

        if WinAreas[i, 2] == 0:
            WinAreas[i, 3] = WinAreas[i, 0] / surfAreas[0, 0]
            surfAreas[0, 1] = surfAreas[0, 1] - WinAreas[i, 0]
        elif WinAreas[i, 2] == 1:
            WinAreas[i, 3] = WinAreas[i, 0] / surfAreas[1, 0]
            surfAreas[1, 1] = surfAreas[1, 1] - WinAreas[i, 0]
        elif WinAreas[i, 2] == 2:
            WinAreas[i, 3] = WinAreas[i, 0] / surfAreas[2, 0]
            surfAreas[2, 1] = surfAreas[2, 1] - WinAreas[i, 0]
        elif WinAreas[i, 2] == 3:
            WinAreas[i, 3] = WinAreas[i, 0] / surfAreas[3, 0]
            surfAreas[3, 1] = surfAreas[3, 1] - WinAreas[i, 0]
        elif WinAreas[i, 2] == 4:
            WinAreas[i, 3] = WinAreas[i, 0] / surfAreas[4, 0]
            surfAreas[4, 1] = surfAreas[4, 1] - WinAreas[i, 0]
        elif WinAreas[i, 2] == 5:
            WinAreas[i, 3] = WinAreas[i, 0] / surfAreas[5, 0]
            surfAreas[5, 1] = surfAreas[5, 1] - WinAreas[i, 0]

    # -----calculates the area percentage of each wall after substract the windows areas-

    for i in range(surfaces.shape[0]):
        surfAreas[i, 2] = surfAreas[i, 1] / surfAreas[i, 0]

    # -----Multiplies the window temperature by the area percentage of the window. this value \
    # -----will be used to ponderate the final temperature of the surface

    for i in range(0, tWindows.shape[0]):
        tPerWin[:, i] = tWindows[i] * WinAreas[i, 3]

    # -----Multiplies the wall temperature by the area percentage of the wall. this value \
    # -----will be used to ponderate the final temperature of the surface (including the windows)

    for i in range(0, tSurfaces.shape[0]):
        tPerWall[:, i] = tSurfaces[i] * surfAreas[i, 2]

    # ------Add up the proportional temperatures of walls and windows to find a final surface temperature.

    for i in range(winVect.shape[0]):
        if WinAreas[i, 2] == 0:
            tPerWall[:, 0] = tPerWall[:, 0] + tPerWin[:, i]
        elif WinAreas[i, 2] == 1:
            tPerWall[:, 1] = tPerWall[:, 1] + tPerWin[:, i]
        elif WinAreas[i, 2] == 2:
            tPerWall[:, 2] = tPerWall[:, 2] + tPerWin[:, i]
        elif WinAreas[i, 2] == 3:
            tPerWall[:, 3] = tPerWall[:, 3] + tPerWin[:, i]
        elif WinAreas[i, 2] == 4:
            tPerWall[:, 4] = tPerWall[:, 4] + tPerWin[:, i]
        elif WinAreas[i, 2] == 5:
            tPerWall[:, 5] = tPerWall[:, 5] + tPerWin[:, i]

    return tPerWall


# ------------------------------------------------------------------------------------
# Mean radiant temperature calculation Simplified-----------
# ------------------------------------------------------------------------------------


def tempRadSimp0(tWalls, l, w, h):
    # -------Areas
    a0 = w * l
    a1 = h * l
    a2 = h * w
    a3 = h * l
    a4 = h * w
    a5 = w * l

    factNom = (a1 + a3) * (0.5 * l) + (a2 + a4) * (0.5 * w) + (a0 + a5) * (0.5 * h)

    # -------area factors
    f0 = (a0 * (h - 0.6)) / factNom
    f1 = (a1 * (0.5 * l)) / factNom
    f2 = (a2 * (0.5 * w)) / factNom
    f3 = (a3 * (0.5 * l)) / factNom
    f4 = (a4 * (0.5 * w)) / factNom
    f5 = (a5 * 0.6) / factNom

    # ------radiative temperature

    tRad = (
        f0 * tWalls[0]
        + f1 * tWalls[1]
        + f2 * tWalls[2]
        + f3 * tWalls[3]
        + f4 * tWalls[4]
        + f5 * tWalls[5]
    )

    return tRad


def tempRadSimp1(tWalls, l, w, h):
    # -------Areas
    a0 = w * l
    a1 = h * l
    a2 = h * w
    a3 = h * l
    a4 = h * w
    a5 = w * l

    factNom = (a1 + a3) * (0.5 * l) + (a2 + a4) * (0.5 * w) + (a0 + a5) * (0.5 * h)

    # -------area factors
    f0 = (a0 * (h - 0.6)) / factNom
    f1 = (a1 * (0.75 * l)) / factNom
    f2 = (a2 * (0.5 * w)) / factNom
    f3 = (a3 * (0.25 * l)) / factNom
    f4 = (a4 * (0.5 * w)) / factNom
    f5 = (a5 * 0.6) / factNom

    # ------radiative temperature

    tRad = (
        f0 * tWalls[0]
        + f1 * tWalls[1]
        + f2 * tWalls[2]
        + f3 * tWalls[3]
        + f4 * tWalls[4]
        + f5 * tWalls[5]
    )

    return tRad


def tempRadSimp2(tWalls, l, w, h):
    # -------Areas
    a0 = w * l
    a1 = h * l
    a2 = h * w
    a3 = h * l
    a4 = h * w
    a5 = w * l

    factNom = (a1 + a3) * (0.5 * l) + (a2 + a4) * (0.5 * w) + (a0 + a5) * (0.5 * h)

    # -------area factors
    f0 = (a0 * (h - 0.6)) / factNom
    f1 = (a1 * (0.5 * l)) / factNom
    f2 = (a2 * (0.75 * w)) / factNom
    f3 = (a3 * (0.5 * l)) / factNom
    f4 = (a4 * (0.25 * w)) / factNom
    f5 = (a5 * 0.6) / factNom

    # ------radiative temperature

    tRad = (
        f0 * tWalls[0]
        + f1 * tWalls[1]
        + f2 * tWalls[2]
        + f3 * tWalls[3]
        + f4 * tWalls[4]
        + f5 * tWalls[5]
    )

    return tRad


def tempRadSimp3(tWalls, l, w, h):
    # -------Areas
    a0 = w * l
    a1 = h * l
    a2 = h * w
    a3 = h * l
    a4 = h * w
    a5 = w * l

    factNom = (a1 + a3) * (0.5 * l) + (a2 + a4) * (0.5 * w) + (a0 + a5) * (0.5 * h)

    # -------area factors
    f0 = (a0 * (h - 0.6)) / factNom
    f1 = (a1 * (0.25 * l)) / factNom
    f2 = (a2 * (0.5 * w)) / factNom
    f3 = (a3 * (0.75 * l)) / factNom
    f4 = (a4 * (0.5 * w)) / factNom
    f5 = (a5 * 0.6) / factNom

    # ------radiative temperature

    tRad = (
        f0 * tWalls[0]
        + f1 * tWalls[1]
        + f2 * tWalls[2]
        + f3 * tWalls[3]
        + f4 * tWalls[4]
        + f5 * tWalls[5]
    )

    return tRad


def tempRadSimp4(tWalls, l, w, h):
    # -------Areas
    a0 = w * l
    a1 = h * l
    a2 = h * w
    a3 = h * l
    a4 = h * w
    a5 = w * l

    factNom = (a1 + a3) * (0.5 * l) + (a2 + a4) * (0.5 * w) + (a0 + a5) * (0.5 * h)

    # -------area factors
    f0 = (a0 * (h - 0.6)) / factNom
    f1 = (a1 * (0.5 * l)) / factNom
    f2 = (a2 * (0.25 * w)) / factNom
    f3 = (a3 * (0.5 * l)) / factNom
    f4 = (a4 * (0.75 * w)) / factNom
    f5 = (a5 * 0.6) / factNom

    # ------radiative temperature

    tRad = (
        f0 * tWalls[0]
        + f1 * tWalls[1]
        + f2 * tWalls[2]
        + f3 * tWalls[3]
        + f4 * tWalls[4]
        + f5 * tWalls[5]
    )

    return tRad


def tempRadSimp(twalls, l, w, h, ts, te, tn, tw, tf, tc):

    tr0 = tempRadSimp0(twalls, l - te - tw, w - ts - tn, h - tf - tc)
    tr1 = tempRadSimp1(twalls, l - te - tw, w - ts - tn, h - tf - tc)
    tr2 = tempRadSimp2(twalls, l - te - tw, w - ts - tn, h - tf - tc)
    tr3 = tempRadSimp3(twalls, l - te - tw, w - ts - tn, h - tf - tc)
    tr4 = tempRadSimp4(twalls, l - te - tw, w - ts - tn, h - tf - tc)

    return (tr0, tr1, tr2, tr3, tr4)


# ------------------------------------------------------------------------------------
# -----Relative humidity calculations -----------
# ------------------------------------------------------------------------------------

# -----Calculates the humidity ratio W using the Relative Humidity, the temperature and the atmospheric pressure (ASHRAE HANDBOOK, 2009)


def humFromRHumTemp(RH, T, P):
    # ----Saturation vapour pressure from ASHRAE
    C8 = -5.8002206e3
    C9 = 1.3914993
    C10 = -4.8640239e-2
    C11 = 4.1764768e-5
    C12 = -1.4452093e-8
    C13 = 6.5459673
    PWS = exp(C8 / T + C9 + C10 * T + C11 * T ** 2 + C12 * T ** 3 + C13 * log(T))
    # ----Vapour pressure
    PW = RH * PWS / 100
    # ----Specific humidity
    W = 0.62198 * PW / (P - PW)
    return W


# -----Calculates the density of moist air [Kg/m^3] based on temperature and relative Humidity (ASHRAE HANDBOOK, 2009)


def rho(temp, hr):
    pres = 101325
    hum = humFromRHumTemp(hr, temp, pres)
    rho = (pres / (1000 * 0.287042 * temp * (1.0 + 1.607858 * hum))) * (1 + hum)
    return rho


# -----Used for verification: calculates  the density of moist air [Kg/m^3] based on temperature and relative Humidity ( https://www.omnicalculator.com/physics/air-density)


def rho_test(temp, rh):
    temp = temp - 273.3
    satVap = 6.1078 * 10 ** (7.5 * temp / (temp + 237.3))
    pv = satVap * (rh / 100)
    pd = 1013.3 - pv
    rho = ((pd * 100) / (287.058 * (temp + 273.15))) + (
        (pv * 100) / (461.495 * (temp * 273.15))
    )
    return rho


# -----Calculates the internal relative humidity , based on meteorological data and FENER results, using (ASHRAE HANDBOOK, 2009) psychrometric equations.


def latEnergyBalance(
    volume,
    infSch,
    airExch,
    latIntHeatFlx,
    rhIn,
    tempInAir,
    relHumOut,
    tempOutAir,
    supAirMass,
    supAirWatFrac,
    timeStep,
):

    # -------------------------Initial Param
    lv = 2.26e6  # ----Specific latent heat of vaporization in J/Kg
    rho_vec = rho(tempInAir, rhIn)
    watFracIn = psych.X_pTphi(101325.0, tempInAir, rhIn)
    watFracOut = psych.X_pTphi(101325.0, tempOutAir, relHumOut / 100.0)
    # -------------------------Main equation
    deltaFrac_vec = watFracIn + (timeStep) / (rho_vec * lv * volume) * (
        latIntHeatFlx
        + infSch * airExch * volume / 3600.0 * rho_vec * lv * (watFracOut - watFracIn)
        + supAirMass * lv * (supAirWatFrac - watFracIn)
    )
    # ------------------------RH Calculation
    relHumIn = min(psych.phi_pTX(101325.0, tempInAir, deltaFrac_vec) * 100, 100.0)

    return relHumIn


# ------------------------------------------------------------------------------------
# -----pmv calculations-----------
# ------------------------------------------------------------------------------------


def pmvInitialize(
    temperatureAir,
    temperatureRad,
    relativeHumidity,
    clo,
    metaRate,
    externalWork,
    airVelocity,
):

    temp_in_k = temperatureAir + 273.15
    temp_rad_k = temperatureRad + 273.15
    temperatureAir = atleast_2d(temperatureAir)
    temperatureRad = atleast_2d(temperatureRad)
    # Return saturation pressure of water as a function of temperature T
    # in the range of 273.16 to 373.16 K.
    # Equation from:
    # Buildings.Media.PerfectGases.MoistAirUnsaturated.saturationPressureLiquidTetens
    # is much more accurate over the range from 0 to 50 °C
    # P=0.61078*exp[(17.2*T)/((T+237.4)], T:C and P:KPa """

    # ~ saturated_vapor_pressure=((611.657*exp((17.2799-4102.99/(temp_in_k-35.719))))/1000)

    saturated_vapor_pressure = (
        exp(16.6536 - 4030.183 / (temperatureAir + 235)) * 1000.0
    )  # Pa

    # Return water vapor pressure s as a function of:
    # relative humidity in % and saturation pressure of water in Kpa
    # Equation from:
    # ASHRAE Handbook-Fundamentals 2009
    # P=RH*saturated_vapor_pressure/100, P:KPa
    water_vapor_pressure = saturated_vapor_pressure * (relativeHumidity / 100.0)

    # thermal insulation of the clothing in m^2K/W
    icl = 0.155 * clo

    # metabolic rate, external work rate and internal heat production in W/m^2
    m = metaRate * 58.15
    w = externalWork * 58.15
    mw = m - w

    # clothing area factor:ASHRAE Handbook
    if icl <= 0.078:
        fcl = 1 + 1.29 * icl
    else:
        fcl = 1.05 + 0.645 * icl

    # heat trasnfer coeff by trasnfer convection
    hcf = 12.1 * sqrt(airVelocity)

    # Calculate surface temperaure of clothing by iteration

    tcla = (
        temp_in_k + (35.5 - temperatureAir) / (3.5 * (6.45 * icl + 0.1))
    ) - 273  # First attemp
    tcla = atleast_2d(tcla)
    a = 35.7 - 0.028 * mw  # float
    b = 39.6e-9 * fcl  # float
    xi = tcla
    xf = xi
    y = zeros([len(tcla)])
    hcn = zeros([len(tcla)])
    hc = zeros([len(tcla)])
    tcl = zeros([len(tcla)])

    for i in range(len(tcla)):

        xi = zeros([150])
        xf = zeros([150])
        c = zeros([150])

        d = zeros([150])
        xi[0] = tcla[i]
        xf[0] = tcla[i]

        for j in range(1, 150):
            xi[j] = (xi[j - 1] + xf[j - 1]) / 2
            hcn[i] = 2.38 * abs(xi[j] - temperatureAir[i]) ** 0.25

            if hcf > hcn[i]:

                hc[i] = hcf
            else:
                hc[i] = hcn[i]

            c[j] = (xi[j] + 273.0) ** 4.0 - (temperatureRad[i] + 273.0) ** 4.0

            d[j] = fcl * hc[i] * (xi[j] - temperatureAir[i])

            xf[j] = a - icl * (b * c[j] + d[j])

            if abs(xf[j] - xi[j]) < 0.000000001:

                y[i] = xf[j]
                break
        tcl[i] = y[i]

    # HEAT LOSS COMPONENTS -----------------------------------------------------

    # heat loss trough skin
    hl1 = 3.05 * 0.001 * (5733 - 6.99 * mw - (water_vapor_pressure))

    # Heat loss by sweating
    if mw > 58.15:
        hl2 = 0.42 * (mw - 58.15)
    else:
        hl2 = 0

    # latent respiration heat loss
    hl3 = 1.7e-5 * m * (5867 - (water_vapor_pressure))

    # Dry respiration heat loss
    hl4 = 0.0014 * m * (34 - temperatureAir)

    # heat loss by radiation
    hl5 = 3.96e-8 * fcl * ((tcl + 273.15) ** 4 - (temp_rad_k) ** 4)
    # ~ hl5=fcl*4.7*0.95*(tcl-temperatureRad)

    # Heat loss by convection
    hl6 = fcl * hc * (tcl - temperatureAir)

    # FINAL CALCULATION PMV AND PPD -------------------------------------------

    # thermal sensation coeff.
    ts = 0.303 * exp(-0.036 * m) + 0.028

    # PMV: Predicted mean vote
    pmv = ts * (mw - hl1 - hl2 - hl3 - hl4 - hl5 - hl6)

    # predited percentage dissat
    ppd = 100.0 - 95.0 * exp(-0.03353 * pmv ** 4 - 0.2179 * pmv ** 2)

    return (pmv, ppd)


def pmvSimp(
    temperatureAir,
    tRad,
    relativeHumidity,
    clo=1.0,
    metaRate=1.0,
    externalWork=0.0,
    airVelocity=0.1,
):

    tcIndicators0 = pmvInitialize(
        temperatureAir,
        tRad[0],
        relativeHumidity,
        clo,
        metaRate,
        externalWork,
        airVelocity,
    )
    tcIndicators1 = pmvInitialize(
        temperatureAir,
        tRad[1],
        relativeHumidity,
        clo,
        metaRate,
        externalWork,
        airVelocity,
    )
    tcIndicators2 = pmvInitialize(
        temperatureAir,
        tRad[2],
        relativeHumidity,
        clo,
        metaRate,
        externalWork,
        airVelocity,
    )
    tcIndicators3 = pmvInitialize(
        temperatureAir,
        tRad[3],
        relativeHumidity,
        clo,
        metaRate,
        externalWork,
        airVelocity,
    )
    tcIndicators4 = pmvInitialize(
        temperatureAir,
        tRad[4],
        relativeHumidity,
        clo,
        metaRate,
        externalWork,
        airVelocity,
    )

    return (
        tcIndicators0[0],
        tcIndicators1[0],
        tcIndicators2[0],
        tcIndicators3[0],
        tcIndicators4[0],
    )


def ppdSimp(
    temperatureAir,
    tRad,
    relativeHumidity,
    clo=1.0,
    metaRate=1.0,
    externalWork=0.0,
    airVelocity=0.1,
):

    tcIndicators0 = pmvInitialize(
        temperatureAir,
        tRad[0],
        relativeHumidity,
        clo,
        metaRate,
        externalWork,
        airVelocity,
    )
    tcIndicators1 = pmvInitialize(
        temperatureAir,
        tRad[1],
        relativeHumidity,
        clo,
        metaRate,
        externalWork,
        airVelocity,
    )
    tcIndicators2 = pmvInitialize(
        temperatureAir,
        tRad[2],
        relativeHumidity,
        clo,
        metaRate,
        externalWork,
        airVelocity,
    )
    tcIndicators3 = pmvInitialize(
        temperatureAir,
        tRad[3],
        relativeHumidity,
        clo,
        metaRate,
        externalWork,
        airVelocity,
    )
    tcIndicators4 = pmvInitialize(
        temperatureAir,
        tRad[4],
        relativeHumidity,
        clo,
        metaRate,
        externalWork,
        airVelocity,
    )

    return (
        tcIndicators0[1],
        tcIndicators1[1],
        tcIndicators2[1],
        tcIndicators3[1],
        tcIndicators4[1],
    )
