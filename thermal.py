#!/usr/bin/env python
# -*- coding: utf-8 -*-


from numpy import zeros, matrix
import param as p
import calcs
from scipy.linalg import solve


# ====================================================
#
# Fraunhofer-Institut für Solare Energiesysteme
# Thermal Systems and Buildings
# Heidenhofstr. 2, 79110 Freiburg, Germany
# http://www.ise.fraunhofer.de
#
# Author: Bruno Bueno
# Team Passive Systems and Daylighting
# Phone: +49(0)761 4588 5377
# bruno.bueno@ise.fraunhofer.de
#
# ====================================================


def genRHTC(
    numWin,
    numFrame,
    numSurf,
    vfWinFrame,
    vfFrameWin,
    vfWinSurf,
    vfSurfWin,
    vfFrameSurf,
    vfSurfFrame,
    vfWinWin,
    vfFrameFrame,
    vfSurfSurf,
    win,
    frame,
    surf,
):
    # Radiation exchange between two diffuse gray surfaces
    # linearization of the Stephan-Boltzman equation
    # secondary reflections are neglected (what is reflected from what is reflected)
    # q12 = F12*(eps1*Eb1 - alp1*G2), where Eb1 is black-body emission of surf1 and G2 is the radiation incident on surf1 from surface 2
    # absortivity = emissivity, alp = eps
    # q12 = F12*eps1*(Eb1 - eps2*Eb2 - (1-eps2)*Eb1) # the last term is what is reflected in surf2 from what is emitted from the whole enclosure,
    # assuming that the average temperature of all surfaces is similar to the temperature of surf1.
    # q12 = F12*eps1*eps2*(Eb1 - Eb2)
    RHTCwinFrame = zeros([numWin, numFrame])
    RHTCframeWin = zeros([numFrame, numWin])
    RHTCwinSurf = zeros([numWin, numSurf])
    RHTCsurfWin = zeros([numSurf, numWin])
    RHTCframeSurf = zeros([numFrame, numSurf])
    RHTCsurfFrame = zeros([numSurf, numFrame])
    RHTCwinWin = zeros([numWin, numWin])
    RHTCframeFrame = zeros([numFrame, numFrame])
    RHTCsurfSurf = zeros([numSurf, numSurf])
    for i in range(numWin):
        for j in range(numFrame):
            k = i * numFrame + j
            RHTCwinFrame[i, j] = (
                vfWinFrame[k]
                * win[i].intEmis
                * frame[j].intEmis
                * 4
                * p.sigma
                * (win[i].intTemp / 2.0 + frame[j].intTemp / 2.0) ** 3.0
            )
    for i in range(numFrame):
        for j in range(numWin):
            k = i * numWin + j
            RHTCframeWin[i, j] = (
                vfFrameWin[k]
                * frame[i].intEmis
                * win[j].intEmis
                * 4
                * p.sigma
                * (win[j].intTemp / 2.0 + frame[i].intTemp / 2.0) ** 3.0
            )
    for i in range(numWin):
        for j in range(numSurf):
            k = i * numSurf + j
            RHTCwinSurf[i, j] = (
                vfWinSurf[k]
                * win[i].intEmis
                * surf[j].intEmis
                * 4
                * p.sigma
                * (win[i].intTemp / 2.0 + surf[j].intTemp / 2.0) ** 3.0
            )
    for i in range(numSurf):
        for j in range(numWin):
            k = i * numWin + j
            RHTCsurfWin[i, j] = (
                vfSurfWin[k]
                * surf[i].intEmis
                * win[j].intEmis
                * 4
                * p.sigma
                * (win[j].intTemp / 2.0 + surf[i].intTemp / 2.0) ** 3.0
            )
    for i in range(numFrame):
        for j in range(numSurf):
            k = i * numSurf + j
            RHTCframeSurf[i, j] = (
                vfFrameSurf[k]
                * frame[i].intEmis
                * surf[j].intEmis
                * 4
                * p.sigma
                * (frame[i].intTemp / 2.0 + surf[j].intTemp / 2.0) ** 3.0
            )
    for i in range(numSurf):
        for j in range(numFrame):
            k = i * numFrame + j
            RHTCsurfFrame[i, j] = (
                vfSurfFrame[k]
                * surf[i].intEmis
                * frame[j].intEmis
                * 4
                * p.sigma
                * (frame[j].intTemp / 2.0 + surf[i].intTemp / 2.0) ** 3.0
            )
    for i in range(numWin):
        for j in range(numWin):
            k = i * numWin + j
            RHTCwinWin[i, j] = (
                vfWinWin[k]
                * win[i].intEmis
                * win[j].intEmis
                * 4
                * p.sigma
                * (win[i].intTemp / 2.0 + win[j].intTemp / 2.0) ** 3.0
            )
    for i in range(numFrame):
        for j in range(numFrame):
            k = i * numFrame + j
            RHTCframeFrame[i, j] = (
                vfFrameFrame[k]
                * frame[i].intEmis
                * frame[j].intEmis
                * 4
                * p.sigma
                * (frame[i].intTemp / 2.0 + frame[j].intTemp / 2.0) ** 3.0
            )
    for i in range(numSurf):
        for j in range(numSurf):
            k = i * numSurf + j
            RHTCsurfSurf[i, j] = (
                vfSurfSurf[k]
                * surf[i].intEmis
                * surf[j].intEmis
                * 4
                * p.sigma
                * (surf[i].intTemp / 2.0 + surf[j].intTemp / 2.0) ** 3.0
            )
    return (
        RHTCwinFrame,
        RHTCframeWin,
        RHTCwinSurf,
        RHTCsurfWin,
        RHTCframeSurf,
        RHTCsurfFrame,
        RHTCwinWin,
        RHTCframeFrame,
        RHTCsurfSurf,
    )


def genCHTC(numWin, numFrame, numSurf, tempInAir, win, frame, surf):
    CHTCwin = zeros(numWin)
    CHTCframe = zeros(numFrame)
    CHTCsurf = zeros(numSurf)
    for i in range(numWin):
        CHTCwin[i] = calcs.winCHTC(tempInAir, win[i].intTemp, win[i].area, win[i].tilt)
    for j in range(numFrame):
        CHTCframe[j] = calcs.surfCHTC(tempInAir, frame[j].intTemp, frame[j].tilt)
    for k in range(numSurf):
        CHTCsurf[k] = calcs.surfCHTC(tempInAir, surf[k].intTemp, surf[k].tilt)
    # ----------------------------------------------------
    return CHTCwin, CHTCframe, CHTCsurf


def genIntHeatGains(
    floorArea,
    powerLight,
    radFracLight,
    lightSch,
    powerEquip,
    radFracEquip,
    equipSch,
    powerPeople,
    radFracPeople,
    occSch,
):
    # ----------------------------------------------------
    # convective heat fluxes from internal heat sources to surfaces (W)
    convIntHeatFlx = (
        floorArea * powerLight * (1 - radFracLight) * lightSch
        + floorArea * powerEquip * (1 - radFracEquip) * equipSch
        + floorArea * powerPeople * (1 - radFracPeople) * occSch
    )
    # ----------------------------------------------------
    # radiant heat fluxes from internal heat sources to surfaces (W)
    radIntHeatFlx = (
        floorArea * powerLight * radFracLight * lightSch
        + floorArea * powerEquip * radFracEquip * equipSch
        + floorArea * powerPeople * radFracPeople * occSch
    )
    # ----------------------------------------------------
    return convIntHeatFlx, radIntHeatFlx


def bldEnergyDemand(
    tempInAir,
    numWin,
    numFrame,
    numSurf,
    areaWin,
    areaFrame,
    areaSurf,
    win,
    frame,
    surf,
    CHTCwin,
    CHTCframe,
    CHTCsurf,
    convIntHeatFlx,
    infSch,
    airExch,
    volume,
    tempOutAir,
):
    # -------------------------------------
    # building energy demand
    elementsContribution = 0.0
    for i in range(numWin):
        elementsContribution = elementsContribution + areaWin[i] * CHTCwin[i] * (
            win[i].intTemp - tempInAir
        )
    for i in range(numFrame):
        elementsContribution = elementsContribution + areaFrame[i] * CHTCframe[i] * (
            frame[i].intTemp - tempInAir
        )
    for i in range(numSurf):
        elementsContribution = elementsContribution + areaSurf[i] * CHTCsurf[i] * (
            surf[i].intTemp - tempInAir
        )
    energyDemand = (
        -elementsContribution
        - convIntHeatFlx
        - infSch
        * airExch
        * volume
        / 3600.0
        * calcs.rho(tempInAir)
        * p.cp
        * (tempOutAir - tempInAir)
    )
    # ----------------------------------------------------
    return energyDemand


def invert(nz, a, c):
    #  a(*,0) lower diagonal (Ai,i-1)
    #  a(*,1) principal diagonal (Ai,i)
    #  a(*,2) upper diagonal (Ai,i+1)
    # --------------------------------------------------------------------------
    x = zeros(nz)
    for i in range(nz - 2, -1, -1):
        c[i] = c[i] - a[i, 2] * c[i + 1] / a[i + 1, 1]
        a[i, 1] = a[i, 1] - a[i, 2] * a[i + 1, 0] / a[i + 1, 1]
    for i in range(1, nz):
        c[i] = c[i] - a[i, 0] * c[i - 1] / a[i - 1, 1]
    for i in range(nz):
        x[i] = c[i] / a[i, 1]
    # --------------------------------------------------------------------------
    return x


def bldEnergyBalanceTFConstTemp(
    numWin,
    numFrame,
    numSurf,
    win,
    frame,
    surf,
    radIntHeatFlx,
    tempInAir,
    tempOutAir,
    skyTemp,
    CHTCwin,
    areaWin,
    CHTCframe,
    areaFrame,
    CHTCsurf,
    areaSurf,
    irrWinAbsExt,
    irrWinAbsInt,
    RHTCwinWin,
    RHTCwinFrame,
    RHTCwinSurf,
    RHTCframeFrame,
    RHTCframeWin,
    RHTCframeSurf,
    RHTCsurfSurf,
    RHTCsurfWin,
    RHTCsurfFrame,
    extRHTC,
    extCHTC,
    irrFrameInt,
    irrFrameExt,
    irrSurfInt,
    irrSurfExt,
    bc,
    sumArea,
):
    # ----------------------------------------------------
    sizeMtx = numWin * 2 + numFrame * 2 + numSurf * 2
    a = matrix(zeros([sizeMtx, sizeMtx]))  # [W/K]
    b = matrix(zeros([sizeMtx, 1]))  # [W]
    # -------------------------------------
    # next numWin lines T_{win}, 1-indoor, 2-outdoor
    xcount = 0
    for i in range(numWin):
        # first column T_{in}
        b[xcount, 0] = (
            irrWinAbsInt[i] * areaWin[i]
            + radIntHeatFlx * areaWin[i] / sumArea
            + CHTCwin[i] * areaWin[i] * tempInAir
        )
        b[xcount + 1, 0] = (
            irrWinAbsExt[i] * areaWin[i]
            + areaWin[i] * win[i].svf * win[i].extEmis * extRHTC * skyTemp
            + areaWin[i] * extCHTC * tempOutAir
        )
        ycount = 0
        # next numWin columns T_{win}, 1-indoor, 2-outdoor
        for j in range(numWin):
            if i == j:
                a[xcount, ycount] = (
                    CHTCwin[i] * areaWin[i]
                    + areaWin[i] / win[i].res
                    + areaWin[i]
                    * (
                        sum(RHTCwinFrame[i, :])
                        + sum(RHTCwinSurf[i, :])
                        + sum(RHTCwinWin[i, :])
                    )
                )
                a[xcount, ycount + 1] = -areaWin[i] / win[i].res
                a[xcount + 1, ycount + 1] = (
                    areaWin[i] / win[i].res
                    + areaWin[i] * win[i].svf * win[i].extEmis * extRHTC
                    + areaWin[i] * extCHTC
                )
                a[xcount + 1, ycount] = -areaWin[i] / win[i].res
            else:
                a[xcount, ycount] = -RHTCwinWin[i, j] * areaWin[i]
            ycount = ycount + 2
        # next numFrame columns T_{frame}, 1-indoor, 2-outdoor
        for j in range(numFrame):
            a[xcount, ycount] = -RHTCwinFrame[i, j] * areaWin[i]
            ycount = ycount + 2
        # next numSurf columns
        for j in range(numSurf):
            a[xcount, ycount] = -RHTCwinSurf[i, j] * areaWin[i]
            ycount = ycount + 2
        xcount = xcount + 2
    # -------------------------------------
    # next numFrame lines T_{frame}, 1-indoor, 2-outdoor
    for i in range(numFrame):
        # first column T_{in}
        b[xcount, 0] = (
            irrFrameInt[i] * (1.0 - frame[i].intAlb) * areaFrame[i]
            + radIntHeatFlx * areaFrame[i] / sumArea
            + CHTCframe[i] * areaFrame[i] * tempInAir
        )
        b[xcount + 1, 0] = (
            irrFrameExt[i] * (1.0 - frame[i].extAlb) * areaFrame[i]
            + areaFrame[i] * frame[i].svf * frame[i].extEmis * extRHTC * skyTemp
            + areaFrame[i] * extCHTC * tempOutAir
        )
        ycount = 0
        # next numWin columns T_{win}, 1-indoor, 2-outdoor
        for j in range(numWin):
            a[xcount, ycount] = -RHTCframeWin[i, j] * areaFrame[i]
            ycount = ycount + 2
        # next numFrame columns T_{frame}, 1-indoor, 2-outdoor
        for j in range(numFrame):
            if i == j:
                a[xcount, ycount] = (
                    CHTCframe[i] * areaFrame[i]
                    + areaFrame[i] / frame[i].res
                    + areaFrame[i]
                    * (
                        sum(RHTCframeWin[i, :])
                        + sum(RHTCframeSurf[i, :])
                        + sum(RHTCframeFrame[i, :])
                    )
                )
                a[xcount, ycount + 1] = -areaFrame[i] / frame[i].res
                a[xcount + 1, ycount + 1] = (
                    areaFrame[i] / frame[i].res
                    + areaFrame[i] * frame[i].svf * frame[i].extEmis * extRHTC
                    + areaFrame[i] * extCHTC
                )
                a[xcount + 1, ycount] = -areaFrame[i] / frame[i].res
            else:
                a[xcount, ycount] = -RHTCframeFrame[i, j] * areaFrame[i]
            ycount = ycount + 2
        # next numSurf columns
        for j in range(numSurf):
            a[xcount, ycount] = -RHTCframeSurf[i, j] * areaFrame[i]
            ycount = ycount + 2
        xcount = xcount + 2
    # -------------------------------------
    # next numSurf lines
    for i in range(numSurf):
        # first column T_{in}
        b[xcount, 0] = (
            irrSurfInt[i] * areaSurf[i]
            + radIntHeatFlx * areaSurf[i] / sumArea
            + CHTCsurf[i] * areaSurf[i] * tempInAir
        )
        for j in range(1, surf[i].nc):
            b[xcount, 0] = (
                b[xcount, 0]
                + surf[i].b[j] * surf[i].tse[j - 1] * areaSurf[i]
                - surf[i].c[j] * surf[i].tsi[j - 1] * areaSurf[i]
            )
        for j in range(1, surf[i].nd):
            b[xcount, 0] = (
                b[xcount, 0] - surf[i].d[j] * surf[i].qsi[j - 1] * areaSurf[i]
            )
        if bc[i] == 1:
            b[xcount + 1, 0] = (
                irrSurfExt[i] * (1.0 - surf[i].extAlb) * areaSurf[i]
                + areaSurf[i] * surf[i].svf * surf[i].extEmis * extRHTC * skyTemp
                + areaSurf[i] * extCHTC * tempOutAir
            )
            for j in range(1, surf[i].nc):
                b[xcount + 1, 0] = (
                    b[xcount + 1, 0]
                    - surf[i].a[j] * surf[i].tse[j - 1] * areaSurf[i]
                    + surf[i].b[j] * surf[i].tsi[j - 1] * areaSurf[i]
                )
            for j in range(1, surf[i].nd):
                b[xcount + 1, 0] = (
                    b[xcount + 1, 0] - surf[i].d[j] * surf[i].qse[j - 1] * areaSurf[i]
                )
        else:
            a[xcount + 1, 0] = a[xcount, 0]
            b[xcount + 1, 0] = b[xcount, 0]
        ycount = 0
        # next numWin columns T_{win}, 1-indoor, 2-outdoor
        for j in range(numWin):
            a[xcount, ycount] = -RHTCsurfWin[i, j] * areaSurf[i]
            if bc[i] < 1:
                a[xcount + 1, ycount] = a[xcount, ycount]
            ycount = ycount + 2
        # next numFrame columns T_{frame}, 1-indoor, 2-outdoor
        for j in range(numFrame):
            a[xcount, ycount] = -RHTCsurfFrame[i, j] * areaSurf[i]
            if bc[i] < 1:
                a[xcount + 1, ycount] = a[xcount, ycount]
            ycount = ycount + 2
        # next numSurf columns T_{surf},
        for j in range(numSurf):
            if i == j:
                a[xcount, ycount] = (
                    CHTCsurf[i] * areaSurf[i]
                    + areaSurf[i]
                    * (
                        sum(RHTCsurfWin[i, :])
                        + sum(RHTCsurfSurf[i, :])
                        + sum(RHTCsurfFrame[i, :])
                    )
                    + areaSurf[i] * surf[i].c[0]
                )
                a[xcount, ycount + 1] = -areaSurf[i] * surf[i].b[0]
                if bc[i] == 1:
                    a[xcount + 1, ycount + 1] = (
                        areaSurf[i] * surf[i].svf * surf[i].extEmis * extRHTC
                        + areaSurf[i] * extCHTC
                        + areaSurf[i] * surf[i].a[0]
                    )
                else:
                    a[xcount + 1, ycount + 1] = (
                        CHTCsurf[i] * areaSurf[i]
                        + areaSurf[i]
                        * (
                            sum(RHTCsurfWin[i, :])
                            + sum(RHTCsurfSurf[i, :])
                            + sum(RHTCsurfFrame[i, :])
                        )
                        + areaSurf[i] * surf[i].a[0]
                    )
                a[xcount + 1, ycount] = -areaSurf[i] * surf[i].b[0]
            else:
                a[xcount, ycount] = -RHTCsurfSurf[i, j] * areaSurf[i]
                if bc[i] < 1:
                    a[xcount + 1, ycount] = a[xcount, ycount]
            ycount = ycount + 2
        xcount = xcount + 2
    # -------------------------------------
    x = solve(a, b)
    # ----------------------------------------------------
    return x


def update(
    x,
    ycount,
    numWin,
    numFrame,
    numSurf,
    win,
    frame,
    surf,
    areaWin,
    areaFrame,
    areaSurf,
    sumArea,
    calorim,
    conIndex,
    CHTCwin,
    CHTCsurf,
    tempInAir,
    radIntHeatFlx,
    RHTCsurfWin,
    RHTCsurfFrame,
    RHTCsurfSurf,
    RHTCwinWin,
    RHTCwinFrame,
    RHTCwinSurf,
):
    # -------------------------------------
    for i in range(numWin):
        # next numWin columns T_{win}, 1-indoor, 2-outdoor
        win[i].intTemp = x[ycount]
        win[i].extTemp = x[ycount + 1]
        ycount = ycount + 2
    for i in range(numFrame):
        # next numFrame columns T_{frame}, 1-indoor, 2-outdoor
        frame[i].intTemp = x[ycount]
        frame[i].extTemp = x[ycount + 1]
        ycount = ycount + 2
    for i in range(numSurf):
        # next numSurf columns
        surf[i].intTemp = x[ycount]
        surf[i].extTemp = x[ycount + 1]
        ycount = ycount + 2
    # -------------------------------------
    # update surface temperatures
    for i in range(numSurf):
        for j in range(surf[i].ns - 1, 0, -1):
            surf[i].tse[j] = surf[i].tse[j - 1]
            surf[i].tsi[j] = surf[i].tsi[j - 1]
        surf[i].tse[0] = surf[i].extTemp
        surf[i].tsi[0] = surf[i].intTemp
    # -------------------------------
    # update window mean temperature and virtual resistance
    for i in range(numWin):
        tMean = (win[i].intTemp + win[i].extTemp) / 2.0
        win[i].res = (
            1.0 / (calorim[conIndex[i]] * (-0.8395 + 0.0065 * tMean)) - 0.04 - 0.13
        )
    # ----------------------------------------------------
    # calculate conduction heat fluxes
    conHeatFluxInt = zeros(numSurf)
    conHeatFluxExt = zeros(numSurf)
    for i in range(numSurf):
        for j in range(surf[i].nc):
            conHeatFluxExt[i] = (
                conHeatFluxExt[i]
                - surf[i].a[j] * surf[i].tse[j]
                + surf[i].b[j] * surf[i].tsi[j]
            )
            conHeatFluxInt[i] = (
                conHeatFluxInt[i]
                + surf[i].b[j] * surf[i].tse[j]
                - surf[i].c[j] * surf[i].tsi[j]
            )
        for j in range(1, surf[i].nd):
            conHeatFluxExt[i] = conHeatFluxExt[i] - surf[i].d[j] * surf[i].qse[j - 1]
            conHeatFluxInt[i] = conHeatFluxInt[i] - surf[i].d[j] * surf[i].qsi[j - 1]
        for j in range(surf[i].ns - 1, 0, -1):
            surf[i].qse[j] = surf[i].qse[j - 1]
            surf[i].qsi[j] = surf[i].qsi[j - 1]
        surf[i].qse[0] = conHeatFluxExt[i]
        surf[i].qsi[0] = conHeatFluxInt[i]
    # ----------------------------------------------------
    # output for verification
    infrSurfInt = zeros(numSurf)
    convSurfInt = zeros(numSurf)
    heatSurfInt = zeros(numSurf)
    for i in range(numSurf):
        convSurfInt[i] = areaSurf[i] * CHTCsurf[i] * (tempInAir - surf[i].intTemp)
        heatSurfInt[i] = radIntHeatFlx * areaSurf[i] / sumArea
        for j in range(numWin):
            for k in range(numFrame):
                for l in range(numSurf):
                    infrSurfInt[i] = infrSurfInt[i] + areaSurf[i] * (
                        RHTCsurfWin[i, j] * (win[j].intTemp - surf[i].intTemp)
                        + RHTCsurfFrame[i, k] * (frame[k].intTemp - surf[i].intTemp)
                        + RHTCsurfSurf[i, l] * (surf[l].intTemp - surf[i].intTemp)
                    )
    # ---------------------------
    infrWinInt = zeros(numWin)
    convWinInt = zeros(numWin)
    heatWinInt = zeros(numWin)
    for i in range(numWin):
        convWinInt[i] = areaWin[i] * CHTCwin[i] * (tempInAir - win[i].intTemp)
        heatWinInt[i] = radIntHeatFlx * areaWin[i] / sumArea
        for j in range(numWin):
            for k in range(numFrame):
                for l in range(numSurf):
                    infrWinInt[i] = infrWinInt[i] + areaWin[i] * (
                        RHTCwinWin[i, j] * (win[j].intTemp - win[i].intTemp)
                        + RHTCwinFrame[i, k] * (frame[k].intTemp - win[i].intTemp)
                        + RHTCwinSurf[i, l] * (surf[l].intTemp - win[i].intTemp)
                    )
    # -------------------------------------
    return (
        win,
        frame,
        surf,
        conHeatFluxInt,
        conHeatFluxExt,
        infrSurfInt,
        convSurfInt,
        heatSurfInt,
        infrWinInt,
        convWinInt,
        heatWinInt,
    )


def bldEnergyBalanceTF(
    numWin,
    numFrame,
    numSurf,
    win,
    frame,
    surf,
    volume,
    infSch,
    airExch,
    convIntHeatFlx,
    radIntHeatFlx,
    tempInAir,
    tempOutAir,
    skyTemp,
    CHTCwin,
    areaWin,
    CHTCframe,
    areaFrame,
    CHTCsurf,
    areaSurf,
    irrWinAbsExt,
    irrWinAbsInt,
    RHTCwinWin,
    RHTCwinFrame,
    RHTCwinSurf,
    RHTCframeFrame,
    RHTCframeWin,
    RHTCframeSurf,
    RHTCsurfSurf,
    RHTCsurfWin,
    RHTCsurfFrame,
    extRHTC,
    extCHTC,
    irrFrameInt,
    irrFrameExt,
    irrSurfInt,
    irrSurfExt,
    bc,
    sumArea,
):
    # ----------------------------------------------------
    # air density
    rho = calcs.rho(tempInAir)
    # ----------------------------------------------------
    sizeMtx = 1 + numWin * 2 + numFrame * 2 + numSurf * 2
    a = matrix(zeros([sizeMtx, sizeMtx]))  # [W/K]
    b = matrix(zeros([sizeMtx, 1]))  # [W]
    # -------------------------------------
    # first line T_{in}
    a[0, 0] = (
        rho * p.cp * volume / p.timeStep
        + infSch * airExch * volume / 3600.0 * rho * p.cp
        + sum(CHTCwin * areaWin)
        + sum(CHTCframe * areaFrame)
        + sum(CHTCsurf * areaSurf)
    )
    b[0, 0] = (
        convIntHeatFlx
        + rho * p.cp * volume / p.timeStep * tempInAir
        + infSch * airExch * volume / 3600.0 * rho * p.cp * tempOutAir
    )
    ycount = 1
    for i in range(numWin):
        # next numWin columns T_{win}, 1-indoor, 2-outdoor
        a[0, ycount] = -CHTCwin[i] * areaWin[i]
        ycount = ycount + 2
    for i in range(numFrame):
        # next numFrame columns T_{frame}, 1-indoor, 2-outdoor
        a[0, ycount] = -CHTCframe[i] * areaFrame[i]
        ycount = ycount + 2
    for i in range(numSurf):
        # next numSurf columns
        a[0, ycount] = -CHTCsurf[i] * areaSurf[i]
        ycount = ycount + 2
    # -------------------------------------
    # next numWin lines T_{win}, 1-indoor, 2-outdoor
    xcount = 1
    for i in range(numWin):
        # first column T_{in}
        a[xcount, 0] = -CHTCwin[i] * areaWin[i]
        b[xcount, 0] = (
            irrWinAbsInt[i] * areaWin[i] + radIntHeatFlx * areaWin[i] / sumArea
        )
        b[xcount + 1, 0] = (
            irrWinAbsExt[i] * areaWin[i]
            + areaWin[i] * win[i].svf * win[i].extEmis * extRHTC * skyTemp
            + areaWin[i] * extCHTC * tempOutAir
        )
        ycount = 1
        # next numWin columns T_{win}, 1-indoor, 2-outdoor
        for j in range(numWin):
            if i == j:
                a[xcount, ycount] = (
                    CHTCwin[i] * areaWin[i]
                    + areaWin[i] / win[i].res
                    + areaWin[i]
                    * (
                        sum(RHTCwinFrame[i, :])
                        + sum(RHTCwinSurf[i, :])
                        + sum(RHTCwinWin[i, :])
                    )
                )
                a[xcount, ycount + 1] = -areaWin[i] / win[i].res
                a[xcount + 1, ycount + 1] = (
                    areaWin[i] / win[i].res
                    + areaWin[i] * win[i].svf * win[i].extEmis * extRHTC
                    + areaWin[i] * extCHTC
                )
                a[xcount + 1, ycount] = -areaWin[i] / win[i].res
            else:
                a[xcount, ycount] = -RHTCwinWin[i, j] * areaWin[i]
            ycount = ycount + 2
        # next numFrame columns T_{frame}, 1-indoor, 2-outdoor
        for j in range(numFrame):
            a[xcount, ycount] = -RHTCwinFrame[i, j] * areaWin[i]
            ycount = ycount + 2
        # next numSurf columns
        for j in range(numSurf):
            a[xcount, ycount] = -RHTCwinSurf[i, j] * areaWin[i]
            ycount = ycount + 2
        xcount = xcount + 2
    # -------------------------------------
    # next numFrame lines T_{frame}, 1-indoor, 2-outdoor
    for i in range(numFrame):
        # first column T_{in}
        a[xcount, 0] = -CHTCframe[i] * areaFrame[i]
        b[xcount, 0] = (
            irrFrameInt[i] * (1.0 - frame[i].intAlb) * areaFrame[i]
            + radIntHeatFlx * areaFrame[i] / sumArea
        )
        b[xcount + 1, 0] = (
            irrFrameExt[i] * (1.0 - frame[i].extAlb) * areaFrame[i]
            + areaFrame[i] * frame[i].svf * frame[i].extEmis * extRHTC * skyTemp
            + areaFrame[i] * extCHTC * tempOutAir
        )
        ycount = 1
        # next numWin columns T_{win}, 1-indoor, 2-outdoor
        for j in range(numWin):
            a[xcount, ycount] = -RHTCframeWin[i, j] * areaFrame[i]
            ycount = ycount + 2
        # next numFrame columns T_{frame}, 1-indoor, 2-outdoor
        for j in range(numFrame):
            if i == j:
                a[xcount, ycount] = (
                    CHTCframe[i] * areaFrame[i]
                    + areaFrame[i] / frame[i].res
                    + areaFrame[i]
                    * (
                        sum(RHTCframeWin[i, :])
                        + sum(RHTCframeSurf[i, :])
                        + sum(RHTCframeFrame[i, :])
                    )
                )
                a[xcount, ycount + 1] = -areaFrame[i] / frame[i].res
                a[xcount + 1, ycount + 1] = (
                    areaFrame[i] / frame[i].res
                    + areaFrame[i] * frame[i].svf * frame[i].extEmis * extRHTC
                    + areaFrame[i] * extCHTC
                )
                a[xcount + 1, ycount] = -areaFrame[i] / frame[i].res
            else:
                a[xcount, ycount] = -RHTCframeFrame[i, j] * areaFrame[i]
            ycount = ycount + 2
        # next numSurf columns
        for j in range(numSurf):
            a[xcount, ycount] = -RHTCframeSurf[i, j] * areaFrame[i]
            ycount = ycount + 2
        xcount = xcount + 2
    # -------------------------------------
    # next numSurf lines
    for i in range(numSurf):
        # first column T_{in}
        a[xcount, 0] = -CHTCsurf[i] * areaSurf[i]
        b[xcount, 0] = (
            irrSurfInt[i] * areaSurf[i] + radIntHeatFlx * areaSurf[i] / sumArea
        )
        for j in range(1, surf[i].nc):
            b[xcount, 0] = (
                b[xcount, 0]
                + surf[i].b[j] * surf[i].tse[j - 1] * areaSurf[i]
                - surf[i].c[j] * surf[i].tsi[j - 1] * areaSurf[i]
            )
        for j in range(1, surf[i].nd):
            b[xcount, 0] = (
                b[xcount, 0] - surf[i].d[j] * surf[i].qsi[j - 1] * areaSurf[i]
            )
        if bc[i] == 1:
            b[xcount + 1, 0] = (
                irrSurfExt[i] * (1.0 - surf[i].extAlb) * areaSurf[i]
                + areaSurf[i] * surf[i].svf * surf[i].extEmis * extRHTC * skyTemp
                + areaSurf[i] * extCHTC * tempOutAir
            )
            for j in range(1, surf[i].nc):
                b[xcount + 1, 0] = (
                    b[xcount + 1, 0]
                    - surf[i].a[j] * surf[i].tse[j - 1] * areaSurf[i]
                    + surf[i].b[j] * surf[i].tsi[j - 1] * areaSurf[i]
                )
            for j in range(1, surf[i].nd):
                b[xcount + 1, 0] = (
                    b[xcount + 1, 0] - surf[i].d[j] * surf[i].qse[j - 1] * areaSurf[i]
                )
        else:
            a[xcount + 1, 0] = a[xcount, 0]
            b[xcount + 1, 0] = b[xcount, 0]
        ycount = 1
        # next numWin columns T_{win}, 1-indoor, 2-outdoor
        for j in range(numWin):
            a[xcount, ycount] = -RHTCsurfWin[i, j] * areaSurf[i]
            if bc[i] < 1:
                a[xcount + 1, ycount] = a[xcount, ycount]
            ycount = ycount + 2
        # next numFrame columns T_{frame}, 1-indoor, 2-outdoor
        for j in range(numFrame):
            a[xcount, ycount] = -RHTCsurfFrame[i, j] * areaSurf[i]
            if bc[i] < 1:
                a[xcount + 1, ycount] = a[xcount, ycount]
            ycount = ycount + 2
        # next numSurf columns T_{surf},
        for j in range(numSurf):
            if i == j:
                a[xcount, ycount] = (
                    CHTCsurf[i] * areaSurf[i]
                    + areaSurf[i]
                    * (
                        sum(RHTCsurfWin[i, :])
                        + sum(RHTCsurfSurf[i, :])
                        + sum(RHTCsurfFrame[i, :])
                    )
                    + areaSurf[i] * surf[i].c[0]
                )
                a[xcount, ycount + 1] = -areaSurf[i] * surf[i].b[0]
                if bc[i] == 1:
                    a[xcount + 1, ycount + 1] = (
                        areaSurf[i] * surf[i].svf * surf[i].extEmis * extRHTC
                        + areaSurf[i] * extCHTC
                        + areaSurf[i] * surf[i].a[0]
                    )
                else:
                    a[xcount + 1, ycount + 1] = (
                        CHTCsurf[i] * areaSurf[i]
                        + areaSurf[i]
                        * (
                            sum(RHTCsurfWin[i, :])
                            + sum(RHTCsurfSurf[i, :])
                            + sum(RHTCsurfFrame[i, :])
                        )
                        + areaSurf[i] * surf[i].a[0]
                    )
                a[xcount + 1, ycount] = -areaSurf[i] * surf[i].b[0]
            else:
                a[xcount, ycount] = -RHTCsurfSurf[i, j] * areaSurf[i]
                if bc[i] < 1:
                    a[xcount + 1, ycount] = a[xcount, ycount]
            ycount = ycount + 2
        xcount = xcount + 2
    # -------------------------------------
    x = solve(a, b)
    # -------------------------------------
    return x


def main(
    numWin,
    numFrame,
    numSurf,
    win,
    frame,
    surf,
    volume,
    infSch,
    airExch,
    convIntHeatFlx,
    radIntHeatFlx,
    tempInAir,
    tempOutAir,
    skyTemp,
    CHTCwin,
    areaWin,
    CHTCframe,
    areaFrame,
    CHTCsurf,
    areaSurf,
    irrWinAbsExt,
    irrWinAbsInt,
    RHTCwinWin,
    RHTCwinFrame,
    RHTCwinSurf,
    RHTCframeFrame,
    RHTCframeWin,
    RHTCframeSurf,
    RHTCsurfSurf,
    RHTCsurfWin,
    RHTCsurfFrame,
    extRHTC,
    extCHTC,
    irrFrameInt,
    irrFrameExt,
    irrSurfInt,
    irrSurfExt,
    bc,
    calorim,
    conIndex,
    heatSpSch,
    coolSpSch,
):
    # -------------------------------------
    sumArea = sum(areaWin) + sum(areaFrame) + sum(areaSurf)
    # -------------------------------------
    x = bldEnergyBalanceTF(
        numWin,
        numFrame,
        numSurf,
        win,
        frame,
        surf,
        volume,
        infSch,
        airExch,
        convIntHeatFlx,
        radIntHeatFlx,
        tempInAir,
        tempOutAir,
        skyTemp,
        CHTCwin,
        areaWin,
        CHTCframe,
        areaFrame,
        CHTCsurf,
        areaSurf,
        irrWinAbsExt,
        irrWinAbsInt,
        RHTCwinWin,
        RHTCwinFrame,
        RHTCwinSurf,
        RHTCframeFrame,
        RHTCframeWin,
        RHTCframeSurf,
        RHTCsurfSurf,
        RHTCsurfWin,
        RHTCsurfFrame,
        extRHTC,
        extCHTC,
        irrFrameInt,
        irrFrameExt,
        irrSurfInt,
        irrSurfExt,
        bc,
        sumArea,
    )
    # -------------------------------------
    tempInAir = x[0]
    ycount = 1
    if heatSpSch - tempInAir > 0.0:
        tempInAir = heatSpSch
        energyDemandFlag = True
    elif tempInAir - coolSpSch > 0.0:
        tempInAir = coolSpSch
        energyDemandFlag = True
    else:
        energyDemandFlag = False
    # -------------------------------------
    if energyDemandFlag:
        x = bldEnergyBalanceTFConstTemp(
            numWin,
            numFrame,
            numSurf,
            win,
            frame,
            surf,
            radIntHeatFlx,
            tempInAir,
            tempOutAir,
            skyTemp,
            CHTCwin,
            areaWin,
            CHTCframe,
            areaFrame,
            CHTCsurf,
            areaSurf,
            irrWinAbsExt,
            irrWinAbsInt,
            RHTCwinWin,
            RHTCwinFrame,
            RHTCwinSurf,
            RHTCframeFrame,
            RHTCframeWin,
            RHTCframeSurf,
            RHTCsurfSurf,
            RHTCsurfWin,
            RHTCsurfFrame,
            extRHTC,
            extCHTC,
            irrFrameInt,
            irrFrameExt,
            irrSurfInt,
            irrSurfExt,
            bc,
            sumArea,
        )
        ycount = 0
    # -------------------------------------
    (
        win,
        frame,
        surf,
        conHeatFluxInt,
        conHeatFluxExt,
        infrSurfInt,
        convSurfInt,
        heatSurfInt,
        infrWinInt,
        convWinInt,
        heatWinInt,
    ) = update(
        x,
        ycount,
        numWin,
        numFrame,
        numSurf,
        win,
        frame,
        surf,
        areaWin,
        areaFrame,
        areaSurf,
        sumArea,
        calorim,
        conIndex,
        CHTCwin,
        CHTCsurf,
        tempInAir,
        radIntHeatFlx,
        RHTCsurfWin,
        RHTCsurfFrame,
        RHTCsurfSurf,
        RHTCwinWin,
        RHTCwinFrame,
        RHTCwinSurf,
    )
    # -------------------------------------
    if energyDemandFlag:
        energyDemand = bldEnergyDemand(
            tempInAir,
            numWin,
            numFrame,
            numSurf,
            areaWin,
            areaFrame,
            areaSurf,
            win,
            frame,
            surf,
            CHTCwin,
            CHTCframe,
            CHTCsurf,
            convIntHeatFlx,
            infSch,
            airExch,
            volume,
            tempOutAir,
        )
    else:
        energyDemand = 0.0
    # -------------------------------------
    return (
        energyDemand,
        tempInAir,
        win,
        frame,
        surf,
        x,
        infrSurfInt,
        convSurfInt,
        heatSurfInt,
        infrWinInt,
        convWinInt,
        heatWinInt,
    )
