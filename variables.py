#!/usr/bin/env python
# -*- coding: utf-8 -*-

from numpy import atleast_1d, genfromtxt, zeros, arange, atleast_2d
import ioRoutines
import calcs


def conWinWithoutGap(conWinVect):
    maxNumPane = 4
    conWin = zeros(int(2 * conWinVect[0] + 1))
    conWin[0] = conWinVect[0]
    count = 1
    for k in range(2 * maxNumPane - 1):
        if conWinVect[2 * k + 1] > 0:
            conWin[count] = conWinVect[2 * k + 1]
            conWin[count + 1] = conWinVect[2 * k + 2]
            count = count + 2
    return conWin


def dmxExt(workDir):
    print("Exterior illuminance daylight matrix")
    illExtDmx = zeros([1, 2306])
    illExtDmx[0, :] = genfromtxt(
        "%sillExt.dmx" % workDir, usecols=arange(0, 2306 * 3, 3)
    )
    # ----------------------------------------------------
    return illExtDmx


def dmxVis(numWin, workDir):
    print("Daylight matrix")
    dmx = zeros([numWin, 145, 2306])
    for i in range(numWin):
        dmx[i, :, :] = genfromtxt(
            "%swin_%i.dmx" % (workDir, i), usecols=arange(0, 2306 * 3, 3)
        )
    # ----------------------------------------------------
    return dmx


def vmxIll(numWin, numSensor, workDir):
    print("Illuminance view matrix")
    illVmx = zeros([numWin, numSensor, 145])
    for i in range(numWin):
        illVmx[i, :, :] = genfromtxt(
            "%sill_%i.vmx" % (workDir, i), usecols=arange(0, 145 * 3, 3)
        )
    # ----------------------------------------------------
    return illVmx


def dmxVisH(numWin, workDir):
    print("Daylight matrix")
    dmx = zeros([numWin, 145, 2306])
    for i in range(numWin):
        dmx[i, :, :] = genfromtxt(
            "%swin_%i.dmx" % (workDir, i),
            skip_header=11,
            usecols=arange(0, 2306 * 3, 3),
        )
    # ----------------------------------------------------
    return dmx


def vmxIllH(numWin, numSensor, workDir):
    print("Illuminance view matrix")
    illVmx = zeros([numWin, numSensor, 145])
    for i in range(numWin):
        illVmx[i, :, :] = genfromtxt(
            "%sill_%i.vmx" % (workDir, i), skip_header=10, usecols=arange(0, 145 * 3, 3)
        )
    # ----------------------------------------------------
    return illVmx


def tmxVis(numConWin, bsdfVis_path):
    print("Visual transmission matrix")
    propagationMtx = calcs.klemsIni()[0]
    tmxv = zeros([numConWin, 145, 145])
    for j in range(numConWin):
        mtx = ioRoutines.xml(bsdfVis_path[j], "Visible", "Transmission Front")
        tmxv[j, :, :] = propagationMtx.dot(mtx)
    # ----------------------------------------------------
    return tmxv


def calorim(numConWin, calorim_path):
    calorim = zeros([numConWin, 146])
    for j in range(numConWin):
        calorim[j, :] = genfromtxt(calorim_path[j], skip_header=1, delimiter=",")
    # ----------------------------------------------------
    return calorim


def tmxSol(numConWin, bsdfSol_path):
    print("Solar transmission matrix")
    propagationMtx = calcs.klemsIni()[0]
    tmxs = zeros([numConWin, 145, 145])
    rmxs = zeros([numConWin, 145, 145])
    for j in range(numConWin):
        mtx = ioRoutines.xml(bsdfSol_path[j], "Solar", "Transmission Front")
        tmxs[j, :, :] = propagationMtx * mtx
        mtx = ioRoutines.xml(bsdfSol_path[j], "Solar", "Reflection Front")
        rmxs[j, :, :] = propagationMtx * mtx
    # ----------------------------------------------------
    return tmxs, rmxs


def smxVis(meteo_path, iniHour, totHour, dayHour, sunAltitude):
    print("Visual sky matrix")
    path, name = ioRoutines.findPathName(meteo_path, "/", ".epw")
    visSmx = ioRoutines.skyMtx(
        "%s%sVis.smx" % (path, name), iniHour, totHour, dayHour, sunAltitude
    )
    # ----------------------------------------------------
    return visSmx


def vmxGlare(numWin, numSensorVert, workDir):
    # ----------------------------------------------------
    illVertVmx = zeros([numWin, numSensorVert, 145])
    for i in range(numWin):
        illVertVmx[i, :, :] = genfromtxt(
            "%sillVert_%i.vmx" % (workDir, i), usecols=arange(0, 145 * 3, 3)
        )
    # ----------------------------------------------------
    return illVertVmx


def winBack(numWin, workDir):
    numSensorWin = atleast_1d(genfromtxt("%snumSensorWin.dat" % workDir, dtype="int"))
    irrWinBackVmx = zeros([numWin, max(numSensorWin), 145])
    for i in range(numWin):
        aux_irrWinBackVmx = genfromtxt(
            "%sirrWinBack_%i.vmx" % (workDir, i), usecols=arange(0, 145 * 3, 3)
        )
        for j in range(numSensorWin[i]):
            irrWinBackVmx[i, j, :] = aux_irrWinBackVmx[j, :]
    return irrWinBackVmx, numSensorWin


def smxSol(meteo_path, iniHour, totHour, dayHour, sunAltitude):
    # ----------------------------------------------------
    #        start_time = time.time()
    print("Solar sky matrix")
    path, name = ioRoutines.findPathName(meteo_path, "/", ".epw")
    # sunSmx = ioRoutines.skyMtx('%s%sSun.smx' % (path,name),iniHour,totHour,dayHour,sunAltitude)
    solSmx = ioRoutines.skyMtx(
        "%s%sSol.smx" % (path, name), iniHour, totHour, dayHour, sunAltitude
    )
    #        print('Time span (s): ',time.time()-start_time)
    #        start_time = time.time()
    # ----------------------------------------------------
    return solSmx


def tpmMtxTherm(numWin, numFrame, numSurf, numConWin, bc, workDir, viewFactor_path):
    # ----------------------------------------------------
    # view factors
    vfFrameWin = atleast_1d(
        genfromtxt("%svfFrameWindows.dat" % viewFactor_path, delimiter=" ", usecols=2)
    )
    vfWinFrame = atleast_1d(
        genfromtxt("%svfWindowsFrame.dat" % viewFactor_path, delimiter=" ", usecols=2)
    )
    vfFrameSurf = atleast_1d(
        genfromtxt("%svfFrameSurface.dat" % viewFactor_path, delimiter=" ", usecols=2)
    )
    vfSurfFrame = atleast_1d(
        genfromtxt("%svfSurfaceFrame.dat" % viewFactor_path, delimiter=" ", usecols=2)
    )
    vfSurfWin = atleast_1d(
        genfromtxt("%svfSurfaceWindows.dat" % viewFactor_path, delimiter=" ", usecols=2)
    )
    vfWinSurf = atleast_1d(
        genfromtxt("%svfWindowsSurface.dat" % viewFactor_path, delimiter=" ", usecols=2)
    )
    vfWinWin = atleast_1d(
        genfromtxt("%svfWindowsWindows.dat" % viewFactor_path, delimiter=" ", usecols=2)
    )
    vfSurfSurf = atleast_1d(
        genfromtxt("%svfSurfaceSurface.dat" % viewFactor_path, delimiter=" ", usecols=2)
    )
    vfFrameFrame = atleast_1d(
        genfromtxt("%svfFrameFrame.dat" % viewFactor_path, delimiter=" ", usecols=2)
    )

    # ----------------------------------------------------
    numSensorSurf = atleast_1d(genfromtxt("%snumSensorSurf.dat" % workDir, dtype="int"))
    numSensorFrame = atleast_1d(
        genfromtxt("%snumSensorFrame.dat" % workDir, dtype="int")
    )
    numSensorWin = atleast_1d(genfromtxt("%snumSensorWin.dat" % workDir, dtype="int"))
    numSensorRev = atleast_1d(genfromtxt("%snumSensorRev.dat" % workDir, dtype="int"))
    # ----------------------------------------------------
    print("Irradiance daylight and view matrices")
    irrSurfExtDmx = zeros([numSurf, max(numSensorSurf), 2306])
    irrFrameExtDmx = zeros([numFrame, max(numSensorFrame), 2306])
    irrSurfFrontVmx = zeros([numWin, numSurf, max(numSensorSurf), 145])
    irrSurfBackVmx = zeros([numWin, numSurf, max(numSensorSurf), 145])
    irrFrameFrontVmx = zeros([numWin, numSurf, max(numSensorFrame), 145])
    irrFrameBackVmx = zeros([numWin, numSurf, max(numSensorFrame), 145])
    irrWinFrontVmx = zeros([numWin, max(numSensorWin), 145])
    irrWinBackVmx = zeros([numWin, max(numSensorWin), 145])
    irrRevFrontVmx = zeros([numWin, numSurf, max(numSensorRev), 145])
    irrRevBackVmx = zeros([numWin, numSurf, max(numSensorRev), 145])
    #        dmx = zeros([numWin,145,2306])
    for k in range(numSurf):
        if bc[k]:
            aux_irrSurfExtDmx = atleast_2d(
                genfromtxt(
                    "%sirrSurfExt_%i.dmx" % (workDir, k), usecols=arange(0, 2306 * 3, 3)
                )
            )
            for j in range(numSensorSurf[k]):
                irrSurfExtDmx[k, j, :] = aux_irrSurfExtDmx[j, :]
    for k in range(numFrame):
        aux_irrFrameExtDmx = genfromtxt(
            "%sirrFrameExt_%i.dmx" % (workDir, k), usecols=arange(0, 2306 * 3, 3)
        )
        for j in range(numSensorFrame[k]):
            irrFrameExtDmx[k, j, :] = aux_irrFrameExtDmx[j, :]
    for i in range(numWin):
        #           dmx[i,:,:] = genfromtxt('%swin_%i.dmx'%(workDir,i),usecols=arange(0,2306*3,3))
        aux_irrWinFrontVmx = genfromtxt(
            "%sirrWinFront_%i.vmx" % (workDir, i), usecols=arange(0, 145 * 3, 3)
        )
        aux_irrWinBackVmx = genfromtxt(
            "%sirrWinBack_%i.vmx" % (workDir, i), usecols=arange(0, 145 * 3, 3)
        )
        for j in range(numSensorWin[i]):
            irrWinFrontVmx[i, j, :] = aux_irrWinFrontVmx[j, :]
            irrWinBackVmx[i, j, :] = aux_irrWinBackVmx[j, :]
        for k in range(numSurf):
            aux_irrSurfFrontVmx = atleast_2d(
                genfromtxt(
                    "%sirrSurfFront_%i_%i.vmx" % (workDir, i, k),
                    usecols=arange(0, 145 * 3, 3),
                )
            )
            aux_irrSurfBackVmx = atleast_2d(
                genfromtxt(
                    "%sirrSurfBack_%i_%i.vmx" % (workDir, i, k),
                    usecols=arange(0, 145 * 3, 3),
                )
            )
            for j in range(numSensorSurf[k]):
                irrSurfFrontVmx[i, k, j, :] = aux_irrSurfFrontVmx[j, :]
                irrSurfBackVmx[i, k, j, :] = aux_irrSurfBackVmx[j, :]
        for k in range(numFrame):
            aux_irrFrameFrontVmx = genfromtxt(
                "%sirrFrameFront_%i_%i.vmx" % (workDir, i, k),
                usecols=arange(0, 145 * 3, 3),
            )
            aux_irrFrameBackVmx = genfromtxt(
                "%sirrFrameBack_%i_%i.vmx" % (workDir, i, k),
                usecols=arange(0, 145 * 3, 3),
            )
            aux_irrRevFrontVmx = genfromtxt(
                "%sirrRevFront_%i_%i.vmx" % (workDir, i, k),
                usecols=arange(0, 145 * 3, 3),
            )
            aux_irrRevBackVmx = genfromtxt(
                "%sirrRevBack_%i_%i.vmx" % (workDir, i, k),
                usecols=arange(0, 145 * 3, 3),
            )
            for j in range(numSensorFrame[k]):
                irrFrameFrontVmx[i, k, j, :] = aux_irrFrameFrontVmx[j, :]
                irrFrameBackVmx[i, k, j, :] = aux_irrFrameBackVmx[j, :]
            for j in range(numSensorRev[k]):
                irrRevFrontVmx[i, k, j, :] = aux_irrRevFrontVmx[j, :]
                irrRevBackVmx[i, k, j, :] = aux_irrRevBackVmx[j, :]
    #        print('Time span (s): ',time.time()-start_time)
    # ----------------------------------------------------
    return (
        irrSurfExtDmx,
        irrFrameExtDmx,
        irrSurfFrontVmx,
        irrSurfBackVmx,
        irrFrameFrontVmx,
        irrFrameBackVmx,
        irrWinFrontVmx,
        irrWinBackVmx,
        irrRevFrontVmx,
        irrRevBackVmx,
        numSensorSurf,
        numSensorFrame,
        numSensorRev,
        numSensorWin,
        vfFrameWin,
        vfWinFrame,
        vfFrameSurf,
        vfSurfFrame,
        vfSurfWin,
        vfWinSurf,
        vfWinWin,
        vfFrameFrame,
        vfSurfSurf,
    )
