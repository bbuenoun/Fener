#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function

import traceback
import argument_parser
import config_parser
import os
import sys
from numpy import (
    savetxt,
    transpose,
    ndarray,
    ones,
    zeros,
    radians,
    mean,
    mod,
    fromiter,
    genfromtxt,
    atleast_2d,
    matrix,
)
import genMtx
import variables
import geometry
import solar
import param as p
import thermal
import ioRoutines
import time
import calcs
import elements
import control
import klems
import thermalComfort
from collections import defaultdict
import multiprocessing as mp

# ====================================================
#
# Fraunhofer Institute for Solar Energie Systems ISE
# Group Solar Building Envelopes
# Heidenhofstr. 2, 79110 Freiburg, Germany
# http://www.ise.fraunhofer.de
#
# Author: Bruno Bueno
# bruno.bueno@ise.fraunhofer.de
#
# ====================================================


def main(opts, config):

    # ----------------------------------------------------
    if opts.meteo:
        print("Building sky matrix...")
        genMtx.sky(
            opts.c,
            config.meteo_path,
            config.lat,
            config.lon,
            config.tzone,
            config.grndAlb,
        )
    # ----------------------------------------------------
    print("Initializing...")
    inobis = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
    iniHour = (inobis[config.iniMonth - 1] + config.iniDay - 1) * 24
    finHour = (inobis[config.endMonth - 1] + config.endDay) * 24
    totHour = finHour - iniHour
    # ----------------------------------------------------
    mm = genfromtxt(
        config.meteo_path,
        skip_header=8 + iniHour,
        skip_footer=8760 - finHour,
        delimiter=",",
        usecols=1,
    )
    dd = genfromtxt(
        config.meteo_path,
        skip_header=8 + iniHour,
        skip_footer=8760 - finHour,
        delimiter=",",
        usecols=2,
    )
    hh = (
        genfromtxt(
            config.meteo_path,
            skip_header=8 + iniHour,
            skip_footer=8760 - finHour,
            delimiter=",",
            usecols=3,
        )
        - 0.5
    )
    irrdirn = genfromtxt(
        config.meteo_path,
        skip_header=8 + iniHour,
        skip_footer=8760 - finHour,
        delimiter=",",
        usecols=14,
    )
    irrdiffh = genfromtxt(
        config.meteo_path,
        skip_header=8 + iniHour,
        skip_footer=8760 - finHour,
        delimiter=",",
        usecols=15,
    )
    tempOutAir = (
        genfromtxt(
            config.meteo_path,
            skip_header=8 + iniHour,
            skip_footer=8760 - finHour,
            delimiter=",",
            usecols=6,
        )
        + 273.15
    )  # K
    path, name = ioRoutines.findPathName(config.meteo_path, "/", ".epw")
    sunAltitude = genfromtxt(
        "%s%s_sunAng.dat" % (path, name),
        skip_header=1 + iniHour,
        skip_footer=8760 - finHour,
        delimiter=",",
        usecols=0,
    )  # deg
    sunAzimuth = genfromtxt(
        "%s%s_sunAng.dat" % (path, name),
        skip_header=1 + iniHour,
        skip_footer=8760 - finHour,
        delimiter=",",
        usecols=1,
    )  # deg
    dayHour = sum(x > 0 for x in sunAltitude)
    # ----------------------------------------------------
    if opts.shoeBox:
        geometry.shoeBox(
            config.surf_path,
            config.rotAng,
            config.length,
            config.width,
            config.height,
            config.thickSouth,
            config.thickEast,
            config.thickNorth,
            config.thickWest,
            config.thickCeiling,
            config.thickFloor,
            config.albWall,
            config.albCeiling,
            config.albFloor,
            config.bcSouth,
            config.bcEast,
            config.bcNorth,
            config.bcWest,
            config.bcCeiling,
            config.bcFloor,
            config.conSouth,
            config.conEast,
            config.conNorth,
            config.conWest,
            config.conCeiling,
            config.conFloor,
        )
    surfVect = atleast_2d(genfromtxt(config.surf_path, skip_header=1, delimiter=","))
    winVect = atleast_2d(genfromtxt(config.win_path, skip_header=1, delimiter=","))
    frameVect = atleast_2d(genfromtxt(config.frame_path, skip_header=1, delimiter=","))
    numWin = winVect.shape[0]
    numFrame = frameVect.shape[0]
    numSurf = surfVect.shape[0]
    # ----------------------------------------------------
    if opts.mask:
        obstMask = genfromtxt(config.obstMask_path, skip_header=1, delimiter=",")
    # ----------------------------------------------------
    if opts.lightSch or opts.therm:
        lightSchYear = genfromtxt(config.lightSch_path, delimiter=",")
        lightSch = genfromtxt(
            config.lightSch_path,
            skip_header=iniHour,
            skip_footer=8760 - finHour,
            delimiter=",",
        )
    if opts.lightSch:
        lightCntrl = atleast_2d(
            genfromtxt(config.lightCntrl_path, skip_header=1, delimiter=",")
        )
    # ----------------------------------------------------
    if opts.mtxCntrl or opts.refeedCntrl:
        thermFlag = False
        glareFlag = False
        if opts.therm:
            thermFlag = True
        if opts.glare or opts.glareSimpl or opts.glareFull:
            glareFlag = True
        cntrlMtx = elements.cntrlMtx(
            config.cntrlMtx_path, opts.daylight, glareFlag, thermFlag, numWin
        )
    # ----------------------------------------------------
    if opts.schCntrl:
        shadingSch = genfromtxt(
            config.shadingSch_path,
            skip_header=iniHour,
            skip_footer=8760 - finHour,
            delimiter=",",
        )
    # ----------------------------------------------------
    if opts.therm:
        # ------------------------------------------------
        skyTemp = (
            genfromtxt(
                config.meteo_path,
                skip_header=8 + iniHour,
                skip_footer=8760 - finHour,
                delimiter=",",
                usecols=12,
            )
            / p.sigma
        ) ** (0.25)
        uWind = genfromtxt(
            config.meteo_path,
            skip_header=8 + iniHour,
            skip_footer=8760 - finHour,
            delimiter=",",
            usecols=21,
        )
        extCHTC = 2.62 + 7.72 * uWind ** (0.5)
        extRHTC = 4 * p.sigma * (0.5 * skyTemp + 0.5 * tempOutAir) ** 3
        matOpaq = atleast_2d(
            genfromtxt(config.matOpaq_path, skip_header=1, delimiter=",")
        )
        conOpaq = atleast_2d(
            genfromtxt(config.constOpaq_path, skip_header=1, delimiter=",")
        )
        infSch = genfromtxt(
            config.infSch_path,
            skip_header=iniHour,
            skip_footer=8760 - finHour,
            delimiter=",",
        )
        occSch = genfromtxt(
            config.occSch_path,
            skip_header=iniHour,
            skip_footer=8760 - finHour,
            delimiter=",",
        )
        equipSch = genfromtxt(
            config.equipSch_path,
            skip_header=iniHour,
            skip_footer=8760 - finHour,
            delimiter=",",
        )
        heatSpSch = (
            genfromtxt(
                config.heatSpSch_path,
                skip_header=iniHour,
                skip_footer=8760 - finHour,
                delimiter=",",
            )
            + 273.15
        )
        coolSpSch = (
            genfromtxt(
                config.coolSpSch_path,
                skip_header=iniHour,
                skip_footer=8760 - finHour,
                delimiter=",",
            )
            + 273.15
        )
    # ----------------------------------------------------
    if opts.thermComf:
        relHumOut = genfromtxt(
            config.meteo_path,
            skip_header=8 + iniHour,
            skip_footer=8760 - finHour,
            delimiter=",",
            usecols=8,
        )  # RH%
        for i in range(len(relHumOut)):
            if relHumOut[i] > 100:
                relHumOut[i] = 100
    # -----------------------------------------------------
    if opts.klems:
        matGlz = atleast_2d(
            genfromtxt(config.matGlz_path, skip_header=1, delimiter=",")
        )
        conWinVect = atleast_2d(
            genfromtxt(config.constWin_path, skip_header=1, delimiter=",")
        )
        # ----------------------------------------------------
    if opts.klems:
        print("Computing Klems BSDF...")
        if config.numBsdfLay > 0:
            klems.main(
                conWinVect,
                matGlz,
                config.bsdfSys_path,
                config.absFront_path,
                config.absBack_path,
                config.numBsdfLay,
                config.bsdfLay_file,
            )
        else:
            klems.main(
                conWinVect,
                matGlz,
                config.bsdfSys_path,
                config.absFront_path,
                config.absBack_path,
            )
    # ----------------------------------------------------
    if opts.daylight or opts.glare or opts.glareSimpl or opts.df:
        tmxv = variables.tmxVis(config.numConWin, config.bsdfSys_path)
    if opts.therm:
        tmxs, rmxs = variables.tmxSol(config.numConWin, config.bsdfSys_path)
    if opts.therm:
        calorim = variables.calorim(config.numConWin, config.calorim_path)
    # ----------------------------------------------------
    conIndex = zeros([numWin, totHour], int)
    for i in range(numWin):
        conIndex[i, 0] = int(winVect[i, 7])
    # ----------------------------------------------------
    surf = ndarray((numSurf,), dtype=object)
    for i in range(numSurf):
        if opts.therm:
            surf[i] = elements.surf(
                    surfVect[i, :], frameVect, i, 1, conOpaq, matOpaq, config.iniTemp
            )
        else:
            surf[i] = elements.surf(surfVect[i, :], frameVect, i)
    bc = fromiter((c.bc for c in surf), float)
    areaSurf = fromiter((c.area for c in surf), float)
    areaRev = fromiter((c.areaRev for c in surf), float)
    # ----------------------------------------------------
    frame = ndarray((numFrame,), dtype=object)
    for i in range(numFrame):
        if opts.therm:
            frame[i] = elements.frame(
                frameVect[i, :], winVect, surf, i, 1, config.iniTemp
            )
        else:
            frame[i] = elements.frame(frameVect[i, :], winVect, surf, i)
    frameSurf = fromiter((c.surf for c in frame), int)
    areaFrame = fromiter((c.area for c in frame), float)
    # ----------------------------------------------------
    win = ndarray((numWin,), dtype=object)
    for i in range(numWin):
        if opts.therm:
            win[i] = elements.win(
                    winVect[i, :],
                    frame,
                    thermFlag=True,
                    iniTemp=config.iniTemp,
                    tmxs=tmxs,
                    rmxs=rmxs,
                    calorim=calorim,
            )
            if opts.cutOff or opts.mtxCntrl or opts.refeedCntrl or opts.schCntrl:
                print("WARNING: the number of layers cannot be dynamically changed")
        else:
            win[i] = elements.win(winVect[i, :], frame)
    areaWin = fromiter((c.area for c in win), float)
    # ----------------------------------------------------
    if opts.geo:
        print("Generating geometry...")
        geometry.fileDel("%s*" % config.workDir_path)
        geometry.zone(surf, frame, win, config.workDir_path)
        # generate obstruction mask
        if opts.mask:
            geometry.genmask(
                surf[config.floor],
                config.height,
                obstMask,
                config.workDir_path,
                opts.mask,
            )
        elif opts.outside:
            ioRoutines.shell(
                "cp -f %s %soutside.rad" % (config.outside_path, config.workDir_path)
            )
        else:
            open("%soutside.rad" % config.workDir_path, "a").close()
        # show geometry
        if opts.display:
            windowRad = "%swindow_0.rad\t" % config.workDir_path
            for i in range(numWin - 1):
                windowRad += "%swindow_%i.rad\t" % (config.workDir_path, i + 1)
            ioRoutines.shell(
                "objline %sroom.rad %s %soutside.rad | x11meta"
                % (config.workDir_path, windowRad, config.workDir_path)
            )
        # daylight input files
        if opts.daylight:
            if opts.grid:
                geometry.photocells(
                    config.numPhotosensX,
                    config.numPhotosensY,
                    config.photosensHeight,
                    config.gridXOffset,
                    config.gridYOffset,
                    config.thickSouth,
                    config.thickEast,
                    config.thickNorth,
                    config.thickWest,
                    config.thickFloor,
                    surf[config.floor],
                    config.workDir_path,
                )
            else:
                geometry.genPtsFile(
                    config.illuPts_path,
                    "%silluSensor.pts" % config.workDir_path,
                    config.rotAng,
                )
        if opts.df:
            geometry.extIlluSensor(
                config.workDir_path,
                surf[config.floor].length / 2.0,
                surf[config.floor].height / 2.0,
                10.0,
            )
        # glare input files
        if opts.glare or opts.glareSimpl or opts.glareFull:
            illuVertPts = geometry.genPtsFile(
                config.illuVertPts_path,
                "%silluVertSensor.pts" % config.workDir_path,
                config.rotAng,
            )
            geometry.genViewFile(illuVertPts, config.workDir_path)
        # therm input files
        if opts.therm:
            geometry.irradianceSensorPts(
                p.irrResSurf,
                p.irrResFrame,
                p.irrResWin,
                p.irrResRev,
                surf,
                frame,
                win,
                config.workDir_path,
            )
            print("Calculating view factors...")
            if opts.shoeBox:
                geometry.shoeBoxVF(
                    config.workDir_path,
                    config.length,
                    config.width,
                    config.height,
                    config.surf_path,
                    config.frame_path,
                    config.win_path,
                )
            else:
                geometry.viewfactors(opts.c, config.workDir_path)
        elif opts.daylight:
            geometry.winSensorPts(p.irrResWin, surf, frame, win, config.workDir_path)
    # ----------------------------------------------------
    if opts.geo:
        print("Building three-phase method matrices...")
        if (
            opts.daylight
            or opts.glare
            or opts.glareSimpl
            or opts.therm
            or opts.df
        ):
            genMtx.dayH(opts.c, numWin, config.workDir_path)
        if opts.daylight or opts.df:
            genMtx.illViewH(opts.c, numWin, config.workDir_path)
        if opts.df:
            genMtx.extIll(opts.c, config.workDir_path)
        if opts.glare or opts.glareSimpl:
            genMtx.glareView(opts.c, numWin, config.workDir_path)
        if opts.therm:
            genMtx.irrView(opts.c, numWin, numFrame, numSurf, config.workDir_path)
            genMtx.ext(opts.c, numSurf, numFrame, bc, config.workDir_path)
        elif opts.daylight:
            genMtx.irrWinBackView(opts.c, numWin, config.workDir_path)
    # ----------------------------------------------------
    # read matrices
    if opts.daylight or opts.glare or opts.glareSimpl:
        visSmx = variables.smxVis(
            config.meteo_path, iniHour, totHour, dayHour, sunAltitude
        )
    if (
        opts.daylight
        or opts.glare
        or opts.glareSimpl
        or opts.df
        or opts.therm
    ):
        dmx = variables.dmxVisH(numWin, config.workDir_path)
    if opts.daylight or opts.df or opts.rad == "radTpm":
        if opts.grid:
            numSensor = config.numPhotosensX * config.numPhotosensY
        else:
            numSensor = atleast_2d(
                genfromtxt(config.illuPts_path, delimiter=" ")
            ).shape[0]
    if opts.daylight or opts.df:
        illVmx = variables.vmxIllH(numWin, numSensor, config.workDir_path)
    if opts.df:
        illExtDmx = variables.dmxExt(config.workDir_path)
    if opts.glare or opts.glareSimpl or opts.glareFull:
        illuVertPts = atleast_2d(
            genfromtxt("%silluVertSensor.pts" % config.workDir_path, delimiter=" ")
        )
        numSensorVert = illuVertPts.shape[0]
        illVertVmx = variables.vmxGlare(numWin, numSensorVert, config.workDir_path)
    else:
        numSensorVert = 0
    if opts.therm:
        solSmx = variables.smxSol(
            config.meteo_path, iniHour, totHour, dayHour, sunAltitude
        )
        (
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
        ) = variables.tpmMtxTherm(
            numWin,
            numFrame,
            numSurf,
            config.numConWin,
            bc,
            config.workDir_path,
            config.workDir_path,
        )
    elif opts.daylight:
        irrWinBackVmx, numSensorWin = variables.winBack(numWin, config.workDir_path)
    # ----------------------------------------------------
    # initialize outputs
    if opts.daylight or opts.df or opts.mtxCntrl or opts.refeedCntrl:
        ill = matrix(zeros([numSensor, totHour]))
        illTrans = zeros(totHour)
        if opts.da:
            goodHour = zeros(numSensor)
            da = zeros(numSensor)
            totWorkHour = 0
            directSunHour = zeros(numSensor)
    if opts.rad == "radTpm":
        illRadTpm = matrix(zeros([numSensor, totHour]))
    if opts.df:
        df = matrix(zeros([numSensor, totHour]))
    if opts.lightSch or opts.therm:
        lightCons = zeros(totHour)
    if (
        opts.glare
        or opts.glareSimpl
        or opts.mtxCntrl
        or opts.refeedCntrl
    ):
        illVert = matrix(zeros([numSensorVert, totHour]))
        dgp = zeros([numSensorVert, totHour])
    if opts.glareFull:
        dgpFull = zeros([numSensorVert, totHour])
    if opts.therm or opts.mtxCntrl or opts.refeedCntrl:
        solTrans = zeros(totHour)
        convIntHeatFlx = zeros(totHour)
        radIntHeatFlx = zeros(totHour)
        irrSurfExt = zeros([numSurf, totHour])
        irrSurfInt = zeros([numSurf, totHour])
        infrSurfInt = zeros([numSurf, totHour])
        convSurfInt = zeros([numSurf, totHour])
        heatSurfInt = zeros([numSurf, totHour])
        tempSurfInt = zeros([numSurf, totHour])
        tempSurfExt = zeros([numSurf, totHour])
        irrFrameExt = zeros([numFrame, totHour])
        irrFrameInt = zeros([numFrame, totHour])
        irrWinAbsExt = zeros([numWin, totHour])
        irrWinAbsInt = zeros([numWin, totHour])
        irrWinInt = zeros([numWin, totHour])
        infrWinInt = zeros([numWin, totHour])
        convWinInt = zeros([numWin, totHour])
        heatWinInt = zeros([numWin, totHour])
        tempWinInt = zeros([numWin, totHour])
        tempWinExt = zeros([numWin, totHour])
        tempInAir = ones(totHour) * config.iniTemp
        energyDemand = zeros(totHour)
        tempInAirTemp = tempInAir[0]
        effGValue = zeros([numWin, totHour])
    if opts.thermComf:
        pmv = zeros([5, totHour])
        ppd = zeros([5, totHour])
        RHIn = ones(totHour) * 75
        RHInTemp = RHIn[0]
        tWalls = zeros([numSurf, totHour])
        tRad = zeros([5, totHour])
    if opts.mtxCntrl or opts.refeedCntrl:
        tempOut24 = ones(int(totHour / 24)) * tempOutAir[0]
        tempOut24temp = 0.0
        itd = 0
    if opts.optStateCntrl: 
        spHorIllu = 1000
        spVerIllu = 2000
    # ----------------------------------------------------
    # daylight factor calculation
    if opts.df:
        df = solar.df(
            ill[:, 0],
            illExtDmx,
            illVmx,
            tmxv,
            conIndex[:, 0],
            dmx,
            numWin,
            config.lat,
            config.lon,
            config.tzone,
            config.workDir_path,
        )
    # ----------------------------------------------------
    # RadianceClassic calculation with BSDF material
    if opts.rad == "radBsdf":
        for ihour in range(len(opts.hour)):
            rad.classicBsdf(
                config.iniMonth,
                config.iniDay + opts.skipDays,
                opts.hour[ihour],
                irrdirn[int(opts.skipDays * 24 + opts.hour[ihour])],
                irrdiffh[int(opts.skipDays * 24 + opts.hour[ihour])],
                opts.c,
                opts.ab,
                opts.ad,
                opts.ar,
                opts.ass,
                opts.aa,
                conIndex[:, 0],
                numWin,
                config.lat,
                config.lon,
                config.tzone,
                config.workDir_path,
                config.output_path,
                config.bsdfSys_path,
            )
    # ----------------------------------------------------
    # daylight hour counter
    dh = 0
    start_time = time.time()
    print("Running dynamic simulation...")
    for it in range(totHour):
        # print time
        if (
            opts.daylight
            or opts.glare
            or opts.glareSimpl
            or opts.glareFull
            or opts.therm
        ):
            if hh[it] > 23.0:
                print(int(mm[it]), int(dd[it]))
        # --------------------------------------------------------
        # matrix shade control algorithm
        if (opts.mtxCntrl or opts.refeedCntrl) and it > 0:
            # incidence angle on window 0, assuming vertical
            incAng = calcs.incAngle(
                radians(sunAltitude[it]),
                radians(sunAzimuth[it]),
                radians(-win[0].rz),
                1,
            )  # radians
            # profile angle on window 0, assuming vertical and horizontal slats
            profAng = calcs.profAngle(
                radians(sunAltitude[it]),
                radians(sunAzimuth[it]),
                radians(-win[0].rz),
                1,
            )  # radians
            # exterior radiation on "south facade"
            if sunAltitude[it] > 0.0:
                irrSurfExtSouth = solar.irrFadExt(
                    numSurf, numSensorSurf, irrSurfExtDmx, transpose(solSmx[dh, :]), bc
                )[1]
            else:
                irrSurfExtSouth = 0.0
            if opts.mtxCntrl and not opts.optStateCntrl:
                conIndex[:, it] = control.mtxCntrl(
                    cntrlMtx,
                    config.conDef,
                    occSch[it],
                    sunAltitude[it],
                    mean(ill[:, it - 1]),
                    max(dgp[:, it - 1]),
                    tempInAir[it - 1],
                    irrSurfExtSouth,
                    tempOut24[itd],
                    incAng,
                    profAng,
                    conIndex[:, it - 1],
                )
            if opts.refeedCntrl and dh < visSmx.shape[0] and not opts.cutOff:
                # control algorithm
                conIndex[:, it] = control.refeed(
                    cntrlMtx,
                    config.conDef,
                    occSch[it],
                    sunAltitude[it],
                    ill[:, it - 1],
                    dgp[:, it - 1],
                    tempInAir[it - 1],
                    irrSurfExtSouth,
                    tempOut24[itd],
                    incAng,
                    profAng,
                    conIndex[:, it - 1],
                    illVmx,
                    tmxv,
                    dmx,
                    transpose(visSmx[dh, :]),
                    numWin,
                    illVert[:, it - 1],
                    illVertVmx,
                    numSensor,
                    numSensorVert,
                )
            if opts.refeedCntrl and opts.cutOff and dh < visSmx.shape[0]:
                # control algorithm
                conIndex[:, it] = control.refeedCutOff(
                    cntrlMtx,
                    config.conDef,
                    occSch[it],
                    sunAltitude[it],
                    ill[:, it - 1],
                    dgp[:, it - 1],
                    tempInAir[it - 1],
                    irrSurfExtSouth,
                    tempOut24[itd],
                    incAng,
                    profAng,
                    conIndex[:, it - 1],
                    illVmx,
                    tmxv,
                    dmx,
                    transpose(visSmx[dh, :]),
                    numWin,
                    illVert[:, it - 1],
                    illVertVmx,
                    numSensor,
                    numSensorVert,
                    sunAzimuth[it],
                    -win[0].rz,
                    irrdirn[it],
                    config.slatWidth,
                    config.slatDist,
                    config.minAng,
                    config.maxAng,
                    config.stepAng,
                    config.conMinAng,
                    config.conMaxAng,
                )
            # Average outside temperature over 24 hours
            tempOut24temp = tempOut24temp + tempOutAir[it] / 24.0
            if mod(it, 24) == 0:
                itd = itd + 1
                tempOut24[itd] = tempOut24temp
                tempOut24temp = 0.0
            # --------------------------------------------------------
            # optimum state control algorithm
            if opts.mtxCntrl and opts.optStateCntrl and dh < visSmx.shape[0]:
                conIndex[:, it] = control.mtxCntrlOptState(
                    cntrlMtx,
                    config.conDef,
                    occSch[it],
                    sunAltitude[it],
                    mean(ill[:, it - 1]),
                    max(dgp[:, it - 1]),
                    tempInAir[it - 1],
                    irrSurfExtSouth,
                    tempOut24[itd],
                    incAng,
                    profAng,
                    conIndex[:, it-1],
                    config.numConWin, 
                    illVmx, 
                    illVertVmx,
                    tmxv, 
                    dmx, 
                    transpose(visSmx[dh, :]), 
                    numWin,
                    spHorIllu,
                    spVerIllu,
                    numSensor,
                    numSensorVert,
                    config.numWinLowerPartition,
                    config.numWinUpperPartition,
                    config.cntrlOpt
                    )
        # --------------------------------------------------------
        # cut-off control algorithm
        else:
            if opts.cutOff and dh < visSmx.shape[0]:
                conIndex[:, it] = control.simpleCutOff(
                    sunAltitude[it],
                    sunAzimuth[it],
                    -win[0].rz,
                    irrdirn[it],
                    config.slatWidth,
                    config.slatDist,
                    config.minAng,
                    config.maxAng,
                    config.stepAng,
                    config.conMinAng,
                    config.conMaxAng,
                )
        # --------------------------------------------------------
        # scheduled shade control algorithm
        if opts.schCntrl:
            for i in range(numWin):
                conIndex[i, it] = shadingSch[it]
        # --------------------------------------------------------
        # daytime
        if sunAltitude[it] > 0.0:
            if opts.daylight or opts.glare or opts.glareSimpl:
                skyVis = visSmx[dh, :]
            if opts.therm:
                skySol = solSmx[dh, :]
            # illuminance
            if opts.daylight:
                ill[:, it] = solar.illum(
                    ill[:, it],
                    illVmx,
                    tmxv,
                    conIndex[:, it],
                    dmx,
                    transpose(skyVis),
                    numWin,
                )
                illTrans[it] = solar.illumTrans(
                    numWin,
                    numSensorWin,
                    areaWin,
                    irrWinBackVmx,
                    tmxv,
                    conIndex[:, it],
                    dmx,
                    transpose(skyVis),
                    illTrans[it],
                )
            if opts.rad == "radTpm":
                illRadTpm[:, it] = rad.multiply(
                    opts.c,
                    illRadTpm[:, it],
                    numWin,
                    numSensor,
                    config.workDir_path,
                    config.bsdfSys_path,
                    conIndex[:, it],
                    mm[it],
                    dd[it],
                    hh[it],
                    irrdirn[it],
                    irrdiffh[it],
                    config.lat,
                    config.lon,
                    config.tzone,
                )
            # glare calculation
            # ~ start_time_glare = time.time()
            if opts.glare:
                illVert[:, it], dgp[:, it] = solar.glare(
                    opts.c,
                    illVert[:, it],
                    dgp[:, it],
                    illVertVmx,
                    tmxv,
                    conIndex[:, it],
                    dmx,
                    transpose(skyVis),
                    numWin,
                    mm[it],
                    dd[it],
                    hh[it],
                    irrdirn[it],
                    irrdiffh[it],
                    config.lat,
                    config.lon,
                    config.tzone,
                    illuVertPts,
                    numSensorVert,
                    config.winRad_path,
                    config.workDir_path,
                    config.input_path,
                )
            if opts.glareSimpl:
                illVert[:, it], dgp[:, it] = solar.glareSimpl(
                    illVert[:, it],
                    dgp[:, it],
                    illVertVmx,
                    tmxv,
                    conIndex[:, it],
                    dmx,
                    transpose(skyVis),
                    numWin,
                    numSensorVert,
                )
            if opts.glareFull:
                dgpFull[:, it] = solar.glareFull(
                    opts.c,
                    dgpFull[:, it],
                    conIndex[:, it],
                    numWin,
                    mm[it],
                    dd[it],
                    hh[it],
                    irrdirn[it],
                    irrdiffh[it],
                    config.lat,
                    config.lon,
                    config.tzone,
                    illuVertPts,
                    numSensorVert,
                    config.winRad_path,
                    config.workDir_path,
                    config.input_path,
                )
            # ~ print('Time glare',time.time()-start_time_glare)
            # exterior and interior solar radiation calculations
            if opts.therm:
                (
                    solTrans[it],
                    irrSurfExt[:, it],
                    irrSurfInt[:, it],
                    irrFrameExt[:, it],
                    irrFrameInt[:, it],
                    irrWinInt[:, it],
                ) = solar.therm(
                    transpose(skySol),
                    win,
                    numWin,
                    areaWin,
                    numSensorSurf,
                    numSensorFrame,
                    numSensorWin,
                    numSensorRev,
                    numSurf,
                    numFrame,
                    frameSurf,
                    areaSurf,
                    areaRev,
                    bc,
                    irrSurfExtDmx,
                    irrFrameExtDmx,
                    dmx,
                    irrWinFrontVmx,
                    irrWinBackVmx,
                    irrSurfFrontVmx,
                    irrSurfBackVmx,
                    irrFrameFrontVmx,
                    irrFrameBackVmx,
                    irrRevFrontVmx,
                    irrRevBackVmx,
                    tmxs,
                    conIndex[:, it],
                    irrSurfInt[:, it],
                    irrSurfExt[:, it],
                    irrFrameInt[:, it],
                    irrFrameExt[:, it],
                    solTrans[it],
                )
            if opts.therm:
                irrWinAbsExt[:, it], irrWinAbsInt[:, it] = solar.calorimWinAbs(
                    numWin, win, dmx, transpose(skySol)
                )
                effGValue[:, it] = solar.effGVal(
                    numWin, win, dmx, transpose(skySol), calorim, conIndex[:, it]
                )
            dh = dh + 1
        # ----------------------------------------------------
        # daylight authonomy and annual (direct) sunlight exposure
        if opts.da:
            goodHour, totWorkHour = solar.daGoodHour(
                ill[:, it],
                goodHour,
                config.minIllum,
                numSensor,
                config.startHour,
                config.endHour,
                hh[it],
                mm[it],
                totWorkHour,
            )
            directSunHour = solar.directSunlightHour(
                ill[:, it], directSunHour, numSensor
            )
        # ----------------------------------------------------
        # lighting power calculation
        if opts.lightSch:
            lightSch[it], lightCons[it] = solar.lightCntrl(
                 lightCntrl,
                 lightSch[it],
                 config.floorArea,
                 config.powerLight,
                 ill[:, it],
            )
            lightSchYear[iniHour + it] = lightSch[it]
        elif opts.therm:
            lightCons[it] = config.floorArea * config.powerLight * lightSch[it]
        # ----------------------------------------------------
        # Thermal simulation
        if opts.therm:
            # update window construction
            if opts.cutOff or opts.mtxCntrl or opts.refeedCntrl or opts.schCntrl:
                for i in range(numWin):
                    win[i].kuhn11(
                        tmxs[conIndex[i, it], :, :],
                        rmxs[conIndex[i, it], :, :],
                        calorim[conIndex[i, it], :],
                    )
            # internal heat gains
            convIntHeatFlx[it], radIntHeatFlx[it] = thermal.genIntHeatGains(
                config.floorArea,
                config.powerLight,
                config.radFracLight,
                lightSch[it],
                config.powerEquip,
                config.radFracEquip,
                equipSch[it],
                config.powerPeople,
                config.radFracPeople,
                occSch[it],
            )
            for itt in range(int(3600 / p.timeStep)):
                # radiant heat transfer coefficients
                (
                    RHTCwinFrame,
                    RHTCframeWin,
                    RHTCwinSurf,
                    RHTCsurfWin,
                    RHTCframeSurf,
                    RHTCsurfFrame,
                    RHTCwinWin,
                    RHTCframeFrame,
                    RHTCsurfSurf,
                ) = thermal.genRHTC(
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
                )
                # convective heat transfer coefficients
                CHTCwin, CHTCframe, CHTCsurf = thermal.genCHTC(
                    numWin, numFrame, numSurf, tempInAirTemp, win, frame, surf
                )
                # Transfer functions for transient heat condution and window heat transfer Kuhn2011
                (
                    energyDemand[it],
                    tempInAirTemp,
                    win,
                    frame,
                    surf,
                    x,
                    infrSurfInt[:, it],
                    convSurfInt[:, it],
                    heatSurfInt[:, it],
                    infrWinInt[:, it],
                    convWinInt[:, it],
                    heatWinInt[:, it],
                ) = thermal.main(
                    numWin,
                    numFrame,
                    numSurf,
                    win,
                    frame,
                    surf,
                    config.volume,
                    infSch[it],
                    config.airExch,
                    convIntHeatFlx[it],
                    radIntHeatFlx[it],
                    tempInAirTemp,
                    tempOutAir[it],
                    skyTemp[it],
                    CHTCwin,
                    areaWin,
                    CHTCframe,
                    areaFrame,
                    CHTCsurf,
                    areaSurf,
                    irrWinAbsExt[:, it],
                    irrWinAbsInt[:, it],
                    RHTCwinWin,
                    RHTCwinFrame,
                    RHTCwinSurf,
                    RHTCframeFrame,
                    RHTCframeWin,
                    RHTCframeSurf,
                    RHTCsurfSurf,
                    RHTCsurfWin,
                    RHTCsurfFrame,
                    extRHTC[it],
                    extCHTC[it],
                    irrFrameInt[:, it],
                    irrFrameExt[:, it],
                    irrSurfInt[:, it],
                    irrSurfExt[:, it],
                    bc,
                    calorim[:, 0],
                    conIndex[:, it],
                    heatSpSch[it],
                    coolSpSch[it],
                )
                # Thermal comfort calculations
                if opts.thermComf:
                    RHInTemp = thermalComfort.latEnergyBalance(
                        config.volume,
                        infSch[it],
                        config.airExch,
                        config.powerPeople,
                        RHInTemp / 100.0,
                        tempInAirTemp,
                        relHumOut[it],
                        tempOutAir[it],
                        0,
                        0,
                        p.timeStep,
                    )
            # -------------------------------
            tempInAir[it] = tempInAirTemp
            tempWinInt[:, it] = fromiter((c.intTemp for c in win), float)
            tempWinExt[:, it] = fromiter((c.extTemp for c in win), float)
            tempSurfInt[:, it] = fromiter((c.intTemp for c in surf), float)
            tempSurfExt[:, it] = fromiter((c.extTemp for c in surf), float)
            # -----------------------------------------------------
            # thermal comfort simulation PMV and PPD: ISO7730
            if opts.thermComf:
                RHIn[it] = RHInTemp
                # Ponderated surface temperatures windows/opaque
                tWalls[:, it] = thermalComfort.winAreas(
                    winVect, frameVect, surfVect, tempWinInt[:, it], tempSurfInt[:, it]
                )
                # Mean radiant temperature for every occupant position (5 positions)
                tRad[:, it] = thermalComfort.tempRadSimp(
                    tWalls[:, it],
                    config.length,
                    config.width,
                    config.height,
                    config.thickSouth,
                    config.thickEast,
                    config.thickNorth,
                    config.thickWest,
                    config.thickFloor,
                    config.thickCeiling,
                )
                # Predicted mean value
                pmv[:, it] = thermalComfort.pmvSimp(
                    tempInAir[it] - 273.15, tRad[:, it] - 273.15, RHIn[it]
                )
                # Predicted percentage of dissatisfied
                ppd[:, it] = thermalComfort.ppdSimp(
                    tempInAir[it] - 273.15, tRad[:, it] - 273.15, RHIn[it]
                )
    # ----------------------------------------------------
    # daylight autonomy and annual (direct) sunlight exposure
    if opts.da:
        da = solar.daCalc(goodHour, numSensor, totWorkHour)
    # ----------------------------------------------------
    # ----------------------------------------------------
    # ---------------------------------------------------------
    # write outputs
    print("Time span (s): ", time.time() - start_time)
    print("Writing outputs...")
    savetxt("%s/hh.out" % config.output_path, hh, delimiter=",", fmt="%1.1f")
    if opts.daylight:
        savetxt(
            "%s/ill.out" % config.output_path,
            transpose(ill),
            delimiter=",",
            fmt="%1.1f",
        )
        savetxt(
            "%s/illTrans.out" % config.output_path, illTrans, delimiter=",", fmt="%1.1f"
        )  # lux
    if opts.rad == "radTpm":
        savetxt(
            "%s/illRadTpm.out" % config.output_path,
            transpose(illRadTpm),
            delimiter=",",
            fmt="%1.1f",
        )
    if opts.lightSch or opts.therm:
        savetxt(
            "%s/lightSch.out" % config.output_path,
            transpose(lightSch),
            delimiter=",",
            fmt="%1.2f",
        )
        savetxt(
            "%s/lights.dat" % config.output_path,
            transpose(lightSchYear),
            delimiter=",",
            fmt="%1.2f",
        )
        savetxt(
            "%s/lightCons.out" % config.output_path,
            transpose(lightCons),
            delimiter=",",
            fmt="%1.2f",
        )  # W
    if opts.cutOff or opts.mtxCntrl or opts.refeedCntrl or opts.schCntrl:
        savetxt(
            "%s/conIndex.out" % config.output_path,
            transpose(conIndex),
            delimiter=",",
            fmt="%i",
        )
    if opts.df:
        savetxt(
            "%s/df.out" % config.output_path, transpose(df), delimiter=",", fmt="%1.5f"
        )
    if opts.da:
        savetxt(
            "%s/da.out" % config.output_path, transpose(da), delimiter=",", fmt="%1.1f"
        )
        savetxt(
            "%s/directSunHour.out" % config.output_path,
            transpose(directSunHour),
            delimiter=",",
            fmt="%1.1f",
        )
    if opts.glare or opts.glareSimpl:
        savetxt(
            "%s/illVert.out" % config.output_path,
            transpose(illVert),
            delimiter=",",
            fmt="%1.1f",
        )
        savetxt(
            "%s/dgp.out" % config.output_path,
            transpose(dgp),
            delimiter=",",
            fmt="%1.2f",
        )
    if opts.glareFull:
        savetxt(
            "%s/dgpFull.out" % config.output_path,
            transpose(dgpFull),
            delimiter=",",
            fmt="%1.2f",
        )
    if opts.therm:
        savetxt(
            "%s/uWind.out" % config.output_path, uWind, delimiter=",", fmt="%1.1f"
        )  # m/s
        savetxt(
            "%s/extCHTC.out" % config.output_path, extCHTC, delimiter=",", fmt="%1.1f"
        )  # W/m2K
        savetxt(
            "%s/tempOutAir.out" % config.output_path,
            tempOutAir,
            delimiter=",",
            fmt="%1.1f",
        )  # K
        savetxt(
            "%s/solTrans.out" % config.output_path, solTrans, delimiter=",", fmt="%1.1f"
        )  # W
        savetxt(
            "%s/convIntHeatFlx.out" % config.output_path,
            convIntHeatFlx,
            delimiter=",",
            fmt="%1.1f",
        )  # W
        savetxt(
            "%s/radIntHeatFlx.out" % config.output_path,
            radIntHeatFlx,
            delimiter=",",
            fmt="%1.1f",
        )  # W
        savetxt(
            "%s/irrSurfExt.out" % config.output_path,
            transpose(irrSurfExt),
            delimiter=",",
            fmt="%1.1f",
        )  # W/m2
        savetxt(
            "%s/irrSurfInt.out" % config.output_path,
            transpose(irrSurfInt),
            delimiter=",",
            fmt="%1.1f",
        )  # W/m2
        savetxt(
            "%s/infrSurfInt.out" % config.output_path,
            transpose(infrSurfInt),
            delimiter=",",
            fmt="%1.1f",
        )  # W
        savetxt(
            "%s/convSurfInt.out" % config.output_path,
            transpose(convSurfInt),
            delimiter=",",
            fmt="%1.1f",
        )  # W
        savetxt(
            "%s/heatSurfInt.out" % config.output_path,
            transpose(heatSurfInt),
            delimiter=",",
            fmt="%1.1f",
        )  # W
        savetxt(
            "%s/tempSurfInt.out" % config.output_path,
            transpose(tempSurfInt),
            delimiter=",",
            fmt="%1.1f",
        )  # K
        savetxt(
            "%s/tempSurfExt.out" % config.output_path,
            transpose(tempSurfExt),
            delimiter=",",
            fmt="%1.1f",
        )  # W/m2
        savetxt(
            "%s/irrFrameExt.out" % config.output_path,
            transpose(irrFrameExt),
            delimiter=",",
            fmt="%1.1f",
        )  # W/m2
        savetxt(
            "%s/irrFrameInt.out" % config.output_path,
            transpose(irrFrameInt),
            delimiter=",",
            fmt="%1.1f",
        )  # W/m2
        savetxt(
            "%s/irrWinAbsExt.out" % config.output_path,
            transpose(irrWinAbsInt),
            delimiter=",",
            fmt="%1.1f",
        )  # W/m2
        savetxt(
            "%s/irrWinAbsInt.out" % config.output_path,
            transpose(irrWinAbsInt),
            delimiter=",",
            fmt="%1.1f",
        )  # W/m2
        savetxt(
            "%s/infrWinInt.out" % config.output_path,
            transpose(infrWinInt),
            delimiter=",",
            fmt="%1.1f",
        )  # W
        savetxt(
            "%s/convWinInt.out" % config.output_path,
            transpose(convWinInt),
            delimiter=",",
            fmt="%1.1f",
        )  # W
        savetxt(
            "%s/heatWinInt.out" % config.output_path,
            transpose(heatWinInt),
            delimiter=",",
            fmt="%1.1f",
        )  # W
        savetxt(
            "%s/tempWinInt.out" % config.output_path,
            transpose(tempWinInt),
            delimiter=",",
            fmt="%1.1f",
        )  # K
        savetxt(
            "%s/tempWinExt.out" % config.output_path,
            transpose(tempWinExt),
            delimiter=",",
            fmt="%1.1f",
        )  # K
        savetxt(
            "%s/effGValue.out" % config.output_path,
            transpose(effGValue),
            delimiter=",",
            fmt="%1.4f",
        )  # -
        savetxt(
            "%s/tempInAir.out" % config.output_path,
            tempInAir,
            delimiter=",",
            fmt="%1.1f",
        )  # K
        savetxt(
            "%s/tempOutAir.out" % config.output_path,
            tempOutAir,
            delimiter=",",
            fmt="%1.1f",
        )  # K
        savetxt(
            "%s/energyDemand.out" % config.output_path,
            energyDemand,
            delimiter=",",
            fmt="%1.1f",
        )  # W
    if opts.thermComf:
        savetxt(
            "%s/pmv.out" % config.output_path,
            transpose(pmv),
            delimiter=",",
            fmt="%1.2f",
        )  # pmv
        savetxt(
            "%s/ppd.out" % config.output_path,
            transpose(ppd),
            delimiter=",",
            fmt="%1.0f",
        )  # ppd
        savetxt(
            "%s/tWalls.out" % config.output_path,
            transpose(tWalls),
            delimiter=",",
            fmt="%1.1f",
        )  # K
        savetxt(
            "%s/tRad.out" % config.output_path,
            transpose(tRad),
            delimiter=",",
            fmt="%1.1f",
        )  # K
        savetxt(
            "%s/RHOut.out" % config.output_path, relHumOut, delimiter=",", fmt="%1.1f"
        )  # %
        savetxt(
            "%s/RHIn.out" % config.output_path, RHIn, delimiter=",", fmt="%1.1f"
        )  # %
    if opts.therm:
        ioRoutines.writeKuhn11Abs(
            config.numConWin, win[0], tmxs, rmxs, calorim, config.output_path
        )
    # ----------------------------------------------------


def create_config(_path):
    with open(_path, "w") as _file:
        print("[PATHS]", file=_file)
        print("meteo_path = your_file", file=_file)
        print("matOpaque = your_file", file=_file)
        print("[VARIABLES]", file=_file)
        print("a = your_value", file=_file)


def _check_options(options):
    assert not (
        (options.mtxCntrl or options.refeedCntrl)
        and (
            (not options.therm)
            or (not options.glare and not options.glareSimpl)
            or not options.daylight
        )
    ), "MtxCntrl must be used combined with a thermal, daylight and glare simulation"

    assert not (
        options.thermComf and not options.therm
    ), "The option thermComf must be used in combination with the option therm"


def _get_config_path(options):
    if opts.config[0] == "/":
        return opts.config
    else:
        return os.getcwd() + "/" + opts.config


def _create_non_existing_directories(config):
    if hasattr(config, "output_path") and not os.path.exists(config.output_path):
        os.makedirs(os.path.join(config.output_path, "figures"))
    if hasattr(config, "workDir_path") and not os.path.exists(config.workDir_path):
        os.makedirs(config.workDir_path)


if __name__ == "__main__":
    try:
        opts = argument_parser.parse()
        _check_options(opts)
        config_path = _get_config_path(opts)
        assert os.path.isfile(config_path), (
            "Config file '%s' not found" % config_path.split("/")[-1]
        )
        config = config_parser.parse(config_path, opts)
        _create_non_existing_directories(config)
        main(opts, config)
    except Exception as error:
        print(
            "".join(
                traceback.format_exception(
                    etype=type(error), value=error, tb=error.__traceback__
                )
            )
        )
        sys.exit("The following error occurred: %s" % error)
