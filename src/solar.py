#!/usr/bin/env python
# -*- coding: utf-8 -*-


from numpy import atleast_2d, genfromtxt, sin, pi, zeros, exp, mean, transpose, ones
import calcs
from ioRoutines import shell
import time
import geometry


def illum(ill, illVmx, tmxv, conIndex, dmx, skyVis, numWin):
    # illuminance calculation
    for i in range(numWin):
        ill = ill + illVmx[i, :, :].dot(
            tmxv[conIndex[i], :, :].dot(dmx[i, :, :].dot(skyVis))
        )
    ill = ill * 179
    return ill


def lightCntrl(lightCntrl, lightSch, floorArea, powerLight, ill):
    # lighting power control
    lightOutput = 0.0
    for i in range(lightCntrl.shape[0]):
        aux = calcs.photocellLightControl(
            lightCntrl[i, 2],
            ill[int(lightCntrl[i, 0])],
            lightCntrl[i, 3],
            lightCntrl[i, 4],
        )
        lightOutput = lightOutput + aux * lightCntrl[i, 1]
    lightSch = lightSch * lightOutput
    lightCons = floorArea * powerLight * lightSch
    # --------------------------------------------------------
    return lightSch, lightCons


def df(ill, illExtDmx, illVmx, tmxv, conIndex, dmx, numWin, lat, lon, tzone, workDir):
    # Generate CIE sky vector
    shell(
        "gensky 6 21 12 -a %f -o %f -m %f -c | genskyvec > %scie.skv"
        % (lat, lon, tzone, workDir)
    )
    skyCie = atleast_2d(
        genfromtxt("%scie.skv" % workDir, usecols=0, skip_header=7, delimiter="\t")
    )
    # interior horizontal illuminance
    illInt = illum(ill, illVmx, tmxv, conIndex, dmx, transpose(skyCie), numWin)
    # exterior horizontal illuminace
    illExt = illExtDmx.dot(transpose(skyCie)) * 179
    # daylight factor
    df = illInt / illExt * 100
    # df = illInt
    # --------------------------------------------------------
    return df


def cie(mm, dd, hh, lat, lon, tzone, irrdirn, irrdiffh, sunAltitude, workDir):
    # Generate CIE sky vector
    shell(
        "gensky %i %i %1.1f -a %f -o %f -m %f +s -g 0.2 -B %f -R %f | genskyvec > %scie.skv"
        % (
            mm,
            dd,
            hh,
            lat,
            lon,
            tzone,
            irrdiffh,
            max(irrdirn * sin(sunAltitude * pi / 180), 0),
            workDir,
        )
    )
    skyCie = atleast_2d(
        genfromtxt("%scie.skv" % workDir, usecols=0, skip_header=7, delimiter="\t")
    )
    return skyCie


def daGoodHour(
    ill,
    goodHour,
    min_illu,
    numSensor,
    start_hour,
    end_hour,
    current_hour,
    current_month,
    totWorkHour,
):
    if start_hour < current_hour < end_hour:
        # if current_month<7 or current_month>8:
        totWorkHour += 1
        for i in range(numSensor):
            if ill[i] > min_illu:
                goodHour[i] += 1
    return goodHour, totWorkHour


def daCalc(goodHour, numSensor, totWorkHour):
    out = zeros(numSensor)
    for i in range(numSensor):
        out[i] = 100 * goodHour[i] / totWorkHour
    return out


def directSunlightHour(ill, directSunHour, numSensor):
    for i in range(numSensor):
        if ill[i] > 1000:  # 1000(lux)
            directSunHour[i] += 1
    return directSunHour


# def glare(nCore,illVert,dgp,illVertVmx,tmxv,conIndex,dmx,skyVis,numWin,mm,dd,hh,irrdirn,irrdiffh,lat,lon,tzone,illuVertPts,numSensorVert,winRad_path,workDir):
#   # glare calculation
#   # vertical illuminance
#   start_time = time.time()
#   for i in range(numWin):
#      illVert = illVert + illVertVmx[i,:,:].dot(tmxv[conIndex[i],:,:].dot(dmx[i,:,:].dot(skyVis)))
#   illVert = illVert*179
#   octFlag = 1
#   # calculation of lowlight correction
#   lowlight = zeros(numSensorVert)
#   for k in range(numSensorVert):
#      if illVert[k]<1000:
#         lowlight[k]=1.0*exp(0.024*illVert[k]-4)/(1+exp(0.024*illVert[k]-4))
#      else:
#         lowlight[k]=1.0
#      if illVert[k]>14310.0:
#         dgp[k] = 1.
#      else:
#         if irrdirn>25:
#            if octFlag:
#               geometry.winRad(numWin,workDir,conIndex)
#               # generate an octree
#               shell("gendaylit %i %i %1.1f -W %1.1f %1.1f -a %f -o %f -m %f -w > %sgendaylit.sky" % (mm,dd,hh,irrdirn,irrdiffh,lat,lon,tzone,workDir))
#               genRtraceOct = "oconv %sroom.mat %sroom.rad"%(workDir,workDir)
#               for i in range(numWin):
#                  genRtraceOct = genRtraceOct + " %s %swindowRad_%i.rad" % (winRad_path[conIndex[i]],workDir,i)
#               genRtraceOct = genRtraceOct+" %sgendaylit.sky %ssky_prev.rad > %srtrace.oct" % (workDir,workDir,workDir)
#               shell(genRtraceOct)
#               octFlag = 0
#            # evalglare
#            #---------------------------
#            #output = shell("vwrays -ff -vf %sview.vf -x 800 -y 800 | rtrace -n %i -ab 0 -dj 0 -dr 2 -lr 4 -lw 0.000001 -st 0.15 -ffc `vwrays -vf %sview.vf -x 800 -y 800 -d` %srtrace.oct | evalglare -i %.3f -1 -vta -vh 180 -vv 180 -vd %1.3f %1.3f %1.3f -vu 0 0 1" % (workDir,nCore,workDir,workDir,illVert[k],illuVertPts[k,3],illuVertPts[k,4],illuVertPts[k,5]),True)
#            #---------------------------
#            #shell("vwrays -ff -vf %sview.vf -x 800 -y 800 | rtrace -n %i -ab 0 -dj 0 -dr 2 -lr 4 -lw 0.000001 -st 0.15 -ffc `vwrays -vf %sview.vf -x 800 -y 800 -d` %srtrace.oct > %srtrace.hdr" % (workDir,nCore,workDir,workDir,workDir))
#            #output = shell("evalglare -i %.3f -1 -vta -vh 180 -vv 180 -vd %1.3f %1.3f %1.3f -vu 0 0 1 %srtrace.hdr" % (mean(illVert),illuVertPts[k,3],illuVertPts[k,4],illuVertPts[k,5],workDir),True)
#            #---------------------------
#            output = shell("rpict -x 800 -y 800 -ab 0 -dj 0 -dr 2 -lr 4 -lw 0.000001  -st 0.15 -ps 0 -pj 0 -vta -vv 180 -vh 180 -vu 0 0 1 -vs 0 -vl 0 -vp %1.3f %1.3f %1.3f -vd %1.3f %1.3f %1.3f %srtrace.oct | evalglare -b 8000 -i %1.3f -D -y -1" % (illuVertPts[k,0],illuVertPts[k,1],illuVertPts[k,2],illuVertPts[k,3],illuVertPts[k,4],illuVertPts[k,5],workDir,illVert[k]),True)
#            #---------------------------
#            dgp[k] = float(output[0])
#         else:
#            dgp[k] = lowlight[k]*(6.22e-5*illVert[k]+0.184)
#      if dgp[k] > 1.: dgp[k] = 1.0
#   # --------------------------------------------------------
#   print('Time span (s): ',time.time()-start_time)
#   # --------------------------------------------------------
#   return illVert,dgp


def glare(
    nCore,
    illVert,
    dgp,
    illVertVmx,
    tmxv,
    conIndex,
    dmx,
    skyVis,
    numWin,
    mm,
    dd,
    hh,
    irrdirn,
    irrdiffh,
    lat,
    lon,
    tzone,
    illuVertPts,
    numSensorVert,
    winRad_path,
    workDir,
    inputDir,
):
    # glare calculation
    # vertical illuminance
    start_time = time.time()
    for i in range(numWin):
        illVert = illVert + illVertVmx[i, :, :].dot(
            tmxv[conIndex[i], :, :].dot(dmx[i, :, :].dot(skyVis))
        )
    illVert = illVert * 179
    octFlag = True
    # calculation of lowlight correction
    lowlight = zeros(numSensorVert)
    for k in range(numSensorVert):
        if illVert[k] < 1000:
            lowlight[k] = (
                1.0 * exp(0.024 * illVert[k] - 4) / (1 + exp(0.024 * illVert[k] - 4))
            )
        else:
            lowlight[k] = 1.0
        if illVert[k] > 14310.0:
            dgp[k] = 1.0
        else:
            if irrdirn > 25:
                if octFlag:
                    geometry.winRad(numWin, workDir, inputDir, conIndex)
                    # generate an octree
                    shell(
                        "gendaylit %i %i %1.1f -W %1.1f %1.1f -a %f -o %f -m %f -w > %sgendaylit.sky"
                        % (mm, dd, hh, irrdirn, irrdiffh, lat, lon, tzone, inputDir)
                    )
                    genRtraceOct = "oconv %sroom.mat %sroom.rad" % (workDir, workDir)
                    for i in range(numWin):
                        genRtraceOct = genRtraceOct + " %s %swindowRad_%i.rad" % (
                            winRad_path[conIndex[i]],
                            inputDir,
                            i,
                        )
                    genRtraceOct = (
                        genRtraceOct
                        + " %sgendaylit.sky %ssky_prev.rad > %srtrace.oct"
                        % (inputDir, workDir, inputDir)
                    )
                    shell(genRtraceOct)
                    octFlag = False
                # evalglare
                # ---------------------------
                # output = shell("vwrays -ff -vf %sview.vf -x 800 -y 800 | rtrace -n %i -ab 0 -dj 0 -dr 2 -lr 4 -lw 0.000001 -st 0.15 -ffc `vwrays -vf %sview.vf -x 800 -y 800 -d` %srtrace.oct | evalglare -i %.3f -1 -vta -vh 180 -vv 180 -vd %1.3f %1.3f %1.3f -vu 0 0 1" % (workDir,nCore,workDir,workDir,illVert[k],illuVertPts[k,3],illuVertPts[k,4],illuVertPts[k,5]),True)
                # ---------------------------
                # shell("vwrays -ff -vf %sview.vf -x 800 -y 800 | rtrace -n %i -ab 0 -dj 0 -dr 2 -lr 4 -lw 0.000001 -st 0.15 -ffc `vwrays -vf %sview.vf -x 800 -y 800 -d` %srtrace.oct > %srtrace.hdr" % (workDir,nCore,workDir,workDir,workDir))
                # output = shell("evalglare -i %.3f -1 -vta -vh 180 -vv 180 -vd %1.3f %1.3f %1.3f -vu 0 0 1 %srtrace.hdr" % (mean(illVert),illuVertPts[k,3],illuVertPts[k,4],illuVertPts[k,5],workDir),True)
                # ---------------------------
                # ~ output = shell("rpict -x 800 -y 800 -ab 1 -dj 0 -dr 2 -lr 4 -lw 0.000001  -st 0.15 -ps 0 -pj 0 -vta -vv 180 -vh 180 -vu 0 0 1 -vs 0 -vl 0 -vp %1.3f %1.3f %1.3f -vd %1.3f %1.3f %1.3f %srtrace.oct | evalglare -b 8000 -i %1.3f -D -y -1" % (illuVertPts[k,0],illuVertPts[k,1],illuVertPts[k,2],illuVertPts[k,3],illuVertPts[k,4],illuVertPts[k,5],inputDir,illVert[k]),True)
                # ~ #---------------------------
                output = shell(
                    "rpict -x 800 -y 800 -ab 0 -dj 0 -dr 2 -lr 4 -lw 0.000001  -st 0.15 -ps 0 -pj 0 -vta -vv 180 -vh 180 -vu 0 0 1 -vs 0 -vl 0 -vp %1.3f %1.3f %1.3f -vd %1.3f %1.3f %1.3f %srtrace.oct | evalglare -b 8000 -i %1.3f -D -y -1"
                    % (
                        illuVertPts[k, 0],
                        illuVertPts[k, 1],
                        illuVertPts[k, 2],
                        illuVertPts[k, 3],
                        illuVertPts[k, 4],
                        illuVertPts[k, 5],
                        inputDir,
                        illVert[k],
                    ),
                    True,
                )
                # ---------------------------
                dgp[k] = float(output)
            else:
                dgp[k] = lowlight[k] * (6.22e-5 * illVert[k] + 0.184)
        if dgp[k] > 1.0:
            dgp[k] = 1.0
    # --------------------------------------------------------
    print("Time span (s): ", time.time() - start_time)
    # --------------------------------------------------------
    return illVert, dgp


def glareFull(
    nCore,
    dgp,
    conIndex,
    numWin,
    mm,
    dd,
    hh,
    irrdirn,
    irrdiffh,
    lat,
    lon,
    tzone,
    illuVertPts,
    numSensorVert,
    winRad_path,
    workDir,
    inputDir,
):
    # glare calculation
    # vertical illuminance
    start_time = time.time()
    octFlag = True
    # calculation of lowlight correction
    for k in range(numSensorVert):
        if octFlag:
            geometry.winRad(numWin, workDir, inputDir, conIndex)
            # generate an octree
            shell(
                "gendaylit %i %i %1.1f -W %1.1f %1.1f -a %f -o %f -m %f -w > %sgendaylit.sky"
                % (mm, dd, hh, irrdirn, irrdiffh, lat, lon, tzone, inputDir)
            )
            genRtraceOct = "oconv %sroom.mat %sroom.rad" % (workDir, workDir)
            for i in range(numWin):
                genRtraceOct = genRtraceOct + " %s %swindowRad_%i.rad" % (
                    winRad_path[conIndex[i]],
                    inputDir,
                    i,
                )
            genRtraceOct = (
                genRtraceOct
                + " %sgendaylit.sky %ssky_prev.rad > %srtrace.oct"
                % (inputDir, workDir, inputDir)
            )
            shell(genRtraceOct)
            octFlag = False
        # evalglare
        # ---------------------------
        output = shell(
            "vwrays -ff -vf %sview.vf -x 800 -y 800 | rtrace -n %i -ab 5 -ad 512 -ar 128 -as 256 -aa 0.15 -dc 1 -ffc `vwrays -vf %sview.vf -x 800 -y 800 -d` %srtrace.oct | evalglare -vta -vh 180 -vv 180 -vd %1.3f %1.3f %1.3f -vu 0 0 1"
            % (
                workDir,
                nCore,
                inputDir,
                inputDir,
                illuVertPts[k, 3],
                illuVertPts[k, 4],
                illuVertPts[k, 5],
            ),
            True,
        )
        # ---------------------------
        dgp[k] = float(output[0].split(" ")[1])
        if dgp[k] > 1.0:
            dgp[k] = 1.0
    # --------------------------------------------------------
    print("Time span (s): ", time.time() - start_time)
    # --------------------------------------------------------
    return dgp


def glareSimpl(
    illVert, dgp, illVertVmx, tmxv, conIndex, dmx, skyVis, numWin, numSensorVert
):
    # Simplified DGP method (Wienold 2009)
    # vertical illuminance
    # start_time = time.time()
    for i in range(numWin):
        illVert = illVert + illVertVmx[i, :, :].dot(
            tmxv[conIndex[i], :, :].dot(dmx[i, :, :].dot(skyVis))
        )
    illVert = illVert * 179
    for k in range(numSensorVert):
        dgp[k] = min(6.22e-5 * mean(illVert[k]) + 0.184, 1.0)
    # --------------------------------------------------------
    return illVert, dgp


def glare_M(
    nCore,
    illVert,
    dgp,
    illVertVmx,
    tmxv,
    conIndex,
    dmx,
    skyVis,
    numWin,
    mm,
    dd,
    hh,
    irrdirn,
    irrdiffh,
    lat,
    lon,
    tzone,
    illuVertPts,
    numSensorVert,
    winRad_path,
    workDir,
    inputDir,
):
    # glare calculation
    # vertical illuminance
    dgpFlag = False
    for i in range(numWin):
        illVert = illVert + illVertVmx[i, :, :].dot(
            tmxv[conIndex[i], :, :].dot(dmx[i, :, :].dot(skyVis))
        )
    illVert = illVert * 179
    # octFlag = 1
    # calculation of lowlight correction
    lowlight = zeros(numSensorVert)
    for k in range(numSensorVert):
        if illVert[k] < 1000:
            lowlight[k] = (
                1.0 * exp(0.024 * illVert[k] - 4) / (1 + exp(0.024 * illVert[k] - 4))
            )
        else:
            lowlight[k] = 1.0
        if illVert[k] > 14310.0:
            dgp[k] = 1.0
        else:
            if irrdirn > 25:
                # if octFlag:
                # geometry.winRad(numWin,workDir,conIndex)
                dgpFlag = True
            else:
                dgp[k] = lowlight[k] * (6.22e-5 * illVert[k] + 0.184)
        if dgp[k] > 1.0:
            dgp[k] = 1.0
    # --------------------------------------------------------
    return illVert, dgp, dgpFlag


def workerstar(arg):
    return do_work(*arg)


def do_work(
    nCore,
    illVert,
    dgp,
    illVertVmx,
    tmxv,
    conIndex,
    dmx,
    skyVis,
    numWin,
    mm,
    dd,
    hh,
    currenthh,
    irrdirn,
    irrdiffh,
    lat,
    lon,
    tzone,
    illuVertPts,
    numSensorVert,
    winRad_path,
    workDir,
    inputDir,
):
    # geometry.winRad(numWin,workDir,conIndex)
    shell(
        "gendaylit %i %i %1.1f -W %1.1f %1.1f -a %f -o %f -m %f -w > %sgendaylit_%i.sky"
        % (mm, dd, hh, irrdirn, irrdiffh, lat, lon, tzone, inputDir, currenthh)
    )
    genRtraceOct = "oconv %sroom.mat %sroom.rad " % (workDir, workDir)
    for i in range(numWin):
        genRtraceOct = genRtraceOct + " %s %swindowRad_%i.rad " % (
            winRad_path[conIndex[i]],
            inputDir,
            i,
        )
    genRtraceOct = (
        genRtraceOct
        + " %sgendaylit_%i.sky %ssky_prev.rad > %srtrace_%i.oct"
        % (inputDir, currenthh, workDir, inputDir, currenthh)
    )
    shell(genRtraceOct)
    for k in range(numSensorVert):
        output = shell(
            "rpict -x 800 -y 800 -ab 0 -dj 0 -dr 2 -lr 4 -lw 0.000001  -st 0.15 -ps 0 -pj 0 -vta -vv 180 -vh 180 -vu 0 0 1 -vs 0 -vl 0 -vp %1.3f %1.3f %1.3f -vd %1.3f %1.3f %1.3f %srtrace_%i.oct | evalglare -b 8000 -i %1.3f -D -y -1"
            % (
                illuVertPts[k, 0],
                illuVertPts[k, 1],
                illuVertPts[k, 2],
                illuVertPts[k, 3],
                illuVertPts[k, 4],
                illuVertPts[k, 5],
                inputDir,
                currenthh,
                illVert[k],
            ),
            True,
        )
        # ---------------------------
        dgp[k] = float(output[0])
    return dgp


def isoWinAbs(numWin, numConWin, numPaneConWin, conWin, conIndex, dmx, skySol):
    #  window layer solar absorption
    irrWinLay = zeros([numWin, max(numPaneConWin)])
    for i in range(numWin):
        for k in range(conWin[conIndex[i]].numPane):
            irrWinLay[i, k] = (
                transpose(conWin[conIndex[i]].absFront[k, 0, :])
                .dot(dmx[i, :, :].dot(skySol))
                .reshape(-1)
                / 145
                * pi
            )
    return irrWinLay


def calorimWinAbs(numWin, win, dmx, skySol):
    #  window layer solar absorption
    irrWinAbsExt = zeros(numWin)
    irrWinAbsInt = zeros(numWin)
    for i in range(numWin):
        irrWinAbsInt[i] = (
            transpose(win[i].absLayInt).dot(dmx[i, :, :].dot(skySol)).reshape(-1)
            / 145
            * pi
        )
        irrWinAbsExt[i] = (
            transpose(win[i].absLayExt).dot(dmx[i, :, :].dot(skySol)).reshape(-1)
            / 145
            * pi
        )
    return irrWinAbsExt, irrWinAbsInt


def effGVal(numWin, win, dmx, skySol, calorim, conIndex):
    #  effective g-value
    effGValue = zeros(numWin)
    for i in range(numWin):
        gValue = zeros(145)
        unitVec = ones(145)
        for j in range(145):
            gValue[j] = calorim[conIndex[i], j + 1]
        effGValue[i] = transpose(gValue).dot(dmx[i, :, :].dot(skySol)).reshape(
            -1
        ) / transpose(unitVec).dot(dmx[i, :, :].dot(skySol)).reshape(-1)
    return effGValue


def irrFadExt(numSurf, numSensorSurf, irrSurfExtDmx, skySol, bc):
    irrSurfExtTmp = zeros([numSurf, max(numSensorSurf)])
    irrSurfExt = zeros(numSurf)
    # exterior solar radiation
    for k in range(numSurf):
        if bc[k]:
            irrSurfExtTmp[k, :] = irrSurfExtDmx[k, :, :].dot(skySol).reshape(-1)
            for i in range(numSensorSurf[k]):
                irrSurfExt[k] = irrSurfExt[k] + irrSurfExtTmp[k, i] / numSensorSurf[k]
    return irrSurfExt


def illumTrans(
    numWin, numSensorWin, areaWin, irrWinBackVmx, tmxv, conIndex, dmx, skyVis, illTrans
):
    irrWinBack = zeros([numWin, max(numSensorWin)])
    for i in range(numWin):
        irrWinBack[i, :] = irrWinBack[i, :] + irrWinBackVmx[i, :, :].dot(
            tmxv[conIndex[i], :, :].dot(dmx[i, :, :].dot(skyVis))
        ).reshape(-1)
        for j in range(numSensorWin[i]):
            illTrans = illTrans + irrWinBack[i, j] / numSensorWin[i]
    return illTrans * 179


def therm(
    skySol,
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
    conIndex,
    irrSurfInt,
    irrSurfExt,
    irrFrameInt,
    irrFrameExt,
    solTrans,
):
    # initialize local variables
    irrSurfExtTmp = zeros([numSurf, max(numSensorSurf)])
    irrSurfExt = zeros(numSurf)

    irrSurfFront = zeros([numSurf, max(numSensorSurf)])
    irrSurfBack = zeros([numSurf, max(numSensorSurf)])

    irrFrameExtTmp = zeros([numFrame, max(numSensorFrame)])
    irrFrameExt = zeros(numFrame)

    irrFrameFront = zeros([numFrame, max(numSensorFrame)])
    # irrFrameBack = zeros([numFrame,max(numSensorFrame)])

    irrWinFront = zeros([numWin, max(numSensorWin)])
    irrWinBack = zeros([numWin, max(numSensorWin)])

    #   irrWinAbsExt = zeros(numWin)
    #   irrWinAbsInt = zeros(numWin)
    irrWinInt = zeros(numWin)

    irrRevFront = zeros([numFrame, max(numSensorRev)])
    irrRevBack = zeros([numFrame, max(numSensorRev)])

    # exterior solar radiation
    for k in range(numSurf):
        if bc[k]:
            irrSurfExtTmp[k, :] = irrSurfExtDmx[k, :, :].dot(skySol).reshape(-1)
            for i in range(numSensorSurf[k]):
                irrSurfExt[k] = irrSurfExt[k] + irrSurfExtTmp[k, i] / numSensorSurf[k]
    for k in range(numFrame):
        irrFrameExtTmp[k, :] = irrFrameExtDmx[k, :, :].dot(skySol).reshape(-1)
        # print(irrFrameExtTmp[k,:])
        for i in range(numSensorFrame[k]):
            irrFrameExt[k] = irrFrameExt[k] + irrFrameExtTmp[k, i] / numSensorFrame[k]
    # interior solar radiation
    for i in range(numWin):
        #      irrWinAbsInt[i] = transpose(win[i].absLayInt).dot(dmx[i,:,:].dot(skySol)).reshape(-1)
        #      irrWinAbsExt[i] = transpose(win[i].absLayExt).dot(dmx[i,:,:].dot(skySol)).reshape(-1)
        irrWinFront[i, :] = irrWinFront[i, :] + irrWinFrontVmx[i, :, :].dot(
            tmxs[conIndex[i], :, :].dot(dmx[i, :, :].dot(skySol))
        ).reshape(-1)
        irrWinBack[i, :] = irrWinBack[i, :] + irrWinBackVmx[i, :, :].dot(
            tmxs[conIndex[i], :, :].dot(dmx[i, :, :].dot(skySol))
        ).reshape(-1)
        for j in range(numSensorWin[i]):
            solTrans = solTrans + irrWinBack[i, j] / numSensorWin[i] * areaWin[i]
            irrWinInt[i] = irrWinInt[i] + irrWinFront[i, j] / numSensorWin[i]
        for k in range(numSurf):
            irrSurfFront[k, :] = irrSurfFront[k, :] + irrSurfFrontVmx[i, k, :, :].dot(
                tmxs[conIndex[i], :, :].dot(dmx[i, :, :].dot(skySol))
            ).reshape(-1)
            irrSurfBack[k, :] = irrSurfBack[k, :] + irrSurfBackVmx[i, k, :, :].dot(
                tmxs[conIndex[i], :, :].dot(dmx[i, :, :].dot(skySol))
            ).reshape(-1)
        for k in range(numFrame):
            irrRevFront[k, :] = irrRevFront[k, :] + irrRevFrontVmx[i, k, :, :].dot(
                tmxs[conIndex[i], :, :].dot(dmx[i, :, :].dot(skySol))
            ).reshape(-1)
            irrRevBack[k, :] = irrRevBack[k, :] + irrRevBackVmx[i, k, :, :].dot(
                tmxs[conIndex[i], :, :].dot(dmx[i, :, :].dot(skySol))
            ).reshape(-1)
            irrFrameFront[k, :] = irrFrameFront[k, :] + irrFrameFrontVmx[
                i, k, :, :
            ].dot(tmxs[conIndex[i], :, :].dot(dmx[i, :, :].dot(skySol))).reshape(-1)
            # irrFrameBack[k,:] = irrFrameBack[k,:] + irrFrameBackVmx[i,k,:,:].dot(tmxs[conIndex[i],:,:].dot(dmx[i,:,:].dot(skySol))).reshape(-1)
    for ks in range(numSurf):
        for j in range(numSensorSurf[ks]):
            irrSurfInt[ks] = (
                irrSurfInt[ks]
                + (irrSurfFront[ks, j] - irrSurfBack[ks, j]) / numSensorSurf[ks]
            )
            # irrSurfInt[ks] = irrSurfInt[ks] + irrSurfFront[ks,j]/numSensorSurf[ks]
        for kf in range(numFrame):
            if frameSurf[kf] == ks:
                for j in range(numSensorRev[kf]):
                    irrSurfInt[ks] = (
                        irrSurfInt[ks]
                        + (irrRevFront[kf, j] - irrRevBack[kf, j])
                        / numSensorRev[kf]
                        * areaRev[kf]
                        / areaSurf[ks]
                    )
                    # irrSurfInt[ks] = irrSurfInt[ks] + irrRevFront[kf,j]/numSensorRev[kf]*areaRev[kf]/areaSurf[ks]
    for kf in range(numFrame):
        for j in range(numSensorFrame[kf]):
            # irrFrameInt[kf] = irrFrameInt[kf] + (irrFrameFront[kf,j]-irrFrameBack[kf,j])/numSensorFrame[kf]
            irrFrameInt[kf] = (
                irrFrameInt[kf] + irrFrameFront[kf, j] / numSensorFrame[kf]
            )
    # --------------------------------------------------------
    return solTrans, irrSurfExt, irrSurfInt, irrFrameExt, irrFrameInt, irrWinInt
