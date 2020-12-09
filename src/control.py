#!/usr/bin/env python
# -*- coding: utf-8 -*-


from numpy import pi, degrees, arccos, mean, radians, zeros, arcsin, cos, matrix
import calcs
import solar


def takeClosest(num,collection):
   return min(collection,key=lambda x:abs(x-num))


def bestIlluminance(numConWin, illVmx, illVertVmx, tmxv, dmx, skyVis, numWin,spHorIllu,spVerIllu,
        numSensor,numSensorVert,numWinLowerPartition,numWinUpperPartition):
   illHor = zeros(numConWin*numConWin)
   illVer = zeros(numConWin*numConWin)
   conIndexLowerPartition = []
   conIndexUpperPartition = []
   conIndexTemp = zeros(numWinLowerPartition+numWinUpperPartition,int)
   for i in range(numConWin):
      for j in range(numConWin):
         for n in range(numWinLowerPartition):
            conIndexTemp[n] = i
         for m in range(numWinUpperPartition):
            conIndexTemp[numWinLowerPartition+m] = j
         k = i*numConWin+j
         # Horizontal illuminance
         illHorAux = matrix(zeros([numSensor, 1]))
         illVerAux = matrix(zeros([numSensorVert, 1]))
         illHorAux = solar.illum(
                illHorAux, 
                illVmx, 
                tmxv, 
                conIndexTemp, 
                dmx,  
                skyVis, 
                numWin
            )
         # Vertical illuminance
         illVerAux = solar.illum(
                illVerAux, 
                illVertVmx, 
                tmxv, 
                conIndexTemp, 
                dmx,  
                skyVis, 
                numWin
            )
         illHor[k] = mean(illHorAux)
         illVer[k] = max(illVerAux)
         if illHor[k] > spHorIllu and illVer[k] < spVerIllu:
            conIndexLowerPartition.append(i)
            conIndexUpperPartition.append(j)
   return conIndexLowerPartition,conIndexUpperPartition

def optState(conIndexPrev,numConWin, illVmx, illVertVmx, tmxv, dmx, skyVis, numWin,spHorIllu,spVerIllu,
   numSensor,numSensorVert,numWinLowerPartition,numWinUpperPartition,cntrlOpt):
    conIndex = zeros(len(conIndexPrev),int)
    if (numWinLowerPartition+numWinUpperPartition) != len(conIndex):
       raise TypeError("The number of windows in partitions do not agree with the total number of windows")
    conIndexLowerPartition = []
    while not conIndexLowerPartition:
        conIndexLowerPartition,conIndexUpperPartition = bestIlluminance(
            numConWin, 
            illVmx, 
            illVertVmx, 
            tmxv, 
            dmx, 
            skyVis, 
            numWin,
            spHorIllu,
            spVerIllu,
            numSensor,
            numSensorVert,
            numWinLowerPartition,
            numWinUpperPartition,
            )
        spHorIllu = spHorIllu -100
        spVerIllu = spVerIllu + 100
    # ----------------------------------------
    if cntrlOpt == 0:
        chosenState = takeClosest(conIndexPrev[0],conIndexUpperPartition)
    # ----------------------------------------
    elif cntrlOpt == 1:
        chosenState = min(conIndexUpperPartition)
    # ----------------------------------------
    elif cntrlOpt == 2:
        chosenState = max(conIndexUpperPartition)
    # ----------------------------------------
    if chosenState in conIndexUpperPartition:
        chosenState_index = conIndexUpperPartition.index(chosenState)
    for n in range(numWinLowerPartition):
        conIndex[n] = conIndexLowerPartition[chosenState_index]
    for m in range(numWinUpperPartition):
        conIndex[numWinLowerPartition+m] = conIndexUpperPartition[chosenState_index]   
    print(conIndex,spHorIllu,spVerIllu)
    return conIndex  


def cutOffAngle(profAng, incAng, width, distance, vertical, slatHor):
    # Cut Off Angle
    # Angles in radians
    if vertical:
        if slatHor:
            cutOffAng = arcsin(distance * cos(profAng) / width) - profAng
        else:
            if profAng <= arccos(width / (width + distance)) and profAng >= -(
                arccos(width / (width + distance))
            ):
                cutOffAng = pi / 2.0
            elif profAng >= pi - abs(arccos(width / (width + distance))):
                cutOffAng = 0.0
            elif profAng <= -pi + abs(arccos(width / (width + distance))):
                cutOffAng = pi
            else:
                cutOffAng = arcsin((width + distance) * cos(profAng) / width) - profAng
    else:
        cutOffAng = arcsin((width + distance) * cos(incAng) / width) - incAng - pi / 2.0
    return cutOffAng


# ~ def mtxCntrl_prev(cntrlMtx,ill,dgp,tempInAir,illSp,thermSp,occ,sunAltitude,conDef):
# ~ if occ == 0. or sunAltitude<=0.:
# ~ conIndex = conDef
# ~ else:
# ~ if ill>illSp and tempInAir<thermSp and dgp < 0.45:
# ~ conIndex = cntrlMtx[0,3]
# ~ if ill<=illSp and tempInAir<thermSp and dgp < 0.45:
# ~ conIndex = cntrlMtx[1,3]
# ~ if ill>illSp and tempInAir>=thermSp and dgp < 0.45:
# ~ conIndex = cntrlMtx[2,3]
# ~ if ill>illSp and tempInAir<thermSp and dgp >= 0.45:
# ~ conIndex = cntrlMtx[3,3]
# ~ if ill<=illSp and tempInAir>=thermSp and dgp < 0.45:
# ~ conIndex = cntrlMtx[4,3]
# ~ if ill>illSp and tempInAir>=thermSp and dgp >= 0.45:
# ~ conIndex = cntrlMtx[5,3]
# ~ if ill<=illSp and tempInAir<thermSp and dgp >= 0.45:
# ~ conIndex = cntrlMtx[6,3]
# ~ if ill<=illSp and tempInAir>=thermSp and dgp >= 0.45:
# ~ conIndex = cntrlMtx[7,3]
# ~ return conIndex


def mtxCntrl(
    cntrlMtx,
    conDef,
    occ,
    sunAltitude,
    ill,
    dgp,
    tempInAir,
    irrSurfExtSouth,
    tempOut24,
    incAng,
    profAng,
    conIndexPrev,
):
    # -------------------------------------------------------
    stateFlag = zeros(len(cntrlMtx.mtx[0, :]))
    for j in range(len(cntrlMtx.mtx[0, :])):
        # -------------------------------------------------------
        # assign values to variables
        if cntrlMtx.var[j] == 0:  # occupation
            cntrlMtx.value[j] = occ
        if cntrlMtx.var[j] == 1:  # daytime
            cntrlMtx.value[j] = sunAltitude
        if cntrlMtx.var[j] == 2:  # average workplane illuminance
            cntrlMtx.value[j] = ill
        if cntrlMtx.var[j] == 3:  # dgp
            cntrlMtx.value[j] = dgp
        if cntrlMtx.var[j] == 4:  # indoor air temperature
            cntrlMtx.value[j] = tempInAir
        if cntrlMtx.var[j] == 5:  # average 24-hour temperature
            cntrlMtx.value[j] = tempOut24
        if cntrlMtx.var[j] == 6:  # exterior irradiance on south facade
            cntrlMtx.value[j] = irrSurfExtSouth
        if cntrlMtx.var[j] == 7:  # incidence angle
            cntrlMtx.value[j] = degrees(incAng)
        if cntrlMtx.var[j] == 8:  # profile angle
            cntrlMtx.value[j] = degrees(profAng)
        # -------------------------------------------------------
        # generate the state line
        if cntrlMtx.value[j] >= cntrlMtx.sp[j]:
            stateFlag[j] = 1
    # -------------------------------------------------------
    # find the right state line in the state matrix
    conIndex = zeros(len(conIndexPrev),int)
    for i in range(len(cntrlMtx.mtx[:, 0])):
        lineFlag = False
        for j in range(len(cntrlMtx.mtx[0, :])):
            if stateFlag[j] == cntrlMtx.mtx[i, j]:
                lineFlag = True
            else:
                lineFlag = False
                break
        if lineFlag:
            if len(cntrlMtx.state[0, :]) == 1:
                for k in range(len(conIndexPrev)):
                    conIndex[k] = cntrlMtx.state[i]
                break
            elif len(cntrlMtx.state[0, :]) > 1:
                for k in range(len(conIndexPrev)):
                    conIndex[k] = cntrlMtx.state[i, k]
                break
            else:
                print("ERROR: Check control definition")
        else:
            for k in range(len(conIndexPrev)):
                conIndex[k] = conDef
        # ~ if lineFlag:
        # ~ for k in range(len(conIndex)):
        # ~ conIndex[k] = cntrlMtx.state[i,k]
        # ~ break
        # ~ else:
        # ~ for k in range(len(conIndex)):
        # ~ conIndex[k] = conDef
    # -------------------------------------------------------
    # ~ print ('conIndex',conIndex)
    return conIndex


def mtxCntrlOptState(cntrlMtx,
        conDef,
        occ,
        sunAltitude,
        ill,
        dgp,
        tempInAir,
        irrSurfExtSouth,
        tempOut24,
        incAng,
        profAng,
        conIndexPrev,
        numConWin, 
        illVmx, 
        illVertVmx, 
        tmxv, 
        dmx, 
        skyVis, 
        numWin,
        spHorIllu,
        spVerIllu,
        numSensor,
        numSensorVert,
        numWinLowerPartition,
        numWinUpperPartition,
        cntrlOpt
    ):
    # -------------------------------------------------------
    conIndex = mtxCntrl(
        cntrlMtx,
        conDef,
        occ,
        sunAltitude,
        ill,
        dgp,
        tempInAir,
        irrSurfExtSouth,
        tempOut24,
        incAng,
        profAng,
        conIndexPrev,
    )
    for k in range(len(conIndex)):
        if conIndex[k] == 90:
            conIndex = optState(
                conIndexPrev,
                numConWin, 
                illVmx,
                illVertVmx,
                tmxv,
                dmx, 
                skyVis, 
                numWin,
                spHorIllu,
                spVerIllu,
                numSensor,
                numSensorVert,
                numWinLowerPartition,
                numWinUpperPartition,
                cntrlOpt
                )
            break
    return conIndex


def blindAngle(minAng, maxAng, stepAng, cutOffAng, slatHor):
    # Slat Angle
    # Angles in degrees
    if slatHor:
        for angle in range(int(minAng), int(maxAng) + int(stepAng), int(stepAng)):
            if cutOffAng <= minAng:
                blindAng = minAng
                break
            elif cutOffAng >= maxAng:
                blindAng = maxAng
                break
            elif cutOffAng < angle:
                blindAng = angle
                break
    else:
        for angle in range(int(maxAng), int(minAng) - int(stepAng), -int(stepAng)):
            if cutOffAng <= minAng:
                blindAng = minAng
                break
            elif cutOffAng >= maxAng:
                blindAng = maxAng
                break
            else:
                if cutOffAng > 90:
                    if cutOffAng >= angle:
                        blindAng = angle
                        break
                else:
                    if cutOffAng > angle:
                        blindAng = angle + stepAng
                        break
                    elif cutOffAng == angle:
                        blindAng = angle
                        break
    return blindAng


def cutOff(
    vertical,
    sunAltitude,
    sunAzimuth,
    surfAzimuth,
    irrdirn,
    slatHor,
    slatWidth,
    slatDist,
    minAng,
    maxAng,
    stepAng,
    conOpen,
    conMinAng,
    conMaxAng,
):
    # cut-off
    if vertical:
        incAng = calcs.incAngle(
            radians(sunAltitude), radians(sunAzimuth), radians(surfAzimuth), vertical
        )  # radians
        profAng = calcs.profAngle(
            radians(sunAltitude), radians(sunAzimuth), radians(surfAzimuth), slatHor
        )  # radians
        if slatHor:
            # when the sun does not hit the window, conIndex is set to "open"
            if incAng >= pi / 2.0 or incAng <= -pi / 2.0 or irrdirn < 25.0:
                conIndex = conOpen
            else:
                # angles in degrees
                cutOffAng = cutOffAngle(
                    profAng, incAng, slatWidth, slatDist, vertical, slatHor
                )
                blindAng = blindAngle(
                    minAng, maxAng, stepAng, degrees(cutOffAng), slatHor
                )
                conIndex = conMinAng + int((blindAng - minAng) / stepAng)
                if conIndex > conMaxAng:
                    conIndex = conMaxAng
        else:

            print(degrees(profAng))
            # when the sun does not hit the window, conIndex is set to "open"
            if irrdirn < 25.0:
                conIndex = conMaxAng
            # West facade_ extreme azimuth angles
            elif profAng < -pi:
                print("west")
                if profAng + 2.0 * pi <= -arccos(slatWidth / (slatWidth + slatDist)):
                    print("w1")
                    cutOffAng = pi - cutOffAngle(
                        -profAng - 2.0 * pi,
                        incAng,
                        slatWidth,
                        slatDist,
                        vertical,
                        slatHor,
                    )
                    blindAng = blindAngle(
                        minAng, maxAng, stepAng, degrees(cutOffAng), slatHor
                    )
                    conIndex = conMinAng + int((blindAng - minAng) / stepAng)
                    print(degrees(cutOffAng))
                elif profAng + 2.0 * pi >= -arccos(
                    slatWidth / (slatWidth + slatDist)
                ) and profAng + 2.0 * pi <= arccos(slatWidth / (slatWidth + slatDist)):
                    print("w2")
                    conIndex = (conMaxAng - conMinAng) / 2.0
                elif profAng + 2.0 * pi >= arccos(slatWidth / (slatWidth + slatDist)):
                    print("w3")
                    cutOffAng = cutOffAngle(
                        profAng + 2.0 * pi,
                        incAng,
                        slatWidth,
                        slatDist,
                        vertical,
                        slatHor,
                    )
                    blindAng = blindAngle(
                        minAng, maxAng, stepAng, degrees(cutOffAng), slatHor
                    )
                    conIndex = conMinAng + int((blindAng - minAng) / stepAng)
                    print(degrees(cutOffAng))
                print(degrees(cutOffAng))
            # East facade_extreme azimuth angles
            elif profAng > pi:
                print("east")
                if profAng - 2.0 * pi < -arccos(slatWidth / (slatWidth + slatDist)):
                    print("e1")
                    cutOffAng = pi - cutOffAngle(
                        -profAng + 2.0 * pi,
                        incAng,
                        slatWidth,
                        slatDist,
                        vertical,
                        slatHor,
                    )
                    blindAng = blindAngle(
                        minAng, maxAng, stepAng, degrees(cutOffAng), slatHor
                    )
                    conIndex = conMinAng + int((blindAng - minAng) / stepAng)
                    print(degrees(cutOffAng))
                elif profAng - 2.0 * pi >= -arccos(
                    slatWidth / (slatWidth + slatDist)
                ) and profAng - 2.0 * pi <= arccos(slatWidth / (slatWidth + slatDist)):
                    print("e2")
                    conIndex = (conMaxAng - conMinAng) / 2.0
                elif profAng - 2.0 * pi >= arccos(slatWidth / (slatWidth + slatDist)):
                    print("e3")
                    cutOffAng = cutOffAngle(
                        profAng - 2.0 * pi,
                        incAng,
                        slatWidth,
                        slatDist,
                        vertical,
                        slatHor,
                    )
                    blindAng = blindAngle(
                        minAng, maxAng, stepAng, degrees(cutOffAng), slatHor
                    )
                    conIndex = conMinAng + int((blindAng - minAng) / stepAng)
                    print(degrees(cutOffAng))
            # South facade_general case
            else:
                print("south")
                if profAng < -arccos(slatWidth / (slatWidth + slatDist)):
                    print("s1")
                    cutOffAng = pi - cutOffAngle(
                        -profAng, incAng, slatWidth, slatDist, vertical, slatHor
                    )
                    blindAng = blindAngle(
                        minAng, maxAng, stepAng, degrees(cutOffAng), slatHor
                    )
                    conIndex = conMinAng + int((blindAng - minAng) / stepAng)
                    print(degrees(cutOffAng))
                elif profAng >= -arccos(
                    slatWidth / (slatWidth + slatDist)
                ) and profAng <= arccos(slatWidth / (slatWidth + slatDist)):
                    print("s2")
                    conIndex = (conMaxAng - conMinAng) / 2.0
                elif profAng > arccos(slatWidth / (slatWidth + slatDist)):
                    print("s3")
                    cutOffAng = cutOffAngle(
                        profAng, incAng, slatWidth, slatDist, vertical, slatHor
                    )
                    blindAng = blindAngle(
                        minAng, maxAng, stepAng, degrees(cutOffAng), slatHor
                    )
                    conIndex = conMinAng + int((blindAng - minAng) / stepAng)
                    print(degrees(cutOffAng))
    else:
        raise TypeError("Cut-off algorithm is only implemented for vertical windows")
    # -------------------------------------------------------
    print("conIndex", conIndex)
    return conIndex


def refeed(
    cntrlMtx,
    conDef,
    occ,
    sunAltitude,
    ill,
    dgp,
    tempInAir,
    irrSurfExtSouth,
    tempOut24,
    incAng,
    profAng,
    conIndexPrev,
    illVmx,
    tmxv,
    dmx,
    skyVis,
    numWin,
    illVert,
    illVertVmx,
    numSensor,
    numSensorVert,
):
    # -------------------------------------------------------
    conIndex = mtxCntrl(
        cntrlMtx,
        conDef,
        occ,
        sunAltitude,
        mean(ill),
        max(dgp),
        tempInAir,
        irrSurfExtSouth,
        tempOut24,
        incAng,
        profAng,
        conIndexPrev,
    )
    # ----------------------------------------------------
    # shade closes
    if conIndex[0] > conIndexPrev[0]:
        print(
            "shade closes",
            conIndex,
            conIndexPrev,
            "ill",
            mean(ill),
            "dgp",
            max(dgp),
            "temp",
            tempInAir,
        )
        if sunAltitude > 0.0:
            illPred = zeros(numSensor)
            illPred = solar.illum(illPred, illVmx, tmxv, conIndex, dmx, skyVis, numWin)
            conIndexPred = mtxCntrl(
                cntrlMtx,
                conDef,
                occ,
                sunAltitude,
                mean(illPred),
                max(dgp),
                tempInAir - 1,
                irrSurfExtSouth,
                tempOut24,
                incAng,
                profAng,
                conIndex,
            )
            # after actuation shade would open again
            print("ill", mean(illPred), "dgp", max(dgp), "temp", tempInAir - 1)
            if conIndexPred[0] < conIndex[0]:
                print(
                    "after actuation shade would open again",
                    "pred",
                    conIndexPred,
                    conIndex,
                )
                conIndex = conIndexPrev
    # shade opens
    elif conIndex[0] < conIndexPrev[0]:
        print(
            "shade opens",
            conIndex,
            conIndexPrev,
            "ill",
            mean(ill),
            "dgp",
            max(dgp),
            "temp",
            tempInAir,
        )
        if sunAltitude > 0.0:
            illVertPred = zeros(numSensorVert)  # matrix(zeros([numSensorVert,1]))
            dgpPred = zeros(numSensorVert)
            illVertPred, dgpPred = solar.glareSimpl(
                illVertPred,
                dgpPred,
                illVertVmx,
                tmxv,
                conIndex,
                dmx,
                skyVis,
                numWin,
                numSensorVert,
            )
            conIndexPred = mtxCntrl(
                cntrlMtx,
                conDef,
                occ,
                sunAltitude,
                mean(ill),
                max(dgpPred),
                tempInAir + 1,
                irrSurfExtSouth,
                tempOut24,
                incAng,
                profAng,
                conIndex,
            )
            # after actuation shade would close again
            print("ill", mean(ill), "dgp", max(dgpPred), "temp", tempInAir + 1)
            if conIndexPred[0] > conIndex[0]:
                print(
                    "after actuation shade would close again",
                    "pred",
                    conIndexPred,
                    conIndex,
                )
                conIndex = conIndexPrev
    # ~ print('conIndex',conIndex)
    # ~ print('------------------------------')
    # -------------------------------------------------------
    return conIndex


def simpleCutOff(
    sunAltitude,
    sunAzimuth,
    surfAzimuth,
    irrdirn,
    slatWidth,
    slatDist,
    minAng,
    maxAng,
    stepAng,
    conMinAng,
    conMaxAng,
):
    # cut-off
    incAng = calcs.incAngle(
        radians(sunAltitude), radians(sunAzimuth), radians(surfAzimuth), 1
    )
    profAng = calcs.profAngle(
        radians(sunAltitude), radians(sunAzimuth), radians(surfAzimuth), 1
    )
    # ~ print ('incAng',degrees(incAng))
    # angles in degrees
    cutOffAng = cutOffAngle(profAng, incAng, slatWidth, slatDist, 1, 1)
    blindAng = blindAngle(minAng, maxAng, stepAng, degrees(cutOffAng), 1)
    conIndex = conMinAng + int((blindAng - minAng) / stepAng)
    if conIndex > conMaxAng:
        conIndex = conMaxAng
    print("cutOffAng", cutOffAng, "blindAng", blindAng, "conIndex", conIndex)
    return conIndex


def refeedCutOff(
    cntrlMtx,
    conDef,
    occ,
    sunAltitude,
    ill,
    dgp,
    tempInAir,
    irrSurfExtSouth,
    tempOut24,
    incAng,
    profAng,
    conIndexPrev,
    illVmx,
    tmxv,
    dmx,
    skyVis,
    numWin,
    illVert,
    illVertVmx,
    numSensor,
    numSensorVert,
    sunAzimuth,
    surfAzimuth,
    irrdirn,
    slatWidth,
    slatDist,
    minAng,
    maxAng,
    stepAng,
    conMinAng,
    conMaxAng,
):
    # -------------------------------------------------------
    conIndex = mtxCntrl(
        cntrlMtx,
        conDef,
        occ,
        sunAltitude,
        mean(ill),
        max(dgp),
        tempInAir,
        irrSurfExtSouth,
        tempOut24,
        incAng,
        profAng,
        conIndexPrev,
    )
    for k in range(len(conIndex)):
        if conIndex[k] == 99:
            conIndex[k] = simpleCutOff(
                sunAltitude,
                sunAzimuth,
                surfAzimuth,
                irrdirn,
                slatWidth,
                slatDist,
                minAng,
                maxAng,
                stepAng,
                conMinAng,
                conMaxAng,
            )
    # ----------------------------------------------------
    # shade closes
    if conIndex[0] > conIndexPrev[0]:
        print(
            "shade closes",
            conIndex,
            conIndexPrev,
            "ill",
            mean(ill),
            "dgp",
            max(dgp),
            "temp",
            tempInAir,
        )
        if sunAltitude > 0.0:
            illPred = zeros(numSensor)
            illPred = solar.illum(illPred, illVmx, tmxv, conIndex, dmx, skyVis, numWin)
            conIndexPred = mtxCntrl(
                cntrlMtx,
                conDef,
                occ,
                sunAltitude,
                mean(illPred),
                max(dgp),
                tempInAir - 1,
                irrSurfExtSouth,
                tempOut24,
                incAng,
                profAng,
                conIndex,
            )
            for k in range(len(conIndex)):
                if conIndexPred[k] == 99:
                    conIndexPred[k] = simpleCutOff(
                        sunAltitude,
                        sunAzimuth,
                        surfAzimuth,
                        irrdirn,
                        slatWidth,
                        slatDist,
                        minAng,
                        maxAng,
                        stepAng,
                        conMinAng,
                        conMaxAng,
                    )
            # after actuation shade would open again
            print("ill", mean(illPred), "dgp", max(dgp), "temp", tempInAir - 1)
            if conIndexPred[0] < conIndex[0]:
                print(
                    "after actuation shade would open again",
                    "pred",
                    conIndexPred,
                    conIndex,
                )
                conIndex = conIndexPrev
    # shade opens
    elif conIndex[0] < conIndexPrev[0]:
        print(
            "shade opens",
            conIndex,
            conIndexPrev,
            "ill",
            mean(ill),
            "dgp",
            max(dgp),
            "temp",
            tempInAir,
        )
        if sunAltitude > 0.0:
            illVertPred = zeros(numSensorVert)  # matrix(zeros([numSensorVert,1]))
            dgpPred = zeros(numSensorVert)
            illVertPred, dgpPred = solar.glareSimpl(
                illVertPred,
                dgpPred,
                illVertVmx,
                tmxv,
                conIndex,
                dmx,
                skyVis,
                numWin,
                numSensorVert,
            )
            conIndexPred = mtxCntrl(
                cntrlMtx,
                conDef,
                occ,
                sunAltitude,
                mean(ill),
                max(dgpPred),
                tempInAir + 1,
                irrSurfExtSouth,
                tempOut24,
                incAng,
                profAng,
                conIndex,
            )
            for k in range(len(conIndex)):
                if conIndexPred[k] == 99:
                    conIndexPred[k] = simpleCutOff(
                        sunAltitude,
                        sunAzimuth,
                        surfAzimuth,
                        irrdirn,
                        slatWidth,
                        slatDist,
                        minAng,
                        maxAng,
                        stepAng,
                        conMinAng,
                        conMaxAng,
                    )
            # after actuation shade would close again
            print("ill", mean(ill), "dgp", max(dgpPred), "temp", tempInAir + 1)
            if conIndexPred[0] > conIndex[0]:
                print(
                    "after actuation shade would close again",
                    "pred",
                    conIndexPred,
                    conIndex,
                )
                conIndex = conIndexPrev
    print("conIndex", conIndex)
    print("------------------------------")
    # -------------------------------------------------------
    return conIndex


# ~ def refeedEC(cntrlMtx,conDef,occ,sunAltitude,ill,dgp,tempInAir,irrSurfExtSouth,\
# ~ tempOut24,incAng,profAng,conIndexPrev,\
# ~ illVmx,tmxv,dmx,skyVis,numWin,\
# ~ illVert,illVertVmx,numSensor,numSensorVert):
# ~ #-------------------------------------------------------
# ~ conIndexAux = mtxCntrl(cntrlMtx,conDef,occ,sunAltitude,mean(ill),max(dgp),tempInAir,irrSurfExtSouth,tempOut24,incAng,profAng)
# ~ conIndex = zeros(numWin,int)
# ~ if conIndexAux == 80:
# ~ conIndex[0] = 2
# ~ conIndex[1] = 0
# ~ conIndex[2] = 0
# ~ conIndex[3] = 2
# ~ # ----------------------------------------------------
# ~ elif conIndexAux == 81:
# ~ conIndex[0] = 2
# ~ conIndex[1] = 1
# ~ conIndex[2] = 1
# ~ conIndex[3] = 2
# ~ else:
# ~ for i in range(numWin):
# ~ conIndex[i] = conIndexAux
# ~ # ----------------------------------------------------
# ~ # shade closes
# ~ if conIndex[0] > conIndexPrev[0]:
# ~ print('shade closes',conIndex,conIndexPrev,'ill',mean(ill),'dgp',max(dgp),'temp',tempInAir)
# ~ if sunAltitude > 0.:
# ~ illPred = zeros(numSensor)
# ~ illPred = solar.illum(illPred,illVmx,tmxv,conIndex,dmx,skyVis,numWin)
# ~ conIndexPred = mtxCntrl(cntrlMtx,conDef,occ,sunAltitude,mean(illPred),max(dgp),tempInAir-1,irrSurfExtSouth,tempOut24,incAng,profAng)
# ~ if conIndexPred == 80 or conIndexPred == 81: conIndexPred = 2
# ~ # after actuation shade would open again
# ~ print('ill',mean(illPred),'dgp',max(dgp),'temp',tempInAir-1)
# ~ if conIndexPred < conIndex[0]:
# ~ print('after actuation shade would open again','pred',conIndexPred,conIndex)
# ~ conIndex = conIndexPrev
# ~ # shade opens
# ~ elif conIndex[0] < conIndexPrev[0]:
# ~ print('shade opens',conIndex,conIndexPrev,'ill',mean(ill),'dgp',max(dgp),'temp',tempInAir)
# ~ if sunAltitude > 0.:
# ~ illVertPred = zeros(numSensorVert) #matrix(zeros([numSensorVert,1]))
# ~ dgpPred = zeros(numSensorVert)
# ~ illVertPred,dgpPred = solar.glareSimpl(illVertPred,dgpPred,illVertVmx,tmxv,conIndex,dmx,skyVis,numWin,numSensorVert)
# ~ conIndexPred = mtxCntrl(cntrlMtx,conDef,occ,sunAltitude,mean(ill),max(dgpPred),tempInAir+1,irrSurfExtSouth,tempOut24,incAng,profAng)
# ~ if conIndexPred == 80 or conIndexPred == 81: conIndexPred = 2
# ~ # after actuation shade would close again
# ~ print('ill',mean(ill),'dgp',max(dgpPred),'temp',tempInAir+1)
# ~ if conIndexPred > conIndex[0]:
# ~ print('after actuation shade would close again','pred',conIndexPred,conIndex)
# ~ conIndex = conIndexPrev
# ~ print('conIndex',conIndex)
# ~ print('------------------------------')
# ~ #-------------------------------------------------------
# ~ return conIndex
