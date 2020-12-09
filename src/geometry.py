#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function

import glob
import os
import fileinput
import itertools
import genMtx
import calcs
import param as p
from numpy import math, cos, sin, arange, tan, pi, atleast_2d, genfromtxt, zeros


def shoeBox(
    surf_path,
    rotAng,
    length,
    width,
    height,
    thickSouth,
    thickEast,
    thickNorth,
    thickWest,
    thickCeiling,
    thickFloor,
    albWall,
    albCeiling,
    albFloor,
    bcSouth,
    bcEast,
    bcNorth,
    bcWest,
    bcCeiling,
    bcFloor,
    conSouth,
    conEast,
    conNorth,
    conWest,
    conCeiling,
    conFloor,
):
    # ----------------------------------------------------
    with open("%s" % surf_path, "w") as fw:
        # ----------------------------------------------------
        fw.write(
            "length (m), height/width(for horizontal surfaces) (m), thickness (m), tx (m), ty (m), tz (m), rx (deg), ry (deg), rz (deg), ExtBoundaryCond {0-interior,1-exterior}, svf, exterior albedo, interior albedo, construction\n"
        )
        fw.write(
            "%1.2f, %1.2f, %1.2f, 0.0, %1.2f, 0.0, 90.0, 0.0, %1.2f, %i, 0.0, %1.2f, %1.2f, %i\n"
            % (
                length,
                width,
                thickFloor,
                width,
                rotAng,
                bcFloor,
                albFloor,
                albFloor,
                conFloor,
            )
        )
        fw.write(
            "%1.2f, %1.2f, %1.2f, 0.0, 0.0, 0.0, 0.0, 0.0, %1.2f, %i, 0.5, %1.2f, %1.2f, %i\n"
            % (length, height, thickSouth, rotAng, bcSouth, albWall, albWall, conSouth)
        )
        fw.write(
            "%1.2f, %1.2f, %1.2f, 0.0, %1.2f, 0.0, 0.0, 0.0, %1.2f, %i, 0.5, %1.2f, %1.2f, %i\n"
            % (
                width,
                height,
                thickEast,
                -length,
                rotAng + 90.0,
                bcEast,
                albWall,
                albWall,
                conEast,
            )
        )
        fw.write(
            "%1.2f, %1.2f, %1.2f, %1.2f, %1.2f, 0.0, 0.0, 0.0, %1.2f, %i, 0.5, %1.2f, %1.2f, %i\n"
            % (
                length,
                height,
                thickNorth,
                -length,
                -width,
                rotAng + 180,
                bcNorth,
                albWall,
                albWall,
                conNorth,
            )
        )
        fw.write(
            "%1.2f, %1.2f, %1.2f, %1.2f, 0.0, 0.0, 0.0, 0.0, %1.2f, %i, 0.5, %1.2f, %1.2f, %i\n"
            % (
                width,
                height,
                thickWest,
                -width,
                rotAng - 90.0,
                bcWest,
                albWall,
                albWall,
                conWest,
            )
        )
        fw.write(
            "%1.2f, %1.2f, %1.2f, 0.0, 0.0, %1.2f, -90.0, 0.0, %1.2f, %i, 1.0, %1.2f, %1.2f, %i\n"
            % (
                length,
                width,
                thickCeiling,
                height,
                rotAng,
                bcCeiling,
                albCeiling,
                albCeiling,
                conCeiling,
            )
        )


def writeRad(fileName):
    with open(fileName, "r") as fw:
        return fw.readlines()


def divide_polygons(element, polygons, y, polygonType):

    result = []

    # when polygon is generating for surface
    if polygonType == "surface":
        a1 = element.xOffset
        a2 = element.xOffset + element.length
        b1 = element.zOffset
        b2 = element.zOffset + element.height
    # when polygon is generating for frame
    else:
        a1 = polygonType.xOffset + element.xOffset
        a2 = polygonType.xOffset + element.xOffset + element.length
        b1 = polygonType.zOffset + element.zOffset
        b2 = polygonType.zOffset + element.zOffset + element.height

    for i1 in range(len(polygons)):
        x1 = polygons[i1][0][0]
        x2 = polygons[i1][2][0]
        z1 = polygons[i1][0][2]
        z2 = polygons[i1][2][2]
        # x1 = polygons[i1][2][0]
        # x2 = polygons[i1][0][0]
        # z1 = polygons[i1][2][2]
        # z2 = polygons[i1][0][2]

        if a2 <= x1 or a1 >= x2:
            result.append(polygons[i1])
            continue
        if b2 <= z1 or b1 >= z2:
            result.append(polygons[i1])
            continue

        if a1 <= x1:
            if a2 >= x2:
                if b1 <= z1:
                    if b2 >= z2:
                        raise ValueError("error in function divide_polygons")
                    else:
                        result.append(
                            [[x1, y, b2], [x1, y, z2], [x2, y, z2], [x2, y, b2]]
                        )
                else:
                    if b2 >= z2:
                        result.append(
                            [[x1, y, z1], [x1, y, b1], [x2, y, b1], [x2, y, z1]]
                        )
                    else:
                        result.append(
                            [[x1, y, z1], [x1, y, b1], [x2, y, b1], [x2, y, z1]]
                        )
                        result.append(
                            [[x1, y, b2], [x1, y, z2], [x2, y, z2], [x2, y, b2]]
                        )
            else:
                if b1 <= z1:
                    if b2 >= z2:
                        result.append(
                            [[a2, y, z1], [a2, y, z2], [x2, y, z2], [x2, y, z1]]
                        )
                    else:
                        result.append(
                            [[x1, y, b2], [x1, y, z2], [a2, y, z2], [a2, y, b2]]
                        )
                        result.append(
                            [[a2, y, z1], [a2, y, z2], [x2, y, z2], [x2, y, z1]]
                        )
                else:
                    if b2 >= z2:
                        result.append(
                            [[x1, y, z1], [x1, y, b1], [a2, y, b1], [a2, y, z1]]
                        )
                        result.append(
                            [[a2, y, z1], [a2, y, z2], [x2, y, z2], [x2, y, z1]]
                        )
                    else:
                        result.append(
                            [[x1, y, b2], [x1, y, z2], [a2, y, z2], [a2, y, b2]]
                        )
                        result.append(
                            [[x1, y, z1], [x1, y, b1], [a2, y, b1], [a2, y, z1]]
                        )
                        result.append(
                            [[a2, y, z1], [a2, y, z2], [x2, y, z2], [x2, y, z1]]
                        )
        else:
            if a2 >= x2:
                if b1 <= z1:
                    if b2 >= z2:
                        result.append(
                            [[x1, y, z1], [x1, y, z2], [a1, y, z2], [a1, y, z1]]
                        )
                    else:
                        result.append(
                            [[x1, y, z1], [x1, y, z2], [a1, y, z2], [a1, y, z1]]
                        )
                        result.append(
                            [[a1, y, b2], [a1, y, z2], [x2, y, z2], [x2, y, b2]]
                        )
                else:
                    if b2 >= z2:
                        result.append(
                            [[x1, y, z1], [x1, y, z2], [a1, y, z2], [a1, y, z1]]
                        )
                        result.append(
                            [[a1, y, z1], [a1, y, b1], [x2, y, b1], [x2, y, z1]]
                        )
                    else:
                        result.append(
                            [[x1, y, z1], [x1, y, z2], [a1, y, z2], [a1, y, z1]]
                        )
                        result.append(
                            [[a1, y, b2], [a1, y, z2], [x2, y, z2], [x2, y, b2]]
                        )
                        result.append(
                            [[a1, y, z1], [a1, y, b1], [x2, y, b1], [x2, y, z1]]
                        )
            else:
                if b1 <= z1:
                    if b2 >= z2:
                        result.append(
                            [[x1, y, z1], [x1, y, z2], [a1, y, z2], [a1, y, z1]]
                        )
                        result.append(
                            [[a2, y, z1], [a2, y, z2], [x2, y, z2], [x2, y, z2]]
                        )
                    else:
                        result.append(
                            [[x1, y, z1], [x1, y, z2], [a1, y, z2], [a1, y, z1]]
                        )
                        result.append(
                            [[a1, y, b2], [a1, y, z2], [a2, y, z2], [a2, y, b2]]
                        )
                        result.append(
                            [[a2, y, z1], [a2, y, z2], [x2, y, z2], [x2, y, z1]]
                        )
                else:
                    if b2 >= z2:
                        result.append(
                            [[x1, y, z1], [x1, y, z2], [a1, y, z2], [a1, y, z1]]
                        )
                        result.append(
                            [[a1, y, z1], [a1, y, b1], [a2, y, b1], [a2, y, z1]]
                        )
                        result.append(
                            [[a2, y, z1], [a2, y, z2], [x2, y, z2], [x2, y, z1]]
                        )
                    else:
                        result.append(
                            [[x1, y, z1], [x1, y, z2], [a1, y, z2], [a1, y, z1]]
                        )
                        result.append(
                            [[a1, y, b2], [a1, y, z2], [a2, y, z2], [a2, y, b2]]
                        )
                        result.append(
                            [[a1, y, z1], [a1, y, b1], [a2, y, b1], [a2, y, z1]]
                        )
                        result.append(
                            [[a2, y, z1], [a2, y, z2], [x2, y, z2], [x2, y, z1]]
                        )

    return result


def layerFrontBack(
    name, outElement, fwo, fwi, mat, iniElement, numElement, inElement, surf
):
    # y_off : 0 or self.thick/ thickL-->0 / "outElement" is a subject from class outElement
    for position in (0, 1):
        if position == 0:
            y_off = 0
            name = "Outer"
            fw = fwo
        else:
            y_off = outElement.thick
            name = "Inner"
            fw = fwi

        # when polygon is generating for surface
        if outElement == surf:
            p1 = [0, y_off, 0]
            p2 = [0, y_off, outElement.height]
            p3 = [outElement.length, y_off, outElement.height]
            p4 = [outElement.length, y_off, 0]
            # p1 = [outElement.length,y_off,outElement.height]
            # p2 = [0,y_off,outElement.height]
            # p3 = [0,y_off,0]
            # p4 = [outElement.length,y_off,0]

            polygons = [[p1, p2, p3, p4]]
            for i1 in range(iniElement, iniElement + numElement):
                polygons = divide_polygons(inElement[i1], polygons, y_off, "surface")
        # when polygon is generating for frame
        else:
            p1 = [outElement.xOffset, y_off + outElement.reveal, outElement.zOffset]
            p2 = [
                outElement.xOffset,
                y_off + outElement.reveal,
                outElement.zOffset + outElement.height,
            ]
            p3 = [
                outElement.xOffset + outElement.length,
                y_off + outElement.reveal,
                outElement.zOffset + outElement.height,
            ]
            p4 = [
                outElement.xOffset + outElement.length,
                y_off + outElement.reveal,
                outElement.zOffset,
            ]
            # p1 = [outElement.xOffset+outElement.length,y_off+outElement.reveal,outElement.zOffset+outElement.height]
            # p2 = [outElement.xOffset,y_off+outElement.reveal,outElement.zOffset+outElement.height]
            # p3 = [outElement.xOffset,y_off+outElement.reveal,outElement.zOffset]
            # p4 = [outElement.xOffset+outElement.length,y_off+outElement.reveal,outElement.zOffset]

            polygons = [[p1, p2, p3, p4]]
            for i1 in range(iniElement, iniElement + numElement):
                polygons = divide_polygons(
                    inElement[i1], polygons, y_off + outElement.reveal, outElement
                )

        for i1 in polygons:
            x1, y1, z1 = surf.orgRot(i1[0][0], i1[0][1], i1[0][2])
            x2, y2, z2 = surf.orgRot(i1[1][0], i1[1][1], i1[1][2])
            x3, y3, z3 = surf.orgRot(i1[2][0], i1[2][1], i1[2][2])
            x4, y4, z4 = surf.orgRot(i1[3][0], i1[3][1], i1[3][2])

            fw.write("%s polygon %s_%s\n" % (mat, name, position))
            fw.write("0\n")
            fw.write("0\n")
            fw.write("12\n")
            fw.write("%1.3f %1.3f %1.3f\n" % (x1, y1, z1))
            fw.write("%1.3f %1.3f %1.3f\n" % (x2, y2, z2))
            fw.write("%1.3f %1.3f %1.3f\n" % (x3, y3, z3))
            fw.write("%1.3f %1.3f %1.3f\n" % (x4, y4, z4))
            fw.write("\n")


def layerSide(name, outElement, fw, mat, y, iniElement, numElement, inElement, surf):
    # when polygon is generating for surface
    if outElement == surf:
        tb = [0, outElement.length]
        rl = [0, outElement.height]
        for i in range(numElement):
            tb = tb + [
                inElement[iniElement + i].xOffset,
                inElement[iniElement + i].xOffset + inElement[iniElement + i].length,
            ]
            rl = rl + [
                inElement[iniElement + i].zOffset,
                inElement[iniElement + i].zOffset + inElement[iniElement + i].height,
            ]
        outerY = 0
        innerY = y
    # when polygon is generating for frame
    else:
        tb = [outElement.xOffset, outElement.xOffset + outElement.length]
        rl = [outElement.zOffset, outElement.zOffset + outElement.height]
        for i in range(numElement):
            tb = tb + [
                outElement.xOffset + inElement[iniElement + i].xOffset,
                outElement.xOffset
                + inElement[iniElement + i].xOffset
                + inElement[iniElement + i].length,
            ]
            rl = rl + [
                outElement.zOffset + inElement[iniElement + i].zOffset,
                outElement.zOffset
                + inElement[iniElement + i].zOffset
                + inElement[iniElement + i].height,
            ]
        outerY = outElement.reveal
        innerY = outElement.reveal + y
    for j in range(0, len(tb), 2):
        for k in range(j, j + 2):
            x1, y1, z1 = surf.orgRot(tb[j], outerY, rl[k])
            x2, y2, z2 = surf.orgRot(tb[j], innerY, rl[k])
            x3, y3, z3 = surf.orgRot(tb[j + 1], innerY, rl[k])
            x4, y4, z4 = surf.orgRot(tb[j + 1], outerY, rl[k])
            # x1,y1,z1 = surf.orgRot(tb[j+1],innerY,rl[k])
            # x2,y2,z2 = surf.orgRot(tb[j],innerY,rl[k])
            # x3,y3,z3 = surf.orgRot(tb[j],outerY,rl[k])
            # x4,y4,z4 = surf.orgRot(tb[j+1],outerY,rl[k])

            fw.write("%s polygon %sSide_horizontal_%i\n" % (mat, name, k))
            fw.write("0\n")
            fw.write("0\n")
            fw.write("12\n")
            fw.write("%1.3f %1.3f %1.3f\n" % (x1, y1, z1))
            fw.write("%1.3f %1.3f %1.3f\n" % (x2, y2, z2))
            fw.write("%1.3f %1.3f %1.3f\n" % (x3, y3, z3))
            fw.write("%1.3f %1.3f %1.3f\n" % (x4, y4, z4))
            fw.write("\n")
    for j in range(0, len(tb), 2):
        for k in range(j, j + 2):
            x1, y1, z1 = surf.orgRot(tb[k], outerY, rl[j])
            x2, y2, z2 = surf.orgRot(tb[k], innerY, rl[j])
            x3, y3, z3 = surf.orgRot(tb[k], innerY, rl[j + 1])
            x4, y4, z4 = surf.orgRot(tb[k], outerY, rl[j + 1])

            # x1,y1,z1 = surf.orgRot(tb[k],innerY,rl[j+1])
            # x2,y2,z2 = surf.orgRot(tb[k],innerY,rl[j])
            # x3,y3,z3 = surf.orgRot(tb[k],outerY,rl[j])
            # x4,y4,z4 = surf.orgRot(tb[k],outerY,rl[j+1])

            fw.write("%s polygon %sSide_vertical_%i\n" % (mat, name, k))
            fw.write("0\n")
            fw.write("0\n")
            fw.write("12\n")
            fw.write("%1.3f %1.3f %1.3f\n" % (x1, y1, z1))
            fw.write("%1.3f %1.3f %1.3f\n" % (x2, y2, z2))
            fw.write("%1.3f %1.3f %1.3f\n" % (x3, y3, z3))
            fw.write("%1.3f %1.3f %1.3f\n" % (x4, y4, z4))
            fw.write("\n")


def zone(surf, frame, win, workDir):
    genSky("sky.rad", "void", workDir)
    genSky("sky_prev.rad", "skyfunc", workDir)
    genMaterial(surf, frame, workDir)

    numWin = len(win)
    numSurf = len(surf)
    numFrame = len(frame)
    # Building Geometry
    with open("%sroom.rad" % workDir, "w") as fr:
        countFrame = numFrame
        countWin = numWin
        # ---------------------------------------------------------------------------------------------------
        # start Surface Loop
        for i in range(numSurf):
            numFrameInSurf = 0
            # ------------------------------------------------------------------------------------------------
            # start Frame Loop
            for j in range(numFrame):
                if frame[j].surf == i:
                    numFrameInSurf += 1
                    numWinInFrame = 0
                    # --------------------------------------------------------------------------------------------
                    # start Window Loop
                    for k in range(numWin):
                        if win[k].frame == j:
                            numWinInFrame += 1
                            with open("%swindow_%i.rad" % (workDir, k), "w") as fw:
                                # with open('%swindowRadVis_%i.rad' %(workDir,k), 'w') as fwr:
                                win[k].genPolygon(
                                    fw, "win_mat", "window_%i" % k, surf[i], frame[j]
                                )
                                # win[k].genPolygon(fwr,'glazSys_mat','window_%i'%k,surf[i],frame[j])
                    # finish Window Loop
                    # --------------------------------------------------------------------------------------------
                    if numWinInFrame == 0:
                        raise TypeError("frame without window!")
                    else:
                        with open(
                            "%sframeOuter%i.rad" % (workDir, j), "w"
                        ) as ffo, open(
                            "%sframeInner%i.rad" % (workDir, j), "w"
                        ) as ffi, open(
                            "%sframeSide%i.rad" % (workDir, j), "w"
                        ) as ffs:
                            layerFrontBack(
                                "frame",
                                frame[j],
                                ffo,
                                ffi,
                                "frame_mat_%i" % j,
                                numWin - countWin,
                                numWinInFrame,
                                win,
                                surf[i],
                            )
                            layerSide(
                                "frame",
                                frame[j],
                                ffs,
                                "frame_mat_%i" % j,
                                frame[j].thick,
                                numWin - countWin,
                                numWinInFrame,
                                win,
                                surf[i],
                            )
                        fr.write(
                            "!xform %s/%sframeOuter%i.rad\n" % (os.getcwd(), workDir, j)
                        )
                        fr.write(
                            "!xform %s/%sframeInner%i.rad\n" % (os.getcwd(), workDir, j)
                        )
                        fr.write(
                            "!xform %s/%sframeSide%i.rad\n" % (os.getcwd(), workDir, j)
                        )
                    countWin = countWin - numWinInFrame
            # finish Frame Loop
            # ------------------------------------------------------------------------------------------------
            with open("%ssurfaceOuter%i.rad" % (workDir, i), "w") as fso, open(
                "%ssurfaceInner%i.rad" % (workDir, i), "w"
            ) as fsi, open("%ssurfaceSide%i.rad" % (workDir, i), "w") as fss:
                layerFrontBack(
                    "surf",
                    surf[i],
                    fso,
                    fsi,
                    "surf_mat_%i" % i,
                    numFrame - countFrame,
                    numFrameInSurf,
                    frame,
                    surf[i],
                )
                layerSide(
                    "surf",
                    surf[i],
                    fss,
                    "surf_mat_%i" % i,
                    surf[i].thick,
                    numFrame - countFrame,
                    numFrameInSurf,
                    frame,
                    surf[i],
                )
            fr.write("!xform %s/%ssurfaceOuter%i.rad\n" % (os.getcwd(), workDir, i))
            fr.write("!xform %s/%ssurfaceInner%i.rad\n" % (os.getcwd(), workDir, i))
            fr.write("!xform %s/%ssurfaceSide%i.rad\n" % (os.getcwd(), workDir, i))
            countFrame = countFrame - numFrameInSurf
    # finish Surface Loop
    # ---------------------------------------------------------------------------------------------------


def winRadBsdf(numWin, workDir, bsdfVis_path, conIndex):
    for i in range(numWin):
        with open("%swindowBSDF_%i.rad" % (workDir, i), "w") as o:
            o.write("void BSDF win_bsdf\n")
            o.write("6 thick %s ux uy uz func.cal\n" % bsdfVis_path[conIndex[i]])
            o.write("0\n")
            o.write("0\n")
            o.write("\n")
            k = 0
            with open("%swindow_%i.rad" % (workDir, i)) as windowFile:
                for line in windowFile:
                    line = line.replace("win_mat", "win_bsdf")
                    if k > 1:
                        o.write(line)
                    k = k + 1
    # ----------------------------------------------------


def winRad(numWin, workDir, inputDir, conIndex):
    # ----------------------------------------------------
    for i in range(numWin):
        with open("%swindowRad_%i.rad" % (inputDir, i), "w") as o:
            k = 0
            with open("%swindow_%i.rad" % (workDir, i)) as windowFile:
                for line in windowFile:
                    line = line.replace("win_mat", "winRad_mat_%i" % conIndex[i])
                    if k > 1:
                        o.write(line)
                    k = k + 1
    # ----------------------------------------------------


def photocells(
    xnum,
    ynum,
    sensorHeight,
    gridXOffset,
    gridYOffset,
    thickSouth,
    thickEast,
    thickNorth,
    thickWest,
    thickFloor,
    surf,
    workDir,
):
    with open("%silluSensor.pts" % workDir, "w") as fw:
        pi = math.pi / 180
        if xnum - 1 > 0:
            xlen = (surf.length - thickEast - thickWest - gridXOffset * 2) / (xnum - 1)
        else:
            print(
                "ERROR: Grid of sensors must contain more than one sensor in the x-axis"
            )
        if ynum - 1 > 0:
            ylen = (surf.height - thickSouth - thickNorth - gridYOffset * 2) / (
                ynum - 1
            )
        else:
            print(
                "ERROR: Grid of sensors must contain more than one sensor in the y-axis"
            )
        xOffset = thickWest + gridXOffset
        yOffset = thickSouth + gridYOffset
        for j in range(ynum):
            for i in range(xnum):
                fw.write(
                    "%1.3f %1.3f %1.3f 0 0 1\n"
                    % (
                        (xOffset + i * xlen) * cos(pi * surf.rz)
                        - (yOffset + j * ylen) * sin(pi * surf.rz),
                        (xOffset + i * xlen) * sin(pi * surf.rz)
                        + (yOffset + j * ylen) * cos(pi * surf.rz),
                        sensorHeight,
                    )
                )
    # -------------------------------------------------------


def extIlluSensor(workDir, xOffset, yOffset, height):
    with open("%sextIlluSensor.pts" % workDir, "w") as fw:
        fw.write("%1.3f\t%1.3f\t%1.3f\t0\t0\t1\n" % (xOffset, yOffset, height))


def winSensorPts(resolutWin, surf, frame, win, workDir):
    with open("%snumSensorWin.dat" % workDir, "w") as fw:
        for ns in range(len(surf)):
            xFrontDir, yFrontDir, zFrontDir = surf[ns].dirRot(0, 1, 0)
            xBackDir, yBackDir, zBackDir = surf[ns].dirRot(0, -1, 0)
            frameInSurf = []
            for nf in range(len(frame)):
                if frame[nf].surf == ns:
                    frameInSurf += [nf]
                    winInFrame = []
                    for nw in range(len(win)):
                        if win[nw].frame == nf:
                            # WINDOWS
                            with open(
                                "%sintWinFrontSensor_%i.pts" % (workDir, nw), "w"
                            ) as fwf, open(
                                "%sintWinBackSensor_%i.pts" % (workDir, nw), "w"
                            ) as fwb:
                                numSensorWin = 0
                                numSensorWinLength = int(win[nw].length / resolutWin)
                                numSensorWinHeight = int(win[nw].height / resolutWin)
                                locWinLength = arange(
                                    1.0 / numSensorWinLength / 2,
                                    1,
                                    1.0 / numSensorWinLength,
                                )
                                locWinHeight = arange(
                                    1.0 / numSensorWinHeight / 2,
                                    1,
                                    1.0 / numSensorWinHeight,
                                )
                                winInFrame += [nw]
                                for i in range(int(numSensorWinLength)):
                                    for j in range(int(numSensorWinHeight)):
                                        txWin = (
                                            frame[nf].xOffset
                                            + win[nw].xOffset
                                            + locWinLength[i] * win[nw].length
                                        )
                                        tyWin = (
                                            frame[nf].reveal + win[nw].reveal + 0.001
                                        )
                                        tzWin = (
                                            frame[nf].zOffset
                                            + win[nw].zOffset
                                            + locWinHeight[j] * win[nw].height
                                        )
                                        xorgWin, yorgWin, zorgWin = surf[ns].orgRot(
                                            txWin, tyWin, tzWin
                                        )
                                        fwf.write(
                                            "%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n"
                                            % (
                                                xorgWin,
                                                yorgWin,
                                                zorgWin,
                                                xFrontDir,
                                                yFrontDir,
                                                zFrontDir,
                                            )
                                        )
                                        fwb.write(
                                            "%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n"
                                            % (
                                                xorgWin,
                                                yorgWin,
                                                zorgWin,
                                                xBackDir,
                                                yBackDir,
                                                zBackDir,
                                            )
                                        )
                                        numSensorWin += 1
                                fw.write("%i\t" % numSensorWin)


def irradianceSensorPts(
    resolutSurf, resolutFrame, resolutWin, resolutRev, surf, frame, win, workDir
):
    with open("%snumSensorSurf.dat" % workDir, "w") as fs, open(
        "%snumSensorFrame.dat" % workDir, "w"
    ) as ff, open("%snumSensorWin.dat" % workDir, "w") as fw, open(
        "%snumSensorRev.dat" % workDir, "w"
    ) as fr:
        for ns in range(len(surf)):
            numSensorSurf = 0
            xFrontDir, yFrontDir, zFrontDir = surf[ns].dirRot(0, 1, 0)
            xBackDir, yBackDir, zBackDir = surf[ns].dirRot(0, -1, 0)
            frameInSurf = []
            for nf in range(len(frame)):
                if frame[nf].surf == ns:
                    frameInSurf += [nf]
                    winInFrame = []
                    for nw in range(len(win)):
                        if win[nw].frame == nf:
                            # WINDOWS
                            with open(
                                "%sintWinFrontSensor_%i.pts" % (workDir, nw), "w"
                            ) as fwf, open(
                                "%sintWinBackSensor_%i.pts" % (workDir, nw), "w"
                            ) as fwb:
                                numSensorWin = 0
                                numSensorWinLength = int(win[nw].length / resolutWin)
                                numSensorWinHeight = int(win[nw].height / resolutWin)
                                locWinLength = arange(
                                    1.0 / numSensorWinLength / 2,
                                    1,
                                    1.0 / numSensorWinLength,
                                )
                                locWinHeight = arange(
                                    1.0 / numSensorWinHeight / 2,
                                    1,
                                    1.0 / numSensorWinHeight,
                                )
                                winInFrame += [nw]
                                for i in range(int(numSensorWinLength)):
                                    for j in range(int(numSensorWinHeight)):
                                        txWin = (
                                            frame[nf].xOffset
                                            + win[nw].xOffset
                                            + locWinLength[i] * win[nw].length
                                        )
                                        tyWin = (
                                            frame[nf].reveal + win[nw].reveal + 0.001
                                        )
                                        tzWin = (
                                            frame[nf].zOffset
                                            + win[nw].zOffset
                                            + locWinHeight[j] * win[nw].height
                                        )
                                        xorgWin, yorgWin, zorgWin = surf[ns].orgRot(
                                            txWin, tyWin, tzWin
                                        )
                                        fwf.write(
                                            "%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n"
                                            % (
                                                xorgWin,
                                                yorgWin,
                                                zorgWin,
                                                xFrontDir,
                                                yFrontDir,
                                                zFrontDir,
                                            )
                                        )
                                        fwb.write(
                                            "%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n"
                                            % (
                                                xorgWin,
                                                yorgWin,
                                                zorgWin,
                                                xBackDir,
                                                yBackDir,
                                                zBackDir,
                                            )
                                        )
                                        numSensorWin += 1
                                fw.write("%i\t" % numSensorWin)

                    # FRAMES
                    with open(
                        "%sintFrameFrontSensor_%i.pts" % (workDir, nf), "w"
                    ) as fff, open(
                        "%sintFrameBackSensor_%i.pts" % (workDir, nf), "w"
                    ) as ffb, open(
                        "%sextFrameSensor_%i.pts" % (workDir, nf), "w"
                    ) as ffEx:
                        numSensorFrame = 0
                        numSensorFrameLength = int(frame[nf].length / resolutFrame)
                        numSensorFrameHeight = int(frame[nf].height / resolutFrame)
                        locFrameLength = arange(
                            1.0 / numSensorFrameLength / 2,
                            1,
                            1.0 / numSensorFrameLength,
                        )
                        locFrameHeight = arange(
                            1.0 / numSensorFrameHeight / 2,
                            1,
                            1.0 / numSensorFrameHeight,
                        )
                        for i in range(int(numSensorFrameLength)):
                            for j in range(int(numSensorFrameHeight)):
                                txFrame = (
                                    0.001
                                    + frame[nf].xOffset
                                    + locFrameLength[i] * frame[nf].length
                                )
                                tyFrame = frame[nf].reveal + frame[nf].thick + 0.001
                                tyFrameEx = frame[nf].reveal - 0.001
                                tzFrame = (
                                    0.001
                                    + frame[nf].zOffset
                                    + locFrameHeight[j] * frame[nf].height
                                )
                                sensorInWin_flag = False
                                for z in range(len(winInFrame)):
                                    winXoffset = (
                                        frame[nf].xOffset + win[winInFrame[z]].xOffset
                                    )
                                    winZoffset = (
                                        frame[nf].zOffset + win[winInFrame[z]].zOffset
                                    )
                                    if (
                                        txFrame > winXoffset - 0.001
                                        and txFrame
                                        < winXoffset + win[winInFrame[z]].length + 0.001
                                        and tzFrame > winZoffset - 0.001
                                        and tzFrame
                                        < winZoffset + win[winInFrame[z]].height + 0.001
                                    ):
                                        sensorInWin_flag = True
                                if not sensorInWin_flag:
                                    xorgFrame, yorgFrame, zorgFrame = surf[ns].orgRot(
                                        txFrame, tyFrame, tzFrame
                                    )
                                    xorgFrameEx, yorgFrameEx, zorgFrameEx = surf[
                                        ns
                                    ].orgRot(txFrame, tyFrameEx, tzFrame)
                                    fff.write(
                                        "%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n"
                                        % (
                                            xorgFrame,
                                            yorgFrame,
                                            zorgFrame,
                                            xFrontDir,
                                            yFrontDir,
                                            zFrontDir,
                                        )
                                    )
                                    ffb.write(
                                        "%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n"
                                        % (
                                            xorgFrame,
                                            yorgFrame,
                                            zorgFrame,
                                            xBackDir,
                                            yBackDir,
                                            zBackDir,
                                        )
                                    )
                                    ffEx.write(
                                        "%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n"
                                        % (
                                            xorgFrameEx,
                                            yorgFrameEx,
                                            zorgFrameEx,
                                            xBackDir,
                                            yBackDir,
                                            zBackDir,
                                        )
                                    )
                                    numSensorFrame += 1
                        ff.write("%i\t" % numSensorFrame)

                    # REVEALS
                    intRevDepth = surf[ns].thick - frame[nf].reveal - frame[nf].thick
                    if intRevDepth > resolutRev:
                        with open(
                            "%sintRevFrontSensor_%i.pts" % (workDir, nf), "w"
                        ) as frf, open(
                            "%sintRevBackSensor_%i.pts" % (workDir, nf), "w"
                        ) as frb:
                            numSensorRev = 0
                            numSensorRevLength = int(frame[nf].length / resolutRev)
                            numSensorRevHeight = int(frame[nf].height / resolutRev)
                            numSensorRevDepth = max(int(intRevDepth / resolutRev), 1)
                            locRevLength = arange(
                                1.0 / numSensorRevLength / 2,
                                1,
                                1.0 / numSensorRevLength,
                            )
                            locRevHeight = arange(
                                1.0 / numSensorRevHeight / 2,
                                1,
                                1.0 / numSensorRevHeight,
                            )
                            locRevDepth = arange(
                                1.0 / numSensorRevDepth / 2, 1, 1.0 / numSensorRevDepth
                            )
                            for i in range(int(numSensorRevLength)):
                                for j in range(int(numSensorRevDepth)):
                                    tx = (
                                        frame[nf].xOffset
                                        + locRevLength[i] * frame[nf].length
                                    )
                                    ty = surf[ns].thick - locRevDepth[j] * intRevDepth
                                    tzDown = frame[nf].zOffset + 0.001
                                    tzUp = frame[nf].zOffset + frame[nf].height - 0.001
                                    xFrontDirRev, yFrontDirRev, zFrontDirRev = surf[
                                        ns
                                    ].dirRot(0, 0, 1)
                                    xBackDirRev, yBackDirRev, zBackDirRev = surf[
                                        ns
                                    ].dirRot(0, 0, -1)
                                    xorg, yorg, zorg = surf[ns].orgRot(tx, ty, tzDown)
                                    frf.write(
                                        "%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n"
                                        % (
                                            xorg,
                                            yorg,
                                            zorg,
                                            xFrontDirRev,
                                            yFrontDirRev,
                                            zFrontDirRev,
                                        )
                                    )
                                    frb.write(
                                        "%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n"
                                        % (
                                            xorg,
                                            yorg,
                                            zorg,
                                            xBackDirRev,
                                            yBackDirRev,
                                            zBackDirRev,
                                        )
                                    )
                                    numSensorRev += 1
                                    xorg, yorg, zorg = surf[ns].orgRot(
                                        tx, ty, tzUp
                                    )  # dirrection is different from 'tz=tzDown' case
                                    frf.write(
                                        "%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n"
                                        % (
                                            xorg,
                                            yorg,
                                            zorg,
                                            xBackDirRev,
                                            yBackDirRev,
                                            zBackDirRev,
                                        )
                                    )
                                    frb.write(
                                        "%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n"
                                        % (
                                            xorg,
                                            yorg,
                                            zorg,
                                            xFrontDirRev,
                                            yFrontDirRev,
                                            zFrontDirRev,
                                        )
                                    )
                                    numSensorRev += 1
                            for j in range(int(numSensorRevHeight)):
                                for j in range(int(numSensorRevDepth)):
                                    txLeft = frame[nf].xOffset + 0.001
                                    txRight = (
                                        frame[nf].xOffset + frame[nf].length - 0.001
                                    )
                                    ty = surf[ns].thick - locRevDepth[j] * intRevDepth
                                    tz = (
                                        frame[nf].zOffset
                                        + locRevHeight[j] * frame[nf].height
                                    )
                                    xFrontDirRev, yFrontDirRev, zFrontDirRev = surf[
                                        ns
                                    ].dirRot(1, 0, 0)
                                    xBackDirRev, yBackDirRev, zBackDirRev = surf[
                                        ns
                                    ].dirRot(-1, 0, 0)
                                    xorg, yorg, zorg = surf[ns].orgRot(txLeft, ty, tz)
                                    frf.write(
                                        "%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n"
                                        % (
                                            xorg,
                                            yorg,
                                            zorg,
                                            xFrontDirRev,
                                            yFrontDirRev,
                                            zFrontDirRev,
                                        )
                                    )
                                    frb.write(
                                        "%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n"
                                        % (
                                            xorg,
                                            yorg,
                                            zorg,
                                            xBackDirRev,
                                            yBackDirRev,
                                            zBackDirRev,
                                        )
                                    )
                                    numSensorRev += 1
                                    xorg, yorg, zorg = surf[ns].orgRot(
                                        txRight, ty, tz
                                    )  # dirrection is different from 'tx=txLeft' case
                                    frf.write(
                                        "%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n"
                                        % (
                                            xorg,
                                            yorg,
                                            zorg,
                                            xBackDirRev,
                                            yBackDirRev,
                                            zBackDirRev,
                                        )
                                    )
                                    frb.write(
                                        "%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n"
                                        % (
                                            xorg,
                                            yorg,
                                            zorg,
                                            xFrontDirRev,
                                            yFrontDirRev,
                                            zFrontDirRev,
                                        )
                                    )
                                    numSensorRev += 1
                            fr.write("%i\t" % numSensorRev)
                    else:
                        raise TypeError(
                            "Reveal is too narrow for the irradiance resolution given in the config file"
                        )

            # SURFACES (including external sensor)
            with open("%sintSurfFrontSensor_%i.pts" % (workDir, ns), "w") as fsf, open(
                "%sintSurfBackSensor_%i.pts" % (workDir, ns), "w"
            ) as fsb, open("%sextSurfSensor_%i.pts" % (workDir, ns), "w") as fsEx:
                numSensorSurf = 0
                numSensorSurfLength = max(
                    int((surf[ns].length - 2 * p.irrSensorXOffset) / resolutSurf), 1
                )
                numSensorSurfHeight = max(
                    int((surf[ns].height - 2 * p.irrSensorZOffset) / resolutSurf), 1
                )
                locSurfLength = arange(
                    1.0 / numSensorSurfLength / 2, 1, 1.0 / numSensorSurfLength
                )
                locSurfHeight = arange(
                    1.0 / numSensorSurfHeight / 2, 1, 1.0 / numSensorSurfHeight
                )
                for i in range(int(numSensorSurfLength)):
                    for j in range(int(numSensorSurfHeight)):
                        txSurf = p.irrSensorXOffset + locSurfLength[i] * (
                            surf[ns].length - 2 * p.irrSensorXOffset
                        )
                        tySurf = surf[ns].thick + 0.001
                        tySurfEx = -0.001
                        tzSurf = p.irrSensorZOffset + locSurfHeight[j] * (
                            surf[ns].height - 2 * p.irrSensorZOffset
                        )
                        sensorInFrame = 0
                        if len(frameInSurf) > 0:
                            for z in range(len(frameInSurf)):
                                if (
                                    txSurf > (frame[frameInSurf[z]].xOffset)
                                    and txSurf
                                    < (
                                        frame[frameInSurf[z]].xOffset
                                        + frame[frameInSurf[z]].length
                                    )
                                    and tzSurf > frame[frameInSurf[z]].zOffset
                                    and tzSurf
                                    < frame[frameInSurf[z]].height
                                    + frame[frameInSurf[z]].zOffset
                                ):
                                    sensorInFrame = 1
                        if sensorInFrame == 0:
                            xorgSurf, yorgSurf, zorgSurf = surf[ns].orgRot(
                                txSurf, tySurf, tzSurf
                            )
                            xorgSurfEx, yorgSurfEx, zorgSurfEx = surf[ns].orgRot(
                                txSurf, tySurfEx, tzSurf
                            )
                            fsf.write(
                                "%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n"
                                % (
                                    xorgSurf,
                                    yorgSurf,
                                    zorgSurf,
                                    xFrontDir,
                                    yFrontDir,
                                    zFrontDir,
                                )
                            )
                            fsb.write(
                                "%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n"
                                % (
                                    xorgSurf,
                                    yorgSurf,
                                    zorgSurf,
                                    xBackDir,
                                    yBackDir,
                                    zBackDir,
                                )
                            )
                            fsEx.write(
                                "%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n"
                                % (
                                    xorgSurfEx,
                                    yorgSurfEx,
                                    zorgSurfEx,
                                    xBackDir,
                                    yBackDir,
                                    zBackDir,
                                )
                            )
                            numSensorSurf += 1
                if numSensorSurf == 0:
                    numSensorSurf = 1
                    xorgSurf, yorgSurf, zorgSurf = surf[ns].orgRot(
                        0.0, surf[ns].thick + 0.001, 0.0
                    )
                    xorgSurfEx, yorgSurfEx, zorgSurfEx = surf[ns].orgRot(
                        0.0, -0.001, 0.0
                    )
                    fsf.write(
                        "%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n"
                        % (
                            xorgSurf,
                            yorgSurf,
                            zorgSurf,
                            xFrontDir,
                            yFrontDir,
                            zFrontDir,
                        )
                    )
                    fsb.write(
                        "%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n"
                        % (xorgSurf, yorgSurf, zorgSurf, xBackDir, yBackDir, zBackDir)
                    )
                    fsEx.write(
                        "%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n"
                        % (
                            xorgSurfEx,
                            yorgSurfEx,
                            zorgSurfEx,
                            xBackDir,
                            yBackDir,
                            zBackDir,
                        )
                    )
                fs.write("%i\t" % numSensorSurf)


def viewfactors(nCores, workDir):

    surface_files = sorted(glob.glob("%ssurfaceInner*" % workDir))
    frame_files = sorted(glob.glob("%sframeInner*" % workDir))
    window_files = sorted(glob.glob("%swindow_*" % workDir))

    header = "#@rfluxmtx h=u"
    for i1 in [surface_files + frame_files + window_files]:
        for line in fileinput.input(i1, inplace=True):
            if fileinput.isfirstline() and line.split(" ")[0] != "#@rfluxmtx h=u":
                print(header)
            print(line, end="")

    dict1 = {"Surface": surface_files, "Frame": frame_files, "Windows": window_files}

    for i2 in itertools.combinations_with_replacement(dict1, 2):

        with open("%svf%s%s.dat" % (workDir, i2[0], i2[1]), "w") as vff:
            for i1 in itertools.product(dict1[i2[0]], dict1[i2[1]]):
                tmp = genMtx.comp_vf(nCores, i1[0], i1[1], workDir)

                tmp2 = tmp[:-1].split("\t")
                if len(tmp2) != 4:
                    print("warning: wrong format of view factor")
                print(i1[0], i1[1], tmp2[0], file=vff)

        for i2 in itertools.combinations(dict1, 2):

            with open("%svf%s%s.dat" % (workDir, i2[1], i2[0]), "w") as vff:
                for i1 in itertools.product(dict1[i2[1]], dict1[i2[0]]):
                    tmp = genMtx.comp_vf(nCores, i1[0], i1[1], workDir)
                    tmp2 = tmp[:-1].split("\t")
                    if len(tmp2) != 4:
                        print("warning: wrong format of view factor")
                    print(i1[0], i1[1], tmp2[0], file=vff)

        header = "#@rfluxmtx h=u"
        for i1 in [surface_files + frame_files + window_files]:
            for line in fileinput.input(i1, inplace=True):
                if (
                    not fileinput.isfirstline()
                    and line.split(" ")[0] != "#@rfluxmtx h=u"
                ):
                    print(line, end="")


def fileDel(fileName):
    files = glob.glob(fileName)
    for f in files:
        os.remove(f)


def genSky(fileName, modifier, workDir):
    with open("%s%s" % (workDir, fileName), "w") as fw:
        fw.write("#@rfluxmtx h=u u=Y")
        fw.write("\n")
        fw.write("%s glow ground_glow\n" % modifier)
        fw.write("0\n")
        fw.write("0\n")
        fw.write("4  1 1 1 0\n")
        fw.write("\n")
        fw.write("ground_glow source ground\n")
        fw.write("0\n")
        fw.write("0\n")
        fw.write("4  0 0 -1 180\n")
        fw.write("\n")
        fw.write("#@rfluxmtx h=r4 u=Y")
        fw.write("\n")
        fw.write("%s glow sky_glow\n" % modifier)
        fw.write("0\n")
        fw.write("0\n")
        fw.write("4  1 1 1 0\n")
        fw.write("\n")
        fw.write("sky_glow source sky\n")
        fw.write("0\n")
        fw.write("0\n")
        fw.write("4  0 0 1 180\n")


# def genGround(floorsurf,buildHeight,grndAlb,workDir):
#    centerX = floorsurf.length/2
#    centerY = floorsurf.height/2
#    radius = (centerX + centerY)/2 + 10*buildHeight
#    with open('%sground.rad'%workDir,'w') as fw:
#       fw.write('void plastic ground_mat\n')
#       fw.write('0\n')
#       fw.write('0\n')
#       fw.write('5  %1.2f %1.2f %1.2f  0 0\n'%(grndAlb,grndAlb,grndAlb))
#       fw.write('\n')
#       fw.write('ground_mat ring ground\n')
#       fw.write('0\n')
#       fw.write('0\n')
#       fw.write('8  %1.3f %1.3f -0.3  0 0 1  0  %1.3f \n' %(centerX,centerY,radius))
#       fw.write('\n')


def genmask(floorsurf, buildHeight, obstMask, workDir, mask):
    pi = math.pi / 180
    centerX = floorsurf.length / 2
    centerY = floorsurf.height / 2
    radius = (centerX + centerY) / 2 + 10 * buildHeight
    with open("%soutside.rad" % workDir, "a") as fw:
        fw.write("void plastic mask_mat\n")
        fw.write("0\n")
        fw.write("0\n")
        fw.write("5  %1.2f %1.2f %1.2f 0 0\n" % (obstMask[0], obstMask[0], obstMask[0]))
        fw.write("\n")
        for i in range(len(obstMask) - 2):
            x1 = centerX + radius * sin(pi * 30 * (i))
            y1 = centerY + -radius * cos(pi * 30 * (i))
            z1 = buildHeight / 2 + 10 * buildHeight * tan(pi * obstMask[i + 1])
            x2 = centerX + radius * sin(pi * 30 * (i + 1))
            y2 = centerY + -radius * cos(pi * 30 * (i + 1))
            z2 = buildHeight / 2 + 10 * buildHeight * tan(pi * obstMask[i + 2])
            fw.write("mask_mat polygon mask%i\n" % i)
            fw.write("0\n")
            fw.write("0\n")
            fw.write("12\n")
            fw.write("%1.3f %1.3f 0\n" % (x1, y1))
            fw.write("%1.3f %1.3f 0 \n" % (x2, y2))
            fw.write("%1.3f %1.3f %1.3f\n" % (x2, y2, z2))
            fw.write("%1.3f %1.3f %1.3f\n" % (x1, y1, z1))
            fw.write("\n")


def genMaterial(surf, frame, workDir):
    with open("%sroom.mat" % workDir, "w") as fw:
        # -----------------------------------------------------------------------
        for ns in range(len(surf)):
            surf[ns].genMatSurf(fw, ns)
        for nf in range(len(frame)):
            frame[nf].genMatFrame(fw, nf)
        # -----------------------------------------------------------------------
        fw.write("void glow win_mat\n")
        fw.write("0\n")
        fw.write("0\n")
        fw.write("4 1 1 1 0\n")
        # -----------------------------------------------------------------------
        #    fw.write('\n')
        #    fw.write('void BSDF win_bsdf\n')
        #    fw.write('6 thick %s ux uy uz func.cal\n'%bsdf_path)
        #    fw.write('0\n')
        #    fw.write('0\n')
        # ------------------------------------------------------------------------
        # reflFront,reflBack,p,qRefl,trans,qTrans= matGlzSys[4],matGlzSys[5],matGlzSys[8],matGlzSys[12],matGlzSys[3],matGlzSys[11]
        # fw.write('void BRTDfunc glazSys_mat\n')
        # fw.write('10\n')
        # fw.write('if(Rdot,rf(%1.3f,%1.0f,%1.2f),rf(%1.3f,%1.0f,%1.2f)) \n' % (reflFront,p,qRefl,reflBack,p,qRefl) )
        # fw.write('if(Rdot,rf(%1.3f,%1.0f,%1.2f),rf(%1.3f,%1.0f,%1.2f)) \n' % (reflFront,p,qRefl,reflBack,p,qRefl) )
        # fw.write('if(Rdot,rf(%1.3f,%1.0f,%1.2f),rf(%1.3f,%1.0f,%1.2f)) \n' % (reflFront,p,qRefl,reflBack,p,qRefl) )
        # fw.write('tau(%1.3f,%1.0f,%1.2f) \n' % (trans,p,qTrans) )
        # fw.write('tau(%1.3f,%1.0f,%1.2f) \n' % (trans,p,qTrans) )
        # fw.write('tau(%1.3f,%1.0f,%1.2f) \n' % (trans,p,qTrans) )
        # fw.write('        0 0 0\n')
        # fw.write('        trans_refl.cal\n')
        # fw.write('0\n')
        # fw.write('9\n')
        # fw.write('        0 0 0\n')
        # fw.write('        0 0 0\n')
        # fw.write('        0 0 0\n')
        # fw.write('\n')


def genPtsFile(pts_path, newPts_path, rotAng):
    rz = rotAng * pi / 180.0
    illuPtsSouth = atleast_2d(genfromtxt(pts_path, delimiter=" "))
    illuPts = zeros([illuPtsSouth.shape[0], 6])
    with open(newPts_path, "w") as fw:
        for i in range(illuPtsSouth.shape[0]):
            illuPts[i, 0] = illuPtsSouth[i, 0] * cos(rz) - illuPtsSouth[i, 1] * sin(rz)
            illuPts[i, 1] = illuPtsSouth[i, 0] * sin(rz) + illuPtsSouth[i, 1] * cos(rz)
            illuPts[i, 2] = illuPtsSouth[i, 2]
            illuPts[i, 3] = illuPtsSouth[i, 3] * cos(rz) - illuPtsSouth[i, 4] * sin(rz)
            illuPts[i, 4] = illuPtsSouth[i, 3] * sin(rz) + illuPtsSouth[i, 4] * cos(rz)
            illuPts[i, 5] = illuPtsSouth[i, 5]
            fw.write(
                "%1.3f %1.3f %1.3f %1.3f %1.3f %1.3f\n"
                % (
                    illuPts[i, 0],
                    illuPts[i, 1],
                    illuPts[i, 2],
                    illuPts[i, 3],
                    illuPts[i, 4],
                    illuPts[i, 5],
                )
            )
    return illuPts


def genViewFile(pts, workDir):
    with open("%sview.vf" % workDir, "w") as fw:
        for i in range(len(pts)):
            tx, ty, tz = pts[i, 0], pts[i, 1], pts[i, 2]
            rx, ry, rz = pts[i, 3], pts[i, 4], pts[i, 5]
            fw.write(
                "rview -vta -vp %1.3f %1.3f %1.3f -vd %1.3f %1.3f %1.3f -vv 180 -vh 180 -vs 0 -vl 0 -vu 0 0 1\n"
                % (tx, ty, tz, rx, ry, rz)
            )


def shoeBoxVF(input_path, length, width, height, surf_path, frame_path, win_path):
    # ----------------------------------------------------
    # view factors
    viewFact = zeros([6, 6])
    viewFact[0, 5] = calcs.viewFactorOppRectangles(length, width, height)
    viewFact[1, 3] = calcs.viewFactorOppRectangles(length, height, width)
    viewFact[2, 4] = calcs.viewFactorOppRectangles(width, height, length)
    viewFact[5, 0] = viewFact[0, 5]
    viewFact[3, 1] = viewFact[1, 3]
    viewFact[4, 2] = viewFact[2, 4]
    viewFact[1, 0] = calcs.viewFactorAdjRectangles(width, length, height)
    viewFact[0, 1] = calcs.viewFactorAdjRectangles(height, length, width)
    viewFact[1, 5] = viewFact[1, 0]
    viewFact[5, 1] = viewFact[0, 1]
    viewFact[3, 5] = viewFact[1, 0]
    viewFact[5, 3] = viewFact[0, 1]
    viewFact[3, 0] = viewFact[1, 0]
    viewFact[0, 3] = viewFact[0, 1]
    viewFact[1, 2] = calcs.viewFactorAdjRectangles(width, height, length)
    viewFact[2, 1] = calcs.viewFactorAdjRectangles(length, height, width)
    viewFact[1, 4] = viewFact[1, 2]
    viewFact[4, 1] = viewFact[2, 1]
    viewFact[3, 4] = viewFact[1, 2]
    viewFact[4, 3] = viewFact[2, 1]
    viewFact[3, 2] = viewFact[1, 2]
    viewFact[2, 3] = viewFact[2, 1]
    viewFact[2, 0] = calcs.viewFactorAdjRectangles(length, width, height)
    viewFact[0, 2] = calcs.viewFactorAdjRectangles(height, width, length)
    viewFact[2, 5] = viewFact[2, 0]
    viewFact[5, 2] = viewFact[0, 2]
    viewFact[4, 5] = viewFact[2, 0]
    viewFact[5, 4] = viewFact[0, 2]
    viewFact[4, 0] = viewFact[2, 0]
    viewFact[0, 4] = viewFact[0, 2]
    # ----------------------------------------------------
    surfVect = atleast_2d(genfromtxt(surf_path, skip_header=1, delimiter=","))
    winVect = atleast_2d(genfromtxt(win_path, skip_header=1, delimiter=","))
    frameVect = atleast_2d(genfromtxt(frame_path, skip_header=1, delimiter=","))
    # ----------------------------------------------------
    numWin = winVect.shape[0]
    numFrame = frameVect.shape[0]
    numSurf = surfVect.shape[0]
    # ----------------------------------------------------
    areaSurf = zeros(numSurf)
    areaFrame = zeros(numFrame)
    areaWin = zeros(numWin)
    for l in range(numWin):
        areaWin[l] = winVect[l, 0] * winVect[l, 1]
    for k in range(numFrame):
        areaFrame[k] = frameVect[k, 0] * frameVect[k, 1]
        for l in range(numWin):
            if winVect[l, 2] == k:
                areaFrame[k] = areaFrame[k] - winVect[l, 0] * winVect[l, 1]
    for i in range(numSurf):
        areaSurf[i] = surfVect[i, 0] * surfVect[i, 1]
        for k in range(numFrame):
            if frameVect[k, 3] == i:
                areaSurf[i] = areaSurf[i] - frameVect[k, 0] * frameVect[k, 1]
    # ----------------------------------------------------
    with open("%svfSurfaceSurface.dat" % input_path, "w") as fss, open(
        "%svfSurfaceWindows.dat" % input_path, "w"
    ) as fsw, open("%svfSurfaceFrame.dat" % input_path, "w") as fsf, open(
        "%svfWindowsWindows.dat" % input_path, "w"
    ) as fww, open(
        "%svfWindowsSurface.dat" % input_path, "w"
    ) as fws, open(
        "%svfWindowsFrame.dat" % input_path, "w"
    ) as fwf, open(
        "%svfFrameWindows.dat" % input_path, "w"
    ) as ffw, open(
        "%svfFrameSurface.dat" % input_path, "w"
    ) as ffs, open(
        "%svfFrameFrame.dat" % input_path, "w"
    ) as fff:
        for i in range(numSurf):
            for j in range(numSurf):
                fss.write(
                    "surf%i surf%i %1.6f\n"
                    % (
                        i,
                        j,
                        viewFact[i, j]
                        * areaSurf[j]
                        / (surfVect[j, 0] * surfVect[j, 1]),
                    )
                )
                for k in range(numFrame):
                    if frameVect[k, 3] == j:
                        fsf.write(
                            "surf%i frame%i %1.6f\n"
                            % (
                                i,
                                k,
                                viewFact[i, j]
                                * areaFrame[k]
                                / (surfVect[j, 0] * surfVect[j, 1]),
                            )
                        )
                        for l in range(numWin):
                            if winVect[l, 2] == k:
                                fsw.write(
                                    "surf%i win%i %1.6f\n"
                                    % (
                                        i,
                                        l,
                                        viewFact[i, j]
                                        * areaWin[l]
                                        / (surfVect[j, 0] * surfVect[j, 1]),
                                    )
                                )
        for i in range(numSurf):
            for h in range(numFrame):
                if frameVect[h, 3] == i:
                    for j in range(numSurf):
                        ffs.write(
                            "frame%i surf%i %1.6f\n"
                            % (
                                h,
                                j,
                                viewFact[i, j]
                                * areaSurf[j]
                                / (surfVect[j, 0] * surfVect[j, 1]),
                            )
                        )
                        for k in range(numFrame):
                            if frameVect[k, 3] == j:
                                fff.write(
                                    "frame%i frame%i %1.6f\n"
                                    % (
                                        h,
                                        k,
                                        viewFact[i, j]
                                        * areaFrame[k]
                                        / (surfVect[j, 0] * surfVect[j, 1]),
                                    )
                                )
                                for l in range(numWin):
                                    if winVect[l, 2] == k:
                                        ffw.write(
                                            "frame%i win%i %1.6f\n"
                                            % (
                                                h,
                                                l,
                                                viewFact[i, j]
                                                * areaWin[l]
                                                / (surfVect[j, 0] * surfVect[j, 1]),
                                            )
                                        )
        for i in range(numSurf):
            for m in range(numFrame):
                if frameVect[m, 3] == i:
                    for h in range(numWin):
                        if winVect[h, 2] == m:
                            for j in range(numSurf):
                                fws.write(
                                    "win%i surf%i %1.6f\n"
                                    % (
                                        h,
                                        j,
                                        viewFact[i, j]
                                        * areaSurf[j]
                                        / (surfVect[j, 0] * surfVect[j, 1]),
                                    )
                                )
                                for k in range(numFrame):
                                    if frameVect[k, 3] == j:
                                        fwf.write(
                                            "win%i frame%i %1.6f\n"
                                            % (
                                                h,
                                                k,
                                                viewFact[i, j]
                                                * areaFrame[k]
                                                / (surfVect[j, 0] * surfVect[j, 1]),
                                            )
                                        )
                                        for l in range(numWin):
                                            if winVect[l, 2] == k:
                                                fww.write(
                                                    "win%i win%i %1.6f\n"
                                                    % (
                                                        h,
                                                        l,
                                                        viewFact[i, j]
                                                        * areaWin[l]
                                                        / (
                                                            surfVect[j, 0]
                                                            * surfVect[j, 1]
                                                        ),
                                                    )
                                                )
