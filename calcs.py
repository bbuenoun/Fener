# -*- coding: utf-8 -*-

import numpy as np
from numpy import (
    matrix,
    arctan,
    interp,
    arccos,
    arctan2,
    zeros,
    ones,
    cos,
    tan,
    sin,
    pi,
    radians,
    arange,
    genfromtxt,
)
import calcs
import ioRoutines
import param as p
from scipy.interpolate import interp2d
from math import log, exp, sqrt, degrees, atan, asin, acos

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


def roosModel(z, t0, r0f, r0b, p, q):
    # ------------------------------------------------------------------------
    # Roos model parameters
    a = 8.0
    b = 0.25 / q
    c = 1 - a - b
    alpha = 5.2 + 0.7 * q
    beta = 2.0
    gamma = (5.26 + 0.06 * p) + (0.73 + 0.04 * p) * q
    # ------------------------------------------------------------------------
    # roos transmission model
    tau = t0 * (
        1.0 - a * (z / 90.0) ** alpha - b * (z / 90.0) ** beta - c * (z / 90.0) ** gamma
    )
    # ------------------------------------------------------------------------
    # reflection model (Kuhn 2006)
    a0 = 1 - r0f - t0
    if z < 75:
        rf = 1.0 - tau - a0
    else:
        rf = 1.0 - tau - a0 * (90.0 - z) / 15.0
    a0 = 1 - r0b - t0
    if z < 75:
        rb = 1.0 - tau - a0
    else:
        rb = 1.0 - tau - a0 * (90.0 - z) / 15.0
    return tau, rf, rb


# ========================================================================
# ========================================================================


def fullRoosModel(polarAngleM, g0, t0v, r0fv, r0bv, t0s, r0fs, r0bs, p, q):
    # ------------------------------------------------------------------------
    # Roos model parameters
    a = 8.0
    b = 0.25 / q
    c = 1 - a - b
    alpha = 5.2 + 0.7 * q
    beta = 2.0
    gamma = (5.26 + 0.06 * p) + (0.73 + 0.04 * p) * q
    # ----------------------------------------------------
    # glass properties
    g = zeros(len(polarAngleM))
    tauv = zeros(len(polarAngleM))
    rfv = zeros(len(polarAngleM))
    rbv = zeros(len(polarAngleM))
    taus = zeros(len(polarAngleM))
    rfs = zeros(len(polarAngleM))
    rbs = zeros(len(polarAngleM))
    a0fv = 1 - r0fv - t0v
    a0bv = 1 - r0bv - t0v
    a0fs = 1 - r0fs - t0s
    a0bs = 1 - r0bs - t0s
    for j in range(len(polarAngleM)):
        # roos transmission model
        z = degrees(polarAngleM[j]) / 90.0
        tauv[j] = t0v * (1.0 - a * z ** alpha - b * z ** beta - c * z ** gamma)
        taus[j] = t0s * (1.0 - a * z ** alpha - b * z ** beta - c * z ** gamma)
        g[j] = g0 * (1.0 - a * z ** alpha - b * z ** beta - c * z ** gamma)
        # ------------------------------------------------------------------------
        # reflection model (Kuhn 2006)
        if degrees(polarAngleM[j]) < 75:
            rfv[j] = 1.0 - tauv[j] - a0fv
            rbv[j] = 1.0 - tauv[j] - a0bv
            rfs[j] = 1.0 - taus[j] - a0fs
            rbs[j] = 1.0 - taus[j] - a0bs
        else:
            rfv[j] = 1.0 - tauv[j] - a0fv * (90.0 - degrees(polarAngleM[j])) / 15.0
            rbv[j] = 1.0 - tauv[j] - a0bv * (90.0 - degrees(polarAngleM[j])) / 15.0
            rfs[j] = 1.0 - taus[j] - a0fs * (90.0 - degrees(polarAngleM[j])) / 15.0
            rbs[j] = 1.0 - taus[j] - a0bs * (90.0 - degrees(polarAngleM[j])) / 15.0
    # ----------------------------------------------------
    return g, tauv, rfv, rbv, taus, rfs, rbs


# ========================================================================
# ========================================================================


def twoLayerSys(
    propagationMtx,
    firstLayerTransMtx,
    firstLayerReflFrontMtx,
    firstLayerReflBackMtx,
    scndLayerTransMtx,
    scndLayerReflFrontMtx,
    scndLayerReflBackMtx,
):
    # ------------------------------------------------------------------------
    auxMtx = (
        np.identity(p.nKlem)
        - propagationMtx
        * firstLayerReflBackMtx
        * propagationMtx
        * scndLayerReflFrontMtx
    )
    transMtx = scndLayerTransMtx * auxMtx.I * propagationMtx * firstLayerTransMtx
    auxMtx = (
        np.identity(p.nKlem)
        - propagationMtx
        * scndLayerReflFrontMtx
        * propagationMtx
        * firstLayerReflBackMtx
    )
    reflFrontMtx = (
        firstLayerReflFrontMtx
        + firstLayerTransMtx
        * auxMtx.I
        * propagationMtx
        * scndLayerReflFrontMtx
        * propagationMtx
        * firstLayerTransMtx
    )
    auxMtx = (
        np.identity(p.nKlem)
        - propagationMtx
        * firstLayerReflBackMtx
        * propagationMtx
        * scndLayerReflFrontMtx
    )
    reflBackMtx = (
        scndLayerReflBackMtx
        + scndLayerTransMtx
        * auxMtx.I
        * propagationMtx
        * firstLayerReflBackMtx
        * propagationMtx
        * scndLayerTransMtx
    )
    # ------------------------------------------------------------------------
    return transMtx, reflFrontMtx, reflBackMtx


# ========================================================================
# ========================================================================


def klemsIni():
    # initial parameters
    polarAngle = [
        0,
        5.0,
        15.0,
        25.0,
        35.0,
        45.0,
        55.0,
        65.0,
        75.0,
    ]  # polar angles of the boundaries [deg]
    azimuthDivision = [1, 8, 16, 20, 24, 24, 24, 16, 12]  # azimuth divisions
    polarAngleR = zeros(len(polarAngle))  # polar angles of the boundaries [rad]
    polarAngleM = zeros(len(polarAngle))  # polar angles of the midpoints [rad]
    for i in range(len(polarAngle) - 2):
        polarAngleM[i + 1] = 0.5 * (polarAngle[i + 1] + polarAngle[i + 2]) * pi / 180.0
        polarAngleR[i + 1] = polarAngle[i + 1] * pi / 180.0
    polarAngleM[len(polarAngle) - 1] = 82.5 * pi / 180.0
    polarAngleR[len(polarAngle) - 1] = polarAngle[len(polarAngle) - 1] * pi / 180.0
    #    #------------------------------------------------------------------------
    #    polarAngle = [0,10.,20.,30.,40.,50.,60.,70.,82.5]  # polar angles of the boundaries [deg]
    #    azimuthDivision = [1,8,16,20,24,24,24,16,12]           # azimuth divisions
    #    polarAngleR = zeros(len(polarAngle))               # polar angles of the boundaries [rad]
    #    polarAngleM = zeros(len(polarAngle))               # polar angles of the midpoints [rad]
    #    for i in range(len(polarAngle)-2):
    #           polarAngleM[i+1] = 0.5*(polarAngle[i+1]+polarAngle[i+2])*pi/180.
    #           polarAngleR[i+1] = polarAngle[i+1]*pi/180.
    #    polarAngleM[len(polarAngle)-1] = 86.25*pi/180.
    #    polarAngleR[len(polarAngle)-1] = polarAngle[len(polarAngle)-1]*pi/180.
    # ------------------------------------------------------------------------
    # compute propagation matrix
    sumAzimuthDivision = 0
    propagationMtx = matrix(zeros([p.nKlem, p.nKlem]))
    propagationMtx[0, 0] = (
        2
        * pi
        / azimuthDivision[0]
        * (-cos(polarAngleR[1]) + cos(polarAngleR[0]))
        * cos(polarAngleM[0])
    )
    for i in range(1, len(polarAngle) - 1):
        sumAzimuthDivision = sumAzimuthDivision + azimuthDivision[i - 1]
        for j in range(azimuthDivision[i]):
            k = sumAzimuthDivision + j
            propagationMtx[k, k] = (
                2
                * pi
                / azimuthDivision[i]
                * (-cos(polarAngleR[i + 1]) + cos(polarAngleR[i]))
                * cos(polarAngleM[i])
            )
    sumAzimuthDivision = sumAzimuthDivision + azimuthDivision[len(polarAngle) - 2]
    for j in range(azimuthDivision[len(polarAngle) - 1]):
        k = sumAzimuthDivision + j
        propagationMtx[k, k] = (
            2
            * pi
            / azimuthDivision[len(polarAngle) - 1]
            * (-cos(pi / 2) + cos(polarAngleR[len(polarAngle) - 1]))
            * cos(polarAngleM[len(polarAngle) - 1])
        )
    # ------------------------------------------------------------------------
    return propagationMtx, polarAngleM, azimuthDivision


# ========================================================================
# ========================================================================


def klemsSolidAngles():
    # initial parameters
    polarAngle = [
        0,
        5.0,
        15.0,
        25.0,
        35.0,
        45.0,
        55.0,
        65.0,
        75.0,
    ]  # polar angles of the boundaries [deg]
    azimuthDivision = [1, 8, 16, 20, 24, 24, 24, 16, 12]  # azimuth divisions
    polarAngleR = zeros(len(polarAngle))  # polar angles of the boundaries [rad]
    polarAngleM = zeros(len(polarAngle))  # polar angles of the midpoints [rad]
    for i in range(len(polarAngle) - 2):
        polarAngleM[i + 1] = 0.5 * (polarAngle[i + 1] + polarAngle[i + 2]) * pi / 180.0
        polarAngleR[i + 1] = polarAngle[i + 1] * pi / 180.0
    polarAngleM[len(polarAngle) - 1] = 82.5 * pi / 180.0
    polarAngleR[len(polarAngle) - 1] = polarAngle[len(polarAngle) - 1] * pi / 180.0
    # ------------------------------------------------------------------------
    # compute propagation matrix
    sumAzimuthDivision = 0
    propagationMtx = matrix(zeros([p.nKlem, p.nKlem]))
    propagationMtx[0, 0] = (
        2 * pi / azimuthDivision[0] * (-cos(polarAngleR[1]) + cos(polarAngleR[0]))
    )
    for i in range(1, len(polarAngle) - 1):
        sumAzimuthDivision = sumAzimuthDivision + azimuthDivision[i - 1]
        for j in range(azimuthDivision[i]):
            k = sumAzimuthDivision + j
            propagationMtx[k, k] = (
                2
                * pi
                / azimuthDivision[i]
                * (-cos(polarAngleR[i + 1]) + cos(polarAngleR[i]))
            )
    sumAzimuthDivision = sumAzimuthDivision + azimuthDivision[len(polarAngle) - 2]
    for j in range(azimuthDivision[len(polarAngle) - 1]):
        k = sumAzimuthDivision + j
        propagationMtx[k, k] = (
            2
            * pi
            / azimuthDivision[len(polarAngle) - 1]
            * (-cos(pi / 2) + cos(polarAngleR[len(polarAngle) - 1]))
        )
    # ------------------------------------------------------------------------
    return propagationMtx


def klemsSolidAnglesBis():
    # initial parameters
    polarAngle = [
        0,
        5.0,
        15.0,
        25.0,
        35.0,
        45.0,
        55.0,
        65.0,
        75.0,
    ]  # polar angles of the boundaries [deg]
    azimuthDivision = [1, 8, 16, 20, 24, 24, 24, 16, 12]  # azimuth divisions
    polarAngleR = zeros(len(polarAngle))  # polar angles of the boundaries [rad]
    polarAngleM = zeros(len(polarAngle))  # polar angles of the midpoints [rad]
    for i in range(len(polarAngle) - 2):
        polarAngleM[i + 1] = 0.5 * (polarAngle[i + 1] + polarAngle[i + 2]) * pi / 180.0
        polarAngleR[i + 1] = polarAngle[i + 1] * pi / 180.0
    polarAngleM[len(polarAngle) - 1] = 82.5 * pi / 180.0
    polarAngleR[len(polarAngle) - 1] = polarAngle[len(polarAngle) - 1] * pi / 180.0
    # ------------------------------------------------------------------------
    # compute solid angles
    sumAzimuthDivision = 0
    solidAngles = zeros(p.nKlem)
    solidAngles[0] = (
        2 * pi / azimuthDivision[0] * (-cos(polarAngleR[1]) + cos(polarAngleR[0]))
    )
    for i in range(1, len(polarAngle) - 1):
        sumAzimuthDivision = sumAzimuthDivision + azimuthDivision[i - 1]
        for j in range(azimuthDivision[i]):
            k = sumAzimuthDivision + j
            solidAngles[k] = (
                2
                * pi
                / azimuthDivision[i]
                * (-cos(polarAngleR[i + 1]) + cos(polarAngleR[i]))
            )
    sumAzimuthDivision = sumAzimuthDivision + azimuthDivision[len(polarAngle) - 2]
    for j in range(azimuthDivision[len(polarAngle) - 1]):
        k = sumAzimuthDivision + j
        solidAngles[k] = (
            2
            * pi
            / azimuthDivision[len(polarAngle) - 1]
            * (-cos(pi / 2) + cos(polarAngleR[len(polarAngle) - 1]))
        )
    # ------------------------------------------------------------------------
    return solidAngles


# ========================================================================
# ========================================================================


def biDirShade(shadTrans, shadRefl):
    thrdLayerTransMtx = matrix(ones([p.nKlem, p.nKlem])) * shadTrans / pi
    thrdLayerReflFrontMtx = matrix(ones([p.nKlem, p.nKlem])) * shadRefl / pi
    thrdLayerReflBackMtx = matrix(ones([p.nKlem, p.nKlem])) * shadRefl / pi
    return thrdLayerTransMtx, thrdLayerReflFrontMtx, thrdLayerReflBackMtx


# ========================================================================
# ========================================================================


def biDirGlazing(trans, reflFront, reflBack, qRoos):
    # run Klems model initialization
    propagationMtx, polarAngleM, azimuthDivision = calcs.klemsIni()
    # ------------------------------------------------------------------------
    # glass properties
    layerSolTrans = zeros(len(polarAngleM))
    layerSolReflFront = zeros(len(polarAngleM))
    layerSolReflBack = zeros(len(polarAngleM))
    for j in range(len(polarAngleM)):
        layerSolTrans[j], layerSolReflFront[j], layerSolReflBack[j] = calcs.roosModel(
            polarAngleM[j] * 180 / pi, trans, reflFront, reflBack, 1, qRoos
        )
    # ------------------------------------------------------------------------
    # compute bi-directional transmittance and reflectance matrices
    layerTransMtx = matrix(zeros([p.nKlem, p.nKlem]))
    layerReflFrontMtx = matrix(zeros([p.nKlem, p.nKlem]))
    layerReflBackMtx = matrix(zeros([p.nKlem, p.nKlem]))
    layerTransMtx[0, 0] = layerSolTrans[0] / propagationMtx[0, 0]
    layerReflFrontMtx[0, 0] = layerSolReflFront[0] / propagationMtx[0, 0]
    layerReflBackMtx[0, 0] = layerSolReflBack[0] / propagationMtx[0, 0]
    sumAzimuthDivision = 0
    for i in range(1, len(polarAngleM)):
        sumAzimuthDivision = sumAzimuthDivision + azimuthDivision[i - 1]
        for j in range(azimuthDivision[i]):
            k = sumAzimuthDivision + j
            layerTransMtx[k, k] = layerSolTrans[i] / propagationMtx[k, k]
            layerReflFrontMtx[k, k] = layerSolReflFront[i] / propagationMtx[k, k]
            layerReflBackMtx[k, k] = layerSolReflBack[i] / propagationMtx[k, k]
    # ------------------------------------------------------------------------
    return layerTransMtx, layerReflFrontMtx, layerReflBackMtx


# ========================================================================
# ========================================================================


def biDirBlinds(tauFile, rhoFile):
    # run Klems model initialization
    propagationMtx, polarAngleM, azimuthDivision = calcs.klemsIni()
    # ------------------------------------------------------------------------
    # read optical data as a function of 'profile angle'
    tauProp = genfromtxt(tauFile, skip_header=0, delimiter="\t")
    rhoProp = genfromtxt(rhoFile, skip_header=0, delimiter="\t")
    # ------------------------------------------------------------------------
    # compute bi-directional transmittance and reflectance matrices
    layerTransMtx = matrix(zeros([p.nKlem, p.nKlem]))
    layerReflFrontMtx = matrix(zeros([p.nKlem, p.nKlem]))
    layerReflBackMtx = matrix(zeros([p.nKlem, p.nKlem]))
    for n in range(p.nKlem):
        layerTransMtx[0, n] = tauProp[17, 2] / pi
        layerReflFrontMtx[0, n] = rhoProp[17, 2] / pi
        layerReflBackMtx[0, n] = rhoProp[17, 2] / pi
    sumAzimuthDivision = 0
    for i in range(1, len(polarAngleM)):
        sumAzimuthDivision = sumAzimuthDivision + azimuthDivision[i - 1]
        for j in range(azimuthDivision[i]):
            k = sumAzimuthDivision + j
            profAngle = arctan(
                sin(2 * pi * j / azimuthDivision[i]) * tan(polarAngleM[i])
            )
            transPatch = calcs.interpOptProp(profAngle, tauProp)
            reflPatch = calcs.interpOptProp(profAngle, rhoProp)
            for n in range(p.nKlem):
                layerTransMtx[k, n] = transPatch / pi
                layerReflFrontMtx[k, n] = reflPatch / pi
                layerReflBackMtx[k, n] = reflPatch / pi
    # ------------------------------------------------------------------------
    return layerTransMtx, layerReflFrontMtx, layerReflBackMtx


# ========================================================================
# ========================================================================


def interpOptProp(profAngle, optProp):
    profAngleDeg = profAngle * 180 / pi
    # ------------------------------------------------------------------------
    if profAngleDeg == -90.0:
        propPatch = optProp[len(optProp[:, 2]) - 1, 2]
    else:
        for i in range(len(optProp[:, 0]) - 1):
            if profAngleDeg <= optProp[i, 0] and profAngleDeg > optProp[i + 1, 0]:
                propPatch = (profAngleDeg - optProp[i + 1, 0]) / (
                    optProp[i, 0] - optProp[i + 1, 0]
                ) * (optProp[i, 2] - optProp[i + 1, 2]) + optProp[i + 1, 2]
    # ------------------------------------------------------------------------
    return propPatch


# ========================================================================
# ========================================================================


def genBSDF(nLayers, transMtx, reflFrontMtx, reflBackMtx):
    # ------------------------------------------------------------------------
    # run Klems model initialization
    propagationMtx, polarAngleM, azimuthDivision = calcs.klemsIni()
    # ------------------------------------------------------------------------
    # read bi-directional transmittance and reflectance matrices for each layer

    firstLayerTransMtx = transMtx[0]
    firstLayerReflFrontMtx = reflFrontMtx[0]
    firstLayerReflBackMtx = reflBackMtx[0]

    if nLayers > 1:
        scndLayerTransMtx = transMtx[1]
        scndLayerReflFrontMtx = reflFrontMtx[1]
        scndLayerReflBackMtx = reflBackMtx[1]

    if nLayers > 2:
        thrdLayerTransMtx = transMtx[2]
        thrdLayerReflFrontMtx = reflFrontMtx[2]
        thrdLayerReflBackMtx = reflBackMtx[2]

        if nLayers > 3:
            fourLayerTransMtx = transMtx[3]
            fourLayerReflFrontMtx = reflFrontMtx[3]
            fourLayerReflBackMtx = reflBackMtx[3]
    # ------------------------------------------------------------------------
    # compute directional-hemispherical transmittance, reflectance, and absorptance
    auxVec = ones(p.nKlem)
    firstLayerTransDirHem = matrix(zeros([1, p.nKlem]))
    firstLayerReflFrontDirHem = matrix(zeros([1, p.nKlem]))
    firstLayerReflBackDirHem = matrix(zeros([1, p.nKlem]))
    firstLayerAbsFrontDirHem = matrix(zeros([1, p.nKlem]))
    firstLayerAbsBackDirHem = matrix(zeros([1, p.nKlem]))
    firstLayerTransDirHem = auxVec * propagationMtx * firstLayerTransMtx
    firstLayerReflFrontDirHem = auxVec * propagationMtx * firstLayerReflFrontMtx
    firstLayerReflBackDirHem = auxVec * propagationMtx * firstLayerReflBackMtx
    firstLayerAbsFrontDirHem = (
        auxVec - firstLayerTransDirHem - firstLayerReflFrontDirHem
    )
    firstLayerAbsBackDirHem = auxVec - firstLayerTransDirHem - firstLayerReflBackDirHem
    if nLayers > 1:
        scndLayerTransDirHem = matrix(zeros([1, p.nKlem]))
        scndLayerReflFrontDirHem = matrix(zeros([1, p.nKlem]))
        scndLayerReflBackDirHem = matrix(zeros([1, p.nKlem]))
        scndLayerAbsFrontDirHem = matrix(zeros([1, p.nKlem]))
        scndLayerAbsBackDirHem = matrix(zeros([1, p.nKlem]))
        scndLayerTransDirHem = auxVec * propagationMtx * scndLayerTransMtx
        scndLayerReflFrontDirHem = auxVec * propagationMtx * scndLayerReflFrontMtx
        scndLayerReflBackDirHem = auxVec * propagationMtx * scndLayerReflBackMtx
        scndLayerAbsFrontDirHem = (
            auxVec - scndLayerTransDirHem - scndLayerReflFrontDirHem
        )
        scndLayerAbsBackDirHem = auxVec - scndLayerTransDirHem - scndLayerReflBackDirHem
    if nLayers > 2:
        thrdLayerTransDirHem = matrix(zeros([1, p.nKlem]))
        thrdLayerReflFrontDirHem = matrix(zeros([1, p.nKlem]))
        thrdLayerReflBackDirHem = matrix(zeros([1, p.nKlem]))
        thrdLayerAbsFrontDirHem = matrix(zeros([1, p.nKlem]))
        thrdLayerAbsBackDirHem = matrix(zeros([1, p.nKlem]))
        thrdLayerTransDirHem = auxVec * propagationMtx * thrdLayerTransMtx
        thrdLayerReflFrontDirHem = auxVec * propagationMtx * thrdLayerReflFrontMtx
        thrdLayerReflBackDirHem = auxVec * propagationMtx * thrdLayerReflBackMtx
        thrdLayerAbsFrontDirHem = (
            auxVec - thrdLayerTransDirHem - thrdLayerReflFrontDirHem
        )
        thrdLayerAbsBackDirHem = auxVec - thrdLayerTransDirHem - thrdLayerReflBackDirHem

        if nLayers > 3:
            fourLayerTransDirHem = matrix(zeros([1, p.nKlem]))
            fourLayerReflFrontDirHem = matrix(zeros([1, p.nKlem]))
            fourLayerReflBackDirHem = matrix(zeros([1, p.nKlem]))
            fourLayerAbsFrontDirHem = matrix(zeros([1, p.nKlem]))
            fourLayerAbsBackDirHem = matrix(zeros([1, p.nKlem]))
            fourLayerTransDirHem = auxVec * propagationMtx * fourLayerTransMtx
            fourLayerReflFrontDirHem = auxVec * propagationMtx * fourLayerReflFrontMtx
            fourLayerReflBackDirHem = auxVec * propagationMtx * fourLayerReflBackMtx
            fourLayerAbsFrontDirHem = (
                auxVec - fourLayerTransDirHem - fourLayerReflFrontDirHem
            )
            fourLayerAbsBackDirHem = (
                auxVec - fourLayerTransDirHem - fourLayerReflBackDirHem
            )
    # ------------------------------------------------------------------------
    # compute multi-layer transmittance and reflectance matrices
    if nLayers > 1:
        (
            firstScndTransMtx,
            firstScndReflFrontMtx,
            firstScndReflBackMtx,
        ) = calcs.twoLayerSys(
            propagationMtx,
            firstLayerTransMtx,
            firstLayerReflFrontMtx,
            firstLayerReflBackMtx,
            scndLayerTransMtx,
            scndLayerReflFrontMtx,
            scndLayerReflBackMtx,
        )
    if nLayers > 2:
        scndThrdTransMtx, scndThrdReflFrontMtx, scndThrdReflBackMtx = calcs.twoLayerSys(
            propagationMtx,
            scndLayerTransMtx,
            scndLayerReflFrontMtx,
            scndLayerReflBackMtx,
            thrdLayerTransMtx,
            thrdLayerReflFrontMtx,
            thrdLayerReflBackMtx,
        )

        (
            firstScndThrdTransMtx,
            firstScndThrdReflFrontMtx,
            firstScndThrdReflBackMtx,
        ) = calcs.twoLayerSys(
            propagationMtx,
            firstScndTransMtx,
            firstScndReflFrontMtx,
            firstScndReflBackMtx,
            thrdLayerTransMtx,
            thrdLayerReflFrontMtx,
            thrdLayerReflBackMtx,
        )
        if nLayers > 3:
            (
                scndThrdFourTransMtx,
                scndThrdFourReflFrontMtx,
                scndThrdFourReflBackMtx,
            ) = calcs.twoLayerSys(
                propagationMtx,
                scndThrdTransMtx,
                scndThrdReflFrontMtx,
                scndThrdReflBackMtx,
                fourLayerTransMtx,
                fourLayerReflFrontMtx,
                fourLayerReflBackMtx,
            )

            (
                thrdFourTransMtx,
                thrdFourReflFrontMtx,
                thrdFourReflBackMtx,
            ) = calcs.twoLayerSys(
                propagationMtx,
                thrdLayerTransMtx,
                thrdLayerReflFrontMtx,
                thrdLayerReflBackMtx,
                fourLayerTransMtx,
                fourLayerReflFrontMtx,
                fourLayerReflBackMtx,
            )

            (
                firstScndThrdFourTransMtx,
                firstScndThrdFourReflFrontMtx,
                firstScndThrdFourReflBackMtx,
            ) = calcs.twoLayerSys(
                propagationMtx,
                firstScndThrdTransMtx,
                firstScndThrdReflFrontMtx,
                firstScndThrdReflBackMtx,
                fourLayerTransMtx,
                fourLayerReflFrontMtx,
                fourLayerReflBackMtx,
            )
    # ------------------------------------------------------------------------
    # compute layer/system absorptance
    firstLayerSysAbsFrontVec = matrix(zeros([1, p.nKlem]))
    firstLayerSysAbsBackVec = matrix(zeros([1, p.nKlem]))
    scndLayerSysAbsFrontVec = matrix(zeros([1, p.nKlem]))
    scndLayerSysAbsBackVec = matrix(zeros([1, p.nKlem]))
    if nLayers == 1:
        firstLayerSysAbsFrontVec = firstLayerAbsFrontDirHem
        firstLayerSysAbsBackVec = firstLayerAbsBackDirHem
    elif nLayers == 2:
        auxMtx = (
            np.identity(p.nKlem)
            - propagationMtx
            * scndLayerReflFrontMtx
            * propagationMtx
            * firstLayerReflBackMtx
        )
        firstLayerSysAbsFrontVec = (
            firstLayerAbsFrontDirHem
            + firstLayerAbsBackDirHem
            * auxMtx.I
            * propagationMtx
            * scndLayerReflFrontMtx
            * propagationMtx
            * firstLayerTransMtx
        )
        firstLayerSysAbsBackVec = (
            firstLayerAbsBackDirHem * auxMtx.I * propagationMtx * scndLayerTransMtx
        )

        auxMtx = (
            np.identity(p.nKlem)
            - propagationMtx
            * firstLayerReflBackMtx
            * propagationMtx
            * scndLayerReflFrontMtx
        )
        scndLayerSysAbsFrontVec = (
            scndLayerAbsFrontDirHem * auxMtx.I * propagationMtx * firstLayerTransMtx
        )
        scndLayerSysAbsBackVec = (
            scndLayerAbsFrontDirHem
            * auxMtx.I
            * propagationMtx
            * firstLayerReflBackMtx
            * propagationMtx
            * scndLayerTransMtx
            + scndLayerAbsBackDirHem
        )
    elif nLayers == 3:
        auxMtx = (
            np.identity(p.nKlem)
            - propagationMtx
            * scndThrdReflFrontMtx
            * propagationMtx
            * firstLayerReflBackMtx
        )
        firstLayerSysAbsFrontVec = (
            firstLayerAbsFrontDirHem
            + firstLayerAbsBackDirHem
            * auxMtx.I
            * propagationMtx
            * scndThrdReflFrontMtx
            * propagationMtx
            * firstLayerTransMtx
        )
        firstLayerSysAbsBackVec = (
            firstLayerAbsBackDirHem * auxMtx.I * propagationMtx * scndThrdTransMtx
        )

        auxMtx1 = (
            np.identity(p.nKlem)
            - propagationMtx
            * firstLayerReflBackMtx
            * propagationMtx
            * scndThrdReflFrontMtx
        )
        auxMtx2 = (
            np.identity(p.nKlem)
            - propagationMtx
            * thrdLayerReflFrontMtx
            * propagationMtx
            * firstScndReflBackMtx
        )
        scndLayerSysAbsFrontVec = (
            scndLayerAbsFrontDirHem * auxMtx1.I * propagationMtx * firstLayerTransMtx
            + scndLayerAbsBackDirHem
            * auxMtx2.I
            * propagationMtx
            * thrdLayerReflFrontMtx
            * propagationMtx
            * firstScndTransMtx
        )
        scndLayerSysAbsBackVec = (
            scndLayerAbsFrontDirHem
            * auxMtx1.I
            * propagationMtx
            * firstLayerReflBackMtx
            * propagationMtx
            * scndThrdTransMtx
            + scndLayerAbsBackDirHem * auxMtx2.I * propagationMtx * thrdLayerTransMtx
        )

        thrdLayerSysAbsFrontVec = matrix(zeros([1, p.nKlem]))
        thrdLayerSysAbsBackVec = matrix(zeros([1, p.nKlem]))

        auxMtx = (
            np.identity(p.nKlem)
            - propagationMtx
            * firstScndReflBackMtx
            * propagationMtx
            * thrdLayerReflFrontMtx
        )
        thrdLayerSysAbsFrontVec = (
            thrdLayerAbsFrontDirHem * auxMtx.I * propagationMtx * firstScndTransMtx
        )
        thrdLayerSysAbsBackVec = (
            thrdLayerAbsBackDirHem
            + thrdLayerAbsFrontDirHem
            * auxMtx.I
            * propagationMtx
            * firstScndReflBackMtx
            * propagationMtx
            * thrdLayerTransMtx
        )

    elif nLayers == 4:
        auxMtx = (
            np.identity(p.nKlem)
            - propagationMtx
            * scndThrdFourReflFrontMtx
            * propagationMtx
            * firstLayerReflBackMtx
        )
        firstLayerSysAbsFrontVec = (
            firstLayerAbsFrontDirHem
            + firstLayerAbsBackDirHem
            * auxMtx.I
            * propagationMtx
            * scndThrdFourReflFrontMtx
            * propagationMtx
            * firstLayerTransMtx
        )
        firstLayerSysAbsBackVec = (
            firstLayerAbsBackDirHem * auxMtx.I * propagationMtx * scndThrdFourTransMtx
        )

        auxMtx1 = (
            np.identity(p.nKlem)
            - propagationMtx
            * firstLayerReflBackMtx
            * propagationMtx
            * scndThrdFourReflFrontMtx
        )
        auxMtx2 = (
            np.identity(p.nKlem)
            - propagationMtx
            * thrdFourReflFrontMtx
            * propagationMtx
            * firstScndReflBackMtx
        )
        scndLayerSysAbsFrontVec = (
            scndLayerAbsFrontDirHem * auxMtx1.I * propagationMtx * firstLayerTransMtx
            + scndLayerAbsBackDirHem
            * auxMtx2.I
            * propagationMtx
            * thrdFourReflFrontMtx
            * propagationMtx
            * firstScndTransMtx
        )
        scndLayerSysAbsBackVec = (
            scndLayerAbsFrontDirHem
            * auxMtx1.I
            * propagationMtx
            * firstLayerReflBackMtx
            * propagationMtx
            * scndThrdFourTransMtx
            + scndLayerAbsBackDirHem * auxMtx2.I * propagationMtx * thrdFourTransMtx
        )

        thrdLayerSysAbsFrontVec = matrix(zeros([1, p.nKlem]))
        thrdLayerSysAbsBackVec = matrix(zeros([1, p.nKlem]))
        auxMtx1 = (
            np.identity(p.nKlem)
            - propagationMtx
            * firstScndReflBackMtx
            * propagationMtx
            * thrdFourReflFrontMtx
        )
        auxMtx2 = (
            np.identity(p.nKlem)
            - propagationMtx
            * fourLayerReflFrontMtx
            * propagationMtx
            * firstScndThrdReflBackMtx
        )
        thrdLayerSysAbsFrontVec = (
            thrdLayerAbsFrontDirHem * auxMtx1.I * propagationMtx * firstScndTransMtx
            + thrdLayerAbsBackDirHem
            * auxMtx2.I
            * propagationMtx
            * fourLayerReflFrontMtx
            * propagationMtx
            * firstScndThrdTransMtx
        )
        thrdLayerSysAbsBackVec = (
            thrdLayerAbsFrontDirHem
            * auxMtx1.I
            * propagationMtx
            * firstScndReflBackMtx
            * propagationMtx
            * thrdFourTransMtx
            + thrdLayerAbsBackDirHem * auxMtx2.I * propagationMtx * fourLayerTransMtx
        )

        auxMtx = (
            np.identity(p.nKlem)
            - propagationMtx
            * fourLayerReflFrontMtx
            * propagationMtx
            * firstScndThrdReflBackMtx
        )
        fourLayerSysAbsFrontVec = (
            fourLayerAbsFrontDirHem * auxMtx.I * propagationMtx * firstScndThrdTransMtx
        )
        fourLayerSysAbsBackVec = (
            fourLayerAbsBackDirHem
            + fourLayerAbsFrontDirHem
            * auxMtx.I
            * propagationMtx
            * firstScndThrdReflBackMtx
            * propagationMtx
            * fourLayerTransMtx
        )
    else:
        print("ERROR: maximum four-layer system implemented")
    # ------------------------------------------------------------------------
    # output
    if nLayers == 1:
        return (
            firstLayerTransMtx,
            firstLayerReflFrontMtx,
            firstLayerReflBackMtx,
            [firstLayerSysAbsFrontVec],
            [firstLayerSysAbsBackVec],
        )
    elif nLayers == 2:
        return (
            firstScndTransMtx,
            firstScndReflFrontMtx,
            firstScndReflBackMtx,
            [firstLayerSysAbsFrontVec, scndLayerSysAbsFrontVec],
            [firstLayerSysAbsBackVec, scndLayerSysAbsBackVec],
        )
    elif nLayers == 3:
        return (
            firstScndThrdTransMtx,
            firstScndThrdReflFrontMtx,
            firstScndThrdReflBackMtx,
            [
                firstLayerSysAbsFrontVec,
                scndLayerSysAbsFrontVec,
                thrdLayerSysAbsFrontVec,
            ],
            [firstLayerSysAbsBackVec, scndLayerSysAbsBackVec, thrdLayerSysAbsBackVec],
        )
    elif nLayers == 4:
        return (
            firstScndThrdFourTransMtx,
            firstScndThrdFourReflFrontMtx,
            firstScndThrdFourReflBackMtx,
            [
                firstLayerSysAbsFrontVec,
                scndLayerSysAbsFrontVec,
                thrdLayerSysAbsFrontVec,
                fourLayerSysAbsFrontVec,
            ],
            [
                firstLayerSysAbsBackVec,
                scndLayerSysAbsBackVec,
                thrdLayerSysAbsBackVec,
                fourLayerSysAbsBackVec,
            ],
        )
    else:
        print("ERROR: maximum four-layer system implemented")


# ========================================================================
# ========================================================================


def photocellLightControl(illuSp, illuDaylight, control, param):
    # ------------------------------------------------------------------------
    fracLightOutput = max(0, (illuSp - illuDaylight) / illuSp)
    # ------------------------------------------------------------------------
    # continuous dimming control {param-minimum light output fraction}
    if control < 0.5:
        if fracLightOutput < param:
            fracPowerInput = 0
        elif fracLightOutput >= param and fracLightOutput <= 1:
            fracPowerInput = (fracLightOutput - param) / (1 - param)
        else:
            print("ERROR: wrong lighting control parameters")
    # ------------------------------------------------------------------------
    # stepped control {param-number of steps}
    elif control > 0.5:
        if fracLightOutput == 0:
            fracPowerInput = 0
        elif fracLightOutput > 0 and fracLightOutput < 1:
            fracPowerInput = (int(param * fracLightOutput) + 1) / param
        elif fracLightOutput == 1:
            fracPowerInput = 1
        else:
            print("ERROR: wrong lighting control parameters")
    else:
        print("ERROR: wrong lighting control parameters")
    # ------------------------------------------------------------------------
    return fracPowerInput


# ========================================================================
# ========================================================================


def schVector(simInp, sch):
    schMultiplier = zeros(8904)
    it = 0
    while it < 8904:
        if simInp[4] < 6:
            for d in range(6 - simInp[4]):
                for h in range(int(sch[1, 0])):
                    schMultiplier[it] = sch[1, 1]
                    it = it + 1
                for hstp in range(int(sch[0, 0]) - 1):
                    for h in range(int(sch[hstp + 2, 0]) - int(sch[hstp + 1, 0])):
                        schMultiplier[it] = sch[hstp + 2, 1]
                        it = it + 1
            for d in range(2):
                for h in range(int(sch[1 + sch[0, 0], 0])):
                    schMultiplier[it] = sch[1 + sch[0, 0], 1]
                    it = it + 1
                for hstp in range(int(sch[0, 1]) - 1):
                    for h in range(
                        int(sch[sch[0, 0] + hstp + 2, 0])
                        - int(sch[sch[0, 0] + hstp + 1, 0])
                    ):
                        schMultiplier[it] = sch[sch[0, 0] + hstp + 2, 1]
                        it = it + 1
            for d in range(simInp[4] - 1):
                for h in range(int(sch[1, 0])):
                    schMultiplier[it] = sch[1, 1]
                    it = it + 1
                for hstp in range(int(sch[0, 0]) - 1):
                    for h in range(int(sch[hstp + 2, 0]) - int(sch[hstp + 1, 0])):
                        schMultiplier[it] = sch[hstp + 2, 1]
                        it = it + 1
        else:
            for d in range(8 - simInp[4]):
                for h in range(int(sch[1 + sch[0, 0], 0])):
                    schMultiplier[it] = sch[1 + sch[0, 0], 1]
                    it = it + 1
                for hstp in range(int(sch[0, 1]) - 1):
                    for h in range(
                        int(sch[sch[0, 0] + hstp + 2, 0])
                        - int(sch[sch[0, 0] + hstp + 1, 0])
                    ):
                        schMultiplier[it] = sch[sch[0, 0] + hstp + 2, 1]
                        it = it + 1
            for d in range(5):
                for h in range(int(sch[1, 0])):
                    schMultiplier[it] = sch[1, 1]
                    it = it + 1
                for hstp in range(int(sch[0, 0]) - 1):
                    for h in range(int(sch[hstp + 2, 0]) - int(sch[hstp + 1, 0])):
                        schMultiplier[it] = sch[hstp + 2, 1]
                        it = it + 1
            for d in range(simInp[4] - 6):
                for h in range(int(sch[1 + sch[0, 0], 0])):
                    schMultiplier[it] = sch[1 + sch[0, 0], 1]
                    it = it + 1
                for hstp in range(int(sch[0, 1]) - 1):
                    for h in range(
                        int(sch[sch[0, 0] + hstp + 2, 0])
                        - int(sch[sch[0, 0] + hstp + 1, 0])
                    ):
                        schMultiplier[it] = sch[sch[0, 0] + hstp + 2, 1]
                        it = it + 1

    # schMultiplierOut = zeros(8904)
    # for i in range(8760):
    #       schMultiplierOut[i] = schMultiplier[i]
    return schMultiplier


# ========================================================================
# ========================================================================


def shadeCntrl(cntrl, sensor):
    for i in range(int(cntrl[2])):
        if sensor > cntrl[3 + i * 2]:
            cntrlFlag = int(cntrl[4 + i * 2])
    return cntrlFlag


# ========================================================================
# ========================================================================


def rho(temp):
    pres = 101325
    hum = calcs.humFromRHumTemp(50, temp, pres)
    rho = pres / (1000 * 0.287042 * temp * (1.0 + 1.607858 * hum))
    return rho


# ========================================================================
# ========================================================================


def humFromRHumTemp(RH, T, P):
    # Saturation vapour pressure from ASHRAE
    C8 = -5.8002206e3
    C9 = 1.3914993
    C10 = -4.8640239e-2
    C11 = 4.1764768e-5
    C12 = -1.4452093e-8
    C13 = 6.5459673
    PWS = exp(C8 / T + C9 + C10 * T + C11 * T ** 2 + C12 * T ** 3 + C13 * log(T))
    # Vapour pressure
    PW = RH * PWS / 100
    # Specific humidity
    W = 0.62198 * PW / (P - PW)
    return W


# ========================================================================
# ========================================================================


def viewFactorAdjRectangles(a, b, c):
    # Mills 99, p.543
    # a - other side of 2
    # b - common side
    # c - other side of 1
    L = c / b
    N = a / b
    F12 = (
        1.0
        / pi
        / L
        * (
            L * arctan(1.0 / L)
            + N * arctan(1.0 / N)
            - (N ** 2.0 + L ** 2.0) ** 0.5 * arctan((N ** 2.0 + L ** 2.0) ** (-0.5))
            + 0.25
            * log(
                ((1.0 + L ** 2.0) * (1.0 + N ** 2.0) / (1.0 + L ** 2.0 + N ** 2.0))
                * (
                    L ** 2.0
                    * (1.0 + L ** 2.0 + N ** 2.0)
                    / (1.0 + L ** 2.0)
                    / (L ** 2.0 + N ** 2.0)
                )
                ** (L ** 2.0)
                * (
                    N ** 2.0
                    * (1.0 + L ** 2.0 + N ** 2.0)
                    / (1.0 + N ** 2.0)
                    / (L ** 2.0 + N ** 2.0)
                )
                ** (N ** 2.0)
            )
        )
    )
    return F12


# ========================================================================
# ========================================================================


def viewFactorOppRectangles(a, b, c):
    # Mills 99, p.543
    # a - length
    # b - width
    # c - distance between rectangles
    X = b / c
    Y = a / c
    F12 = (
        2.0
        / pi
        / X
        / Y
        * (
            log(
                ((1.0 + X ** 2.0) * (1.0 + Y ** 2.0) / (1.0 + X ** 2.0 + Y ** 2.0))
                ** 0.5
            )
            - X * arctan(X)
            - Y * arctan(Y)
            + X * (1.0 + Y ** 2.0) ** 0.5 * arctan(X / (1.0 + Y ** 2.0) ** 0.5)
            + Y * (1.0 + X ** 2.0) ** 0.5 * arctan(Y / (1.0 + X ** 2.0) ** 0.5)
        )
    )
    return F12


# ========================================================================
# ========================================================================


def chtcGlazCav(widthGap, ventGap, frac, temp1, temp2, areaWin, vertical):
    # ----------------------------------------------------
    if ventGap < 1:
        global temp
        temp = (temp1 + temp2) / 2
        if frac[0] == 1 and frac[1] == 0:
            densGap, visGap, condGap, cpGap, mMGap = calcs.gasProp(0, temp)
        elif frac[0] == 0 and frac[1] == 1:
            densGap, visGap, condGap, cpGap, mMGap = calcs.gasProp(1, temp)
        else:
            densGap, visGap, condGap, cpGap = calcs.gasMixProp(frac, temp)
        # ----------------------------------------------------
        rayleigh = (
            densGap ** 2.0
            * widthGap ** 3.0
            * p.g
            * cpGap
            * abs(temp1 - temp2)
            / visGap
            / condGap
            / temp
        )
        # ----------------------------------------------------
        # Vertical cavity
        if vertical:
            # According to ISO15099
            if rayleigh > 5e4:
                nu1 = 0.0673838 * rayleigh ** (1.0 / 3.0)
            elif rayleigh > 1e4:
                nu1 = 0.028154 * rayleigh ** 0.4134
            else:
                nu1 = 1 + 1.7596678e-10 * rayleigh ** 2.2984755
            aspRat = sqrt(areaWin) / widthGap
            nu2 = 0.242 * (rayleigh / aspRat) ** 0.272
            nusselt = max(nu1, nu2)
        # According to Mills
        #   if rayleigh>0:
        #      nu1 = 0.0605*rayleigh**(1./3.)
        #      nu2 = (1+(0.104*rayleigh**0.293/(1+(6310/rayleigh)**1.36))**3)**(1./3.)
        #      nu3 = 0.242*(rayleigh/sqrt(areaWin)*widthGap)**0.272
        #      nusselt = max(nu1,nu2)
        #      nusselt = max(nusselt,nu3)
        #   else:
        #      nusselt = 0.0
        # ----------------------------------------------------
        # Horizontal cavity
        else:
            aux1 = 0.0
            if rayleigh > 0.01:
                aux1 = ((1 - 1708.0 / rayleigh) + abs((1 - 1708.0 / rayleigh))) / 2.0
            aux2 = (
                ((rayleigh / 5830.0) ** (1.0 / 3.0) - 1)
                + abs(((rayleigh / 5830.0) ** (1.0 / 3.0) - 1.0))
            ) / 2.0
            nusselt = 1.0 + 1.44 * aux1 + aux2
        # ----------------------------------------------------
        # convective heat transfer coefficent
        chtc = nusselt * condGap / widthGap
    else:
        chtc = 25.0
    # ----------------------------------------------------
    return chtc


# ========================================================================
# ========================================================================


def gasProp(gas, temp):
    if gas == 0:  # air
        aCon = 2.873e-3
        bCon = 7.760e-5
        aVis = 3.723e-6
        bVis = 4.94e-8
        aCp = 1002.737
        bCp = 1.2324e-2
        mM = 0.02897
    if gas == 1:  # argon
        aCon = 2.285e-3
        bCon = 5.149e-5
        aVis = 3.379e-6
        bVis = 6.451e-8
        aCp = 521.9285
        bCp = 0
        mM = 0.039948
    densGap = 101300 * mM / p.r / temp
    visGap = aVis + bVis * temp
    condGap = aCon + bCon * temp
    cpGap = aCp + bCp * temp
    return densGap, visGap, condGap, cpGap, mM


# ========================================================================
# ========================================================================


def gasMixProp(frac, temp):
    numGas = 2
    dens = zeros(numGas)
    vis = zeros(numGas)
    cond = zeros(numGas)
    cp = zeros(numGas)
    mM = zeros(numGas)
    dens[0], vis[0], cond[0], cp[0], mM[0] = calcs.gasProp(0, temp)  # air
    dens[1], vis[1], cond[1], cp[1], mM[1] = calcs.gasProp(1, temp)  # argon
    mMGap = 0
    cpGap = 0
    for i in range(numGas):
        mMGap = mMGap + frac[i] * mM[i]
        cpGap = cpGap + frac[i] * mM[i] * cp[i]
    cpGap = cpGap / mMGap
    densGap = 101300 * mMGap / p.r / temp
    phi_01 = (
        (1 + (vis[0] / vis[1]) ** 0.5 * (mM[1] / mM[0]) ** 0.25) ** 2
        / 2
        / sqrt(2)
        / (1 + mM[0] / mM[1]) ** 0.5
    )
    phi_10 = (
        (1 + (vis[1] / vis[0]) ** 0.5 * (mM[0] / mM[1]) ** 0.25) ** 2
        / 2
        / sqrt(2)
        / (1 + mM[1] / mM[0]) ** 0.5
    )
    visGap = vis[1] / (1 + phi_10 * frac[0] / frac[1]) + vis[0] / (
        1 + phi_01 * frac[1] / frac[0]
    )
    cond_p = zeros(numGas)
    cond_pp = zeros(numGas)
    for i in range(numGas):
        cond_p[i] = 15 / 4 * p.r / mM[i] * vis[i]
        cond_pp[i] = cond[i] - cond_p[i]
    psi_01 = (
        (1 + (cond[0] / cond[1]) ** 0.5 * (mM[1] / mM[0]) ** 0.25) ** 2
        / 2
        / sqrt(2)
        / (1 + mM[0] / mM[1]) ** 0.5
        * (1 + 2.41 * (mM[0] - mM[1]) * (mM[0] - 0.142 * mM[1]) / (mM[0] + mM[1]) ** 2)
    )
    psi_10 = (
        (1 + (cond[1] / cond[0]) ** 0.5 * (mM[0] / mM[1]) ** 0.25) ** 2
        / 2
        / sqrt(2)
        / (1 + mM[1] / mM[0]) ** 0.5
        * (1 + 2.41 * (mM[1] - mM[0]) * (mM[1] - 0.142 * mM[0]) / (mM[1] + mM[0]) ** 2)
    )
    phi_01 = (
        (1 + (cond[0] / cond[1]) ** 0.5 * (mM[1] / mM[0]) ** 0.25) ** 2
        / 2
        / sqrt(2)
        / (1 + mM[0] / mM[1]) ** 0.5
    )
    phi_10 = (
        (1 + (cond[1] / cond[0]) ** 0.5 * (mM[0] / mM[1]) ** 0.25) ** 2
        / 2
        / sqrt(2)
        / (1 + mM[1] / mM[0]) ** 0.5
    )
    condGap_p = cond_p[1] / (1 + psi_10 * frac[0] / frac[1]) + cond_p[0] / (
        1 + psi_01 * frac[1] / frac[0]
    )
    condGap_pp = cond_pp[1] / (1 + phi_10 * frac[0] / frac[1]) + cond_pp[0] / (
        1 + phi_01 * frac[1] / frac[0]
    )
    condGap = condGap_p + condGap_pp
    return densGap, visGap, condGap, cpGap


# ========================================================================
# ========================================================================


def klemsPatch(sunAltitude, sunAzimuth, surfTilt, surfAzimuth):
    # ----------------------------------------------------
    # Calculates the Klems patch for a certain sun position, surface azimuth and facade tilt angle
    # Convention: Sun and surface azimuth is zero at south, positive towards west.
    # Surf tilt: 0 degrees= vertical, +90 degrees= horizontal
    # BB and ThP 2014
    # ----------------------------------------------------
    # The angles for the two rotation (around Z axis and then the new Y axis);
    gamma = -surfAzimuth
    beta = surfTilt
    # Position of the sun in the sun cartesian coordinate system (X=South, Y=West, Z=Upward), (3,1) matrix
    vectSun = matrix(
        [
            [cos(sunAzimuth) * cos(sunAltitude)],
            [-sin(sunAzimuth) * cos(sunAltitude)],
            [sin(sunAltitude)],
        ]
    )
    # rotation matrix on the z-axis,(3,3) matrix
    rotMtx1 = matrix(
        [[cos(gamma), sin(gamma), 0], [-sin(gamma), cos(gamma), 0], [0, 0, 1]]
    )
    # rotation matrix on the y'-axis, (3,3) matrix
    rotMtx2 = matrix([[cos(beta), 0, -sin(beta)], [0, 1, 0], [sin(beta), 0, cos(beta)]])
    # sun position in cartesian klems coordinates, (3,1) matrix
    vectSurf = rotMtx2 * (rotMtx1 * vectSun)
    # sun positions in spherical klems coordinates, (3,1) matrix
    thetaKlems = arccos(vectSurf[2, 0])  # arcos(Z)
    phiKlems = arctan2(vectSurf[0, 0], vectSurf[1, 0])  # arctan(X/Y)
    if phiKlems < 0:
        phiKlems = 2 * pi + phiKlems
    # patch selection
    patch = searchPatch(phiKlems * 180.0 / pi, thetaKlems * 180.0 / pi)
    # ----------------------------------------------------
    return patch


# ========================================================================
# ========================================================================


def winCHTC(tempInAir, intTemp, area, tilt):
    # convection heat transfer coefficients - internal side ISO15099
    # tilt angle of the window in degrees
    # tilt = 0 is horizontal glazing, heat flow upwards
    # tilt = 90 is vertical glazing, heat flow upwards
    # tilt = 180 is horizontal glazing, heat flow downwards
    minCHTC = 0.1
    if tempInAir < intTemp:
        tilt = 180 - tilt
    dens, vis, cond, cpGap, mMGap = calcs.gasProp(
        0, tempInAir + 0.25 * (intTemp - tempInAir)
    )
    height = sqrt(area)
    rayleigh = (
        dens ** 2
        * height ** 3
        * p.g
        * cpGap
        * abs(intTemp - tempInAir)
        / vis
        / cond
        / (tempInAir + 0.25 * (intTemp - tempInAir) + 273.15)
    )
    if tilt >= 0.0 and tilt < 15:
        nusselt = 0.13 * rayleigh ** (1.0 / 3.0)
    elif tilt >= 15 and tilt <= 90:
        rayleigh_cv = 2.5e5 * (exp(0.72 * tilt) / sin(tilt * pi / 180)) ** (1.0 / 5.0)
        if rayleigh_cv >= rayleigh:
            nusselt = 0.56 * (rayleigh * sin(tilt * pi / 180)) ** (1.0 / 4.0)
        else:
            nusselt = 0.13 * (
                rayleigh ** (1.0 / 3.0) - rayleigh_cv ** (1.0 / 3.0)
            ) + 0.56 * (rayleigh_cv * sin(tilt * pi / 180)) ** (1.0 / 4.0)
    elif tilt > 90 and tilt <= 179:
        nusselt = 0.56 * (rayleigh * sin(tilt * pi / 180)) ** (1.0 / 4.0)
    elif tilt > 179 and tilt <= 180:
        nusselt = 0.58 * (rayleigh * sin(tilt * pi / 180)) ** (1.0 / 5.0)
    # ----------------------------------------------------
    return max(nusselt * cond / height, minCHTC)


# ========================================================================
# ========================================================================


def surfCHTC(tempInAir, intTemp, tilt):
    # convection heat transfer coefficients !!! check tilt
    if tilt < 10:
        if intTemp > tempInAir:
            CHTC = max(0.76 * (abs(intTemp - tempInAir) ** (1.0 / 3.0)), 0.1)  # 0.948
        else:
            CHTC = max(1.51 * (abs(intTemp - tempInAir) ** (1.0 / 3.0)), 0.1)  # 4.040
    elif tilt > 170:
        if intTemp >= tempInAir:
            CHTC = max(1.51 * (abs(intTemp - tempInAir) ** (1.0 / 3.0)), 0.1)  # 4.040
        else:
            CHTC = max(0.76 * (abs(intTemp - tempInAir) ** (1.0 / 3.0)), 0.1)  # 0.948
    else:
        CHTC = max(1.31 * (abs(intTemp - tempInAir) ** (1.0 / 3.0)), 0.1)  # 3.076
    return CHTC


# def surfCHTC(tempInAir,intTemp,tilt):
#   # convection heat transfer coefficients !!! check tilt
#   if tilt<10:
#      if intTemp > tempInAir:
#         CHTC = 0.948
#      else:
#         CHTC = 4.040
#   elif tilt>170:
#      if intTemp >= tempInAir:
#         CHTC = 4.040
#      else:
#         CHTC = 0.948
#   else:
#      CHTC = 3.076
#   return CHTC

# ========================================================================
# ========================================================================


def transientConduction(
    dt, nLayer, layerTemp, volHeat, thermalCond, width, flx, bc, temp2
):
    # ----------------------------------------------------
    # implicit coefficient
    fimp = 0.5
    # explicit coefficient
    fexp = 0.5
    # mean thermal conductivity over distance between 2 layers
    tcp = zeros(nLayer + 1)
    # thermal capacity times layer depth
    hcp = zeros(nLayer)
    # lower, main, and upper diagonals
    za = zeros([nLayer, 3])
    # RHS
    zy = zeros(nLayer)
    # --------------------------------------------------------------------------
    hcp[0] = volHeat[0] * width[0]
    tcp[0] = 2.0 * thermalCond[0] / width[0]
    tcp[nLayer] = 2.0 * thermalCond[nLayer - 1] / width[nLayer - 1]
    for j in range(1, nLayer):
        tcp[j] = 2.0 / (width[j - 1] / thermalCond[j - 1] + width[j] / thermalCond[j])
        hcp[j] = volHeat[j] * width[j]
    # --------------------------------------------------------------------------
    za[0, 0] = 0.0
    za[0, 1] = hcp[0] / dt + fimp * tcp[1]
    za[0, 2] = -fimp * tcp[1]
    zy[0] = (
        hcp[0] / dt * layerTemp[0]
        - fexp * tcp[1] * (layerTemp[0] - layerTemp[1])
        + flx[0]
    )
    # --------------------------------------------------------------------------
    for j in range(1, nLayer - 1):
        za[j, 0] = fimp * (-tcp[j])
        za[j, 1] = hcp[j] / dt + fimp * (tcp[j] + tcp[j + 1])
        za[j, 2] = fimp * (-tcp[j + 1])
        zy[j] = (
            hcp[j] / dt * layerTemp[j]
            + fexp
            * (
                tcp[j] * layerTemp[j - 1]
                - tcp[j] * layerTemp[j]
                - tcp[j + 1] * layerTemp[j]
                + tcp[j + 1] * layerTemp[j + 1]
            )
            + flx[j]
        )
    # --------------------------------------------------------------------------
    if bc == 1:  # het flux
        za[nLayer - 1, 0] = fimp * (-tcp[nLayer - 1])
        za[nLayer - 1, 1] = hcp[nLayer - 1] / dt + fimp * tcp[nLayer - 1]
        za[nLayer - 1, 2] = 0.0
        zy[nLayer - 1] = (
            hcp[nLayer - 1] / dt * layerTemp[nLayer - 1]
            + fexp * tcp[nLayer - 1] * (layerTemp[nLayer - 2] - layerTemp[nLayer - 1])
            + flx[nLayer - 1]
        )
    elif bc == 0:  # deep-temperature
        za[nLayer - 1, 0] = 0
        za[nLayer - 1, 1] = 1
        za[nLayer - 1, 2] = 0.0
        zy[nLayer - 1] = temp2
    else:
        print("ERROR: check input parameters for TransientConduction routine")
    # --------------------------------------------------------------------------
    layerTemp = calcs.invert(nLayer, za, zy)
    # --------------------------------------------------------------------------
    extLayerTemp = flx[0] / tcp[0] + layerTemp[0]
    intLayerTemp = flx[nLayer - 1] / tcp[nLayer] + layerTemp[nLayer - 1]
    # --------------------------------------------------------------------------
    return layerTemp, extLayerTemp, intLayerTemp


# ========================================================================
# ========================================================================


def transientConductionWin(
    dt,
    nLayer,
    layerTemp,
    volHeat,
    thermalCond,
    width,
    intEmissWin,
    extEmissWin,
    widthGap,
    ventGap,
    fracGap,
    areaWin,
    flx,
    vertical,
    bc,
    temp2,
):
    # ----------------------------------------------------
    # implicit coefficient
    fimp = 0.5
    # explicit coefficient
    fexp = 0.5
    # mean thermal conductivity over distance between 2 layers
    tcp = zeros(nLayer)
    # thermal capacity times layer depth
    hcp = zeros(nLayer)
    # lower, main, and upper diagonals
    za = zeros([nLayer, 3])
    # RHS
    zy = zeros(nLayer)
    # resistance
    res = zeros(nLayer + 1)
    # --------------------------------------------------------------------------
    res[0] = width[0] / 2 / thermalCond[0]
    res[nLayer] = width[nLayer - 1] / 2 / thermalCond[nLayer - 1]
    hcp[0] = volHeat[0] * width[0]
    for j in range(1, nLayer):
        chtc = calcs.chtcGlazCav(
            widthGap[j - 1],
            ventGap[j - 1],
            fracGap[j - 1, :],
            layerTemp[j - 1],
            layerTemp[j],
            areaWin,
            vertical,
        )
        rhtc = (
            4
            * p.sigma
            * ((layerTemp[j - 1] + layerTemp[j]) / 2) ** 3
            / (1 / intEmissWin[j - 1] + 1 / extEmissWin[j] - 1)
        )
        res[j] = (
            width[j - 1] / 2 / thermalCond[j - 1]
            + width[j] / 2 / thermalCond[j]
            + 1 / (chtc + rhtc)
        )
        tcp[j] = 1 / res[j]
        hcp[j] = volHeat[j] * width[j]
    # --------------------------------------------------------------------------
    za[0, 0] = 0.0
    za[0, 1] = hcp[0] / dt + fimp * tcp[1]
    za[0, 2] = -fimp * tcp[1]
    zy[0] = (
        hcp[0] / dt * layerTemp[0]
        - fexp * tcp[1] * (layerTemp[0] - layerTemp[1])
        + flx[0]
    )
    # --------------------------------------------------------------------------
    for j in range(1, nLayer - 1):
        za[j, 0] = fimp * (-tcp[j])
        za[j, 1] = hcp[j] / dt + fimp * (tcp[j] + tcp[j + 1])
        za[j, 2] = fimp * (-tcp[j + 1])
        zy[j] = (
            hcp[j] / dt * layerTemp[j]
            + fexp
            * (
                tcp[j] * layerTemp[j - 1]
                - tcp[j] * layerTemp[j]
                - tcp[j + 1] * layerTemp[j]
                + tcp[j + 1] * layerTemp[j + 1]
            )
            + flx[j]
        )
    # --------------------------------------------------------------------------
    if bc == 1:  # het flux
        za[nLayer - 1, 0] = fimp * (-tcp[nLayer - 1])
        za[nLayer - 1, 1] = hcp[nLayer - 1] / dt + fimp * tcp[nLayer - 1]
        za[nLayer - 1, 2] = 0.0
        zy[nLayer - 1] = (
            hcp[nLayer - 1] / dt * layerTemp[nLayer - 1]
            + fexp * tcp[nLayer - 1] * (layerTemp[nLayer - 2] - layerTemp[nLayer - 1])
            + flx[nLayer - 1]
        )
    elif bc == 0:  # deep-temperature
        za[nLayer - 1, 0] = 0
        za[nLayer - 1, 1] = 1
        za[nLayer - 1, 2] = 0.0
        zy[nLayer - 1] = temp2
    else:
        print("ERROR: check input parameters for TransientConduction routine")
    # --------------------------------------------------------------------------
    layerTemp = calcs.invert(nLayer, za, zy)
    # --------------------------------------------------------------------------
    extLayerTemp = layerTemp[0]
    intLayerTemp = layerTemp[nLayer - 1]
    # --------------------------------------------------------------------------
    return layerTemp, extLayerTemp, intLayerTemp


# ========================================================================
# ========================================================================


def invert(nz, a, c):
    #  a(*,1) lower diagonal (Ai,i-1)
    #  a(*,2) principal diagonal (Ai,i)
    #  a(*,3) upper diagonal (Ai,i+1)
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


# ========================================================================
# ========================================================================


def incAngle(sunAltitude, sunAzimuth, surfAzimuth, vertical):
    # angles in radians
    if vertical:
        incAng = arccos(cos(sunAltitude) * cos(sunAzimuth - surfAzimuth))
        # if incAng > pi/2. or incAng < -pi/2.:
        #   incAng = pi/2.
    else:
        incAng = pi / 2.0 - sunAltitude
    # ----------------------------------------------------------------------
    return incAng


# ========================================================================
# ========================================================================


def profAngle(sunAltitude, sunAzimuth, surfAzimuth, slatHor):
    # only valid for vertical surfaces
    # angles in radians
    # if the slats are horizontal
    if slatHor:
        profAng = arctan(tan(sunAltitude) / cos(sunAzimuth - surfAzimuth))
    # if the slats are vertical
    else:
        #      profAng = arctan(tan(sunAzimuth-surfAzimuth)/cos(sunAltitude))
        profAng = sunAzimuth - surfAzimuth
        # profAng = -sunAzimuth+surfAzimuth
    # ----------------------------------------------------------------------
    return profAng


# ========================================================================
# ========================================================================


def klemsAng():
    # ----------------------------------------------------
    # calculate the theta and phi angles of each patch (degrees)
    # This function is for incident hemisphere.
    # Theta is the polar angle. Phi is the azimuth angle.
    # Convention: theta =0 when the point is at the summit of the sphere, and theta=90 when
    # the point is on the surface
    # Convention: if the surface would be perfectly vertical and facing south, phi=0 would
    # be facing east, phi=90 facing upwards, phi=180 facing west, phi=270 facing downward.
    polarAngle = [
        0,
        5.0,
        15.0,
        25.0,
        35.0,
        45.0,
        55.0,
        65.0,
        75.0,
        90.0,
    ]  # polar angles of the boundaries [deg]
    azimuthDivision = [1, 8, 16, 20, 24, 24, 24, 16, 12]  # azimuth divisions
    patch = 0
    phi = zeros(145)
    theta = zeros(145)
    sumAzimuthDivision = 1
    for i in range(1, len(polarAngle) - 1):
        azimuthStep = 360 / azimuthDivision[i]
        for count in range(azimuthDivision[i]):
            patch = sumAzimuthDivision + count
            theta[patch] = polarAngle[i] / 2.0 + polarAngle[i + 1] / 2.0
            phi[patch] = 180 + azimuthStep * count
            if phi[patch] >= 360:
                phi[patch] = phi[patch] - 360
        sumAzimuthDivision = sumAzimuthDivision + azimuthDivision[i]
    return theta, phi


# ========================================================================
# ========================================================================


def klemsProf():
    # ----------------------------------------------------
    # calculate the profile angle of each patch (radians)
    theta, phi = klemsAng()
    profAng = zeros(len(theta))
    for patch in range(1, len(theta)):
        if phi[patch] > 270:
            phiRad = radians(phi[patch] - 360)
        elif phi[patch] > 90:
            phiRad = radians(180 - phi[patch])
        else:
            phiRad = radians(phi[patch])
        profAng[patch] = arctan(
            tan(radians(theta[patch]))
            * tan(phiRad)
            / (tan(radians(theta[patch])) ** 2 + tan(phiRad) ** 2) ** 0.5
        )
    #       print(patch,theta[patch],phi[patch],degrees(profAng[patch]))
    return profAng


# ========================================================================
# ========================================================================


def diffPropGlass(angles, prop):
    # ------------------------------------------------------------------------
    a = [0.0170, 0.1294, 0.2241, 0.2588, 0.2241, 0.1294, 0.0170]
    diffProp = a[0] * prop[0]
    for i in range(1, 7):
        diffProp = diffProp + a[i] * (interp(15 * i, angles, prop))
    # ------------------------------------------------------------------------
    return diffProp


# ========================================================================
# ========================================================================


def diffPropBlinds(angles, prop):
    # ------------------------------------------------------------------------
    b = [0.1304, 0.1261, 0.1130, 0.0923, 0.0653, 0.0338, 0.0043]
    diffProp = b[0] * interp(0, angles, prop)
    for i in range(1, 7):
        diffProp = diffProp + b[i] * (
            interp(15 * i, angles, prop) + interp(-15 * i, angles, prop)
        )
    # ------------------------------------------------------------------------
    return diffProp


# ========================================================================
# ========================================================================


def kuhn2006_intBlinds(
    UGlass,
    gGlass,
    tauGlass,
    tauGlassDif,
    rhoGlass,
    rhoGlassDif,
    tauBlind,
    tauBlindDif,
    rhoBlind,
    rhoBlindDif,
    kapp,
):
    tauTot = zeros(145)
    rhoTot = zeros(145)
    gTot_1 = zeros(145)
    gTot_2 = zeros(145)
    absBlind = zeros(145)
    del2 = 18  # W m-2 K-1
    delInt = 1.0 / (1.0 / UGlass + 1.0 / del2)
    # -----------------------------------------------------
    absGlassDif = 1 - tauGlassDif - rhoGlassDif
    for i in range(145):
        absBlind[i] = 1 - tauBlind[i] - rhoBlind[i]
        # -----------------------------------------------------
        tauTot[i] = tauGlass[i] * tauBlind[i] + tauGlass[i] * rhoBlind[
            i
        ] * rhoGlassDif * tauBlindDif / (1 - rhoBlindDif * rhoGlassDif)
        # -----------------------------------------------------
        gTot_2[i] = (
            gGlass[i]
            - tauGlass[i] * tauGlassDif * rhoBlind[i] / (1 - rhoBlindDif * rhoGlassDif)
            - kapp
            * tauGlass[i]
            * absGlassDif
            * rhoBlind[i]
            / (1 - rhoBlindDif * rhoGlassDif)
            - gGlass[i] * absBlind[i] * delInt / del2
        )
        # -----------------------------------------------------
        rhoTot[i] = rhoGlass[i] + tauGlass[i] * tauGlassDif * rhoBlind[i] / (
            1 - rhoBlindDif * rhoGlassDif
        )
        # -----------------------------------------------------
        gTot_1[i] = (
            tauTot[i]
            + gGlass[i]
            - tauGlass[i]
            + (1 - kapp)
            * tauGlass[i]
            * absGlassDif
            * rhoBlind[i]
            / (1 - rhoBlindDif * rhoGlassDif)
            + tauGlass[i]
            * absBlind[i]
            * (1 - delInt / del2)
            / (1 - rhoBlindDif * rhoGlassDif)
        )

    ##--------------extras---------------------------------

    # ~ tauetot=tauTot[i]
    # ~ qigzg= gGlass[i] - tauGlass[i]
    # ~ reflected_absorbed=(1-kapp)*tauGlass[i]*absGlassDif*rhoBlind[i]/(1-rhoBlindDif*rhoGlassDif)
    # ~ extra_gains=tauGlass[i]*absBlind[i]*(1-delInt/del2)/(1-rhoBlindDif*rhoGlassDif)
    # ~
    # ~ back_exterior=(tauGlass[i]*tauGlassDif*rhoBlind[i]/(1-rhoBlindDif*rhoGlassDif))
    # ~ reflected_absorbed=(kapp*tauGlass[i]*absGlassDif*rhoBlind[i]/(1-rhoBlindDif*rhoGlassDif))
    # ~ extra_losses= (gGlass[i]*absBlind[i]*delInt/del2)
    # ~ a=gTot_2[i]-gTot_1[i]
    # ~
    # ~ print('difference',a)
    # ~ if i==0:
    # ~ print('percentages',gTot_1[i],tauetot,qigzg,reflected_absorbed,extra_gains)
    # ~ temp=[gTot_1[i],tauetot,qigzg,reflected_absorbed,extra_gains]
    # ~ print('percentages',temp)
    # ~ savetxt("/net/p/600/groupdrives/TSB/Aktuell/620_eeb/621_sbe/Aktuell/Berechnungsprogramme/Fener/fener180419/config/camilo/thesis/estadisticas/temp.out", temp, fmt='%1.6f')

    ## -------------------------------------------
    return tauTot, rhoTot, gTot_1


# ========================================================================
# ========================================================================


def kuhn2006_extBlinds(
    UGlass, gGlass, tauGlass, rhoGlassDif, tauBlind, rhoBlind, rhoBlindDif
):
    tauTot = zeros(145)
    # rhoTot = zeros(145)
    gTot = zeros(145)
    absBlind = zeros(145)
    del1 = 6.0  # W m-2 K-1
    del2 = 18.0  # W m-2 K-1
    delExt = 1.0 / (1.0 / UGlass + 1.0 / del1 + 1.0 / del2)
    # -----------------------------------------------------
    for i in range(145):
        absBlind[i] = 1 - tauBlind[i] - rhoBlind[i]
        # -----------------------------------------------------
        gTot[i] = (
            tauBlind[i] * gGlass[i]
            + absBlind[i] * delExt / del2
            + tauBlind[i] * (1 - gGlass[i]) * delExt / del1
        )
        # -----------------------------------------------------
        tauTot[i] = tauGlass[i] * tauBlind[i] / (1 - rhoBlindDif * rhoGlassDif)
    # -------------------------------------------
    return tauTot, gTot


# ========================================================================
# ========================================================================


def propKlemsProfAng(angles, prop):
    # ----------------------------------------------------
    # g-value interpolation in klems angles
    propKlems = []
    anglesRad = radians(angles)
    # calculate profile angle of each Klems patch
    profAngRad = calcs.klemsProf()
    # interpolate g-values
    for i in range(len(profAngRad)):
        flag = False
        for j in range(1, len(anglesRad)):
            if (
                profAngRad[i] >= anglesRad[j - 1]
                and profAngRad[i] < anglesRad[j]
                and not flag
            ):
                anglesRadRange = [anglesRad[j - 1], anglesRad[j]]
                anglesPropRange = [prop[j - 1], prop[j]]
                propKlems.append(interp(profAngRad[i], anglesRadRange, anglesPropRange))
                flag = True
        if not flag:
            print("ERROR: check routine")
    return propKlems


# ========================================================================
# ========================================================================


def propKlemsIncAng(angles, prop):
    # Klems model initialization
    propagationMtx, polarAngleM, azimuthDivision = calcs.klemsIni()
    # ----------------------------------------------------
    # g-value interpolation in klems angles
    propKlems = []
    anglesRad = radians(angles)
    g = interp(polarAngleM, anglesRad, prop, left=prop[0], right=0.0)
    k = 0
    for i in range(len(polarAngleM)):
        for j in range(azimuthDivision[i]):
            if k < 144:
                propKlems.append(g[i])
            else:
                propKlems.append(g[i])
            k = k + 1
    return propKlems


# ========================================================================
# ========================================================================


def kuhn11(tmxs, rmxs, calorim):
    # Generation of a virtual two-layer virtual system based on Kuhn et al 2011
    resExt = 0.13
    resInt = 0.04
    solidAngleMtx = calcs.klemsSolidAngles()
    res = 1 / calorim[0] - resExt - resInt
    absLayInt = zeros(145)
    absLayExt = zeros(145)
    transTot = zeros(145)
    reflTot = zeros(145)
    for i in range(145):
        gValue = calorim[i + 1]
        transTot[i] = sum(tmxs[i, :])
        reflTot[i] = sum(rmxs[i, :])
        absTot = 1 - transTot[i] - reflTot[i]
        solHeatIn = gValue - transTot[i]
        absLayInt[i] = (solHeatIn * (resExt + res + resInt) - absTot * resExt) / res
        absLayExt[i] = (
            -solHeatIn * (resExt + res + resInt) + absTot * (resExt + res)
        ) / res
        # weight by the solid angle and apply irradiance/radiance factor pi
        absLayInt[i] = absLayInt[i] * solidAngleMtx[i, i] / (2 * pi) * pi
        absLayExt[i] = absLayExt[i] * solidAngleMtx[i, i] / (2 * pi) * pi
    return res, absLayInt, absLayExt


# ========================================================================
# ========================================================================


def vectOrder(vect):
    # organizes a vector by the order of its values

    vectOrd = vect
    n = len(vect)
    indexOrd = arange(n)
    for i in range(n):
        for j in range(n):
            if vectOrd[i] < vectOrd[j]:
                aux = vectOrd[j]
                aux2 = indexOrd[j]
                vectOrd[j] = vectOrd[i]
                indexOrd[j] = indexOrd[i]
                vectOrd[i] = aux
                indexOrd[i] = aux2
    return vectOrd, indexOrd


# ========================================================================
# ========================================================================


def buildBSDFfromTauTot(fileName, tauTotV, rhoTotV, tauTotS, rhoTotS):
    # -----------------------------------------------------------------------------------
    # Build a BSDF from tauTot
    # -----------------------------------------------------------------------------------
    transMtxV = matrix(zeros([145, 145]))
    reflFrontMtxV = matrix(zeros([145, 145]))
    reflBackMtxV = matrix(zeros([145, 145]))
    transMtxS = matrix(zeros([145, 145]))
    reflFrontMtxS = matrix(zeros([145, 145]))
    reflBackMtxS = matrix(zeros([145, 145]))
    tmxs = matrix(zeros([145, 145]))
    rmxs = matrix(zeros([145, 145]))
    propagationMtx = calcs.klemsIni()[0]
    for i in range(145):
        transMtxV[:, i] = tauTotV[i] / propagationMtx[i, i] / 145
        reflFrontMtxV[:, i] = rhoTotV[i] / propagationMtx[i, i] / 145
        reflBackMtxV[:, i] = rhoTotV[i] / propagationMtx[i, i] / 145
        transMtxS[:, i] = tauTotS[i] / propagationMtx[i, i] / 145
        reflFrontMtxS[:, i] = rhoTotS[i] / propagationMtx[i, i] / 145
        reflBackMtxS[:, i] = rhoTotS[i] / propagationMtx[i, i] / 145
        tmxs[:, i] = tauTotS[i] / 145
        rmxs[:, i] = rhoTotS[i] / 145

    # write bsdf file
    ioRoutines.writeXml(
        "%s" % fileName,
        transMtxV,
        reflFrontMtxV,
        transMtxV,
        reflBackMtxV,
        transMtxS,
        reflFrontMtxS,
        transMtxS,
        reflBackMtxS,
    )

    return (
        transMtxV,
        reflFrontMtxV,
        reflBackMtxV,
        transMtxS,
        reflFrontMtxS,
        reflBackMtxS,
    )


# ========================================================================
# ========================================================================


# ! THIS ROUTINE MUST BE REVIEWED
def AzE_from_phitheta(phi, theta):
    # calculates azimuth angle (Az) and altitude angle (E) [according to EN14500 convention]
    # from phi and theta angles of Klems patches
    # Theta is the polar angle. Phi is the azimuth angle.
    # Convention: theta =0 when the point is at the summit of the sphere, and theta=90 when
    # the point is on the surface
    # Convention: if the surface would be perfectly vertical and facing south, phi=0 would
    # be facing east, phi=90 facing upwards, phi=180 facing west, phi=270 facing downward.
    # Units: radiants
    if phi == 0 and theta == 0:
        return 0, 0
    tan_Az = -cos(phi) * tan(
        theta
    )  # doesn't work if cos phi =0 and tan theta is +- inf
    sin_E = sin(phi) * sin(theta)  # doesn't work if sin phi =0 and tan theta is +- inf
    return atan(tan_Az), asin(sin_E)


def klems_spherical_distance_interpolation_Azimuth_Elevation(Az, E, values):
    # calculate values at klems patches from values at pairs of
    # elevation(E) and azimuth angles (Az) [according to EN14500 convention]
    # through 2d interpolation
    # Units: degrees
    theta_list, phi_list = calcs.klemsAng()
    klems = zeros(145)
    f = interp2d(E, Az, values, "quintic")
    for i in range(145):
        Az_i, E_i = AzE_from_phitheta(phi_list[i] * pi / 180, theta_list[i] * pi / 180)
        klems[i] = f(E_i * 180 / pi, Az_i * 180 / pi)
    return klems


def findSpecPatch(patch):
    theta_list, phi_list = klemsAng()
    theta = theta_list[patch]
    phi = phi_list[patch] + 180.0
    if phi >= 360.0:
        phi = phi - 360.0
    specPatch = searchPatch(phi, theta)
    return specPatch


def findAzimuthPatches(patch):
    # extracts those klems patches with similar azimuth (+-7.5 deg) to
    # the azimuth of the input 'patch'. It also extracts the theta and phi
    # angles of those patches. theta is negative if the patch is on the
    # other 'azimuth side' of input patch.
    patches = []
    theta_list, phi_list = klemsAng()
    phi_s = phi_list[patch]
    for i in range(145):
        phi = phi_list[i]
        if (
            (phi <= phi_s + 7.5 and phi >= phi_s - 7.5)
            or (phi <= phi_s + 7.5 + 180 and phi >= phi_s - 7.5 + 180)
            or (phi <= phi_s + 7.5 - 180 and phi >= phi_s - 7.5 - 180)
        ):
            patches.append(i)
    # addition of the first patch
    flag = True
    for j in patches:
        if j == 0:
            flag = False
    if flag:
        patches.append(0)
    # theta and phi lists of selected patches.
    theta_patches = zeros(len(patches))
    phi_patches = zeros(len(patches))
    for i in range(len(patches)):
        theta_patches[i] = theta_list[patches[i]]
        phi_patches[i] = phi_list[patches[i]]
        if (
            (phi_s < 180 and phi_patches[i] > 180)
            or (phi_patches[i] < 180 and phi_s > 180)
            or (phi_s == 0 and phi_patches[i] == 180)
            or (phi_patches[i] == 0 and phi_s == 180)
        ):
            theta_patches[i] = -theta_patches[i]
    return patches, theta_patches, phi_patches


# ~ def klemsLinearDistribution(specPatch,value,theta,theta_max,theta_min):
# ~ theta_list,phi_list= klemsAng()
# ~ theta_s = theta_list[specPatch]
# ~ interpValue = 0
# ~ if theta > theta_s:
# ~ interpValue = interp(theta,[theta_s,theta_max],[value,0])
# ~ else:
# ~ interpValue = interp(theta,[theta_min,theta_s],[0,value])
# ~ return interpValue


def phitheta_from_AzE(Az, E):
    # this routine has been developed and tested in the context of the arkol project 24.07.2018
    # calculates Klems phi and theta angles from
    # azimuth angle (Az) and altitude angle (E) [according to EN14500 convention]
    # Valid only for a south-facing facade
    # Klems convention: Theta is the polar angle. Phi is the azimuth angle.
    # theta =0 when the point is at the summit of the sphere, and theta=90 when
    # the point is on the surface
    # Convention: if the surface would be perfectly vertical and facing south, phi=0 would
    # be facing west, phi=90 facing downwards, phi=180 facing east, phi=270 facing upwards.
    # Units: radiants
    if abs(Az) < 0.0001:
        theta = abs(E)
        if E >= 0.0001:
            phi = radians(270)
        else:
            phi = radians(90)

    else:
        theta = acos(cos(abs(Az)) * cos(abs(E)))
        if E >= 0:
            if Az >= 0.0001:
                phi = -atan(tan(E) / sin(Az))
            else:
                phi = atan(tan(E) / sin(abs(Az))) + pi
        else:
            if Az >= 0.0001:
                phi = atan(tan(abs(E)) / sin(Az))
            else:
                phi = -atan(tan(abs(E)) / sin(abs(Az))) + pi
        if phi < 0:
            phi += 2 * pi
    return phi, theta


def searchPatch(phi, theta):
    # this routine has been developed and tested in the context of the arkol project 24.07.2018
    # Delivers the Klems patch when provided with the phi and theta angles (degrees)
    # This function is for incident hemisphere.
    # Theta is the polar angle. Phi is the azimuth angle.
    # Convention: theta =0 when the point is at the summit of the sphere, and theta=90 when
    # the point is on the surface
    # Convention: if the surface would be perfectly vertical and facing south, phi=0 would
    # be facing west, phi=90 facing downwards, phi=180 facing east, phi=270 facing upwards.
    # ----------------------------------------------------
    if phi > 360 or phi < 0 or theta > 90 or theta < 0:
        print("ERROR: WRONG ENTRY (theta,phi)", theta, phi)
    patch = -1
    polarAngle = [
        0,
        5.0,
        15.0,
        25.0,
        35.0,
        45.0,
        55.0,
        65.0,
        75.0,
        90.0,
    ]  # polar angles of the boundaries [deg]
    azimuthDivision = [1, 8, 16, 20, 24, 24, 24, 16, 12]  # azimuth divisions
    sumAzimuthDivision = 1
    if theta >= polarAngle[0] and theta < polarAngle[1]:
        patch = 0
    for i in range(2, len(polarAngle)):
        if theta >= polarAngle[i - 1] and theta < polarAngle[i]:
            azimuthStep = 360 / azimuthDivision[i - 1]
            for count in range(azimuthDivision[i - 1]):
                lowLimit = azimuthStep / 2.0 + azimuthStep * (count - 1)
                highLimit = azimuthStep / 2.0 + azimuthStep * count
                if lowLimit >= 360:
                    lowLimit = lowLimit - 360.0
                if highLimit >= 360:
                    highLimit = highLimit - 360.0
                if (
                    (phi >= lowLimit and phi < highLimit)
                    or (highLimit < lowLimit and phi < highLimit)
                    or (highLimit < lowLimit and phi >= lowLimit)
                ):
                    patch = sumAzimuthDivision + count
        sumAzimuthDivision = sumAzimuthDivision + azimuthDivision[i - 1]
    # ----------------------------------------------------
    return patch
