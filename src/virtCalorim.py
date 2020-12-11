#!/usr/bin/env python
# -*- coding: utf-8 -*-


from numpy import zeros, transpose, ones
import iso
import ioRoutines
import calcs


def unitWinAbs(numConWin, numPaneConWin, conWin, patchIndex):
    #  window layer solar absorption
    unitRadVect = zeros(145)
    unitRadVect[patchIndex] = 1
    irrWinLay = zeros([numConWin, int(max(numPaneConWin))])
    for i in range(numConWin):
        for k in range(conWin[i].numPane):
            irrWinLay[i, k] = (
                transpose(conWin[i].absFront[k, 0, :]).dot(unitRadVect).reshape(-1)
            )
    return irrWinLay


def gValue(
    numConWin,
    numFrame,
    numSurf,
    numPaneConWin,
    irrWinLay,
    win,
    frame,
    surf,
    tmxs,
    patchIndex,
):
    # window heat transfer according to ISO15099
    # boundary conditions
    tempInAir_u = 21 + 273.15
    tempOutAir_u = -18.0 + 273.15
    tempInAir_s = 24 + 273.15
    tempOutAir_s = 32.0 + 273.15  # 32.0+273.15
    extCHTC = 25.0  # 25.
    intCHTC = 7.7 * ones(numConWin)
    irrad = 680.0  # 680
    # with solar radiation
    irrWinLayAux = zeros([1, int(max(numPaneConWin))])
    qInt_w = zeros(numConWin)
    for i in range(numConWin):
        for k in range(int(max(numPaneConWin))):
            irrWinLayAux[0, k] = irrWinLay[i, k] * irrad
        win_w = [win[0]]
        win_w = iso.compute(
            win_w,
            surf,
            frame,
            tempOutAir_s,
            tempInAir_s,
            extCHTC,
            intCHTC,
            0.0,
            0.0,
            0.0,
            0.0,
            irrWinLayAux,
            dynamicFlag=False,
        )
        qInt_w[i] = -win_w[
            0
        ].qInt  # win.qInt refers to heat flux from indoor to surface
    # without solar radiation
    #    U-value boundary conditions
    uValue = zeros(numConWin)
    for i in range(numConWin):
        win_wo_u = [win[0]]
        win_wo_u = iso.compute(
            win_wo_u,
            surf,
            frame,
            tempOutAir_u,
            tempInAir_u,
            extCHTC,
            intCHTC,
            0.0,
            0.0,
            0.0,
            0.0,
            zeros([1, int(max(numPaneConWin))]),
            dynamicFlag=False,
        )
        uValue[i] = win_wo_u[0].uValue
    #    solar boundary conditions
    qInt_wo_s = zeros(numConWin)
    for i in range(numConWin):
        win_wo_s = [win[0]]
        win_wo_s = iso.compute(
            win_wo_s,
            surf,
            frame,
            tempOutAir_s,
            tempInAir_s,
            extCHTC,
            intCHTC,
            0.0,
            0.0,
            0.0,
            0.0,
            zeros([1, int(max(numPaneConWin))]),
            dynamicFlag=False,
        )
        qInt_wo_s[i] = -win_wo_s[0].qInt
    # calculate angular g-value
    calorim = zeros(numConWin)
    for i in range(numConWin):
        calorim[i] = sum(tmxs[i, :, patchIndex]) + (qInt_w[i] - qInt_wo_s[i]) / irrad
    # -----------------------------------------------------
    return calorim, uValue


def main(
    numConWin,
    numFrame,
    numSurf,
    numPaneConWin,
    conWin,
    win,
    frame,
    surf,
    tmxs,
    calorim_path,
):
    calorim = zeros([numConWin, 146])
    calorim[:, 0] = gValue(
        numConWin,
        numFrame,
        numSurf,
        numPaneConWin,
        zeros([int(numConWin), int(max(numPaneConWin))]),
        win,
        frame,
        surf,
        tmxs,
        0,
    )[1]
    for patchIndex in range(145):
        irrWinLay = unitWinAbs(numConWin, numPaneConWin, conWin, patchIndex)
        calorim[:, patchIndex + 1] = gValue(
            numConWin,
            numFrame,
            numSurf,
            numPaneConWin,
            irrWinLay,
            win,
            frame,
            surf,
            tmxs,
            patchIndex,
        )[0]
    ioRoutines.writeDat(calorim_path, numConWin, calorim)
    print("Calorimetric file(s) saved in: %s" % calorim_path)
    print("----------------")
    for i in range(numConWin):
        print("U-value", calorim[i, 0])
        print(
            "g(0)",
            calorim[i, calcs.searchPatch(270, 0) + 1],
            "t(0)",
            sum(tmxs[i, :, calcs.searchPatch(270, 0)]),
        )
        print(
            "g(60)",
            calorim[i, calcs.searchPatch(270, 60) + 1],
            "t(60)",
            sum(tmxs[i, :, calcs.searchPatch(270, 60)]),
        )
        print("----------------")
