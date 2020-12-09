#!/usr/bin/env python
# -*- coding: utf-8 -*-


from numpy import zeros, savetxt
import ioRoutines
import calcs
from variables import conWinWithoutGap
import param as p

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


# transmission matrix, BSDF
def main(
    conWinVect,
    matGlz,
    bsdfSys_path,
    absFront_path,
    absBack_path,
    numBsdfLay=0,
    bsdfLay_file=0,
):
    for j in range(conWinVect.shape[0]):
        conWin = conWinWithoutGap(conWinVect[j, :])
        layerTransMtxV = zeros([int(conWin[0]), p.nKlem, p.nKlem])
        layerReflFrontMtxV = zeros([int(conWin[0]), p.nKlem, p.nKlem])
        layerReflBackMtxV = zeros([int(conWin[0]), p.nKlem, p.nKlem])
        layerTransMtxS = zeros([int(conWin[0]), p.nKlem, p.nKlem])
        layerReflFrontMtxS = zeros([int(conWin[0]), p.nKlem, p.nKlem])
        layerReflBackMtxS = zeros([int(conWin[0]), p.nKlem, p.nKlem])
        for k in range(int(conWin[0])):
            if conWin[2 * k + 1] == 1:
                (
                    layerTransMtxV[k, :, :],
                    layerReflFrontMtxV[k, :, :],
                    layerReflBackMtxV[k, :, :],
                ) = calcs.biDirGlazing(
                    matGlz[int(conWin[2 * k + 2]), 4],
                    matGlz[int(conWin[2 * k + 2]), 5],
                    matGlz[int(conWin[2 * k + 2]), 6],
                    matGlz[int(conWin[2 * k + 2]), 11],
                )
                (
                    layerTransMtxS[k, :, :],
                    layerReflFrontMtxS[k, :, :],
                    layerReflBackMtxS[k, :, :],
                ) = calcs.biDirGlazing(
                    matGlz[int(conWin[2 * k + 2]), 1],
                    matGlz[int(conWin[2 * k + 2]), 2],
                    matGlz[int(conWin[2 * k + 2]), 3],
                    matGlz[int(conWin[2 * k + 2]), 11],
                )
            elif conWin[2 * k + 1] == 4:
                if numBsdfLay == 0:
                    print(
                        "ERROR: Your window contruction contains BSDF layers. Check the input variables numBsdfLay and bsdfLay_file"
                    )
                layerTransMtxV[k, :, :] = ioRoutines.xml(
                    "%s" % (bsdfLay_file[int(conWin[2 * k + 2])]),
                    "Visible",
                    "Transmission Front",
                )
                layerReflFrontMtxV[k, :, :] = ioRoutines.xml(
                    "%s" % (bsdfLay_file[int(conWin[2 * k + 2])]),
                    "Visible",
                    "Reflection Front",
                )
                layerReflBackMtxV[k, :, :] = ioRoutines.xml(
                    "%s" % (bsdfLay_file[int(conWin[2 * k + 2])]),
                    "Visible",
                    "Reflection Back",
                )
                layerTransMtxS[k, :, :] = ioRoutines.xml(
                    "%s" % (bsdfLay_file[int(conWin[2 * k + 2])]),
                    "Solar",
                    "Transmission Front",
                )
                layerReflFrontMtxS[k, :, :] = ioRoutines.xml(
                    "%s" % (bsdfLay_file[int(conWin[2 * k + 2])]),
                    "Solar",
                    "Reflection Front",
                )
                layerReflBackMtxS[k, :, :] = ioRoutines.xml(
                    "%s" % (bsdfLay_file[int(conWin[2 * k + 2])]),
                    "Solar",
                    "Reflection Back",
                )
        transMtxV, reflFrontMtxV, reflBackMtxV, absFrontV, absBackV = calcs.genBSDF(
            conWin[0], layerTransMtxV, layerReflFrontMtxV, layerReflBackMtxV
        )
        transMtxS, reflFrontMtxS, reflBackMtxS, absFront, absBack = calcs.genBSDF(
            conWin[0], layerTransMtxS, layerReflFrontMtxS, layerReflBackMtxS
        )
        ioRoutines.writeXml(
            "%s" % (bsdfSys_path[j]),
            transMtxV,
            reflFrontMtxV,
            transMtxV,
            reflBackMtxV,
            transMtxS,
            reflFrontMtxS,
            transMtxS,
            reflBackMtxS,
        )
        for k in range(int(conWin[0])):
            savetxt("%s" % (absFront_path[j, k]), absFront[k], delimiter=",")
            savetxt("%s" % (absBack_path[j, k]), absBack[k], delimiter=",")
        print("BSDF file(s) saved in: %s" % bsdfSys_path)
        print(
            "Layer absortivities file(s) saved in: %s and %s"
            % (absFront_path, absBack_path)
        )
