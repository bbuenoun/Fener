#!/usr/bin/env python
# -*- coding: utf-8 -*-


from numpy import zeros, ones, cos, sin, math, atleast_2d, genfromtxt, arccos, pi

# import sys
import ftrans
import param


def convertWidthLayer(widthReal):
    # ----------------------------------------------------
    # widthReal – thickness of a layer in reality (e.g. brick 0.3)
    # widthModel – thickness of a layer in the model (e.g.widthModel_1=0.05,widthModel_2=0.1,etc)
    widthModel = []
    if widthReal <= 0.05:
        widthModel.append(widthReal)
    elif widthReal <= 0.1:
        widthModel.append(widthReal / 2)
        widthModel.append(widthReal / 2)
    else:
        for n in range(2, 10):
            if widthReal > 0.1 * (n - 1) and widthReal <= 0.1 * n:
                widthModel.append(widthReal / (2 * n))
                for i in range(1, n):
                    widthModel.append(widthReal / n)
                widthModel.append(widthReal / (2 * n))
    # ----------------------------------------------------
    return widthModel


class surf(object):
    def __init__(
        self,
        vect,
        frameVect,
        surfID,
        thermFlag=False,
        const=0,
        mat=0,
        iniTemp=0,
        finDiffFlag=False,
        simTimeStep=0,
    ):
        self.length = vect[0]
        self.height = vect[1]
        self.thick = vect[2]
        self.tx = vect[3]
        self.ty = vect[4]
        self.tz = vect[5]
        self.rx = vect[6]
        self.ry = vect[7]
        self.rz = vect[8]
        self.bc = int(vect[9])
        # ---------------------------------------------
        self.tilt = arccos(-sin(self.rx * pi / 180)) * 180 / pi
        self.area = vect[0] * vect[1]
        self.areaRev = 0
        for i in range(frameVect.shape[0]):
            if frameVect[i, 3] == surfID:
                self.area = self.area - frameVect[i, 0] * frameVect[i, 1]
                self.areaRev = (
                    self.areaRev
                    + frameVect[i, 0] * (self.thick - frameVect[i, 6])
                    + frameVect[i, 1] * (self.thick - frameVect[i, 6])
                )
        if self.area < 0.0:
            raise TypeError("check geometry")
        # ---------------------------------------------
        self.intAlb = vect[12]
        self.extAlb = vect[11]
        if thermFlag:
            self.svf = vect[10]
            self.extEmis = 0.9
            self.intEmis = 0.9
            self.extTemp = iniTemp
            self.intTemp = iniTemp
            # ---------------------------------------------
            numLayerReal = int(const[int(vect[13]), 0])
            layers = 0
            for i in range(numLayerReal):
                aux = convertWidthLayer(
                    mat[int(const[int(vect[13]), numLayerReal - i]), 0]
                )
                layers = layers + len(aux)
            self.numLayer = layers
            thickLay = zeros(self.numLayer)
            conduct = zeros(self.numLayer)
            volHeat = zeros(self.numLayer)
            layers = 0
            for i in range(numLayerReal):
                aux = convertWidthLayer(
                    mat[int(const[int(vect[13]), numLayerReal - i]), 0]
                )
                for j in range(len(aux)):
                    thickLay[j + layers] = aux[j]
                    conduct[j + layers] = mat[
                        int(const[int(vect[13]), numLayerReal - i]), 1
                    ]
                    volHeat[j + layers] = (
                        mat[int(const[int(vect[13]), numLayerReal - i]), 3]
                        * mat[int(const[int(vect[13]), numLayerReal - i]), 2]
                    )
                layers = layers + len(aux)
            # ---------------------------------------------
            if finDiffFlag:
                self.res = zeros(self.numLayer)
                self.cap = zeros(self.numLayer + 1)
                self.temp = ones(self.numLayer + 1) * iniTemp
                self.cap[0] = thickLay[0] / 2.0 * volHeat[0]
                self.res[0] = thickLay[0] / conduct[0]
                for i in range(1, self.numLayer):
                    self.cap[i] = (
                        thickLay[i - 1] / 2.0 * volHeat[i - 1]
                        + thickLay[i] / 2.0 * volHeat[i]
                    )
                    self.res[i] = thickLay[i] / conduct[i]
                self.cap[self.numLayer] = (
                    thickLay[self.numLayer - 1] / 2.0 * volHeat[self.numLayer - 1]
                )
            # ---------------------------------------------
            # transfer functions
            if simTimeStep > 0:
                timeStep = simTimeStep
            else:
                timeStep = param.timeStep
            (
                self.u,
                self.a,
                self.c,
                self.b,
                self.d,
                self.nc,
                self.nd,
                self.ns,
            ) = ftrans.ztran(
                layers, thickLay[::-1], conduct[::-1], volHeat[::-1], timeStep
            )
            self.qse = zeros(self.ns)
            self.qsi = zeros(self.ns)
            self.tse = ones(self.ns) * iniTemp
            self.tsi = ones(self.ns) * iniTemp

    def genbox(self, fw, surfFile):
        fw.write(
            "!genbox surf_mat %s %1.3f %1.3f %1.3f |xform -rx %1.3f -t %1.3f %1.3f %1.3f -rz %1.3f\n"
            % (
                surfFile,
                self.length,
                self.thick,
                self.height,
                self.rx,
                self.tx,
                self.ty,
                self.tz,
                self.rz,
            )
        )
        fw.write("\n")

    def genMatSurf(self, fw, ns):
        fw.write("void plastic surf_mat_%i\n" % ns)
        fw.write("0\n")
        fw.write("0\n")
        fw.write("5  %1.3f %1.3f %1.3f 0 0\n" % (self.intAlb, self.intAlb, self.intAlb))
        fw.write("\n")

    def extEnergyBal(self, extRHTC, skyTemp, extCHTC, tempOutAir, irrSurfExt):
        if self.bc:
            self.extFlux = (
                self.svf * self.emis[0] * extRHTC * (skyTemp - self.extTemp)
                + extCHTC * (tempOutAir - self.extTemp)
                + self.solAbs[0] * irrSurfExt
            )

    def intEnergyBal(self, RHTCsurf, tempRad, CHTCsurf, tempInAir, irrSurfInt):
        self.intFlux = (
            self.emis[self.numLayer - 1] * RHTCsurf * (tempRad - self.intTemp)
            + CHTCsurf * (tempInAir - self.intTemp)
            + self.solAbs[self.numLayer - 1] * irrSurfInt
        )

    def dirRot(
        self, tx, ty, tz
    ):  # rotation for direction(unit) vector, following right hand rule
        pi = math.pi / 180
        # xorg= cos(pi*(-self.rz))*tx - sin(pi*(-self.rz))*cos(pi*(-self.rx))*ty + sin(pi*(-self.rz))*sin(pi*(-self.rx))*tz
        # yorg= sin(pi*(-self.rz))*tx + cos(pi*(-self.rz))*cos(pi*(-self.rx))*ty - cos(pi*(-self.rz))*sin(pi*(-self.rx))*tz
        # zorg= sin(pi*(-self.rx))*ty + cos(pi*(-self.rx))*tz
        xorg = (
            cos(pi * (-self.rz)) * tx
            + sin(pi * (-self.rz)) * cos(pi * (-self.rx)) * ty
            + sin(pi * (-self.rz)) * sin(pi * (-self.rx)) * tz
        )
        yorg = (
            -sin(pi * (-self.rz)) * tx
            + cos(pi * (-self.rz)) * cos(pi * (-self.rx)) * ty
            + cos(pi * (-self.rz)) * sin(pi * (-self.rx)) * tz
        )
        zorg = -sin(pi * (-self.rx)) * ty + cos(pi * (-self.rx)) * tz

        return xorg, yorg, zorg

    def orgRot(self, tx, ty, tz):  # rotation for coordinate
        pi = math.pi / 180
        # rotate with x-axis
        x = tx
        y = cos(pi * (self.rx)) * ty - sin(pi * (self.rx)) * tz
        z = sin(pi * (self.rx)) * ty + cos(pi * (self.rx)) * tz
        # move coordinates
        x += self.tx
        y += self.ty
        z += self.tz
        # rotate with z-axis
        xorg = cos(pi * (self.rz)) * x - sin(pi * (self.rz)) * y
        yorg = sin(pi * (self.rz)) * x + cos(pi * (self.rz)) * y
        zorg = z
        return xorg, yorg, zorg


class frame(object):
    def __init__(self, vect, winVect, surf, frameID, thermFlag=False, iniTemp=0):
        self.length = vect[0]
        self.height = vect[1]
        self.thick = vect[2]
        self.surf = vect[3]
        self.xOffset = vect[4]
        self.zOffset = vect[5]
        self.reveal = vect[6]
        self.res = 1.0 / vect[7] - 0.13 - 0.04
        self.extAlb = vect[8]
        self.intAlb = vect[9]
        self.svf = vect[10]
        self.tilt = surf[int(self.surf)].tilt
        self.rz = surf[int(self.surf)].rz
        self.rx = surf[int(self.surf)].rx
        self.ry = surf[int(self.surf)].ry
        self.area = vect[0] * vect[1]
        for i in range(winVect.shape[0]):
            if winVect[i, 2] == frameID:
                self.area = self.area - winVect[i, 0] * winVect[i, 1]
        if thermFlag:
            self.extTemp = iniTemp
            self.intTemp = iniTemp
            # ---------------------------------------------
            self.extEmis = 0.9
            self.intEmis = 0.9
            # ---------------------------------------------

    def genMatFrame(self, fw, nf):
        fw.write("void plastic frame_mat_%i\n" % nf)
        fw.write("0\n")
        fw.write("0\n")
        fw.write("5  %1.3f %1.3f %1.3f 0 0\n" % (self.intAlb, self.intAlb, self.intAlb))
        fw.write("\n")


class win(object):
    def __init__(
        self,
        vect,
        frame,
        thermFlag=False,
        iniTemp=0,
        tmxs=0,
        rmxs=0,
        calorim=0,
        isoFlag=False,
        conWinVect=0,
        matGas=0,
        matGlz=0,
        matBSDF=0,
    ):
        self.length = vect[0]
        self.height = vect[1]
        self.frame = vect[2]
        self.xOffset = vect[3]
        self.zOffset = vect[4]
        self.reveal = vect[5]
        self.svf = vect[6]
        self.area = vect[0] * vect[1]
        self.tilt = frame[int(self.frame)].tilt
        self.rz = frame[int(self.frame)].rz
        self.con = vect[7]
        if frame[int(self.frame)].rx == 0.0 and frame[int(self.frame)].ry == 0.0:
            self.vertical = 1
        else:
            self.vertical = 0
        # ---------------------------------------------
        if thermFlag:
            # ---------------------------------------------
            self.extEmis = 0.9
            self.intEmis = 0.9
            # ---------------------------------------------
            self.extTemp = iniTemp
            self.intTemp = iniTemp
            # ---------------------------------------------
            # TK's black box model
            self.kuhn11(
                tmxs[int(vect[7]), :, :],
                rmxs[int(vect[7]), :, :],
                calorim[int(vect[7]), :],
            )
        # ---------------------------------------------
        if isoFlag:
            self.gapTemp = ones(int(conWinVect[int(vect[7]), 0])) * iniTemp
            self.paneTemp = ones(int(conWinVect[int(vect[7]), 0])) * iniTemp
            self.iso(conWinVect[int(vect[7]), :], matGas, matGlz, matBSDF)
            # ---------------------------------------------
            self.extTemp = iniTemp
            self.intTemp = iniTemp
            # secondary heat transferred to the room
            self.qInt = 0.0
            self.uValue = 0.0

    def kuhn11(self, tmxs, rmxs, calorim):
        # Generation of a virtual two-layer virtual system based on Kuhn et al 2011
        resInt = 0.13
        resExt = 0.04
        # leeward resExt = 0.065, windward resExt = 0.045, according to TARP (EP EngineeringRef)
        # solidAngleMtx = calcs.klemsSolidAngles()
        self.res = 1.0 / (calorim[0] * (-0.8395 + 0.0065 * 297.15)) - resExt - resInt
        self.res = self.res
        # ~ self.res = 1/calorim[0]-resExt-resInt
        self.absLayInt = zeros(145)
        self.absLayExt = zeros(145)
        transTot = zeros(145)
        reflTot = zeros(145)
        for i in range(145):
            gValue = calorim[i + 1]
            transTot[i] = sum(tmxs[:, i])
            reflTot[i] = sum(rmxs[:, i])
            absTot = 1 - transTot[i] - reflTot[i]
            solHeatIn = gValue - transTot[i]
            self.absLayInt[i] = (
                solHeatIn * (resExt + self.res + resInt) - absTot * resExt
            ) / self.res
            self.absLayExt[i] = (
                -solHeatIn * (resExt + self.res + resInt) + absTot * (resExt + self.res)
            ) / self.res
            self.absLayInt[i] = self.absLayInt[i]  # ELIMINAR!!!!
            self.absLayExt[i] = absTot - self.absLayInt[i]  # ELIMINAR !!!
            # weight by the solid angle and apply irradiance/radiance factor pi
            # ~ self.absLayInt[i] = self.absLayInt[i]*solidAngleMtx[i,i]/(2*pi)*pi
            # ~ self.absLayExt[i] = self.absLayExt[i]*solidAngleMtx[i,i]/(2*pi)*pi
            # ~ self.absLayInt[i] = self.absLayInt[i]/145*pi
            # ~ self.absLayExt[i] = self.absLayExt[i]/145*pi
        # ~ print('g0',calorim[1],'g60',calorim[calcs.searchPatch(0,60)+1])
        # ~ print('absTot0',1-transTot[0]-reflTot[0],'absTot60',1-transTot[calcs.searchPatch(0,60)]-reflTot[calcs.searchPatch(0,60)])
        # ~ print('transTot0',transTot[0],'transTot60',transTot[calcs.searchPatch(0,60)])
        # ~ print('absLayInt0',self.absLayInt[0],'absLayInt60',self.absLayInt[calcs.searchPatch(0,60)])
        # ~ print('absLayExt0',self.absLayExt[0],'absLayExt60',self.absLayExt[calcs.searchPatch(0,60)])
        # ~ tout = 32
        # ~ hout = 53
        # ~ resExtTest = 1/hout
        # ~ tmean = (tout+25)/2+273.15
        # ~ uvalueTest = calorim[0]*(-0.8395+0.0065*tmean)
        # ~ print('Uvalue',calorim[0],uvalueTest)
        # ~ resTest = 1/uvalueTest-resExtTest-resInt
        # ~ print('Rs',self.res,resTest)
        # ~ qiTest = (resExtTest*self.absLayExt[0]+(resExtTest+resTest)*self.absLayInt[0])*uvalueTest
        # ~ qi000 = (resExt*self.absLayExt[0]+(resExt+self.res)*self.absLayInt[0])*calorim[0]
        # ~ print('g0',calorim[1],transTot[0]+qiTest)
        # ~ exit()

    def iso(self, conWin, matGas, matGlz, matBSDF):
        # ---------------------------------------------
        # maxNumPane = 4
        self.numPane = int(conWin[0])
        self.gapThick = ones(self.numPane) * -77777
        self.airFrac = ones(self.numPane) * -77777
        self.argFrac = ones(self.numPane) * -77777
        self.paneThick = zeros(self.numPane)
        self.paneCond = zeros(self.numPane)
        # self.volHeat =  zeros(self.numPane)
        self.frontEmis = zeros(self.numPane)
        self.backEmis = zeros(self.numPane)
        # self.topOpen =  zeros(self.numPane)
        # self.bottomOpen =  zeros(self.numPane)
        # self.leftOpen =  zeros(self.numPane)
        # self.rightOpen =  zeros(self.numPane)
        # self.frontOpen =  zeros(self.numPane)
        paneCount = 0
        paneType = int(conWin[1])
        paneValue = int(conWin[2])
        if paneType == 1:
            self.paneThick[paneCount] = matGlz[paneValue, 0]
            self.paneCond[paneCount] = matGlz[paneValue, 10]
            # self.volHeat[paneCount] = matGlz[paneValue,12]
            self.frontEmis[paneCount] = matGlz[paneValue, 8]
            self.backEmis[paneCount] = matGlz[paneValue, 9]
        if paneType == 4:
            self.paneThick[paneCount] = matBSDF[paneValue, 0]
            self.paneCond[paneCount] = matBSDF[paneValue, 1]
            # self.volHeat[paneCount] = matBSDF[paneValue,9]
            self.frontEmis[paneCount] = matBSDF[paneValue, 2]
            self.backEmis[paneCount] = matBSDF[paneValue, 3]
            # self.topOpen[paneCount] = matBSDF[paneValue,4]
            # self.bottomOpen[paneCount] = matBSDF[paneValue,5]
            # self.leftOpen[paneCount] = matBSDF[paneValue,6]
            # self.rightOpen[paneCount] = matBSDF[paneValue,7]
            # self.frontOpen[paneCount] = matBSDF[paneValue,8]
        for i in range(1, 2 * self.numPane - 1):
            paneType = int(conWin[2 * i + 1])
            paneValue = int(conWin[2 * i + 2])
            if paneType == 0:
                paneCount = paneCount + 1
                self.gapThick[paneCount] = matGas[paneValue, 0]
                self.airFrac[paneCount] = matGas[paneValue, 1]
                self.argFrac[paneCount] = matGas[paneValue, 2]
            else:
                if paneType == 1:
                    self.paneThick[paneCount] = matGlz[paneValue, 0]
                    self.paneCond[paneCount] = matGlz[paneValue, 10]
                    # self.volHeat[paneCount] = matGlz[paneValue,12]
                    self.frontEmis[paneCount] = matGlz[paneValue, 8]
                    self.backEmis[paneCount] = matGlz[paneValue, 9]
                if paneType == 4:
                    self.paneThick[paneCount] = matBSDF[paneValue, 0]
                    self.paneCond[paneCount] = matBSDF[paneValue, 1]
                    # self.volHeat[paneCount] = matBSDF[paneValue,9]
                    self.frontEmis[paneCount] = matBSDF[paneValue, 2]
                    self.backEmis[paneCount] = matBSDF[paneValue, 3]
                    # self.topOpen[paneCount] = matBSDF[paneValue,4]
                    # self.bottomOpen[paneCount] = matBSDF[paneValue,5]
                    # self.leftOpen[paneCount] = matBSDF[paneValue,6]
                    # self.rightOpen[paneCount] = matBSDF[paneValue,7]
                    # self.frontOpen[paneCount] = matBSDF[paneValue,8]
        # ---------------------------------------------
        self.extEmis = self.frontEmis[0]
        self.intEmis = self.backEmis[self.numPane - 1]

    def genPolygon(self, fw, mat, name, surf, frame):
        ux, uy, uz = surf.dirRot(0, 0, 1)
        # x1,y1,z1 = surf.orgRot(frame.xOffset+self.xOffset+self.length,frame.reveal+self.reveal,frame.zOffset+self.height+self.zOffset)
        # x2,y2,z2 = surf.orgRot(frame.xOffset+self.xOffset,frame.reveal+self.reveal,frame.zOffset+self.height+self.zOffset)
        # x3,y3,z3 = surf.orgRot(frame.xOffset+self.xOffset,frame.reveal+self.reveal,frame.zOffset+self.zOffset)
        # x4,y4,z4 = surf.orgRot(frame.xOffset+self.xOffset+self.length,frame.reveal+self.reveal,frame.zOffset+self.zOffset)
        x1, y1, z1 = surf.orgRot(
            frame.xOffset + self.xOffset,
            frame.reveal + self.reveal,
            frame.zOffset + self.zOffset,
        )
        x2, y2, z2 = surf.orgRot(
            frame.xOffset + self.xOffset,
            frame.reveal + self.reveal,
            frame.zOffset + self.height + self.zOffset,
        )
        x3, y3, z3 = surf.orgRot(
            frame.xOffset + self.xOffset + self.length,
            frame.reveal + self.reveal,
            frame.zOffset + self.height + self.zOffset,
        )
        x4, y4, z4 = surf.orgRot(
            frame.xOffset + self.xOffset + self.length,
            frame.reveal + self.reveal,
            frame.zOffset + self.zOffset,
        )
        fw.write("#@rfluxmtx h=kf u=%1.3f,%1.3f,%1.3f\n" % (ux, uy, uz))
        fw.write("\n")
        fw.write("%s polygon %s\n" % (mat, name))
        fw.write("0\n")
        fw.write("0\n")
        fw.write("12  %1.3f %1.3f %1.3f \n" % (x1, y1, z1))
        fw.write("    %1.3f %1.3f %1.3f \n" % (x2, y2, z2))
        fw.write("    %1.3f %1.3f %1.3f \n" % (x3, y3, z3))
        fw.write("    %1.3f %1.3f %1.3f \n" % (x4, y4, z4))
        fw.write("\n")

    def extEnergyBal(self, extRHTC, skyTemp, extCHTC, tempOutAir, conWin):
        self.extFlux = self.svf * conWin[self.con].frontEmis[0] * extRHTC * (
            skyTemp - self.extTemp
        ) + extCHTC * (tempOutAir - self.extTemp)

    def intEnergyBal(self, RHTCwin, tempRad, CHTCwin, tempInAir, irrWinInt, conWin):
        self.intFlux = (
            conWin[self.con].backEmis[conWin[self.con].numPane - 1]
            * RHTCwin
            * (tempRad - self.intTemp)
            + CHTCwin * (tempInAir - self.intTemp)
            + conWin[self.con].absBack[conWin[self.con].numPane - 1, 0] * irrWinInt
        )


class conWin(object):
    def __init__(self, vect, i, input_path, absFront_path, absBack_path):
        # maxNumPane = 4
        self.numPane = int(vect[0])
        self.absFront = zeros([self.numPane, 1, 145])
        self.absBack = zeros([self.numPane, 1, 145])
        for k in range(self.numPane):
            self.absFront[k, :, :] = atleast_2d(
                genfromtxt("%s" % absFront_path[i, k], delimiter=",")
            )
            self.absBack[k, :, :] = atleast_2d(
                genfromtxt("%s" % absBack_path[i, k], delimiter=",")
            )

    def updateTemp(self, tempOutAir, tempInAir):
        self.temp[0] = tempOutAir
        self.temp[self.numPane - 1] = tempInAir


class cntrlMtx(object):
    def __init__(self, cntrlMtx_path, daylightFlag, glareFlag, thermFlag, numWin):
        with open(cntrlMtx_path) as f:
            lines = f.readlines()

        aux1 = lines[1].split(",")
        aux2 = lines[2].split(",")
        aux3 = lines[3].split(",")
        # number of conditions
        numCond = len(lines) - 3
        # number of variables
        numVar = len(aux1)
        # number of states
        numState = len(aux3) - numVar
        if numState != 1 and numState != numWin:
            raise TypeError(
                "Mismatch between control matrix states and number of windows"
            )
            # sys.exit()
        # control variables
        self.var = zeros(numVar, int)
        # control setpoints
        self.value = zeros(numVar)
        self.sp = zeros(numVar)
        for j in range(numVar):
            self.var[j] = int(aux1[j])
            if (
                (self.var[j] == 2 and not daylightFlag)
                or (self.var[j] == 3 and not glareFlag)
                or (self.var[j] > 3 and not thermFlag)
            ):
                raise TypeError(
                    "Mismatch between control matrix variables and simulation options"
                )
                # sys.exit()
            self.sp[j] = float(aux2[j])
        # matrix of conditions
        self.mtx = zeros([numCond, numVar])
        # vector of states
        self.state = zeros([numCond, numState])
        for i in range(3, len(lines)):
            aux = lines[i].split(",")
            for j in range(len(aux) - numState):
                self.mtx[i - 3, j] = float(aux[j])
            for k in range(numState):
                self.state[i - 3, k] = int(aux[len(aux) - numState + k])
