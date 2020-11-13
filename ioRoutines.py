#!/usr/bin/env python
# -*- coding: utf-8 -*-

from subprocess import run, PIPE
from numpy import matrix, zeros
from xml.etree import ElementTree
from lxml import etree
import re


def shell(cmd, outputFlag=False):
    print(cmd)
    stdout = PIPE if outputFlag else None
    completedProcess = run(cmd, stdout=stdout, stderr=PIPE, shell=True, check=True)
    if outputFlag:
        return completedProcess.stdout.decode("utf-8")


def shell_output(cmd):
    return shell(cmd, True)


def findPathName(s, first, last):
    try:
        start = s.rindex(first) + len(first)
        end = s.rindex(last, start)
        return s[:start], s[start:end]
    except ValueError:
        return ""


def skyMtx(file_name, iniHour, totHour, dayHour, sunAltitude):
    with open(file_name, "r") as _file:
        M = list()
        for s in _file:
            M.append(s)

    mtx = matrix(zeros([dayHour, 2306]))
    for j in range(2306):
        dh = 0
        for i in range(8760):
            k = j * 8761 + i
            if i >= iniHour and i < iniHour + totHour:
                if sunAltitude[i - iniHour] > 0.0:
                    mtx[dh, j] = float(M[k].split()[0])
                    dh = dh + 1
    return mtx


def xmlBlank(file_name, wLength, wLengthDataDir):
    with open(file_name, "rt") as f:
        tree = ElementTree.parse(f)
    flagV = False
    flagT = False
    mtx = zeros([145, 145])
    for node in tree.iter():
        if node.text == wLength:
            flagV = True
        if node.text == wLengthDataDir:
            flagT = True
        if node.text and flagT and flagV and node.tag.endswith("ScatteringData"):
            a = node.text.split(" ")  # string to array
            for i in range(145):  # array to matrix
                for j in range(145):
                    k = i * 145 + j
                    mtx[i, j] = float(a[k + 1])
            flagT, flagV = False, False
        elif node.tag.endswith("ScatteringData"):
            flagT, flagV = False, False
    return mtx


def xmlTab(file_name, wLength, wLengthDataDir):
    with open(file_name, "rt") as f:
        tree = ElementTree.parse(f)
    flagV = False
    flagT = False
    mtx = zeros([145, 145])
    for node in tree.iter():
        if node.text == wLength:
            flagV = True
        if node.text == wLengthDataDir:
            flagT = True
        if node.text and flagT and flagV and node.tag.endswith("ScatteringData"):
            a = node.text.split("\t")  # string to array
            for i in range(145):  # array to matrix
                for j in range(145):
                    k = i * 145 + j
                    mtx[i, j] = float(a[k + 1])
            flagT, flagV = False, False
        elif node.tag.endswith("ScatteringData"):
            flagT, flagV = False, False
    return mtx


def xml(file_name, wLength, wLengthDataDir):
    tree = etree.parse(file_name)
    root = tree.getroot()
    for wavelengthData in root.findall('.//{*}WavelengthData'):
        for wavelength in wavelengthData.findall('{*}Wavelength'):
            if wavelength.text == wLength:
                for wavelength_data_direction in wavelengthData.findall('.//{*}WavelengthDataDirection'):
                    if wavelength_data_direction.text == wLengthDataDir:
                         for scattering_data in wavelengthData.findall('.//{*}ScatteringData'):
                             floats = []
                             try:
                                for x in scattering_data.text.split():
                                    floats.append(x.replace(',',''))
                             except:
                                for x in scattering_data.text.split(','):
                                    floats.append(x)
    mtx = zeros([145, 145])
    for i in range(145):
        for j in range(145):
            k = i * 145 + j
            mtx[i, j] = floats[k]
    return mtx


#~ def xml(file_name, wLength, wLengthDataDir):
    #~ with open(file_name, "rt") as f:
        #~ tree = ElementTree.parse(f)
    #~ flagV = False
    #~ flagT = False
    #~ mtx = zeros([145, 145])
    #~ for node in tree.iter():
        #~ if node.text == wLength:
            #~ flagV = True
        #~ if node.text == wLengthDataDir:
            #~ flagT = True
        #~ if node.text and flagT and flagV and node.tag.endswith("ScatteringData"):
            #~ a = node.text.split(",")  # string to array
            #~ for i in range(145):  # array to matrix
                #~ for j in range(145):
                    #~ k = i * 145 + j
                    #~ mtx[i, j] = a[k]
            #~ flagT, flagV = False, False
        #~ elif node.tag.endswith("ScatteringData"):
            #~ flagT, flagV = False, False
    #~ return mtx


def writeDat(file_name, numConWin, calorim):
    for i in range(numConWin):
        with open("%s" % (file_name[i]), "w") as fw:
            fw.write(
                "# uValue, gValue(klemsPatch=0), gValue(klemsPatch=1), gValue(klemsPatch=2), gValue(klemsPatch=3),... gValue(klemsPatch=144)\n"
            )
            for j in range(len(calorim[i, :]) - 1):
                fw.write("%1.4f, " % calorim[i, j])
            fw.write("%1.4f " % calorim[i, len(calorim[i, :]) - 1])


def writeKuhn11Abs(numConWin, win, tmxs, rmxs, calorim, output_path):
    winAux = win
    for i in range(numConWin):
        winAux.kuhn11(tmxs[i, :, :], rmxs[i, :, :], calorim[i, :])
        with open("%sbb_absFront_%i.vct" % (output_path, i), "w") as ff, open(
            "%sbb_absBack_%i.vct" % (output_path, i), "w"
        ) as fb:
            for j in range(144):
                ff.write("%1.4f, " % winAux.absLayExt[j])
                fb.write("%1.4f, " % winAux.absLayInt[j])
            ff.write("%1.4f " % winAux.absLayExt[144])
            fb.write("%1.4f " % winAux.absLayInt[144])


def writeXml(
    file_name,
    transFrontMtxV,
    refFrontMtxV,
    transBackMtxV,
    refBackMtxV,
    transFrontMtxS,
    refFrontMtxS,
    transBackMtxS,
    refBackMtxS,
):
    length = 145
    with open(file_name, "w") as fw:
        fw.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fw.write('<WindowElement xmlns="http://www.ise.fraunhofer.de">\n')
        fw.write(
            "<!-- File produced by Fener, Fraunhofer-Institut für Solare Energiesysteme -->\n"
        )
        fw.write("<WindowElementType>System</WindowElementType>\n")
        fw.write("<FileType>BSDF</FileType>\n")
        fw.write("<Optical>\n")
        fw.write("<Layer>\n")
        fw.write("	<Material>\n")
        fw.write("		<Name>Name</Name>\n")
        fw.write("		<Manufacturer>Manufacturer</Manufacturer>\n")
        fw.write('		<Thickness unit="meter"></Thickness>\n')
        fw.write('		<Width unit="meter"></Width>\n')
        fw.write('		<Height unit="meter"></Height>\n')
        fw.write("		<DeviceType>Other</DeviceType>\n")
        fw.write("	</Material>\n")
        fw.write("	<DataDefinition>\n")
        fw.write("		<IncidentDataStructure>Columns</IncidentDataStructure>\n")
        fw.write("		<AngleBasis>\n")
        fw.write("			<AngleBasisName>LBNL/Klems Full</AngleBasisName>\n")
        fw.write("			<AngleBasisBlock>\n")
        fw.write("			<Theta>0</Theta>\n")
        fw.write("			<nPhis>1</nPhis>\n")
        fw.write("			<ThetaBounds>\n")
        fw.write("				<LowerTheta>0</LowerTheta>\n")
        fw.write("				<UpperTheta>5</UpperTheta>\n")
        fw.write("			</ThetaBounds>\n")
        fw.write("			</AngleBasisBlock>\n")
        fw.write("			<AngleBasisBlock>\n")
        fw.write("			<Theta>10</Theta>\n")
        fw.write("			<nPhis>8</nPhis>\n")
        fw.write("			<ThetaBounds>\n")
        fw.write("				<LowerTheta>5</LowerTheta>\n")
        fw.write("				<UpperTheta>15</UpperTheta>\n")
        fw.write("			</ThetaBounds>\n")
        fw.write("			</AngleBasisBlock>\n")
        fw.write("			<AngleBasisBlock>\n")
        fw.write("			<Theta>20</Theta>\n")
        fw.write("			<nPhis>16</nPhis>\n")
        fw.write("			<ThetaBounds>\n")
        fw.write("				<LowerTheta>15</LowerTheta>\n")
        fw.write("				<UpperTheta>25</UpperTheta>\n")
        fw.write("			</ThetaBounds>\n")
        fw.write("			</AngleBasisBlock>\n")
        fw.write("			<AngleBasisBlock>\n")
        fw.write("			<Theta>30</Theta>\n")
        fw.write("			<nPhis>20</nPhis>\n")
        fw.write("			<ThetaBounds>\n")
        fw.write("				<LowerTheta>25</LowerTheta>\n")
        fw.write("				<UpperTheta>35</UpperTheta>\n")
        fw.write("			</ThetaBounds>\n")
        fw.write("			</AngleBasisBlock>\n")
        fw.write("			<AngleBasisBlock>\n")
        fw.write("			<Theta>40</Theta>\n")
        fw.write("			<nPhis>24</nPhis>\n")
        fw.write("			<ThetaBounds>\n")
        fw.write("				<LowerTheta>35</LowerTheta>\n")
        fw.write("				<UpperTheta>45</UpperTheta>\n")
        fw.write("			</ThetaBounds>\n")
        fw.write("			</AngleBasisBlock>\n")
        fw.write("			<AngleBasisBlock>\n")
        fw.write("			<Theta>50</Theta>\n")
        fw.write("			<nPhis>24</nPhis>\n")
        fw.write("			<ThetaBounds>\n")
        fw.write("				<LowerTheta>45</LowerTheta>\n")
        fw.write("				<UpperTheta>55</UpperTheta>\n")
        fw.write("			</ThetaBounds>\n")
        fw.write("			</AngleBasisBlock>\n")
        fw.write("			<AngleBasisBlock>\n")
        fw.write("			<Theta>60</Theta>\n")
        fw.write("			<nPhis>24</nPhis>\n")
        fw.write("			<ThetaBounds>\n")
        fw.write("				<LowerTheta>55</LowerTheta>\n")
        fw.write("				<UpperTheta>65</UpperTheta>\n")
        fw.write("			</ThetaBounds>\n")
        fw.write("			</AngleBasisBlock>\n")
        fw.write("			<AngleBasisBlock>\n")
        fw.write("			<Theta>70</Theta>\n")
        fw.write("			<nPhis>16</nPhis>\n")
        fw.write("			<ThetaBounds>\n")
        fw.write("				<LowerTheta>65</LowerTheta>\n")
        fw.write("				<UpperTheta>75</UpperTheta>\n")
        fw.write("			</ThetaBounds>\n")
        fw.write("			</AngleBasisBlock>\n")
        fw.write("			<AngleBasisBlock>\n")
        fw.write("			<Theta>82.5</Theta>\n")
        fw.write("			<nPhis>12</nPhis>\n")
        fw.write("			<ThetaBounds>\n")
        fw.write("				<LowerTheta>75</LowerTheta>\n")
        fw.write("				<UpperTheta>90</UpperTheta>\n")
        fw.write("			</ThetaBounds>\n")
        fw.write("			</AngleBasisBlock>\n")
        fw.write("		</AngleBasis>\n")
        fw.write("	</DataDefinition>\n")
        fw.write("	<WavelengthData>\n")
        fw.write("		<LayerNumber>System</LayerNumber>\n")
        fw.write('		<Wavelength unit="Integral">Visible</Wavelength>\n')
        fw.write("		<SourceSpectrum>CIE Illuminant D65 1nm.ssp</SourceSpectrum>\n")
        fw.write("		<DetectorSpectrum>ASTM E308 1931 Y.dsp</DetectorSpectrum>\n")
        fw.write("		<WavelengthDataBlock>\n")
        fw.write(
            "			<WavelengthDataDirection>Transmission Front</WavelengthDataDirection>\n"
        )
        fw.write("			<ColumnAngleBasis>LBNL/Klems Full</ColumnAngleBasis>\n")
        fw.write("			<RowAngleBasis>LBNL/Klems Full</RowAngleBasis>\n")
        fw.write("			<ScatteringDataType>BTDF</ScatteringDataType>\n")
        fw.write("			<ScatteringData>\n")
        for i in range(length):
            for j in range(length):
                fw.write("%f,\n" % transFrontMtxV[i, j])
            fw.write("\n")
        fw.write("			</ScatteringData>\n")
        fw.write("		</WavelengthDataBlock>\n")
        fw.write("	</WavelengthData>\n")
        fw.write("	<WavelengthData>\n")
        fw.write("		<LayerNumber>System</LayerNumber>\n")
        fw.write('		<Wavelength unit="Integral">Visible</Wavelength>\n')
        fw.write("		<SourceSpectrum>CIE Illuminant D65 1nm.ssp</SourceSpectrum>\n")
        fw.write("		<DetectorSpectrum>ASTM E308 1931 Y.dsp</DetectorSpectrum>\n")
        fw.write("		<WavelengthDataBlock>\n")
        fw.write(
            "			<WavelengthDataDirection>Reflection Front</WavelengthDataDirection>\n"
        )
        fw.write("			<ColumnAngleBasis>LBNL/Klems Full</ColumnAngleBasis>\n")
        fw.write("			<RowAngleBasis>LBNL/Klems Full</RowAngleBasis>\n")
        fw.write("			<ScatteringDataType>BTDF</ScatteringDataType>\n")
        fw.write("			<ScatteringData>\n")
        for i in range(length):
            for j in range(length):
                fw.write("%f,\n" % refFrontMtxV[i, j])
            fw.write("\n")
        fw.write("			</ScatteringData>\n")
        fw.write("		</WavelengthDataBlock>\n")
        fw.write("	</WavelengthData>\n")
        fw.write("	<WavelengthData>\n")
        fw.write("		<LayerNumber>System</LayerNumber>\n")
        fw.write('		<Wavelength unit="Integral">Visible</Wavelength>\n')
        fw.write("		<SourceSpectrum>CIE Illuminant D65 1nm.ssp</SourceSpectrum>\n")
        fw.write("		<DetectorSpectrum>ASTM E308 1931 Y.dsp</DetectorSpectrum>\n")
        fw.write("		<WavelengthDataBlock>\n")
        fw.write(
            "			<WavelengthDataDirection>Transmission Back</WavelengthDataDirection>\n"
        )
        fw.write("			<ColumnAngleBasis>LBNL/Klems Full</ColumnAngleBasis>\n")
        fw.write("			<RowAngleBasis>LBNL/Klems Full</RowAngleBasis>\n")
        fw.write("			<ScatteringDataType>BTDF</ScatteringDataType>\n")
        fw.write("			<ScatteringData>\n")
        for i in range(length):
            for j in range(length):
                fw.write("%f,\n" % transBackMtxV[i, j])
            fw.write("\n")
        fw.write("			</ScatteringData>\n")
        fw.write("		</WavelengthDataBlock>\n")
        fw.write("	</WavelengthData>\n")
        fw.write("	<WavelengthData>\n")
        fw.write("		<LayerNumber>System</LayerNumber>\n")
        fw.write('		<Wavelength unit="Integral">Visible</Wavelength>\n')
        fw.write("		<SourceSpectrum>CIE Illuminant D65 1nm.ssp</SourceSpectrum>\n")
        fw.write("		<DetectorSpectrum>ASTM E308 1931 Y.dsp</DetectorSpectrum>\n")
        fw.write("		<WavelengthDataBlock>\n")
        fw.write(
            "			<WavelengthDataDirection>Reflection Back</WavelengthDataDirection>\n"
        )
        fw.write("			<ColumnAngleBasis>LBNL/Klems Full</ColumnAngleBasis>\n")
        fw.write("			<RowAngleBasis>LBNL/Klems Full</RowAngleBasis>\n")
        fw.write("			<ScatteringDataType>BTDF</ScatteringDataType>\n")
        fw.write("			<ScatteringData>\n")
        for i in range(length):
            for j in range(length):
                fw.write("%f,\n" % refBackMtxV[i, j])
            fw.write("\n")
        fw.write("			</ScatteringData>\n")
        fw.write("		</WavelengthDataBlock>\n")
        fw.write("	</WavelengthData>\n")

        fw.write("	<WavelengthData>\n")
        fw.write("		<LayerNumber>System</LayerNumber>\n")
        fw.write('		<Wavelength unit="Integral">Solar</Wavelength>\n')
        fw.write("		<SourceSpectrum>CIE Illuminant D65 1nm.ssp</SourceSpectrum>\n")
        fw.write("		<DetectorSpectrum>ASTM E308 1931 Y.dsp</DetectorSpectrum>\n")
        fw.write("		<WavelengthDataBlock>\n")
        fw.write(
            "			<WavelengthDataDirection>Transmission Front</WavelengthDataDirection>\n"
        )
        fw.write("			<ColumnAngleBasis>LBNL/Klems Full</ColumnAngleBasis>\n")
        fw.write("			<RowAngleBasis>LBNL/Klems Full</RowAngleBasis>\n")
        fw.write("			<ScatteringDataType>BTDF</ScatteringDataType>\n")
        fw.write("			<ScatteringData>\n")
        for i in range(length):
            for j in range(length):
                fw.write("%f,\n" % transFrontMtxS[i, j])
            fw.write("\n")
        fw.write("			</ScatteringData>\n")
        fw.write("		</WavelengthDataBlock>\n")
        fw.write("	</WavelengthData>\n")
        fw.write("	<WavelengthData>\n")
        fw.write("		<LayerNumber>System</LayerNumber>\n")
        fw.write('		<Wavelength unit="Integral">Solar</Wavelength>\n')
        fw.write("		<SourceSpectrum>CIE Illuminant D65 1nm.ssp</SourceSpectrum>\n")
        fw.write("		<DetectorSpectrum>ASTM E308 1931 Y.dsp</DetectorSpectrum>\n")
        fw.write("		<WavelengthDataBlock>\n")
        fw.write(
            "			<WavelengthDataDirection>Reflection Front</WavelengthDataDirection>\n"
        )
        fw.write("			<ColumnAngleBasis>LBNL/Klems Full</ColumnAngleBasis>\n")
        fw.write("			<RowAngleBasis>LBNL/Klems Full</RowAngleBasis>\n")
        fw.write("			<ScatteringDataType>BTDF</ScatteringDataType>\n")
        fw.write("			<ScatteringData>\n")
        for i in range(length):
            for j in range(length):
                fw.write("%f,\n" % refFrontMtxS[i, j])
            fw.write("\n")
        fw.write("			</ScatteringData>\n")
        fw.write("		</WavelengthDataBlock>\n")
        fw.write("	</WavelengthData>\n")
        fw.write("	<WavelengthData>\n")
        fw.write("		<LayerNumber>System</LayerNumber>\n")
        fw.write('		<Wavelength unit="Integral">Solar</Wavelength>\n')
        fw.write("		<SourceSpectrum>CIE Illuminant D65 1nm.ssp</SourceSpectrum>\n")
        fw.write("		<DetectorSpectrum>ASTM E308 1931 Y.dsp</DetectorSpectrum>\n")
        fw.write("		<WavelengthDataBlock>\n")
        fw.write(
            "			<WavelengthDataDirection>Transmission Back</WavelengthDataDirection>\n"
        )
        fw.write("			<ColumnAngleBasis>LBNL/Klems Full</ColumnAngleBasis>\n")
        fw.write("			<RowAngleBasis>LBNL/Klems Full</RowAngleBasis>\n")
        fw.write("			<ScatteringDataType>BTDF</ScatteringDataType>\n")
        fw.write("			<ScatteringData>\n")
        for i in range(length):
            for j in range(length):
                fw.write("%f,\n" % transBackMtxS[i, j])
            fw.write("\n")
        fw.write("			</ScatteringData>\n")
        fw.write("		</WavelengthDataBlock>\n")
        fw.write("	</WavelengthData>\n")
        fw.write("	<WavelengthData>\n")
        fw.write("		<LayerNumber>System</LayerNumber>\n")
        fw.write('		<Wavelength unit="Integral">Solar</Wavelength>\n')
        fw.write("		<SourceSpectrum>CIE Illuminant D65 1nm.ssp</SourceSpectrum>\n")
        fw.write("		<DetectorSpectrum>ASTM E308 1931 Y.dsp</DetectorSpectrum>\n")
        fw.write("		<WavelengthDataBlock>\n")
        fw.write(
            "			<WavelengthDataDirection>Reflection Back</WavelengthDataDirection>\n"
        )
        fw.write("			<ColumnAngleBasis>LBNL/Klems Full</ColumnAngleBasis>\n")
        fw.write("			<RowAngleBasis>LBNL/Klems Full</RowAngleBasis>\n")
        fw.write("			<ScatteringDataType>BTDF</ScatteringDataType>\n")
        fw.write("			<ScatteringData>\n")
        for i in range(length):
            for j in range(length):
                fw.write("%f,\n" % refBackMtxS[i, j])
            fw.write("\n")
        fw.write("			</ScatteringData>\n")
        fw.write("		</WavelengthDataBlock>\n")
        fw.write("	</WavelengthData>\n")
        fw.write("</Layer>\n")
        fw.write("</Optical>\n")
        fw.write("</WindowElement>\n")
