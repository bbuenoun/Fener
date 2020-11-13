#!/usr/bin/env python
# -*- coding: utf-8 -*-


import ioRoutines
import multiprocessing
import os
from typing import Tuple

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

# TODO Replace ioRoutines.findPathName by os.path.split and os.path.splitext
# TODO Do not hardcode the separator "/"
def _directory_path_and_file_name(path: str) -> Tuple[str, str]:
    return ioRoutines.findPathName(path, "/", os.path.splitext(path)[1])


def _sky_path(weather_file_path: str, suffix: str, extension: str) -> str:
    return "{0}{1}{suffix}{extension}".format(
        *_directory_path_and_file_name(weather_file_path),
        suffix=suffix,
        extension=extension,
    )


def weather_sky_path(weather_file_path: str) -> str:
    return _sky_path(weather_file_path, suffix="", extension=".wea")


def sol_sky_path(weather_file_path: str) -> str:
    return _sky_path(weather_file_path, "Sol", ".smx")


def vis_sky_path(weather_file_path: str) -> str:
    return _sky_path(weather_file_path, "Vis", ".smx")


def sun_sky_path(weather_file_path: str) -> str:
    return _sky_path(weather_file_path, "Sun", ".smx")


def dif_sky_path(weather_file_path: str) -> str:
    return _sky_path(weather_file_path, "Dif", ".smx")


def sun_ang_sky_path(weather_file_path: str) -> str:
    return _sky_path(weather_file_path, "_sunAng", ".dat")


# ----------------------------------------------------
# sky matrix
def sky(cpu_count, weather_file_path, lat, lon, tzone, grndAlb, bipv=False):
    ioRoutines.shell(
        "epw2wea %s %s" % (weather_file_path, weather_sky_path(weather_file_path))
    )
    # ioRoutines.shell("gendaymtx -m 4 -O1 -d -h %s%s.wea > %s%sSun.smx" % (path,name,path,name))
    genday_tasks = [
        "gendaymtx -m 4 -O1 -h -g %1.2f %1.2f %1.2f %s > %s"
        % (
            grndAlb,
            grndAlb,
            grndAlb,
            weather_sky_path(weather_file_path),
            sol_sky_path(weather_file_path),
        ),
        "gendaymtx -m 4 -h -g %1.2f %1.2f %1.2f %s > %s"
        % (
            grndAlb,
            grndAlb,
            grndAlb,
            weather_sky_path(weather_file_path),
            vis_sky_path(weather_file_path),
        ),
    ]
    if bipv:
        genday_tasks.extend(
            [
                "gendaymtx -m 4 -O1 -d -h -g %1.2f %1.2f %1.2f %s > %s"
                % (
                    grndAlb,
                    grndAlb,
                    grndAlb,
                    weather_sky_path(weather_file_path),
                    sun_sky_path(weather_file_path),
                ),
                "gendaymtx -m 4 -O1 -s -h -g %1.2f %1.2f %1.2f %s > %s"
                % (
                    grndAlb,
                    grndAlb,
                    grndAlb,
                    weather_sky_path(weather_file_path),
                    dif_sky_path(weather_file_path),
                ),
            ]
        )
    # For an explanation on the differences between `map*` and `imap*` see
    # https://stackoverflow.com/questions/26520781/multiprocessing-pool-whats-the-difference-between-map-async-and-imap/26521507#26521507
    with multiprocessing.Pool(cpu_count) as pool:
        pool.map(ioRoutines.shell, genday_tasks)
    nobis = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    gensky_tasks = [
        "gensky %i %i %1.1f -a %f -o %f -m %f"
        % (mm + 1, dd + 1, hh + 0.5, lat, lon, tzone)
        for mm in range(12)
        for dd in range(nobis[mm])
        for hh in range(24)
    ]
    with multiprocessing.Pool(cpu_count) as pool:
        coordsun_temps = pool.map(ioRoutines.shell_output, gensky_tasks)
    with open(sun_ang_sky_path(weather_file_path), "w") as fw:
        fw.write("%s, %s\n" % ("Sun Altitude (deg)", "Sun Azimuth (deg)"))
        for coordsun_temp in coordsun_temps:
            fw.write(
                "%s, %s\n" % (coordsun_temp.split()[21], coordsun_temp.split()[22])
            )


# ----------------------------------------------------
# daylight coefficient matrix (NOT OPERATIVE)
def dc(opts):
    ioRoutines.shell(
        "rfluxmtx -n %i -I+ -ab 12 -ad 10000 -lw 1e-4 < points.txt - skies/sky.rad materials/testroom.mat objects/testroom.rad > dc.dmx"
        % (opts.c)
    )


# ----------------------------------------------------
# daylight matrices
def day(nCore, numWin, workDir):
    for i in range(numWin):
        ioRoutines.shell(
            "rfluxmtx -n %i -c 20000 -h %swindow_%i.rad %ssky.rad -w %sroom.mat %sroom.rad %soutside.rad> %swin_%i.dmx"
            % (nCore, workDir, i, workDir, workDir, workDir, workDir, workDir, i)
        )


# ----------------------------------------------------
# daylight matrices
def dayH(nCore, numWin, workDir):
    for i in range(numWin):
        ioRoutines.shell(
            "rfluxmtx -n %i -c 20000 %swindow_%i.rad %ssky.rad -w %sroom.mat %sroom.rad %soutside.rad> %swin_%i.dmx"
            % (nCore, workDir, i, workDir, workDir, workDir, workDir, workDir, i)
        )


# ----------------------------------------------------
# calculate view factors
def comp_vf(nCores, sender, receiver, workDir):
    return ioRoutines.shell(
        """rfluxmtx -n %i -c 20000 -h "!xform -I %s" %s -w %sroom.mat"""
        % (nCores, sender, receiver, workDir),
        outputFlag=True,
    )


# ----------------------------------------------------
# illuminance view matrix
def illView(nCore, numWin, workDir):
    for i in range(numWin):
        ioRoutines.shell(
            "rfluxmtx -n %i -o %sill_%i.vmx -I+ -ab 12 -ad 50000 -lw 2e-5 -h - %swindow_%i.rad %sroom.mat %sroom.rad < %silluSensor.pts"
            % (nCore, workDir, i, workDir, i, workDir, workDir, workDir)
        )


# ----------------------------------------------------
# illuminance view matrix
def illViewH(nCore, numWin, workDir):
    for i in range(numWin):
        ioRoutines.shell(
            "rfluxmtx -n %i -o %sill_%i.vmx -I+ -ab 12 -ad 50000 -lw 2e-5 - %swindow_%i.rad %sroom.mat %sroom.rad < %silluSensor.pts"
            % (nCore, workDir, i, workDir, i, workDir, workDir, workDir)
        )


# ----------------------------------------------------
# glare view matrix
def glareView(nCore, numWin, workDir):
    for i in range(numWin):
        ioRoutines.shell(
            "rfluxmtx -n %i -o %sillVert_%i.vmx -I+ -ab 12 -ad 50000 -lw 2e-5 -h - %swindow_%i.rad %sroom.mat %sroom.rad < %silluVertSensor.pts"
            % (nCore, workDir, i, workDir, i, workDir, workDir, workDir)
        )


def irrWinBackView(nCore, numWin, workDir):
    for i in range(numWin):
        ioRoutines.shell(
            "rfluxmtx -n %i -o %sirrWinBack_%i.vmx -I+ -ab 12 -ad 50000 -lw 2e-5 -h - %swindow_%i.rad %sroom.mat %sroom.rad < %sintWinBackSensor_%i.pts"
            % (nCore, workDir, i, workDir, i, workDir, workDir, workDir, i)
        )


# ----------------------------------------------------
# thermal view matrices
def irrView(nCore, numWin, numFrame, numSurf, workDir):
    for i in range(numWin):
        ioRoutines.shell(
            "rfluxmtx -n %i -o %sirrWinFront_%i.vmx -I+ -ab 12 -ad 50000 -lw 2e-5 -h - %swindow_%i.rad %sroom.mat %sroom.rad < %sintWinFrontSensor_%i.pts"
            % (nCore, workDir, i, workDir, i, workDir, workDir, workDir, i)
        )
        ioRoutines.shell(
            "rfluxmtx -n %i -o %sirrWinBack_%i.vmx -I+ -ab 12 -ad 50000 -lw 2e-5 -h - %swindow_%i.rad %sroom.mat %sroom.rad < %sintWinBackSensor_%i.pts"
            % (nCore, workDir, i, workDir, i, workDir, workDir, workDir, i)
        )
        for j in range(numSurf):
            ioRoutines.shell(
                "rfluxmtx -n %i -o %sirrSurfFront_%i_%i.vmx -I+ -ab 12 -ad 50000 -lw 2e-5 -h - %swindow_%i.rad %sroom.mat %sroom.rad < %sintSurfFrontSensor_%i.pts"
                % (nCore, workDir, i, j, workDir, i, workDir, workDir, workDir, j)
            )
            ioRoutines.shell(
                "rfluxmtx -n %i -o %sirrSurfBack_%i_%i.vmx -I+ -ab 12 -ad 50000 -lw 2e-5 -h - %swindow_%i.rad %sroom.mat %sroom.rad < %sintSurfBackSensor_%i.pts"
                % (nCore, workDir, i, j, workDir, i, workDir, workDir, workDir, j)
            )
        for j in range(numFrame):
            ioRoutines.shell(
                "rfluxmtx -n %i -o %sirrFrameFront_%i_%i.vmx -I+ -ab 12 -ad 50000 -lw 2e-5 -h - %swindow_%i.rad %sroom.mat %sroom.rad < %sintFrameFrontSensor_%i.pts"
                % (nCore, workDir, i, j, workDir, i, workDir, workDir, workDir, j)
            )
            ioRoutines.shell(
                "rfluxmtx -n %i -o %sirrFrameBack_%i_%i.vmx -I+ -ab 12 -ad 50000 -lw 2e-5 -h - %swindow_%i.rad %sroom.mat %sroom.rad < %sintFrameBackSensor_%i.pts"
                % (nCore, workDir, i, j, workDir, i, workDir, workDir, workDir, j)
            )
            ioRoutines.shell(
                "rfluxmtx -n %i -o %sirrRevFront_%i_%i.vmx -I+ -ab 12 -ad 50000 -lw 2e-5 -h - %swindow_%i.rad %sroom.mat %sroom.rad < %sintRevFrontSensor_%i.pts"
                % (nCore, workDir, i, j, workDir, i, workDir, workDir, workDir, j)
            )
            ioRoutines.shell(
                "rfluxmtx -n %i -o %sirrRevBack_%i_%i.vmx -I+ -ab 12 -ad 50000 -lw 2e-5 -h - %swindow_%i.rad %sroom.mat %sroom.rad < %sintRevBackSensor_%i.pts"
                % (nCore, workDir, i, j, workDir, i, workDir, workDir, workDir, j)
            )


# ----------------------------------------------------
# daylight matrix for outdoor calculations
def ext(nCore, numSurf, numFrame, bc, workDir):
    for i in range(numSurf):
        if bc[i]:  # exterior boundary conditions
            ioRoutines.shell(
                "rfluxmtx -n %i -o %sirrSurfExt_%i.dmx -I+ -ab 12 -ad 50000 -lw 2e-5 -h - %ssky.rad -w %sroom.mat %sroom.rad %soutside.rad < %sextSurfSensor_%i.pts"
                % (nCore, workDir, i, workDir, workDir, workDir, workDir, workDir, i)
            )
            # ioRoutines.shell("rfluxmtx -n %i -I+ -ab 12 -ad 50000 -lw 2e-5  < radiance/extSurfSensor_%i.pts - radiance/sky.rad radiance/room.mat radiance/room.rad > mtx/irrExt_%i.dmx" % (nCore,i,i), True)

    for i in range(numFrame):
        ioRoutines.shell(
            "rfluxmtx -n %i -o %sirrFrameExt_%i.dmx -I+ -ab 12 -ad 50000 -lw 2e-5 -h - %ssky.rad -w %sroom.mat %sroom.rad %soutside.rad < %sextFrameSensor_%i.pts"
            % (nCore, workDir, i, workDir, workDir, workDir, workDir, workDir, i)
        )


# ----------------------------------------------------
# daylight matrix for outdoor horizontal illuminance
def extIll(nCore, workDir):
    ioRoutines.shell(
        "rfluxmtx -n %i -o %sillExt.dmx -I+ -ab 12 -ad 50000 -lw 2e-5 -h - %ssky.rad -w %sroom.mat %sroom.rad %soutside.rad < %sextIlluSensor.pts"
        % (nCore, workDir, workDir, workDir, workDir, workDir, workDir)
    )


# ----------------------------------------------------
# daylight matrix for bipv calculations
def bipv(nCore, workDir):
    ioRoutines.shell(
        "rfluxmtx -n %i -o %sbipv.dmx -I+ -ab 12 -ad 50000 -lw 2e-5 -h - %ssky.rad -w %sroom.mat %sroom.rad %soutside.rad < %sbipvSensor.pts"
        % (nCore, workDir, workDir, workDir, workDir, workDir, workDir)
    )
