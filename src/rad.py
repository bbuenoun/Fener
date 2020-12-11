#!/usr/bin/env python
# -*- coding: utf-8 -*-


import ioRoutines
import geometry
import multiprocessing


def classicBsdf(
    month,
    day,
    hour,
    irrdirn,
    irrdiffh,
    cores,
    ab,
    ad,
    ar,
    ass,
    aa,
    conIndex,
    numWin,
    lat,
    lon,
    tzone,
    workDir_path,
    output_path,
    bsdfVis_path,
):
    geometry.winRadBsdf(numWin, workDir_path, bsdfVis_path, conIndex)
    ioRoutines.shell(
        "gendaylit %i %i %f -a %1.2f -o %1.2f -m %1.2f -W %1.2f %1.2f > %sgendaylit.rad"
        % (month, day, hour, lat, lon, tzone, irrdirn, irrdiffh, workDir_path)
    )
    radWin = ""
    for i in range(numWin):
        radWin = radWin + "%swindowBSDF_%i.rad " % (workDir_path, i)
    ioRoutines.shell(
        "oconv -f %sgendaylit.rad %ssky_prev.rad %sroom.mat %sroom.rad %s > %sclassicBsdf.oct"
        % (workDir_path, workDir_path, workDir_path, workDir_path, radWin, workDir_path)
    )
    ioRoutines.shell(
        "cat %silluSensor.pts | rtrace -h -I -n %i -ab %i -ad %i -ar %i -as %i -aa %1.2f %sclassicBsdf.oct | rcalc -e '$1=47.4*$1+120*$2+11.6*$3' > %sillRadBsdf_%1.1f.dat"
        % (workDir_path, cores, ab, ad, ar, ass, aa, workDir_path, output_path, hour)
    )


def multiply(
    cpu_count,
    ill,
    numWin,
    numSensor,
    workDir_path,
    bsdfVis_path,
    conIndex,
    month,
    day,
    hour,
    irrdirn,
    irrdiffh,
    lat,
    lon,
    tzone,
):
    # illuminance calculation
    head = 9
    ioRoutines.shell(
        "gendaylit %i %i %f -a %1.2f -o %1.2f -m %1.2f -W %1.2f %1.2f | genskyvec -m 4 -c 1 1 1 > %stempVis.skv"
        % (month, day, hour, lat, lon, tzone, irrdirn, irrdiffh, workDir_path)
    )
    tasks = [
        "dctimestep %sill_0.vmx %s %swin_0.dmx %stempVis.skv | rcalc -e '$1=179*($1*0.265+$2*0.670+$3*0.065)'"
        % (workDir_path, bsdfVis_path[conIndex[i]], workDir_path, workDir_path)
        for i in range(numWin)
    ]
    with multiprocessing.Pool(cpu_count) as pool:
        outRaws = pool.map(ioRoutines.shell_output, tasks)
    for outRaw in outRaws:
        out = outRaw.split("\n")
        for j in range(numSensor):
            ill[j] = ill[j] + float(out[j + head])
    return ill
