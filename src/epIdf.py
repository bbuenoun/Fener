#!/usr/bin/env python
# -*- coding: utf-8 -*-


from numpy import mod
import os
import ioRoutines
import tempfile


def matrix2D(fw, name, mat, length1, length2):
    length = length1 * length2
    fw.write("Matrix:TwoDimension, %s, " % name)
    fw.write("%1.0f, " % length1)
    fw.write("%1.0f, " % length2)
    for i in range(length1):
        for j in range(length2):
            k = i * length2 + j
            if mod(k, 40) == 0:
                fw.write("\n    ")
            if k < length - 1:
                fw.write("%1.5f, " % mat[i, j])
            else:
                fw.write("%1.5f;\n" % mat[i, j])
                fw.write("\n")


def southSurface(fw, name, bc, length, width, height, con):
    fw.write("BuildingSurface:Detailed,\n")
    fw.write("    %s,       !- Name\n" % name)
    fw.write("    WALL,           !- Surface Type\n")
    fw.write("    conOpaq_%i,       !- Construction Name\n" % con)
    fw.write("    Zone,           !- Zone Name\n")
    if bc == 1:
        fw.write("    Outdoors,       !- Outside Boundary Condition\n")
        fw.write("    ,               !- Outside Boundary Condition Object\n")
        fw.write("    SunExposed,     !- Sun Exposure\n")
        fw.write("    WindExposed,    !- Wind Exposure\n")
    elif bc == 0:
        fw.write("    Surface,       !- Outside Boundary Condition\n")
        fw.write("    %s,             !- Outside Boundary Condition Object\n" % name)
        fw.write("    NoSun,     !- Sun Exposure\n")
        fw.write("    NoWind,    !- Wind Exposure\n")
    fw.write("    0.5,            !- View Factor to Ground\n")
    fw.write("    4,              !- Number of Vertices\n")
    fw.write("    %1.3f,%1.3f,%1.3f, \n" % (0.000, 0.000, height))
    fw.write("    %1.3f,%1.3f,%1.3f, \n" % (0.000, 0.000, 0.000))
    fw.write("    %1.3f,%1.3f,%1.3f, \n" % (length, 0.000, 0.000))
    fw.write("    %1.3f,%1.3f,%1.3f; \n" % (length, 0.000, height))
    fw.write("\n")


def eastSurface(fw, name, bc, length, width, height, con):
    fw.write("BuildingSurface:Detailed,\n")
    fw.write("    %s,       !- Name\n" % name)
    fw.write("    WALL,           !- Surface Type\n")
    fw.write("    conOpaq_%i,       !- Construction Name\n" % con)
    fw.write("    Zone,           !- Zone Name\n")
    if bc == 1:
        fw.write("    Outdoors,       !- Outside Boundary Condition\n")
        fw.write("    ,               !- Outside Boundary Condition Object\n")
        fw.write("    SunExposed,     !- Sun Exposure\n")
        fw.write("    WindExposed,    !- Wind Exposure\n")
    elif bc == 0:
        fw.write("    Surface,       !- Outside Boundary Condition\n")
        fw.write("    %s,             !- Outside Boundary Condition Object\n" % name)
        fw.write("    NoSun,     !- Sun Exposure\n")
        fw.write("    NoWind,    !- Wind Exposure\n")
    fw.write("    0.5,            !- View Factor to Ground\n")
    fw.write("    4,              !- Number of Vertices\n")
    fw.write("    %1.3f,%1.3f,%1.3f, \n" % (length, 0.000, height))
    fw.write("    %1.3f,%1.3f,%1.3f, \n" % (length, 0.000, 0.000))
    fw.write("    %1.3f,%1.3f,%1.3f, \n" % (length, width, 0.000))
    fw.write("    %1.3f,%1.3f,%1.3f; \n" % (length, width, height))
    fw.write("\n")


def northSurface(fw, name, bc, length, width, height, con):
    fw.write("BuildingSurface:Detailed,\n")
    fw.write("    %s,       !- Name\n" % name)
    fw.write("    WALL,           !- Surface Type\n")
    fw.write("    conOpaq_%i,       !- Construction Name\n" % con)
    fw.write("    Zone,           !- Zone Name\n")
    if bc == 1:
        fw.write("    Outdoors,       !- Outside Boundary Condition\n")
        fw.write("    ,               !- Outside Boundary Condition Object\n")
        fw.write("    SunExposed,     !- Sun Exposure\n")
        fw.write("    WindExposed,    !- Wind Exposure\n")
    elif bc == 0:
        fw.write("    Surface,       !- Outside Boundary Condition\n")
        fw.write("    %s,             !- Outside Boundary Condition Object\n" % name)
        fw.write("    NoSun,     !- Sun Exposure\n")
        fw.write("    NoWind,    !- Wind Exposure\n")
    fw.write("    0.5,            !- View Factor to Ground\n")
    fw.write("    4,              !- Number of Vertices\n")
    fw.write("    %1.3f,%1.3f,%1.3f, \n" % (length, width, height))
    fw.write("    %1.3f,%1.3f,%1.3f, \n" % (length, width, 0.000))
    fw.write("    %1.3f,%1.3f,%1.3f, \n" % (0.000, width, 0.000))
    fw.write("    %1.3f,%1.3f,%1.3f; \n" % (0.000, width, height))
    fw.write("\n")


def westSurface(fw, name, bc, length, width, height, con):
    fw.write("BuildingSurface:Detailed,\n")
    fw.write("    %s,       !- Name\n" % name)
    fw.write("    WALL,           !- Surface Type\n")
    fw.write("    conOpaq_%i,       !- Construction Name\n" % con)
    fw.write("    Zone,           !- Zone Name\n")
    if bc == 1:
        fw.write("    Outdoors,       !- Outside Boundary Condition\n")
        fw.write("    ,               !- Outside Boundary Condition Object\n")
        fw.write("    SunExposed,     !- Sun Exposure\n")
        fw.write("    WindExposed,    !- Wind Exposure\n")
    elif bc == 0:
        fw.write("    Surface,       !- Outside Boundary Condition\n")
        fw.write("    %s,             !- Outside Boundary Condition Object\n" % name)
        fw.write("    NoSun,     !- Sun Exposure\n")
        fw.write("    NoWind,    !- Wind Exposure\n")
    fw.write("    0.5,            !- View Factor to Ground\n")
    fw.write("    4,              !- Number of Vertices\n")
    fw.write("    %1.3f,%1.3f,%1.3f, \n" % (0.000, width, height))
    fw.write("    %1.3f,%1.3f,%1.3f, \n" % (0.000, width, 0.000))
    fw.write("    %1.3f,%1.3f,%1.3f, \n" % (0.000, 0.000, 0.000))
    fw.write("    %1.3f,%1.3f,%1.3f; \n" % (0.000, 0.000, height))
    fw.write("\n")


def ceilingSurface(fw, name, bc, length, width, height, con):
    fw.write("BuildingSurface:Detailed,\n")
    fw.write("    %s,       !- Name\n" % name)
    fw.write("    CEILING,           !- Surface Type\n")
    fw.write("    conOpaq_%i,       !- Construction Name\n" % con)
    fw.write("    Zone,           !- Zone Name\n")
    if bc == 1:
        fw.write("    Outdoors,       !- Outside Boundary Condition\n")
        fw.write("    ,               !- Outside Boundary Condition Object\n")
        fw.write("    SunExposed,     !- Sun Exposure\n")
        fw.write("    WindExposed,    !- Wind Exposure\n")
    elif bc == 0:
        fw.write("    Surface,       !- Outside Boundary Condition\n")
        fw.write("    %s,             !- Outside Boundary Condition Object\n" % name)
        fw.write("    NoSun,     !- Sun Exposure\n")
        fw.write("    NoWind,    !- Wind Exposure\n")
    fw.write("    0.5,            !- View Factor to Ground\n")
    fw.write("    4,              !- Number of Vertices\n")
    fw.write("    %1.3f,%1.3f,%1.3f, \n" % (0.000, width, height))
    fw.write("    %1.3f,%1.3f,%1.3f, \n" % (0.000, 0.000, height))
    fw.write("    %1.3f,%1.3f,%1.3f, \n" % (length, 0.000, height))
    fw.write("    %1.3f,%1.3f,%1.3f; \n" % (length, width, height))
    fw.write("\n")


def floorSurface(fw, name, bc, length, width, height, con):
    fw.write("BuildingSurface:Detailed,\n")
    fw.write("    %s,       !- Name\n" % name)
    fw.write("    FLOOR,           !- Surface Type\n")
    fw.write("    conOpaq_%i,       !- Construction Name\n" % con)
    fw.write("    Zone,           !- Zone Name\n")
    if bc == 1:
        fw.write("    Outdoors,       !- Outside Boundary Condition\n")
        fw.write("    ,               !- Outside Boundary Condition Object\n")
        fw.write("    SunExposed,     !- Sun Exposure\n")
        fw.write("    WindExposed,    !- Wind Exposure\n")
    elif bc == 0:
        fw.write("    Surface,       !- Outside Boundary Condition\n")
        fw.write("    %s,             !- Outside Boundary Condition Object\n" % name)
        fw.write("    NoSun,     !- Sun Exposure\n")
        fw.write("    NoWind,    !- Wind Exposure\n")
    fw.write("    0.5,            !- View Factor to Ground\n")
    fw.write("    4,              !- Number of Vertices\n")
    fw.write("    %1.3f,%1.3f,%1.3f, \n" % (0.000, 0.000, 0.000))
    fw.write("    %1.3f,%1.3f,%1.3f, \n" % (0.000, width, 0.000))
    fw.write("    %1.3f,%1.3f,%1.3f, \n" % (length, width, 0.000))
    fw.write("    %1.3f,%1.3f,%1.3f; \n" % (length, 0.000, 0.000))
    fw.write("\n")


def internalMass(fw, surf):
    fw.write("InternalMass,\n")
    fw.write("    InternalMass,        !- Surface Name\n")
    fw.write("    Mass,                !- Construction Name\n")
    fw.write("    Zone,                !- Zone\n")
    fw.write("    %1.2f;           !- Total area exposed to Zone {m2}\n" % surf)
    fw.write("\n")


def shadeSurf(fw, name, shdSurf):
    fw.write("Shading:Site,\n")
    fw.write("    %s,     !- Name\n" % name)
    fw.write("    %1.2f,      !- Azimuth Angle {deg}\n" % shdSurf[0])
    fw.write("    90,        !- Tilt Angle {deg}\n")
    fw.write("    %1.2f,    !- Starting X Coordinate {m}\n" % shdSurf[1])
    fw.write("    %1.2f,    !- Starting Y Coordinate {m}\n" % shdSurf[2])
    fw.write("    %1.2f,      !- Starting Z Coordinate {m}\n" % shdSurf[3])
    fw.write("    %1.2f,     !- Length {m}\n" % shdSurf[4])
    fw.write("    %1.2f;      !- Height {m}\n" % shdSurf[5])
    fw.write("\n")


def schedule(fw, name, sch):
    fw.write("Schedule:Compact,\n")
    fw.write("    %s,     !- Name\n" % name)
    fw.write("    Any Number,         !- Schedule Type Limits Name\n")
    fw.write("    Through: 12/31,    !- Field 1\n")
    fw.write(
        "    For: WeekDays SummerDesignDay WinterDesignDay CustomDay1 CustomDay2 Holiday, !- Field 2\n"
    )
    for i in range(int(sch[0, 0])):
        fw.write("    Until: %1.0f:00,%1.4f, \n" % (sch[i + 1, 0], sch[i + 1, 1]))
    fw.write("    For: Weekends, !- Field 2\n")
    for i in range(int(sch[0, 1]) - 1):
        fw.write(
            "    Until: %1.0f:00,%1.4f, \n"
            % (sch[int(sch[0, 0]) + i + 1, 0], sch[int(sch[0, 0]) + i + 1, 1])
        )
    fw.write(
        "    Until: %1.0f:00,%1.4f; \n"
        % (sch[int(sch[0, 0] + sch[0, 1]), 0], sch[int(sch[0, 0] + sch[0, 1]), 1])
    )
    fw.write("\n")


def shadeCntrl(fw, j, cntrl, numWin):
    # print "WARNING: Shade control is currently not operative in EnergyPlus due to an unsolved bug of the program"
    # if int(cntrl[0]) == 0:
    #     surf = "SouthWall"
    # elif int(cntrl[0]) == 1:
    #     surf = "EastWall"
    # elif int(cntrl[0]) == 2:
    #     surf = "NorthWall"
    # elif int(cntrl[0]) == 3:
    #     surf = "WestWall"
    # ------------------------------------------------------------------------
    for i in range(int(cntrl[3])):
        fw.write("EnergyManagementSystem:ConstructionIndexVariable,\n")
        fw.write("    conWin_%i,    \n" % int(cntrl[i * 2 + 5]))
        fw.write("    conWin_%i;    \n" % int(cntrl[i * 2 + 5]))
        fw.write("\n")
    # ------------------------------------------------------------------------
    fw.write("EnergyManagementSystem:Sensor,\n")
    fw.write("    Sensor%i,       \n" % j)
    fw.write("    ShadeControl%i,             \n" % j)
    fw.write("    Schedule Value;             \n")
    fw.write("\n")
    for k in range(numWin[int(cntrl[0])]):
        # ------------------------------------------------------------------------
        fw.write("EnergyManagementSystem:ProgramCallingManager,\n")
        fw.write("    Manager%i%i,       \n" % (j, k))
        fw.write("    BeginTimestepBeforePredictor, \n")
        fw.write("    Control%i%i;       \n" % (j, k))
        fw.write("\n")
        # ------------------------------------------------------------------------
        fw.write("EnergyManagementSystem:Actuator,\n")
        fw.write("    Actuator%i%i,       \n" % (j, k))
        if int(cntrl[0]) == 0:
            fw.write("    SouthWindow%i,    \n" % k)
        if int(cntrl[0]) == 1:
            fw.write("    EastWindow%i,    \n" % k)
        if int(cntrl[0]) == 2:
            fw.write("    NorthWindow%i,    \n" % k)
        if int(cntrl[0]) == 3:
            fw.write("    WestWindow%i,    \n" % k)
        if int(cntrl[0]) == 4:
            fw.write("    CeilingWindow%i,    \n" % k)
        fw.write("    Surface,          \n")
        fw.write("    Construction State; \n")
        fw.write("\n")
        # ------------------------------------------------------------------------
        fw.write("EnergyManagementSystem:Program,\n")
        fw.write("    Control%i%i,                 \n" % (j, k))
        fw.write("    IF Sensor%i < 1,     \n" % j)
        fw.write("    SET Actuator%i%i = ConWin0, \n" % (j, k))
        fw.write("    ELSEIF Sensor%i < 2,     \n" % j)
        fw.write("    SET Actuator%i%i = ConWin1, \n" % (j, k))
        fw.write("    ELSEIF Sensor%i < 3,     \n" % j)
        fw.write("    SET Actuator%i%i = ConWin2, \n" % (j, k))
        fw.write("    ELSEIF Sensor%i < 4,     \n" % j)
        fw.write("    SET Actuator%i%i = ConWin3, \n" % (j, k))
        fw.write("    ELSEIF Sensor%i < 5,     \n" % j)
        fw.write("    SET Actuator%i%i = ConWin4, \n" % (j, k))
        fw.write("    ELSEIF Sensor%i < 6,     \n" % j)
        fw.write("    SET Actuator%i%i = ConWin5, \n" % (j, k))
        fw.write("    ELSEIF Sensor%i < 7,     \n" % j)
        fw.write("    SET Actuator%i%i = ConWin6, \n" % (j, k))
        fw.write("    ELSE ,                     \n")
        fw.write("    SET Actuator%i%i = ConWin0, \n" % (j, k))
        fw.write("    ENDIF;                     \n")
        fw.write("\n")
    # ------------------------------------------------------------------------
    fw.write("Schedule:File,\n")
    fw.write("    ShadeControl%i    ,       \n" % j)
    fw.write("    Any Number,          \n")
    fw.write("    winConSch.csv,          \n")
    fw.write("    %i,          \n" % int(cntrl[0] + 2))
    fw.write("    1,          \n")
    fw.write("    8760,          \n")
    fw.write("    Comma;          \n")
    fw.write("\n")


# def shadeCntrl(fw,j,cntrl):
#    if int(cntrl[0])==0:
#        surf = 'SouthWall'
#    elif int(cntrl[0])==1:
#        surf = 'EastWall'
#    elif int(cntrl[0])==2:
#        surf = 'NorthWall'
#    elif int(cntrl[0])==3:
#        surf = 'WestWall'
#    #------------------------------------------------------------------------
#    if int(cntrl[1])==0:
#        sensor = 'Surface Outside Face Solar Radiation Heat Gain Rate per Area'
#        #print "WARNING: Shade control is currently not operative in EnergyPlus due to an unsolved bug of the program"
#    elif int(cntrl[1])==1:
#        sensor = 'Glare'
#        print "WARNING: No glare control implemented for EneryPlus simulation"
#    else:
#        print "ERROR: Only two control type implemented"
#    #------------------------------------------------------------------------
#    for i in range(int(cntrl[2])):
#            fw.write('EnergyManagementSystem:ConstructionIndexVariable,\n')
#            fw.write('    ConWin%i,    \n' % int(cntrl[i*2+4]))
#            fw.write('    ConWin%i;    \n' % int(cntrl[i*2+4]))
#                fw.write('\n')
#    #------------------------------------------------------------------------
#    fw.write('EnergyManagementSystem:Sensor,\n')
#    fw.write('    Sensor%i,       \n' % j)
#    fw.write('    %s,             \n' % surf)
#    fw.write('    %s;             \n' % sensor)
#    fw.write('\n')
#    #------------------------------------------------------------------------
#    fw.write('EnergyManagementSystem:Actuator,\n')
#    fw.write('    Actuator%i,       \n' % j)
#    fw.write('    SouthWindow%i,    \n' % j)
#    fw.write('    Surface,          \n')
#    fw.write('    Construction State; \n')
#    fw.write('\n')
#    #------------------------------------------------------------------------
#    fw.write('EnergyManagementSystem:ProgramCallingManager,\n')
#    fw.write('    Manager%i,       \n' % j)
#    fw.write('    BeginTimestepBeforePredictor, \n')
#    fw.write('    Control%i;       \n' % j)
#    fw.write('\n')
#    #------------------------------------------------------------------------
#    fw.write('EnergyManagementSystem:Program,\n')
#    fw.write('    Control%i,                 \n' % j)
#    fw.write('    IF Sensor%i <= %1.2f,     \n' % (j,cntrl[5]))
#    fw.write('    SET Actuator%i = ConWin%i, \n' % (j,cntrl[4]))
#    for i in range(int(cntrl[2])-2):
#            fw.write('    ELSEIF Sensor%i <= %1.2f,     \n' % (j,cntrl[i*2+7]))
#            fw.write('    SET Actuator%i = ConWin%i, \n' % (j,cntrl[i*2+6]))
#    fw.write('    ELSE ,                     \n')
#    fw.write('    SET Actuator%i = ConWin%i, \n' % (j,cntrl[int(cntrl[2])+4]))
#    fw.write('    ENDIF;                     \n')
#    fw.write('\n')


def southWindow(fw, name, win, frame, thickFloor, thickWest):
    fw.write("FenestrationSurface:Detailed,\n")
    fw.write("    %s,   !- Name\n" % name)
    fw.write("    Window,         !- Surface Type\n")
    fw.write("    conWin_%i,   !- Construction Name\n" % int(win.con))
    fw.write("    SouthWall,      !- Building Surface Name\n")
    fw.write("    ,               !- Outside Boundary Condition Object\n")
    fw.write("    0.5,            !- View Factor to Ground\n")
    fw.write("    ,               !- Shading Control Name\n")
    fw.write("    %s,   !- Frame and Divider Name\n" % name)
    fw.write("    1,              !- Multiplier\n")
    fw.write("    4,              !- Number of Vertices\n")
    fw.write(
        "    %1.3f,%1.3f,%1.3f, \n"
        % (
            (frame.xOffset + win.xOffset) - thickWest,
            0.000,
            win.height + (frame.zOffset + win.zOffset) - thickFloor,
        )
    )
    fw.write(
        "    %1.3f,%1.3f,%1.3f, \n"
        % (
            (frame.xOffset + win.xOffset) - thickWest,
            0.000,
            (frame.zOffset + win.zOffset) - thickFloor,
        )
    )
    fw.write(
        "    %1.3f,%1.3f,%1.3f, \n"
        % (
            (frame.xOffset + win.xOffset) + win.length - thickWest,
            0.000,
            (frame.zOffset + win.zOffset) - thickFloor,
        )
    )
    fw.write(
        "    %1.3f,%1.3f,%1.3f; \n"
        % (
            (frame.xOffset + win.xOffset) + win.length - thickWest,
            0.000,
            win.height + (frame.zOffset + win.zOffset) - thickFloor,
        )
    )
    fw.write("\n")


def eastWindow(fw, name, win, frame, lenth, thickFloor, thickSouth):
    fw.write("FenestrationSurface:Detailed,\n")
    fw.write("    %s,   !- Name\n" % name)
    fw.write("    Window,         !- Surface Type\n")
    fw.write("    conWin_%i,   !- Construction Name\n" % int(win.con))
    fw.write("    EastWall,      !- Building Surface Name\n")
    fw.write("    ,               !- Outside Boundary Condition Object\n")
    fw.write("    0.5,            !- View Factor to Ground\n")
    fw.write("    ,               !- Shading Control Name\n")
    fw.write("    %s,   !- Frame and Divider Name\n" % name)
    fw.write("    1,              !- Multiplier\n")
    fw.write("    4,              !- Number of Vertices\n")
    fw.write(
        "    %1.3f,%1.3f,%1.3f, \n"
        % (
            lenth,
            (frame.xOffset + win.xOffset) - thickSouth,
            win.height + (frame.zOffset + win.zOffset) - thickFloor,
        )
    )
    fw.write(
        "    %1.3f,%1.3f,%1.3f, \n"
        % (
            lenth,
            (frame.xOffset + win.xOffset) - thickSouth,
            (frame.zOffset + win.zOffset) - thickFloor,
        )
    )
    fw.write(
        "    %1.3f,%1.3f,%1.3f, \n"
        % (
            lenth,
            (frame.xOffset + win.xOffset) + win.length - thickSouth,
            (frame.zOffset + win.zOffset) - thickFloor,
        )
    )
    fw.write(
        "    %1.3f,%1.3f,%1.3f; \n"
        % (
            lenth,
            (frame.xOffset + win.xOffset) + win.length - thickSouth,
            win.height + (frame.zOffset + win.zOffset) - thickFloor,
        )
    )
    fw.write("\n")


def northWindow(fw, name, win, frame, width, thickFloor, rightOffset):
    fw.write("FenestrationSurface:Detailed,\n")
    fw.write("    %s,   !- Name\n" % name)
    fw.write("    Window,         !- Surface Type\n")
    fw.write("    conWin_%i,   !- Construction Name\n" % int(win.con))
    fw.write("    NorthWall,      !- Building Surface Name\n")
    fw.write("    ,               !- Outside Boundary Condition Object\n")
    fw.write("    0.5,            !- View Factor to Ground\n")
    fw.write("    ,               !- Shading Control Name\n")
    fw.write("    %s,   !- Frame and Divider Name\n" % name)
    fw.write("    1,              !- Multiplier\n")
    fw.write("    4,              !- Number of Vertices\n")
    fw.write(
        "    %1.3f,%1.3f,%1.3f, \n"
        % (
            rightOffset + win.length,
            width,
            win.height + (frame.zOffset + win.zOffset) - thickFloor,
        )
    )
    fw.write(
        "    %1.3f,%1.3f,%1.3f, \n"
        % (rightOffset + win.length, width, (frame.zOffset + win.zOffset) - thickFloor)
    )
    fw.write(
        "    %1.3f,%1.3f,%1.3f, \n"
        % (rightOffset, width, (frame.zOffset + win.zOffset) - thickFloor)
    )
    fw.write(
        "    %1.3f,%1.3f,%1.3f; \n"
        % (rightOffset, width, win.height + (frame.zOffset + win.zOffset) - thickFloor)
    )
    fw.write("\n")


def westWindow(fw, name, win, frame, thickFloor, rightOffset):
    fw.write("FenestrationSurface:Detailed,\n")
    fw.write("    %s,   !- Name\n" % name)
    fw.write("    Window,         !- Surface Type\n")
    fw.write("    conWin_%i,   !- Construction Name\n" % int(win.con))
    fw.write("    WestWall,      !- Building Surface Name\n")
    fw.write("    ,               !- Outside Boundary Condition Object\n")
    fw.write("    0.5,            !- View Factor to Ground\n")
    fw.write("    ,               !- Shading Control Name\n")
    fw.write("    %s,   !- Frame and Divider Name\n" % name)
    fw.write("    1,              !- Multiplier\n")
    fw.write("    4,              !- Number of Vertices\n")
    fw.write(
        "    %1.3f,%1.3f,%1.3f, \n"
        % (
            0.000,
            rightOffset + win.length,
            win.height + (frame.zOffset + win.zOffset) - thickFloor,
        )
    )
    fw.write(
        "    %1.3f,%1.3f,%1.3f, \n"
        % (0.000, rightOffset + win.length, (frame.zOffset + win.zOffset) - thickFloor)
    )
    fw.write(
        "    %1.3f,%1.3f,%1.3f, \n"
        % (0.000, rightOffset, (frame.zOffset + win.zOffset) - thickFloor)
    )
    fw.write(
        "    %1.3f,%1.3f,%1.3f; \n"
        % (0.000, rightOffset, win.height + (frame.zOffset + win.zOffset) - thickFloor)
    )
    fw.write("\n")


def ceilingWindow(fw, name, win, frame, height, thickWest, thickSouth):
    fw.write("FenestrationSurface:Detailed,\n")
    fw.write("    %s,   !- Name\n" % name)
    fw.write("    Window,         !- Surface Type\n")
    fw.write("    conWin_%i,   !- Construction Name\n" % int(win.con))
    fw.write("    Ceiling,      !- Building Surface Name\n")
    fw.write("    ,               !- Outside Boundary Condition Object\n")
    fw.write("    0.0,            !- View Factor to Ground\n")
    fw.write("    ,               !- Shading Control Name\n")
    fw.write("    %s,   !- Frame and Divider Name\n" % name)
    fw.write("    1,              !- Multiplier\n")
    fw.write("    4,              !- Number of Vertices\n")
    fw.write(
        "    %1.3f,%1.3f,%1.3f, \n"
        % (
            (frame.xOffset + win.xOffset) - thickWest,
            (frame.zOffset + win.zOffset) + win.height - thickSouth,
            height,
        )
    )
    fw.write(
        "    %1.3f,%1.3f,%1.3f, \n"
        % (
            (frame.xOffset + win.xOffset) - thickWest,
            (frame.zOffset + win.zOffset) - thickSouth,
            height,
        )
    )
    fw.write(
        "    %1.3f,%1.3f,%1.3f, \n"
        % (
            (frame.xOffset + win.xOffset) + win.length - thickWest,
            (frame.zOffset + win.zOffset) - thickSouth,
            height,
        )
    )
    fw.write(
        "    %1.3f,%1.3f,%1.3f; \n"
        % (
            (frame.xOffset + win.xOffset) + win.length - thickWest,
            (frame.zOffset + win.zOffset) + win.height - thickSouth,
            height,
        )
    )
    fw.write("\n")


def addFrame(fw, name, frame, frameWidth, thick):
    # if win[5]>0:
    #   print "WARNING: Simulation of frames is currently not operative in EnergyPlus due to an unsolved bug of the program"
    fw.write("WindowProperty:FrameAndDivider,\n")
    fw.write("    %s,   ! Frame/Divider Name\n" % name)
    fw.write("    %1.3f,         ! Frame Width\n" % frameWidth)
    fw.write("    0.,         ! Frame Outside Projection\n")
    fw.write("    0.,         ! Frame Inside Projection\n")
    fw.write("    %1.3f,         ! Frame Conductance\n" % (1.0 / frame.res))
    fw.write(
        "    1.0,         ! Ratio of Frame-Edge Glass Conductance to Center-Of-Glass Conductance\n"
    )
    fw.write("    %1.3f,         ! Frame Solar Absorptance\n" % (1.0 - frame.extAlb))
    fw.write("    %1.3f,         ! Frame Visible Absorptance\n" % (1.0 - frame.extAlb))
    fw.write("    %1.3f,         ! Frame Thermal Emissivity\n" % frame.extEmis)
    fw.write("    DividedLite,    ! Divider Type\n")
    fw.write("    0.,         ! Divider Width\n")
    fw.write("    0.,         ! Number of Horizontal Dividers\n")
    fw.write("    0.,         ! Number of Vertical Dividers\n")
    fw.write("    0.,         ! Divider Outside Projection\n")
    fw.write("    0.,         ! Divider Inside Projection\n")
    fw.write("    0.,         ! Divider Conductance\n")
    fw.write(
        "    1.,         ! Ratio of Divider-Edge Glass Conductance to Center-Of-Glass Conductance\n"
    )
    fw.write("    0.5,         ! Divider Solar Absorptance\n")
    fw.write("    0.5,         ! Divider Visible Absorptance\n")
    fw.write("    0.8,         ! Divider Thermal Emissivity\n")
    fw.write("    0.5,         ! Outside Reveal Solar Absorptance\n")
    fw.write("    %1.3f,         ! Inside Sill Depth (m)\n" % (thick - frame.reveal))
    fw.write("    0.5,         ! Inside Sill Solar Absorptance\n")
    fw.write("    %1.3f,         ! Inside Reveal Depth (m)\n" % (thick - frame.reveal))
    fw.write("    0.5;         ! Inside Reveal Solar Absorptance\n")
    fw.write("\n")


def overhang(fw, name, depth):
    fw.write("Shading:Overhang,\n")
    fw.write("    Overhang%s, !- Name\n" % name)
    fw.write("    %s, !- Window or Door Name\n" % name)
    fw.write("    0, !- Height above Window or Door {m}\n")
    fw.write("    90, !- Tilt Angle from Window/Door {deg}\n")
    fw.write("    0, !- Left extension from Window/Door Width {m}\n")
    fw.write("    0, !- Right extension from Window/Door Width {m}\n")
    fw.write("    %1.3f; !- Depth {m}\n" % depth)
    fw.write("\n")


def fin(fw, name, depth):
    fw.write("Shading:Fin:Projection,\n")
    fw.write("    Fin%s, !- Name\n" % name)
    fw.write("    %s, !- Window or Door Name\n" % name)
    fw.write("    0, !- Left Extension from Window/Door {m}\n")
    fw.write("    0, !- Left Distance Above Top of Window {m}\n")
    fw.write("    0, !- Left Distance Below Bottom of Window {m}\n")
    fw.write("    90, !- Left Tilt Angle from Window/Door {deg}\n")
    fw.write("    %1.3f, !- Left Fin Depth {m}\n" % depth)
    fw.write("    0, !- Right Extension from Window/Door {m}\n")
    fw.write("    0, !- Right Distance Above Top of Window {m}\n")
    fw.write("    0, !- Right Distance Below Bottom of Window {m}\n")
    fw.write("    90, !- Right Tilt Angle from Window/Door {deg}\n")
    fw.write("    %1.3f; !- Right Fin Depth {m}\n" % depth)
    fw.write("\n")


def constructionWin(fw, name, win):
    fw.write("Construction:ComplexFenestrationState,\n")
    fw.write("    %s,            ! Name\n" % name)
    fw.write("    LBNLWindow,         !- basis type\n")
    fw.write("    None,               !- basis symmetry type\n")
    fw.write("    ThermParam,         !- window thermal model\n")
    fw.write("    CFS_Glz_Basis,      !- basis matrix name\n")
    fw.write("    Tfsol_%s,      !- Tfsol\n" % name)
    fw.write("    Rbsol_%s,      !- Rbsol\n" % name)
    fw.write("    Tfvis_%s,      !- Tfvis\n" % name)
    fw.write("    Tbvis_%s,      !- Tbvis\n" % name)
    if win[1] == 1:
        fw.write("    matGlz_%i,            !- layer 1 name\n" % win[2])
    elif win[1] == 2:
        fw.write("    MatShd%i,            !- layer 1 name\n" % win[2])
    elif win[1] == 3:
        fw.write("    MatBlnd%i,           !- layer 1 name\n" % win[2])
    elif win[1] == 4:
        fw.write("    matBSDF_%i,           !- layer 1 name\n" % win[2])
    fw.write("    fAbs_%s_Lay_0,  !- fAbs\n" % name)
    if win[0] == 1:
        fw.write("    bAbs_%s_Lay_0;  !- bAbs\n" % name)
    else:
        fw.write("    bAbs_%s_Lay_0,  !- bAbs\n" % name)
        fw.write("    GapMatGas_%i, \n" % win[4])
        fw.write("    , \n")
        fw.write("    , \n")
        if win[5] == 1:
            fw.write("    matGlz_%i,            \n" % win[6])
        elif win[5] == 2:
            fw.write("    MatShd%i,            \n" % win[6])
        elif win[5] == 3:
            fw.write("    MatBlnd%i,           \n" % win[6])
        elif win[5] == 4:
            fw.write("    matBSDF_%i,           \n" % win[6])
        fw.write("    fAbs_%s_Lay_1,  !- fAbs\n" % name)
        if win[0] == 2:
            fw.write("    bAbs_%s_Lay_1;  !- bAbs\n" % name)
        else:
            fw.write("    bAbs_%s_Lay_1,  !- bAbs\n" % name)
            fw.write("    GapMatGas_%i, \n" % win[8])
            fw.write("    , \n")
            fw.write("    , \n")
            if win[9] == 1:
                fw.write("    matGlz_%i,            \n" % win[10])
            elif win[9] == 2:
                fw.write("    MatShd%i,            \n" % win[10])
            elif win[9] == 3:
                fw.write("    MatBlnd%i,           \n" % win[10])
            elif win[9] == 4:
                fw.write("    matBSDF_%i,           \n" % win[10])
            fw.write("    fAbs_%s_Lay_2,  !- fAbs\n" % name)
            if win[0] == 3:
                fw.write("    bAbs_%s_Lay_2;  !- bAbs\n" % name)
            else:
                fw.write("    bAbs_%s_Lay_2,  !- bAbs\n" % name)
                fw.write("    GapMatGas_%i, \n" % win[12])
                fw.write("    , \n")
                fw.write("    , \n")
                if win[13] == 1:
                    fw.write("    matGlz_%i,            \n" % win[14])
                elif win[13] == 2:
                    fw.write("    MatShd%i,            \n" % win[14])
                elif win[13] == 3:
                    fw.write("    MatBlnd%i,           \n" % win[14])
                elif win[13] == 4:
                    fw.write("    matBSDF_%i,           \n" % win[14])
                fw.write("    fAbs_%s_Lay_3,  !- fAbs\n" % name)
                fw.write("    bAbs_%s_Lay_3;  !- bAbs\n" % name)

    fw.write("\n")


def material(fw, name, mat):
    fw.write("Material,\n")
    fw.write("    %s,     !- Name\n" % name)
    fw.write("    Rough,     !- Roughness\n")
    fw.write("    %1.4f, \n" % mat[0])
    fw.write("    %1.4f, \n" % mat[1])
    fw.write("    %1.4f, \n" % mat[2])
    fw.write("    %1.4f, \n" % mat[3])
    fw.write("    0.9,\n")
    fw.write("    0.5,\n")
    fw.write("    0.5;\n")
    fw.write("\n")


def construction(fw, name, con):
    fw.write("Construction,\n")
    fw.write("    %s,     !- Name\n" % name)
    for i in range(int(con[0]) - 1):
        fw.write("    mat_%i, \n" % con[i + 1])
    fw.write("    mat_%i; \n" % con[int(con[0])])
    fw.write("\n")


def materialGlz(fw, name, mat):
    fw.write("WindowMaterial:Glazing,\n")
    fw.write("    %s,     !- Name\n" % name)
    fw.write("    BSDF,     !- Optical Data Type\n")
    fw.write("    ,         \n")
    for i in range(11):
        fw.write("    %1.4f, \n" % mat[i])
    fw.write("    ,         \n")
    fw.write("    ,         \n")
    fw.write("    7.2e10,   !-Youngs modulus\n")
    fw.write("    0.22;     !-Poissons ratio\n")
    fw.write("\n")


def materialGas(fw, name, mat):
    if mat[1] > 0.999:
        fw.write("WindowMaterial:Gas,\n")
        fw.write("    %s,    !- Name        \n" % name)
        fw.write("    Air,   !- Gas Type       \n")
        fw.write("    %1.4f;     !- thickness\n" % mat[0])
    elif mat[2] > 0.999:
        fw.write("WindowMaterial:Gas,\n")
        fw.write("    %s,    !- Name        \n" % name)
        fw.write("    Argon,   !- Gas Type       \n")
        fw.write("    %1.4f;     !- thickness\n" % mat[0])
    else:
        fw.write("WindowMaterial:GasMixture,\n")
        fw.write("    %s,    !- Name\n" % name)
        fw.write("    %1.4f,     !- thickness\n" % mat[0])
        fw.write("    2,          !- Number of Gases in the mixture\n")
        fw.write("    Air,        !- Gas Name\n")
        fw.write("    %1.4f,     !- Fraction\n" % mat[1])
        fw.write("    Argon,      !- Gas Name\n")
        fw.write("    %1.4f;     !- Fraction\n" % mat[2])
    fw.write("\n")
    fw.write("WindowMaterial:Gap,\n")
    fw.write("    Gap%s,    !- gap name: Air (10) / Argon (90) Mix\n" % name)
    fw.write("    %1.4f,     !- thickness\n" % mat[0])
    fw.write("    %s,    !- Gas (or Gas Mixture) name\n" % name)
    fw.write("    101325.,    !- pressure\n")
    fw.write("    ;           !- deflection state\n")
    fw.write("\n")


def materialBSDF(fw, name, mat):
    fw.write("WindowMaterial:ComplexShade,\n")
    fw.write("    %s,    !- Name\n" % name)
    fw.write("    BSDF,       !- Optical Data Type\n")
    fw.write("    %1.4f,     !- thickness\n" % mat[0])
    fw.write("    %1.4f,     !- layer conductivity\n" % mat[1])
    fw.write("    0.0,        !- IR transmittance\n")
    fw.write("    %1.4f,     !- front emissivity\n" % mat[2])
    fw.write("    %1.4f,     !- back emissivity\n" % mat[3])
    fw.write("    %1.4f,     !- top opening multiplier\n" % mat[4])
    fw.write("    %1.4f,     !- bottom opening multiplier\n" % mat[5])
    fw.write("    %1.4f,     !- left side opening multiplier\n" % mat[6])
    fw.write("    %1.4f,     !- right side opening multiplier\n" % mat[7])
    fw.write("    %1.4f;     !- front opening multiplier\n" % mat[8])
    fw.write("\n")


def gen(
    ep_path,
    rotAng,
    volume,
    floorArea,
    iniMonth,
    iniDay,
    endMonth,
    endDay,
    iniDayWeek,
    airExch,
    powerLight,
    radFracLight,
    powerEquip,
    radFracEquip,
    powerPeople,
    radFracPeople,
    output_path,
    infSch_path,
    lightSch_path,
    equipSch_path,
    occSch_path,
    coolSpSch_path,
    heatSpSch_path,
    matOpaq,
    conOpaq,
    matGlz,
    matGas,
    matBSDF,
    conWinVect,
    conWin,
    bsdfSys_path,
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
    numWin,
    win,
    numFrame,
    frame,
    frameWidth,
):
    # ----------------------------------------------------
    ioRoutines.shell("cp -f %s %sinfiltration.dat" % (infSch_path, ep_path))
    ioRoutines.shell("cp -f %s %sequipment.dat" % (equipSch_path, ep_path))
    ioRoutines.shell("cp -f %s %soccupation.dat" % (occSch_path, ep_path))
    ioRoutines.shell("cp -f %s %scoolSetpoint.dat" % (coolSpSch_path, ep_path))
    ioRoutines.shell("cp -f %s %sheatSetpoint.dat" % (heatSpSch_path, ep_path))
    ioRoutines.shell("cp -f %slights.dat %slights.dat" % (output_path, ep_path))
    # ----------------------------------------------------
    # module creation by substitution
    with open("%swebport.idf" % ep_path, "w") as o:
        with open("template.idf") as template_file:
            for line in template_file:
                line = line.replace("@ORIENT", str(-rotAng))
                line = line.replace("@HEIGHT", str(volume / floorArea))
                line = line.replace("@VOL", str(volume))
                line = line.replace("@FLOORAREA", str(floorArea))
                line = line.replace("@ACH", str(airExch))
                line = line.replace("@BEGMONTH", str(iniMonth))
                line = line.replace("@BEGDAY", str(iniDay))
                line = line.replace("@ENDMONTH", str(endMonth))
                line = line.replace("@ENDDAY", str(endDay))
                line = line.replace("@LIGHTLEVEL", str(powerLight))
                line = line.replace("@LIGHTRAD", str(radFracLight))
                line = line.replace("@EQUIPLEVEL", str(powerEquip))
                line = line.replace("@EQUIPRAD", str(radFracEquip))
                line = line.replace("@NUMPEOPLE", str(powerPeople / 0.531 / 120.0))
                line = line.replace("@PEOPLERAD", str(radFracPeople))
                line = line.replace("@MASSAREA", str(0.0))
                o.write(line)
    # ----------------------------------------------------
    # module creation by object inception
    with open("%swebport.idf" % ep_path, "a") as fw:
        # Opaque surfaces
        fw.write("\n")
        fw.write("!================================\n")
        fw.write("!         Opaque Surfaces        \n")
        fw.write("!================================\n")
        fw.write("\n")
        southSurface(
            fw,
            "SouthWall",
            bcSouth,
            length - thickWest - thickEast,
            width - thickNorth - thickSouth,
            height - thickCeiling - thickFloor,
            conSouth,
        )
        eastSurface(
            fw,
            "EastWall",
            bcEast,
            length - thickWest - thickEast,
            width - thickNorth - thickSouth,
            height - thickCeiling - thickFloor,
            conEast,
        )
        northSurface(
            fw,
            "NorthWall",
            bcNorth,
            length - thickWest - thickEast,
            width - thickNorth - thickSouth,
            height - thickCeiling - thickFloor,
            conNorth,
        )
        westSurface(
            fw,
            "WestWall",
            bcWest,
            length - thickWest - thickEast,
            width - thickNorth - thickSouth,
            height - thickCeiling - thickFloor,
            conWest,
        )
        ceilingSurface(
            fw,
            "Ceiling",
            bcCeiling,
            length - thickWest - thickEast,
            width - thickNorth - thickSouth,
            height - thickCeiling - thickFloor,
            conCeiling,
        )
        floorSurface(
            fw,
            "Floor",
            bcFloor,
            length - thickWest - thickEast,
            width - thickNorth - thickSouth,
            height - thickCeiling - thickFloor,
            conFloor,
        )
        # Windows
        fw.write("\n")
        fw.write("!================================\n")
        fw.write("!         Windows        \n")
        fw.write("!================================\n")
        fw.write("\n")
        for i in range(numWin):
            for j in range(numFrame):
                if win[i].frame == j and frame[j].surf == 1:
                    southWindow(
                        fw,
                        "SouthWindow_%i" % i,
                        win[i],
                        frame[j],
                        thickFloor,
                        thickWest,
                    )
                    addFrame(fw, "SouthWindow_%i" % i, frame[j], frameWidth, thickSouth)
                if win[i].frame == j and frame[j].surf == 2:
                    eastWindow(
                        fw,
                        "EastWindow_%i" % i,
                        win[i],
                        frame[j],
                        length - thickWest - thickEast,
                        thickFloor,
                        thickSouth,
                    )
                    addFrame(fw, "SouthWindow_%i" % i, frame[j], frameWidth, thickEast)
                if win[i].frame == j and frame[j].surf == 3:
                    northWindow(
                        fw,
                        "NorthWindow_%i" % i,
                        win[i],
                        frame[j],
                        width - thickNorth - thickSouth,
                        thickFloor,
                        length - thickWest - win[i].xOffset - win[i].length,
                    )
                    addFrame(fw, "SouthWindow_%i" % i, frame[j], frameWidth, thickNorth)
                if win[i].frame == j and frame[j].surf == 4:
                    westWindow(
                        fw,
                        "WestWindow_%i" % i,
                        win[i],
                        frame[j],
                        thickFloor,
                        width - thickSouth - win[i].xOffset - win[i].length,
                    )
                    addFrame(fw, "SouthWindow_%i" % i, frame[j], frameWidth, thickWest)
                if win[i].frame == j and frame[j].surf == 5:
                    ceilingWindow(
                        fw,
                        "CeilingWindow_%i" % i,
                        win[i],
                        frame[j],
                        height - thickCeiling - thickFloor,
                        thickWest,
                        thickSouth,
                    )
                    addFrame(
                        fw, "SouthWindow_%i" % i, frame[j], frameWidth, thickCeiling
                    )
        # Materials
        fw.write("\n")
        fw.write("!================================\n")
        fw.write("!         Materials        \n")
        fw.write("!================================\n")
        fw.write("\n")
        for i in range(matOpaq.shape[0]):
            material(fw, "mat_%i" % i, matOpaq[i, :])
        # Glazing material
        for i in range(matGlz.shape[0]):
            materialGlz(fw, "matGlz_%i" % i, matGlz[i, :])
        # Gap material
        for i in range(matGas.shape[0]):
            materialGas(fw, "matGas_%i" % i, matGas[i, :])
        # BSDF material
        for i in range(matBSDF.shape[0]):
            materialBSDF(fw, "matBSDF_%i" % i, matBSDF[i, :])
        # Opaque constructions
        fw.write("\n")
        fw.write("!================================\n")
        fw.write("!         Constructions        \n")
        fw.write("!================================\n")
        fw.write("\n")
        for i in range(conOpaq.shape[0]):
            construction(fw, "conOpaq_%i" % i, conOpaq[i, :])
        # Window construction
        for i in range(conWinVect.shape[0]):
            constructionWin(fw, "conWin_%i" % i, conWinVect[i, :])
            TransMtxV = ioRoutines.xml(bsdfSys_path[i], "Visible", "Transmission Front")
            TransMtxS = ioRoutines.xml(bsdfSys_path[i], "Solar", "Transmission Front")
            ReflBackMtxS = ioRoutines.xml(bsdfSys_path[i], "Solar", "Reflection Back")
            matrix2D(fw, "Tfsol_conWin_%i" % i, TransMtxS, 145, 145)
            matrix2D(fw, "Rbsol_conWin_%i" % i, ReflBackMtxS, 145, 145)
            matrix2D(fw, "Tfvis_conWin_%i" % i, TransMtxV, 145, 145)
            matrix2D(fw, "Tbvis_conWin_%i" % i, TransMtxV, 145, 145)
            for k in range(conWin[i].numPane):
                matrix2D(
                    fw,
                    "fAbs_conWin_%i_Lay_%i" % (i, k),
                    conWin[i].absFront[k, :, :],
                    1,
                    145,
                )
                matrix2D(
                    fw,
                    "bAbs_conWin_%i_Lay_%i" % (i, k),
                    conWin[i].absBack[k, :, :],
                    1,
                    145,
                )


def run(ep_path, meteo_path):
    weatherFile = os.path.join(os.getcwd(), meteo_path)
    with tempfile.TemporaryDirectory() as temp_directory_path:
        ioRoutines.shell(
            "cp -rf %s %s" % (os.path.join(ep_path, "*"), temp_directory_path)
        )
        ioRoutines.shell(
            "cp -f %s %s"
            % (weatherFile, os.path.join(temp_directory_path, "weather.epw"))
        )
        ioRoutines.shell(
            "runenergyplus %s %s"
            % (
                os.path.join(temp_directory_path, "webport.idf"),
                os.path.join(temp_directory_path, "weather.epw"),
            )
        )
        ioRoutines.shell(
            "cp -rf %s %s" % (os.path.join(temp_directory_path, "Output"), ep_path)
        )
    with open("%sOutput/webport.err" % ep_path) as f:
        for line in f:
            print(line)
