#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse
import sys
import os
from numpy import zeros, genfromtxt, savetxt


def schVector(sch, iniDayWeek):
    schMultiplier = zeros(8904)
    it = 0
    while it < 8904:
        if iniDayWeek < 6:
            for d in range(6 - iniDayWeek):
                for h in range(int(sch[1, 0])):
                    schMultiplier[it] = sch[1, 1]
                    it = it + 1
                for hstp in range(int(sch[0, 0]) - 1):
                    for h in range(int(sch[hstp + 2, 0]) - int(sch[hstp + 1, 0])):
                        schMultiplier[it] = sch[hstp + 2, 1]
                        it = it + 1
            for d in range(2):
                for h in range(int(sch[int(1 + sch[0, 0]), 0])):
                    schMultiplier[it] = sch[int(1 + sch[0, 0]), 1]
                    it = it + 1
                for hstp in range(int(sch[0, 1]) - 1):
                    for h in range(
                        int(sch[int(sch[0, 0] + hstp + 2), 0])
                        - int(sch[int(sch[0, 0] + hstp + 1), 0])
                    ):
                        schMultiplier[it] = sch[int(sch[0, 0] + hstp + 2), 1]
                        it = it + 1
            for d in range(iniDayWeek - 1):
                for h in range(int(sch[1, 0])):
                    schMultiplier[it] = sch[1, 1]
                    it = it + 1
                for hstp in range(int(sch[0, 0]) - 1):
                    for h in range(int(sch[int(hstp + 2), 0]) - int(sch[int(hstp + 1), 0])):
                        schMultiplier[it] = sch[int(hstp + 2), 1]
                        it = it + 1
        else:
            for d in range(8 - iniDayWeek):
                for h in range(int(sch[int(1 + sch[0, 0]), 0])):
                    schMultiplier[it] = sch[int(1 + sch[0, 0]), 1]
                    it = it + 1
                for hstp in range(int(sch[0, 1]) - 1):
                    for h in range(
                        int(sch[int(sch[0, 0] + hstp + 2), 0])
                        - int(sch[int(sch[0, 0] + hstp + 1), 0])
                    ):
                        schMultiplier[it] = sch[int(sch[0, 0] + hstp + 2), 1]
                        it = it + 1
            for d in range(5):
                for h in range(int(sch[1, 0])):
                    schMultiplier[it] = sch[1, 1]
                    it = it + 1
                for hstp in range(int(sch[0, 0]) - 1):
                    for h in range(int(sch[int(hstp + 2), 0]) - int(sch[int(hstp + 1), 0])):
                        schMultiplier[it] = sch[int(hstp + 2), 1]
                        it = it + 1
            for d in range(iniDayWeek - 6):
                for h in range(int(sch[int(1 + sch[0, 0]), 0])):
                    schMultiplier[it] = sch[int(1 + sch[0, 0]), 1]
                    it = it + 1
                for hstp in range(int(sch[0, 1]) - 1):
                    for h in range(
                        int(sch[int(sch[0, 0] + hstp + 2), 0])
                        - int(sch[int(sch[0, 0] + hstp + 1), 0])
                    ):
                        schMultiplier[it] = sch[int(sch[0, 0] + hstp + 2), 1]
                        it = it + 1

    schMultiplierOut = zeros(8760)
    for i in range(8760):
        schMultiplierOut[i] = schMultiplier[i]
    return schMultiplierOut


def configSchedule(obj, num, k, kmax, line):
    myline = line.split(",")
    if num[0] < 1:
        for j in range(2):
            num[j] = float(myline[j])
        kmax = num[0] + num[1] + 1
        obj = zeros([kmax, 2])
        k = 0
    if num[0] > 0 and k < kmax:
        for j in range(2):
            obj[k, j] = float(myline[j])
        k = k + 1
    return obj, num, k, kmax


def file_exists(_path):
    _out = True
    try:
        open(_path, "r").close()
    except IOError:
        _out = False
    return _out


def main(compact_path, annual_path):
    sch = genfromtxt(compact_path, skip_header=1, delimiter=",")
    schMultiplier = schVector(sch, opts.w)
    savetxt(annual_path, schMultiplier, fmt="%1.1f")
    print("schedule generated")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="generates an annual hourly schedule file from a compact definition"
    )
    parser.add_argument("compact_file", help="compact defintion file")
    parser.add_argument("annual_file", help="annual schedule file")
    parser.add_argument(
        "-w",
        action="store",
        type=int,
        default=1,
        help="day of week for first of January {1-Monday, 2-Tuesday}",
    )

    opts = parser.parse_args()

    if opts.compact_file[0] == "/":
        compact_path = opts.compact_file
    else:
        compact_path = os.getcwd() + "/" + opts.compact_file

    if opts.annual_file[0] == "/":
        annual_path = opts.annual_file
    else:
        annual_path = os.getcwd() + "/" + opts.annual_file

    if not file_exists(compact_path):
        print("file not found", file=sys.stderr)
    else:
        main(compact_path, annual_path)
