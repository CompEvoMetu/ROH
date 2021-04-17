#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

f = open(sys.argv[1])
pstart = 0
for i in f:
    line=i.rstrip().split()
    pend = line[1]
    print("1", pstart, pend, "_".join(map(str, [1, pstart, pend])))
    print("1", line[1], line[2], "_".join(map(str, [1,line[1], line[2]])))
    pstart = line[2]
if int(pstart) < 249250621:
    pend=249250621
    print("1", pstart, pend, "_".join(map(str, [1, pstart, pend])))

f.close()
