#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 26 18:06:59 2017

@author: jolyon
"""

from mpmath import mp
import time
from math import sin, pi

mp.dps = 100

start = time.time()

for _ in range(1000000):
    mp.sin(mp.pi)
 #   sin(pi)

end = time.time()

print("Finished in", round(end - start, 4), "s")
