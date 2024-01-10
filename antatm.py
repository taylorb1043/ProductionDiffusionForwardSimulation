#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Syntax: pressure = antatm(elevation)

Units: elevation in m; pressure in hPa. Accepts vector inputs.

This function converts elevation to atmospheric pressure according to a best-fit relationship for Antarctic stations from:

Radok, Allison, Wendler, 1996. Atmospheric pressure over the interior of Antarctica. Antarctic Science 0 p. 209. 
 
Written by Greg Balco -- UW Cosmogenic Nuclide Lab
balcs@u.washington.edu
First version, March, 2001
Checked October, 2005
Part of the CRONUS-Earth online calculators: 
     http://hess.ess.washington.edu/math

Copyright 2001-2007, University of Washington
All rights reserved
Developed in part with funding from the National Science Foundation.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License, version 2, as published by the Free Software Foundation (www.fsf.org).

Translated into Python and edited by: Taylor Bourikas
Contact: tbourika@purdue.edu
Last modified: 10.26.23
"""
import math
import numpy as np

out_1 = 989.1 * math.exp(z/(-7588))

#Make it a row vector
if np.size(out_1, 1) > np.size(out_1, 2):
    out = np.transpose(out_1)
else:
    out = out_1