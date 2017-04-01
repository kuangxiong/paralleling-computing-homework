#!/usr/bin/python
# -*- coding: UTF-8 -*-
from numpy import *
#import numpy as np
a = [[7,1,1,1],[1,9,1,1],[1,1,11,1], [1, 1, 1,13]]
a = array(a)
b = [0, 1,2,3]
b = array(b)
#solve linear equation ax=b
x = linalg.solve(a, b)
print x
