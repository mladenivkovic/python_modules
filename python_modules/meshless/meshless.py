#!/usr/bin/env python3

#===========================================
#
# A module containing common routines for 
# the meshless effective area visualisiation
# with 2d datasets
#
#===========================================


from meshlessio import read_file
from kernels import W, kernels
from particles import find_index
import settings


