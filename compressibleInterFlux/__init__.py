#!/usr/bin/env python

#--------------------------------------------------------------------------------------
## pythonFlu - Python wrapping for OpenFOAM C++ API
## Copyright (C) 2010- Alexey Petrov
## Copyright (C) 2009-2010 Pebble Bed Modular Reactor (Pty) Limited (PBMR)
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 
## See http://sourceforge.net/projects/pythonflu
##
## Author : Alexey PETROV
##


#----------------------------------------------------------------------------
from Foam import FOAM_VERSION
if FOAM_VERSION( "==", "010701" ):
    from compressibleInterFlux.r1_7_1 import *
    pass


#--------------------------------------------------------------------------------------
from Foam import FOAM_VERSION
if FOAM_VERSION( "==", "020000" ) or FOAM_VERSION( "==", "020001" ):
    from compressibleInterFlux.r2_0_0 import *
    pass


#--------------------------------------------------------------------------------------
from Foam import FOAM_VERSION
if FOAM_VERSION( ">=", "020100" ):
    from compressibleInterFlux.r2_1_0 import *
    pass


#--------------------------------------------------------------------------------------
def entry_point():
    try:
       engine = main_standalone
       pass
    except NameError:
       print
       print "There is no implementation of the current OpenFOAM version"
       print
       import os; os._exit( os.EX_OK )
       pass
    
    import sys; argv = sys.argv
    return engine( len( argv ), argv )

#--------------------------------------------------------------------------------------
if __name__ == "__main__" :
    entry_point()
    pass
    
    
#--------------------------------------------------------------------------------------
