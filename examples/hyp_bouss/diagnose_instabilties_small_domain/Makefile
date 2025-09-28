# Makefile for Clawpack code in this directory.
# This version only sets the TEMP files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(TEMPCLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = geoclaw                  # Clawpack package to use
EXE = xgeoclaw                 # Executable to create
SETRUN_FILE = setrun.py        # File containing function to make data
OUTDIR = _output               # Directory for output
SETPLOT_FILE = setplot.py      # File containing function to set plots
PLOTDIR = _plots               # Directory for plots

# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
FFLAGS ?= 

# ---------------------------------
# package sources for this program:
# ---------------------------------

# GEOLIB = $(CLAW)/geoclaw/src/2d/shallow
# include $(GEOLIB)/Makefile.geoclaw


# Add Boussinesq Makefile
include $(TEMPBOUSSLIB2D)/Makefile.geoclaw_bouss

# ---------------------------------------
# package sources specifically to exclude
# (i.e. if a custom replacement source 
#  under a different name is provided)
# ---------------------------------------

EXCLUDE_MODULES = \

EXCLUDE_SOURCES = \

# ----------------------------------------
# List of custom sources for this program:
# ----------------------------------------


MODULES = \
	lu_solver.f90\
	bouss_module.f90 \

SOURCES = \
	qinit.f90\
	b4step2.f90\
	setaux.f90\
	filval.f90\
	
#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

# Construct the topography data
.PHONY: topo all
topo:
	python maketopo.py
bouss_regions:
	python makeboussregions.py

all: 
	$(MAKE) topo
	$(MAKE) bouss_regions
	$(MAKE) .plots
	$(MAKE) .htmls

