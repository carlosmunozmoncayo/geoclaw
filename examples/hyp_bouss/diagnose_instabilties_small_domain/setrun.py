"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""

from __future__ import absolute_import
from __future__ import print_function
import os
import numpy as np

try:
    CLAW = os.environ['CLAW']
except:
    raise Exception("*** Must first set CLAW enviornment variable")

# Scratch directory for storing topo and dtopo files:
scratch_dir = os.path.join(CLAW, 'geoclaw', 'scratch')

#System to be solved
Bouss_Hyp_rel = -1
if Bouss_Hyp_rel==-1:
    numeqn = 5
    numwaves = 5
elif Bouss_Hyp_rel==-2:
    numeqn = 7
    numwaves = 7
elif Bouss_Hyp_rel==0:
    numeqn = 5
    numwaves = 5



#------------------------------
def setrun(claw_pkg='geoclaw'):
#------------------------------

    """
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "geoclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData

    """

    from clawpack.clawutil import data

    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)


    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------
    
    #probdata = rundata.new_UserData(name='probdata',fname='setprob.data')


    #------------------------------------------------------------------
    # GeoClaw specific parameters:
    #------------------------------------------------------------------
    rundata = setgeo(rundata)

    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #   (or to amr2ez.data for AMR)
    #------------------------------------------------------------------
    clawdata = rundata.clawdata  # initialized when rundata instantiated


    # Set single grid parameters first.
    # See below for AMR parameters.


    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim

    # Lower and upper edge of computational domain:
    clawdata.lower[0] = -92.0#-91.5#-120.0      # west longitude
    clawdata.upper[0] = -90.0#-90.0#-60.0       # east longitude

    clawdata.lower[1] = -2.0#-60.0       # south latitude
    clawdata.upper[1] = 0.0         # north latitude



    # Number of grid cells: Coarsest grid
    clawdata.num_cells[0] = 40#500
    clawdata.num_cells[1] = 40#500

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = numeqn

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 5#First 3 taken

    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 2

    
    
    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0


    # Restart from checkpoint file of a previous run?
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in 
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False              # True to restart from prior results
    clawdata.restart_file = 'fort.chk00096'  # File to use for restart data

    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

    clawdata.output_style = 1

    if clawdata.output_style==1:
        # Output nout frames at equally spaced times up to tfinal:
        clawdata.num_output_times = 180*4
        clawdata.tfinal = 10000.*4.
        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list of output times.
        clawdata.output_times = [1.0]

    elif clawdata.output_style == 3:
        # Output every iout timesteps with a total of ntot time steps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 10
        # clawdata.output_step_interval = 4
        # clawdata.total_steps = 3077
        clawdata.output_t0 = True
        

    clawdata.output_format = 'ascii'      # 'ascii' or 'binary' 

    clawdata.output_q_components = 'all'   # need all
    clawdata.output_aux_components = 'all'  # eta=h+B is in q
    clawdata.output_aux_onlyonce = False    # output aux arrays each frame



    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 1



    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = True#True

    # Initial time step for variable dt.
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = 0.2

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99

    # Desired Courant number if variable dt used, and max to allow without
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.45
    clawdata.cfl_max = 0.5

    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 10000000

    #Setting aux array
    clawdata.num_aux = 5


    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2
    
    # Use dimensional splitting? (not yet available for AMR)
    clawdata.dimensional_split = 'unsplit'
    
    # For unsplit method, transverse_waves can be 
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 'all'

    # Number of waves in the Riemann solution:
    clawdata.num_waves = numwaves
    
    # List of limiters to use for each wave family:  
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'mc'       ==> MC limiter
    #   4 or 'vanleer'  ==> van Leer
    clawdata.limiter =clawdata.num_waves*['minmod']

    clawdata.use_fwaves = True    # True ==> use f-wave version of algorithms
    
    # Source terms splitting:
    #   src_split == 0 or 'none'    ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used, 
    #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = 'godunov'#'godunov'


    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2

    # Choice of BCs at xlower and xupper:
    #   0 => user specified (must modify bcN.f to use this option)
    #   1 => extrapolation (non-reflecting outflow)
    #   2 => periodic (must specify this at both boundaries)
    #   3 => solid wall for systems where q(2) is normal velocity

    clawdata.bc_lower[0] = 'extrap'
    clawdata.bc_upper[0] = 'extrap'

    clawdata.bc_lower[1] = 'extrap'
    clawdata.bc_upper[1] = 'extrap'



    # --------------
    # Checkpointing:
    # --------------

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    clawdata.checkpt_style = 0

    if clawdata.checkpt_style == 0:
        # Do not checkpoint at all
        pass

    elif np.abs(clawdata.checkpt_style) == 1:
        # Checkpoint only at tfinal.
        pass

    elif np.abs(clawdata.checkpt_style) == 2:
        # Specify a list of checkpoint times.  
        clawdata.checkpt_times = [0.1,0.15]

    elif np.abs(clawdata.checkpt_style) == 3:
        # Checkpoint every checkpt_interval timesteps (on Level 1)
        # and at the final time.
        clawdata.checkpt_interval = 5


    # ---------------
    # AMR parameters:
    # ---------------
    amrdata = rundata.amrdata

    # max number of refinement levels:
    amrdata.amr_levels_max = 1

    #Forcing a large number of cells on a single grid (just for debugging, comment out)
    amrdata.max1d = 1000# amrdata.memsize = 2*1048575

    # List of refinement ratios at each level (length at least mxnest-1)
    amrdata.refinement_ratios_x = [2,4,4,4,2]
    amrdata.refinement_ratios_y = [2,4,4,4,2]
    amrdata.refinement_ratios_t = [2,4,4,4,2]


    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length maux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).

    amrdata.aux_type = ['center','capacity','yleft','center','center']


    # Flag using refinement routine flag2refine rather than richardson error
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag_richardson_tol = 0.002  # Richardson tolerance
    amrdata.flag2refine = True

    # steps to take on each level L between regriddings of level L+1:
    amrdata.regrid_interval = 3

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    amrdata.regrid_buffer_width  = 3

    # clustering alg. cutoff for (# flagged pts) / (total # of cells refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged cells)
    amrdata.clustering_cutoff = 0.700000

    # print info about each regridding up to this level:
    amrdata.verbosity_regrid = 0  

    #  ----- For developers ----- 
    # Toggle debugging print statements:
    amrdata.dprint = False      # print domain flags
    amrdata.eprint = False      # print err est flags
    amrdata.edebug = False      # even more err est flags
    amrdata.gprint = False      # grid bisection/clustering
    amrdata.nprint = False      # proper nesting output
    amrdata.pprint = False      # proj. of tagged points
    amrdata.rprint = False      # print regridding summary
    amrdata.sprint = False      # space/memory output
    amrdata.tprint = True       # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting
    
    # More AMR parameters can be set -- see the defaults in pyclaw/data.py

    # ---------------
    # Regions:
    # ---------------
    rundata.regiondata.regions = []
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]

    if 0:
        # Allow only level 1 as default everywhere:
        rundata.regiondata.regions.append([1, 1, 0., 1e9, -180, 180, -90, 90])

        # Force refinement around earthquake source region for first hour:
        rundata.regiondata.regions.append([3, 3, 0., 3600., -85,-72,-38,-25])

        # Allow up to level 3 in northeastern part of domain:
        rundata.regiondata.regions.append([1, 3, 0., 1.e9, -90,-60,-30,0])

    # ---------------
    # Gauges:
    # ---------------
    rundata.gaugedata.gauges = []
    # for gauges append lines of the form  [gaugeno, x, y, t1, t2]
    #rundata.gaugedata.gauges.append([32412, -86.392, -17.975, 0., 1.e10])
    

    return rundata
    # end of function setrun
    # ----------------------


#-------------------
def setgeo(rundata):
#-------------------
    """
    Set GeoClaw specific runtime parameters.
    For documentation see ....
    """

    try:
        geo_data = rundata.geo_data
    except:
        print("*** Error, this rundata has no geo_data attribute")
        raise AttributeError("Missing geo_data attribute")
       
    # == Physics ==
    geo_data.gravity = 9.81
    geo_data.coordinate_system = 2
    geo_data.earth_radius = 6367.5e3

    # == Forcing Options
    geo_data.coriolis_forcing = False

    # == Algorithm and Initial Conditions ==
    geo_data.sea_level = 0.0
    geo_data.dry_tolerance = 1.e-3
    geo_data.friction_forcing = False
    geo_data.manning_coefficient =.025
    geo_data.friction_depth = 1e6

    #Boussinesq data
    from clawpack.geoclaw.data import BoussData, BoussFlagRegionData
    rundata.add_data(BoussData(),'bouss_data')
    
    rundata.bouss_data.bouss_min_level = 1    # coarsest level to apply bouss (not used for EDC)
    rundata.bouss_data.bouss_max_level = 10   # finest level to apply bouss   (not used for EDC)
    rundata.bouss_data.bouss_solver = 3       # 1=GMRES, 2=Pardiso, 3=PETSc   (not used for EDC)
    rundata.bouss_data.bouss_tstart = 0.0      # time to switch from SWE       (not used for EDC)
    rundata.bouss_data.bouss_tfinal = 1.e10 

    rundata.bouss_data.bouss_equations = Bouss_Hyp_rel    # -1=EDC_HypRel, -2=BBBD_HypRel, 0=SWE, 1=MS, 2=SGN
    rundata.bouss_data.bouss_trans_low = 10#1.e3  # depth to switch to SWE
    rundata.bouss_data.bouss_trans_up = 100#1.e3
    rundata.bouss_data.bouss_EDC_c_sq = 4.e3*geo_data.gravity*3.# EDC wave speed
    rundata.bouss_data.bouss_EDC_gamma = 3./2.  # EDC gamma (2 for Saint-Marie's sytem)
    rundata.bouss_data.bouss_csq_index = 3  # index of c^2 in aux array (Python indexing)
    rundata.bouss_data.bouss_decay_rate_index = 4  # index of decay rate in aux array (Python indexing)
    rundata.bouss_data.bouss_transition_type = 0  #0=depth-based, 1=distance_based
    rundata.bouss_data.bouss_transition_type_fun = 1  # depth to start transition

    

    from clawpack.amrclaw.data import FlagRegion
    # nearshore_continental = FlagRegion(num_dim=2)
    # nearshore_continental.name = 'RR_nearshore_continental'
    # nearshore_continental.minlevel = 1
    # nearshore_continental.maxlevel = 10
    # nearshore_continental.t1 = 0.0
    # nearshore_continental.t2 = 1e9
    # nearshore_continental.spatial_region_type = 2 
    # nearshore_continental.spatial_region_file = \
    #         os.path.abspath('rr_nearshore_continental_SouthAmerica.data')
    # bouss_flag_regions.append(nearshore_continental)

    # nearshore_ocean = FlagRegion(num_dim=2)
    # nearshore_ocean.name = 'RR_nearshore_ocean'
    # nearshore_ocean.minlevel = 1
    # nearshore_ocean.maxlevel = 10
    # nearshore_ocean.t1 = 0.0
    # nearshore_ocean.t2 = 1e9
    # nearshore_ocean.spatial_region_type = 2 
    # nearshore_ocean.spatial_region_file = \
    #         os.path.abspath('rr_nearshore_ocean_SouthAmerica.data')
    # bouss_flag_regions.append(nearshore_ocean)

    thisdir = os.path.dirname(__file__)
    #Region for forced refinement
    flagregions = rundata.flagregiondata.flagregions
    initial_conditions = FlagRegion(num_dim=2)
    initial_conditions.name = 'RR_initial_conditions'
    initial_conditions.minlevel = 3
    initial_conditions.maxlevel = 3
    initial_conditions.t1 = 0.0
    initial_conditions.t2 = 100.0
    initial_conditions.spatial_region_type = 2
    initial_conditions.spatial_region_file = os.path.join(thisdir,"flag2refine_rrs/rr1.data")
    flagregions.append(initial_conditions)

    #Set Bouss regions
    rundata.add_data(BoussFlagRegionData(),'bouss_regions_data')
    bouss_flag_regions = rundata.bouss_regions_data.flagregions

    
    culstered_rrs_path = os.path.join(thisdir,"clustered_rrs")
    for file in os.listdir(culstered_rrs_path):
        if file.endswith(".data"):# and file.startswith("rr"):
            rr =  FlagRegion(num_dim=2)
        #Strip the extension
        rr.name = os.path.splitext(file)[0]
        rr.minlevel = 1
        rr.maxlevel = 10
        rr.t1 = 0.0
        rr.t2 = 1.e9
        rr.spatial_region_type = 2
        rr.spatial_region_file = \
            os.path.abspath(os.path.join(culstered_rrs_path,file))
        bouss_flag_regions.append(rr)
 



    # Refinement settings
    refinement_data = rundata.refinement_data
    refinement_data.variable_dt_refinement_ratios = True
    refinement_data.wave_tolerance = 0.05

    # == settopo.data values ==
    topo_data = rundata.topo_data
    # for topography, append lines of the form
    #    [topotype, fname]
    # topo_path = 'etopo_no_dry.asc'
    topo_path = os.path.join(scratch_dir, 'etopo10min120W60W60S0S.asc')
    # topo_path = 'etopo10min120W60W60S0S_smoothed.tt2'
    topo_data.topofiles.append([2, topo_path])

    # == setdtopo.data values ==
    dtopo_data = rundata.dtopo_data
    # for moving topography, append lines of the form :   (<= 1 allowed for now!)
    #   [topotype, fname]
    # dtopo_path = os.path.join(scratch_dir, 'dtopo_usgs100227.tt3')
    # dtopo_data.dtopofiles.append([3,dtopo_path])
    # dtopo_data.dt_max_dtopo = 0.2


    # == setqinit.data values ==
    rundata.qinit_data.qinit_type = 0
    rundata.qinit_data.qinitfiles = []
    # for qinit perturbations, append lines of the form: (<= 1 allowed for now!)
    #   [minlev, maxlev, fname]

    return rundata
    # end of function setgeo
    # ----------------------



if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    from clawpack.geoclaw import kmltools

    rundata = setrun(*sys.argv[1:])
    rundata.write()

    kmltools.make_input_data_kmls(rundata)
