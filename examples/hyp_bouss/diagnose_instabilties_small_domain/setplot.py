
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt


#--------------------------
def setplot(plotdata=None):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from clawpack.visclaw import colormaps, geoplot
    from clawpack.visclaw.data import ClawPlotData

    from numpy import linspace

    if plotdata is None:
        plotdata = ClawPlotData()

    plotdata.clearfigures()  # clear any old figures,axes,items data


    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:

    def addgauges(current_data):
        from clawpack.visclaw import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos='all', format_string='ko', add_labels=True)
    
    def fixup(current_data):
        import pylab
        #addgauges(current_data)
        t = current_data.t
        t = t / 3600.  # hours
        pylab.title('Surface at %4.2f hours' % t, fontsize=15)
        #pylab.xticks(fontsize=15)
        #pylab.yticks(fontsize=15)

    def fixup2(current_data):
        import pylab
        #addgauges(current_data)
        t = current_data.t
        t = t / 3600.  # hours
        pylab.title('System used at %4.2f hours' % t, fontsize=15)
        #pylab.xticks(fontsize=15)
        #pylab.yticks(fontsize=15)

    def fixup2(current_data):
        import pylab
        #addgauges(current_data)
        t = current_data.t
        t = t / 3600.  # hours
        pylab.title(r'Velocity $|u| + |v|$'+'at %4.2f hours' % t, fontsize=15)
        #pylab.xticks(fontsize=15)
        #pylab.yticks(fontsize=15)


    def decayrate(current_data):
        return current_data.aux[4,:,:]
    
    def abs_momentum(current_data):
        h = current_data.q[0,:,:]
        u = np.where(h>1.e-3, current_data.q[1,:,:]/h, 0.0)
        v = np.where(h>1.e-3, current_data.q[2,:,:]/h, 0.0)
        return np.abs(u) + np.abs(v)

    #-----------------------------------------
    # Figure for surface
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.scaled = True
    plotaxes.afteraxes = fixup

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    #plotitem.plot_var = geoplot.surface
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -0.1#-1.e-4#-5.e-9
    plotitem.pcolor_cmax = 0.1#1.e-4#5.e-9
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0,0,0]
    plotitem.patchedges_show = [0,0,0,0,0]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = 5*[False]#[0,0,0,0,0]
    plotitem.patchedges_show = 5*[False]# [0,0,0,0,0]
    plotaxes.xlimits = [-92.0,-90.0]#[-120,-60]
    plotaxes.ylimits = [-2.0,0]#[-60,0]

    # #Auxiliary var 2
    plotfigure = plotdata.new_plotfigure(name='System', figno=1)
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'System'
    plotaxes.scaled = True
    plotaxes.afteraxes = fixup2
    
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = decayrate
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 1.0
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = 5*[False]#[0,0,0,0,0]
    plotitem.patchedges_show = 5*[False]#[0]*5#[1,1,1,1,1]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = 5*[False]#[0]*5#[1,0,0,0,0]
    plotitem.patchedges_show = 5*[False]#[0]*5

    plotaxes.xlimits = [-92,-90]#[-120,-60]
    plotaxes.ylimits = [-2,0]#[-60,0]
    import os
    from clawpack.amrclaw import region_tools
    thisdir = os.path.dirname(__file__)
    def boundary_ruled_rectangle(current_data):
        #Iterate over all files in this dir that end with .data
        #and start with Ruled_Rectangle
        culstered_rrs_path = os.path.join(thisdir,"clustered_rrs")
        for file in os.listdir(culstered_rrs_path):
            if file.endswith(".data"):# and file.startswith("rr"):
                #Read the file
                rr = region_tools.RuledRectangle()
                rr.read(os.path.join(culstered_rrs_path,file))
                #Get the vertices of the rectangle
                x,y = rr.vertices()
                #Plot with green color and dashed lines
                plt.plot(x,y,'purple',linewidth=0.5)
    plotaxes.afteraxes = boundary_ruled_rectangle


    # Velocities
    plotfigure = plotdata.new_plotfigure(name='momentum', figno=2)
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Velocity |u| + |v|'
    plotaxes.scaled = True
    plotaxes.afteraxes = fixup2
    
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = abs_momentum
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 0.1
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = 5*[False]#[0,0,0,0,0]
    plotitem.patchedges_show = [0]*5#[1,1,1,1,1]




    # #-----------------------------------------
    # # Figures for gauges
    # #-----------------------------------------
    # plotfigure = plotdata.new_plotfigure(name='Surface at gauges', figno=300, \
    #                 type='each_gauge')
    # plotfigure.clf_each_gauge = True

    # # Set up for axes in this figure:
    # plotaxes = plotfigure.new_plotaxes()
    # plotaxes.xlimits = 'auto'
    # plotaxes.ylimits = 'auto'
    # plotaxes.title = 'Surface'

    # # Plot surface as blue curve:
    # plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    # plotitem.plot_var = 5
    # plotitem.plotstyle = 'b-'


    # def add_zeroline(current_data):
    #     from pylab import plot, legend, xticks, floor, axis, xlabel
    #     t = current_data.t 
    #     gaugeno = current_data.gaugeno

    #     plot(t, 0*t, 'k')
    #     n = int(floor(t.max()/3600.) + 2)
    #     xticks([3600*i for i in range(n)], ['%i' % i for i in range(n)])
    #     xlabel('time (hours)')

    # plotaxes.afteraxes = add_zeroline



    #-----------------------------------------
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = [15*i for i in range(0,41)]
    #'all'          # list of frames to print
    plotdata.print_gaugenos = 'all'          # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

