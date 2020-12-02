#!/usr/bin/env python

# Copyright 2019-2020 Juliane Mai - juliane.mai(at)uwaterloo.ca
#
# License
# This file is part of Juliane Mai's personal code library.
#
# Juliane Mai's personal code library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Juliane Mai's personal code library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public Licensefstop
# along with Juliane Mai's personal code library.  If not, see <http://www.gnu.org/licenses/>.
#
# run with:
#     run figure_2.py -t pdf -p figure_2

from __future__ import print_function

"""

Plots Salmon River catchment and annual precipitation and temperature

History
-------
Written,  JM, Jul 2019
"""

# -------------------------------------------------------------------------
# Command line arguments - if script
#

# Comment|Uncomment - Begin
if __name__ == '__main__':

    import argparse
    import numpy as np

    plotname    = ''
    outtype     = ''
    usetex      = False
    serif       = False
    doabc       = True
    nsets       = 10            # number of Sobol sequences
    nboot       = 1             # Set to 1 for single run of SI and STI calculation
    
    parser   = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                      description='''Benchmark example to test Sensitivity Analysis for models with multiple process options.''')
    parser.add_argument('-p', '--plotname', action='store',
                        default=plotname, dest='plotname', metavar='plotname',
                        help='Name of plot output file for types pdf, html or d3, '
                        'and name basis for type png (default: '+__file__[0:__file__.rfind(".")]+').')
    parser.add_argument('-s', '--serif', action='store_true', default=serif, dest="serif",
                    help="Use serif font; default sans serif.")
    parser.add_argument('-t', '--type', action='store',
                        default=outtype, dest='outtype', metavar='outtype',
                        help='Output type is pdf, png, html, or d3 (default: open screen windows).')
    parser.add_argument('-u', '--usetex', action='store_true', default=usetex, dest="usetex",
                        help="Use LaTeX to render text in pdf, png and html.")
    parser.add_argument('-b', '--nboot', action='store',
                        default=nboot, dest='nboot', metavar='nboot',
                        help='Number of bootstrap samples (default: nboot=10).')
    parser.add_argument('-n', '--nsets', action='store',
                        default=nsets, dest='nsets', metavar='nsets',
                        help='Number of sensitivity samples (default: nsets=10).')

    args     = parser.parse_args()
    plotname = args.plotname
    outtype  = args.outtype
    serif    = args.serif
    usetex   = args.usetex
    nboot    = np.int(args.nboot)
    nsets    = np.int(args.nsets)

    del parser, args
    # Comment|Uncomment - End

    # -----------------------
    # add subolder scripts/lib to search path
    # -----------------------
    import sys
    import os 
    dir_path = os.path.dirname(os.path.realpath(__file__))
    sys.path.append(dir_path+'/../lib')


# -------------------------------------------------------------------------
# Function definition - if function
#

    # Check input
    outtype = outtype.lower()
    outtypes = ['', 'pdf', 'png', 'html', 'd3']
    if outtype not in outtypes:
        print('\nError: output type must be in ', outtypes)
        import sys
        sys.exit()

    import numpy as np
    import xarray as xr
    import time
    import datetime
    import color                           # in lib/
    from   position        import position # in lib/
    from   autostring      import astr     # in lib/
    from   abc2plot        import abc2plot # in lib/
    from   fread           import fread    # in lib/
    from   str2tex         import str2tex  # in lib/

    
    t1 = time.time()

    if (outtype == 'd3'):
        try:
            import mpld3
        except:
            print("No mpld3 found. Use output type html instead of d3.")
            outtype = 'html'


    # -------------------------------------------------------------------------
    # Setup
    #
    dowhite    = False  # True: black background, False: white background
    title      = False   # True: title on plots, False: no plot titles
    textbox    = False  # if true: additional information is set as text box within plot
    textbox_x  = 0.95
    textbox_y  = 0.85

    # -------------------------------------------------------------------------
    # Setup Calculations
    #
    if dowhite:
        fgcolor = 'white'
        bgcolor = 'black'
    else:
        fgcolor = 'black'
        bgcolor = 'white'

    # colors
    cols1 = color.get_brewer('YlOrRd9', rgb=True)
    cols1 = color.get_brewer( 'WhiteYellowOrangeRed',rgb=True)[30:]
    cols1 = color.get_brewer( 'dark_rainbow_256',rgb=True)   # blue to red

    cols2 = color.get_brewer('YlOrRd9', rgb=True)[::-1]
    cols2 = color.get_brewer( 'WhiteYellowOrangeRed',rgb=True)[30:][::-1]
    cols2 = color.get_brewer( 'dark_rainbow_256',rgb=True)[::-1]  # red to blue

    cols3 = [cols2[0],cols2[95],cols2[-1]]  # red, yellow, blue
    cols3 = [color.colours('gray'),cols2[0],color.colours('white')]  # gray red white

    # -------------------------------------------------------------------------
    # Read latlon info from file
    # -------------------------------------------------------------------------
    latlon = fread("../examples/data_in/data_obs/08KC001/salmon_latlon.dat",skip=1)

    # -------------------------------------------------------------------------
    # Read forcings from file
    # -------------------------------------------------------------------------
    forcing = fread("../examples/data_in/data_obs/08KC001/Salmon-River-Near-Prince-George_meteo_daily.rvt", skip=4, comment=":")
    ntime     = np.shape(forcing)[0]
    ntime_doy = 365
    forc_time = np.array([ datetime.datetime(1954,1,1,0,0) + datetime.timedelta(itime) for itime in range(ntime) ])
    leap      = np.array([ True if (forc_time[ii].month == 2 and forc_time[ii].day == 29) else False for ii in range(ntime) ])

    # RAINFALL   SNOWFALL   TEMP_DAILY_MIN   TEMP_DAILY_MAX   PET
    rain      = forcing[~leap,0]
    snow      = forcing[~leap,1]
    tmin      = forcing[~leap,2]
    tmax      = forcing[~leap,3]
    pet       = forcing[~leap,4]
    pre       = forcing[~leap,0]+forcing[~leap,1]
    tavg      = (tmin+tmax)/2.0

    # average over day of year
    forc_time = forc_time[~leap]
    pre_doy  = np.median(np.reshape(pre[:],  [int(ntime/ntime_doy),ntime_doy]),axis=0)
    snow_doy = np.median(np.reshape(snow[:], [int(ntime/ntime_doy),ntime_doy]),axis=0)
    rain_doy = np.median(np.reshape(rain[:], [int(ntime/ntime_doy),ntime_doy]),axis=0)
    tavg_doy = np.mean(np.reshape(tavg[:],   [int(ntime/ntime_doy),ntime_doy]),axis=0)

    # mean annual precip
    years = [ tt.year for tt in forc_time ]
    uyears = np.unique(years)
    annual_precip = []
    for uu in uyears:
        idx = np.where(years == uu)
        annual_precip.append(np.sum(pre[idx]))
    annual_temp = np.mean(tavg_doy)
    # print("annual_precip      = ", annual_precip)
    print()
    print('---------------------------------------------')
    print("mean annual temp   = ", np.mean(annual_temp))
    print("mean annual precip = ", np.mean(annual_precip))
    print('---------------------------------------------')
    print()

    # mean monthly temperature and average total precip per month
    years  = [ tt.year for tt in forc_time ]
    months = [ tt.month for tt in forc_time ]
    uyears = np.unique(years)
    month_prec = []
    month_snow = []
    month_rain = []
    month_temp = []
    for uu in uyears:
        tmp_prec = []
        tmp_snow = []
        tmp_rain = []
        tmp_temp = []
        for mm in np.arange(1,13):
            idx = np.where((years == uu) & (months == mm))
            tmp_prec.append(np.sum(pre[idx]))
            tmp_snow.append(np.sum(snow[idx]))
            tmp_rain.append(np.sum(rain[idx]))
            tmp_temp.append(np.mean(tavg[idx]))
        month_prec.append(tmp_prec)
        month_snow.append(tmp_snow)
        month_rain.append(tmp_rain)
        month_temp.append(tmp_temp)
    month_prec = np.mean(np.array(month_prec),axis=0)
    month_snow = np.mean(np.array(month_snow),axis=0)
    month_rain = np.mean(np.array(month_rain),axis=0)
    month_temp = np.mean(np.array(month_temp),axis=0)

    # -------------------------------------------------------------------------
    # Colors
    # -------------------------------------------------------------------------
    ocean_color = (151/256., 183/256., 224/256.)
    infil_color = color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[20]
    quick_color = color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[55]
    evapo_color = color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[80]
    basef_color = color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[105]
    snowb_color = color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[130]
    convs_color = color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[155]
    convd_color = color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[180]
    potme_color = color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[205]
    perco_color = color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[230]
    soilm_color = color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[255]
    
    # -------------------------------------------------------------------------
    # Plotting of results
    # -------------------------------------------------------------------------
    # Main plot
    ncol        = 2           # number columns
    nrow        = 4           # number of rows
    textsize    = 10          # standard text size
    dxabc       = 0.03          # % of (max-min) shift to the right from left y-axis for a,b,c,... labels
    dyabc       = 0.92          # % of (max-min) shift up from lower x-axis for a,b,c,... labels
    dxsig       = 1.23        # % of (max-min) shift to the right from left y-axis for signature
    dysig       = -0.075      # % of (max-min) shift up from lower x-axis for signature
    dxtit       = 0           # % of (max-min) shift to the right from left y-axis for title
    dytit       = 1.2         # % of (max-min) shift up from lower x-axis for title
    hspace      = 0.16        # x-space between subplots
    vspace      = 0.06        # y-space between subplots

    lwidth      = 0.5         # linewidth
    elwidth     = 0.5         # errorbar line width
    alwidth     = 1.0         # axis line width
    glwidth     = 0.5         # grid line width
    msize       = 8.0         # marker size
    mwidth      = 0.0         # marker edge width
    mcol1       = '0.7'       # primary marker colour
    mcol2       = '0.0'       # secondary
    mcol3       = '0.0'       # third
    mcols       = color.colours(['blue','green','yellow','orange','red','darkgray','darkblue','black','darkgreen','gray'])
    lcol0       = color.colours('black')    # line colour
    lcol1       = color.colours('blue')     # line colour
    lcol2       = color.colours('green')    # line colour
    lcol3       = color.colours('yellow')   # line colour
    lcols       = color.colours(['black','blue','green','yellow'])
    markers     = ['o','v','s','^']

    # Legend
    llxbbox     = 0.98        # x-anchor legend bounding box
    llybbox     = 0.98        # y-anchor legend bounding box
    llrspace    = 0.          # spacing between rows in legend
    llcspace    = 1.0         # spacing between columns in legend
    llhtextpad  = 0.4         # the pad between the legend handle and text
    llhlength   = 1.5         # the length of the legend handles
    frameon     = False       # if True, draw a frame around the legend. If None, use rc
      
    import matplotlib as mpl
    import matplotlib.patches as patches
    from matplotlib.patches import Rectangle, Circle, Polygon
    from mpl_toolkits.basemap import Basemap
    mpl.use('TkAgg')
    
    if (outtype == 'pdf'):
        mpl.use('PDF') # set directly after import matplotlib
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        # Customize: http://matplotlib.sourceforge.net/users/customizing.html
        mpl.rc('ps', papersize='a4', usedistiller='xpdf') # ps2pdf
        # mpl.rc('figure', figsize=(8.27,11.69)) # a4 portrait
        mpl.rc('figure', figsize=(7.48,9.06)) # WRR maximal figure size
        if usetex:
            mpl.rc('text', usetex=True)
            if not serif:
                #   r'\usepackage{helvet}',                             # use Helvetica
                mpl.rcParams['text.latex.preamble'] = [
                    r'\usepackage[math,lf,mathtabular,footnotefigures]{MyriadPro}', # use MyriadPro font
                    r'\renewcommand{\familydefault}{\sfdefault}',       # normal text font is sans serif
                    r'\figureversion{lining,tabular}',
                    r'\usepackage{wasysym}',                            # for permil symbol (load after MyriadPro)
                    ]
            else:
                mpl.rcParams['text.latex.preamble'] = [
                    r'\usepackage{wasysym}'                     # for permil symbol
                    ]
        else:
            if serif:
                mpl.rcParams['font.family']     = 'serif'
                mpl.rcParams['font.sans-serif'] = 'Times'
            else:
                mpl.rcParams['font.family']     = 'sans-serif'
                mpl.rcParams['font.sans-serif'] = 'Arial'       # Arial, Verdana
    elif (outtype == 'png') or (outtype == 'html') or (outtype == 'd3'):
        mpl.use('Agg') # set directly after import matplotlib
        import matplotlib.pyplot as plt
        # mpl.rc('figure', figsize=(8.27,11.69)) # a4 portrait
        mpl.rc('figure', figsize=(7.48,9.06)) # WRR maximal figure size
        if usetex:
            mpl.rc('text', usetex=True)
            if not serif:
                #   r'\usepackage{helvet}',                             # use Helvetica
                mpl.rcParams['text.latex.preamble'] = [
                    r'\usepackage[math,lf,mathtabular,footnotefigures]{MyriadPro}', # use MyriadPro font
                    r'\renewcommand{\familydefault}{\sfdefault}',       # normal text font is sans serif
                    r'\figureversion{lining,tabular}',
                    r'\usepackage{wasysym}',                            # for permil symbol (load after MyriadPro)
                    ]
            else:
                mpl.rcParams['text.latex.preamble'] = [
                    r'\usepackage{wasysym}'                     # for permil symbol
                    ]
        else:
            if serif:
                mpl.rcParams['font.family']     = 'serif'
                mpl.rcParams['font.sans-serif'] = 'Times'
            else:
                mpl.rcParams['font.family']     = 'sans-serif'
                mpl.rcParams['font.sans-serif'] = 'Arial'       # Arial, Verdana
        mpl.rc('savefig', dpi=dpi, format='png')
    else:
        import matplotlib.pyplot as plt
        # mpl.rc('figure', figsize=(4./5.*8.27,4./5.*11.69)) # a4 portrait
        mpl.rc('figure', figsize=(7.48,9.06)) # WRR maximal figure size
    mpl.rc('text.latex') #, unicode=True)
    mpl.rc('font', size=textsize)
    mpl.rc('path', simplify=False) # do not remove
    # print(mpl.rcParams)
    mpl.rc('axes', linewidth=alwidth, edgecolor=fgcolor, facecolor=bgcolor, labelcolor=fgcolor)
    mpl.rc('figure', edgecolor=bgcolor, facecolor='grey')
    mpl.rc('grid', color=fgcolor)
    mpl.rc('lines', linewidth=lwidth, color=fgcolor)
    mpl.rc('patch', edgecolor=fgcolor)
    mpl.rc('savefig', edgecolor=bgcolor, facecolor=bgcolor)
    mpl.rc('patch', edgecolor=fgcolor)
    mpl.rc('text', color=fgcolor)
    mpl.rc('xtick', color=fgcolor)
    mpl.rc('ytick', color=fgcolor)

    if (outtype == 'pdf'):
        pdffile = plotname+'.pdf'
        print('Plot PDF ', pdffile)
        pdf_pages = PdfPages(pdffile)
    elif (outtype == 'png'):
        print('Plot PNG ', plotname)
    else:
        print('Plot X')

    t1  = time.time()
    ifig = 0

    figsize = mpl.rcParams['figure.figsize']
    mpl.rcParams['axes.linewidth'] = lwidth

    
    ifig = 0
    
    ifig += 1
    iplot = 0
    print('Plot - Fig ', ifig)
    fig = plt.figure(ifig)

    # -----------------------
    # plot
    # -----------------------
    ylim = [-0.1, 0.6]

    # -------------
    # Salmon River watershed
    # -------------
    iplot += 1
    #                           [left, bottom, width, height]
    sub = fig.add_axes(np.array([ 0.025 ,  0.745 ,  0.3475,  0.155 ]))

    # Map: Canada - Lake Erie
    llcrnrlon =  -130.0 #-81.25 #-85.5
    urcrnrlon =  -113.0 #-77.0
    llcrnrlat =   48.0 #39.5
    urcrnrlat =   56.0
    lat_1     =   52.0 # 42.0  # first  "equator"
    lat_2     =   52.0 # 42.0  # second "equator"
    lat_0     =   52.0 # 42.5  # center of the map
    lon_0     =  -121.5 #-79.00 #-82.0  # center of the map
    # m = Basemap(projection='lcc',
    #             llcrnrlon=-80, llcrnrlat=43, urcrnrlon=-75, urcrnrlat=47,
    #             lon_0=-77.5, lat_0=43, 
    #             lat_1=44, lat_2=44, 
    #             resolution='i') # Lambert conformal
    m = Basemap(projection='lcc', area_thresh=10000.,
                llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon, llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat,
                lat_1=lat_1, lat_2=lat_2, lat_0=lat_0, lon_0=lon_0,
                resolution=None) # Lambert conformal

    # draw parallels and meridians.
    # labels: [left, right, top, bottom]
    if iplot%ncol == 1:  # plot in 1st column
        m.drawparallels(np.arange(-80.,81.,5.),  labels=[1,0,0,0], dashes=[1,1], linewidth=0.25, color='0.5')
    elif iplot%ncol == 0:  # plot in last column
        m.drawparallels(np.arange(-80.,81.,5.),  labels=[0,1,0,0], dashes=[1,1], linewidth=0.25, color='0.5')
    else:
        m.drawparallels(np.arange(-80.,81.,5.),  labels=[0,0,0,0], dashes=[1,1], linewidth=0.25, color='0.5')
    m.drawmeridians(np.arange(-180.,181.,6.),labels=[0,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5')
    # draw cooastlines and countries
    #m.drawcoastlines(linewidth=0.3)
    #m.drawmapboundary(fill_color=ocean_color, linewidth=0.3)
    #m.drawcountries(color='black', linewidth=0.3)
    # m.fillcontinents(lake_color=ocean_color)  # color='white', 
    m.shadedrelief()
    # m.drawlsmask()

    # convert coordinates to map coordinates
    xpt,ypt = m(latlon[:,1],latlon[:,0])
    mlatlon = np.transpose(np.array([xpt,ypt])) #zip(xpt,ypt)

    # Catchments
    facecolor = perco_color
    sub.add_patch(Polygon(mlatlon, facecolor=facecolor, edgecolor='none', alpha=1.0))

    # some reference points
    nref = 4
    ref_latlon = np.ones([nref,2]) * -9999.0
    ref_names  = []
    ref_latlon[0,0] = 49.2827 ; ref_latlon[0,1] = -123.1207  ; ref_names.append( 'Vancouver'     )   # Vancouver
    ref_latlon[1,0] = 53.9171 ; ref_latlon[1,1] = -122.7497  ; ref_names.append( 'Prince George' )   # Prince George
    ref_latlon[2,0] = 52.8737 ; ref_latlon[2,1] = -118.0814  ; ref_names.append( 'Jasper'        )   # Jasper
    ref_latlon[3,0] = 51.1784 ; ref_latlon[3,1] = -115.5708  ; ref_names.append( 'Banff '        )   # Banff

    for iref in range(nref):
            xpt1, ypt1 = m(ref_latlon[iref,1],ref_latlon[iref,0])
            sub.plot(xpt1, ypt1,
                     #linewidth=lwidth,
                     #color=color,
                     marker='o',
                     markeredgecolor='black',
                     markerfacecolor='black',
                     markersize=3.0, #msize,
                     markeredgewidth=mwidth)
            ha = ['right', 'left']
            va = ['bottom','top']

            if iref in [0,1,2]:
                ha = 'left'
                dx = 0.3
            else:
                ha = 'right'
                dx = -0.3
            x2, y2 = m(ref_latlon[iref,1]+dx,ref_latlon[iref,0])
            sub.annotate(ref_names[iref],
                xy=(xpt1, ypt1),  xycoords='data',
                xytext=(x2, y2), textcoords='data',
                fontsize='small',
                verticalalignment='center',horizontalalignment=ha#,
                # arrowprops=dict(arrowstyle="->",relpos=(1.0,0.5),linewidth=0.4)
                )

    # Title
    #sub.set_title(str2tex('Salmon River watershed',usetex=usetex),fontsize=textsize)

    # Fake subplot for numbering
    if doabc:
        pos = position(nrow,ncol,iplot,hspace=hspace,vspace=vspace)
        lsub = fig.add_axes([0.29, 0.88, 0.15, 0.0155])

        lsub.set_xlim([0,1])
        lsub.set_ylim([0,1])

        # subplot numbering
        abc2plot(lsub, dxabc, dyabc, iplot, lower=False, parenthesis='none',
                     bold=True, large=True,
                     mathrm=True, usetex=usetex,
                     horizontalalignment='left', verticalalignment='top')

        lsub.set_title('')
        lsub.set_xlabel('')
        lsub.set_ylabel('')
        lsub.set_xticks([])
        lsub.set_yticks([])
        lsub.set_axis_off()

    # -------------
    # Annual precipitation and temperature
    # -------------
    iplot += 1
    #                           [left, bottom, width, height]
    sub = fig.add_axes(np.array([ 0.4125,  0.745 ,  0.4475,  0.155 ]))

    # make color a bit darker for reviewer (use color of "water" in panels C and D)
    quick_color_dark = (5/256.,113/256.,176/256.)

    #ppre   = sub.plot( np.arange(ntime_doy), pre_doy,  color = quick_color)
    #psnow  = sub.plot( np.arange(ntime_doy), snow_doy,  color = infil_color)
    #prain  = sub.plot( np.arange(ntime_doy), rain_doy,  color = quick_color)
    width = 0.9
    psnow  = sub.bar( np.arange(12), month_snow, width, color = quick_color, bottom = month_rain)
    prain  = sub.bar( np.arange(12), month_rain, width, color = quick_color_dark)

    sub2 = sub.twinx()
    # ptavg = sub2.plot( np.arange(ntime_doy), tavg_doy, color = perco_color)
    ptavg = sub2.plot( np.arange(12), month_temp, color = perco_color, linewidth = 2*lwidth)

    print("month_snow = ",month_snow)
    print("month_rain = ",month_rain)
    print("month_prec = ",month_prec)
    print("month_temp = ",month_temp)

    ylim = [0,82]
    sub.set_ylim(ylim)
    ylim = [-15,18]
    sub2.set_ylim(ylim)

    # nopt = len(procopt)
    #sub.set_xticks(np.arange(12), ['J','F','M','A','M','J','J','A','S','O','N','D'])
    plt.xticks(np.arange(12), ['J','F','M','A','M','J','J','A','S','O','N','D'])
    
    # plt.title(str2tex('Sensitivities of Process Options',usetex=usetex))
    sub.set_xlabel(str2tex("Month",usetex=usetex), color='black')
    sub.set_ylabel(str2tex("Avg. Monthly Precip. [mm]",usetex=usetex), color=quick_color_dark)
    sub2.set_ylabel(str2tex('Avg. Monthly Temp. [$^\circ C$]',usetex=usetex), color=perco_color)

    # Create custom artists
    #      (left, bottom), width, height
    boxSTi_1  = patches.Rectangle( (0.045, 0.92), 0.02, 0.05, color = quick_color,      alpha=0.8, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_2  = patches.Rectangle( (0.045, 0.84), 0.02, 0.05, color = quick_color_dark, alpha=0.8, fill  = True, transform=sub.transAxes, clip_on=False )
    line      = patches.Rectangle( (0.045, 0.78), 0.02, 0.00, color = perco_color,      alpha=1.0, fill  = True, transform=sub.transAxes, clip_on=False )
    sub.add_patch(boxSTi_1)
    sub.add_patch(boxSTi_2)
    sub.add_patch(line)
    sub.text(0.1, 0.94, str2tex("Snow",usetex=usetex),         fontsize='small', horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)
    sub.text(0.1, 0.86, str2tex("Rain",usetex=usetex),         fontsize='small', horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)
    sub.text(0.1, 0.78, str2tex("Temperature",usetex=usetex),  fontsize='small', horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)
   
    abc2plot(sub,0.93,0.96,iplot,bold=True,usetex=usetex,mathrm=True, large=True, parenthesis='none',verticalalignment='top')

    # -------------
    # Salmon River watershed (soil)
    # -------------
    iplot += 1
    #                           [left, bottom, width, height]
    lleft   = -0.0615
    bbottom = 0.435
    wwidth  = 0.5 #0.3475
    hheight = 0.25 #0.155
    sub = fig.add_axes(np.array([ lleft, bbottom, wwidth, hheight ]))

    # Map: Canada - Lake Erie
    llcrnrlon =  -124.1 #-81.25 #-85.5
    urcrnrlon =  -122.5 #-77.0
    llcrnrlat =   53.85 #39.5
    urcrnrlat =   55.1
    lat_1     =   (llcrnrlat+urcrnrlat)/2.0 # 42.0  # first  "equator"
    lat_2     =   (llcrnrlat+urcrnrlat)/2.0 # 42.0  # second "equator"
    lat_0     =   (llcrnrlat+urcrnrlat)/2.0 # 42.5  # center of the map
    lon_0     =   (llcrnrlon+urcrnrlon)/2.0 #-79.00 #-82.0  # center of the map
    # m = Basemap(projection='lcc',
    #             llcrnrlon=-80, llcrnrlat=43, urcrnrlon=-75, urcrnrlat=47,
    #             lon_0=-77.5, lat_0=43, 
    #             lat_1=44, lat_2=44, 
    #             resolution='i') # Lambert conformal
    m = Basemap(projection='lcc', #area_thresh=10000.,
                llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon, llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat,
                lat_1=lat_1, lat_2=lat_2, lat_0=lat_0, lon_0=lon_0,
                resolution='i') # Lambert conformal

    # draw parallels and meridians.
    # labels: [left, right, top, bottom]
    if iplot%ncol == 1:  # plot in 1st column
        m.drawparallels(np.arange(-80.,81.,1.),  labels=[1,0,0,0], dashes=[1,1], linewidth=0.25, color='0.5')
    elif iplot%ncol == 0:  # plot in last column
        m.drawparallels(np.arange(-80.,81.,1.),  labels=[0,1,0,0], dashes=[1,1], linewidth=0.25, color='0.5')
    else:
        m.drawparallels(np.arange(-80.,81.,1.),  labels=[0,0,0,0], dashes=[1,1], linewidth=0.25, color='0.5')
    m.drawmeridians(np.arange(-180.,181.,1.),labels=[0,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5')

    # draw cooastlines and countries
    # m.drawcoastlines(linewidth=0.3)
    m.drawmapboundary(fill_color=ocean_color, linewidth=0.3)
    m.drawcountries(color='black', linewidth=0.3)
    m.fillcontinents(color='white', lake_color=ocean_color)

    # scalebar
    length = 40.0 # km
    m.drawmapscale(llcrnrlon*0.45+urcrnrlon*0.55, llcrnrlat*0.08+urcrnrlat*0.92, lon_0, lat_0, length, barstyle='fancy', fontsize='small')

    # ------------------------------
    # read data from NetCDF
    # ------------------------------
    ds       = xr.open_dataset('../examples/data_in/data_obs/08KC001/salmon_soil.nc')
    lon      = ds['lon'].data        # 1D or 2D field
    lat      = ds['lat'].data        # 1D or 2D field
    var      = ds['soil'][0].data

    # 1D lats and lons
    if (len(np.shape(lon)) == 1):
        nlat = lat.shape[0]
        nlon = lon.shape[0]
        lon       = np.array([ lon for ilat in range(nlat) ])
        lon[:,-1] = 359.999999   # this is a hack
        lat = np.transpose(np.array([ lat for ilon in range(nlon) ]))
    elif (len(np.shape(lon)) == 2):
        nlat = lat.shape[0]
        nlon = lon.shape[1]
    else:
        raise ValueError('plot_CASPAR_data: lon and lat has to be either 1D or 2D')

    # boundaries between lats and lons
    lonh = np.empty((nlat+1,nlon+1), dtype=np.float)
    lath = np.empty((nlat+1,nlon+1), dtype=np.float)

    tmp = [ [ (lat[ii+1,jj+1]-lat[ii,jj])/2 for jj in range(nlon-1) ] + [ (lat[ii+1,nlon-1]-lat[ii,nlon-2])/2 ] for ii in range(nlat-1) ]
    dlat = np.array(tmp + [ tmp[-1] ])

    tmp = [ [ (lon[ii+1,jj+1]-lon[ii,jj])/2 for jj in range(nlon-1) ] + [ (lon[ii+1,nlon-1]-lon[ii,nlon-2])/2 ] for ii in range(nlat-1) ]
    dlon = np.array(tmp + [ tmp[-1] ])

    lonh[0:nlat,0:nlon] = lon - dlon
    lath[0:nlat,0:nlon] = lat - dlat

    # make lat and lon one column and row wider such that all
    lonh[nlat,0:nlon] = lonh[nlat-1,0:nlon] + (lonh[nlat-1,0:nlon] - lonh[nlat-2,0:nlon])  
    lath[nlat,0:nlon] = lath[nlat-1,0:nlon] + (lath[nlat-1,0:nlon] - lath[nlat-2,0:nlon])  
    lonh[0:nlat,nlon] = lonh[0:nlat,nlon-1] + (lonh[0:nlat,nlon-1] - lonh[0:nlat,nlon-2])  
    lath[0:nlat,nlon] = lath[0:nlat,nlon-1] + (lath[0:nlat,nlon-1] - lath[0:nlat,nlon-2])  
    lonh[nlat,nlon]   = lonh[nlat-1,nlon-1] + (lonh[nlat-1,nlon-1] - lonh[nlat-2,nlon-2])  
    lath[nlat,nlon]   = lath[nlat-1,nlon-1] + (lath[nlat-1,nlon-1] - lath[nlat-2,nlon-2])

    # geo-referenced
    xx,  yy  = m(lon,lat)
    xxh, yyh = m(lonh,lath)
    zz = var[:,:]  
    var_minval = np.nanmin(var)
    var_maxval = np.nanmax(var)

    # colors
    # RGB (166,97,26 )		1 - Loam
	# RGB (223,194,125 ) 	2 - Sandy clay loam
	# RGB (128,205,193 )	3 - Loamy sand
	# RGB (5,113,176)		5 - Water
    cc = [(166/256.,97/256.,26/256.),(223/256.,194/256.,125/256.),(128/256.,205/256.,193/256.),(0.2,0.2,0.2),(5/256.,113/256.,176/256.)]
    cmap = mpl.colors.ListedColormap(cc)

    # plot
    variable_plot = m.pcolor(xxh, yyh, zz, cmap=cmap, zorder = 100, alpha=0.9)    

    # convert coordinates to map coordinates
    xpt,ypt = m(latlon[:,1],latlon[:,0])
    mlatlon = np.transpose(np.array([xpt,ypt])) # zip(xpt,ypt)

    # Catchments
    facecolor = perco_color
    sub.add_patch(Polygon(mlatlon, facecolor='None', linewidth=lwidth, edgecolor='black', alpha=1.0, zorder=200))

    # some reference points
    nref = 4
    ref_latlon = np.ones([nref,2]) * -9999.0
    ref_names  = []
    ref_latlon[0,0] = 49.2827 ; ref_latlon[0,1] = -123.1207  ; ref_names.append( 'Vancouver'     )   # Vancouver
    ref_latlon[1,0] = 53.9171 ; ref_latlon[1,1] = -122.7497  ; ref_names.append( 'Prince George' )   # Prince George
    ref_latlon[2,0] = 52.8737 ; ref_latlon[2,1] = -118.0814  ; ref_names.append( 'Jasper'        )   # Jasper
    ref_latlon[3,0] = 51.1784 ; ref_latlon[3,1] = -115.5708  ; ref_names.append( 'Banff '        )   # Banff

    for iref in range(nref):
            xpt1, ypt1 = m(ref_latlon[iref,1],ref_latlon[iref,0])
            sub.plot(xpt1, ypt1,
                     #linewidth=lwidth,
                     #color=color,
                     marker='o',
                     markeredgecolor='black',
                     markerfacecolor='black',
                     markersize=3.0, #msize,
                     markeredgewidth=mwidth)
            ha = ['right', 'left']
            va = ['bottom','top']

            if iref in [0,2]:
                ha = 'left'
                dx = 0.1
            else:
                ha = 'right'
                dx = -0.1
            x2, y2 = m(ref_latlon[iref,1]+dx,ref_latlon[iref,0])
            sub.annotate(ref_names[iref],
                xy=(xpt1, ypt1),  xycoords='data',
                xytext=(x2, y2), textcoords='data',
                fontsize='small',
                verticalalignment='center',horizontalalignment=ha#,
                # arrowprops=dict(arrowstyle="->",relpos=(1.0,0.5),linewidth=0.4)
                )

    # Title
    #sub.set_title(str2tex('Salmon River watershed',usetex=usetex),fontsize=textsize)

    # percentages
    category_vals = np.unique(var)[~np.isnan(np.unique(var))]
    count_all = np.count_nonzero(~np.isnan(var))
    perc = []
    for icat in category_vals:
        count_cat = np.shape(np.where(var==icat)[0])[0]
        perc.append(count_cat*1.0/count_all*100.)

    print("percentages categories: ",perc)

    # Create custom artists
    #               RGB (166,97,26 )	1 - Loam
	#               RGB (223,194,125 ) 	2 - Sandy clay loam
	#               RGB (128,205,193 )	3 - Loamy sand
	#               RGB (5,113,176)		5 - Water
    #      (left, bottom), width, height
    boxSTi_1  = patches.Rectangle( (1.045, 0.92), 0.05, 0.05, color = cc[0], alpha=0.6, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_2  = patches.Rectangle( (1.045, 0.81), 0.05, 0.05, color = cc[1], alpha=0.6, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_3  = patches.Rectangle( (1.045, 0.70), 0.05, 0.05, color = cc[2], alpha=0.6, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_4  = patches.Rectangle( (1.045, 0.62), 0.05, 0.05, color = cc[4], alpha=0.6, fill  = True, transform=sub.transAxes, clip_on=False )
    sub.add_patch(boxSTi_1)
    sub.add_patch(boxSTi_2)
    sub.add_patch(boxSTi_3)
    sub.add_patch(boxSTi_4)
    sub.text(1.13, 0.94, str2tex("Loam ("+astr(perc[0],prec=0)+"%)",usetex=usetex),            fontsize='small', horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)
    sub.text(1.13, 0.83, str2tex("Sandy Clay\nLoam ("+astr(perc[1],prec=0)+"%)",usetex=usetex), fontsize='small', horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)
    sub.text(1.13, 0.72, str2tex("Loamy Sand ("+astr(perc[2],prec=0)+"%)",usetex=usetex),      fontsize='small', horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)
    sub.text(1.13, 0.64, str2tex("Water ("+astr(perc[3],prec=0)+"%)",usetex=usetex),           fontsize='small', horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)

    # Fake subplot for numbering
    if doabc:
        pos = position(nrow,ncol,iplot,hspace=hspace,vspace=vspace)
        lsub = fig.add_axes([lleft+wwidth-0.165, bbottom+hheight-0.098, 0.1, 0.1])

        lsub.set_xlim([0,1])
        lsub.set_ylim([0,1])

        # subplot numbering
        abc2plot(lsub, dxabc, dyabc, iplot, lower=False, parenthesis='none',
                     bold=True, large=True,
                     mathrm=True, usetex=usetex,
                     horizontalalignment='left', verticalalignment='top')

        lsub.set_title('')
        lsub.set_xlabel('')
        lsub.set_ylabel('')
        lsub.set_xticks([])
        lsub.set_yticks([])
        lsub.set_axis_off()

    sol = var

    # -------------
    # Salmon River watershed (landuse)
    # -------------
    iplot += 1
    #                           [left, bottom, width, height]
    lleft   = 0.38
    bbottom = 0.435
    wwidth  = 0.5 #0.3475
    hheight = 0.25 #0.155
    sub = fig.add_axes(np.array([ lleft, bbottom, wwidth, hheight ]))

    # Map: Canada - Lake Erie
    llcrnrlon =  -124.1 #-81.25 #-85.5
    urcrnrlon =  -122.5 #-77.0
    llcrnrlat =   53.85 #39.5
    urcrnrlat =   55.1
    lat_1     =   (llcrnrlat+urcrnrlat)/2.0 # 42.0  # first  "equator"
    lat_2     =   (llcrnrlat+urcrnrlat)/2.0 # 42.0  # second "equator"
    lat_0     =   (llcrnrlat+urcrnrlat)/2.0 # 42.5  # center of the map
    lon_0     =   (llcrnrlon+urcrnrlon)/2.0 #-79.00 #-82.0  # center of the map
    # m = Basemap(projection='lcc',
    #             llcrnrlon=-80, llcrnrlat=43, urcrnrlon=-75, urcrnrlat=47,
    #             lon_0=-77.5, lat_0=43, 
    #             lat_1=44, lat_2=44, 
    #             resolution='i') # Lambert conformal
    m = Basemap(projection='lcc', #area_thresh=10000.,
                llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon, llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat,
                lat_1=lat_1, lat_2=lat_2, lat_0=lat_0, lon_0=lon_0,
                resolution='i') # Lambert conformal

    # draw parallels and meridians.
    # labels: [left, right, top, bottom]
    m.drawparallels(np.arange(-80.,81.,1.),  labels=[1,0,0,0], dashes=[1,1], linewidth=0.25, color='0.5')
    m.drawmeridians(np.arange(-180.,181.,1.),labels=[0,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5')

    # draw cooastlines and countries
    # m.drawcoastlines(linewidth=0.3)
    m.drawmapboundary(fill_color=ocean_color, linewidth=0.3)
    m.drawcountries(color='black', linewidth=0.3)
    m.fillcontinents(color='white', lake_color=ocean_color)

    # scalebar
    length = 40.0 # km
    m.drawmapscale(llcrnrlon*0.45+urcrnrlon*0.55, llcrnrlat*0.08+urcrnrlat*0.92, lon_0, lat_0, length, barstyle='fancy', fontsize='small')

    # ------------------------------
    # read data from NetCDF
    # ------------------------------
    ds       = xr.open_dataset('../examples/data_in/data_obs/08KC001/salmon_landuse.nc')
    lon      = ds['lon'].data        # 1D or 2D field
    lat      = ds['lat'].data        # 1D or 2D field
    var      = ds['landuse'][0].data

    # 1D lats and lons
    if (len(np.shape(lon)) == 1):
        nlat = lat.shape[0]
        nlon = lon.shape[0]
        lon       = np.array([ lon for ilat in range(nlat) ])
        lon[:,-1] = 359.999999   # this is a hack
        lat = np.transpose(np.array([ lat for ilon in range(nlon) ]))
    elif (len(np.shape(lon)) == 2):
        nlat = lat.shape[0]
        nlon = lon.shape[1]
    else:
        raise ValueError('plot_CASPAR_data: lon and lat has to be either 1D or 2D')

    # boundaries between lats and lons
    lonh = np.empty((nlat+1,nlon+1), dtype=np.float)
    lath = np.empty((nlat+1,nlon+1), dtype=np.float)

    tmp = [ [ (lat[ii+1,jj+1]-lat[ii,jj])/2 for jj in range(nlon-1) ] + [ (lat[ii+1,nlon-1]-lat[ii,nlon-2])/2 ] for ii in range(nlat-1) ]
    dlat = np.array(tmp + [ tmp[-1] ])

    tmp = [ [ (lon[ii+1,jj+1]-lon[ii,jj])/2 for jj in range(nlon-1) ] + [ (lon[ii+1,nlon-1]-lon[ii,nlon-2])/2 ] for ii in range(nlat-1) ]
    dlon = np.array(tmp + [ tmp[-1] ])

    lonh[0:nlat,0:nlon] = lon - dlon
    lath[0:nlat,0:nlon] = lat - dlat

    # make lat and lon one column and row wider such that all
    lonh[nlat,0:nlon] = lonh[nlat-1,0:nlon] + (lonh[nlat-1,0:nlon] - lonh[nlat-2,0:nlon])  
    lath[nlat,0:nlon] = lath[nlat-1,0:nlon] + (lath[nlat-1,0:nlon] - lath[nlat-2,0:nlon])  
    lonh[0:nlat,nlon] = lonh[0:nlat,nlon-1] + (lonh[0:nlat,nlon-1] - lonh[0:nlat,nlon-2])  
    lath[0:nlat,nlon] = lath[0:nlat,nlon-1] + (lath[0:nlat,nlon-1] - lath[0:nlat,nlon-2])  
    lonh[nlat,nlon]   = lonh[nlat-1,nlon-1] + (lonh[nlat-1,nlon-1] - lonh[nlat-2,nlon-2])  
    lath[nlat,nlon]   = lath[nlat-1,nlon-1] + (lath[nlat-1,nlon-1] - lath[nlat-2,nlon-2])

    # geo-referenced
    xx,  yy  = m(lon,lat)
    xxh, yyh = m(lonh,lath)
    zz = var[:,:]  
    var_minval = np.nanmin(var)
    var_maxval = np.nanmax(var)

    # colors
	# RGB (44,162,95)		 1 - Evergreen Needleleaf Forests
    # RGB (153,216,201)		 5 - Mixed Forest
    # RGB (229,245,249)		10 - Grassland
	# RGB (5,113,176)		11 - Water    
    cc = [(44/256.,162/256.,95/256.),(0.2,0.2,0.2),(0.2,0.2,0.2),(0.2,0.2,0.2),(153/256.,216/256.,201/256.),
          (0.2,0.2,0.2),(0.2,0.2,0.2),(0.2,0.2,0.2),(0.2,0.2,0.2),(229/256.,245/256.,249/256.),(5/256.,113/256.,176/256.)]
    cmap = mpl.colors.ListedColormap(cc)

    # plot
    variable_plot = m.pcolor(xxh, yyh, zz, cmap=cmap, zorder = 100, alpha=0.9)  

    # convert coordinates to map coordinates
    xpt,ypt = m(latlon[:,1],latlon[:,0])
    mlatlon = np.transpose(np.array([xpt,ypt])) # zip(xpt,ypt)

    # Catchments
    facecolor = perco_color
    sub.add_patch(Polygon(mlatlon, facecolor='None', linewidth=lwidth, edgecolor='black', alpha=1.0))

    # some reference points
    nref = 4
    ref_latlon = np.ones([nref,2]) * -9999.0
    ref_names  = []
    ref_latlon[0,0] = 49.2827 ; ref_latlon[0,1] = -123.1207  ; ref_names.append( 'Vancouver'     )   # Vancouver
    ref_latlon[1,0] = 53.9171 ; ref_latlon[1,1] = -122.7497  ; ref_names.append( 'Prince George' )   # Prince George
    ref_latlon[2,0] = 52.8737 ; ref_latlon[2,1] = -118.0814  ; ref_names.append( 'Jasper'        )   # Jasper
    ref_latlon[3,0] = 51.1784 ; ref_latlon[3,1] = -115.5708  ; ref_names.append( 'Banff '        )   # Banff

    for iref in range(nref):
            xpt1, ypt1 = m(ref_latlon[iref,1],ref_latlon[iref,0])
            sub.plot(xpt1, ypt1,
                     #linewidth=lwidth,
                     #color=color,
                     marker='o',
                     markeredgecolor='black',
                     markerfacecolor='black',
                     markersize=3.0, #msize,
                     markeredgewidth=mwidth)
            ha = ['right', 'left']
            va = ['bottom','top']

            if iref in [0,2]:
                ha = 'left'
                dx = 0.1
            else:
                ha = 'right'
                dx = -0.1
            x2, y2 = m(ref_latlon[iref,1]+dx,ref_latlon[iref,0])
            sub.annotate(ref_names[iref],
                xy=(xpt1, ypt1),  xycoords='data',
                xytext=(x2, y2), textcoords='data',
                fontsize='small',
                verticalalignment='center',horizontalalignment=ha#,
                # arrowprops=dict(arrowstyle="->",relpos=(1.0,0.5),linewidth=0.4)
                )

    # Title
    #sub.set_title(str2tex('Salmon River watershed',usetex=usetex),fontsize=textsize)

    # percentages
    category_vals = np.unique(var)[~np.isnan(np.unique(var))]
    count_all = np.count_nonzero(~np.isnan(var))
    perc = []
    for icat in category_vals:
        count_cat = np.shape(np.where(var==icat)[0])[0]
        perc.append(count_cat*1.0/count_all*100.)

    print("percentages categories: ",perc)

    # Create custom artists
    #               RGB (44,162,95)		 1 - Evergreen Needleleaf Forests
    #               RGB (153,216,201)	 5 - Mixed Forest
    #               RGB (229,245,249)	10 - Grassland
	#               RGB (5,113,176)		11 - Water   
    #      (left, bottom), width, height
    boxSTi_1  = patches.Rectangle( (1.045, 0.92), 0.05, 0.05, color = cc[0],  alpha=0.6, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_2  = patches.Rectangle( (1.045, 0.81), 0.05, 0.05, color = cc[4],  alpha=0.6, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_3  = patches.Rectangle( (1.045, 0.73), 0.05, 0.05, color = cc[9],  alpha=0.6, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_4  = patches.Rectangle( (1.045, 0.65), 0.05, 0.05, color = cc[10], alpha=0.6, fill  = True, transform=sub.transAxes, clip_on=False )
    sub.add_patch(boxSTi_1)
    sub.add_patch(boxSTi_2)
    sub.add_patch(boxSTi_3)
    sub.add_patch(boxSTi_4)
    sub.text(1.13, 0.94, str2tex("Evergreen Needleleaf \nForest ("+astr(perc[0],prec=0)+"%)",usetex=usetex),
                 fontsize='small', horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)
    sub.text(1.13, 0.83, str2tex("Mixed Forest ("+astr(perc[1],prec=0)+"%)",usetex=usetex),
                 fontsize='small', horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)
    sub.text(1.13, 0.75, str2tex("Grassland ("+astr(perc[2],prec=0)+"%)",usetex=usetex),
                 fontsize='small', horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)
    sub.text(1.13, 0.67, str2tex("Water ("+astr(perc[3],prec=0)+"%)",usetex=usetex),
                 fontsize='small', horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)

    # Fake subplot for numbering
    if doabc:
        pos = position(nrow,ncol,iplot,hspace=hspace,vspace=vspace)
        lsub = fig.add_axes([lleft+wwidth-0.165, bbottom+hheight-0.098, 0.1, 0.1])

        lsub.set_xlim([0,1])
        lsub.set_ylim([0,1])

        # subplot numbering
        abc2plot(lsub, dxabc, dyabc, iplot, lower=False, parenthesis='none',
                     bold=True, large=True,
                     mathrm=True, usetex=usetex,
                     horizontalalignment='left', verticalalignment='top')

        lsub.set_title('')
        lsub.set_xlabel('')
        lsub.set_ylabel('')
        lsub.set_xticks([])
        lsub.set_yticks([])
        lsub.set_axis_off()


    if (outtype == 'pdf'):
        pdf_pages.savefig(fig)
        plt.close(fig)
    elif (outtype == 'png'):
        pngfile = pngbase+"{0:04d}".format(ifig)+".png"
        fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
        plt.close(fig)

    # --------------------------------------
    # Finish
    # --------------------------------------
    if (outtype == 'pdf'):
        pdf_pages.close()
    elif (outtype == 'png'):
        pass
    else:
        plt.show()

    
    t2  = time.time()
    str = '  Time plot [m]: '+astr((t2-t1)/60.,1) if (t2-t1)>60. else '  Time plot [s]: '+astr(t2-t1,0)
    print(str)


