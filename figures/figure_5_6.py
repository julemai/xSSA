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
#     run figure_5_6.py -t pdf -p figure_5_6 -n 10
#     python figure_5_6.py -t pdf -p 08KC001  -i ../data_out/08KC001/results_nsets2.pkl -n 2   -o "pkl"
#     python figure_5_6.py -t pdf -p 08KC001  -i ../data_out/08KC001/results_nsets2.nc  -n 2   -o "nc"

from __future__ import print_function

"""

Plots results of RAVEN sensitivity analysis with multiple process options

History
-------
Written,  JM, Jun 2019
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
    nsets       = 100           # number of Sobol sequences
    nboot       = 1             # Set to 1 for single run of SI and STI calculation
    variable    = 'Q'           # model output variable
    inputfile   = None          # default "results_nsets<nsets>_snow.pkl"
    intype      = None
    
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
    parser.add_argument('-v', '--variable', action='store',
                        default=variable, dest='variable', metavar='variable',
                        help='Model output variable name. Must match key in dictionaries in Pickle. (default: variable="Q").')
    parser.add_argument('-i', '--inputfile', action='store',
                        default=inputfile, dest='inputfile', metavar='inputfile',
                        help="Name of file that contains all model runs and Sobol' indexes. (default: results_nsets<nsets>_snow.pkl).")
    parser.add_argument('-o', '--intype', action='store',
                        default=intype, dest='intype', metavar='intype',
                        help='Type of inputfile: needs to be pkl, msgpack, json, or nc (default: nc).')
    

    args       = parser.parse_args()
    plotname   = args.plotname
    outtype    = args.outtype
    serif      = args.serif
    usetex     = args.usetex
    nboot      = np.int(args.nboot)
    nsets      = np.int(args.nsets)
    variable   = args.variable
    inputfile  = args.inputfile
    intype     = args.intype

    if inputfile is None:
        inputfile = "results_nsets"+str(nsets)+"_snow.pkl"

    if outtype is None:
        outtype = 'pkl'

    if inputfile.split('.')[-1] != intype:
        raise ValueError('Type of input file (option -o) does not match file ending of input file (option -i).')

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
    import copy
    import time
    #import re
    #import os
    import datetime
    #from   raven_templates import RVI, RVT, RVP, RVH, RVC
    #from   raven_common    import writeString, makeDirectories
    #from   pathlib2        import Path
    #import subprocess
    #import shutil
    #from   fread           import fread

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
    # Read results
    # -------------------------------------------------------------------------
    if intype == 'pkl':
        import pickle
        setup         = pickle.load( open( inputfile, "rb" ) )
        sobol_indexes = setup['sobol_indexes']
        
    elif intype == 'nc':
        import netCDF4 as nc4
        nc4_in = nc4.Dataset(inputfile, "r", format="NETCDF4")

        # ---------------------
        # sensitivity indexes: sobol_indexes['paras']['msi'][variable]
        # ---------------------
        sobol_indexes = {}   

        for analysis_type in ['paras', 'process_options', 'processes']:
            
            tmp2 = {}

            # aggregated sensitivity indexes
            for sensi_index_type in ['msi', 'msti', 'wsi', 'wsti']:
                tmp = {}
                ncvar_name = sensi_index_type+'_'+analysis_type.split('_')[-1]      # msi_paras, msi_options, msi_processes
                tmp[variable] = nc4_in.groups[variable].variables[ncvar_name][:]
                tmp2[sensi_index_type] = tmp

            # time-dependent sensitivity indexes
            if analysis_type == 'processes':
                for sensi_index_type in ['sti']:
                    tmp = {}
                    ncvar_name = sensi_index_type+'_'+analysis_type.split('_')[-1]      # msi_paras, msi_options, msi_processes
                    tmp[variable] = nc4_in.groups[variable].variables[ncvar_name][:]
                    tmp2[sensi_index_type] = tmp

            sobol_indexes[analysis_type] = tmp2

        # ---------------------
        # model outouts to derive timesetp weights: setup['f_a'][variable]
        # ---------------------
        setup = {}

        for ncvar_name in ['f_a', 'f_b']:
            tmp = {}
            tmp[variable] = nc4_in.groups[variable].variables[ncvar_name][:]   
            setup[ncvar_name] = tmp      
        
        nc4_in.close()

    # -------------------------------------------------------------------------
    # Colors
    # -------------------------------------------------------------------------
    infil_color = color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[20]   #  [20]
    quick_color = color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[43]   #  [55] 
    evapo_color = color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[66]   #  [80] 
    basef_color = color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[89]   #  [105]
    snowb_color = color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[112]  #  [130]
    convs_color = color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[135]  #  [155]
    convd_color = color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[157]  #  [180]
    potme_color = color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[179]  #  [205]
    perco_color = color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[201]  #  [230]
    rspar_color = color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[223]  #  [230]
    rscor_color = color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[255]  #  [230]
    soilm_color = (0.7,0.7,0.7) #color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[255]
    
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
    from matplotlib.lines import Line2D
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

    unequal_second_row = 0.12
    ifig = 0
    
    # -------------------------------------------------
    # variance-weighted mean Sobol' indexes
    # -------------------------------------------------
    
    ifig += 1
    iplot = 0
    print('Plot - Fig ', ifig)
    fig = plt.figure(ifig)

    # -----------------------
    # plot
    # -----------------------
    ylim = [-0.1, 0.6]

    # [left, bottom, width, height]
    pos_a      = [0.100,0.1750,0.255,0.7750]
    pos_b      = [0.605,0.4825,0.255,0.4675]
    pos_c      = [0.605,0.1750,0.255,0.2275]
    pos_legend = [0.180,0.0200,0.610,0.09]

    # -------------
    # Parameter sensitivities
    # -------------
    iplot += 1

    # [left, bottom, width, height]
    # sub = fig.add_axes(position(nrow, 1, iplot, hspace=hspace/2, vspace=vspace))
    # print("(A) position was: ",position(nrow, 1, iplot, hspace=hspace/2, vspace=vspace))
    
    sub = fig.add_axes(pos_a)
    

    active_snowproc = "+".join(inputfile.split('.')[0].split('_')[4:-1])
    if (active_snowproc == "HMETS_SIMPLE_HBV"):
        paras = [      '$x_{1}$',  '$x_{2}$',  '$x_{3}$',  '$x_{4}$',  '$x_{5}$',  '$x_{6}$',  '$x_{7}$',  '$x_{8}$',  '$x_{9}$',  '$x_{10}$',
                       '$x_{11}$', '$x_{12}$', '$x_{13}$', '$x_{14}$', '$x_{15}$', '$x_{16}$', '$x_{17}$', '$x_{18}$', '$x_{19}$', '$x_{20}$',
                       '$x_{21}$', '$x_{22}$', '$x_{23}$', '$x_{24}$', '$x_{25}$', '$x_{26}$', '$x_{27}$', '$x_{28}$', '$x_{29}$', '$x_{30}$',
                       '$x_{31}$', '$x_{32}$', '$x_{33}$', '$x_{34}$', '$x_{35}$',
                       '$r_{1}$', '$r_{2}$', '$r_{3}$', '$r_{4}$', '$r_{5}$', '$r_{6}$', '$r_{7}$', '$r_{8}$']
    elif (active_snowproc == "HMETS_HBV"):
        paras = [      '$x_{1}$',  '$x_{2}$',  '$x_{3}$',  '$x_{4}$',  '$x_{5}$',  '$x_{6}$',  '$x_{7}$',  '$x_{8}$',  '$x_{9}$',  '$x_{10}$',
                       '$x_{11}$', '$x_{12}$', '$x_{13}$', '$x_{14}$', '$x_{15}$', '$x_{16}$', '$x_{17}$', '$x_{18}$', '$x_{19}$', '$x_{20}$',
                       '$x_{21}$', '$x_{22}$', '$x_{23}$', '$x_{24}$', '$x_{25}$', '$x_{26}$', '$x_{27}$', '$x_{28}$', '$x_{29}$', '$x_{30}$',
                       '$x_{31}$', '$x_{32}$', '$x_{33}$', '$x_{34}$', '$x_{35}$',
                       '$r_{1}$', '$r_{2}$', '$r_{3}$', '$r_{4}$', '$r_{5}$', '$r_{6}$', '$r_{7}$']
    elif (active_snowproc == "HMETS" or active_snowproc == "SIMPLE" or active_snowproc == "HBV"):
        paras = [      '$x_{1}$',  '$x_{2}$',  '$x_{3}$',  '$x_{4}$',  '$x_{5}$',  '$x_{6}$',  '$x_{7}$',  '$x_{8}$',  '$x_{9}$',  '$x_{10}$',
                       '$x_{11}$', '$x_{12}$', '$x_{13}$', '$x_{14}$', '$x_{15}$', '$x_{16}$', '$x_{17}$', '$x_{18}$', '$x_{19}$', '$x_{20}$',
                       '$x_{21}$', '$x_{22}$', '$x_{23}$', '$x_{24}$', '$x_{25}$', '$x_{26}$', '$x_{27}$', '$x_{28}$', '$x_{29}$', '$x_{30}$',
                       '$x_{31}$', '$x_{32}$', '$x_{33}$', '$x_{34}$', '$x_{35}$',
                       '$r_{1}$', '$r_{2}$', '$r_{3}$', '$r_{4}$', '$r_{5}$', '$r_{6}$']
    else:
        # all
        paras = [      '$x_{1}$',  '$x_{2}$',  '$x_{3}$',  '$x_{4}$',  '$x_{5}$',  '$x_{6}$',  '$x_{7}$',  '$x_{8}$',  '$x_{9}$',  '$x_{10}$',
                       '$x_{11}$', '$x_{12}$', '$x_{13}$', '$x_{14}$', '$x_{15}$', '$x_{16}$', '$x_{17}$', '$x_{18}$', '$x_{19}$', '$x_{20}$',
                       '$x_{21}$', '$x_{22}$', '$x_{23}$', '$x_{24}$', '$x_{25}$', '$x_{26}$', '$x_{27}$', '$x_{28}$', '$x_{29}$', '$x_{30}$',
                       '$x_{31}$', '$x_{32}$', '$x_{33}$', '$x_{34}$', '$x_{35}$',
                       '$r_{1}$', '$r_{2}$', '$r_{3}$', '$r_{4}$', '$r_{5}$', '$r_{6}$', '$r_{7}$', '$r_{8}$']
    paras = [ str2tex(ii,usetex=usetex) for ii in paras ]

    keys = sobol_indexes['paras']['wsi'].keys()
    ikey = variable
    if not( variable in keys):
        print("")
        print("Variables in pickle: ",keys)
        print("Variable given:      ",variable)
        raise ValueError("Variable given is not available in Pickle!")

    if (ikey == 'Q'):
        ppp = [4,5,7,15,17,18,23,24,25,26,28,30,33]
        print("sel parameters: ",["x_{"+str(ipp+1)+"}" for ipp in ppp])
        perc = np.sum(sobol_indexes['paras']['wsti'][ikey][ppp])/np.sum(sobol_indexes['paras']['wsti'][ikey][0:35])*100.
        print("Overall sensitivity of sel parameters (w/o ri): ",perc,"%")
        perc = np.sum(sobol_indexes['paras']['wsti'][ikey][ppp])/np.sum(sobol_indexes['paras']['wsti'][ikey][:])*100.
        print("Overall sensitivity of sel parameters (w/  ri): ",perc,"%")
    if (ikey == 'infiltration'):
        ppp = [4,5,7,17,18,23,24,25,26,28,30,33]
        print("sel parameters: ",["x_{"+str(ipp+1)+"}" for ipp in ppp])
        perc = np.sum(sobol_indexes['paras']['wsti'][ikey][ppp])/np.sum(sobol_indexes['paras']['wsti'][ikey][0:35])*100.
        print("Overall sensitivity of sel parameters (w/o ri): ",perc,"%")
        perc = np.sum(sobol_indexes['paras']['wsti'][ikey][ppp])/np.sum(sobol_indexes['paras']['wsti'][ikey][:])*100.
        print("Overall sensitivity of sel parameters (w/  ri): ",perc,"%")
        
    
    npara = np.shape(sobol_indexes['paras']['wsi'][ikey])[0]
    mark1 = sub.barh(np.arange(npara), sobol_indexes['paras']['wsti'][ikey], align='center', alpha=0.3,
                    color=(infil_color, infil_color, infil_color,
                           quick_color, quick_color, quick_color, quick_color,
                           evapo_color, evapo_color, evapo_color,
                           basef_color, basef_color,
                           snowb_color, snowb_color, snowb_color, snowb_color, snowb_color, snowb_color, snowb_color,
                           convs_color,convs_color,
                           convd_color,convd_color,
                           potme_color,potme_color,potme_color,potme_color,
                           perco_color,
                           soilm_color,soilm_color,
                           rspar_color, rspar_color,     # x31, x32
                           rscor_color, rscor_color,     # x33, x34
                           perco_color,                  # x35
                           # weights generating random numbers
                           infil_color, infil_color,
                           quick_color, quick_color,
                           evapo_color,
                           basef_color,
                           snowb_color, snowb_color) ) #, snowb_color) )    # STi wmean
    mark2 = sub.barh(np.arange(npara), sobol_indexes['paras']['wsi'][ikey], align='center', alpha=1.0,
                    color=(infil_color, infil_color, infil_color,
                           quick_color, quick_color, quick_color, quick_color,
                           evapo_color, evapo_color, evapo_color,
                           basef_color, basef_color,
                           snowb_color, snowb_color, snowb_color, snowb_color, snowb_color, snowb_color, snowb_color,
                           convs_color,convs_color,
                           convd_color,convd_color,
                           potme_color,potme_color,potme_color,potme_color,
                           perco_color,
                           soilm_color,soilm_color,
                           rspar_color, rspar_color,     # x31, x32
                           rscor_color, rscor_color,     # x33, x34
                           perco_color,                  # x35
                           # weights generating random numbers
                           infil_color, infil_color,
                           quick_color, quick_color,
                           evapo_color,
                           basef_color,
                           snowb_color, snowb_color) ) #, snowb_color) )    # Si  wmean

    sub.set_ylim([-1,npara])
    sub.set_xlim(ylim)

    npara = len(paras)
    plt.yticks(np.arange(npara), paras,rotation=0,fontsize='medium')
    
    plt.title(str2tex('Sensitivities of\n Model Parameters',usetex=usetex))
    plt.xlabel(str2tex("Sobol' Index",usetex=usetex))

    sub.invert_yaxis()  # labels read top-to-bottom 

    abc2plot(sub,0.95,0.99,iplot,bold=True,usetex=usetex,mathrm=True, large=True, parenthesis='none',verticalalignment='top',horizontalalignment='right')

    # -------------
    # Process option sensitivities
    # -------------
    iplot += 2
    #                  [left, bottom, width, height] 
    # sub = fig.add_axes(position(nrow, ncol, iplot, hspace=hspace/2, vspace=vspace)+[0.0,0.0,unequal_second_row,0.0])
    # print("(B) position was: ",position(nrow, ncol, iplot, hspace=hspace/2, vspace=vspace)+[0.0,0.0,unequal_second_row,0.0])
    
    # [left, bottom, width, height]
    # sub = fig.add_axes(position(nrow, 1, iplot, hspace=hspace/2, vspace=vspace))
    sub = fig.add_axes(pos_b)

    active_snowproc = "_".join(inputfile.split('.')[0].split('_')[4:-1])
    if (active_snowproc == "HMETS_SIMPLE_HBV"):
        procopt = ['INF_HMETS $M_1$', 'INF_VIC_ARNO $M_2$', 'INF_HBV $M_3$',
                   'BASE_LINEAR_ANALYTIC $N_1$', 'BASE_VIC $N_2$', 'BASE_TOPMODEL $N_3$',
                   'SOILEVAP_ALL $O_1$', 'SOILEVAP_TOPMODEL $O_2$',
                   'BASE_LINEAR_ANALYTIC $P_1$', 'BASE_POWER_LAW $P_2$',
                   'SNOBAL_HMETS $Q_1$', 'SNOBAL_SIMPLE_MELT $Q_2$', 'SNOBAL_HBV $Q_3$',
                   'CONVOL_GAMMA $R_1$', 'CONVOL_GAMMA_2 $S_1$', 'POTMELT_HMETS $T_1$', 'PERC_LINEAR $U_1$',
                   'RAINSNOW_HBV $V_1$', 'RAINSNOW_CORRECTION $W_1$',
                   'Infiltration weight gen. $r_{1}$', 'Infiltration weight gen. $r_{2}$',
                   'Quickflow weight gen. $r_{3}$', 'Quickflow weight gen. $r_{4}$',
                   'Evaporation weight gen. $r_{5}$',
                   'Baseflow weight gen. $r_{6}$',
                   'Snow Balance weight gen. $r_{7}$', 'Snow Balance weight gen. $r_{8}$']
        ccc = (snowb_color, snowb_color, snowb_color)
    elif (active_snowproc == "HMETS_HBV"):
        procopt = ['INF_HMETS', 'INF_VIC_ARNO', 'INF_HBV', 'BASE_LINEAR_ANALYTIC', 'BASE_VIC', 'BASE_TOPMODEL',
                   'SOILEVAP_ALL', 'SOILEVAP_TOPMODEL', 'BASE_LINEAR_ANALYTIC', 'BASE_POWER_LAW', 'SNOBAL_HMETS', 'SNOBAL_HBV',
                   'CONVOL_GAMMA', 'CONVOL_GAMMA_2', 'POTMELT_HMETS', 'PERC_LINEAR',
                   'RAINSNOW_HBV $V_1$', 'RAINSNOW_CORRECTION $W_1$',
                   'Infiltration weight gen. $r_{1}$', 'Infiltration weight gen. $r_{2}$', 'Quickflow weight gen. $r_{3}$', 'Quickflow weight gen. $r_{4}$',
                   'Evaporation weight gen. $r_{5}$', 'Baseflow weight gen. $r_{6}$', 'Snow Balance weight gen. $r_{7}$']
        ccc = (snowb_color, snowb_color)
    elif (active_snowproc == "HMETS" or active_snowproc == "SIMPLE" or active_snowproc == "HBV"):
        procopt = ['INF_HMETS', 'INF_VIC_ARNO', 'INF_HBV', 'BASE_LINEAR_ANALYTIC', 'BASE_VIC', 'BASE_TOPMODEL',
                   'SOILEVAP_ALL', 'SOILEVAP_TOPMODEL', 'BASE_LINEAR_ANALYTIC', 'BASE_POWER_LAW', 'SNOBAL_'+active_snowproc,
                   'CONVOL_GAMMA', 'CONVOL_GAMMA_2', 'POTMELT_HMETS', 'PERC_LINEAR',
                   'RAINSNOW_HBV $V_1$', 'RAINSNOW_CORRECTION $W_1$',
                   'Infiltration weight gen. $r_{1}$', 'Infiltration weight gen. $r_{2}$', 'Quickflow weight gen. $r_{3}$', 'Quickflow weight gen. $r_{4}$',
                   'Evaporation weight gen. $r_{5}$', 'Baseflow weight gen. $r_{6}$']
        ccc = (snowb_color)
        ccc = (ccc,)
    else:
        # all
        procopt = ['INF_HMETS $M_1$', 'INF_VIC_ARNO $M_2$', 'INF_HBV $M_3$',
                   'BASE_LINEAR_ANALYTIC $N_1$', 'BASE_VIC $N_2$', 'BASE_TOPMODEL $N_3$',
                   'SOILEVAP_ALL $O_1$', 'SOILEVAP_TOPMODEL $O_2$',
                   'BASE_LINEAR_ANALYTIC $P_1$', 'BASE_POWER_LAW $P_2$',
                   'SNOBAL_HMETS $Q_1$', 'SNOBAL_SIMPLE_MELT $Q_2$', 'SNOBAL_HBV $Q_3$',
                   'CONVOL_GAMMA $R_1$', 'CONVOL_GAMMA_2 $S_1$', 'POTMELT_HMETS $T_1$', 'PERC_LINEAR $U_1$',
                   'RAINSNOW_HBV $V_1$', 'RAINSNOW_CORRECTION $W_1$',
                   'Infiltration weight gen. $r_{1}$', 'Infiltration weight gen. $r_{2}$', 'Quickflow weight gen. $r_{3}$', 'Quickflow weight gen. $r_{4}$',
                   'Evaporation weight gen. $r_{5}$', 'Baseflow weight gen. $r_{6}$', 'Snow Balance weight gen. $r_{7}$', 'Snow Balance weight gen. $r_{8}$']
        ccc = (snowb_color, snowb_color, snowb_color)

    keys = sobol_indexes['process_options']['wsi'].keys()
    ikey = variable
    if not( variable in keys):
        print("")
        print("Variables in pickle: ",keys)
        print("Variable given:      ",variable)
        raise ValueError("Variable given is not available in Pickle!")

    if ikey == 'Q':
        iprocopt = 'P_1'
        sublist = [s for s in procopt if iprocopt in s]
        iiii = procopt.index(sublist[0])
        print(sublist[0]+' :: ST^w_{'+iprocopt+'} = ',sobol_indexes['process_options']['wsti'][ikey][iiii])

        iprocopt = 'P_2'
        sublist = [s for s in procopt if iprocopt in s]
        iiii = procopt.index(sublist[0])
        print(sublist[0]+' :: ST^w_{'+iprocopt+'} = ',sobol_indexes['process_options']['wsti'][ikey][iiii])

        iprocopt = 'R_1'
        sublist = [s for s in procopt if iprocopt in s]
        iiii = procopt.index(sublist[0])
        print(sublist[0]+' :: ST^w_{'+iprocopt+'} = ',sobol_indexes['process_options']['wsti'][ikey][iiii])

        iprocopt = 'S_1'
        sublist = [s for s in procopt if iprocopt in s]
        iiii = procopt.index(sublist[0])
        print(sublist[0]+' :: ST^w_{'+iprocopt+'} = ',sobol_indexes['process_options']['wsti'][ikey][iiii])

    procopt = [ str2tex(ii,usetex=usetex) for ii in procopt ]
    
    nopt = np.shape(sobol_indexes['process_options']['wsi'][ikey])[0]
    mark1 = sub.barh(np.arange(nopt), sobol_indexes['process_options']['wsti'][ikey], align='center', alpha=0.6,
                    color=(infil_color, infil_color, infil_color,
                           quick_color, quick_color, quick_color,
                           evapo_color, evapo_color,
                           basef_color, basef_color)+ccc+(
                           convs_color, convd_color, potme_color, perco_color,
                           rspar_color, rscor_color,
                           infil_color, infil_color,
                           quick_color, quick_color,
                           evapo_color,
                           basef_color,
                           snowb_color, snowb_color, snowb_color))    # STi wmean
    mark2 = sub.barh(np.arange(nopt), sobol_indexes['process_options']['wsi'][ikey], align='center', alpha=1.0,
                    color=(infil_color, infil_color, infil_color,
                           quick_color, quick_color, quick_color,
                           evapo_color, evapo_color,
                           basef_color, basef_color)+ccc+(
                           convs_color, convd_color, potme_color, perco_color,
                           rspar_color, rscor_color,
                           infil_color, infil_color,
                           quick_color, quick_color,
                           evapo_color,
                           basef_color,
                           snowb_color, snowb_color, snowb_color))        # Si  wmean

    sub.set_ylim([-1,nopt])
    sub.set_xlim(ylim)

    nopt = len(procopt)
    plt.yticks(np.arange(nopt), procopt,rotation=0,fontsize='medium')
    
    plt.title(str2tex('Sensitivities of\n Process Options',usetex=usetex))
    # plt.xlabel(str2tex("Sobol' Index",usetex=usetex))

    sub.invert_yaxis()  # labels read top-to-bottom 

    abc2plot(sub,0.95,0.985,iplot-1,bold=True,usetex=usetex,mathrm=True, large=True, parenthesis='none',verticalalignment='top',horizontalalignment='right')

    # -------------
    # Process sensitivities
    # -------------
    iplot += 1
    #                  [left, bottom, width, height] 
    # sub = fig.add_axes(position(nrow, ncol, iplot, hspace=hspace/2, vspace=vspace)+[unequal_second_row,0.0,-unequal_second_row,0.0])
    # print("(C) position was: ",position(nrow, ncol, iplot, hspace=hspace/2, vspace=vspace)+[unequal_second_row,0.0,-unequal_second_row,0.0])

    # [left, bottom, width, height]
    # sub = fig.add_axes(position(nrow, 1, iplot, hspace=hspace/2, vspace=vspace))
    sub = fig.add_axes(pos_c)

    processes = ['Infiltration $M$','Quickflow $N$','Evaporation $O$','Baseflow $P$','Snow Balance $Q$', 'Convolution (srfc runoff) $R$',
                 'Convolution (dlyd runoff) $S$', 'Potential Melt $T$', 'Percolation $U$',
                 'Rain-Snow Partitioning $V$', 'Precipitation Correction $W$']  # , 'Soil model'
    processes = [ str2tex(ii,usetex=usetex) for ii in processes ]

    keys = sobol_indexes['processes']['wsi'].keys()
    ikey = variable
    if not( variable in keys):
        print("")
        print("Variables in pickle: ",keys)
        print("Variable given:      ",variable)
        raise ValueError("Variable given is not available in Pickle!")

    if (ikey == 'Q'):
        ppp = [1,4]
        print("Surface:: sel processes: ",[chr(ord('M')+ipp) for ipp in ppp])
        perc = np.sum(sobol_indexes['processes']['wsti'][ikey][ppp])/np.sum(sobol_indexes['processes']['wsti'][ikey][:])*100.
        print("Overall sensitivity of sel processes: ",perc,"%")

        ppp = [3,5,6]
        print("Deep:: sel processes: ",[chr(ord('M')+ipp) for ipp in ppp])
        perc = np.sum(sobol_indexes['processes']['wsti'][ikey][ppp])/np.sum(sobol_indexes['processes']['wsti'][ikey][:])*100.
        print("Overall sensitivity of sel processes: ",perc,"%")

        ppp = [0,2,8]
        print("Soil:: sel processes: ",[chr(ord('M')+ipp) for ipp in ppp])
        perc = np.sum(sobol_indexes['processes']['wsti'][ikey][ppp])/np.sum(sobol_indexes['processes']['wsti'][ikey][:])*100.
        print("Overall sensitivity of sel processes: ",perc,"%")

        ppp = [7,9,10]
        print("Input:: sel processes: ",[chr(ord('M')+ipp) for ipp in ppp])
        perc = np.sum(sobol_indexes['processes']['wsti'][ikey][ppp])/np.sum(sobol_indexes['processes']['wsti'][ikey][:])*100.
        print("Overall sensitivity of sel processes: ",perc,"%")
    
    nproc = np.shape(sobol_indexes['processes']['wsi'][ikey])[0]    
    mark1 = sub.barh(    np.arange(nproc), sobol_indexes['processes']['wsti'][ikey], align='center', alpha=0.6,
                        color=(infil_color, quick_color, evapo_color, basef_color, snowb_color,
                               convs_color, convd_color, potme_color, perco_color,
                               rspar_color, rscor_color)) #, soilm_color))    # STi wmean
    mark2 = sub.barh(    np.arange(nproc), sobol_indexes['processes']['wsi'][ikey], align='center', alpha=1.0,
                        color=(infil_color, quick_color, evapo_color, basef_color, snowb_color,
                               convs_color, convd_color, potme_color, perco_color,
                               rspar_color, rscor_color)) #, soilm_color))    # Si  wmean

    sub.set_ylim([-1,nproc])
    sub.set_xlim(ylim)
    
    plt.yticks(np.arange(nproc), processes,rotation=0,fontsize='medium')
    plt.title(str2tex('Sensitivities of\n Processes',usetex=usetex))
    plt.xlabel(str2tex("Sobol' Index",usetex=usetex))

    sub.invert_yaxis()  # labels read top-to-bottom

    abc2plot(sub,0.95,0.97,iplot-1,bold=True,usetex=usetex,mathrm=True, large=True, parenthesis='none',verticalalignment='top',horizontalalignment='right')


    # -------------
    # Legend
    # -------------
    #                  [left, bottom, width, height] 
    # sub = fig.add_axes(position(nrow, ncol, iplot, hspace=hspace/2, vspace=vspace)+[unequal_second_row/2,-0.04,-unequal_second_row,0.02])
    # print('Legend position was: ',position(nrow, ncol, iplot, hspace=hspace/2, vspace=vspace)+[unequal_second_row/2,-0.04,-unequal_second_row,0.02])
    sub = fig.add_axes(pos_legend)

    sub.axis('off')
    # Create custom artists
    #      (left, bottom), width, height
    boxSi_1   = patches.Rectangle( (0.00, 0.90), 0.03, 0.12, color = infil_color, alpha=1.0, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSi_2   = patches.Rectangle( (0.04, 0.90), 0.03, 0.12, color = quick_color, alpha=1.0, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSi_3   = patches.Rectangle( (0.08, 0.90), 0.03, 0.12, color = evapo_color, alpha=1.0, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSi_4   = patches.Rectangle( (0.12, 0.90), 0.03, 0.12, color = basef_color, alpha=1.0, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSi_5   = patches.Rectangle( (0.16, 0.90), 0.03, 0.12, color = snowb_color, alpha=1.0, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSi_6   = patches.Rectangle( (0.20, 0.90), 0.03, 0.12, color = convs_color, alpha=1.0, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSi_7   = patches.Rectangle( (0.24, 0.90), 0.03, 0.12, color = convd_color, alpha=1.0, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSi_8   = patches.Rectangle( (0.28, 0.90), 0.03, 0.12, color = potme_color, alpha=1.0, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSi_9   = patches.Rectangle( (0.32, 0.90), 0.03, 0.12, color = perco_color, alpha=1.0, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSi_10  = patches.Rectangle( (0.36, 0.90), 0.03, 0.12, color = rspar_color, alpha=1.0, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSi_11  = patches.Rectangle( (0.40, 0.90), 0.03, 0.12, color = rscor_color, alpha=1.0, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSi_12  = patches.Rectangle( (0.44, 0.90), 0.03, 0.12, color = soilm_color, alpha=1.0, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_1  = patches.Rectangle( (0.00, 0.64), 0.03, 0.12, color = infil_color, alpha=0.4, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_2  = patches.Rectangle( (0.04, 0.64), 0.03, 0.12, color = quick_color, alpha=0.4, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_3  = patches.Rectangle( (0.08, 0.64), 0.03, 0.12, color = evapo_color, alpha=0.4, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_4  = patches.Rectangle( (0.12, 0.64), 0.03, 0.12, color = basef_color, alpha=0.4, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_5  = patches.Rectangle( (0.16, 0.64), 0.03, 0.12, color = snowb_color, alpha=0.4, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_6  = patches.Rectangle( (0.20, 0.64), 0.03, 0.12, color = convs_color, alpha=0.4, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_7  = patches.Rectangle( (0.24, 0.64), 0.03, 0.12, color = convd_color, alpha=0.4, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_8  = patches.Rectangle( (0.28, 0.64), 0.03, 0.12, color = potme_color, alpha=0.4, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_9  = patches.Rectangle( (0.32, 0.64), 0.03, 0.12, color = perco_color, alpha=0.4, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_10 = patches.Rectangle( (0.36, 0.64), 0.03, 0.12, color = rspar_color, alpha=0.4, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_11 = patches.Rectangle( (0.40, 0.64), 0.03, 0.12, color = rscor_color, alpha=0.4, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_12 = patches.Rectangle( (0.44, 0.64), 0.03, 0.12, color = soilm_color, alpha=0.4, fill  = True, transform=sub.transAxes, clip_on=False )
    sub.add_patch(boxSi_1)  ;  sub.add_patch(boxSTi_1)
    sub.add_patch(boxSi_2)  ;  sub.add_patch(boxSTi_2)
    sub.add_patch(boxSi_3)  ;  sub.add_patch(boxSTi_3)
    sub.add_patch(boxSi_4)  ;  sub.add_patch(boxSTi_4)
    sub.add_patch(boxSi_5)  ;  sub.add_patch(boxSTi_5)
    sub.add_patch(boxSi_6)  ;  sub.add_patch(boxSTi_6)
    sub.add_patch(boxSi_7)  ;  sub.add_patch(boxSTi_7)
    sub.add_patch(boxSi_8)  ;  sub.add_patch(boxSTi_8)
    sub.add_patch(boxSi_9)  ;  sub.add_patch(boxSTi_9)
    sub.add_patch(boxSi_10) ;  sub.add_patch(boxSTi_10)
    sub.add_patch(boxSi_11) ;  sub.add_patch(boxSTi_11)
    sub.add_patch(boxSi_12) ;  sub.add_patch(boxSTi_12)
    
    sub.text(0.50, 0.92, str2tex("Sobol' main effect $\overline{S_i^w}$",usetex=usetex), horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)
    sub.text(0.50, 0.66, str2tex("Sobol' total effect $\overline{ST_i^w}$",usetex=usetex), horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)


    if (outtype == 'pdf'):
        pdf_pages.savefig(fig)
        plt.close(fig)
    elif (outtype == 'png'):
        pngfile = pngbase+"{0:04d}".format(ifig)+".png"
        fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
        plt.close(fig)



    # -------------------------------------------------
    # Sobol' indexes of all processes over time (stacked)
    # -------------------------------------------------
    
    ifig += 1
    iplot = 0
    print('Plot - Fig ', ifig)
    fig = plt.figure(ifig)

    # -----------------------
    # plot
    # -----------------------
    ylim = [0.0, 1.0]

    # -------------
    # Parameter sensitivities
    # -------------
    iplot += 1
    #                  [left, bottom, width, height] 
    sub = fig.add_axes(position(nrow, 1, iplot, hspace=hspace/2, vspace=vspace)-[0.02,0,0,0])

    keys = sobol_indexes['processes']['sti'].keys()
    ikey = variable
    if not( variable in keys):
        print("")
        print("Variables in pickle: ",keys)
        print("Variable given:      ",variable)
        raise ValueError("Variable given is not available in Pickle!")

    ntime = np.shape(sobol_indexes['processes']['sti'][ikey])[0]
    nproc = np.shape(sobol_indexes['processes']['sti'][ikey])[1]
    ntime_doy = 365

    # determine leap days and discrad them
    start_date = datetime.datetime(1991,1,1,0,0)
    times = np.array([start_date+datetime.timedelta(ii) for ii in range(ntime)])
    leap = np.array([ True if (times[ii].month == 2 and times[ii].day == 29) else False for ii in range(ntime) ])
    # leap[-1] = True  # by accident last day is 2011-01-01 instaed of 2010-12-31

    # weights used for STi in 'Mai1999'
    varA    = np.var(setup['f_a'][ikey][~leap,:],axis=1)
    denomA  = 1./np.sum(varA)
    weights = varA*denomA # (ntime, nproc)

    # average hydrograph
    median_Q = np.median(setup['f_a'][ikey][~leap,:],axis=1)

    # reshape such that (ntime,nproc) --> (nyears, 365, nproc)
    tmp_sobol    = copy.deepcopy(sobol_indexes['processes']['sti'][ikey][~leap,:])     # (ntime, nproc)
    tmp_weights  = copy.deepcopy(weights)  # (ntime)
    tmp_median_Q = copy.deepcopy(median_Q)   # (ntime)
    sobol_doy    = np.ones([int(ntime/ntime_doy),ntime_doy,nproc]) * -9999.0
    weights_doy  = np.ones([int(ntime/ntime_doy),ntime_doy]) * -9999.0
    median_Q_doy = np.ones([int(ntime/ntime_doy),ntime_doy]) * -9999.0
    for iproc in range(nproc):
        sobol_doy[:,:,iproc]   = np.reshape(tmp_sobol[:,iproc],  [int(ntime/ntime_doy),ntime_doy])
        # sobol_doy[:,:,iproc] = np.reshape(tmp_sobol[:,iproc],  [int(ntime/ntime_doy),ntime_doy])
        sobol_doy[:,:,iproc]   = np.where(np.isinf(np.reshape(tmp_sobol[:,iproc],  [int(ntime/ntime_doy),ntime_doy])),np.nan,np.reshape(tmp_sobol[:,iproc],  [int(ntime/ntime_doy),ntime_doy]))
    weights_doy[:,:]  = np.reshape(tmp_weights[:],  [int(ntime/ntime_doy),ntime_doy])
    median_Q_doy[:,:] = np.reshape(tmp_median_Q[:], [int(ntime/ntime_doy),ntime_doy])

    # average over years
    sobol_doy_mean    = np.nanmean(sobol_doy,axis=0)
    weights_doy_mean  = np.nanmean(weights_doy,axis=0)
    median_Q_doy_mean = np.nanmean(median_Q_doy,axis=0)

    # to scale Sobol' indexes such that they sum up to 1.0 (over all processes)
    csum   = np.sum(sobol_doy_mean,axis=1)

    # colors for all processes
    colors = [infil_color, quick_color, evapo_color, basef_color, snowb_color, convs_color, convd_color, potme_color, perco_color, rspar_color, rscor_color, soilm_color]

    
    width = 1.0
    for iproc in range(nproc):
        p1 = sub.bar(np.arange(ntime_doy),
                         sobol_doy_mean[:,iproc]/csum,
                         width,
                         color=colors[iproc],
                         bottom=np.sum(sobol_doy_mean[:,0:iproc],axis=1)/csum)

    ff = open('.'.join(inputfile.split('.')[:-1])+'_wSTi_processes.csv',"w")
    # ff.write('# CSV of process sensitivities over time\n')
    # ff.write('# Original file: '+inputfile+'\n')
    ff.write('date, '+', '.join(processes[::-1])+'\n')
    for itime in range(ntime_doy):
        ff.write((start_date+datetime.timedelta(days=itime)).strftime("%Y-%m-%d")+', '+', '.join(astr(sobol_doy_mean[itime,::-1]/csum[itime],prec=5))+'\n')
    ff.close()

    # ff = open('.'.join(inputfile.split('.')[:-1])+'_hydrograph.csv',"w")
    # ff.write('# Mean over 20 years of median hydrograph of all model runs sets A\n')
    # ff.write('# Original file: '+inputfile+'\n')
    # ff.write('date, median_Q_doy_mean \n')
    # for itime in range(ntime_doy):
    #     ff.write((start_date+datetime.timedelta(days=itime)).strftime("%Y-%m-%d")+', '+astr(median_Q_doy_mean[itime],prec=5)+'\n')
    # ff.close()

    ff = open('.'.join(inputfile.split('.')[:-1])+'_wSTi_weights.csv',"w")
    #ff.write('# Mean over 20 years of median hydrograph of all model runs sets A\n')
    #ff.write('# Original file: '+inputfile+'\n')
    ff.write('date, weight \n')
    for itime in range(ntime_doy):
        ff.write((start_date+datetime.timedelta(days=itime)).strftime("%Y-%m-%d")+', '+astr(weights_doy_mean[itime],prec=5)+'\n')
    ff.close()

    day1 = 150
    day2 = 300
    print("Average sensitivity of processes druing summer (DOY "+str(day1)+"-"+str(day2)+"):")
    for iproc in range(nproc):
        print(ikey," :: ",iproc," :: ",np.mean((sobol_doy_mean[:,iproc]/csum)[day1:day2])*100.0,"%")

    sub2 = sub.twinx()
    sub2.plot(np.arange(ntime_doy),weights_doy_mean,color='black',linewidth=lwidth*2)
    sub2.set_ylabel(str2tex('Weight',usetex=usetex), color='black')

    # only for month ticks instead of DOY
    sub3 = sub.twiny()
    #sub3.plot(np.arange(ntime_doy),weights_doy_mean,color='black',linewidth=lwidth*0.)  # second time plot (fake)
    monthlength = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    sub3.set_xlim([0,365])
    sub3.set_xlabel(str2tex('Months',usetex=usetex), color='black')
    sub3.set_xticks(np.cumsum(monthlength)) #, ['J','F','M','A','M','J','J','A','S','O','N','D'])
    sub3.set_xticklabels('')
    sub3.set_xticks((monthlength*1.0/2.)+np.cumsum(np.append(0,monthlength)[0:12]), minor=True)
    sub3.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'], minor=True)
    sub3.tick_params(which='minor',length=0)   # dont show minor ticks; only labels
    #plt.xticks(np.cumsum(np.array([31,28,31,30,31,30,31,31,30,31,30,31])), ['J','F','M','A','M','J','J','A','S','O','N','D'])

    sub.set_xlim([0,ntime_doy])
    sub.set_ylim(ylim)

    #npara = len(paras)
    #plt.xticks(np.arange(npara), paras,rotation=90,fontsize='small')
    
    # sub.set_title(str2tex('Sensitivities of Model Processes',usetex=usetex))
    sub.set_xlabel(str2tex("Day of Year",usetex=usetex))
    sub.set_ylabel(str2tex("(normalized) Total\n Sobol' Index $ST_i$",usetex=usetex))

    # Create custom artists
    #      (left, bottom), width, height
    boxSTi_1  = patches.Rectangle( (0.00, -0.54), 0.02, 0.05, color = infil_color, alpha=0.6, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_2  = patches.Rectangle( (0.00, -0.63), 0.02, 0.05, color = quick_color, alpha=0.6, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_3  = patches.Rectangle( (0.00, -0.72), 0.02, 0.05, color = evapo_color, alpha=0.6, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_4  = patches.Rectangle( (0.00, -0.81), 0.02, 0.05, color = basef_color, alpha=0.6, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_5  = patches.Rectangle( (0.22, -0.54), 0.02, 0.05, color = snowb_color, alpha=0.6, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_6  = patches.Rectangle( (0.22, -0.63), 0.02, 0.05, color = convs_color, alpha=0.6, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_7  = patches.Rectangle( (0.22, -0.72), 0.02, 0.05, color = convd_color, alpha=0.6, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_8  = patches.Rectangle( (0.22, -0.81), 0.02, 0.05, color = potme_color, alpha=0.6, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_9  = patches.Rectangle( (0.62, -0.54), 0.02, 0.05, color = perco_color, alpha=0.6, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_10 = patches.Rectangle( (0.62, -0.63), 0.02, 0.05, color = rspar_color, alpha=0.6, fill  = True, transform=sub.transAxes, clip_on=False )
    boxSTi_11 = patches.Rectangle( (0.62, -0.72), 0.02, 0.05, color = rscor_color, alpha=0.6, fill  = True, transform=sub.transAxes, clip_on=False )
    #boxSTi_10 = patches.Rectangle( (0.22, -0.86), 0.02, 0.05, color = soilm_color, alpha=0.6, fill  = True, transform=sub.transAxes, clip_on=False )
    line      = patches.Rectangle( (0.62, -0.79), 0.02, 0.00, color = 'black',     alpha=1.0, fill  = True, transform=sub.transAxes, clip_on=False )
    sub.add_patch(boxSTi_1)
    sub.add_patch(boxSTi_2)
    sub.add_patch(boxSTi_3)
    sub.add_patch(boxSTi_4)
    sub.add_patch(boxSTi_5)
    sub.add_patch(boxSTi_6)
    sub.add_patch(boxSTi_7)
    sub.add_patch(boxSTi_8)
    sub.add_patch(boxSTi_9)
    sub.add_patch(boxSTi_10)
    sub.add_patch(boxSTi_11)
    #sub.add_patch(boxSTi_12)
    sub.add_patch(line)

    sub.text(0.00, -0.40, str2tex("Processes:",usetex=usetex),                   horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)
    sub.text(0.04, -0.53, str2tex("Infiltration $M$",usetex=usetex),                 fontsize='medium', horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)
    sub.text(0.04, -0.62, str2tex("Quickflow $N$",usetex=usetex),                    fontsize='medium', horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)
    sub.text(0.04, -0.71, str2tex("Evaporation $O$",usetex=usetex),                  fontsize='medium', horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)
    sub.text(0.04, -0.80, str2tex("Baseflow $P$",usetex=usetex),                     fontsize='medium', horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)
    sub.text(0.26, -0.53, str2tex("Snow Balance $Q$",usetex=usetex),                 fontsize='medium', horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)
    sub.text(0.26, -0.62, str2tex("Convolution (Surface Runoff) $R$",usetex=usetex), fontsize='medium', horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)
    sub.text(0.26, -0.71, str2tex("Convolution (Delayed Runoff) $S$",usetex=usetex), fontsize='medium', horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)
    sub.text(0.26, -0.80, str2tex("Potential Melt $T$",usetex=usetex),               fontsize='medium', horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)
    sub.text(0.66, -0.53, str2tex("Percolation $U$",usetex=usetex),                  fontsize='medium', horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)
    sub.text(0.66, -0.62, str2tex("Rain-Snow Partitioning $V$",usetex=usetex),       fontsize='medium', horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)
    sub.text(0.66, -0.71, str2tex("Precipitation Correction $W$",usetex=usetex),     fontsize='medium', horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)
    #sub.text(0.26, -0.85, str2tex("Soil Model",usetex=usetex),                      fontsize='medium', horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)
    sub.text(0.66, -0.80, str2tex("Weight of Timestep",usetex=usetex),               fontsize='medium', horizontalalignment='left', verticalalignment='center', transform=sub.transAxes)
    
    #abc2plot(sub,dxabc/2,dyabc,iplot,bold=True,usetex=usetex,mathrm=True, large=True, parenthesis='none',verticalalignment='top')

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


