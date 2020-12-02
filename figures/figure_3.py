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
#     run figure_3.py -t pdf -p figure_3

from __future__ import print_function

"""

Plots results of benchmark function compared to theoretical indexes.
Either 'shared' or 'disjoint' benchmark model can be chosen.

Results in Pickle file are created by running:
    python __run_benchmark_DVM.py -o "benchmark_model_runs_DVM.pkl"
    python __run_benchmark_xSSA.py -o "benchmark_model_runs_xSSA.pkl"

with 
    do_disjoint_example_DVM  = True    in  __run_benchmark_DVM.py
    do_realistic_example_DVM = True    in  __run_benchmark_DVM.py
    do_disjoint_example  = True        in  __run_benchmark_xSSA.py
    do_shared_example    = True        in  __run_benchmark_xSSA.py

History
-------
Written,  JM, Apr 2020
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
    picklefile  = None          # default "results_shared-benchmark-model_DVM.pkl"
    model       = 'shared'
    
    parser   = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                      description='''Benchmark example to test Sensitivity Analysis for models with multiple process options.''')
    parser.add_argument('-p', '--plotname', action='store',
                        default=plotname, dest='plotname', metavar='plotname',
                        help='Name of plot output file for types pdf, html or d3, '
                        'and name basis for type png (default: '+__file__[0:__file__.rfind(".")]+').')
    parser.add_argument('-t', '--type', action='store',
                        default=outtype, dest='outtype', metavar='outtype',
                        help='Output type is pdf, png, html, or d3 (default: open screen windows).')
    parser.add_argument('-u', '--usetex', action='store_true', default=usetex, dest="usetex",
                        help="Use LaTeX to render text in pdf, png and html.")
    parser.add_argument('-s', '--serif', action='store_true', default=serif, dest="serif",
                    help="Use serif font; default sans serif.")
    parser.add_argument('-i', '--picklefile', action='store',
                        default=picklefile, dest='picklefile', metavar='picklefile',
                        help="Name of Pickle that contains all model runs and Sobol' indexes. (default: results_nsets<nsets>.pkl).")
    parser.add_argument('-m', '--model', action='store',
                        default=model, dest='model', metavar='model',
                        help="Benchmark model. Either 'disjoint' or 'shared'. (default: 'shared').")

    args       = parser.parse_args()
    plotname   = args.plotname
    outtype    = args.outtype
    usetex     = args.usetex
    serif      = args.serif
    picklefile = args.picklefile
    model      = args.model


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
    import pickle
    import time

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
    # Colors
    # -------------------------------------------------------------------------
    colors_benchmark = ( #color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[20], 
                         color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[55],
                         color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[80],
                         color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[105],
                         #color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[130],
                         color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[155],
                         #color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[180],
                         color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[205],
                         #color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[230],
                         #color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[230],
                         color.get_brewer( 'WhiteBlueGreenYellowRed',rgb=True)[255],
                         (0.7,0.7,0.7) )
    
    # -------------------------------------------------------------------------
    # Plotting of results
    # -------------------------------------------------------------------------
    # Main plot
    ncol        = 2           # number columns
    nrow        = 5           # number of rows
    textsize    = 10          # standard text size
    dxabc       = 0.60          # % of (max-min) shift to the right from left y-axis for a,b,c,... labels
    dyabc       = 0.95          # % of (max-min) shift up from lower x-axis for a,b,c,... labels
    dxsig       = 1.23        # % of (max-min) shift to the right from left y-axis for signature
    dysig       = -0.075      # % of (max-min) shift up from lower x-axis for signature
    dxtit       = 0           # % of (max-min) shift to the right from left y-axis for title
    dytit       = 1.2         # % of (max-min) shift up from lower x-axis for title
    hspace      = 0.12        # x-space between subplots
    vspace      = 0.03        # y-space between subplots

    lwidth      = 1.0         # linewidth
    elwidth     = 0.5         # errorbar line width
    alwidth     = 1.0         # axis line width
    glwidth     = 0.5         # grid line width
    msize       = 8.0         # marker size
    mwidth      = 1.0         # marker edge width
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
    llxbbox     = 1.02        # x-anchor legend bounding box
    llybbox     = -0.08        # y-anchor legend bounding box
    llrspace    = 0.          # spacing between rows in legend
    llcspace    = 0.5         # spacing between columns in legend
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

    
    ifig = 0


    ifig += 1
    iplot = 0
    print('Plot - Fig ', ifig)
    fig = plt.figure(ifig)

        
    for model in ['disjoint', 'shared']:

        for method in ['DVM', 'xSSA']:

            print("Model: ",model)
    
            # -------------------------------------------------------------------------
            # Set processes and options
            # -------------------------------------------------------------------------
            if model == 'disjoint':
                paras = ['$x_{1}$', '$x_{2}$', '$x_{3}$', '$x_{4}$', '$x_{5}$', '$x_{6}$', '$x_{7}$', '$r_{1}$', '$r_{2}$', '$r_{3}$', '$r_{4}$' ]
                paras = [ str2tex(ii,usetex=usetex) for ii in paras ]
            
                procopt = ['$A_{1}$', '$A_{2}$', '$B_{1}$', '$B_{2}$', '$B_{3}$', '$C_{1}$', '$C_{2}$', '$r_{1}$', '$r_{2}$', '$r_{3}$', '$r_{4}$' ]
                procopt = [ str2tex(ii,usetex=usetex) for ii in procopt ]
                
                processes = ['$A$', '$B$', '$C$' ]
                processes = [ str2tex(ii,usetex=usetex) for ii in processes ]

            elif model == 'shared':
                paras = ['$x_{1}$', '$x_{2}$', '$x_{3}$', '$x_{4}$', '$x_{5}$', '$x_{6}$', '$x_{7}$', '$r_{1}$', '$r_{2}$', '$r_{3}$', '$r_{4}$' ]
                paras = [ str2tex(ii,usetex=usetex) for ii in paras ]
            
                procopt = ['$D_{1}$', '$D_{2}$', '$E_{1}$', '$E_{2}$', '$E_{3}$', '$F_{1}$', '$F_{2}$', '$r_{1}$', '$r_{2}$', '$r_{3}$', '$r_{4}$' ]
                procopt = [ str2tex(ii,usetex=usetex) for ii in procopt ]
                
                processes = ['$D$', '$E$', '$F$' ]
                processes = [ str2tex(ii,usetex=usetex) for ii in processes ]
            else:
                raise ValueError('Model (option -m) needs to be either "disjoint" or "shared".')

            # -------------------------------------------------------------------------
            # Read results
            # -------------------------------------------------------------------------
            if method == 'DVM' or method == 'xSSA':
                picklefile = "../examples/data_out/benchmark/results_"+model+"-benchmark-model_"+method+".pkl"
            else:
                raise ValueError('Method needs to be either "DVM" or "xSSA".')

                
            import pickle
            setup         = pickle.load( open( picklefile, "rb" ) )
            if model == 'shared':
                numerical_results = setup['numerical_si_sti_shared']
            elif model == 'disjoint':
                numerical_results = setup['numerical_si_sti_disjoint']
            else:
                raise ValueError('Model (option -m) needs to be either "disjoint" or "shared".')
            
            # theoretical sensitivity indexes (converting to same structure as numerical results)
            theoretical_results = {}
            theoretical_results['processes']              = {}
            if model == 'shared':
                theoretical_results['processes']['si']        = setup['theo_si_sti_shared']['wa_A-wb_B-wc_C_process'][0]
                theoretical_results['processes']['sti']       = setup['theo_si_sti_shared']['wa_A-wb_B-wc_C_process'][1]
            elif model == 'disjoint':
                theoretical_results['processes']['si']        = setup['theo_si_sti_disjoint']['wa_A-wb_B-wc_C_process'][0]
                theoretical_results['processes']['sti']       = setup['theo_si_sti_disjoint']['wa_A-wb_B-wc_C_process'][1]
            else:
                raise ValueError('Model (option -m) needs to be either "disjoint" or "shared".')


            


            # -----------------------
            # plot
            # -----------------------
            ylim = [-0.9, 0.4]





            # -------------
            # Process sensitivities (ni constant)
            # -------------
            if method == 'DVM':

                nsets          = np.sort(list(map(int,numerical_results.keys())))
                nsets_parasets = np.sort(list(map(int,numerical_results[str(nsets[0])].keys())))

                si_numerical  = np.array([[ numerical_results[str(iset)][str(jset)]['processes'][ 'si']['out'] for jset in nsets_parasets ] for iset in nsets ])
                sti_numerical = np.array([[ numerical_results[str(iset)][str(jset)]['processes']['sti']['out'] for jset in nsets_parasets ] for iset in nsets ])
                
                jset = 128
                idx_jset = np.where(nsets_parasets==jset)[0][0]

                iplot += 1
                sub = fig.add_axes(position(nrow, ncol, iplot, hspace=hspace/2, vspace=vspace))
                sub.set_xscale('log')
                sub.set_ylim(ylim)
            
                for iiprocesses,iprocesses in enumerate(processes):
                    if iiprocesses < int(np.ceil(len(procopt)/2.)):
                        lstyle = '-'
                    else:
                        lstyle = '--'
                    mark1 = sub.plot(nsets, sti_numerical[:,idx_jset,iiprocesses]-theoretical_results['processes']['sti'][iiprocesses], alpha=1.0,
                                     linewidth=lwidth,
                                     linestyle=lstyle,
                                     color=colors_benchmark[iiprocesses%int(np.ceil(len(procopt)/2.))],
                                     marker='o',
                                     markeredgecolor=colors_benchmark[iiprocesses%int(np.ceil(len(procopt)/2.))],
                                     markerfacecolor='None',
                                     markersize=3.0, #msize,
                                     markeredgewidth=mwidth,
                                     label=iprocesses)    # Sti wmean
                sub.text(0.5,0.95, '$n_i = '+str(jset)+'$',verticalalignment='top',horizontalalignment='center',transform=sub.transAxes)
                plt.xlabel(str2tex("Number of Sobol' Reference Sets $K$",usetex=usetex))
                if iplot%2 == 1:
                    plt.ylabel(str2tex("$ST_{i}^{(appr)} - ST_{i}^{(theo)}$",usetex=usetex))
                else:
                    plt.ylabel(str2tex("",usetex=usetex))
                abc2plot(sub,dxabc,dyabc,iplot,bold=True,usetex=usetex,mathrm=True, large=True, parenthesis='none',verticalalignment='top')
                #plt.title(str2tex("Main Sobol' Index $S_i$",usetex=usetex))

            elif method == 'xSSA':

                nsets          = np.sort(list(map(int,numerical_results.keys())))

                si_numerical  = np.array([ numerical_results[str(iset)]['processes'][ 'si']['out'] for iset in nsets ])
                sti_numerical = np.array([ numerical_results[str(iset)]['processes']['sti']['out'] for iset in nsets ])
                
                iplot += 1
                sub = fig.add_axes(position(nrow, ncol, iplot, hspace=hspace/2, vspace=vspace))
                sub.set_xscale('log')
                sub.set_ylim(ylim)
            
                for iiprocesses,iprocesses in enumerate(processes):
                    if iiprocesses < int(np.ceil(len(procopt)/2.)):
                        lstyle = '-'
                    else:
                        lstyle = '--'
                    mark1 = sub.plot(nsets, sti_numerical[:,iiprocesses]-theoretical_results['processes']['sti'][iiprocesses], alpha=1.0,
                                     linewidth=lwidth,
                                     linestyle=lstyle,
                                     color=colors_benchmark[iiprocesses%int(np.ceil(len(procopt)/2.))],
                                     marker='o',
                                     markeredgecolor=colors_benchmark[iiprocesses%int(np.ceil(len(procopt)/2.))],
                                     markerfacecolor='None',
                                     markersize=3.0, #msize,
                                     markeredgewidth=mwidth,
                                     label=iprocesses)    # Sti wmean
                plt.xlabel(str2tex("Number of Sobol' Reference Sets $K$",usetex=usetex))
                if iplot%2 == 1:
                    plt.ylabel(str2tex("$ST_{i}^{(appr)} - ST_{i}^{(theo)}$",usetex=usetex))
                else:
                    plt.ylabel(str2tex("",usetex=usetex))
                abc2plot(sub,dxabc,dyabc,iplot,bold=True,usetex=usetex,mathrm=True, large=True, parenthesis='none',verticalalignment='top')
                #plt.title(str2tex("Main Sobol' Index $S_i$",usetex=usetex))

            else:
                raise ValueError('Method needs to be either "DVM" or "xSSA".')

            # legend
            ll = plt.legend(frameon=frameon, ncol=4,bbox_to_anchor=(llxbbox,llybbox), loc='lower right',
                           scatterpoints=1, numpoints=1,
                           labelspacing=llrspace, columnspacing=llcspace, handletextpad=llhtextpad, handlelength=llhlength)
            plt.setp(ll.get_texts(), fontsize='small')

            # # remove tick labels
            # ixticks = sub.get_xticks()
            # xnames  = [ r'' for i in ixticks]
            # xticknames = plt.setp(sub, xticklabels=xnames)

            # text to the right
            if iplot == 2:
                sub.text(1.10, 0.5, str2tex("Disjoint-parameter\n benchmark",usetex=usetex),rotation=90, fontsize='large', horizontalalignment='center', verticalalignment='center', transform=sub.transAxes)
            if iplot == 4:
                sub.text(1.10, 0.5, str2tex("Shared-parameter\n benchmark",usetex=usetex),rotation=90, fontsize='large', horizontalalignment='center', verticalalignment='center', transform=sub.transAxes)

            # text on top
            if iplot == 1:
                sub.text(0.5, 1.04, str2tex("Existing Discrete Values Method",usetex=usetex),rotation=0, fontsize='large', horizontalalignment='center', verticalalignment='bottom', transform=sub.transAxes)
            if iplot == 2:
                sub.text(0.5, 1.04, str2tex("Proposed xSSA Method",usetex=usetex),rotation=0, fontsize='large', horizontalalignment='center', verticalalignment='bottom', transform=sub.transAxes)


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


