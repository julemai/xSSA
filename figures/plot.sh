#!/bin/bash
#
# Produces plots for publication:
#
#      Mai, J., Craig, J. R., and Tolson, B. A. (2020).
#      Simultaneously Determining Global Sensitivities of Model Parameters and Model Structure. 
#      Hydrol. Earth Syst. Sci., accepted. 
#      https://doi.org/10.5194/hess-2020-215
#
#
# Copyright 2020 Juliane Mai - juliane.mai(at)uwaterloo.ca
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

# You should have received a copy of the GNU Lesser General Public License
# along with Juliane Mai's personal code library.  If not, see <http://www.gnu.org/licenses/>.
#

# Switches
dotex=1       	   #   	       LaTeX fonts in plots

dofig1=0      	   #   	       Model setups
dofig2=0           #   	       (A) Salmon River catchment, (B) monthly temperature/precipitation, (C) soil, and (D) landuse
dofig3=0           #   	       Simple/Realistic benchmark model: Error between theoretical and numerical values for different budgets of BARONI method
dofig4=1           #   	       Realistic benchmark model: Error between theoretical and numerical values for different budgets
dofig5_6=0         #   	       Results RAVEN sensitivity analysis
dofigS1=0      	   #           HMETS flowchart

verbose=0 # 0: pipe stdout and stderr to /dev/null
          # 1: pipe stdout to /dev/null
          # 2: print everything

# pathes
inputpath='../data'
outpath='pdfs'

# pdf margins
pdfmargins=3

# Treat switches
if [[ ${dotex} -eq 1 ]] ; then
    texit='-u'
else
    texit=''
fi

pipeit=''
if [[ ${verbose} -eq 0 ]] ; then pipeit=' > /dev/null 2>&1' ; fi
if [[ ${verbose} -eq 1 ]] ; then pipeit=' > /dev/null' ; fi

# Figures
if [[ ${dofig1} -eq 1 ]] ; then
    echo ''
    echo 'Figure 1 in progress...'
    # --------------------------------------------------------------------
    # Model setups
    # --------------------------------------------------------------------
    pdflatex ../figures/figure_1_model_setups/model_setups.tex  ${pipeit}
    pdfcrop model_setups.pdf
    mv model_setups-crop.pdf ${outpath}/figure_1.pdf
    rm model_setups.pdf
    rm model_setups.log
    rm model_setups.aux
fi

if [[ ${dofig2} -eq 1 ]] ; then
    echo ''
    echo 'figure 2 in progress...'
    python figure_2.py -t pdf -p 'figure_2' ${texit}  > /dev/null 2>&1
    mv figure_2.pdf ${outpath}/figure_2.pdf
    pdfcrop --margins ${pdfmargins} ${outpath}/figure_2.pdf > /dev/null 2>&1
    mv ${outpath}/figure_2-crop.pdf ${outpath}/figure_2.pdf
fi

if [[ ${dofig3} -eq 1 ]] ; then
    echo ''
    echo 'figure 3 in progress...'
    python figure_3.py -p figure_3 -t pdf ${texit}
    mv figure_3.pdf ${outpath}/figure_3.pdf
    pdfcrop --margins ${pdfmargins} ${outpath}/figure_3.pdf
    mv ${outpath}/figure_3-crop.pdf ${outpath}/figure_3.pdf
fi

if [[ ${dofig4} -eq 1 ]] ; then
    echo ''
    echo 'figure 4 in progress...'
    python figure_4.py -p figure_4 -t pdf -m 'shared' ${texit}
    mv figure_4.pdf ${outpath}/figure_4.pdf
    pdfcrop --margins ${pdfmargins} ${outpath}/figure_4.pdf > /dev/null 2>&1
    mv ${outpath}/figure_4-crop.pdf ${outpath}/figure_4.pdf
fi

if [[ ${dofig5_6} -eq 1 ]] ; then
    echo ''
    echo 'figure 5 and 6 in progress...'

    nsets=1000
    picklefile="../examples/data_out/08KC001/results_nsets${nsets}.nc"
    vars='Q'
    for ivar in ${vars} ; do
	python figure_5_6.py -n ${nsets} -t pdf -p 'figure_5_6' -v ${ivar} -i ${picklefile} -n ${nsets} -o "nc" ${texit}
	mv figure_5_6.pdf ${outpath}/figure_5_6.pdf
	pdfcrop --margins ${pdfmargins} ${outpath}/figure_5_6.pdf
	pdfsplit ${outpath}/figure_5_6-crop.pdf
	#rm ${outpath}/figure_5_6-crop1.pdf 
	mv ${outpath}/figure_5_6-crop1.pdf ${outpath}/figure_5.pdf
	mv ${outpath}/figure_5_6-crop2.pdf ${outpath}/figure_6.pdf
	#rm ${outpath}/figure_5_6-crop4.pdf
	rm ${outpath}/figure_5_6.pdf
	rm ${outpath}/figure_5_6-crop.pdf
    done
fi

if [[ ${dofigS1} -eq 1 ]] ; then
    echo ''
    echo 'Figure S1 in progress...'
    # --------------------------------------------------------------------
    # Model setups
    # --------------------------------------------------------------------
    pdflatex ../figures/figure_S1_HMETS_flowchart/HMETS_flowchart.tex  ${pipeit}
    pdfcrop HMETS_flowchart.pdf
    mv HMETS_flowchart-crop.pdf ${outpath}/figure_S1.pdf
    rm HMETS_flowchart.pdf
    rm HMETS_flowchart.log
    rm HMETS_flowchart.aux
fi

exit 0
