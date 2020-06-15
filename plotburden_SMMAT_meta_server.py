#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import subprocess
import pandas as pd
from pandas import notnull, isnull
import numpy as np
from numpy import log10, append, nan
import numpy as np
import re
import json, requests, asr
import urllib.request
import urllib.parse
import sys
import csv
import seaborn as sns
import random

#Personal libraries
import helper_functions
from helper_functions import *
from callbacks import *
from gene_plotter import *

pheno=sys.argv[1]
gene=sys.argv[2]
condition_string=sys.argv[3]
smmat_set_file=sys.argv[4]
smmat_out_file=sys.argv[5]
co_names=sys.argv[6]
sp_results=sys.argv[7]
vcf=sys.argv[8]
window=sys.argv[9]
output=sys.argv[10]
linkedFeatures=sys.argv[11]
chop=sys.argv[12]
#contdir=sys.argv[13]

import gene_plotter
gene_plotter.linkedFeatures=linkedFeatures
helper_functions.contdir=os.path.dirname(__file__)
# global variables
server = "https://rest.ensembl.org";
helper_functions.server=server;

# Getting Gene coordinates and extend the coordinates with the specified window:
info("Querying Ensembl for coordinates of "+gene+"...")
gc=get_coordinates(gene);
gc.extend(int(window))

#Extract coordinates:

c=gc.chrom
start = gc.start
end = gc.end
gene_start=gc.gstart
gene_end=gc.gend
ensid=gc.gene_id

# Report coordinates:
info("\t\t⇰ Ensembl provided the coordinates "+str(c)+":"+str(gene_start)+"-"+str(gene_end)+" for gene "+gene)
info("\t\t⇰ Plot boundaries: "+ str(c)+":"+str(start)+"-"+str(end))

## Getting variant consequences for all variants in the region
info("Querying Ensembl for SNP consequences and phenotype associations.")
resp=get_rsid_in_region(gc)
#resp.to_csv(gene+".snp.data", index=None, sep=",", quoting=csv.QUOTE_NONNUMERIC);
#resp=pd.read_table("snp.data", sep=",")
resp['pheno'].replace(to_replace="Annotated by HGMD*", value="", inplace=True, regex=True)
resp['pheno'].replace(to_replace="ClinVar.*not specified", value="", inplace=True, regex=True)
resp.loc[isnull(resp.pheno), 'pheno']="none"
resp.pheno=resp.pheno.str.strip()
resp.loc[resp.pheno=="", 'pheno']="none"
resp['ensembl_rs']=resp['rs']
resp.drop('rs', axis=1, inplace=True)
info("\t\t⇰ Ensembl provided", len(resp),"known SNPs, ", len(resp[resp.pheno!="none"]), "have associated phenotypes.")




## Get the single point results. Returns one merged (outer) dataframe with all columns suffixed by the cohort name
## We are just going to use this for annotation purposes
sp = fetch_single_point_meta(gc, sp_results, co_names)
info("Read", len(sp), "lines from single-point analysis.");
sp=pd.merge(sp, resp, on='ps', how='outer')
sp.loc[isnull(sp.ensembl_rs), 'ensembl_rs']="novel"
sp.loc[isnull(sp.consequence), 'consequence']="novel"
sp=sp[notnull(sp.chr)]
sp['ensembl_consequence']=sp['consequence']
sp['chr']=sp['chr'].astype(int)
info("Getting consequences for novel variants...")
sp=get_csq_novel_variants(sp, 'chr', 'ps', 'allele0', 'allele1')



## Get the burden p-values
## Returns a dictionary indexed by the cohort name or "meta"
info("Reading burden P-values...")
#results=read_burden_ps(co_names,smmat_out_file, ensid, pheno, condition_string)
import pickle
with open('results.bin', 'rb') as config_dictionary_file:
	results = pickle.load(config_dictionary_file)

## read all variants in the gene sets including those not in some cohorts
info("Reading variants from gene set...")
variants=read_variants_from_gene_set_SMMAT(ensid, condition_string, smmat_set_file)
info("Read ", variants.count()[0], "variants in burden across all cohorts")
## Now for the plot data
##
cohdat={}
lddat={}
i=0
import pickle
for n in co_names.split(","):
	info("Preparing plot data for cohort "+n)
	rawdat, ld=produce_single_cohort_df(gc, sp_results.split(',')[i], resp, vcf.split(",")[i], smmat_out_file.split(',')[i], smmat_set_file, pheno, condition_string, n, sp, variants)
	cohdat[i]=rawdat
	lddat[i]=ld
#	cohdat[i]=ColumnDataSource(data=dict(ps=rawdat.ps, logsp=-log10(rawdat.p_score), radii=rawdat.radii, alpha=rawdat.alpha, color=rawdat.color, mafcolor=rawdat.mafcolor, weightcolor=rawdat.weightcolor, outcol=rawdat.outcolor, outalpha=rawdat.outalpha, alpha_prevsig=rawdat.alpha_prevsig, snpid=rawdat.rs, rs=rawdat.ensembl_rs, maf=rawdat.maf, csq=rawdat.ensembl_consequence))
	i=i+1

cohdat[i], lddat[i] = produce_meta_df(gc, sp, variants, vcf, co_names)
with open('sp.bin', 'wb') as config_dictionary_file:
	pickle.dump(sp, config_dictionary_file)


# with open('resp.bin', 'wb') as config_dictionary_file:
# 	pickle.dump(resp, config_dictionary_file)
#
# with open('ld.bin', 'wb') as config_dictionary_file:
# 	pickle.dump(lddat, config_dictionary_file)
#
# with open('results.bin', 'wb') as config_dictionary_file:
# 	pickle.dump(results, config_dictionary_file)
# with open('gc.bin', 'wb') as config_dictionary_file:
# 	pickle.dump(gc, config_dictionary_file)
#
# # =================================================================================================
# # =================================================================================================
# # =================================================================================================
#
# import pickle
#
# with open('cohdat.bin', 'rb') as config_dictionary_file:
# 	cohdat = pickle.load(config_dictionary_file)
#
# with open('gc.bin', 'rb') as config_dictionary_file:
# 	gc = pickle.load(config_dictionary_file)
# c=gc.chrom
# start = gc.start
# end = gc.end
# gene_start=gc.gstart
# gene_end=gc.gend
# ensid=gc.gene_id
#
#
# with open('ld.bin', 'rb') as config_dictionary_file:
# 	lddat = pickle.load(config_dictionary_file)
#
# with open('results.bin', 'rb') as config_dictionary_file:
# 	results = pickle.load(config_dictionary_file)
#
# with open('resp.bin', 'rb') as config_dictionary_file:
# 	resp = pickle.load(config_dictionary_file)

## Preparing plot variables
## maxlogp is the min of all ps, that is used to set plot ylim
## cohdat is overwritten to become an array of dataframes containing the SP info for the respective cohorts

i=0
rawdats=[]
maxlogp=0
for n in co_names.split(","):
	rawdats.append(cohdat[i])
	rawdat=cohdat[i]
	if -log10(rawdat.p_score.min())>maxlogp:
		maxlogp=-log10(rawdat.p_score.min())
	cohdat[i]=dict(ps=rawdat.ps, logsp=-log10(rawdat.p_score), radii=rawdat.radii, alpha=rawdat.alpha, color=rawdat.color, mafcolor=rawdat.mafcolor, weightcolor=rawdat.weightcolor, outcol=rawdat.outcolor, outalpha=rawdat.outalpha, alpha_prevsig=rawdat.alpha_prevsig, snpid=rawdat.rs, rs=rawdat.ensembl_rs, maf=rawdat.maf, csq=rawdat.ensembl_consequence)
	i=i+1

## meta-analysis (beware, ALL columns are as is)
rawdats.append(cohdat[i])
rawdat=cohdat[i]
print(rawdat.columns)
if -log10(rawdat["P-valuemeta"].min())>maxlogp:
	maxlogp=-log10(rawdat["P-valuemeta"].min())
cohdat[i]=dict(ps=rawdat.ps, logsp=-log10(rawdat["P-valuemeta"]), radii=rawdat.radii, alpha=rawdat.alpha, color=rawdat.color, mafcolor=rawdat.mafcolor, weightcolor=rawdat.weightcolor, outcol=rawdat.outcolor, outalpha=rawdat.outalpha, alpha_prevsig=rawdat.alpha_prevsig, snpid=rawdat.chr.astype(str)+":"+rawdat.ps.astype(str), rs=rawdat.ensembl_rs, maf=rawdat.maf, csq=rawdat.ensembl_consequence)


## Creating the df containing ALL points coloured by cohort
co_split=co_names.split(",")
palette=sns.color_palette("Set2").as_hex()
random.shuffle(palette)

cter=0
bigdf=None
for data in cohdat:
	df=pd.DataFrame(cohdat[data])
	if(data != (len(cohdat)-1)):
		# set color to cohort-based color
		df['cocolor']=palette[cter]
		#set additional column to cohort name
		df['cohort']=co_split[cter]
	#rbind it
	if(bigdf is None):
		bigdf=df
	else:
		bigdf=bigdf.append(df, ignore_index=True)
	cter=cter+1


## append the p ranges to sp for the m/a segments
for co in co_split:
    if co=="meta":
        col='P-valuemeta'
    else:
        col='p_score'+co
    sp['logp'+co]=-log10(sp[col])


sp['minp']=sp[["logp"+s for s in co_split]].min(axis=1)
sp['maxp']=sp[["logp"+s for s in co_split]].max(axis=1)
sp['segcol']="gray"
sp=pd.merge(sp, variants, on="ps", how='outer')
sp.dropna(subset=['chr'], inplace=True)
#sp.loc[sp.weight.isnull(), 'segcol']="#FFFFFF00"


# with open('cohdat.bin', 'wb') as config_dictionary_file:
# 	pickle.dump(cohdat, config_dictionary_file)
# with open('rawdat.bin', 'wb') as config_dictionary_file:
# 	pickle.dump(rawdats, config_dictionary_file)

rawdat=rawdats[0]
if -log10(min(results.values()))>maxlogp:
	maxlogp=-log10(min(results.values()))

import pickle
plotdat=dict(rawdats=rawdats, rawdat=rawdat, maxlogp=maxlogp, gene=gene, gc=gc, resp=resp, lddat=lddat, sp=sp, cohdat=cohdat, co_split=co_split, results=results)
with open('plotdat.bin', 'wb') as config_dictionary_file:
	pickle.dump(plotdat, config_dictionary_file)


import pickle


info("Loading Bokeh...")
from bokeh.plotting import figure, output_file, show, save, curdoc
from bokeh.layouts import layout, widgetbox, row, column
from bokeh.models.widgets import Button, RadioButtonGroup, Div
from bokeh.models import ColumnDataSource, CustomJS, HoverTool, LabelSet, OpenURL, TapTool, Axis, SaveTool

### Initialising the burden p-value and corresponding segment
gc2=get_coordinates(gene);
logburdenp=-1*log10(results[co_split[0]])
if (max(rawdats[0].loc[rawdats[0].weight.notnull(), 'ps'])-min(rawdats[0].loc[rawdats[0].weight.notnull(), 'ps']) < 500):
	eseg=gc2.end
	sseg=gc2.start
else:
	eseg=max(rawdats[0].loc[rawdats[0].weight.notnull(), 'ps'])
	sseg=min(rawdats[0].loc[rawdats[0].weight.notnull(), 'ps'])
segsource=ColumnDataSource(data=dict(y0=[logburdenp], y1=[logburdenp], x0=[sseg], x1=[eseg], alpha=[1], color=["firebrick"]))

### Initialising the figure
p1=figure(width=1500, x_range=[start, end], tools="box_zoom,tap,xwheel_zoom,reset,save", y_range=[-0.5, maxlogp+0.5])

## Previous signals use the resp dataframe. Putting this code here ensures they are behind everything
resp=resp[notnull(resp.pheno)]
resp=resp[resp.pheno!="none"]
resp['alpha']=0
resp['y']=None
ray_source=ColumnDataSource(data=dict(ps=resp.ps, alpha=resp.alpha, y=resp.y, pheno=resp.pheno))
displayhits = CustomJS(args=dict(source=ray_source), code=displayhits_code)
p1.segment(x0='ps', x1='ps', y0=p1.y_range.start, color="firebrick", y1=p1.y_range.end, alpha='alpha', source=ray_source)
traits=LabelSet(x='ps', y=p1.y_range.end, y_offset=-0.5, text='pheno', level='glyph', text_alpha='alpha', angle=90, angle_units='deg', text_font_size='10pt', text_align='right', text_font_style='italic', source=ray_source)
p1.add_layout(traits)


### Initialising LD data
ld=lddat[0]
ld_source=ColumnDataSource(data=dict(x1=ld.BP_A, x2=ld.BP_B, r2=ld.R2, dp=ld.DP))

### the meta-analysis segments, initialised void
metasegsource=ColumnDataSource(data=dict(ps=[], minp=[], maxp=[], segcol=[]))
p1.segment(x0='ps', x1='ps', y0='minp', y1='maxp', line_color='segcol', source=metasegsource)


## Initialising source for points
source = ColumnDataSource(cohdat[0])
mainplot_points=p1.circle(x='ps', y='logsp', radius='radii', fill_alpha='alpha', fill_color='color', line_color='outcol', line_alpha='outalpha', line_width=6, radius_units='screen', source=source)
p1.xaxis.visible = False

## now that we have the figure, add the segment
p1.segment(y0='y0', y1='y1' , x0='x0', x1='x1', color='color', alpha='alpha', source=segsource, line_width=3)


x0=rawdat.ps[100]
y0=-log10(rawdat.p_score[100])
x1=rawdat.ps[400]
y1=-log10(rawdat.p_score[400])
bzier=ColumnDataSource(data=dict(x0=[], y0=[], x1=[], y1=[], cx0=[], cy0=[], cx1=[], cy1=[], col=[]))
p1.bezier(x0='x0', y0='y0', x1='x1', y1='y1', cx0='cx0', cy0='cy0', cx1='cx1', cy1='cy1', color='col', line_width=2, source=bzier)


## Destined to die: JS callbacks
showhide_sp=CustomJS(args=dict(source=source), code=showhide_sp_code)
changecolor=CustomJS(args=dict(source=source), code=changecolor_code)
hideburden=CustomJS(args=dict(source=segsource), code=hideburden_code)
ld_hover = CustomJS(args=dict(lds=ld_source, rawdat=source), code=ld_hover_code)
signalling=ColumnDataSource(data=dict(way=[0]))
ldbz_hover = CustomJS(args=dict(lds=ld_source, rawdat=source, bezier=bzier, signalling=signalling), code=ldbz_hover_code)
changehover = CustomJS(args=dict(signalling=signalling, rawdat=source, bezier=bzier), code=changehover_code)
testhover=CustomJS(args=dict(source=source), code=hover_test_code)

p1.add_tools(HoverTool(callback=ldbz_hover, tooltips=[("SNPid", "@snpid"), ("RSid", "@rs"), ("MAF", "@maf"), ("consequence", "@csq")], renderers=[mainplot_points]))

taptool = p1.select(type=TapTool)
taptool.callback = OpenURL(url="http://www.ensembl.org/Homo_sapiens/Variation/Explore?db=core;v=@rs;vdb=variation")

co_split.append("meta")
p_source=Div(text="""<strong>Show data for :</strong>""", width=100)
control_source = RadioButtonGroup(labels=co_split, active=0)
#control_source.js_on_change('active1',changecohort)
p_rbg=Div(text="""<strong>Single-point :</strong>""", width=100)
rbg = RadioButtonGroup(labels=["Hide non-burden", "Show all"], active=0)
p_chcolor=Div(text="""<strong>Colouring :</strong>""", width=100)
chcolor = RadioButtonGroup(labels=["None", "MAF", "Weight"], active=0)
p_burden=Div(text="""<strong>Show burden :</strong>""", width=100)
burden = RadioButtonGroup(labels=["Yes", "No"], active=0)
p_ld=Div(text="""<strong>LD behaviour :</strong>""", width=100)
control_ld = RadioButtonGroup(labels=["Highlight", "Fountain"], active=0)
control_ld.js_on_change('active', changehover)
p_signals=Div(text="""<strong>Show Existing associations :</strong>""", width=100)
control_signals = RadioButtonGroup(labels=["No", "Rays"], active=0)
p_meta=Div(text="""<strong>Show all cohorts in meta-analysis :</strong>""", width=100)
control_meta = RadioButtonGroup(labels=["No", "Yes"], active=0, disabled=True)
p1.yaxis[0].major_label_text_font_size = "13pt"
p1.yaxis.axis_label = '-log₁₀(p)'
p1.yaxis.axis_label_text_font_size = "15pt"

#gc.extend(-int(window)) # <- Not needed anymore. gc object contains the gene start and gene end position: gc.gstart and gc.gend



## Callback for the source change
def callback_changesource(arg):
	print(control_source.active)
	print(arg)
	if(arg == (len(co_split)-1)):
		## we are in case of meta-analysis
		control_meta.disabled=False
	else:
		control_meta.disabled=True
	print(co_split[arg])
	if (max(rawdats[arg].loc[rawdats[arg].weight.notnull(), 'ps'])-min(rawdats[arg].loc[rawdats[arg].weight.notnull(), 'ps']) < 500):
		eseg=gc2.end
		sseg=gc2.start
	else:
		eseg=max(rawdats[arg].loc[rawdats[arg].weight.notnull(), 'ps'])
		sseg=min(rawdats[arg].loc[rawdats[arg].weight.notnull(), 'ps'])
	segsource.data=dict(y0=[-1*log10(results[co_split[arg]])], y1=[-1*log10(results[co_split[arg]])], x0=[sseg], x1=[eseg], alpha=[1], color=["firebrick"])
	## This changes LD
	ld=lddat[control_source.active]
	ld_source.data=dict(x1=ld.BP_A, x2=ld.BP_B, r2=ld.R2, dp=ld.DP)
	## BEWARE FOR COMPREHENSION
	## The source is actually called by showsp, which both maintains the current state of the single-point button AND switches the source.
	callback_showsp(rbg.active)

## Callback for the single-point display
def callback_showsp(arg):
	currawdat=rawdats[control_source.active]
	if arg == 0:
		# do not show sp
		print("Disabled single point")
		currawdat.loc[currawdat.weight.isnull(), 'alpha']=0
	else:
		print("Enabled single point")
		currawdat.loc[currawdat.weight.isnull(), 'alpha']=1
	#source.data=cohdat[i]
	if control_source.active==(len(co_split)-1):
		## we are in m/a
		source.data=dict(ps=currawdat.ps, logsp=-log10(currawdat["P-valuemeta"]), radii=currawdat.radii, alpha=currawdat.alpha, color=currawdat.color, mafcolor=currawdat.mafcolor, weightcolor=currawdat.weightcolor, outcol=currawdat.outcolor, outalpha=currawdat.outalpha, alpha_prevsig=currawdat.alpha_prevsig, snpid=currawdat.chr.astype(str)+":"+currawdat.ps.astype(str), rs=currawdat.ensembl_rs, maf=currawdat.maf, csq=currawdat.ensembl_consequence)
	else:
		source.data=dict(ps=currawdat.ps, logsp=-log10(currawdat.p_score), radii=currawdat.radii, alpha=currawdat.alpha, color=currawdat.color, mafcolor=currawdat.mafcolor, weightcolor=currawdat.weightcolor, outcol=currawdat.outcolor, outalpha=currawdat.outalpha, alpha_prevsig=currawdat.alpha_prevsig, snpid=currawdat.rs, rs=currawdat.ensembl_rs, maf=currawdat.maf, csq=currawdat.ensembl_consequence)
	callback_chcolor(chcolor.active)
	callback_burden(burden.active)

## Callback for the association rays
def callback_toggleassoc(arg):
	print("Ray control toggled to "+str(arg))
	if(control_signals.active==0):
		ray_source.data["alpha"]=np.repeat(0, len(ray_source.data["ps"]))
	else:
		ray_source.data["alpha"]=np.repeat(1, len(ray_source.data["ps"]))

## callback for the color change
def callback_chcolor(arg):
	if arg == 0 :
		if control_meta.disabled | control_meta.active==0:
			source.data['color']=np.repeat("#3288bd", len(source.data['ps']))
		else:
			# in this case the source is set to bigdf which has the color
			source.data['color']=source.data['cocolor']
	if arg == 1 :
		## color by Weight
		source.data['color']=source.data['mafcolor']
	if arg == 2 :
		## color by MAF
		source.data['color']=source.data['weightcolor']

## callback to toggle burden
def callback_burden(arg):
	revert=(1-arg)*(1-arg)
	segsource.data['alpha']=[revert]

## test callback ld hovered
def callback_hover(cb_data):
	print(cb_data.Index)

def callback_meta(arg):
	## here we toggle between bigdf or no.
	if(arg==1):
		metasegsource.data=sp
		if rbg.active == 0 :
			bigdf.loc[bigdf.weightcolor=="#939393", 'alpha']=0
		else:
			bigdf.loc[bigdf.weightcolor=="#939393", 'alpha']=1
		source.data=bigdf
		callback_chcolor(chcolor.active)
	else:
		callback_changesource(control_source.active)

## adding callbacks to elements
control_meta.on_click(callback_meta)
control_source.on_click(callback_changesource)
control_signals.on_click(callback_toggleassoc)
rbg.on_click(callback_showsp)
chcolor.on_click(callback_chcolor)
burden.on_click(callback_burden)

#ld_hovertool=HoverTool(callback=callback_hover, renderers=[mainplot_points])
#ld_hovertool=HoverTool(callback=ld_hover, renderers=[mainplot_points])
#p1.add_tools(ld_hovertool)


# p2=draw_genes(gc, window, width=1500, chop=chop)
# p2.x_range=p1.x_range
# p2.xaxis[0].formatter.use_scientific = False
# p2.xaxis[0].major_label_text_font_size = "13pt"
# p2.xaxis.axis_label = 'position on chromosome '+str(rawdat.chr[1].astype(np.int64))
# p2.xaxis.axis_label_text_font_size = "15pt"


bbox=column(row([p_source, control_source]),row([p_rbg, rbg]), row([p_chcolor, chcolor]), row([p_burden, burden]), row([p_ld, control_ld]), row([p_signals, control_signals]), row([p_meta, control_meta]))
#l=layout([[p1, bbox], [p2]])
l=layout([[p1, bbox]])
#p1.output_backend = "svg" #NOT FUNCTIONAL
#p2.output_backend = "svg"
#save(l)
curdoc().add_root(l)
#rawdat.to_csv(output+".csv", index=False)
#ld.to_csv(output+".ld.csv",index=False)
