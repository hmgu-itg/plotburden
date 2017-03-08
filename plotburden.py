#!/usr/bin/env python
import sys
import subprocess
import pandas as pd
from pandas import notnull, isnull
from numpy import log10, append, nan
import numpy as np
import re
import json, requests, asr
import urllib.request
import urllib.parse
import sys
import csv

#Personal libraries
import helper_functions
from helper_functions import *
from callbacks import *
from gene_plotter import *

gene=sys.argv[1]
input_monster=sys.argv[2]
output_monster=sys.argv[3]
sp_results=sys.argv[4]
vcf=sys.argv[5]
window=sys.argv[6]
output=sys.argv[7]




# global variables
server = "https://rest.ensembl.org";
helper_functions.server=server;



# Getting Gene coordinates
info("Querying Ensembl for gene coordinates...")
gc=get_coordinates(gene);
c=gc.chrom
start=gc.start
end=gc.end
info("\t\t⇰ Ensembl provided the coordinates "+str(c)+":"+str(start)+"-"+str(end)+" for gene "+gene)
start-=int(window)
end+=int(window)
gc.start=start
gc.end=end


## Getting variant consequences for all variants in the region
info("Querying Ensembl for SNP consequences and phenotype associations.")
resp=get_rsid_in_region(gc)
#resp.to_csv("snp.data", index=None, sep=",", quoting=csv.QUOTE_NONNUMERIC);
#resp=pd.read_table("~/snp.data", sep=",")
resp['pheno'].replace(to_replace=".*HGMD but.*", value="", inplace=True, regex=True)
info("\t\t⇰ Ensembl provided", len(resp),"known SNPs, ", len(resp[resp.pheno!=""]), "have associated phenotypes.")




## Get the single point results
sp = fetch_single_point(gc, sp_results)
info("Read", len(sp), "lines from single-point analysis.");
sp=pd.merge(sp, resp, on='ps', how='outer')
#rs_y is the ensembl rsid
sp.loc[isnull(sp.rs_y), 'rs_y']="novel"
sp.loc[isnull(sp.consequence), 'consequence']="novel"
sp.loc[isnull(sp.pheno), 'pheno']="none"
sp=sp[notnull(sp.chr)]


## Get the burden p from MONSTER output
results=pd.read_table(output_monster);
burdenp=results.p_MONSTER[0];
logburdenp=-log10(burdenp);
info("Burden p-value is", burdenp, "(", logburdenp, ").");



## Get the weights and variants in burden
variants=read_variants_from_gene_set(gc, input_monster)
info("Read ", variants.count()[0], "variants in burden")
rawdat=pd.merge(sp, variants, on='ps', how='outer')
if rawdat[rawdat.chr.isnull()].ps.size > 0 :
	warn(str(rawdat[rawdat.chr.isnull()].ps.size)+" variants were not found in the single point. They will be removed, but this is not a normal thing, please check your results.")
rawdat.dropna(subset=['chr'], inplace=True)



## Calculate LD
info("Calculating LD...")
task = subprocess.Popen(["/nfs/users/nfs_a/ag15/getld.sh", vcf, "chr"+str(c)+":"+str(start)+"-"+str(end), str(sp.size), str(end-start)], stdout=subprocess.PIPE);
ld=pd.read_table(task.stdout, sep='\s+');
os.remove("plink.log")
os.remove("plink.nosex")



## Defining plot-specific data
info("Defining plot-specific data...")
rawdat['radii']=150
rawdat.loc[rawdat.weight.notnull(), 'radii']=150+1000*rawdat.weight[rawdat.weight.notnull()]/max(rawdat.weight[rawdat.weight.notnull()])
rawdat['alpha']=0
rawdat.loc[rawdat.weight.notnull(), 'alpha']=0.8
rawdat['alpha_prevsig']=0
rawdat.loc[(rawdat.pheno!="none") & (rawdat.alpha>0), 'alpha_prevsig']=1
from bokeh.palettes import PuOr8 as palette
from bokeh.palettes import Viridis8 as palWeight
# Spectral9 Palette : ['#3288bd', '#66c2a5', '#abdda4', '#e6f598', '#ffffbf', '#fee08b', '#fdae61', '#f46d43', '#d53e4f']
palWeight.append("#939393")
rawdat['maf']=[af if af<0.5 else 1-af for af in rawdat.af]
rawdat['mafcolor']=[palette[i] for i in pd.cut(rawdat.maf, [-1, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.6]).cat.codes]
rawdat['color']="#1F77B4"
rawdat['weightcolor']=[palWeight[i] for i in pd.cut(rawdat.weight, 7).cat.codes]
rawdat['outcolor']="#3288bd"
rawdat["outalpha"]=0


# Plotting
info("Loading Bokeh...")
from bokeh.plotting import figure, output_file, show, save
from bokeh.layouts import layout, widgetbox, row, column
from bokeh.models.widgets import Button, RadioButtonGroup, Div
from bokeh.models import ColumnDataSource, CustomJS, HoverTool, LabelSet
output_file(output)
p1=figure(width=1500, x_range=[start, end], tools="box_zoom,tap,xwheel_zoom,reset", y_range=[-0.5, max(append(-log10(rawdat.p_score), logburdenp))+0.5])

ld_source=ColumnDataSource(data=dict(x1=ld.BP_A, x2=ld.BP_B, r2=ld.R2, dp=ld.DP))
source = ColumnDataSource(data=dict(ps=rawdat.ps, logsp=-log10(rawdat.p_score), radii=rawdat.radii, alpha=rawdat.alpha, color=rawdat.color, mafcolor=rawdat.mafcolor, weightcolor=rawdat.weightcolor, outcol=rawdat.outcolor, outalpha=rawdat.outalpha, alpha_prevsig=rawdat.alpha_prevsig))

## TODO
## this is the glyph for a third option in the "current signals"
## that should somehow highlight the points. The below draws an underlying asterisk
## not so good because circle is an actual circle, not a point. If this goes forward the points need to be made points.
## the callbacks also need to be written to update alpha_prevsig in rawdat.
#p1.asterisk(x='ps', y='logsp', size=20, alpha='alpha_prevsig', color="#F0027F", line_width=2, source=source)

p1.circle(x='ps', y='logsp', radius='radii', fill_alpha='alpha', fill_color='color', line_color='outcol', line_alpha='outalpha', line_width=6, source=source)
p1.xaxis.visible = False
segsource=ColumnDataSource(data=dict(y0=[logburdenp], y1=[logburdenp], x0=[min(rawdat.loc[rawdat.weight.notnull(), 'ps'])], x1=[max(rawdat.loc[rawdat.weight.notnull(), 'ps'])], alpha=[1], color=["firebrick"]))
p1.segment(y0='y0', y1='y1' , x0='x0', x1='x1', color='color', alpha='alpha', source=segsource, line_width=3)
x0=rawdat.ps[100]
y0=-log10(rawdat.p_score[100])
x1=rawdat.ps[400]
y1=-log10(rawdat.p_score[400])
bzier=ColumnDataSource(data=dict(x0=[], y0=[], x1=[], y1=[], cx0=[], cy0=[], cx1=[], cy1=[], col=[]))
p1.bezier(x0='x0', y0='y0', x1='x1', y1='y1', cx0='cx0', cy0='cy0', cx1='cx1', cy1='cy1', color='col', line_width=2, source=bzier)

showhide_sp=CustomJS(args=dict(source=source), code=showhide_sp_code)

changecolor=CustomJS(args=dict(source=source), code=changecolor_code)

hideburden=CustomJS(args=dict(source=segsource), code=hideburden_code)

ld_hover = CustomJS(args=dict(lds=ld_source, rawdat=source), code=ld_hover_code)


signalling=ColumnDataSource(data=dict(way=[0]))

ldbz_hover = CustomJS(args=dict(lds=ld_source, rawdat=source, bezier=bzier, signalling=signalling), code=ldbz_hover_code)


changehover = CustomJS(args=dict(signalling=signalling, rawdat=source, bezier=bzier), code=changehover_code)

resp=resp[notnull(resp.pheno)]
resp['alpha']=0
resp['y']=None
ray_source=ColumnDataSource(data=dict(ps=resp.ps, alpha=resp.alpha, y=resp.y, pheno=resp.pheno))
displayhits = CustomJS(args=dict(source=ray_source), code=displayhits_code)



p1.ray(x='ps', y='y', color="firebrick", length=0, angle=90, angle_units='deg', alpha='alpha', source=ray_source)
traits=LabelSet(x='ps', y=p1.y_range.end, y_offset=-0.5, text='pheno', level='glyph', text_alpha='alpha', angle=90, angle_units='deg', text_font_size='10pt', text_align='right', text_font_style='italic', source=ray_source)
p1.add_layout(traits)

p1.add_tools(HoverTool(callback=ldbz_hover, tooltips=None))
p_rbg=Div(text="""<strong>Single-point :</strong>""", width=100)
rbg = RadioButtonGroup(labels=["Hide non-burden", "Show all"], active=0, callback=showhide_sp, name="Hello")
p_chcolor=Div(text="""<strong>Colouring :</strong>""", width=100)
chcolor = RadioButtonGroup(labels=["None", "MAF", "Weight"], active=0, callback=changecolor)
p_burden=Div(text="""<strong>Show burden :</strong>""", width=100)
burden = RadioButtonGroup(labels=["Yes", "No"], active=0, callback=hideburden)
p_ld=Div(text="""<strong>LD behaviour :</strong>""", width=100)
control_ld = RadioButtonGroup(labels=["Highlight", "Fountain"], active=0, callback=changehover)
p_signals=Div(text="""<strong>Show Existing associations :</strong>""", width=100)
control_signals = RadioButtonGroup(labels=["No", "Rays"], active=0, callback=displayhits)

gc.extend(-int(window))
p2=draw_genes(gc, window, width=1500)
p2.x_range=p1.x_range

bbox=column(row([p_rbg, rbg]), row([p_chcolor, chcolor]), row([p_burden, burden]), row([p_ld, control_ld]), row([p_signals, control_signals]))
l=layout([[p1, bbox], [p2]])
save(l)
