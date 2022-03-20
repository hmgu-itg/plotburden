import os
import sys
import json

import numpy as np
import pandas as pd
from bokeh.models import (FixedTicker, FuncTickFormatter, TableColumn, DataTable,
						  ColumnDataSource, CustomJS, HoverTool, LabelSet, OpenURL,
						  TapTool)
from bokeh.models.widgets import RadioButtonGroup, Div
from bokeh.layouts import layout, row, column
from bokeh.plotting import figure, curdoc

from plotburden import helper_functions, gene_plotter
import callbacks as cb


helper_functions.contdir=os.path.dirname(__file__)
# global variables
server = "https://rest.ensembl.org"
helper_functions.server=server


plotdatfile=sys.argv[1]

#plotdat=dict(rawdats=rawdats, rawdat=rawdat, maxlogp=maxlogp, gene=gene, gc=gc, resp=resp, lddat=lddat, sp=sp, cohdat=cohdat, co_split=co_split, results=results)

# with open(plotdatfile, 'rb') as config_dictionary_file:
# 	plotdat = pickle.load(config_dictionary_file)

with open(plotdatfile, 'r') as config_dictionary_file:
	plotdat = json.load(config_dictionary_file)


# rawdats=plotdat['rawdats']
rawdats = [pd.DataFrame(v) for v in plotdat['rawdats']]
# rawdat=plotdat['rawdat']
rawdat = pd.DataFrame(plotdat['rawdat'])

maxlogp=plotdat['maxlogp']
gene=plotdat['gene']
# gc=plotdat['gc']
gc: dict = plotdat['gc']
gc = helper_functions.GeneCoordinates(gc['chrom'], gc['start'], gc['end'], gc['gene_id'], gc['name'])

# resp=plotdat['resp']
resp = pd.DataFrame(plotdat['resp'])
# lddat=plotdat['lddat']
lddat = {int(k): pd.DataFrame(v) for k, v in plotdat['lddat'].items()}

# sp=plotdat['sp']
sp = pd.DataFrame(plotdat['sp'])

# cohdat=plotdat['cohdat']
cohdat = {int(k): pd.DataFrame(v) for k, v in plotdat['cohdat'].items()}

co_split=plotdat['co_split']
results=plotdat['results']
# bigdf=plotdat['bigdf']
bigdf = pd.DataFrame(plotdat['bigdf'])

window=plotdat['window']
chop=plotdat['chop']
pheno=plotdat['pheno']
condition_string=plotdat['condition_string']
linkedFeatures=plotdat['linkedFeatures']
meta=plotdat['meta']

cohort_color=bigdf[["cocolor", "cohort"]]
cohort_color.drop_duplicates(inplace=True)
cohort_color=cohort_color.set_index('cohort').to_dict()
cohort_color=cohort_color['cocolor']



c=gc.chrom
start = gc.start
end = gc.end
gene_start=gc.gstart
gene_end=gc.gend
ensid=gc.gene_id


helper_functions.info("Loading Bokeh...")

### Initialising the burden p-value and corresponding segment

# from pathlib import Path
# g = Path('gene.bin')
# if g.exists():
# 	with open(g, 'rb') as f:
# 		gc2 = pickle.load(f)
# else:
# 	gc2 = helper_functions.get_coordinates(gene)
# 	with open(g, 'wb') as f:
# 		pickle.dump(gc2, f)
gc2 = helper_functions.get_coordinates(gene)

burden_p = results[co_split[0]]
logburdenp = -1 * np.log10(burden_p)
if (max(rawdats[0].loc[rawdats[0].weight.notnull(), 'ps'])-min(rawdats[0].loc[rawdats[0].weight.notnull(), 'ps']) < 500):
	eseg=gc2.end
	sseg=gc2.start
else:
	eseg=max(rawdats[0].loc[rawdats[0].weight.notnull(), 'ps'])
	sseg=min(rawdats[0].loc[rawdats[0].weight.notnull(), 'ps'])
segsource=ColumnDataSource(data=dict(y0=[logburdenp], y1=[logburdenp], x0=[sseg], x1=[eseg], alpha=[1], color=["firebrick"]))

### Initialising the figure
p1=figure(x_range=[start, end], tools="box_zoom,lasso_select,tap,xwheel_zoom,reset,save", y_range=[-0.5, maxlogp+0.5], width=1500)

## Previous signals use the resp dataframe. Putting this code here ensures they are behind everything
resp=resp[pd.notnull(resp.pheno)]
resp=resp[resp.pheno!="none"]
resp['alpha']=0
resp['y']=None
ray_source=ColumnDataSource(data=dict(ps=resp.ps, alpha=resp.alpha, y=resp.y, pheno=resp.pheno))
displayhits = CustomJS(args=dict(source=ray_source), code=cb.displayhits_code)
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
sourceadd=ColumnDataSource(data=dict(ps=[], logsp=[], radii=[], alpha=[], color=[], mafcolor=[], weightcolor=[], outcol=[], outalpha=[], alpha_prevsig=[], snpid=[], rs=[], maf=[], csq=[]))
mainplot_points=p1.circle(x='ps', y='logsp', radius='radii', fill_alpha='alpha', fill_color='color', line_color='outcol', line_alpha='outalpha', line_width=6, radius_units='screen', source=source)
additional_points=p1.circle(x='ps', y='logsp', radius='radii', fill_alpha='alpha', fill_color='color', line_color='outcol', line_alpha='outalpha', line_width=6, radius_units='screen', source=sourceadd)

p1.xaxis.visible = False

## now that we have the figure, add the segment
p1.segment(y0='y0', y1='y1' , x0='x0', x1='x1', color='color', alpha='alpha', source=segsource, line_width=3)


x0=rawdat.ps[100]
y0=rawdat.logp[100]
x1=rawdat.ps[400]
y1=rawdat.logp[400]
bzier=ColumnDataSource(data=dict(x0=[], y0=[], x1=[], y1=[], cx0=[], cy0=[], cx1=[], cy1=[], col=[]))
p1.bezier(x0='x0', y0='y0', x1='x1', y1='y1', cx0='cx0', cy0='cy0', cx1='cx1', cy1='cy1', color='col', line_width=2, source=bzier)


## Destined to die: JS callbacks
showhide_sp=CustomJS(args=dict(source=source), code=cb.showhide_sp_code)
changecolor=CustomJS(args=dict(source=source), code=cb.changecolor_code)
hideburden=CustomJS(args=dict(source=segsource), code=cb.hideburden_code)
ld_hover = CustomJS(args=dict(lds=ld_source, rawdat=source), code=cb.ld_hover_code)
signalling=ColumnDataSource(data=dict(way=[0]))
ldbz_hover = CustomJS(args=dict(lds=ld_source, rawdat=source, bezier=bzier, signalling=signalling), code=cb.ldbz_hover_code)
changehover = CustomJS(args=dict(signalling=signalling, rawdat=source, bezier=bzier), code=cb.changehover_code)
testhover=CustomJS(args=dict(source=source), code=cb.hover_test_code)

p1.add_tools(HoverTool(callback=ldbz_hover, tooltips=[("SNPid", "@snpid"), ("RSid", "@rs"), ("p-value", "@p_value"), ("MAF", "@maf"), ("consequence", "@csq")], renderers=[mainplot_points]))

taptool = p1.select(type=TapTool)
taptool.callback = OpenURL(url="http://www.ensembl.org/Homo_sapiens/Variation/Explore?db=core;v=@rs;vdb=variation")

p_title=Div(text="""<h1>"""+pheno+""" burden in <i>"""+gene+"""</i> ("""+condition_string+""")</h1>""", width=1500, style={'text-align':'center', 'margin':'auto', "width":"100%"})

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
p_click=Div(text="""<strong>On click :</strong>""", width=100)
control_click = RadioButtonGroup(labels=["Ensembl", "Show meta-analysis"], active=0)
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
	burden_p = results[co_split[arg]]
	logburdenp = -1 * np.log10(burden_p)
	segsource.data=dict(y0=[logburdenp], y1=[logburdenp], x0=[sseg], x1=[eseg], alpha=[1], color=["firebrick"])
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
	if meta is True and control_source.active==(len(co_split)-1):
		## we are in m/a
		source.data=dict(ps=currawdat.ps, p_value=currawdat.p_score, logsp=currawdat.logp_meta, radii=currawdat.radii, alpha=currawdat.alpha, color=currawdat.color, mafcolor=currawdat.mafcolor, weightcolor=currawdat.weightcolor, outcol=currawdat.outcolor, outalpha=currawdat.outalpha, alpha_prevsig=currawdat.alpha_prevsig, snpid=currawdat.chr.astype(str)+":"+currawdat.ps.astype(str), rs=currawdat.ensembl_rs, maf=currawdat.maf, csq=currawdat.ensembl_consequence)
	else:
		source.data=dict(ps=currawdat.ps, p_value=currawdat.p_score, logsp=currawdat.logp, radii=currawdat.radii, alpha=currawdat.alpha, color=currawdat.color, mafcolor=currawdat.mafcolor, weightcolor=currawdat.weightcolor, outcol=currawdat.outcolor, outalpha=currawdat.outalpha, alpha_prevsig=currawdat.alpha_prevsig, snpid=currawdat.rs, rs=currawdat.ensembl_rs, maf=currawdat.maf, csq=currawdat.ensembl_consequence)
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
		if (control_meta.disabled) or (control_meta.active==0):
			source.data['color']=np.repeat("#3288bd", len(source.data['ps']))
		else:
			## UNREACHABLE unless control_meta is uncommented in the code
			# in this case the source is set to bigdf which has the color
			#print(source.data)
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
		if rbg.active == 0 :
			usedf=bigdf[bigdf.weightcolor!="#939393"]
			usesp=sp[sp.weight.notnull()]
		else:
			bigdf.loc[bigdf.weightcolor=="#939393", 'alpha']=1
			usedf=bigdf
			usesp=sp
		source.data=usedf
		metasegsource.data=usesp
		callback_chcolor(chcolor.active)
	else:
		callback_changesource(control_source.active)
		metasegsource.data=dict(ps=[], minp=[], maxp=[], segcol=[])

p3 = figure(width=430, title="Forest Plot (enable show meta-analysis and select points)", y_range=[-1, len(co_split)])
def generate_forestplot_df(selected):
	print("selected points:")
	print(selected)
	directions=list(selected.Direction_meta)
	ptdf=[]
	k=0
	p3.title.text="Forest plot for "+str(int(selected.chr))+":"+str(int(selected.ps))
	for n in co_split:
		if(n=="meta"):
			x=selected["Effect_"+n].squeeze()
			y=k
			start=x-selected['StdErr_'+n].squeeze()
			end=x+selected['StdErr_'+n].squeeze()
		else:
			x=abs(selected["beta_"+n].squeeze())
			if(directions[k]=="-"):
				x=-1*x
			y=k
			start=x-selected['se_'+n].squeeze()
			end=x+selected['se_'+n].squeeze()
		color=cohort_color[n]
		ptdf.append([x,y,start,end, color])
		k=k+1
	ptdf=pd.DataFrame(ptdf, columns=["x", "y", "start", "end", "color"])
	return(ptdf)
forestsource=ColumnDataSource(dict(x=[], y=[], start=[], end=[], color=[]))
p3.segment(x0=0, x1=0, y0=-20, y1=20, color="#F4A582", line_width=2, line_dash="dashed")
p3.square(x='x', y='y', color='color', source=forestsource)
p3.segment(x0='start',x1='end', y0='y', y1='y', color='color', source=forestsource)
p3.yaxis.ticker = FixedTicker(ticks=list(range(0, len(co_split))))
namedict = { i : co_split[i] for i in range(0, len(co_split) ) }
p3.yaxis.formatter = FuncTickFormatter(args=dict(namedict=namedict), code="""return(namedict[tick])""")

#callback for when user selects meta points (draws segments, colors, etc)
def callback_pointclicked(attr, old, new):
	df=pd.DataFrame(source.data)
	print("Length of df is "+str(len(df.index)))
	print("New is ")
	print(new)

	selected=df.iloc[new]
	selected=selected[selected.alpha!=0]
	todisplay=bigdf[bigdf.ps.isin(selected.ps)]
	print(todisplay)
	toappend=dict(ps=todisplay.ps, logsp=todisplay.logsp, radii=todisplay.radii, alpha=todisplay.alpha, color=todisplay.cocolor, mafcolor=todisplay.mafcolor, weightcolor=todisplay.weightcolor, outcol=todisplay.weightcolor, outalpha=todisplay.outalpha, alpha_prevsig=todisplay.alpha_prevsig, snpid=todisplay.snpid, rs=todisplay.rs, maf=todisplay.maf, csq=todisplay.csq)
	toappend=pd.DataFrame(toappend)
	toappend.alpha=1

	sourceadd.data=toappend
	df=pd.DataFrame(sourceadd.data)
	print(len(df.index))

	#todisplay=sp[sp.ps.isin(selected.ps)]
	#mss=dict(ps=todisplay.ps, minp=todisplay.minp, maxp=todisplay.minp, segcol=todisplay.segcol)
	df=pd.DataFrame(metasegsource.data)
	print(len(df.index))
	metasegsource.data=sp[sp.ps.isin(selected.ps)]
	df=pd.DataFrame(metasegsource.data)
	print(len(df.index))

	forest=sp[sp.ps.isin(selected.ps)]
	forest=forest.iloc[-1] ## selects the last to display. Is a choice, whatever
	forestsource.data=generate_forestplot_df(forest)

def callback_click(arg):
	if(arg == 0):
		taptool.callback = OpenURL(url="http://www.ensembl.org/Homo_sapiens/Variation/Explore?db=core;v=@rs;vdb=variation")
	else:
		taptool.callback = None
		mainplot_points.data_source.selected.on_change('indices', callback_pointclicked)


## adding callbacks to elements
control_meta.on_click(callback_meta)
control_source.on_click(callback_changesource)
control_signals.on_click(callback_toggleassoc)
rbg.on_click(callback_showsp)
chcolor.on_click(callback_chcolor)
burden.on_click(callback_burden)
control_click.on_click(callback_click)
#ld_hovertool=HoverTool(callback=callback_hover, renderers=[mainplot_points])
#ld_hovertool=HoverTool(callback=ld_hover, renderers=[mainplot_points])
#p1.add_tools(ld_hovertool)

#window=100000
chop=False
gene_plotter.linkedFeatures=linkedFeatures # "Linked_features.bed.gz"

p2 = gene_plotter.draw_genes(gc, window, width=1500, chop=chop)
p2.x_range=p1.x_range
p2.xaxis[0].formatter.use_scientific = False
p2.xaxis[0].major_label_text_font_size = "13pt"
p2.xaxis.axis_label = 'position on chromosome '+str(rawdat.chr[1].astype(np.int64))
p2.xaxis.axis_label_text_font_size = "15pt"


col_signals_table=[
TableColumn(field="chr", title="chr"),
TableColumn(field="ps", title="ps"),
TableColumn(field="a1", title="VCF A1"),
TableColumn(field="a2", title="VCF A2"),
TableColumn(field="ensembl_rs", title="RSid"),
TableColumn(field="ensembl_consequence", title="most severe consequence"),
TableColumn(field="weight", title="weight")
]
k = 0
for n in co_split:
	if(n=="meta"):
		col_signals_table.append(TableColumn(field="P-value_meta", title="P (meta)"))
		col_signals_table.append(TableColumn(field="Effect_meta", title="effect (meta)"))
		col_signals_table.append(TableColumn(field="StdErr_meta", title="S.E. (meta)"))
		col_signals_table.append(TableColumn(field="Allele1_meta", title="A1 (meta)"))
		col_signals_table.append(TableColumn(field="Allele2_meta", title="A2 (meta)"))
		col_signals_table.append(TableColumn(field="Freq1_meta", title="A2 (meta)"))
		col_signals_table.append(TableColumn(field="HetPVal_meta", title="het. P"))
	else:
		col_signals_table.append(TableColumn(field="p_score_"+n, title="P ("+n+")"))
		col_signals_table.append(TableColumn(field="beta_"+n, title="effect ("+n+")"))
		col_signals_table.append(TableColumn(field="se_"+n, title="S.E. ("+n+")"))
		col_signals_table.append(TableColumn(field="allele1_"+n, title="A1 ("+n+")"))
		col_signals_table.append(TableColumn(field="allele0_"+n, title="A2 ("+n+")"))
		col_signals_table.append(TableColumn(field="af_"+n, title="AF ("+n+")"))
	k = k + 1

signals_table=DataTable(source=metasegsource, columns=col_signals_table, width=1900)
p_sep=Div(text="""<h3>Signals (click on "show meta-analysis" and select points):</h3>""", width=1900)

#bbox=column(row([p_source, control_source]),row([p_rbg, rbg]), row([p_chcolor, chcolor]), row([p_burden, burden]), row([p_ld, control_ld]), row([p_signals, control_signals]), row([p_meta, control_meta]))
bbox=column(row([p_source, control_source]),row([p_rbg, rbg]), row([p_chcolor, chcolor]), row([p_burden, burden]), row([p_ld, control_ld]), row([p_signals, control_signals]), row([p_click, control_click]))

# Right hand side Burden p-value display table
burden_p_sep = Div(text="<h3>Burden p-values</h3>", height=50, width=430)
burden_p_data = ColumnDataSource(data = {n: ['{:0.5e}'.format(results[n])] for n in co_split})
burden_p_cols = [TableColumn(field=n, title=n) for n in co_split]
burden_p_table = DataTable(source=burden_p_data, columns=burden_p_cols, height=55, width=430)

l=layout([
[p_title],
[
[p1,p2],
[bbox, burden_p_sep, burden_p_table, p3]
],
[p_sep],
[signals_table]
])
# l=layout([[p1, bbox]])
#p1.output_backend = "svg" #NOT FUNCTIONAL
#p2.output_backend = "svg"
#save(l)
curdoc().add_root(l)
#rawdat.to_csv(output+".csv", index=False)
#ld.to_csv(output+".ld.csv",index=False)
