import json, requests, asr, subprocess, sys, requests, re
import urllib.request
import urllib.parse
import pandas as pd
import numpy as np
from pandas import notnull, isnull
from numpy import log10, append, nan
from numpy import frompyfunc
from mpmath import *


pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

# number by which to multiply the normal number of requests to Ensembl
# Increase if lots of 504 Timeout errors
ENSEMBL_USELESSNESS_COEFFICIENT=3

def info (*strs):
	outstr="[INFO]";
	for string in strs:
		outstr+=" "+str(string)
	print(outstr);
	return;

def warn (*strs):
	outstr="[WARNING]";
	for string in strs:
		outstr+=" "+str(string)
	print(outstr, file=sys.stderr);
	return;

class GeneCoordinates:
	chrom=0
	start=0
	end=0
	gene_id=""
	name=""
	def __init__(self, chrom, start, end, gene_id, name):
		self.chrom = chrom
		self.start = start
		self.end = end
		self.gene_id=gene_id
		self.name=name
		self.gstart = start
		self.gend = end

	def extend(self, margin):
		self.start-=int(margin);
		if self.start < 0: self.start = 0;
		self.end+=int(margin);


def get_csq_novel_variants(e, chrcol, pscol, a1col, a2col):
	global server
	e.loc[(e['ensembl_rs']=="novel") & (e[a1col]==e[a2col]),'ensembl_consequence']='double allele'
	novelsnps=e.loc[(e['ensembl_rs']=="novel") & (e['ensembl_consequence']!='double allele'),]
	info("Sending query for "+str(len(novelsnps.index))+" novel SNPs.")
	csq=pd.DataFrame()
	if novelsnps.empty:
		return e
	pd.options.mode.chained_assignment = None
	novelsnps['query']=novelsnps[chrcol].astype(str)+" "+novelsnps[pscol].astype(int).astype(str)+" . "+novelsnps[a1col].astype(str)+" "+novelsnps[a2col].astype(str)+" . . ."
	n=200
	for i in range(0, len(novelsnps['query']), n):
		request='{ "variants" : ["'+'", "'.join(novelsnps['query'][i:i+n])+'" ] }'
		ext = "/vep/homo_sapiens/region"
		headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
		info("\t\t\t🌐   Querying Ensembl VEP (POST) :"+server+ext)
		r = requests.post(server+ext, headers=headers, data=request)
		if not r.ok:
			print("headers :"+request)
			r.raise_for_status()
			sys.exit()
		jData = json.loads(r.text)
		csq=csq.append(pd.DataFrame(jData))

	# request='{ "variants" : ["'+'", "'.join(novelsnps['query'])+'" ] }'
	# ext = "/vep/homo_sapiens/region"
	# headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
	# info("\t\t\t🌐   Querying Ensembl VEP (POST) :"+server+ext)
	# r = requests.post(server+ext, headers=headers, data=request)

	# if not r.ok:
	#     print("headers :"+request)
	#     r.raise_for_status()
	#     sys.exit()

	# jData = json.loads(r.text)
	# csq=pd.DataFrame(jData)


	for index,row in csq.iterrows():
		e.loc[e['ps']==row['start'],'ensembl_consequence']=row['most_severe_consequence']

	e['ensembl_consequence'].replace('_', ' ')
	return e

def get_coordinates(gene_name):
	'''
	Function to return the genomic coordinates of a gene submitted
	output data: chromosome, start, end, stable ID
	'''

	global server
	ext = "/lookup/symbol/homo_sapiens/%s?expand=0" % (gene_name)

	r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

	if not r.ok:
	  r.raise_for_status()
	  sys.exit()

	decoded = r.json()
	gc=GeneCoordinates(int(decoded["seq_region_name"]), int(decoded["start"]), int(decoded["end"]), decoded["id"], gene_name)
	return(gc)


def get_rsid_in_region(gc):
	c=gc.chrom
	start=gc.start
	end=gc.end
	url = server+'/overlap/region/human/'+str(c)+":"+str(start)+"-"+str(end)+'?feature=variation;content-type=application/json;';
	info("\t\t\t🌐   Querying Ensembl with region "+url)
	response = urllib.request.urlopen(url).read().decode('utf-8')
	jData = json.loads(response)
	snps=pd.DataFrame(jData)
	snps['location']=snps.seq_region_name.map(str)+":"+snps.start.map(str)+"-"+snps.end.map(str)

	url = server+'/phenotype/region/homo_sapiens/'+str(c)+":"+str(start)+"-"+str(end)+'?feature_type=Variation;content-type=application/json;';
	info("\t\t\t🌐   Querying Ensembl with region "+url)
	response = urllib.request.urlopen(url).read().decode('utf-8')
	jData = json.loads(response)
	pheno=pd.DataFrame(jData)
	#print(pheno['phenotype_associations'])
	pheno['pheno']=""
	pheno['location']=""
	for index, variant in pheno.iterrows():
		for assoc in variant.phenotype_associations:
			if assoc['source'] != 'COSMIC':
				variant.pheno=assoc['description'] if (variant.pheno=="") else variant.pheno+";"+assoc['description'];
				variant.location=assoc['location']
	resp=pd.merge(snps, pheno, on='location', how='outer')
	resp.drop(["alleles", "assembly_name", "clinical_significance", "feature_type", "end", "seq_region_name", "phenotype_associations", "strand", "source", "id_y", "location"], axis=1, inplace=True)
	resp.dropna(inplace=True, subset=["id_x"])
	resp.rename(columns = {'start':'ps', 'id_x':'rs', 'consequence_type':'consequence'}, inplace = True)
	return(resp)

def get_rsid_in_region_old(gc):
	c=gc.chrom
	start=gc.start
	end=gc.end
	url = server+'/overlap/region/human/'+str(c)+":"+str(start)+"-"+str(end)+'?feature=variation;content-type=application/json;';
	info("Querying Ensembl with region "+url)
	response = urllib.request.urlopen(url).read().decode('utf-8')
	jData = json.loads(response)
	cat=pd.DataFrame(jData)
	nsplit=ENSEMBL_USELESSNESS_COEFFICIENT*int(cat.count()[0]/1000)
	resp=pd.DataFrame(columns=('rs', 'ps', 'consequence', 'pheno'))
	info("\t\t\t🌐   Performing "+str(nsplit)+" phenotype requests...");
	ext = "/variation/homo_sapiens?phenotypes=1"
	j=0
	for i in np.array_split(cat['id'], nsplit):
		j=j+1
		info("\t\t\t⌛   GWASCAT (POST) at "+server+ext+" ("+str(j)+'/'+str(nsplit)+")")
		header=pd.DataFrame()
		header['ids']=i
		data = header.to_dict(orient="list")
		headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
		data=json.dumps(data)
		data = data.encode('ascii')
		req = urllib.request.Request(server+ext, data=data, headers=headers)
		with urllib.request.urlopen(req) as response:
			r=response.read().decode('utf-8')
			jData = json.loads(r)
			#pp.pprint(jData)
			for rsid in jData:
				ps=jData[rsid]['mappings'][0]['end']
				consequence=jData[rsid]['most_severe_consequence']
				pheno=""
				for k in jData[rsid]['phenotypes']:
					pheno=k['trait'] if (pheno=="") else pheno+";"+k['trait'];
				resp=resp.append({'rs':rsid, 'ps':int(ps), 'consequence':consequence, 'pheno':pheno}, ignore_index=True)
	return(resp)

def fetch_single_point(gc, sp_results):
	c=gc.chrom
	start=gc.start
	end=gc.end
	sp = pd.DataFrame();
	task = subprocess.Popen(["tabix", "-h", sp_results, str(c)+":"+str(start)+"-"+str(end)], stdout=subprocess.PIPE);
	sp=pd.read_table(task.stdout);
	task = subprocess.Popen(["zgrep", "-m1", "chr", sp_results], stdout=subprocess.PIPE);
	sp.columns=task.stdout.read().decode('UTF-8').split();
	return(sp)

def fetch_single_point_meta(gc, sp_results, co_names):
	c=gc.chrom
	start=gc.start
	end=gc.end
	spfiles=sp_results.split(",")
	co_names=co_names.split(",")
	co_names.append("meta")
	if len(spfiles)!=len(co_names):
		sys.exit("cohort and single-point information not coherent (both must be comma-separated)")
	i=0
	import os.path
	retdf=pd.DataFrame()
	for file in spfiles:
		if not os.path.isfile(file):
			sys.exit("Single-point file "+file+" is not reachable.")
		sp = pd.DataFrame();
		# def logmp(x):
		# 	return np.float(log10(mpf(x)))
		# getlogp = frompyfunc(logmp, 1, 1)
		try:
			task = subprocess.Popen(["zgrep", "-m1", "chr", file], stdout=subprocess.PIPE);
			cols = task.stdout.read().decode('UTF-8').replace('#', '').split();
			task = subprocess.Popen(["tabix", file, str(c)+":"+str(start)+"-"+str(end)], stdout=subprocess.PIPE);
			if(co_names[i]=="meta"):
				sp=pd.read_table(task.stdout, header=0, names=cols, dtype={"P-value":np.unicode_});
				sp["logp"]=[-1*float(str(log10(mpf(x)))) for x in sp["P-value"]]
			else:
				sp=pd.read_table(task.stdout, header=0, names=cols, dtype={"p_score":np.unicode_});
				sp["logp"]=[-1*float(str(log10(mpf(x)))) for x in sp["p_score"]]
		except:
			e = sys.exc_info()[0]
			info(e)
			sys.exit("Error while running tabix on single-point file "+file+" is not reachable.")

		if sp.empty:
			sys.exit("Tabixing "+file+"for region "+str(c)+":"+str(start)+"-"+str(end)+" returned 0 rows.")
		if(co_names[i]=="meta"):
			sp.rename(columns={"MarkerName" : "ps"}, inplace=True)
		sp=sp.add_suffix(co_names[i])
		if retdf.empty:
			retdf=sp
			retdf['ps']=retdf["ps"+co_names[i]]
			retdf['chr']=retdf["chr"+co_names[i]]
			retdf['allele0']=retdf["allele0"+co_names[i]]
			retdf['allele1']=retdf["allele1"+co_names[i]]
		else:

			retdf=pd.merge(retdf, sp, left_on="ps", right_on="ps"+co_names[i], how="outer")
			## ALSO THIS WOULD WORK df.loc[df['foo'].isnull(),'foo'] = df['bar']
			retdf.ps.fillna(retdf["ps"+co_names[i]], inplace=True)
			retdf.chr.fillna(retdf["chr"+co_names[i]], inplace=True)
		i=i+1
	return(retdf)

def read_large_results_file(fn, gene,protein, condition_string):
	task=subprocess.Popen(["bzgrep", "-w", gene+"."+condition_string, fn], stdout=subprocess.PIPE)
	results=pd.read_table(task.stdout, header=None, names=["protein","group","n_variants","miss_min","miss_mean","miss_max","freq_min","freq_mean","freq_max","B_score","B_var","B_pval","S_pval","O_pval","O_minp","O_minp.rho","E_pval"]);
	return(results[results.protein==protein])

def read_sc_results_file(fn, gene,pheno, condition_string):
	task=subprocess.Popen(["zgrep", "-w", "^"+gene, fn], stdout=subprocess.PIPE)
	results=pd.read_table(task.stdout, header=None, names=["gene","pheno","condition","symbol","n_variants","miss_min","miss_mean","miss_max","freq_min","freq_mean","freq_max","B_score","B_var","B_pval","S_pval","O_pval","O_minp","O_minp.rho","E_pval"]);
	results=results[(results.pheno==pheno) & (results.condition ==condition_string)]
	return(results.O_minp.iloc[0])

def read_meta_results_file(fn, gene,pheno, condition_string):
	task=subprocess.Popen(["zgrep", "-w", "^"+gene, fn], stdout=subprocess.PIPE)
	results=pd.read_table(task.stdout, header=None, names=["gene","pheno","condition","symbol","n_variants","B_score","B_var","B_pval","S_pval","O_pval","O_minp","O_minp.rho","E_pval"]);
	results=results[(results.pheno==pheno) & (results.condition ==condition_string)]
	return(results.O_minp.iloc[0])

def read_burden_ps(co_names, smmat_out_file, ensid, pheno, condition_string):
	co_names=co_names.split(",")
	smmat_out_file=smmat_out_file.split(",")
	if len(smmat_out_file)!=(len(co_names)+1):
		print(co_names)
		print(smmat_out_file)
		sys.exit("cohort and meta-analysis information not coherent (both must be comma-separated)")
	i=0
	burdp={}
	for file in smmat_out_file:
		if(i==len(co_names)):
			burdp["meta"]=read_meta_results_file(file, ensid, pheno, condition_string)
		else:
			burdp[co_names[i]]=read_sc_results_file(file, ensid, pheno, condition_string)
		i=i+1
	return(burdp)

def read_variants_from_gene_set(gc, input_monster):
	c=gc.chrom
	variants=pd.read_table(input_monster)
	variants.columns=[w.replace('chr'+str(c), '') for w in variants.columns]
	variants.columns=[re.sub("[A-Z]*", '', w) for w in variants.columns]
	variants.drop(variants.columns[[0,1]], axis=1, inplace=True)
	variants=variants.transpose()
	variants['ps']=variants.index
	variants.index=range(variants.count()[0])
	if (len(variants.columns)==2):
		variants.columns=["weight", "ps"]
	else:
		#Sometimes runs have no weights
		variants.columns=["ps"]
		variants['weight']=1
	variants.ps=pd.to_numeric(variants.ps, errors='coerce')
	return(variants)


def read_variants_from_gene_set_SMMAT(gene, condition_string, smmat_set_file):
	variantset=pd.read_table("/storage/hmgu/projects/SCALLOP_WGS_meta/from.farm/HA.HP.OY.variantset.noInf.txt", header=None, names=["set", "chr", "ps", "a1", "a2", "weight"])
	#returns df with columns above
	variantset.drop(['chr'], axis=1, inplace=True)
	return(variantset[(variantset.set==gene+"."+condition_string)])

def produce_meta_df(gc, sp, variants, vcf_files, co_names):
	#      chrMANOLIS rsMANOLIS  psMANOLIS  n_missMANOLIS allele1MANOLIS allele0MANOLIS  afMANOLIS  betaMANOLIS  seMANOLIS  l_remleMANOLIS  l_mleMANOLIS
	#p_waldMANOLIS  p_lrtMANOLIS  p_scoreMANOLIS           ps  chr allele0 allele1  chrPomak      rsPomak      psPomak  n_missPomak allele1Pomak allele0Po
	#mak  afPomak  betaPomak  sePomak  l_remlePomak  l_mlePomak   p_waldPomak    p_lrtPomak  p_scorePomak  chrmeta       psmeta Allele1meta Allele2meta  F
	#req1meta  FreqSEmeta  MinFreqmeta  MaxFreqmeta  Effectmeta  StdErrmeta   P-valuemeta Directionmeta  HetISqmeta  HetChiSqmeta  HetDfmeta  HetPValmeta
	# NMISSTOTALmeta       consequence pheno    ensembl_rs ensembl_consequence
	c=gc.chrom
	start = gc.start
	end = gc.end
	gene_start=gc.gstart
	gene_end=gc.gend
	ensid=gc.gene_id
	rawdat=pd.merge(sp, variants, on='ps', how='outer')
	if rawdat[rawdat.chr.isnull()].ps.size > 0 :
		warn(str(rawdat[rawdat.chr.isnull()].ps.size)+" variants from the gene set were not found in the single point.")
	rawdat.dropna(subset=['chr'], inplace=True)

	info("Calculating LD...")
	info("getld_meta.sh", co_names, vcf_files, "chr"+str(c)+":"+str(start)+"-"+str(end), str(sp.size), str(end-start))
	import os
	import subprocess
	task = subprocess.Popen([contdir+"/getld_meta.sh", co_names, vcf_files, "chr"+str(c)+":"+str(start)+"-"+str(end), str(sp.size), str(end-start)], stdout=subprocess.PIPE);
	print(task.stderr)
	ld=pd.read_table(task.stdout, sep='\s+');
	info("Computed LD between ", str(len(ld.index)), " variant pairs.")

	## Defining plot-specific data
	info("Defining plot-specific data...")
	rawdat['radii']=3
	denom=rawdat.weight[rawdat.weight.notnull()]
	if len(denom):
		denom=max(denom)
	else:
		denom=1
	rawdat.loc[rawdat.weight.notnull(), 'radii']=3+20*rawdat.weight[rawdat.weight.notnull()]/denom
	rawdat['alpha']=0
	rawdat.loc[rawdat.weight.notnull(), 'alpha']=0.8
	rawdat['alpha_prevsig']=0
	rawdat.loc[(rawdat.pheno!="none") & (rawdat.alpha>0), 'alpha_prevsig']=1
	from bokeh.palettes import PuOr8 as palette
	from bokeh.palettes import Viridis8 as palWeight
	# Spectral9 Palette : ['#3288bd', '#66c2a5', '#abdda4', '#e6f598', '#ffffbf', '#fee08b', '#fdae61', '#f46d43', '#d53e4f']
	palWeight=[x for x in palWeight]
	palWeight.append("#939393")
	rawdat['maf']=[af if af<0.5 else 1-af for af in rawdat.Freq1meta]
	rawdat['mafcolor']=[palette[i] for i in pd.cut(rawdat.maf, [-1, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.6]).cat.codes]
	rawdat['color']="#1F77B4"
	rawdat['weightcolor']=[palWeight[i] for i in pd.cut(rawdat.weight, 7).cat.codes]
	rawdat['outcolor']="#3288bd"
	rawdat["outalpha"]=0
	return(rawdat, ld)



def produce_single_cohort_df(gc, sp_results, resp, vcf, smmat_out_file, smmat_set_file, pheno, condition_string, coname, megasp, variants):
	from pandas import notnull, isnull
	from numpy import log10, append, nan
	global contdir
	c=gc.chrom
	start = gc.start
	end = gc.end
	gene_start=gc.gstart
	gene_end=gc.gend
	ensid=gc.gene_id
	## Get the single point results
	sp=megasp.loc[:, megasp.columns.str.endswith(coname) | megasp.columns.isin(['consequence', 'pheno', 'ensembl_rs', 'ensembl_consequence'])]
	sp.columns=sp.columns.str.replace(coname+'$', '')
	sp=sp[sp.chr.notnull()]
	## Get the weights and variants in burden


	rawdat=pd.merge(sp, variants, on='ps', how='outer')
	if rawdat[rawdat.chr.isnull()].ps.size > 0 :
		warn(str(rawdat[rawdat.chr.isnull()].ps.size)+" variants from the gene set were not found in the single point.")
	rawdat.dropna(subset=['chr'], inplace=True)



	## Calculate LD
	info("Calculating LD...")
	info("getld.sh", vcf, "chr"+str(c)+":"+str(start)+"-"+str(end), str(sp.size), str(end-start))
	import os
	import subprocess
	task = subprocess.Popen([contdir+"/getld.sh", vcf, "chr"+str(c)+":"+str(start)+"-"+str(end), str(sp.size), str(end-start)], stdout=subprocess.PIPE);
	ld=pd.read_table(task.stdout, sep='\s+');
	os.remove("plink.log")
	os.remove("plink.nosex")



	## Defining plot-specific data
	info("Defining plot-specific data...")
	rawdat['radii']=3
	denom=rawdat.weight[rawdat.weight.notnull()]
	if len(denom):
		denom=max(denom)
	else:
		denom=1
	rawdat.loc[rawdat.weight.notnull(), 'radii']=3+20*rawdat.weight[rawdat.weight.notnull()]/denom
	rawdat['alpha']=0
	rawdat.loc[rawdat.weight.notnull(), 'alpha']=0.8
	rawdat['alpha_prevsig']=0
	rawdat.loc[(rawdat.pheno!="none") & (rawdat.alpha>0), 'alpha_prevsig']=1
	from bokeh.palettes import PuOr8 as palette
	from bokeh.palettes import Viridis8 as palWeight
	# Spectral9 Palette : ['#3288bd', '#66c2a5', '#abdda4', '#e6f598', '#ffffbf', '#fee08b', '#fdae61', '#f46d43', '#d53e4f']
	palWeight=[x for x in palWeight]
	palWeight.append("#939393")
	rawdat['maf']=[af if af<0.5 else 1-af for af in rawdat.af]
	rawdat['mafcolor']=[palette[i] for i in pd.cut(rawdat.maf, [-1, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.6]).cat.codes]
	rawdat['color']="#1F77B4"
	rawdat['weightcolor']=[palWeight[i] for i in pd.cut(rawdat.weight, 7).cat.codes]
	rawdat['outcolor']="#3288bd"
	rawdat["outalpha"]=0
	return(rawdat, ld)
