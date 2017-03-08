import json, requests, asr, subprocess, sys, requests, re
import urllib.request
import urllib.parse
import pandas as pd
import numpy as np

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

	def extend(self, margin):
		self.start-=int(margin);
		self.end+=int(margin);


def get_coordinates(gene_name):
    '''
    Function to return the genomic coordinates of a gene submitted (GRCh37)
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

def read_variants_from_gene_set(gc, input_monster):
	c=gc.chrom
	variants=pd.read_table(input_monster)
	variants.columns=[w.replace('chr'+str(c)+"_", '') for w in variants.columns]
	variants.columns=[re.sub("_.*", '', w) for w in variants.columns]
	variants.drop(variants.columns[[0,1]], axis=1, inplace=True)
	variants=variants.transpose()
	variants['ps']=variants.index
	variants.index=range(variants.count()[0])
	variants.columns=["weight", "ps"]
	variants.ps=pd.to_numeric(variants.ps, errors='coerce')
	return(variants)