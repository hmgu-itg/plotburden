#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import random
import pickle

import click
import numpy as np
import pandas as pd
import seaborn as sns


import gene_plotter
import helper_functions
from helper_functions import info


__version__ = '0.0.1'

@click.command()
@click.option('-p', '--pheno', type = click.STRING, required=True, help = 'Phenotype name')
@click.option('-g', '--gene',type = click.STRING, required=True, help = 'Gene name')
@click.option('-c', '--condition', 'condition_string', type = click.STRING, required=True, help = 'Condition')
@click.option('-w', '--window', type = click.INT, default = 100_000, show_default=True, help = 'Window size')
@click.option('--variant-set', type = click.Path(exists=True), required=True, help = 'Variant set file')
@click.option('--cohort-name', type = click.STRING, required=True, multiple=True, help = 'Single cohort name')
@click.option('--cohort-rv', type = click.Path(exists=True, dir_okay=False), required=True, multiple=True, help = 'Single cohort rare variant file')
@click.option('--cohort-sp', type = click.Path(exists=True, dir_okay=False), required=True, multiple=True, help = 'Single cohort single-point file')
@click.option('--cohort-vcf', type = click.Path(exists=True, dir_okay=False), required=True, multiple=True, help = 'Single cohort VCF file')
@click.option('-o', '--out', 'output', type = click.Path(exists=False, dir_okay=False), required=True, help = 'Output prefix')
@click.option('-l', '--linked', 'linkedFeatures', type = click.Path(exists=True, dir_okay=False), required=True, help = 'Filepath to LinkedFeatures file')
@click.option('--chop/--no-chop', default=False, show_default=True, help = 'Whether to chop or not')
@click.version_option(__version__)
def cli(pheno, gene, condition_string, window, variant_set, cohort_name, cohort_rv, cohort_sp, cohort_vcf, output, linkedFeatures, chop, **kwargs):
    '''
    Prepares data for plotting
    '''
    cohort_data = list(zip(cohort_name, cohort_rv, cohort_sp, cohort_vcf))

    # smmat_set_file=sys.argv[4]
    # smmat_out_file=sys.argv[5]
    # co_names=sys.argv[6]
    # sp_results=sys.argv[7]
    # vcf=sys.argv[8]


    gene_plotter.linkedFeatures=linkedFeatures
    helper_functions.contdir=os.path.dirname(__file__)

    # Getting Gene coordinates and extend the coordinates with the specified window:
    info(f"Querying Ensembl for coordinates of {gene}...")
    gc = helper_functions.get_coordinates(gene)
    gc.extend(int(window))

    #Extract coordinates:

    c = gc.chrom
    start = gc.start
    end = gc.end
    ensid=gc.gene_id

    # Report coordinates:
    info(f"    ⇰ Ensembl provided the coordinates chr{c}:{start}-{end} for gene {gene}")
    info(f"    ⇰ Plot boundaries: chr{c}:{start}-{end}")

    ## Getting variant consequences for all variants in the region
    info("Querying Ensembl for SNP consequences and phenotype associations.")
    resp = helper_functions.get_rsid_in_region(gc)

    #resp.to_csv(gene+".snp.data", index=None, sep=",", quoting=csv.QUOTE_NONNUMERIC)
    #resp=pd.read_table("snp.data", sep=",")
    resp['pheno'].replace(to_replace="Annotated by HGMD*", value="", inplace=True, regex=True)
    resp['pheno'].replace(to_replace="ClinVar.*not specified", value="", inplace=True, regex=True)
    resp.loc[pd.isnull(resp.pheno), 'pheno']="none"
    resp.pheno=resp.pheno.str.strip()
    resp.loc[resp.pheno=="", 'pheno']="none"
    resp['ensembl_rs']=resp['rs']
    resp.drop('rs', axis=1, inplace=True)
    info(f"    ⇰ Ensembl provided {len(resp)} known SNPs, {len(resp[resp.pheno!='none'])} have associated phenotypes.")
    return



    ## Get the single point results. Returns one merged (outer) dataframe with all columns suffixed by the cohort name
    ## We are just going to use this for annotation purposes
    sp = helper_functions.fetch_single_point_meta(gc, sp_results, co_names)
    info("Read", len(sp), "lines from single-point analysis.")
    sp=pd.merge(sp, resp, on='ps', how='outer')
    sp.loc[pd.isnull(sp.ensembl_rs), 'ensembl_rs']="novel"
    sp.loc[pd.isnull(sp.consequence), 'consequence']="novel"
    sp=sp[pd.notnull(sp.chr)]
    sp['ensembl_consequence']=sp['consequence']
    sp['chr']=sp['chr'].astype(int)
    info("Getting consequences for novel variants...")
    sp = helper_functions.get_csq_novel_variants(sp, 'chr', 'ps', 'allele0', 'allele1')



    ## Get the burden p-values
    ## Returns a dictionary indexed by the cohort name or "meta"
    info("Reading burden P-values...")
    results = helper_functions.read_burden_ps(co_names, smmat_out_file, ensid, pheno, condition_string)
    # import pickle
    # with open('results.bin', 'rb') as config_dictionary_file:
    #     results = pickle.load(config_dictionary_file)

    ## read all variants in the gene sets including those not in some cohorts
    info("Reading variants from gene set...")
    variants = helper_functions.read_variants_from_gene_set_SMMAT(ensid, condition_string, smmat_set_file)
    info("Read ", variants.count()[0], "variants in burden across all cohorts")
    ## Now for the plot data
    ##
    cohdat={}
    lddat={}
    i=0
    import pickle
    for n in co_names.split(","):
        info("Preparing plot data for cohort "+n)
        rawdat, ld=helper_functions.produce_single_cohort_df(gc, sp_results.split(',')[i], resp, vcf.split(",")[i], smmat_out_file.split(',')[i], smmat_set_file, pheno, condition_string, n, sp, variants)
        cohdat[i]=rawdat
        lddat[i]=ld
    #    cohdat[i]=ColumnDataSource(data=dict(ps=rawdat.ps, logsp=-log10(rawdat.p_score), radii=rawdat.radii, alpha=rawdat.alpha, color=rawdat.color, mafcolor=rawdat.mafcolor, weightcolor=rawdat.weightcolor, outcol=rawdat.outcolor, outalpha=rawdat.outalpha, alpha_prevsig=rawdat.alpha_prevsig, snpid=rawdat.rs, rs=rawdat.ensembl_rs, maf=rawdat.maf, csq=rawdat.ensembl_consequence))
        i=i+1

    cohdat[i], lddat[i] = helper_functions.produce_meta_df(gc, sp, variants, vcf, co_names)
    with open('sp.bin', 'wb') as config_dictionary_file:
        pickle.dump(sp, config_dictionary_file)


    # with open('resp.bin', 'wb') as config_dictionary_file:
    #     pickle.dump(resp, config_dictionary_file)
    #
    # with open('ld.bin', 'wb') as config_dictionary_file:
    #     pickle.dump(lddat, config_dictionary_file)
    #
    # with open('results.bin', 'wb') as config_dictionary_file:
    #     pickle.dump(results, config_dictionary_file)
    # with open('gc.bin', 'wb') as config_dictionary_file:
    #     pickle.dump(gc, config_dictionary_file)
    #
    # # =================================================================================================
    # # =================================================================================================
    # # =================================================================================================
    #
    # import pickle
    #
    # with open('cohdat.bin', 'rb') as config_dictionary_file:
    #     cohdat = pickle.load(config_dictionary_file)
    #
    # with open('gc.bin', 'rb') as config_dictionary_file:
    #     gc = pickle.load(config_dictionary_file)
    # c=gc.chrom
    # start = gc.start
    # end = gc.end
    # gene_start=gc.gstart
    # gene_end=gc.gend
    # ensid=gc.gene_id
    #
    #
    # with open('ld.bin', 'rb') as config_dictionary_file:
    #     lddat = pickle.load(config_dictionary_file)
    #
    # with open('results.bin', 'rb') as config_dictionary_file:
    #     results = pickle.load(config_dictionary_file)
    #
    # with open('resp.bin', 'rb') as config_dictionary_file:
    #     resp = pickle.load(config_dictionary_file)

    ## Preparing plot variables
    ## maxlogp is the min of all ps, that is used to set plot ylim
    ## cohdat is overwritten to become an array of dataframes containing the SP info for the respective cohorts

    i=0
    rawdats=[]
    maxlogp=0
    # Single-points
    for n in co_names.split(","):
        rawdats.append(cohdat[i])
        rawdat=cohdat[i]
    #    if -log10(rawdat.p_score.min())>maxlogp:
        if rawdat.logp.max()>maxlogp:
            print(n)
            print(rawdat.columns)
            maxlogp=rawdat.logp.max()
        cohdat[i]=dict(ps=rawdat.ps, p_value=rawdat.p_score, logsp=rawdat.logp, radii=rawdat.radii, alpha=rawdat.alpha, color=rawdat.color, mafcolor=rawdat.mafcolor, weightcolor=rawdat.weightcolor, outcol=rawdat.outcolor, outalpha=rawdat.outalpha, alpha_prevsig=rawdat.alpha_prevsig, snpid=rawdat.rs, rs=rawdat.ensembl_rs, maf=rawdat.maf, csq=rawdat.ensembl_consequence)
        i=i+1

    ## meta-analysis (beware, ALL columns are as is)
    rawdats.append(cohdat[i])
    rawdat=cohdat[i]
    print(rawdat.columns)
    if rawdat.logpmeta.max()>maxlogp:
        maxlogp=rawdat.logpmeta.max()
    rawdat.rename(columns = {'P-valuemeta': 'p_score'}, inplace = True) # WARNING! This code is only valid for METAL output meta-analysis file!
    cohdat[i]=dict(ps=rawdat.ps, p_score=rawdat.p_score, logsp=rawdat.logpmeta, radii=rawdat.radii, alpha=rawdat.alpha, color=rawdat.color, mafcolor=rawdat.mafcolor, weightcolor=rawdat.weightcolor, outcol=rawdat.outcolor, outalpha=rawdat.outalpha, alpha_prevsig=rawdat.alpha_prevsig, snpid=rawdat.chr.astype(str)+":"+rawdat.ps.astype(str), rs=rawdat.ensembl_rs, maf=rawdat.maf, csq=rawdat.ensembl_consequence)


    ## Creating the df containing ALL points coloured by cohort
    co_split=co_names.split(",")
    co_split.append("meta")

    palette=sns.color_palette("Set2").as_hex()
    random.shuffle(palette)
    #palette.insert(0, '#3288bd')
    cter=0
    bigdf=None
    for data in cohdat:
        df=pd.DataFrame(cohdat[data])
        # set color to cohort-based color
        if data==(len(cohdat)-1):
            df['cocolor']='#3288bd'
        else:
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
    # for co in co_split:
    #     if co=="meta":
    #         col='P-valuemeta'
    #     else:
    #         col='p_score'+co
    #     sp['logp'+co]=-log10(sp[col])


    sp['minp']=sp[["logp"+s for s in co_split]].min(axis=1)
    sp['maxp']=sp[["logp"+s for s in co_split]].max(axis=1)
    sp['segcol']="gray"
    sp=pd.merge(sp, variants, on="ps", how='outer')
    sp.dropna(subset=['chr'], inplace=True)
    #sp.loc[sp.weight.isnull(), 'segcol']="#FFFFFF00"


    # with open('cohdat.bin', 'wb') as config_dictionary_file:
    #     pickle.dump(cohdat, config_dictionary_file)
    # with open('rawdat.bin', 'wb') as config_dictionary_file:
    #     pickle.dump(rawdats, config_dictionary_file)

    ## we modified the whole code to have logp calculated within a helper Function
    ## burden logp very unklikely to be that low so we leave logp calculation as is
    rawdat=rawdats[0]
    if -np.log10(min(results.values()))>maxlogp:
        maxlogp=-np.log10(min(results.values()))


    plotdat=dict(rawdats=rawdats, rawdat=rawdat, maxlogp=maxlogp, gene=gene, gc=gc, resp=resp, lddat=lddat, sp=sp, cohdat=cohdat, co_split=co_split, results=results, bigdf=bigdf, window=window, chop=chop, pheno=pheno, condition_string=condition_string, linkedFeatures=linkedFeatures)
    with open(output+'.plotdat.bin', 'wb') as config_dictionary_file:
        pickle.dump(plotdat, config_dictionary_file)


if __name__ == "__main__":
    cli()