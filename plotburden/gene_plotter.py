#!/usr/bin/env python3
import sys
import json

import requests
import pandas as pd
from pybedtools import BedTool
from bokeh.plotting import figure
from bokeh.models import HoverTool, BoxAnnotation, Label

from .helper_functions import info


def draw_genes(gc, window, width=900, height=400, chop = "No"):
    '''
    Based on the gene name this function draws selected and overlapping genomic features on
    a bokeh image object created based on the submitted width and height values.
    '''


    # If we are chomping the gene track, we have to cut back the height:
    if chop == "Yes":
        height = 180

    # Based on the gene name get all selected exons and regulatory features that are associated with this gene:
    df = get_genomic_features(gc)
    info('Selected genomic features extracted.')
    # Get the boundaries of the used features:
    start_region = df.start.min() - int(window)
    start_region = 0 if start_region < 0 else start_region
    end_region = df.end.max() + int(window)
    chromosome = df.chrom.tolist()[0]

    # Return overlapping genomic features:
    overlapping_df = get_overlapping_features(chromosome, start_region, end_region)
    info('Overlapping genomic features extracted.')

    # Get rank of each genomic features:
    gene_rank = get_gene_rank(overlapping_df[overlapping_df.feature_type.isin(["gene"])])
    regulatory_rank = get_gene_rank(overlapping_df[overlapping_df.feature_type.isin(['regulatory'])])

    # Currently that's the lowest y values of the plot:
    y_pos_min = df.y_position.min() - 1.5
    gene_ypos = pd.Series([ y_pos_min - 0.5 *gene_rank[x]
                           for x in overlapping_df[overlapping_df.feature_type != "regulatory"]["name"]],
        index = overlapping_df[overlapping_df.feature_type != "regulatory"]["name"].index)
    y_pos_min = gene_ypos.min() - 1
    regulatory_ypos = pd.Series([ y_pos_min - 0.3 * regulatory_rank[x]
                                 for x in overlapping_df[overlapping_df.feature_type == "regulatory"]["name"]],
        index = overlapping_df[overlapping_df.feature_type == "regulatory"]["name"].index)
    Y_series = gene_ypos.append(regulatory_ypos, ignore_index=False).sort_index()
    overlapping_df["y_position"] = Y_series

    # Adding biotype based color to dataframe:
    overlapping_df["color"] = [get_biotype_color(x) for x in overlapping_df["biotype"]]

    # Adjustig y position according to the chop factor:
    if chop == "Yes":
        Y_min = 8.5
    else:
        Y_min = overlapping_df["y_position"].min()-0.5

    ##
    ## Based on the returned features let's calculate the boundaries of the plots:
    ##
    #tools = [WheelZoomTool(), PanTool(), ResetTool()]
    p = figure(width=width, height=height,
               y_range = (Y_min,11),
               x_range = (start_region, end_region), tools = "xwheel_zoom,xpan,reset,save")
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    p.yaxis.visible = False

    # Drawing highlight for the selected features:
    highlight_box = BoxAnnotation(bottom=8.5, fill_alpha=0.6, fill_color='#FFECD2',level='glyph',)
    p.renderers.extend([highlight_box])

    # Drawing selected features:
    selected_exons = draw_box(p, df[df.feature_type == 'exon'].copy())
    selected_genes = draw_line(p, df[df.feature_type == 'gene'].copy())
    selected_regulatory = draw_box(p, df[df.feature_type == 'regulatory'].copy(), height=0.2)

    # Drawing overlapping features:
    if chop != "Yes":
        overlapping_exons = draw_box(p, overlapping_df[overlapping_df.feature_type == 'exon'].copy())
        overlapping_genes = draw_line(p, overlapping_df[overlapping_df.feature_type == 'gene'].copy())
        overlapping_regulatory = draw_box(p, overlapping_df[overlapping_df.feature_type == 'regulatory'].copy(), height=0.2)

        # Adding object with tooltip:
        hover = HoverTool(renderers = [selected_exons, selected_genes, selected_regulatory,
                          overlapping_exons, overlapping_genes, overlapping_regulatory],
            tooltips=[("Name:", "@name"),("Biotype:", "@biotype")])
        p.add_tools(hover)

        # Adding text:
        overlapping_text = Label(x=65, y=8, x_units='screen', y_units='data',text='Overlapping genomic features')
        p.add_layout(overlapping_text)
    else:
        # Adding object with tooltip:
        hover = HoverTool(renderers = [selected_exons, selected_genes, selected_regulatory],
                         tooltips=[("Name:", "@name"),("Biotype:", "@biotype")])
        p.add_tools(hover)

    # Adding text to the plot:
    selected_text = Label(x=65, y=10.5, x_units='screen', y_units='data',text='Selected genomic features')
    p.add_layout(selected_text)

    return(p)

def get_overlapping_features(chromosome, start, end):
    '''
    This function retrieves a list of overlapping genomic features based on the provided
    coordinates.

    Retruned features: exon, gene, regulatory features.

    Input: chromosome, start, end
    Output: dataframe
        columns: chromosome, start, end, feature type, info (for gene and exon: gene name, or regulatory type)
    '''
    # Initalize request:
    server = "http://rest.ensembl.org"
    ext = "/overlap/region/human/%s:%s-%s?" %(chromosome, start, end)
    features = "feature=exon&feature=transcript&feature=gene&feature=regulatory"

    r = requests.get(server+ext+features, headers={ "Content-Type" : "application/json"})

    if not r.ok:
      r.raise_for_status()
      sys.exit()

    gene_name_mapping = {}
    feature_container = []

    info('List of overlapping features retrieved.')

    # Looping through all overlapping features:
    for feature in  r.json():

        # Extract relevant fields:
        seq_region_name = feature["seq_region_name"]
        start = feature["start"]
        end = feature["end"]
        ft = feature["feature_type"]

        # Extract eature types and transcript and gene IDs:
        if ft == "regulatory":
            infocol = feature["description"]
            #print(feature);
            feature_container.append([seq_region_name, start, end, ft, feature["id"], feature["id"], infocol])
        elif ft == "exon":
            infocol = feature["Parent"]
            feature_container.append([seq_region_name, start, end, ft, infocol,"", ""])
        elif ft == "transcript":
            gene_name_mapping[feature["id"]] = feature["Parent"]
        elif ft == "gene":
            try:
                infocol = feature["external_name"]
            except KeyError:
                infocol = feature["gene_id"]
                
            gene_name_mapping[feature["id"]] = [infocol, feature["biotype"]]
            feature_container.append([seq_region_name, start, end, ft, infocol, feature["id"], feature["biotype"]])

    info('List of overlapping features processed.')

    # Assign gene name for exons:
    for ft in feature_container:
        if ft[3] == "exon":
            g_ID = gene_name_mapping[ft[4]]
            g_name = gene_name_mapping[g_ID][0]
            g_biotype = gene_name_mapping[g_ID][1]
            ft[4] = g_name
            ft[6] = g_biotype

    # Return data frame:
    overlapping_df = pd.DataFrame(data=feature_container,
                                  columns=['chromosome', 'start', 'end', 'feature_type', 'name', 'ID', 'biotype'],
                                  index=list(range(len(feature_container))))

    # Adding extra column that tells in which row a feature should be plotted:
    #rank = get_gene_rank(overlapping_d[overlapping_d.feature_type == "gene"])

    return(overlapping_df)

def get_gene_rank(df):
    '''
    This function ranks the genes from 0 to n to make sure gene with the same
    rank won't overlap. It helps to avoid overlapping genes on the plot.

    Returned values: ranks (dict), pos (array)
    '''

    ranks = {}
    pos = []

    for index, row in df.iterrows():

        # For the first row we initalize the array:
        if len(pos) == 0:
            pos = [row["end"]]
            ranks[row["name"]] = 0
            continue

        # looping through all indices of the array until we found a smaller one:
        flag = 0
        for i in range(len(pos)):
            if pos[i] < row["start"]:
                pos[i] = row["end"]
                ranks[row["name"]] = i
                flag = 1
                break

        if flag == 0:
            pos.append(row["end"])
            ranks[row["name"]] = len(pos) - 1

    return(ranks)

def get_genomic_features(gc):
    '''
    Returns a dataframe read from a logfile provided as input.
    '''
    gene_name=gc.name
    # the extracted features:
    gencode_features = ['exon']
    regulatory_features = ['promoter', 'enhancer', 'TF_binding_site']

    # Genomic feature file:
    #featureFile = '/lustre/scratch113/projects/helic/ds26/project_burden/2016.10.10/Linked_features.bed.gz'
    #featureFile = '/lustre/scratch119/humgen/projects/helic/ds26/project_burden/2016.10.10/Linked_features.bed.gz'
    global linkedFeatures
    featureFile = BedTool(linkedFeatures)

    # Extract all lines overlapping with this gene:
    (chromosome, start, end, gene_ID, name) = (gc.chrom, gc.gstart, gc.gend, gc.gene_id, gc.name)
    info('Coordinates of the queried gene retrieved.')

    # Creating bed formatted text from the submitted genomic coordinates:
    #print(str(chromosome), start, end, name)
    geneBed = BedTool([(str(chromosome), start, end, name)])

    # Overlapping genomic features are selected:
    intersectfeature = geneBed.intersect(featureFile, wb = True, sorted = True).to_dataframe()
    intersectfeature.columns = ['chrom', 'start', 'end', 'name', 'chrom2', 'start2', 'end2',
           'geneID', 'annot']

    # Open and read files:
    bedlines = {"GENCODE" : [], "regulatory" : []}
    for line in intersectfeature[intersectfeature.geneID == gc.gene_id]["annot"]:
        # content = line.strip()

        # Extract json string and load data:
        jsonData = json.loads(line)

        # Creating bedfile:
        if jsonData["source"] == "GENCODE":
            if jsonData["class"] in gencode_features:
                bedlines["GENCODE"].append((jsonData["chr"], jsonData["start"], jsonData["end"], jsonData["source"]))
        else:
            if jsonData["class"] in regulatory_features:
                bedlines["regulatory"].append((jsonData["chr"], jsonData["start"], jsonData["end"],
                                           jsonData["class"], jsonData["regulatory_ID"], jsonData["regulatory_ID"]))

    # If GENOCE features are selected in the set, creating dataframe:
    full_df = pd.DataFrame()
    Y_top = 10
    if len(bedlines["GENCODE"]) > 0:
        a = BedTool(bedlines["GENCODE"]).sort()
        GENCODE_df = a.merge().to_dataframe()

        # Adding extra annotations:
        GENCODE_df["name"] = gene_name
        GENCODE_df["feature_type"] = "exon"
        GENCODE_df.loc[-1] = [chromosome, GENCODE_df.start.min(), GENCODE_df.end.max(), gene_name, "gene"]
        GENCODE_df["ID"] = gene_ID
        GENCODE_df['biotype'] = "protein_coding"

        # Position and the color is fixed:
        GENCODE_df['y_position'] = Y_top
        GENCODE_df['color'] = "LightCoral"

        # Adding dataframe to the full df:
        full_df = GENCODE_df

        # Adjusting Y_top
        Y_top = 9

    # If regulatory features are selected in the set, create dataframe:
    if(len(bedlines['regulatory'])>0):
        a = BedTool(bedlines["regulatory"]).sort()
        regulatory_df = a.to_dataframe().drop_duplicates()

        # Fixing header, adding columns:
        regulatory_df.columns = ['chrom', 'start', 'end', 'biotype', 'ID', 'name']
        regulatory_df["feature_type"] = "regulatory"
        regulatory_df["color"] = [get_biotype_color(x) for x in regulatory_df["biotype"]]

        # Assigning ranks and y positions based on the rank:
        rank_dict = get_gene_rank(regulatory_df[["start", "end", 'name']])
        regulatory_df['y_position'] = [Y_top - 0.5 * rank_dict[x] for x in regulatory_df["name"]]

        # Merging dataframe if exists:
        try:
            full_df = full_df.append(regulatory_df)
        except:
            full_df = regulatory_df
    return(full_df)

def get_biotype_color(biotype):
    '''
    Based on a list of biotypes, this function returns a list colors
    If a biotype is not found, gray color will be assigned
    '''

    # Coloring features based on their biotype:
    biotypeColor = {

        # Gene biotypes:
        'protein_coding' : "lightsalmon",

        'antisense' : '#525252',
        'processed_pseudogene' : '#5E5E5E',
        'unprocessed_pseudogene' : '#696969',
        'transcribed_processed_pseudogene' : "#757575",
        'ribozyme' : "#7F7F7F",
        'sense_intronic' : '#8A8A8A',
        'transcribed_unprocessed_pseudogene' : '#969696',
        'miRNA' : '#A3A3A3',
        'rRNA' : '#ABABAB',
        'snRNA' : '#B5B5B5',
        'misc_RNA' : '#BEBEBE',

        # Regulatory feature type:
        'CTCF binding site' : "#78AB46",
        'Open chromatin region' : "#61B329",
        'Transcription factor binding site' : "#4DBD33",
        'TF_binding_site' : "#4DBD33",
        'Predicted enhancer region' : "#7BCC70",
        'enhancer' : "#3F9E4D",
        'Predicted promoter flanking region' : "#3D9140",
        'Predicted promoter' : '#548B54',
        'promoter' : '#548B54'
    }

    try:
        return biotypeColor[biotype]
    except:
        return "#CDCDCD"

def draw_box(z, df, height = 0.4):
    '''
    Draws a box based on the dataframe
    required fields: y_position, start, end, color
    '''
    df["top"] = df["y_position"] + height/2
    df['bottom'] = df["y_position"] - height/2

    x = z.quad(bottom = 'bottom',
           top = 'top',
           left='start',
           right='end',
           color='color',
          source = df)
    return(x)

def draw_line(z, df, lw = 2):
    '''
    Draws a box based on the dataframe
    required fields: y_position, start, end, color
    '''
    x = z.segment(x0='start',x1='end',y0="y_position",y1="y_position",
              line_width=lw,color="color", source=df)
    return(x)
