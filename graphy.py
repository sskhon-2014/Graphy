import pandas as pd
import itertools
import sys
import os
import spectra
from io import StringIO
from ipywidgets import interact
import ipywidgets as widgets
from IPython.display import display
from ipywidgets import FloatSlider
from IPython.display import clear_output
import pip 
pip.main(['install', 'pysam'])

import pip 
pip.main(['install', 's3fs'])
clear_output()

import clip_tools as ct
from functions import *

updates = "\
3/23/17: Now the tracks are displayed in the order they were selected! <p> \
4/11/17: Added custom section to bedfiles. Covers custom bedfiles, including custom miRNA tracks <p> \
4/27/17: Updated LeftToRight button\
"


###################################################


file_utr3 = "s3:graphy101/genomes/refseq_3utr.bed"
file_refseq = "s3:graphy101/genomes/refseq_names.txt"
file_targetscan = "s3:graphy101/genomes/Targetscan7_fixed.bed"
file_refgenebed = 's3:graphy101/genomes/refGene.bed'
default_track_folder = "s3:graphy101/GRAPHY_TEST/"
# Placeholder defaults
placeholder_file_bedfile = "s3:graphy101/RK53/mmu-miR-29a-3p.bed"
###############################

# Placeholder defaults
placeholder_file_bedfile = "s3:graphy101/RK53/mmu-miR-29a-3p.bed"

bedfiles_allowed = {"Peaks":{"name":"Peaks",
                          "bedfile":"s3:graphy101/FASTUTR/peaks.bed",
                          "bedtype":"bed"},
               "Targetscan":{"name":"Targetscan",
                          "bedfile":file_targetscan,
                          "bedtype":"targetscan"}}



### Lookup 3'UTR by Refseq
df_refseq_3utr = pd.read_table("s3:graphy101/genomes/mm10_refseq_3utr.bed",names=['chrom','start','stop','name','score','strand'])
df_refseq_3utr.name = ["_".join(x.split("_")[0:2]) for x in df_refseq_3utr.name]
df_refseq_3utr = df_refseq_3utr.drop_duplicates("name")
df_refseq_3utr = df_refseq_3utr.set_index("name")

###################################################
from IPython.core.display import Image, display
display(Image('Data/logo.png', width=300, unconfined=False))
###################################################


# Make a genelookup object.
bd = ct.BetweenDict(file_refgenebed)

##############################################

print_plotting_parameters = True

##############################################


legend = True # adds legend box to alignment bam track.
staticaxes = False # keeps axes the same across all bam files. (for normalized data)
LeftToRight = False # Makes all genes read left to right. 
axis_off = False  # Don't show the axis labels.

### FILE OUTPUT OPTIONS ###
#dpi = 300
# figwidth = "scale"  # Can be "scale" or <int>  (15 is a good basic value)
# outputformat = "pdf"
# outputsuffix = ""

# Mostly Untouched Default Variables.



# Set static parameters
limits = "default" # can also be [int_start,int_stop]

##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################



defaultloc = 'chr5:142,903,034-142,906,867'
defaultgeneid = 'Actb'
defaultrefseqid = 'NM_007393'
defaultstrand = "-"

#### Accordian Change File Locations ####
w_default_folder = widgets.Text(
    description='RnaSeq Folder',
    value=default_track_folder,
)
w_default_folder.width="80%"
w_default_folder_valid = widgets.Valid(value=True)
hboxdefaultfolder = widgets.HBox([w_default_folder,w_default_folder_valid])
page1 = widgets.VBox(children=[hboxdefaultfolder])

accord = widgets.Accordion(children=[page1], width="80%")
display(accord)

accord.set_title(1, 'Defaults')


#### Define Widgets ###
lookuptype_widget = widgets.RadioButtons(
    options={'By location':'location','By geneid (3\'UTR only)':"geneid"},
    value='location',
    description='Lookup:',
    disabled=False
)

location_widget = widgets.Text(
    value=defaultloc,
    placeholder='Type something',
    description='Location:',
    disabled=False
)


strand_widget = widgets.RadioButtons(
    options=['+','-'],
    value=defaultstrand,
    description='Strand:',
    disabled=False
)


geneid_widget = widgets.Dropdown(
        options=[defaultgeneid],
        value=defaultgeneid,
        description='GeneID:',
        disabled=False,
        button_style='' # 'success', 'info', 'warning', 'danger' or ''
    )

refseqid_widget = widgets.Dropdown(
        options=[defaultrefseqid],
        value=defaultrefseqid,
        description='RefseqID:',
        disabled=False,
        button_style='' # 'success', 'info', 'warning', 'danger' or ''
    )


## Updating Widgets ##

def update_by_type(*args):
    if lookuptype_widget.value == "location":
        location_widget.value = defaultloc
    if lookuptype_widget.value == "geneid":
        location_widget.value = defaultgeneid

def update_geneid(*args):
    if lookuptype_widget.value == "location":
        chrom,start,stop = ct.easylocation(location_widget.value)
        strand = strand_widget.value
        genelist =  list(get_gene_id(chrom,start,stop,strand,bd,choice='all'))
        if len(genelist)>0:
            geneid_widget.options = genelist
        else:
            geneid_widget.options = ["None"]
    elif lookuptype_widget.value == "geneid":
        geneid_widget.options = [location_widget.value]
        geneid_widget.value = location_widget.value
        

def update_refseqid(*args):
    reflist =  get_refseq_id(file_refseq,geneid_widget.value,choice="all")
    print(reflist)
    if len(reflist)>0:
        refseqid_widget.options = reflist
    else:
        refseqid_widget.options = ["None"]

lookuptype_widget.observe(update_by_type,'value')
location_widget.observe(update_geneid, 'value')
strand_widget.observe(update_geneid, 'value')
geneid_widget.observe(update_refseqid,'value')

display(widgets.HTML(value="<h5>Updates:</h5><p>" + updates))
display(widgets.HTML(value="<h3>Select Region</h3>"))

### Define INTERACT: All interdependent widgets must be in here ###
def printer(lookuptype,loc, strand, geneid,refseqid):
    print("")
interact(printer,
         lookuptype = lookuptype_widget,
         loc=location_widget,
         strand=strand_widget,
         geneid=geneid_widget,
         refseqid = refseqid_widget,
         continuous_update=False)


#Lists the contents of 'GRAPHY_TEST' into list_dir_default_track_folder
import boto3
s3 = boto3.resource('s3')
my_bucket = s3.Bucket('graphy101')
list_dir_default_track_folder = []
for object in my_bucket.objects.all():
    if str(object.key)[0:11] == 'GRAPHY_TEST':
            if str(object.key)[12:] != '':
                list_dir_default_track_folder.append(str(object.key)[12:])

### Select File ###
all_files = list_dir_default_track_folder
file_types = ["bam","bw"]
file_names = []
for f in all_files:
    if f.split(".")[-1] in file_types:
        file_names.append(f)
file_widget = widgets.SelectMultiple(
    options=file_names,
    disabled=False
)

def update_antisense_check(*args):
    if len(file_widget.value)==0:
        w_as_selectors.layout.display="none"
        w_antisense_check.options = []
    if len(file_widget.value)>0:
        w_as_selectors.layout.display=""
        # Add or remove from antisence_check.options
        options = list(w_antisense_check.options)
        for i in file_widget.value:
            if not i in options:
                options.append(i)
        for i in options:
            if not i in file_widget.value:
                options.remove(i)
        w_antisense_check.options = options



file_widget.observe(update_antisense_check,'value')

w_aschecktitle= widgets.HTML("Are any of these tracks antisense?")
# The antisense_check box is also holding the order of the files selected
w_antisense_check = widgets.SelectMultiple(
    options=[]
)
w_as_selectors = widgets.VBox([w_aschecktitle,w_antisense_check])
w_as_selectors.layout.display="none"

### Style Selectors ###
buttonwidth="15%"
buttonstyle='info'
w_labels = widgets.ToggleButton(
    value=False,
    description='BedLabels',
    button_style=buttonstyle, # 'success', 'info', 'warning', 'danger' or ''
    tooltip='Label each element of a bedtrack',
    width=buttonwidth,
)
w_invert = widgets.ToggleButton(
    value=False,
    description='Left to Right',
    button_style=buttonstyle, # 'success', 'info', 'warning', 'danger' or ''
    tooltip="Always display genes from left to right regardless of strand.",
    width=buttonwidth,

)
w_legend = widgets.ToggleButton(
    value=True,
    description='Legend',
    button_style=buttonstyle, # 'success', 'info', 'warning', 'danger' or ''
    tooltip="Show Legend for Seq-Tracks",
    width=buttonwidth,

)
w_stagger = widgets.ToggleButton(
    value=False,
    description='Stagger Bedtracks',
    button_style=buttonstyle, # 'success', 'info', 'warning', 'danger' or ''
    tooltip="Stagger Bedtracks",
    width=buttonwidth,
)
w_refseq = widgets.ToggleButton(
    value=True,
    description='Show Refseq',
    button_style=buttonstyle, # 'success', 'info', 'warning', 'danger' or ''
    tooltip="Refseq track?",
    width=buttonwidth,
)
w_shadeBed = widgets.ToggleButton(
    value=False,
    description='Shade Bed',
    button_style=buttonstyle, # 'success', 'info', 'warning', 'danger' or ''
    tooltip="Highlight bedtrack regions",
    width=buttonwidth,
)
w_staticaxes = widgets.ToggleButton(
    value=True,
    description='Autoscale',
    button_style=buttonstyle, # 'success', 'info', 'warning', 'danger' or ''
    tooltip="Scale the rnaseq tracks?",
    width=buttonwidth,
)
w_boundingbox = widgets.ToggleButton(
    value=True,
    description='Show BoundingBox',
    button_style=buttonstyle, # 'success', 'info', 'warning', 'danger' or ''
    tooltip="Show axes and border box?",
    width=buttonwidth,
)
formattingwidth ="20%"
w_xscale = widgets.FloatText(
    description="x scale (inches)",
    placeholder="Int",
    tooltip="can enter 0 to autoscale",
    width=formattingwidth,
)
w_fontsize = widgets.FloatText(
    description="Font Size",
    placeholder="Int",
    tooltip="can enter 0 to autoscale",
    value=12,
    width=formattingwidth,
)
w_dpi = widgets.FloatText(
    description="DPI",
    placeholder="300",
    value="300",
    width=formattingwidth,
)
w_fileformat = widgets.Text(
    description="File Format",
    placeholder="Pdf,Svg,PNG",
    value="pdf",
    width=formattingwidth,
)
w_filesuffix = widgets.Text(
    description="File Suffix",
    placeholder="None",
    value="",
    width=formattingwidth,
)

w_outputfolder = widgets.Text(
    value="Output/",
    description="Output Folder",
    width="80%",
)
w_outputfolder_valid = widgets.Valid(value=True)
hboxoutputfolder = widgets.HBox([w_outputfolder,w_outputfolder_valid])
def output_folder_check(*args):
    if os.path.isdir(w_outputfolder.value) and w_outputfolder.value[-1] == "/" :
        w_outputfolder_valid.value=True
    else:
        w_outputfolder_valid.value=False
w_outputfolder.observe(output_folder_check,'value')


### Color Selectors ###
w_cp1=widgets.ColorPicker(value="#0080ff",width="100px",concise=False)
w_cp2=widgets.ColorPicker(value="#ff0000",width="100px",concise=False)
w_cp3=widgets.ColorPicker(value="#00ff00",width="100px",concise=False)
w_cp4=widgets.ColorPicker(value="#800080",width="100px",concise=False)



selectors_color = widgets.HBox([w_cp1,w_cp2,w_cp3,w_cp4])
selectors1 = widgets.HBox([w_labels,
                           w_invert,
                           w_legend,
                           w_stagger,
                           w_refseq,
                           w_shadeBed,
                           w_staticaxes,
                           w_boundingbox,
                          ])
selectors2 = widgets.HBox([w_xscale,
                           w_fontsize,
                           w_dpi,
                           w_fileformat,
                           w_filesuffix,
                          ])


### MiRNA and Bedtrack Picker ###
bedtracks_widget = widgets.RadioButtons(
    options=["None"] + list(bedfiles_allowed.keys()) + ["Custom","Other"],
    description='Bedfiles:',
    disabled=False
)
other_bedtracks_widget = widgets.Text(
    value=placeholder_file_bedfile,
    description="Other:",
    width="90%"
)
miRNA_widget = widgets.Select(
    options=["None"],
    value="None",
    description = "miRNA Family",
)


def update_miRNAlist(*args):
    miRNAwrapper.layout.display="none"
    otherbedtrackwrapper.layout.display="none"
    def get_targetscan_families(targetscan_file):
        beddb = pd.read_table(targetscan_file,names=["chrom","start","stop","miRNA"],usecols=[0,1,2,3])
        families = beddb.miRNA.unique()
        f = families.tolist()
        return sorted(f[1:])
    if bedtracks_widget.value=="Targetscan":
        miRNA_widget.description = "miRNA Family"
        miRNAwrapper.layout.display=""
        tsfamilies = get_targetscan_families(bedfiles_allowed[bedtracks_widget.value]["bedfile"])
        miRNA_widget.options = ["ALL"] + tsfamilies
    if bedtracks_widget.value=="Other":
        otherbedtrackwrapper.layout.display=""
    if bedtracks_widget.value in ["None","Peaks"]:
        miRNAwrapper.layout.display="none"
        miRNA_widget.options = ["None"]
        miRNA_widget.value = "None"
    if bedtracks_widget.value=="Custom":
        miRNAwrapper.layout.display=""
        if os.path.isdir(w_default_folder.value) and w_default_folder.value[-1] == "/" :
            w_default_folder_valid.value=True
            track_folder = w_default_folder.value
            all_files = os.listdir(track_folder)
            file_names = []
            file_types = ["bed"]
            for f in all_files:
                if f.split(".")[-1] in file_types:
                    file_names.append(f)
            miRNA_widget.options = file_names
            miRNA_widget.description = "Bedfile"
        else:
            w_default_folder_valid.value=False
            file_names = []
            miRNA_widget.options =[]
# def update_bedfile_other(*args):
#     if bedtracks_widget.value == "Other":


bedtracks_widget.observe(update_miRNAlist,'value')
#otherbedtrackwrapper.observe(update_bedfile_other,'value')

display(widgets.HTML(value="<h3>Files</h3>"))
display(bedtracks_widget)
otherbedtrackwrapper = widgets.VBox([other_bedtracks_widget])
display(otherbedtrackwrapper)
otherbedtrackwrapper.layout.display="none"
# Puts miRNA_widget in a VBox that can have its display toggled. 
miRNAwrapper = widgets.VBox([miRNA_widget])
display(miRNAwrapper)
miRNAwrapper.layout.display="none"

def miRNAUpdater(bedfile,miRNApick):
    print(bedfile)

miRNA_widget.layout.display=""
# interact(miRNAUpdater,
#          bedfile=bedtracks_widget,
#          miRNApick = miRNA_widget,
#          continuous_update=False,
#          font_size=4)
miRNA_widget.layout.width="100%"
file_widget.layout.width="100%"
def update_filelist(*args):
    if os.path.isdir(w_default_folder.value) and w_default_folder.value[-1] == "/" :
        w_default_folder_valid.value=True
        print(w_default_folder.value)
        track_folder = w_default_folder.value
        all_files = os.listdir(track_folder)
        file_names = []
        file_types = ["bam","bw"]
        for f in all_files:
            if f.split(".")[-1] in file_types:
                file_names.append(f)
        file_widget.options = file_names
    else:
        w_default_folder_valid.value=False
        file_names = []
        file_widget.options =[]
w_default_folder.observe(update_filelist,'value')

# Display #
display(widgets.HTML(value="Alignments to Use"))
display(file_widget)
display(w_as_selectors)
display(widgets.HTML(value="<h3>Formatting</h3>"))
display(selectors1)
display(widgets.HTML(value="<h4>Color Picker</h4>"))
display(selectors_color)
display(widgets.HTML(value="<h4>Output Format</h4>"))
display(selectors2)
display(hboxoutputfolder)
display(widgets.HTML(value="<br><br>"))

accord.selected_index=0 #starts the accordian closed.


###### LAUNCH BUTTON #####
button = widgets.Button(description="Graph!")
display(button)
def on_button_clicked(b):
    x = 1#clear_output()
    
    ## SANITY CHECKS ###
    if not len(file_widget.value)>0:
        print("ERROR: It's required to pick a .bam alignment!")
        return
    #####################



    if lookuptype_widget.value=="location":
        chrom,start,stop = ct.easylocation(location_widget.value)
        strand = strand_widget.value
    elif lookuptype_widget.value=="geneid":
        chrom,start,stop,value,strand = df_refseq_3utr.loc[refseqid_widget.value]
    if bedtracks_widget.value == "None":
        bedtrack=False
        bedfile=None
        bedtype=None
        name=None
    elif bedtracks_widget.value == "Peaks":
        bedtrack=True
        bedfile = bedfiles_allowed[bedtracks_widget.value]["bedfile"]
        bedtype = "bed"
        name = "Peaks"
    elif bedtracks_widget.value == "Other":
        bedtrack=True
        bedfile = other_bedtracks_widget.value
        bedtype = "bed"
        name = "Bedtrack"       
    elif bedtracks_widget.value == "Targetscan":
        bedtrack=True
        bedfile = bedfiles_allowed[bedtracks_widget.value]["bedfile"]
        if miRNA_widget.value == "ALL":
            bedtype = bedfiles_allowed[bedtracks_widget.value]["bedtype"]
            name = bedfiles_allowed[bedtracks_widget.value]["name"]
        else:
            bedtype = "targetscan"
            name = miRNA_widget.value
    elif bedtracks_widget.value == "Custom":
        bedtrack=True
        bedfile = w_default_folder.value+miRNA_widget.value
        bedtype = "bed"
        name = miRNA_widget.value
    track_type=[]
    for t in w_antisense_check.options:
        if t in w_antisense_check.value:
            track_type.append("as")
        else:
            track_type.append("s")
    geneid= geneid_widget.value
    refseqid = refseqid_widget.value
    track_names = [f for f in w_antisense_check.options if f.split(".")[-1]=="bam"]
    bigwignames = [f for f in w_antisense_check.options if f.split(".")[-1]=="bw"]
    track_files = [w_default_folder.value + f for f in track_names]
    bigwigfiles = [w_default_folder.value + f for f in bigwignames]
    depths = get_depth_data(track_files,track_names,chrom,start,stop,strand,track_type)
    wig_df_list = get_wig_data(bigwigfiles,bigwignames,chrom,start,stop)
    
    outputformat = w_fileformat.value
    figwidth = w_xscale.value
    dpi = w_dpi.value
    LeftToRight = w_invert.value
    annotate_bed = w_labels.value
    staggerbed = w_stagger.value
    refseqtrack = w_refseq.value
    shade_by_bed = w_shadeBed.value
    outputsuffix = w_filesuffix.value
    staticaxes = w_staticaxes.value
    axis_off= not w_boundingbox.value
    fontsize = w_fontsize.value
    legend = w_legend.value
    figheight = 2.5*len(w_antisense_check.options)
    output_folder = w_outputfolder.value

    color_values = [w_cp1.value,w_cp2.value,w_cp3.value]
    shade = itertools.cycle(color_values)
    colors = itertools.cycle([spectra.html(i).darken(20).hexcode for i in color_values])
    plot(figwidth,figheight,refseqtrack,LeftToRight,strand,depths,
       colors,shade,limits,bedtrack,start,stop,staggerbed,bigwignames,
        wig_df_list,shade_by_bed,output_folder,geneid,outputsuffix,outputformat,dpi,track_names,axis_off,
       legend,staticaxes,bedfile,bedtype,name,chrom,refseqid,annotate_bed,fontsize)

button.on_click(on_button_clicked)


