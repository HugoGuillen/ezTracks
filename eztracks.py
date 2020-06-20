import sys
import os
from os import path,listdir
import pandas as pd
import csv
import subprocess
import argparse
import configparser
from collections import OrderedDict
import matplotlib
import matplotlib.pyplot as plt


#~~~~~~~~HELPERS FOR BUILDING PYGT CONFIG FILE
#~~~Gencode annotation
str_track_gtf="""[x-axis]
where=top

[test gtf]
title={gtf_title}
file={gtf_input}
file_type=gtf 
height={gtf_height}
#~I think this is a bug in gtf track implementation
color=black
#~Labels on annotations
labels=true
fontsize=8
prefered_name=transcript_name
#~Arrows on TSSs (comment line to remove)
style=tssarrow

[spacer]
height = 0.4

"""

str_track_line="""
[lspacer]
file_type=hlines
show_data_range=false
y_values=.5
overlay_previous=share-y
height = 0.8

#[spacer]
#height = 0.1

"""

#~~~Group spacer (more than 1 bed file per group)
str_track_grouplabel="""
[bedspacer {group}]
file = {emptybed}
title = {group}
file_type = bed
gene_rows = 1
height=.8
"""

#~~~Group spacer (General bed formatting)
str_track_bed="""
[test {group}-{title}]
file = {filename}
title = {title}
color = {color}
labels = false
file_type = bed
display = collapsed
gene_rows = 1
arrowhead_included = true
style = UCSC
height=.8

"""

#~~~~~~~~~~~~~~~~~~~~
def load_config(config_ini):
    if not path.exists(config_ini):
        raise FileNotFoundError("Config file not found. Exiting.")
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(config_ini)            
    return config

def load_tracks(config_obj):
    config=config_obj
    tracks = OrderedDict()
    for group in config['tracks']:
        P = config['tracks'][group]
        if not path.exists(P):
            raise FileNotFoundError('ERROR: track group',P,'not found.')
        if path.isfile(P):
            tracks[group] = [(P,group)]
        else:            
            files = [path.join(P,x) for x in listdir(P) if x.endswith('.bed') or x.endswith('.bed.gz')]
            if len(files)==0:
                print('# No bed tracks found in %s. Skipping.',file=sys.stderr)
                continue
            tracks[group] = []
            for p in files:                
                p_name = path.basename(p).replace('.bed.gz','').replace('.bed','')                
                tracks[group].append((p,p_name))                
    return tracks
    
def print_config(config_obj):
    for s in config_obj.sections():
        print('[%s]'%s)
        for k in config_obj[s]:
            print('%s = %s'%(k,config_obj[s][k]))

#~~~~~~~~~~~~~~~~~~~~
def check(config_ini):
    errors = 0
    #~~~Check programs
    try:
        call = 'echo $CONDA_PREFIX;python -V'
        call_output = subprocess.check_output(call,shell=True).decode().strip().replace('\n',': ')
        print('# conda_env:',call_output,file=sys.stderr)
    except:
        print('ERROR: Conda environment not found. Exiting.',file=sys.stderr)
        errors+=1
    try:
        call = 'pyGenomeTracks --version'
        call_output = subprocess.check_output(call,shell=True).decode().strip()
        print('#',call_output,file=sys.stderr)
    except:
        print('ERROR: pyGenomeTracks not found. Exiting.',file=sys.stderr)
        errors+=1
    try:
        call = 'bedtools -version'
        call_output = subprocess.check_output(call,shell=True).decode().strip()
        print('#',call_output,file=sys.stderr)
    except:
        print('ERROR: bedtools not found. Exiting.',file=sys.stderr)
        errors+=1
    
    #~~~Check configuration
    if path.exists(config_ini):
        print('# config_file found at',path.abspath(config_ini),file=sys.stderr)
    else:
        print('ERROR: config_file not found. Exiting.',file=sys.stderr)
        return
    config = load_config(config_ini)
    required_params = ['output_path','input_gtf']
    for param in required_params:
        if param not in config['default']:
            print('ERROR: required parameter %s not found in configuration.'%param,file=sys.stderr)
            errors+=1
    if 'region' not in config['default'] and 'transcript' not in config['default']:
        print('ERROR: you need to specify a region or transcript in the configuration.',file=sys.stderr)
        errors+=1
    required_files = ['input_gtf']
    for  _file in required_files:
        filepath = config['default'][_file]
        if not path.exists(filepath):
            print('ERROR: required %s (%s) not found.'%(_file,filepath),file=sys.stderr)
            errors+=1
    
    #~~~Check track files
    tracks = OrderedDict()
    for group in config['tracks']:
        P = config['tracks'][group]
        if not path.exists(P):
            print('ERROR: track group',P,'not found.',file=sys.stderr)
            errors+=1
            continue
        print('# Found %s: %s'%(group,P),file=sys.stderr)
        if path.isfile(P):
            tracks[group] = [(P,group)]
        else:
            tracks[group] = []
            files = [path.join(P,x) for x in listdir(P) if x.endswith('.bed') or x.endswith('.bed.gz')]
            for p in files:                
                p_name = path.basename(p).replace('.bed.gz','').replace('.bed','')                
                tracks[group].append((p,p_name))
                print('# Found %s-%s: %s'%(group,p_name,p),file=sys.stderr)    

    #~~~End checkup
    print('# Checkup finalized. %d error(s) detected.'%errors,file=sys.stderr)
    if errors>0:
        print('# Please solve all the configuration errors before calling "prepare" and "draw".',file=sys.stderr)    
    
#~~~~~~~~~~~~~~~~~~~~
def prepare(config_ini):
    calls = []
    call_gtf = 'intersectBed -sorted -u -a {input_gtf} -b {query_bed} | sort -k1,1 -k4,4n > {prep_gtf}'    
    call_intersect = 'intersectBed -sorted -u -a {input_bed} -b {query_bed} > {output_bed}'
    config = load_config(config_ini)
    output_path = path.abspath(config['default']['output_path'])
    input_gtf = config['default']['input_gtf']
    prep_path = path.join(output_path,'prep')
    prep_gtf = path.join(prep_path,'input.gtf')
    query_bed = path.join(prep_path,'query.bed')
    folders = [output_path,prep_path]
    process_transcript_relative = False
    for folder in folders:
        if not path.exists(folder):
            os.makedirs(folder)
            print('# Created dir %s'%folder,file=sys.stderr)
    if 'region' in config['default']:
        region = config['default']['region']
        print('# REGION MODE: %s'%region,file=sys.stderr)
        _chr = region.split(':')[0]
        start = str(int(region.split(':')[1].split('-')[0])-1)
        end = region.split(':')[1].split('-')[1]        
        with open(query_bed,'w') as f:            
            f.write('\t'.join([_chr,start,end])+'\n')
            print('# Wrote',query_bed,file=sys.stderr)    
        call = call_gtf.format(query_bed=query_bed,input_gtf=input_gtf,prep_gtf=prep_gtf)
        calls.append(call)
    elif 'transcript' in config['default']:
        transcript = config['default']['transcript']        
        relative_plot = config.getboolean('default','relative',fallback=False)
        if relative_plot:
            feature = 'exon'
            call_intersect = 'intersectBed -sorted -a {input_bed} -b {query_bed} > {output_bed}'
            process_transcript_relative = True
        else:
            feature= 'transcript'
        print('# TRANSCRIPT MODE (relative=%s): %s'%(relative_plot,transcript),file=sys.stderr)
        call_gtf_search = '{compressed}grep {transcript} {gtf} | sort -k1,1 -k4,4n > {prep_gtf}'
        call = call_gtf_search.format(transcript=transcript,gtf=input_gtf,prep_gtf=prep_gtf,
                                     compressed='z' if input_gtf.endswith('.gz') else '')
        calls.append(call)        
        call_bed = 'grep -P "\\t{feature}\\t" {prep_gtf} | awk \'BEGIN{{OFS=FS="\\t"}}{{print $1,$4-1,$5}}\' | sort -k1,1 -k2,2n > {query_bed}'
        call = call_bed.format(feature=feature,prep_gtf=prep_gtf,query_bed=query_bed)
        calls.append(call)
    if not process_transcript_relative:
        TRACKS = load_tracks(config)
        for group,tracks in TRACKS.items():
            track_path = path.join(prep_path,group)
            if not path.exists(track_path):
                os.makedirs(track_path)
                print('# Created dir %s'%track_path,file=sys.stderr)
            for track in tracks:
                output_bed = path.join(track_path,track[1]+'.bed')            
                call = call_intersect.format(query_bed=query_bed,input_bed=track[0],output_bed=output_bed)
                calls.append(call)    
        for call in calls:
            pass
            print(call,file=sys.stderr)        
            result = subprocess.run(call,shell=True,stdout=subprocess.PIPE)
            msg = result.stdout.decode().strip()
            if msg!='':
                print(msg,file=sys.stderr)    

        print('# Done.',file=sys.stderr)
        return
    #REMOVE, PUT EVERYTHING AFTER THE IF ELSE 
    for call in calls:        
        print(call,file=sys.stderr)
        result = subprocess.run(call,shell=True,stdout=subprocess.PIPE)
        msg = result.stdout.decode().strip()
        if msg!='':
            print(msg,file=sys.stderr) 
    #~~~~if process_transcript_relative    
    index_bed = path.join(prep_path,'index.bed')
    vert_bed = path.join(prep_path,'vertical.bed')
    #~~~Prepare index (chr,start,end,intron_length_to_subtract)
    df = pd.read_csv(query_bed,sep='\t',header=None)    
    X = [df.iloc[0][1]]+list(df[2].values)[:-1]
    df[3] = X
    df['length'] = df[2]-df[1]
    df['intron'] = df[1]-df[3]
    df['cumsum'] = df['intron'].cumsum()
    df[[0,1,2,'cumsum']].to_csv(index_bed,sep='\t',header=None,index=None)
    t_length = df['length'].sum()
    #~~~Prepare gtf
    dgt = pd.read_csv(prep_gtf,sep='\t',header=None)
    offset_gtf = dgt[3].min()
    dfexons = dgt[dgt[2]=='exon'].sort_values(by=3)
    dftrans = dgt[dgt[2]=='transcript'].copy()
    dfexons[3] = (dfexons[3]-df['cumsum'].values)-offset_gtf+1
    dfexons[4] = (dfexons[4]-df['cumsum'].values)-offset_gtf+1
    dfexons['bed'] = dfexons[3]-1
    dfexons[[0,'bed',3]].to_csv(vert_bed,sep='\t',index=None,header=None)
    with open(vert_bed,'a') as f:
        f.write('\t'.join([dfexons.iloc[0][0],str(t_length-1),str(t_length)])+'\n')
    dftrans[3] = 1
    dftrans[4] = t_length
    dgt = pd.concat((dftrans,dfexons)).sort_index()
    dgt.to_csv(prep_gtf,sep='\t',index=None,header=None,quoting=csv.QUOTE_NONE)    
    calls= []
    TRACKS = load_tracks(config)
    offset_bed = offset_gtf-1
    call_intersect = 'intersectBed -sorted -wb -a {input_bed} -b {query_bed} '
    call_intersect+= ' | awk \'BEGIN{{OFS=FS="\\t"}}{{printf $1FS$2-$NF-{offset}FS$3-$NF-{offset};for(i=4;i<=NF-4;i++) printf FS$i;print ""}}\' > {output_bed}'    
    for group,tracks in TRACKS.items():
        track_path = path.join(prep_path,group)
        if not path.exists(track_path):
            os.makedirs(track_path)
            print('# Created dir %s'%track_path,file=sys.stderr)
        for track in tracks:
            output_bed = path.join(track_path,track[1]+'.bed')            
            call = call_intersect.format(query_bed=index_bed,input_bed=track[0],output_bed=output_bed,offset=offset_bed)
            calls.append(call)
    for call in calls:        
        print(call,file=sys.stderr)
        result = subprocess.run(call,shell=True,stdout=subprocess.PIPE)
        msg = result.stdout.decode().strip()
        if msg!='':
            print(msg,file=sys.stderr) 
    #~~~
    


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
call_PGT = 'pyGenomeTracks --tracks {output_ini} --region {chr}:{view_start}-{view_end} --fontSize 12 --trackLabelFraction 0.2 --trackLabelHAlign left --width {width} --dpi {dpi} -o {output_img}'

def draw(config_ini):
    config = load_config(config_ini)
    #~~~
    relative_mode = False
    if 'transcript' in config['default'] and config.getboolean('default','relative',fallback=False):
        relative_mode = True
    #~~~
    output_path = path.abspath(config['default']['output_path'])
    prep_path = path.join(output_path,'prep')
    prep_gtf = path.join(prep_path,'input.gtf')
    output_ini = path.join(output_path,'config.ini')    
    output_img = path.join(output_path, config['plot'].get('output_img','my_tracks.pdf'))
    empty_bed = path.join(prep_path,'empty.bed')
    with open(empty_bed,'w') as f:
        pass
    
    #~~~Graphics parameters
    gtf_title = config['plot'].get('gtf_title','Annotation')
    padding = int(config['plot'].get('padding','100'))    
    width = config['plot'].get('width','40')
    dpi = config['plot'].get('dpi','100')
    gtf_height = config['plot'].get('gtf_height','10')
    
    #~~~Get plotting coordinates
    if 'region' in config['default']:
        region = config['default']['region']
        print('# REGION MODE: %s'%region,file=sys.stderr)
        _chr = region.split(':')[0]
        start = int(region.split(':')[1].split('-')[0])-1
        end = int(region.split(':')[1].split('-')[1])
        view_start = start-padding
        view_end = end+padding
    elif 'transcript' in config['default']:
        transcript = config['default']['transcript']        
        relative_plot = config.getboolean('default','relative',fallback=False)
        if not relative_plot or relative_plot:
            df = pd.read_csv(prep_gtf,sep='\t',header=None)
            _chr = df.iloc[0][0]
            start = df[3].min()-1
            end = df[4].max()
            view_start = max(start-padding,0)
            view_end = end+padding
    #~~~ Generate HTML color list
    cmap = plt.get_cmap('tab10')
    colors = [matplotlib.colors.rgb2hex(cmap(i)[:3]) for i in range(cmap.N)]
    
    #~~~Create ini file (GTF)
    with open(output_ini,'w') as f:
        f.write(str_track_gtf.format(gtf_input=prep_gtf,gtf_title=gtf_title,gtf_height=gtf_height))
        if relative_mode:            
            f.write("[vlines]\nfile = {file}\ntype = vlines\n".format(file=path.join(prep_path,'vertical.bed')))
    
    #~~~Add tracks to ini file
    TRACKS = load_tracks(config)
    track_count = 0
    with open(output_ini,'a') as f:
        for idx,(group,tracks) in enumerate(TRACKS.items()):
            color = colors[idx%10]
            track_path = path.join(prep_path,group)
            if len(tracks) == 1: #Only one track in the group
                filename = path.join(track_path,tracks[0][1]+'.bed')
                if os.stat(filename).st_size == 0:
                    continue
                f.write(str_track_bed.format(group=group,color=color,
                                        filename=filename,
                                        title=group))
                f.write(str_track_line)
                track_count+=1
            else:
                f.write(str_track_grouplabel.format(group=group,emptybed=empty_bed))
                for track in tracks:                    
                    filename = path.join(track_path,track[1]+'.bed')                    
                    if os.stat(filename).st_size == 0:
                        continue
                    f.write(str_track_bed.format(group=group,color=color,
                                                 filename=filename,
                                                 title=track[1]))
                    f.write(str_track_line)
                    track_count+=1                    
    call = call_PGT.format(output_ini=output_ini,chr=_chr,view_start=view_start,view_end=view_end,
                         output_img=output_img,width=width,dpi=dpi)    
    result = subprocess.run(call,shell=True,stderr=subprocess.PIPE)
    msg = result.stderr.decode().strip()
    if msg!='':
        print(msg,file=sys.stderr)
    print('# Wrote',output_img,file=sys.stderr)
    print('#######################################################################',file=sys.stderr)
    print('### If you want to customize the image modify %s and then rerun\n%s'%(path.abspath(output_ini),call),file=sys.stderr)
    print('#######################################################################',file=sys.stderr)

    
    
def main():
    parser = argparse.ArgumentParser(description='ezTracks. Plot a single GTF annotation followed by grouped bed files.\n\nHugo Guillen, 2020.')
    
    parser.add_argument('action',choices=['check','prepare','draw'],help="""check: check if config file, environment, and input files are correct.
prepare: preprocess input files into output defined in config file.
draw: generate output image into output defined in config file.""")
    parser.add_argument('configfile',help="path to configuration file. For possible options, please refer to the example located at test_data/config.ini")
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)    
    args = parser.parse_args()    
    if args.action == 'check': 
        check(args.configfile)
    elif args.action == 'prepare':
        prepare(args.configfile)
    elif args.action == 'draw':
        draw(args.configfile)
    else:
        parser.print_help()
        
if __name__ == "__main__":
    main()