#!/usr/bin/python

# -*- coding: utf-8 -*-

# Module metadata variables
__author__ = "Diego Mena Santos"
__credits__ = ["Diego Mena Santos", "Jose Rodriguez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "0.1.0"
__maintainer__ = "Diego Mena Santos"
__email__ = "diego.mena@cnic.es;jmrodriguezc@cnic.es"
__status__ = "Development"


### Import modules

import pandas as pd
import numpy as np
import sys
import os
import configparser
import warnings
import requests
import sys
import argparse
import logging
warnings.filterwarnings("ignore")


def calculate_starts(df):

    """
    Calculates the start position of the feature in the genome
    """

    empieza=[]
    sobra=[]
    for i in range(len(df)):
        if i==0:
            inicio=0
            final=((df['gle'][i]+1)-(df['glb'][i]))%3
            empieza.append(0)
            sobra.append(final)
        else:
            inicio=3-final
            if inicio==3:
                empieza.append(0)
            else:
                empieza.append(inicio)
            final=((df['gle'][i]+1)-(df['glb'][i]+inicio))%3
            sobra.append(final)

    df['start']=empieza

    return df


def prot_location(id):

    """
    Recovers exons genomic positions for a protein from Uniprot 
    """
    url=f"https://www.ebi.ac.uk/proteins/api/coordinates/{id}"
    r = requests.get(url, headers={ "Accept" : "application/json"})

    gl=r.json()
    try:
        result=gl['gnCoordinate'][0]['genomicLocation']
        df_1=pd.DataFrame()
        reverse=False
        if gl['gnCoordinate'][0]['genomicLocation']['reverseStrand']==True:
            reverse=True
            strand_d='-'
        else:
            strand_d='+'

        exon=1

        for i in result['exon']:
            if reverse==False:
                a=pd.Series([i['proteinLocation']['begin']['position'],i['proteinLocation']['end']['position'],i['genomeLocation']['begin']['position'],i['genomeLocation']['end']['position']])
                df_1[f'{exon}']=a
                exon+=1     
            else:
                a=pd.Series([i['proteinLocation']['begin']['position'],i['proteinLocation']['end']['position'],i['genomeLocation']['end']['position'],i['genomeLocation']['begin']['position']])
                df_1[f'{exon}']=a
                exon+=1

        df_1=df_1.transpose()
        df_1.columns=['b','e','glb','gle']

        df_def=calculate_starts(df_1)

        if result['chromosome']=='X':
            chr='X'
        elif result['chromosome']=='Y':
            chr='Y'
        else:
            chr=result['chromosome']

        return df_def,strand_d, gl['sequence'],chr
    
    except KeyError:
        pass


def read_file(path,q_column):

    """
    Read infile
    """
    
    df=pd.read_csv(path,sep='\t',low_memory=False)
    proteins=df[q_column].unique()

    return proteins,df


def establish_feature(row,criteria,chr,strand,attributes,g_col_head,
                      nms_name,score_NM,score_pgm,missing_cleavages_header,
                      dt_name,qf_name,start_pos_prot_head,end_pos_prot_head,
                      dp_name,score_qfs,source_comparison):

    """
    Defines the score, feature, source, chromosome, frame and attributes columns
    """

    if criteria==True:
        row['feature']=row[g_col_head]
        if row['feature']==nms_name:
            if type(row[score_NM])!=np.nan:
                row['score']=row[score_NM]
            else:
                row['score']='.'
        else:
            if type(row[score_pgm])!=np.nan:
                row['score']=row[score_pgm]
            else:
                row['score']='.'
    else:
        if row[missing_cleavages_header]==dt_name:
            row['feature']=qf_name
            positions=row[qf_name].split(':')[1].split('_')
            row[start_pos_prot_head]=int(positions[0])
            try:
                row[end_pos_prot_head]=int(positions[1])
            except ValueError:
                 row[end_pos_prot_head]=int(row[start_pos_prot_head])
        elif row[missing_cleavages_header]==dp_name:
            row['feature']=dp_name
        if type(row[score_qfs])!=np.nan:
            row['score']=row[score_qfs]
        else:
            row['score']='.'

    row['chr']=chr
    row['strand']=strand
    row['frame']='.'
    row['source']=source_comparison

    atts=row.loc[attributes].to_dict()

    atts_str=''

    for i in atts.keys():
            atts_str+=str(i)+'='+str(atts[i])+';'

    row['attributes']=atts_str[:-1]


    return row


def locateinGenome(row,strand,coordinates,mod_name,mod_pos_prot_head,start_pos_prot_head,
                   end_pos_prot_head):
    
    """
    Locate the features in the genome
    """

    if strand=='-':
        factor=-1
        coordinates=coordinates.rename(columns={'gle':'glb','glb':'gle'})
    else:
        factor=1

    if row['feature']==mod_name:

        prot_coordinates=coordinates.loc[(coordinates['b']<=row[mod_pos_prot_head])&(coordinates['e']>=row[mod_pos_prot_head])].reset_index()
        if len(prot_coordinates)==1:
            if prot_coordinates['start'][0]==0:
                rest=0
            elif prot_coordinates['start'][0]==1:
                rest=1
            else:
                rest=-1
            row['glb']=(prot_coordinates['glb'][0]-prot_coordinates['start'][0]*factor-rest*factor)+(3*(row['n']-prot_coordinates['b'][0]))*factor
            row['gle']=row['glb']+2*factor
        
        elif len(prot_coordinates)==2:
            if prot_coordinates.iloc[-1]['start']==1:
                row['glb']=prot_coordinates.iloc[0]['gle']-1*factor
                row['gle']=prot_coordinates.iloc[-1]['glb']
            else:
                row['glb']=prot_coordinates.iloc[0]['gle']
                row['gle']=prot_coordinates.iloc[-1]['glb']+1*factor

    else:
        qf_start=int(row[start_pos_prot_head])
        qf_end=int(row[end_pos_prot_head])
        a=coordinates.loc[(coordinates['b'].isin(range(qf_start,qf_end)))|(coordinates['e'].isin(range(qf_start,qf_end)))].reset_index(drop=True)
        if len(a)==0:
            a=coordinates.loc[(coordinates['b']<=qf_start)&(coordinates['e']>=qf_end)].reset_index(drop=True)
            if a['start'][0]==0:
                rest=0
            elif a['start'][0]==1:
                rest=1
            else:
                rest=-1
            row['glb']=(a['glb'][0]-a['start'][0]*factor-rest*factor)+(3*(row[start_pos_prot_head]-a['b'][0]))*factor
            row['gle']=int((a['glb'][0]-a['start'][0]*factor-rest*factor)+(3*(row[end_pos_prot_head]-a['b'][0]))*factor+2*factor)

        
        else:
            if a['start'][0]==0:
                rest=0
            elif a['start'][0]==1:
                rest=1
            else:
                rest=-1
            row['glb']=(a['glb'][0]-a['start'][0]*factor-rest*factor)+(3*(row[start_pos_prot_head]-a['b'][0]))*factor

            if a.iloc[-1]['start']==0:
                rest=0
            elif a.iloc[-1]['start']==1:
                rest=1
            else:
                rest=-1
            row['gle']=(a.iloc[-1]['glb']-a.iloc[-1]['start']*factor-rest*factor)+(3*(row[end_pos_prot_head]-a.iloc[-1]['b'])*factor)+2*factor

    return row


def create_gff(proteins,df,attributes,g_col_head,nms_name,score_NM,score_pgm,
               missing_cleavages_header,dt_name,qf_name,start_pos_prot_head,end_pos_prot_head,
               dp_name,score_qfs,source_comparison,mod_name,mod_pos_prot_head):
    
    """
    Create the gff file
    """

    gffs=[]

    for i in proteins:

        try:
            coordinates, strand,sequence,chr=prot_location(i)
            df_filt=df[df.q.eq(i)].reset_index(drop=True)
            df_filt_qfs=df_filt.drop_duplicates(['qf','Missing_Cleavages'])
            df_filt[['glb','gle','feature','frame','strand','attributes','chr','score']]=''

            try:
                df_filt=df_filt.apply(lambda y: establish_feature(y,True,chr,strand,attributes,g_col_head,nms_name,score_NM,score_pgm,missing_cleavages_header,dt_name,qf_name,start_pos_prot_head,end_pos_prot_head,dp_name,score_qfs,source_comparison),axis=1)
                df_filt_qfs=df_filt_qfs.apply(lambda y: establish_feature(y,False,chr,strand,attributes,g_col_head,nms_name,score_NM,score_pgm,missing_cleavages_header,dt_name,qf_name,start_pos_prot_head,end_pos_prot_head,dp_name,score_qfs,source_comparison),axis=1)

            except KeyError as e:

                print('Check the headers in the input and the config, something does not match')
                print(e)
                sys.exit()

            
            df_prueba_mods=df_filt.apply(lambda y: locateinGenome(y,strand,coordinates,mod_name,mod_pos_prot_head,start_pos_prot_head,end_pos_prot_head),axis=1)
            df_prueba_qfs=df_filt_qfs.apply(lambda y: locateinGenome(y,strand,coordinates,mod_name,mod_pos_prot_head,start_pos_prot_head,end_pos_prot_head),axis=1)

            df_def=pd.concat([df_prueba_mods,df_prueba_qfs],axis=0).fillna('.')
            
            if strand=='-':
                df_def=df_def.rename(columns={'glb':'gle','gle':'glb'})
                # print('hola')

            df_def=df_def.loc[:,['chr','source','feature','glb','gle','score','strand','frame','attributes']]


            gffs.append(df_def)


        except TypeError:
            # print(i, "has no coordinates, maybe it's an inmunoglobulin")
            logging.info(i+" has no coordinates, maybe it's an inmunoglobulin")
            continue

        
        except KeyError:
            # print(i, 'coordinates do not matcht with yor data')
            logging.warning(i+ ' coordinates do not matcht with yor data')
            continue

    return gffs


def write_file(file,name):

    """
    Write the file
    """

    gff=pd.concat(file,axis=0)
    gff.to_csv(name,sep='\t',header=False,index=False)


def main(infile, q_column_header,attributes,
         g_col_head,nms_name,score_NM,score_pgm,missing_cleavages_header,
         dt_name,qf_name,start_pos_prot_head,end_pos_prot_head,dp_name,
         score_qfs,source_comparison,mod_name,mod_pos_prot_head,output):


    proteins,df=read_file(infile,q_column_header)

    
    logging.info('Locating features in genome')
    print()
    file=create_gff(proteins,df,attributes,g_col_head,nms_name,score_NM,score_pgm,missing_cleavages_header,dt_name,qf_name,start_pos_prot_head,end_pos_prot_head,dp_name,score_qfs,source_comparison,mod_name,mod_pos_prot_head)
    print()
    logging.info('Writting GFF file')
    write_file(file,output)



if __name__ == '__main__':
    
        # parse arguments
    parser = argparse.ArgumentParser(
        description='GFF creator',
        epilog='''
        Example:
            python gff_creator.py -i infile -o output -c config
        ''')
      
    # default PTMaps configuration file
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/Solver.ini")
    parser.add_argument('-i', '--infile', required=True, help='Path to input file')
    parser.add_argument('-o',  '--output', required=True, help='Name of the output file')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')
    args = parser.parse_args()

    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)

    q_column_header=config['Gff_creator']['q_column_header']
    missing_cleavages_header=config['Gff_creator']['missing_cleavages_header']
    start_pos_prot_head=config['Gff_creator']['start_position_in_protein_header']
    end_pos_prot_head=config['Gff_creator']['end_position_in_protein_header']
    mod_pos_prot_head=config['Gff_creator']['mod_position_in_protein_header']
    g_col_head=config['Gff_creator']['g_column_header']
    source_comparison=config['Gff_creator']['source_comparison']
    score_qfs=config['Gff_creator']['score_qfs']
    score_qfs=config['Gff_creator']['score_qfs']
    score_pgm=config['Gff_creator']['score_pgm']
    score_NM=config['Gff_creator']['score_NM']
    attributes=config['Gff_creator']['attributes'].split(',')
    nms_name=config['Gff_creator']['NM']
    mod_name=config['Gff_creator']['Mod']
    qf_name=config['Gff_creator']['qf']
    dp_name=config['Gff_creator']['dp']
    dt_name=config['Gff_creator']['dt']
    log_file = args.infile[:-4] + '_gff_creator_log.txt'


    logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    datefmt='%m/%d/%Y %I:%M:%S %p',
                    handlers=[logging.FileHandler(log_file),
                                logging.StreamHandler()])


    main(args.infile, q_column_header,attributes,
         g_col_head,nms_name,score_NM,score_pgm,missing_cleavages_header,
         dt_name,qf_name,start_pos_prot_head,end_pos_prot_head,dp_name,
         score_qfs,source_comparison,mod_name,mod_pos_prot_head,args.output)
    
    logging.info("end of the script")