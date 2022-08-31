import pandas as pd
import string as s
import random
from random import randint
import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO
from fuzzysearch import find_near_matches
import re
import os
import argparse

def RBD_RBM_matrix(Input_folder, Output_path, Data_matrix_folder):
    # dictionary
    c_index=["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
    c_number=list(range(0,20))
    r_positions=list(range(331,532))
    r_number=list(range(0,201))
    col=dict(zip(c_index,c_number))
    row=dict(zip(r_number,r_positions))
    # read binding score matrix
    wuhan_m=np.load(os.path.join( Data_matrix_folder,"wuhan_m.npy"))
    # read expression score matrix
    wuhan_expm=np.load(os.path.join( Data_matrix_folder,"wuhan_expm.npy"))
    
    # files
    file_list=os.listdir(Input_folder)
    # create df
    df = pd.DataFrame(columns=['file_name', 'S_sequence',"RBD_start","RBD_end","RBD_seq","Len_RBD","RBD_score_sum","RBD_score","RBM_seq","RBM_StartatRBD","RBM_EndatRBD","Len_RBM","RBM_score_sum","RBM_score","RBD_ex_sum","RBD_ex_score","RBM_ex_sum","RBM_ex_score"])
    
    for i in range(len(file_list)):
        if file_list[i].endswith("fasta"):
            path= os.path.join(Input_folder,file_list[i])
            record = SeqIO.read(path, "fasta")
            target_start=re.search(r"ATGTTTGTTTTT", str(record.seq))
            target_end=re.search(r"CATTACACATAA", str(record.seq))
            if target_start is None or target_end is None:
                print(file_list[i],":mismatch s gnen")
                pass
            else:
                s_gene=record[target_start.start():target_end.end()]
                s_protein=s_gene.translate()
                # search RBD
                bs=find_near_matches('NITNL', str(s_protein.seq), max_l_dist=1)
                be=find_near_matches('GPKKST', str(s_protein.seq), max_l_dist=1)
                #search RBM
                xs=find_near_matches('GCVIAWN', str(s_protein.seq), max_l_dist=1)
                xe=find_near_matches('PYRVVV', str(s_protein.seq), max_l_dist=1)
                
                if len(bs)!=1 or len(be)!=1 or len(xs)!=1 or len(xe)!=1:
                    print(file_list[i],":mismatch RBD or RBM")
                    pass
                else:
                    rbd_start = bs[0].start
                    rbd_end= be[0].end
                    #RBD
                    rbd_seq=s_protein[rbd_start:rbd_end]
                    rbd_seq=str(rbd_seq.seq)

                    rbm_start=xs[0].end                     
                    rbm_end= xe[0].start
                    #RBM
                    rbm_seq=s_protein[rbm_start:rbm_end]
                    rbm_seq=str(rbm_seq.seq)
                    # return end position by re will add 1; be careful
                    df.loc[i]= file_list[i], str(s_protein.seq), rbd_start, rbd_end, rbd_seq, len(rbd_seq), "" , "", rbm_seq, rbm_start, rbm_end, len(rbm_seq), "", "", "", "", "", ""

    # reset_index; prevent loc index error 
    df.reset_index(inplace=True)
    # map dict
    for i in range(len(df)):
        rbd_seq=list(df["RBD_seq"][i])
        rbm_seq=list(df["RBM_seq"][i])
        rbd_map=list([*map(col.get, rbd_seq)])
        rbm_map=list([*map(col.get, rbm_seq)])
        # build binding score
        rbd_score=[]
        rbm_score=[]
        # build expression score
        rbdx_score=[]
        rbmx_score=[]
        # map scores
        for j in range(len(rbd_map)):
            bind_rbd=wuhan_m[j,rbd_map[j]]
            rbd_score.append(bind_rbd)
            ex_rbd=wuhan_expm[j,rbd_map[j]]
            rbdx_score.append(ex_rbd)
        rbd_sum=sum(rbd_score)
        rbdx_sum=sum(rbdx_score)
        df["RBD_score_sum"][i], df["RBD_score"][i], df["RBD_ex_sum"][i], df["RBD_ex_score"][i] = rbd_sum, rbd_score, rbdx_sum, rbdx_score
        for k in range(len(rbm_map)):
            rbm_m=wuhan_m[107:176]
            rbmx_m=wuhan_expm[107:176]
            bind_rbm=rbm_m[k,rbm_map[k]]
            rbm_score.append(bind_rbm)
            ex_rbm=rbmx_m[k,rbm_map[k]]
            rbmx_score.append(ex_rbm)      
        rbm_sum=sum(rbm_score)
        rbmx_sum=sum(rbmx_score)
        df["RBM_score_sum"][i], df["RBM_score"][i],df["RBM_ex_sum"][i],df["RBM_ex_score"][i] = rbm_sum, rbm_score, rbmx_sum, rbmx_score
    df=df.drop(["index"],axis=1)
    df.to_csv(os.path.join( Output_path,"summary.csv"),index = False)
    

    
parser = argparse.ArgumentParser(prog='python RBD_RBM_matrix.py',formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i","--Input_folder",help="The path of the folder containing fasta files")
parser.add_argument("-o","--Output_path",help="The path of the output folder")
parser.add_argument("-d","--Data_matrix_folder",help="The path of the folder containing wuhan npy files")
args = parser.parse_args()
 
if __name__=="__main__":
    Input_folder = args.Input_folder
    Output_path = args.Output_path
    Data_matrix_folder = args.Data_matrix_folder
    RBD_RBM_matrix(Input_folder, Output_path, Data_matrix_folder)