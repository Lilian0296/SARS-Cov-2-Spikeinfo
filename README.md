# SARS-Cov-2-Spikeinfo
Calculate the binding and expression scores of RBD and RBM to ACE2

-----------
A surface spike protein regulates the entry of SARS-Cov-2 into host cells.
* RBD: receptor-binding domain of the spike protein of SARS-Cov-2
* RBM: receptor-binding motif of the spike protein of SARS-Cov-2
#### Data matrix was built from the original Wuhan strain. NA was filled by median. 
#### Raw data matrix from the paper: Deep mutational scanning of SARS-CoV-2 receptor binding domain reveals constraints on folding and ACE2 binding



### Notice
1. Input is the path of the folder containing multiple fasta files. (These fasta files cannot be necessary in-frame.)
2. RBD and RBM search employ fuzzy matches with a maximum Levenshtein Distance of 1.
3. Mismatch fasta files will be printed while running the script.


### Running

```
python RBD_RBM_matrix.py -i ./input_folder/ -o ./output_path/ -d ./data_matrix_folder/

# -i, --Input_folder, The path of the folder containing fasta files
# -o, --Output_path, The path of the output folder
# -d, --Data_matrix_folder, The path of the folder containing wuhan npy files

```

### Output
#### Output a summary csv file containing following information:
1. File name
2. S_sequences
3. RBD and RBM sequences
4. The position of RBD start and end 
5. The position of RBM start and end 
6. The length of RBD and RBM
7. The binding score of RBD and RBM to ACE2 
  (RBD_score_sum, RBD_score (for each amino acid position);
   RBM_score_sum, RBM_score (for each amino acid position))
8. The expression score of RBD and RBM
  (RBD_ex_sum, RBD_ex_score (for each amino acid position);
   RBM_ex_sum, RBM_ex_score (for each amino acid position))


### References
----------

[1] Starr, T. N., Greaney, A. J., Hilton, S. K., Ellis, D., Crawford, K. H., Dingens, A. S., ... & Bloom, J. D. (2020). Deep mutational scanning of SARS-CoV-2 receptor binding domain reveals constraints on folding and ACE2 binding. cell, 182(5), 1295-1310
<https://doi.org/10.1016/j.cell.2020.08.012>
