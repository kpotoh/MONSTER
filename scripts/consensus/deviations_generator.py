import numpy as np
import pandas as pd
from Bio import Seq
import sys
import getopt
import os
import json

IUAPC_table = ["A","T","C","G","U","R","Y","K","M","S","W","B","D","H","V","N","-"]

def main(argv):
    # Parse arguments
    try:
        opts, args = getopt.getopt(argv,"h",["seq_input=","log_input=","output="])
    except getopt.GetoptError:
        print(argv)
        print("Error in inputs, please check the list of parameters with -h.")
        sys.exit(2)
    seq_input = ""
    log_input = ""
    output_name = "deviations.json"
    for opt, arg in opts:
        if opt == "-h":
            print("\n")
            print("Help section:\n")
            print("This script will derive the consensus and all the deviations observed in the sequences.\nFor coding sequences it will also produce the expected deviations.")
            print("We recommend to check the gitlab (https://gitlab.epfl.ch/baffou/mutational-spectrum) to be sure of the I/O format.\n\n")
            print("Inputs:")
            print("-seq_input : mandatory parameter indicating the location of the .csv file containing the sequences (It should be the output file of the BLAST script)")
            print("-log_input: mandatory parameter indicating the location of the .csv file containing the sequences (It should be the output file of the BLAST script)")
            print("-output: name of the output file for the deviations. It will be a json file so don't forget the .json extension. (Default: deviations.json)\n\n")
            print("Outputs:")
            print("The script will output a single file. It will be a csv file containing for each original query, the consensus sequence with the observed deviations \n and if possible the expected deviations.\n")
            sys.exit()
        elif opt == "--seq_input":
            seq_input = arg
        elif opt == "--log_input":
            log_input = arg
        elif opt == "--output":
            output_name = arg
    ## Inputs assertions
    if seq_input == "":
        print("The sequences input has not been provided. Please check the help section (-h) to see mandatory parameters.")
        sys.exit(2)
    if log_input == "":
        print("The log input has not been provided. Please check the help section (-h) to see mandatory parameters.")
        sys.exit(2)
    if not os.path.isfile(seq_input):
        print("The sequences input file doesn't exist, please provide an existing file.")
        sys.exit(2)
    if not os.path.isfile(log_input):
        print("The log input file doesn't exist, please provide an existing file.")
        sys.exit(2)
    if os.path.isfile(output_name):
        print(f"The file {output_name} already exists, please delete it or provide another name.")
        sys.exit(2)
    try:
        pd.read_csv(seq_input)
    except:
        print("Cannot open the sequences file, please check that the format is the same as presented in the gitlab (https://gitlab.epfl.ch/baffou/mutational-spectrum)")
        sys.exit(2)
    try:
        pd.read_csv(log_input)
    except:
        print("Cannot open the log file, please check that the format is the same as presented in the gitlab (https://gitlab.epfl.ch/baffou/mutational-spectrum)")
        sys.exit(2)
    # Script execution 
    deviations = generate_deviations(seq_input, log_input)
    with open(output_name, "w") as dev_out_file:
        json.dump(deviations, dev_out_file)
    print(f"Script executed with success! The output file is : {output_name}.")

def generate_deviations(seq_file_name, log_file_name):
    """
    Given sequences, produces a list of all deviations for each of the original queries.
    If the query is a coding sequence, then it will also list the amino acids deviations 
    and the expected deviations.
    
    Parameters:
        seq_file_name : Name of the file containing the sequences. It should be in the same format as the BLAST script output (check gitlab : https://gitlab.epfl.ch/baffou/mutational-spectrum)
        log_file_name : Name of the file containing the log for the sequences. It should be in the same format as the BLAST script output (check gitlab : https://gitlab.epfl.ch/baffou/mutational-spectrum)
    
    Returns:
        output : dictionnary whose keys are the query ids, and the values the consensus for the given query and the deviations.
    """
    sequences_df = pd.read_csv(seq_file_name)
    log_df = pd.read_csv(log_file_name)
    deviations = dict()
    for query_id in log_df.index: # Create deviations list per original queries
        query_seqs = sequences_df[sequences_df["Control_ID"] == query_id] 
        seqs_list = [q for q in query_seqs["Sequence"]]
        seqs_df = dataframe_builder(seqs_list)
        consensus, maf_table = consensus_builder(seqs_df)
        coding = log_df.loc[query_id]["CDS"] == 1
        gen_code = log_df.loc[query_id]["Code"] 
        deviations_obs = create_observed_deviations(query_seqs,consensus,maf_table,coding,gen_code)
        deviations[int(query_id)] = {"Consensus" : consensus,"Observed_Deviations" : deviations_obs}
        if coding:
            expected_deviations = create_expected_deviations(consensus,gen_code)
            deviations[int(query_id)]["Expected_Deviations"] = expected_deviations      
    return deviations

def dataframe_builder(seq_list):
    """
    Given a list of sequences, output a pandas DataFrame with a column for every site
    
    Parameters:
        seq_list : list of sequences (string). Be carefull, every sequences must have the same length
    
    Returns:
        output : pandas Dataframe with each sequence and a column per site
    """
    formated_seq_list = list(map(lambda s : list(s), seq_list))
    return pd.DataFrame(formated_seq_list)

def consensus_builder(seq_df):
    """
    Given a pandas DataFrame containing the sequences, it will create the consensus sequence and the Minor Allele Frequency table
    
    Parameters:
        seq_df : pandas Dataframe with each sequence and a column per site
        
    Returns:
        consensus : consensus sequence(string)
        maf_table : minor allele frequency table(dict) (an entry per site)
    """
    consensus = ""
    maf_table = dict()
    seq_nb = seq_df.shape[0]
    for c in seq_df.columns:
        nucleotide_count = seq_df[c].value_counts()
        consensus += nucleotide_count.index[0] # !!! It may be biased because in case of equality we always take the first index -> randomize if equality?
        nuc_dict = dict()
        for nuc in IUAPC_table: # !!! Need to see how to handle new letters (R,N,...)
            if nuc in nucleotide_count:
                nuc_dict[nuc] = nucleotide_count.loc[nuc]/seq_nb
            else:
                nuc_dict[nuc] = 0.0
        maf_table[c] = nuc_dict
    return consensus, maf_table

def create_observed_deviations(seq_df,consensus,maf_table,coding,genetic_code,neighbours_nb=2):
    """
    It will list all the observed deviations from consensus. For coding sequences, the deviations are in the format : {Seq_ID,Codon_Pos,Nucleotide_Pos,From,To,MAF,Amino_Acid_From,Amino_Acid_To,Synonymous}.
    For non-coding sequences, the deviations are in the format : {Seq_ID,Codon_Pos,Nucleotide_Pos,From,To,MAF} (because we don't have codon alignment).
    where "Seq_ID" is the id of the sequence, "Codon_Pos" is the codon position, "Nucleotide_Pos" is the nucleotide position, "From" is the original value in the consensus, "To" is 
    the value in the deviations, Amino_Acid_From is the original amino acid in the consensus, Amino_Acid_To is the amino acid in the deviation,
    "Synonymous" indicates if the mutation is synonymous or not (boolean).
    We make extensive use of numpy masking and slicing thus if the code isn't clear, check the corresponding gitlab (https://gitlab.epfl.ch/baffou/mutational-spectrum) section.
    
    Parameters:
        seq_df : pandas DataFrame with every sequences in a single column
        consensus : consensus sequence (string)
        maf_table : minor allele frequency table (dict)
        neighbours_nb: number(int) of desired neighbours in the output, to the left and to the right (default : 2)
    
    Returns:
        deviations : list of dictionnaries containing every deviations for the given sequence
    """
    deviations = []
    sequences = np.asarray([list("#"*neighbours_nb+q.casefold()+"#"*neighbours_nb) for q in seq_df["Sequence"]]) # Pad sequences to take care of extremities
    formated_consensus = "#"*neighbours_nb+consensus.casefold()+"#"*neighbours_nb
    deviations_indices = np.where(sequences != list(formated_consensus))
    start_indices = deviations_indices[1]-neighbours_nb # Sub-sequence starting point
    end_indices = deviations_indices[1]+neighbours_nb# Sub-sequence ending point
    columns_indices = np.arange(sequences.shape[1])
    selection_mask = (start_indices[:,None] <= columns_indices) & (end_indices[:,None] >= columns_indices) # Mask used for slicing the sequences
    deviations_array = sequences[deviations_indices[0]][selection_mask].reshape(deviations_indices[0].shape[0],2*neighbours_nb+1) # Take valid sequences and reshape to have constant output
    deviations_subseq = ["".join(s) for s in deviations_array]   
    deviations_subseq = list(map(lambda s : s[:neighbours_nb]+s[neighbours_nb:].capitalize(),deviations_subseq)) # Format in ..xxXxx..
    for i in range(len(deviations_subseq)):
        seq_id = int(deviations_indices[0][i])
        pos_in_seq = int(deviations_indices[1][i]-neighbours_nb) # The "-neighbours_nb" is to take care of padding and its index shifting
        origin_seq = formated_consensus[pos_in_seq:pos_in_seq+2*neighbours_nb+1]
        origin_seq = origin_seq[:neighbours_nb]+origin_seq[neighbours_nb:].capitalize() # Format in ..xxXxx..
        dev_seq = deviations_subseq[i]
        if coding:
            codon_start = int(pos_in_seq/3)
            original_codon = consensus[codon_start*3:(codon_start+1)*3].upper()
            original_codon = original_codon.replace("-","N")
            original_codon = original_codon.replace("#","N")
            new_codon = "".join(sequences[seq_id][neighbours_nb+codon_start*3:neighbours_nb+(codon_start+1)*3]).upper()
            new_codon = new_codon.replace("-","N")
            new_codon = new_codon.replace("#","N")
            aa_from = Seq.translate(original_codon,table=genetic_code)
            aa_to = Seq.translate(new_codon,table=genetic_code)
            synonymous = aa_from == aa_to
            deviations.append({"Seq_ID" : seq_id, "Codon_Pos" : codon_start+1, "Nucleotide_Pos" : pos_in_seq+1, "From" : origin_seq, "To" : dev_seq, "MAF" : maf_table[pos_in_seq][dev_seq[neighbours_nb]], "Amino_Acid_From" : aa_from, "Amino_Acid_To" : aa_to, "Synonymous" : synonymous})
        else:
            deviations.append({"Seq_ID" : seq_id, "Nucleotide_Pos" : pos_in_seq+1, "From" : origin_seq, "To" : dev_seq, "MAF" : maf_table[pos_in_seq][dev_seq[neighbours_nb]]})
    return deviations

def create_expected_deviations(consensus,genetic_code):
    """
    Create the exèected deviations list based on the consensus sequence. The deviations are in the format : {Codon_Pos,Nucleotide_Pos,From,To,Amino_Acid_From,Amino_Acid_To,Synonymous},
    where "Codon_Pos" is the codon position, "Nucleotide_Pos" is the nucleotide position, "From" is the original value in the consensus, "To" is 
    the value in the deviations, Amino_Acid_From is the original codon in the consensus, Amino_Acid_To is the codon in the deviation,
    "Synonymous" indicates if the mutation is synonymous or not (boolean).
    
    Parameters:
        consensus : consensus sequence (string)
        genetic_code : NCBI genetic code id (int) (can be found in https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
    
    Returns:
        deviations : list of dictionnaries containing every deviations from the consensus sequence
    """
    deviations = []
    for i in range(0,len(consensus),3): # for each codon we deviate each one of the three nucleotide to the three other possible ones
        for j in range(3):
            nucleotide_list = ["A","T","C","G"]
            nucleotide_list.remove(consensus[i+j])
            for nuc in nucleotide_list:
                original_codon = consensus[i:i+3]
                original_codon = original_codon.replace("-","N")
                new_codon = consensus[i:i+j]+nuc+consensus[i+j+1:i+3]
                new_codon = new_codon.replace("-","N")
                synonymous = True
                origin_translation = Seq.translate(original_codon,genetic_code)
                new_translation = Seq.translate(new_codon,genetic_code)
                if origin_translation != new_translation:
                    synonymous = False
                origin = consensus.casefold()
                left_index = i+j-2 if i+j-2 >= 0 else 0 # we have to make a index correction when we are at the left extremity (the right one is handled by python)
                origin = origin[left_index:i+j]+origin[i+j:i+j+3].capitalize()
                deviated = consensus.casefold()
                deviated = deviated[left_index:i+j]+nuc+deviated[i+j+1:i+j+3]
                deviation_dict = {"Codon_Pos" : i+1, "Nucleotide_Pos" : i+j+1, "From" : origin, "To" : deviated, "Amino_Acid_From" : origin_translation, "Amino_Acid_To" : new_translation, "Synonymous" : synonymous}
                deviations.append(deviation_dict)
    return deviations

if __name__ == "__main__":
    main(sys.argv[1:])