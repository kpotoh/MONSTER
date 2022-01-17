import numpy as np
import pandas as pd
from Bio import Seq
import sys
import getopt
import os
import json

IUAPC_table = ["A", "T", "C", "G", "U", "R", "Y",
               "K", "M", "S", "W", "B", "D", "H", "V", "N", "-"]


def main(argv):
    """This script will derive all the mutations and consensus for each queries.

    For every given query, it will generate the observed and (when possible)
    the expected deviations. The input file should be in the same format as
    the output of the Step 1. Thus if you use this script in standalone mode,
    we recommend to check the gitlab (https://gitlab.epfl.ch/baffou/mutational-spectrum)
    to be sure of the I/O format. The input files should be in CSV format and
    the output files in JSON format.

    Parameter:

    -seq_input: (String) Mandatory parameter indicating the location of the .csv 
        file containing the sequences (It should be the output file of the Step 1 of Pipeline).
    -log_input: (String) Mandatory parameter indicating the location of the .csv 
        file containing the logs (It should be the output file of the Step 1 of Pipeline).

    Keyword Arguments:

    -output: (String) Output file name, should be a .json file (default deviations.json).
    -neighbours_nb: (int) Number of desired neighbours (on the left
        and on the right) in the output (default : 2).
    -codon: (flag) Indicate if codons should be provided in the output (default : False)
    -h: (flag) Opens this Help section.
    
    Output:

    Single file name as output, containing the deviations for each inital queries.

    """
    # Parse arguments
    try:
        opts, args = getopt.getopt(
            argv, "h", ["seq_input=", "log_input=", "output=", "neighbours_nb=", "codon"])
    except getopt.GetoptError:
        print(argv)
        print("Error in inputs, please check the list of parameters with -h.")
        sys.exit(2)
    seq_input = ""
    log_input = ""
    output_name = "deviations.json"
    neighbours_nb = 2
    codon = False
    for opt, arg in opts:
        if opt == "-h":
            print("\n")
            print("Help section:\n")
            print("This script will derive all the mutations and consensus for each queries.\n")
            print("For every given query, it will generate the observed and (when possible)")
            print("the expected deviations. The input file should be in the same format as")
            print("the output of the Step 1. Thus if you use this script in standalone mode,")
            print("we recommend to check the gitlab (https://gitlab.epfl.ch/baffou/mutational-spectrum)")
            print("to be sure of the I/O format. The input files should be in CSV format and")
            print("the output files in JSON format.\n")
            print("Parameter:\n")
            print("-seq_input: (String) Mandatory parameter indicating the location of the .csv") 
            print("    file containing the sequences (It should be the output file of the Step 1 of Pipeline).")
            print("-log_input: (String) Mandatory parameter indicating the location of the .csv")
            print("    file containing the logs (It should be the output file of the Step 1 of Pipeline).\n")
            print("Keyword Arguments:\n")
            print("-output: (String) Output file name, should be a .json file (default deviations.json).")
            print("-neighbours_nb: (int) Number of desired neighbours (on the left")
            print("    and on the right) in the output (default : 2).")
            print("-codon: (flag) Indicate if codons should be provided in the output (default : False).")
            print("-h: (flag) Opens this Help section.\n")
            print("Output:\n")
            print("Single file name as output, containing the deviations for each inital queries.")
            sys.exit()
        elif opt == "--seq_input":
            seq_input = arg
        elif opt == "--log_input":
            log_input = arg
        elif opt == "--output":
            output_name = arg
        elif opt == "--neighbours_nb":
            neighbours_nb = int(arg)
        elif opt == "--codon":
            codon = True
    # Inputs assertions
    if seq_input == "":
        print("The sequences input has not been provided. Please check the help",end="") 
        print("section (-h) to see mandatory parameters.")
        sys.exit(2)
    if log_input == "":
        print("The log input has not been provided. Please check the help",end="") 
        print("section (-h) to see mandatory parameters.")
        sys.exit(2)
    if not os.path.isfile(seq_input):
        print("The sequences input file doesn't exist, please provide an existing file.")
        sys.exit(2)
    if not os.path.isfile(log_input):
        print("The log input file doesn't exist, please provide an existing file.")
        sys.exit(2)
    if os.path.isfile(output_name):
        print(
            f"The file {output_name} already exists, please delete it or provide another name.")
        sys.exit(2)
    if neighbours_nb < 0 or neighbours_nb > 3:
        print(
            "The number of neighbours is invalid, it should be comprised between 0 and 3.")
        sys.exit(2)
    try:
        pd.read_csv(seq_input)
    except:
        print("Cannot open the sequences file, please check that the format",end="")
        print("is the same as presented in the gitlab",end="") 
        print("(https://gitlab.epfl.ch/baffou/mutational-spectrum)")
        sys.exit(2)
    try:
        pd.read_csv(log_input)
    except:
        print("Cannot open the logs file, please check that the format",end="")
        print("is the same as presented in the gitlab",end="") 
        print("(https://gitlab.epfl.ch/baffou/mutational-spectrum)")
        sys.exit(2)
    # Script execution
    deviations = generate_deviations(
        seq_input, log_input, codon=codon, neighbours_nb=neighbours_nb)
    with open(output_name, "w") as dev_out_file:
        json.dump(deviations, dev_out_file)
    print(f"Step 2 executed with success!\nThe output file is : {output_name}.")


def generate_deviations(seq_file_name, log_file_name, codon=False, neighbours_nb=2):
    """ Generate observed and, if possible, expected deviations.

    Given some sequences, it produces a list of all deviations for each of 
    the original queries. If the query is a coding sequence, then it will
    also list the amino acids deviations and the expected deviations.

    Parameters:

    -seq_file_name: (String) Name of the file containing the sequences. 
        The format should be as specified in our gitlab : 
        https://gitlab.epfl.ch/baffou/mutational-spectrum.
    -log_file_name: (String) Name of the file containing the log for the sequences. 
        The format should be as specified in our gitlab : 
        https://gitlab.epfl.ch/baffou/mutational-spectrum.

    Keywords Arguments:

    -codon: (bool) Indicate if codons should be provided in the output (default : False)
    -neighbours_nb: (int) Number of desired neighbours (on the left
        and on the right) in the output (default : 2).

    Output:
        
    -output : (dict) Contains all deviations for each query. The keys are the 
        query ids, and the values the corresponding consensus and deviations.

    """
    sequences_df = pd.read_csv(seq_file_name)
    sequences_df = sequences_df.set_index("Unnamed: 0")
    log_df = pd.read_csv(log_file_name)
    log_df = log_df.set_index("Unnamed: 0")
    deviations = dict()
    for query_id in log_df.index:  # Create deviations list per original queries
        query_seqs = sequences_df[sequences_df["Control_ID"] == query_id]
        seqs_list = [q for q in query_seqs["Sequence"]]
        seqs_df = dataframe_builder(seqs_list)
        consensus, maf_table = consensus_builder(seqs_df)
        coding = log_df.loc[query_id]["CDS"] == 1
        gen_code = log_df.loc[query_id]["Code"]
        deviations_obs = create_observed_deviations(
        query_seqs, consensus, maf_table, coding, gen_code, codon=codon, neighbours_nb=neighbours_nb)
        deviations[int(query_id)] = {
            "Consensus": consensus, "Observed_Deviations": deviations_obs}
        if coding:
            expected_deviations = create_expected_deviations(
            consensus, gen_code, codon=codon, neighbours_nb=neighbours_nb)
            deviations[int(query_id)
                ]["Expected_Deviations"] = expected_deviations
    return deviations


def create_observed_deviations(seq_df, consensus, maf_table, coding, genetic_code, codon=False, neighbours_nb=2):
    """Derive observed deviations from given sequences and the consensus.
    
    It will list all the observed deviations from consensus. The output format
    is different depending if the original sequence was coding or not. For a 
    detailed output format, please check the gitlab: (https://gitlab.epfl.ch/baffou/mutational-spectrum).
    We make extensive use of numpy masking and slicing thus if the code isn't
    clear, check the corresponding gitlab section.

    Parameters:

    -seq_df: (Dataframe) Contains each sequences with a column per site.
    -consensus: (String) Consensus sequence.
    -maf_table: (dict) Minor allele frequency table (keys are site and values are alleles freq).
    -coding: (bool) Indicates if the original sequence was coding.
    -genetic_code: (int) Id of the genetic coding table to use. 
        The list of the different genetic codes can be found on 
        (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). 

    Keywords Arguments:

    -codon: (bool) Indicate if codons should be provided in the output (default : False)
    -neighbours_nb: (int) Number of desired neighbours (on the left
        and on the right) in the output (default : 2).

    Output:

    -deviations: (list) Contains dictionaries for each deviations for the given sequence.

    """
    deviations = []
    # Pad sequences with '#' symbol to add neighbours to extremities
    sequences = np.asarray([list("#"*neighbours_nb+q.casefold()+"#"*neighbours_nb)
                            for q in seq_df["Sequence"]])  
    formated_consensus = "#"*neighbours_nb+consensus.casefold()+"#"*neighbours_nb
    deviations_indices = np.where(sequences != list(formated_consensus))
    start_indices = deviations_indices[1]-neighbours_nb
    end_indices = deviations_indices[1] + neighbours_nb 
    columns_indices = np.arange(sequences.shape[1])
    selection_mask = (start_indices[:, None] <= columns_indices) & (
        end_indices[:, None] >= columns_indices)
    # Take valid sequences with masking and reshape to have constant size outputs
    deviations_array = sequences[deviations_indices[0]][selection_mask].reshape(
        deviations_indices[0].shape[0], 2*neighbours_nb+1)  
    deviations_subseq = ["".join(s) for s in deviations_array]
    # Format in ..xxXxx..
    deviations_subseq = list(map(
        lambda s: s[:neighbours_nb]+s[neighbours_nb:].capitalize(), deviations_subseq))
    for i in range(len(deviations_subseq)):
        seq_id = int(deviations_indices[0][i])
        # The "-neighbours_nb" is to take care of padding and its index shifting effect
        pos_in_seq = int(deviations_indices[1][i]-neighbours_nb)
        origin_seq = formated_consensus[pos_in_seq:pos_in_seq +
                                        2*neighbours_nb+1]
        # Format in ..xxXxx..
        origin_seq = (origin_seq[:neighbours_nb]
                        + origin_seq[neighbours_nb:].capitalize())
        dev_seq = deviations_subseq[i]
        if coding:
            codon_start = int(pos_in_seq/3)
            original_codon = consensus[codon_start*3:(codon_start+1)*3].upper()
            original_codon = original_codon.replace("-", "N")
            original_codon = original_codon.replace("#", "N")
            new_codon = "".join(
                sequences[seq_id][neighbours_nb+codon_start*3:neighbours_nb+(codon_start+1)*3]).upper()
            new_codon = new_codon.replace("-", "N")
            new_codon = new_codon.replace("#", "N")
            aa_from = Seq.translate(original_codon, table=genetic_code)
            aa_to = Seq.translate(new_codon, table=genetic_code)
            synonymous = aa_from == aa_to
            deviation_dict = {"Seq_ID": seq_id, "Codon_Pos": codon_start+1, "Nucleotide_Pos": pos_in_seq+1,
                                "From": origin_seq, "To": dev_seq, "MAF": maf_table[pos_in_seq][dev_seq[neighbours_nb]],
                                "Amino_Acid_From": aa_from, "Amino_Acid_To": aa_to, "Synonymous": synonymous}
            if codon:
                deviation_dict["Codon_From"] = original_codon
                deviation_dict["Codon_To"] = new_codon
            deviations.append(deviation_dict)
        else:
            deviations.append({"Seq_ID": seq_id, "Nucleotide_Pos": pos_in_seq+1, "From": origin_seq,
                               "To": dev_seq, "MAF": maf_table[pos_in_seq][dev_seq[neighbours_nb]]})
    return deviations


def create_expected_deviations(consensus, genetic_code, codon=False, neighbours_nb=2):
    """Derive expected deviations from given consensus.
    
    It will list all the expected deviations from consensus. For a 
    detailed output format, please check the gitlab: (https://gitlab.epfl.ch/baffou/mutational-spectrum).

    Parameters:

    -consensus: (String) Consensus sequence.
    -genetic_code: (int) Id of the genetic coding table to use. 
        The list of the different genetic codes can be found on 
        (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). 

    Keywords Arguments:

    -codon: (bool) Indicate if codons should be provided in the output (default : False)
    -neighbours_nb: (int) Number of desired neighbours (on the left
        and on the right) in the output (default : 2).

    Output:

    -deviations: (list) Contains dictionaries for each deviations for the given sequence.

    """
    deviations = []
    # For each codon we deviate each one its three nucleotides
    # To the three other possible nucleotides.
    for i in range(0, len(consensus), 3):
        for j in range(3):
            if consensus[i+j] != "-": # We don't want to generate deviations from gaps
                nucleotide_list = ["A", "T", "C", "G"]
                nucleotide_list.remove(consensus[i+j])
                for nuc in nucleotide_list:
                    original_codon = consensus[i:i+3]
                    original_codon = original_codon.replace("-", "N")
                    new_codon = consensus[i:i+j]+nuc+consensus[i+j+1:i+3]
                    new_codon = new_codon.replace("-", "N")
                    origin_translation = Seq.translate(
                        original_codon, genetic_code)
                    new_translation = Seq.translate(new_codon, genetic_code)
                    synonymous = (origin_translation == new_translation)
                    origin = consensus.casefold()
                    # We make an index correction when we are 
                    # At the left extremity (the right one is handled by python)
                    left_index = i+j-neighbours_nb if i+j-neighbours_nb >= 0 else 0
                    origin = (origin[left_index:i+j]
                            + origin[i+j:i+j+neighbours_nb+1].capitalize())
                    deviated = consensus.casefold()
                    deviated = (deviated[left_index:i+j]+nuc
                            + deviated[i+j+1:i+j+neighbours_nb+1])
                    deviation_dict = {"Codon_Pos": i+1, "Nucleotide_Pos": i+j+1, "From": origin, "To": deviated,
                                    "Amino_Acid_From": origin_translation, "Amino_Acid_To": new_translation, 
                                    "Synonymous": synonymous}
                    if codon:
                        deviation_dict["Codon_From"] = original_codon
                        deviation_dict["Codon_To"] = new_codon
                    deviations.append(deviation_dict)
    return deviations


def dataframe_builder(seq_list):
    """Create a DataFrame with a column per site."""
    formated_seq_list = list(map(lambda s: list(s), seq_list))
    return pd.DataFrame(formated_seq_list)


def consensus_builder(seq_df):
    """Create the consensus and the MAF table.
    
    Given a pandas DataFrame containing the sequences, it will create 
    the consensus sequence and the Minor Allele Frequency table.

    Parameters:
        
    -seq_df : (Dataframe) Contains each sequences with a column per site.

    Outputs:
        
    -consensus: (String) Consensus sequence.
    -maf_table: (dict) Minor allele frequency table (keys are site and values are alleles freq).

    """
    consensus = ""
    maf_table = dict()
    seq_nb = seq_df.shape[0]
    for c in seq_df.columns:
        nucleotide_count = seq_df[c].value_counts()
        # !!! It may be biased because in case of equality we always take the first index -> randomize if equality?
        consensus += nucleotide_count.index[0]
        nuc_dict = dict()
        for nuc in IUAPC_table:
            if nuc in nucleotide_count:
                nuc_dict[nuc] = nucleotide_count.loc[nuc]/seq_nb
            else:
                nuc_dict[nuc] = 0.0
        maf_table[c] = nuc_dict
    return consensus, maf_table


if __name__ == "__main__":
    main(sys.argv[1:])
