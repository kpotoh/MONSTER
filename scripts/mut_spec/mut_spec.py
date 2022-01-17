import json
from collections import Counter
import sys
import getopt
import os


def main(argv):
    """This script will generate the mutational spectrum for each queries.

    Based on the observed and (when available) expected deviations, it will create the mutational spectrum
    for every given query. The input file should be in the same format as the output of the Step 2. Thus
    if you use this script in standalone mode, we recommend to check the gitlab (https://gitlab.epfl.ch/baffou/mutational-spectrum)
    to be sure of the I/O format. The files (both input and output) should be in JSON format.

    Parameter:

    -deviations: Mandatory parameter indicating the location of the .json file containing the deviations for each queries 
    (It should be the output file of the Consensus script, i.e. Step 2 of Pipeline)

    Keyword Arguments:

    -output: Output file name, should be a .json file (default mut_spec.json).
    -c: Indicate if context should be taken into account (default=False).

    Output:

    Single file name as output, containing the mutational spectrums for each inital queries.

    """
    try:
        opts, args = getopt.getopt(argv,"hc",["deviations=","output="])
    except getopt.GetoptError:
        print(argv)
        print("Error in inputs, please check the list of parameters with -h.")
        sys.exit(2)
    dev_input = ""
    output_name = "mut_spec.json"
    context = False
    for opt, arg in opts:
        if opt == "-h":
            print("\n")
            print("Help section:\n")
            print("Based on the observed and (when available) expected deviations, it will create the mutational spectrum \nfor every given query. The input file should be in the same format as the output of the Step 2.")
            print("Thus if you use this script in standalone mode, we recommend to check the gitlab (https://gitlab.epfl.ch/baffou/mutational-spectrum)")
            print("to be sure of the I/O format. The files (both input and output) should be in JSON format.")
            print("Arguments:\n")
            print("-deviations: Mandatory parameter indicating the location of the .json file containing the deviations for each queries \n (It should be the output file of the Consensus script, i.e. Step 2 of Pipeline)")
            print("-output: Output file name, should be a .json file.")
            print("-c: Indicate if context should be taken into account")
            print("\n\nOutputs:\n")
            print("The script will output a single file. It will be a .json containing the mutational frequencies.\n")
            sys.exit()
        elif opt == "-c":
            context = True
        elif opt == "--deviations":
            dev_input = arg
        elif opt == "--output":
            output_name = arg

    ## Inputs assertions
    if dev_input == "":
        print("The file input name containing the deviations should be provided (--deviations). Please check the help section (-h) to see mandatory parameters.")
        sys.exit(2)
    if not os.path.isfile(dev_input):
        print("The deviations input file doesn't exist, please provide an existing file.")
        sys.exit(2)
    if os.path.isfile(output_name):
        print(f"The file {output_name} already exists, please delete it or provide another name for the output file.")
        sys.exit(2)
    try:
        with open(dev_input, "r") as input_deviations:
            deviations = json.load(input_deviations)
    except:
        print("Cannot open the file containing the deviations, please check that \nthe format is the same as presented in the gitlab (https://gitlab.epfl.ch/baffou/mutational-spectrum)")
        sys.exit(2)
    if dev_input[-5:] != ".json" or output_name[-5:] != ".json":
        print("The files extensions must be .json, please be sure to provide them in the parameters.")
        sys.exit(2)
    
    # Script execution
    with open(dev_input, "r") as input_deviations:
        deviations = json.load(input_deviations)
    mutational_spectrums = dict()
    filtered_deviations = filter_deviations(deviations)
    for query_id, dev in filtered_deviations.items():
        print(f"Creating mutational spectrum for query {query_id}")
        query_context = context
        if dev["Observed_Deviations"]:
            if len(dev["Observed_Deviations"][0][0]) < 3: # Check that context is available (i.e. we have the neighbours a.a.)
                query_context = False
        mut_spec = create_mutational_spectrum(dev,query_context)
        formated_mut_spec = format_mut_spec(mut_spec)
        mutational_spectrums[query_id] = formated_mut_spec
    with open(output_name, "w") as dev_out_file:
        json.dump(mutational_spectrums, dev_out_file)
    print(f"Script executed with success! The output file is : {output_name}.")


def create_mutational_spectrum(deviations, context=False):
    """Create the mutational spectrum for a given query and its devitations.

    This function creates the mutational spectrum for a single query. It takes all its deviations and
    count the occurences of each type of SNP (with or without context depending on parameters). If the sequence
    is coding, it normalizes the result with the expected deviations occurences. The output format is a dictionnary
    whose keys are the possible mutations and the values their frequencies/amount.

    Parameters:

    -deviations: (dict) Deviations dictionnary for a given sequence.

    Keyword Arguments:

    -coding: (bool) Indicate if context should be taken into account (default=False).

    Output:

    -mut_spec: (dict) Mutational spectrum dictionnary for the given query

    """
    if not deviations["Observed_Deviations"]:
        print("No mutation of interest observed for this query.")
        default_output = dict.fromkeys(context_mutations,0.0) if context else dict.fromkeys(non_context_mutations,0.0)
        return default_output
    coding = "Expected_Deviations" in deviations.keys()
    mut_spec = dict()
    if coding: # Coding sequence case
        obs_deviations, exp_deviations = format_deviations(deviations,coding,context)
        counter_obs = Counter(obs_deviations)
        counter_exp = Counter(exp_deviations)
        mutations_list = context_mutations if context else non_context_mutations
        for dev in mutations_list:
            if counter_obs[dev] == 0:
                mut_spec[dev] = 0
            else:
                if counter_exp == 0:
                    mut_spec[dev] = -1
                else:
                    mut_spec[dev] = counter_obs[dev]/counter_exp[dev]
    else:
        counter = Counter(deviations["Observed_Deviations"])
        deviations_nb = len(deviations["Observed_Deviations"]) # !!! How to handle the case of non-coding mut spec?
        for dev in context_mutations:
            mut_spec[dev] = counter[dev]/deviations_nb
    return mut_spec


def format_mut_spec(mutational_spectrum):
    """Format annotation of mutational spectrum."""
    formated_mut_spec = dict()
    for deviation, freq in mutational_spectrum.items():
        formated_key = deviation[0]+">"+deviation[1] # Output will look like aAa>bBb if context and A>B if no context
        formated_mut_spec[formated_key] = freq
    return formated_mut_spec


def format_deviations(deviations,coding,context=False):
    """Format deviation to only keep nucleotide of interest.

    For the mutational spectrum, we are interested either in SNP (thus keeping only the nucleotide that mutated),
    or in SNP but with context (i.e. with one nucleotide on the left and one on the right). This function filter
    both observed and (if available) expected deviations based on the parameters. The output format is a list
    of tuples : (from,to). 

    Parameters:

    -deviations: (dict) Deviations dictionnary for a given sequence.
    -coding: (bool) Indicate if the sequences is coding or not.

    Keyword Arguments:

    -context: (bool) Indicate if context should be taken into account (default=False).

    Outputs:

    -obs_deviations: (list) Formated observed deviations.
    -exp_deviations: (list) Formated expected deviations (always empty if non-coding).

    """
    mut_pos = int((len(deviations["Observed_Deviations"][0][0])-1)/2) # Take the nucleotide in the middle
    obs_deviations = []
    exp_deviations = []
    if context:
        obs_deviations = [(d[0][mut_pos-1:mut_pos+2],d[1][mut_pos-1:mut_pos+2]) for d in deviations["Observed_Deviations"]]
        if coding :
            exp_deviations = [(d[0][mut_pos-1:mut_pos+2],d[1][mut_pos-1:mut_pos+2]) for d in deviations["Expected_Deviations"]]
    else:
        obs_deviations = [(d[0][mut_pos],d[1][mut_pos]) for d in deviations["Observed_Deviations"]]
        if coding:
            exp_deviations = [(d[0][mut_pos],d[1][mut_pos]) for d in deviations["Expected_Deviations"]]
    
    return obs_deviations, exp_deviations


def filter_deviations(deviations):
    """Filter out non-synonymous deviations.

    If the sequences is a coding, then keep only synonymous deviations. Otherwise keep every deviations.
    Additionally for each query, keep only deviations (i.e. From/To fields).

    Parameters:

    -deviations: dictionnary of all the deviations for the original queries (keys = query_ids and values = deviations and consensus.)

    Output:

    -new_deviations: dictionnary of selected deviations (keys = query_ids and values = deviations)

    """
    new_deviations = dict()
    for query_id,query_dev in deviations.items():
        query_deviations = {}
        if "Expected_Deviations" in query_dev.keys(): # Coding sequence case
            query_deviations["Observed_Deviations"] = [(d["From"],d["To"]) for d in query_dev["Observed_Deviations"] if d["Synonymous"]]
            query_deviations["Expected_Deviations"] = [(d["From"],d["To"]) for d in query_dev["Expected_Deviations"] if d["Synonymous"]]
        else:
            query_deviations["Observed_Deviations"] = [(d["From"],d["To"]) for d in query_dev["Observed_Deviations"]]
        new_deviations[query_id] = query_deviations
    return new_deviations


def context_mutations_gen():
    """Return all SNP mutations with context."""
    bases = ["A","T","C","G"]
    possible_mutations = []
    for context_left in bases:
        for b_from in bases:
            for b_to in bases:
                if b_to != b_from:
                    for context_right in bases:
                        original = context_left.casefold()+b_from+context_right.casefold()
                        mutated = context_left.casefold()+b_to+context_right.casefold()
                        possible_mutations.append((original,mutated))
    return possible_mutations
                    
    
def non_context_mutations_gen():
    """Return all SNP mutations."""
    bases = ["A","T","C","G"]
    possible_mutations = []
    for b_from in bases:
        for b_to in bases:
            if b_from != b_to:
                possible_mutations.append((b_from,b_to))
    return possible_mutations


context_mutations = context_mutations_gen() # List of SNP mutations with context
non_context_mutations = non_context_mutations_gen() # List of SNP mutations

if __name__ == "__main__":
    main(sys.argv[1:])