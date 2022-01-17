import matplotlib.pyplot as plt
import json
import sys
import getopt
import os

def main(argv):
    try:
        opts, args = getopt.getopt(argv,"h",["input="])
    except getopt.GetoptError:
        print(argv)
        print("Error in inputs, please check the list of parameters with -h.")
        sys.exit(2)
    mut_spec_input = ""
    for opt, arg in opts:
        if opt == "-h":
            print("\n")
            print("Help section:\n\n")
            print("This script provide a really simple visualisation of the mutational spectrums of different queries.\n\n")
            print("Paramters:\n")
            print("-input: File name of the mutation spectrums. (Should be a .json file in the same format as the output of Step 2)\n")
            sys.exit()
        elif opt == "--input":
            mut_spec_input = arg
    if mut_spec_input == "":
        print("The file input name containing the mutational spectrums should be provided (--input). Please check the help section (-h) to see mandatory parameters.")
        sys.exit(2)
    ## inputs assertions
    if not os.path.isfile(mut_spec_input):
        print("The deviations input file doesn't exist, please provide an existing file.")
        sys.exit(2)
    try:
        with open(mut_spec_input, "r") as input_spec:
            mut_spec = json.load(input_spec)
    except:
        print("Cannot open the mutational spectrums please check that the format is the same as presented in the gitlab (https://gitlab.epfl.ch/baffou/mutational-spectrum)")
        sys.exit(2)
    if mut_spec_input[-5:] != ".json":
        print("The file extension must be .json, please be sure to provide them in the parameters.")
        sys.exit(2)
    
    # script
    with open(mut_spec_input, "r") as input_spec:
        mut_spec = json.load(input_spec)
    
    fig, axs = plt.subplots(len(mut_spec.keys()), figsize=(16,12))
    i = 0
    for query_id,mutations in mut_spec.items():
        mut_dict = mutations
        axs[i].bar(range(len(mut_dict)), list(mut_dict.values()), align='center')
        axs[i].set_xticks(range(len(mut_dict)))
        axs[i].set_xticklabels(list(mut_dict.keys()), rotation=-90)
        axs[i].set_title(f"Query {query_id}")
        axs[i].set_xlabel("Mutations")
        axs[i].set_ylabel("Amount")
        i += 1
    fig.tight_layout(pad=10.0)
    plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])