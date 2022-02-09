import step_1
import step_1_browser
import step_2
import step_3
import sys
import getopt
import os

docstring = """
This script will apply a pipeline to reconstruct mutational spectrums.

For every given query, it will succesively call different scripts to make 
the BLAST queries, extract the consensus and the observed (and if possible
expected) deviations to finally reconstruct the mutation spectrums. To see 
a detailed description of the pipeline please check the gitlab: 
https://gitlab.epfl.ch/baffou/mutational-spectrum . Note that each step of
the pipeline can be run as a standalone script. Thus you can tweak intermediate
inputs and outputs as you want or even add some extra steps in between. For 
the I/O format of this script but also for the sub-scripts please check the
corresponding section in the gitlab.

Parameter:

-input_file: (String) Name of the control file containing the sequences. Please
    look at the gitlab for the file format.

Keyword Arguments:


-out_folder: (String) Name of the folder where both intermediate and final
    outputs will be stored (default=outputs).
-hits_nb: (int) Max number of blast hits to return per query (default=100).
    It must be one of the following values: [10, 50, 100, 250, 500, 1000, 5000].
-neighbours_nb: (int) Number of desired neighbours (on the left
    and on the right) in the output (default=2).
-hit_size_treshold: (int) Min number of hits for the sequences to be kept (default=10).
-codon: (flag) Indicate if codons should be provided in the output (default=False)
-b: (flag) Indicates if the step 1 should run in browser mode rather than 
    with biopython (default=False).
-v: (flag) Indicates if the browser should be visible (default=False).   
-c: (flag) Indicate if context should be taken into account (default=False).
-h: (flag) Opens this Help section.

Outputs:

The script will output 4 different files. Three of them are intermediate results:
    -sequences.csv: (DataFrame) Contains the results of all BLAST queries.
    -logs.csv: (DataFrame) Contains extra informations on queries.
    -deviations.json: (dict) Contains all observed and expected deviations.
The last one contains the mutational spectrums for each valid initial query: mut_spec.JSON

"""


def main(argv):
    # Parse arguments
    try:
        opts, args = getopt.getopt(
            argv, "hbcv", [
                "input_file=", "out_folder=", "hits_nb=",
                "neighbours_nb=", "hit_size_treshold=", "codon"
            ]
        )
    except getopt.GetoptError:
        print(argv)
        print("Error in inputs, please be sure to provide at least the path to the control file.")
        sys.exit(2)

    browser_mode = False
    context = False
    headless = True
    input_file = ""
    out_folder = "outputs"
    hit_nb = 100
    neighbours_nb = 2
    hit_size_treshold = 10
    codon = False
    for opt, arg in opts:
        if opt == "-h":
            print("Help section:\n", docstring)
            sys.exit()
        elif opt == "-b":
            browser_mode = True
        elif opt == "-c":
            context = True
        elif opt == "-v":
            headless = False
        elif opt == "--input_file":
            input_file = arg
        elif opt == "--out_folder":
            out_folder = arg
        elif opt == "--hits_nb":
            hit_nb = int(arg)
        elif opt == "--neighbours_nb":
            neighbours_nb = int(arg)
        elif opt == "--codon":
            codon = True
        elif opt == "--hit_size_treshold":
            hit_size_treshold = int(arg)
    # Inputs Assertion
    if neighbours_nb < 0 or neighbours_nb > 3:
        print("The number of neighbours is invalid, it should be comprised between 0 and 3.")
        sys.exit(2)
    if hit_nb < 1 or hit_nb > 5000:
        print("The number of hits should be in the interval from 1 to 5000.")
        sys.exit(2)
    if hit_size_treshold < 1 or hit_size_treshold > hit_nb:
        print("The min number of hits should be in the interval from 1 to hit_nb.")
        sys.exit(2)
    if not os.path.isfile(input_file):
        print("The control file doesn't exist, please provide an existing file")
        sys.exit(2)
    if os.path.isdir(out_folder):
        print(f"The folder {out_folder} already exists, please delete it or provide another name")
        _ans = input("Delete folder? (Y/n)\n")
        if _ans.lower() in {"", "y", "yes"}:
            os.rmdir(out_folder)
        else:
            print("Provide another folder or delete passed one manually")
            sys.exit(2)
    try:
        os.mkdir(out_folder)
    except:
        print(f"Impossible to create output dir {out_folder}. Please be sure of your permissions.")
        sys.exit(2)
    # Step 1
    step_1_arguments = ["--input_file", input_file, "--seq_out", out_folder+"/sequences.csv",
                        "--log_out", out_folder+"/logs.csv", "--hits_nb", hit_nb, "--hit_size_treshold", hit_size_treshold]
    if browser_mode:
        if not headless:
            step_1_arguments += ["-v"]
        step_1_browser.main(step_1_arguments)
    else:
        step_1.main(step_1_arguments)
    print("-"*72)
    # Step 2
    step_2_arguments = ["--seq_input", out_folder+"/sequences.csv", "--log_input", out_folder+"/logs.csv",
                        "--output", out_folder+"/deviations.json", "--neighbours_nb", str(neighbours_nb)]
    if codon:
        step_2_arguments += ["--codon"]
    step_2.main(step_2_arguments)
    print("-"*72)
    # Step 3
    step_3_arguments = ["--deviations", out_folder +
                        "/deviations.json", "--output", out_folder+"/mut_spec.json"]
    if context:
        step_3_arguments += ["-c"]
    step_3.main(step_3_arguments)
    print("-"*72)
    print(
        f"Pipeline executed with success! \nOutput files are located in {out_folder}")


if __name__ == "__main__":
    main(sys.argv[1:])
