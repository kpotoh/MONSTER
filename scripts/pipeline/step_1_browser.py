import time
import pathlib
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.chrome.options import Options
from webdriver_manager.chrome import ChromeDriverManager
from selenium import webdriver
import selenium
import os
import pandas as pd
import numpy as np
import sys
import getopt
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Seq
from Bio import SeqIO
from Bio import Entrez

from step_1 import query_extraction

Entrez.email = "your@email.com"
possible_hit_nb = [10, 50, 100, 250, 500, 1000, 5000]


def main(argv):
    """This script will retrieve all the matching sequences for each queries.

    For every given query, it will make a BLAST query (either using tblastn if
    the sequence is coding or blastn otherwise) in order to retrieve every matching
    sequences. The BLAST queries are made with an automated browser which uses the 
    NCBI BLAST portal. The input file should be in the same format as specified in 
    the gitlab: https://gitlab.epfl.ch/baffou/mutational-spectrum . The output files 
    should be in CSV format because we are using DataFrames as output.

    Parameter:

    -input_file: (String) Name of the control file containing the sequences. Please
        look at the gitlab for the file format.

    Keyword Arguments:

    -seq_out: (String) Name of the output file which will contain the sequences. It
        must be a .csv file (default=sequences.csv). 
    -log_out: (String) Name of the output file which will contain the logs. It
        must be a .csv file (default=blast_logs.csv).
    -hits_nb: (int) Max number of blast hits to return per query (default=100).
        It must be one of the following values: [10, 50, 100, 250, 500, 1000, 5000].
    -hit_size_treshold: (int) Min number of hits for the sequences to be kept (default=10).
    -v: (flag) Indicates if the browser should be visible or run in background.
    -h: (flag) Opens this Help section.    

    Outputs:

    Two files named as seq_out and log_out which contains the resulting BLAST hits
    sequences and logs for each query.

    """
    # Parse arguments
    try:
        opts, args = getopt.getopt(
            argv, "hv", ["input_file=", "seq_out=", "log_out=", "hits_nb=", "hit_size_treshold="])
    except getopt.GetoptError:
        print(argv)
        print("Error in inputs, please be sure to provide at least the path to the control file.")
        sys.exit(2)
    input_file = ""
    seq_out = "sequences.csv"
    log_out = "blast_logs.csv"
    headless = True
    hit_index = possible_hit_nb.index(100)
    hit_size_treshold = 10
    for opt, arg in opts:
        if opt == "-h":
            print("help section:\n")
            print("This script will retrieve all the matching sequences for each queries.\n")
            print("For every given query, it will make a BLAST query (either using tblastn if")
            print("the sequence is coding or blastn otherwise) in order to retrieve every matching")
            print("sequences. The BLAST queries are made with an automated browser which uses the") 
            print("NCBI BLAST portal. The input file should be in the same format as specified in") 
            print("the gitlab: https://gitlab.epfl.ch/baffou/mutational-spectrum . The output files") 
            print("should be in CSV format because we are using DataFrames as output.\n")
            print("Parameter:\n")
            print("-input_file: (String) Name of the control file containing the sequences. Please")
            print("    look at the gitlab for the file format.\n")
            print("Keyword Arguments:\n")
            print("-seq_out: (String) Name of the output file which will contain the sequences. It")
            print("    must be a .csv file (default=sequences.csv).") 
            print("-log_out: (String) Name of the output file which will contain the logs. It")
            print("    must be a .csv file (default=blast_logs.csv).")
            print("-hits_nb: (int) Max number of blast hits to return per query (default=100).")
            print("    It must be one of the following values: [10, 50, 100, 250, 500, 1000, 5000].")
            print("-hit_size_treshold: (int) Min number of hits for the sequences to be kept (default=10).")
            print("-v: (flag) Indicates if the browser should be visible or run in background.")
            print("-h: (flag) Opens this Help section.\n")    
            print("Outputs:\n")
            print("Two files named as seq_out and log_out which contains the resulting BLAST hits")
            print("sequences and logs for each query.")
            sys.exit()
        elif opt == "-v":
            headless = False
        elif opt == "--input_file":
            input_file = arg
        elif opt == "--seq_out":
            print(arg)
            seq_out = arg
        elif opt == "--log_out":
            log_out = arg
        elif opt == "--hits_nb":
            if int(arg) not in possible_hit_nb:
                print("The given number of sequences is not in the list of possible choices.",end="") 
                print("Please select among 10, 50, 100, 250, 500, 1000 or 5000.")
                sys.exit(2)
            hit_index = possible_hit_nb.index(int(arg))
        elif opt == "--hit_size_treshold":
            hit_size_treshold = int(arg)
    # Inputs Assertion
    if hit_size_treshold < 1 or hit_size_treshold > possible_hit_nb[hit_index]:
        print("The min number of hits should be in the interval from 1 to hit_nb.")
        sys.exit(2)
    if not os.path.isfile(input_file):
        print("The control file doesn't exist, please provide an existing file")
        sys.exit(2)
    if not seq_out[-4:] == ".csv" or not log_out[-4:] == ".csv":
        print("We use csv as format output, please provide output names with the extension .csv at the end")
        sys.exit(2)
    if os.path.isfile(seq_out):
        print(
            f"The file {seq_out} already exists, please delete it or provide another name")
        sys.exit(2)
    if os.path.isfile(log_out):
        print(
            f"The file {log_out} already exists, please delete it or provide another name")
        sys.exit(2)
    try:
        opening_test = query_extraction(input_file)
    except:
        print("Impossible to parse the control file. Be sure that the format is the same",end="")
        print("as in the github (i.e. tab separated parameters and line separated queries).")
        sys.exit(2)
    # Setup browser
    down_arrow_nb = hit_index-possible_hit_nb.index(100)
    driver = setup_browser(headless)
    # Script execution
    sequences, log = blast_query(input_file, driver, down_arrow_nb, hit_size_treshold)
    pd.DataFrame(sequences).to_csv(seq_out)
    pd.DataFrame(log).T.to_csv(log_out)
    driver.close()
    print("Step 1 executed with success!")
    print(
        f"Output files are : {seq_out} for the sequences and {log_out} for the complementary informations.")


def setup_browser(headless):
    # Setup browser
    chrome_options = Options()
    if headless:
        chrome_options.add_argument("headless")
    chrome_options.add_experimental_option("prefs", {
        "download.default_directory": str(pathlib.Path().resolve()),
        "download.prompt_for_download": False,
        "download.directory_upgrade": True,
        "safebrowsing.enabled": False
    })
    chrome_options.add_argument("--log-level=3")
    driver = webdriver.Chrome(
        ChromeDriverManager().install(), options=chrome_options)
    params = {'behavior': 'allow',
              'downloadPath': str(pathlib.Path().resolve())}
    driver.execute_cdp_cmd('Page.setDownloadBehavior', params)
    return driver


def blast_query(file_name, driver, down_arrow_nb, hit_size_treshold=10):
    """ Make the BLAST queries for each query in the control file.

    It will first parse the control file. Then iteratively make a blast querie,
    parse it and process it, for each initial query sequence. 

    Parameters:
    
    -file_name: (string) Name of the control file containing the sequences.
    -hit_nb: (int) Max number of hits to return.

    Keywords Argument:

    -hit_size_treshold: (int) Min number of hits for the sequences to be kept.

    Outputs:

    -sequence: (DataFrame) Sequences and additional informations resulting 
        from the BLAST hits for each queries.
    -header: (DataFrame) Logs for the different given queries

    """
    queries_df = query_extraction(file_name)
    header = queries_df.to_dict(orient="index")
    queries_parameters = zip(queries_df["CDS"], queries_df["Sequence"],
                            queries_df["Species"], queries_df["Code"],
                            queries_df.index)
    sequences = []
    for param in queries_parameters:
        t_start = time.time()
        print(f"Making query: {param[4]}")
        seq_tuples = []
        if param[0] == "0":
            # Case of proteine coding sequence
            seq_tuples = blastn_query(
                param[1], param[2], param[3], param[4], header, driver, down_arrow_nb)
        else:
            # Case of non-coding sequence
            seq_tuples = tblastn_query(
                param[1], param[2], param[3], param[4], header, driver, down_arrow_nb)
        if len(seq_tuples) >= hit_size_treshold:
            sequences += seq_tuples
        else:
            header.pop(param[4])
        # Use of timers to avoid spamming the server
        t_end = time.time()
        if t_end - t_start < 60:
            time.sleep(60 - (t_end - t_start))
    return sequences, header


def tblastn_query(seq, organism, gen_code, control_id, header_dict, driver, down_arrow_nb):
    """ Make and process a BLAST query over biopython.

    It will make a tblastn query based on the given parameters. It will retrieve
    for each amino acid sequence resulting from hits, the corresponding 
    nucleotide sequence. Then make a quality check and finally outputs a
    list containing all the matching sequences in the good format with 
    additonnal information.

    Parameters:

    -seq: (String) Query sequence.
    -organism: (String) Name of the organism compatible with the NCBI database.
    -gen_code: (int) Id of the genetic coding table to use. 
        The list of the different genetic codes can be found on 
        (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).
    -control_id: (int) Id of the query in the control file.
    -header_dict: (dict) Header used to keep logs.
    -driver: (Selenium driver) Driver used to control the browser.
    -down_arrow_nb: (int) Number of movement to select right number of hits.

    Output:
        
    -valid_sequences: (list) Contains all the matching sequences with extra 
        information.

    """
    valid_sequences = []
    print("Waiting for BLAST output")
    # Nested try to avoid false negative (state impossible when it is not)
    try:
        temp_out_name = make_query("tblastn", seq, organism, driver, down_arrow_nb)
    except OrganismError:
        try:
            print("Organism unfound, trying again.")
            temp_out_name = make_query("tblastn", seq, organism, driver, down_arrow_nb)
        except OrganismError:
            print("The organism provided is not in the NCBI database.")
            return []
        except:
            print("An error occured while making the query, please try again.")
            return []
    except:
        print("An error occured while making the query, please try again.")
        return []
    with open(temp_out_name) as result:
        blast_record = NCBIXML.read(result)
    print("Processing BLAST output")
    codons_sequences = retrieve_codons(blast_record)
    # !!! Use table 1 because cannot specify gen code in browser]we uses the translation to ensure codon alignment
    translated_sequences = list(
        map(lambda a: (a[0], Seq.translate(a[1])), codons_sequences))
    min_start, max_end = find_query_align_position(blast_record)
    blast_sequences = []
    # Retrieve all alignements to compare with the translations
    for alignment in blast_record.alignments:   
        for hsp in alignment.hsps:
            blast_sequences.append(
                (hsp.sbjct, hsp.sbjct_start, hsp.sbjct_end, hsp.query_start, hsp.query_end))
    header_dict[control_id]["Align_Start"] = min_start
    header_dict[control_id]["Align_End"] = max_end
    # Ensure codon alignment, quality check and output format
    for seq_pair in zip(translated_sequences, blast_sequences, codons_sequences):
        if seq_pair[0][1] == seq_pair[1][0]:  
            padded_seq = str(seq_pair[2][1])
            padded_seq = (3*(seq_pair[1][3]-min_start))*"-" + \
                padded_seq+(3*(max_end - seq_pair[1][4]))*"-"
            seq_dict = {"Control_ID": control_id,
                        "Seq_ID": seq_pair[2][0], "Sequence": padded_seq, "Start_Pos": seq_pair[1][1], 
                        "End_Pos": seq_pair[1][2]}
            valid_sequences.append(seq_dict)
    header_dict[control_id]["Hits_Number"] = len(valid_sequences)
    print("Query process finished. Removing temp file.")
    if os.path.exists(temp_out_name):
        os.remove(temp_out_name)
        print("Temp file removed.\n")
    else:
        print("The temp file doesn't exist, script output potentially corrupted. Abortion.")
        sys.exit(2)
    return valid_sequences


def blastn_query(seq, organism, gen_code, control_id, header_dict, driver, down_arrow_nb):
    """ Make and process a BLAST query over biopython.

    It will make a blastn query based on the given parameters. It will retrieve
    for each amino acid sequence resulting from hits, the corresponding 
    nucleotide sequence. Then make a quality check and finally outputs a
    list containing all the matching sequences in the good format with 
    additonnal information.

    Parameters:

    -seq: (String) Query sequence.
    -organism: (String) Name of the organism compatible with the NCBI database.
    -gen_code: (int) Id of the genetic coding table to use. 
        The list of the different genetic codes can be found on 
        (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).
    -control_id: (int) Id of the query in the control file.
    -header_dict: (dict) Header used to keep logs.
    -driver: (Selenium driver) Driver used to control the browser.
    -down_arrow_nb: (int) Number of movement to select right number of hits.

    Output:
        
    -blast_sequences: (list) Contains all the matching sequences with extra 
        information.

    """
    print("Waiting for BLAST output")
    # Nested try to avoid false negative (state impossible when it is not)
    try:
        temp_out_name = make_query("blastn", seq, organism, driver, down_arrow_nb)
    except OrganismError:
        try:
            temp_out_name = make_query("blastn", seq, organism, driver, down_arrow_nb)
        except OrganismError:
            print("The organism provided is not in the NCBI database.")
            return []
        except:
            print("An error occured while making the query, please try again.")
            return []
    except:
        print("An error occured while making the query, please try again.")
        return []
    with open(temp_out_name) as result:
        blast_record = NCBIXML.read(result)
    print("Processing BLAST output")
    blast_sequences = []
    min_start, max_end = find_query_align_position(blast_record)
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            new_sequence = delete_insertions(hsp, seq)
            new_sequence = add_gaps(min_start, max_end, new_sequence, hsp)
            seq_dict = {"Control_ID": control_id, "Seq_ID": alignment.accession,
                        "Sequence": new_sequence, "Start_pos": hsp.sbjct_start, "End_pos": hsp.sbjct_end}
            blast_sequences.append(seq_dict)
    header_dict[control_id]["Align_Start"] = min_start
    header_dict[control_id]["Align_End"] = max_end
    header_dict[control_id]["Hits_Number"] = len(blast_sequences)
    print("Query process finished. Removing temp file.")
    if os.path.exists(temp_out_name):
        os.remove(temp_out_name)
        print("Temp file removed.")
    else:
        print("The temp file doesn't exist, script output potentially corrupted. Abortion.")
        sys.exit(2)
    return blast_sequences


def retrieve_codons(record):
    """ Given a blast record, retrieve every aligned codons sequences."""
    max_entrez_query_nb = 200  # Entrez.fetch() doesn't accept more than 200 inputs
    sequences = []
    id_list = []
    for alignment in record.alignments:
        for hsp in alignment.hsps:
            id_list.append(
                (alignment.accession, hsp.sbjct_start, hsp.sbjct_end))
    # Split the main list into sublists to avoid big chunk of ids
    id_list_split = [id_list[x:x+max_entrez_query_nb]
                     for x in range(0, len(id_list), max_entrez_query_nb)]
    for sub_id_list in id_list_split:
        id_string = ""
        for i in sub_id_list:
            id_string += i[0]+","
        id_string = id_string[:-1]
        # Retrieve the nucleotides sequence
        handle = Entrez.efetch(db="nucleotide", id=id_string,
                               rettype="fasta", retmode="text")
        fasta_file = SeqIO.parse(handle, format="fasta")
        for seq in zip(fasta_file, sub_id_list):
            start = seq[1][1] - 1  # Include the fist nucleotide
            stop = seq[1][2]
            name = seq[1][0]
            sequences.append((name, seq[0].seq[start:stop]))
    return sequences


def find_query_align_position(record):
    """Find the lowest start point and the largest end point of alignments.""" 
    min_query_start = np.inf
    max_query_end = -1
    for alignment in record.alignments:
        for hsp in alignment.hsps:
            if hsp.query_start < min_query_start:
                min_query_start = hsp.query_start
            if hsp.query_end > max_query_end:
                max_query_end = hsp.query_end
    return min_query_start, max_query_end


def add_gaps(min_start, max_end, seq, hsp):
    """Pad the sequence with gaps on the left and right to have fixed length sequences"""
    return "-"*(hsp.query_start-min_start) + seq + "-"*(max_end-hsp.query_end)


def delete_insertions(hsp, original_sequence):
    """Remove every insertion from the BLAST output sequence"""
    insertion_indices = search_gap(hsp.query)
    original_gaps = set(search_gap(original_sequence))
    new_subjct = hsp.sbjct
    indice_corrector = 0
    for i in insertion_indices:
        if i not in original_gaps:
            new_subjct = (new_subjct[:i-indice_corrector]
                + new_subjct[i+1-indice_corrector:])
            # Need a correction as we are deleting elements, thus shifting indices
            indice_corrector += 1
    return new_subjct


def search_gap(seq):
    """Find every gaps and output their indices."""
    search = seq.find("-")
    summer = 0
    insertion_indices = []
    while search != -1:
        summer += search
        insertion_indices.append(summer)
        summer += 1  # Avoid to stay on a gap forever
        search = seq[summer:].find("-")
    return insertion_indices


def make_query(program, sequence, organism, driver, down_arrow_nb):
    """ Make a BLAST query using an automated browser.
    
    It will first select a program between tblastn and blastn. Then make a
    query using the NCBI BLAST portal and the different parameters provided.
    It will finaly download all the sequences resulting from BLAST hits into
    a XML file. If the organism is not in the NCBI database, it will raise
    an Exception.

    Parameters:

    -program: (String) Blast program to use between tblastn and blastn.
    -sequence: (String) Original query sequence.
    -organism: (String) Name of the specie.
    -driver: (Selenium driver) Driver used to control the browser.
    -down_arrow_nb: (int) Number of movement to select right number of hits.

    Output:

    XML file containing the sequences

    """
    url = ""
    name_index = 0
    if program == "tblastn":
        url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=tblastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome"
        name_index = 5
    else:
        url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome"
        name_index = 6
    # Enter query parameter
    driver.get(url)
    time.sleep(1)
    seq_box = driver.find_element_by_id("seq")
    seq_box.clear()
    seq_box.send_keys(sequence)
    org_box = driver.find_element_by_id("qorganism")
    org_box.clear()
    org_box.send_keys(organism)
    time.sleep(1.5)
    org_box.send_keys(Keys.ARROW_DOWN)
    org_box.send_keys(Keys.ENTER)
    opt_box = driver.find_element_by_id("btnDescrOver")
    opt_box.click()
    hit_nb_box = driver.find_element_by_id("NUM_SEQ")
    actual_hit_index = possible_hit_nb.index(
        int(hit_nb_box.get_attribute("value"))) - possible_hit_nb.index(100)
    if actual_hit_index != down_arrow_nb:
        if down_arrow_nb < 0:
            for i in range(-down_arrow_nb):
                hit_nb_box.send_keys(Keys.ARROW_UP)
        else:
            for i in range(down_arrow_nb):
                hit_nb_box.send_keys(Keys.ARROW_DOWN)
    blast_submit = driver.find_elements_by_id("blastButton1")
    blast_submit[0].click()
    # Wait until response
    waiting = True
    while waiting:
        try:
            align_btn = driver.find_element_by_id("btnAlign")
        except:
            time.sleep(10)
            try:
                org_error_box = driver.find_element_by_id("lpgMsg")
            except:
                pass
            else:
                if org_error_box.is_displayed():
                    raise OrganismError
        else:
            waiting = False
    out_name_box = driver.find_element_by_id("brc")
    # Use an offset (4) to select the RID of the BLAST job
    out_name = out_name_box.text.split()[name_index][4:]+"-Alignment.xml"
    align_btn = driver.find_element_by_id("btnAlign")
    align_btn.click()
    time.sleep(1)
    align_select = driver.find_element_by_id("alignViewSelect")
    align_select.send_keys("QueryAnchoredNoIdentities")
    # Wait for alignments
    waiting = True
    while waiting:
        try:
            align_ready = driver.find_element_by_id("alnDwnld_1")
        except:
            time.sleep(2)
        else:
            waiting = False
    # Download output
    download_btn = driver.find_element_by_id("btnDwnldAln")
    download_btn.click()
    time.sleep(1)
    download_xml = driver.find_element_by_id("dwXMLAln")
    download_xml.click()
    # wait to avoid closing browser before having downloaded file
    time.sleep(3)
    return out_name


class OrganismError(Exception):
    """Raised when the organism is not found in the NCBI database."""
    pass


if __name__ == "__main__":
    main(sys.argv[1:])
