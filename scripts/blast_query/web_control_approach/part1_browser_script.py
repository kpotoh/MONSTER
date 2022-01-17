import sys
import getopt
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Seq
from Bio import SeqIO
from Bio import Entrez
Entrez.email = "your@email.com"
import numpy as np
import pandas as pd
import os
import selenium
from selenium import webdriver
from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.keys import Keys
import pathlib
import time

possible_hit_nb = [10,50,100,250,500,1000,5000]

def main(argv):
    # Parse arguments
    try:
        opts, args = getopt.getopt(argv,"hv",["input_file=","seq_out=","log_out=","seq_nb="])
    except getopt.GetoptError:
        print(argv)
        print("Error in inputs, please be sure to provide at least the path to the control file.")
        sys.exit(2)
    input_file = ""
    seq_out = "sequences.csv"
    log_out = "blast_logs.csv"
    headless = False
    hit_index = possible_hit_nb.index(100)
    for opt, arg in opts:
        if opt == "-h":
            print("help section:\n")
            print("This script will make blast queries using the Biopython module and the NCBI API.\n\n")
            print("Inputs:")
            print("-input_file : mandatory parameters indicating the Control File")
            print("-seq_out: name of the output file for the sequences. It will be a csv file so don't forget the .csv extension. (Default : sequences.csv)")
            print("-log_out: name of the output file for the log corresponding to the sequences. It will be a csb file so don't forget the .csv extension. (Default: blast_log.csv)\n\n")
            print("Outputs:")
            print("The script will output two files. One containing the sequences resulting from blast hits. The other contains logs about the queries.")
            print("")
            sys.exit()
        elif opt == "-v":
            headless = True
        elif opt == "--input_file":
            input_file = arg
        elif opt == "--seq_out":
            print(arg)
            seq_out = arg
        elif opt == "--log_out":
            log_out = arg
        elif opt == "--seq_nb":
            if int(arg) not in possible_hit_nb:
                print("The given number of sequences is not in the list of possible choices. Please select among 10, 50, 100, 250, 500, 1000 or 5000.")
                sys.exit(2)
            hit_index = possible_hit_nb.index(int(arg))
    # Inputs Assertion
    if not os.path.isfile(input_file):
        print("The control file doesn't exist, please provide an existing file")
        sys.exit(2)
    if not seq_out[-4:] == ".csv" or not log_out[-4:] == ".csv":
        print("We use csv as format output, please provide output names with the extension .csv at the end")
        sys.exit(2)
    if os.path.isfile(seq_out):
        print(f"The file {seq_out} already exists, please delete it or provide another name")
        sys.exit(2)
    if os.path.isfile(log_out):
        print(f"The file {log_out} already exists, please delete it or provide another name")
        sys.exit(2)
    try:
        opening_test = query_extraction(input_file)
    except :
        print("Impossible to parse the control file. Be sure that the format is the same as in the github (i.e. tab separated parameters and line separated queries).")
        sys.exit(2)
    # Script execution
    down_arrow_nb = hit_index-possible_hit_nb.index(100)
    # Setup browser
    driver = setup_browser(headless)
    sequences, log = blast_query(input_file,driver,down_arrow_nb)
    pd.DataFrame(sequences).to_csv(seq_out) 
    pd.DataFrame(log).T.to_csv(log_out)
    driver.close()
    print(f"Script executed with success! Output files are : {seq_out} for the sequences and {log_out} for the complementary informations")

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
    driver = webdriver.Chrome(ChromeDriverManager().install(),options=chrome_options)
    params = {'behavior' : 'allow', 'downloadPath':str(pathlib.Path().resolve())}
    driver.execute_cdp_cmd('Page.setDownloadBehavior', params)
    return driver

def full_script(file_name,seq_output_name="sequences.csv", log_output_name="blast_log.csv"):
    """
    Given a control file, will output matching nucleotide sequences with BLAST for every given sequences.
    It will assert the inputs, make the blast queries, format and save the output.
    
    Parameters:
        file_name : name(string) of the control file
        seq_output_name : name (with csv extension) of the sequence output file (string, default value : sequences.csv)
        log_output_name : name (with csv extension) of the extra information output file (string, default value : blast_log.csv)
    
    Returns:
        None (but produce two files: the alignments and a log file)
    """
    if not os.path.isfile(file_name):
        print("The control file doesn't exist, please provide an existing file")
        sys.exit(2)
    if not seq_output_name[-4:] == ".csv" or not log_output_name[-4:] == ".csv":
        print("We use csv as format output, please provide output names with the extension .csv at the end")
        sys.exit(2)
    if os.path.isfile(seq_output_name):
        print(f"The file {seq_output_name} already exists, please delete it or provide another name")
        sys.exit(2)
    if os.path.isfile(log_output_name):
        print(f"The file {log_output_name} already exists, please delete it or provide another name")
        sys.exit(2)
    try:
        opening_test = query_extraction(file_name)
    except :
        print("Impossible to parse the control file. Be sure that the format is the same as in the github (i.e. tab separated parameters and line separated queries).")
        sys.exit(2)
    sequences, log = blast_query(file_name)
    pd.DataFrame(sequences).to_csv(seq_output_name) 
    pd.DataFrame(log).T.to_csv(log_output_name)
    print(f"Script executed with success! Output files are : {seq_output_name} for the sequences and {log_output_name} for the complementary informations")
    
def query_extraction(file_name):
    """
    Given a control file with tab separated parameters, it will produce a pandas DataFrame with all queries informations.
    
    Parameter:
        file_name : string that indicates the name of the control file
        
    Return:
        query_df : pandas Dataframe with queries informations
    """
    with open(file_name) as f:
        raw_lines = f.readlines()
    queries = list(map(lambda l : l[:-1].split("\t"),raw_lines)) # parsing the control file
    query_df = pd.DataFrame(data = queries[1:], columns=queries[0])
    query_df["Code"] = pd.to_numeric(query_df["Code"],downcast="integer")
    return query_df

def blast_query(file_name,driver,down_arrow_nb):
    """
    Given a control file, will output matching nucleotide sequences with BLAST for every given sequences
    
    Parameter:
        file_name : string indicating the control file name
        
    Return:
        A list of tuple which contains a query identificator (the index of the query), a nucleotide sequence and an alignement id
    """
    queries_df = query_extraction(file_name)
    header = queries_df.to_dict(orient="index")
    queries_parameters = zip(queries_df["CDS"],queries_df["Sequence"],queries_df["Species"],queries_df["Code"],queries_df.index)
    sequences = []
    for param in queries_parameters:
        print(f"Making query: {param[4]}")
        seq_tuples = []
        if param[0] == "0":
            seq_tuples = blastn_query(param[1],param[2],param[3],param[4],header,driver,down_arrow_nb) # Case of proteine coding sequence
        else:
            seq_tuples = tblastn_query(param[1],param[2],param[3],param[4],header,driver,down_arrow_nb) # Case of non-coding sequence
        sequences += seq_tuples
    return sequences, header

def tblastn_query(seq,organism,gen_code,control_id,header_dict,driver,down_arrow_nb):
    """
    Given parameters, will make a BLAST query for coding sequences and output the codon aligned nucleotide sequence
    together with the id of the sequence.
    
    Parameter:
        parameters : tuple with the sequence (as as string) in first position and the organism in second position
        
    Return:
        A list of tuple which contains every nucleotides sequences for a given sequence with the alignement id
    """
    valid_sequences =[]
    print("Waiting for BLAST output")
    temp_out_name = make_query("tblastn",seq,organism,driver,down_arrow_nb)
    with open(temp_out_name) as result:
        blast_record = NCBIXML.read(result)
    print("Processing BLAST output")
    codons_sequences = retrieve_codons(blast_record)
    translated_sequences = list(map(lambda a: (a[0],Seq.translate(a[1])),codons_sequences)) # !!! [use table 1 because cannot specify gen code in browser]we uses the translation to ensure codon alignment
    min_start, max_end = find_query_align_position(blast_record) # we need to know the overall size of the query sequence to fullfill with gaps
    ###
    blast_sequences = []
    for alignment in blast_record.alignments:   # retrieve all alignements to compare with the translations
        for hsp in alignment.hsps:
            blast_sequences.append((hsp.sbjct,hsp.sbjct_start,hsp.sbjct_end,hsp.query_start,hsp.query_end))
    header_dict[control_id]["align_start"] = min_start
    header_dict[control_id]["align_end"] = max_end 
    ###
    for seq_pair in zip(translated_sequences,blast_sequences,codons_sequences):
        if seq_pair[0][1] == seq_pair[1][0]: # ensure codon alignment and quality check
            padded_seq = str(seq_pair[2][1])
            padded_seq = (3*(seq_pair[1][3]-min_start))*"-"+padded_seq+(3*(max_end - seq_pair[1][4]))*"-"
            seq_dict = {"Control_ID" : control_id ,"Seq_ID" : seq_pair[2][0], "Sequence" : padded_seq, "Start_Pos" : seq_pair[1][1], "End_Pos" : seq_pair[1][2]}
            valid_sequences.append(seq_dict) 
    print("Query process finished. Removing temp file.")
    if os.path.exists(temp_out_name):
        os.remove(temp_out_name)
        print("Temp file removed, going to next query.")
    else:
        print("The temp file doesn't exist, script output potentially corrupted. Abortion.")
        sys.exit(2) 
    return valid_sequences
    
def blastn_query(seq,organism,gen_code,control_id,header_dict,driver,down_arrow_nb):
    """
    Given parameters, will make a BLAST query for non coding sequences and output the non aligned nucleotide sequence
    together with the id of the sequence.
    
    Parameter:
        parameters : tuple with the sequence (as as string) in first position and the organism in second position
        
    Return:
        A list of tuple which contains every nucleotides sequences for a given sequence with the alignement id
    """
    print("Waiting for BLAST output")
    temp_out_name = make_query("blastn",seq,organism,driver,down_arrow_nb)
    with open(temp_out_name) as result:
        blast_record = NCBIXML.read(result)
    print("Processing BLAST output")
    blast_sequences = []
    min_start, max_end = find_query_align_position(blast_record) # we need to know the overall size of the query sequence to fullfill with gaps
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            new_sequence = delete_insertions(hsp,seq)
            new_sequence = add_gaps(min_start,max_end,new_sequence,hsp)
            seq_dict = {"Control_ID" : control_id, "Seq_ID" : alignment.accession, "Sequence" : new_sequence, "Start_pos" : hsp.sbjct_start, "End_pos" : hsp.sbjct_end}
            blast_sequences.append(seq_dict)
    ##
    header_dict[control_id]["align_start"] = min_start
    header_dict[control_id]["align_end"] = max_end
    ##
    print("Query process finished. Removing temp file.")
    if os.path.exists(temp_out_name):
        os.remove(temp_out_name)
        print("Temp file removed, going to next query.")
    else:
        print("The temp file doesn't exist, script output potentially corrupted. Abortion.")
        sys.exit(2) 
    return blast_sequences

def retrieve_codons(record):
    """ 
    Given a blast record, will output every aligned codons sequences
    
    Parameter:
        record : blast XML record (either tblastn or blastn)
    
    Return : 
        A list of tuple which contains every HSP codons sequences together with the alignement id
    """
    max_entrez_query_nb = 200 #The Entrez.fetch function doesn't accept more than 200 input at a time
    sequences = []
    id_list = []
    for alignment in record.alignments:
        for hsp in alignment.hsps:
            id_list.append((alignment.accession,hsp.sbjct_start,hsp.sbjct_end)) # create a tuple with sequence id and start, stop ids in the sequence
    ###
    id_list_split = [id_list[x:x+max_entrez_query_nb] for x in range(0, len(id_list), max_entrez_query_nb)]# split the main list into sublists to avoid big chunk of ids
    for sub_id_list in id_list_split:
        id_string = ""
        for i in sub_id_list:
            id_string+=i[0]+","
        id_string = id_string[:-1]
        handle = Entrez.efetch(db="nucleotide", id=id_string, rettype="fasta", retmode="text") #we retrieve the nucleotides sequence
        fasta_file = SeqIO.parse(handle,format="fasta")
        for seq in zip(fasta_file,sub_id_list):
            start = seq[1][1] - 1 # the fist nucleotide is contained
            stop = seq[1][2] 
            name = seq[1][0]
            sequences.append((name,seq[0].seq[start:stop]))
    return sequences

def find_query_align_position(record):
    """
    Given a BLAST record, will find the lowest starting point and the largest end point  of alignments in the query sequence
    
    Parameters:
        record : Bio.Blast.Record whose alignment points need to be find
        
    Return:
        min_query_start : int indicating the lowest index of starting point
        max_query_end : int indicating the largest index of ending point

    """
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
    """
    Full fill a given sequence with gaps on the left and right to obtain fixed length sequences
    
    Parameters:
        min_start: index(int) of the lowest query starting point
        max_end: index(int) of the largest query ending point
        seq: BLAST output sequences (string)
        hsp: HSP object in Bio.Blast.Record designing the output of the BLAST record for this sequence
        
    Return:
        output : fixed length sequence fullfilled with gaps
    """
    return "-"*(hsp.query_start-min_start) + seq + "-"*(max_end-hsp.query_end)

def delete_insertions(hsp,original_sequence):
    """
    Given a BLAST alignment and the original sequence, retrieve the BLAST output sequence without any insertions
    
    Parameters: 
        hsp : HSP object in Bio.Blast.Record designing the output of the BLAST record for a sequence
        original_sequence : original query sequence (string) for the hsp
        
    Return:
        new_subjct : BLAST output sequence (string) without any insertion
    """
    insertion_indices = search_gap(hsp.query)
    original_gaps = set(search_gap(original_sequence)) #cast in set for faster index research (as they are unique)
    new_subjct = hsp.sbjct
    indice_corrector = 0
    for i in insertion_indices:
        if i not in original_gaps:
            new_subjct = new_subjct[:i-indice_corrector]+new_subjct[i+1-indice_corrector:]
            indice_corrector += 1 # we need a correction as we are deleting elements, thus shifting indices
    return new_subjct

def search_gap(seq):
    """
    Given a sequence, it will find every gaps and output their indices.
    
    Parameters: 
        seq: sequence of interest(string)
    
    Return:
        insertion_indices: list of indices(int) where gaps have been found
    """
    search = seq.find("-")
    summer = 0
    insertion_indices = []
    while search != -1:
        summer += search
        insertion_indices.append(summer)
        summer += 1 # used to avoid to stay on a gap forever
        search = seq[summer:].find("-")
    return insertion_indices

def make_query(program, sequence, organism, driver, down_arrow_nb):
    url = ""
    name_index = 0
    if program == "tblastn":
        url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=tblastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome"
        name_index = 5
    else :
        url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome"
        name_index = 6
    ### make query
    driver.get(url)
    time.sleep(1)
    seq_box = driver.find_element_by_id("seq")
    seq_box.clear()
    seq_box.send_keys(sequence)
    org_box = driver.find_element_by_id("qorganism")
    org_box.clear()
    org_box.send_keys(organism)
    time.sleep(1)
    org_box.send_keys(Keys.ARROW_DOWN) # need to validate input with the tax_id
    org_box.send_keys(Keys.ENTER)
    opt_box = driver.find_element_by_id("btnDescrOver")
    opt_box.click()
    hit_nb_box = driver.find_element_by_id("NUM_SEQ")
    if down_arrow_nb < 0:
        for i in range(-down_arrow_nb):
            hit_nb_box.send_keys(Keys.ARROW_UP)
    else:
        for i in range(down_arrow_nb):
            hit_nb_box.send_keys(Keys.ARROW_DOWN)
    blast_submit = driver.find_elements_by_id("blastButton1")
    blast_submit[0].click()
    ### wait until response
    waiting = True
    while waiting:
        try:
            align_btn = driver.find_element_by_id("btnAlign")
        except:
            time.sleep(10)
        else: 
            waiting = False
    ### download output
    out_name_box = driver.find_element_by_id("brc")
    out_name = out_name_box.text.split()[name_index][4:]+"-Alignment.xml" # we have an offset to select the RID section and then cut the string "RID-"
    print(f"The outname is {out_name}\n")
    align_btn = driver.find_element_by_id("btnAlign")
    align_btn.click()
    time.sleep(1)
    align_select = driver.find_element_by_id("alignViewSelect")
    align_select.send_keys("QueryAnchoredNoIdentities")
    # need to wait while producing the alignments
    waiting = True
    while waiting:
        try:
            align_ready = driver.find_element_by_id("alnDwnld_1")
        except:
            time.sleep(2)
        else: 
            waiting = False
    download_btn = driver.find_element_by_id("btnDwnldAln")
    download_btn.click()
    time.sleep(1)
    download_xml = driver.find_element_by_id("dwXMLAln")
    download_xml.click()
    time.sleep(3) ## wait to avoid closing browser before having downloaded the entire file
    #driver.close()
    return out_name

if __name__ == "__main__":
    main(sys.argv[1:])

