# Mutational Spectrum

The mutational spectrum is a great tool to analyse the structure of mutations. With this project, we try to provide a fast and simple way to create this mutational spectrum. Provided a number (small or large) of sequences into a control file, the script will make BLAST queries in order to find every matching sequences and from there derive the mutational spectrum. 

# Table of Contents
- [Mutational Spectrum](#mutational-spectrum)
- [Table of Contents](#table-of-contents)
- [Environment](#environment)
- [Pipeline Description](#pipeline-description)
  - [Main characteristics of the pipeline:](#main-characteristics-of-the-pipeline)
  - [Pipeline Structure](#pipeline-structure)
    - [Step 1)](#step-1)
    - [Step 2)](#step-2)
    - [Step 3)](#step-3)
- [Requirements](#requirements)
- [Pipeline Parameters](#pipeline-parameters)
- [Input/Output](#inputoutput)
  - [Step 1](#step-1-1)
    - [Input](#input)
    - [Columns Description:](#columns-description)
    - [Input Example](#input-example)
    - [Script Parameters](#script-parameters)
    - [Output](#output)
      - [Sequences output](#sequences-output)
      - [Log output](#log-output)
  - [Step 2](#step-2-1)
    - [Input](#input-1)
    - [Script Parameters](#script-parameters-1)
    - [Output](#output-1)
  - [Step 3](#step-3-1)
    - [Input](#input-2)
    - [Script Parameters:](#script-parameters-2)
    - [Output](#output-2)
- [Automated browser vs Biopython](#automated-browser-vs-biopython)
  - [Browser](#browser)
  - [Biopython](#biopython)
  - [Conclusion](#conclusion)
- [Limitations](#limitations)

# Environment

**Requirements**:
- unix
- installed Chrome bwowser
- python 3.8+

**Setup enviroment**:
```
pip install -r requirements.txt
```

**Tests & Dev**:
```
pip install -r requirements.dev.txt  # only for development
pytest scripts/pipeline/test_pipe.py  # project structure will be modified soon
```


# Pipeline Description


## Main characteristics of the pipeline

* Easy to use python-based pipeline (no supercomputers and unix) for all OSs.
* Phylogeny-free approach allows to work with all species where consensus is used as if ancestral sequence.
* Modular structure: there is a possibility to use any steps of choice with standard inputs/outputs.
* All different genetic codes are available.
* An internal normalization on species-specific nucleotide/codon contents (expected spectra) allows to derive a mutational spectrum, comparable between different species.
* Query-anchored alignment and results give a possibility to use any additional query-specific annotations (RefSeq etc).
* Reports of Flanking nucleotides xxAxx > xxTxx give a possibility to reconstruct deep mutational spectrum with context (1, 3 or 5 bases).
* Flexible threshold of the minor allele frequencies (with default 10%) used to reconstruct spectra. 
* Compatibility of outputs on each step with R and popular R packages.  

## Pipeline Structure

We can split the pipeline in three main steps:  

- Step 1: Parse the control file and make the BLAST queries. Format the output as query anchored alignments (codon alignment in the case of coding sequences).  

- Step 2: Derive consensus and observed deviations from the previous alignments output. In the case of coding sequences, also create the expected deviations, which will be useful for Step 3 when we will create the mutational spectrum.    

- Step 3: Create the Mutational Spectrum.  

Each step is executed with a different script which can be run in standalone mode (provided that you feed them with inputs in correct format). The outputs of each steps are saved, thus you can access the intermediate results/files to play with them or if you need extra informations on your queries.  

### Step 1)

The BLAST queries (either tblastn if coding sequence or blastn otherwise) are made via the BLAST module of the Biopython library by default. You can specify if you prefer to use a browser automated version, c.f below for [pros and cons](#automated-browser-vs-biopython). Once the queries done, we retrieve the alignments and create a pandas DataFrame containing the aligned sequences together with some informations (such as start and end point) for each query. A log file containing each query sequences where BLAST hits has been found is also created. It contains also the alignment positions on the original queries.

### Step 2)

The consensus is derived from the aligned sequences by taking for each position the nucleotide corresponding to the major allele frequency. The observed deviations are every deviations from a sequence compared to the consensus (and **not** the original query). The expected deviations are all possible SNP from consensus but they are only relevant with coding sequences.

### Step 3)

The deviations are first filtered, e.g. in coding case we keep only synonymous mutations, then normalized. The normalization is different depending if the original query sequences is coding or not. If it is coding, we normalize by the number of mutations of the same type in the expected deviations. If it is non-coding, we normalize by the percentage of occurences of the ancestor of the nucleotide that changes. Then we make a normalization over all these number to have frequencies.

# Requirements

The code is in Python 3.8.5 and uses four external libraries:

- numpy (v1.19.2)
- pandas (v1.1.3)
- biopython (v1.78)
- selenium (v3.141.0)

You can use the following command to install the requirements:

`pip install -r requirements.txt`


# Pipeline Parameters

The pipeline takes as parameters:  

- input_file: (String) Name of the control file. The control file should be in the format depicted in the previous section.
- out_folder: (String) Name of the folder where both intermediate and final outputs will be stored (default=outputs).
- hits_nb: (int) Max number of blast hits to return per query (default=100). It must be one of the following values: [10, 50, 100, 250,500, 1000, 5000].
- neighbours_nb: (int) Number of desired neighbours (on the left and on the right) in the output (default=2).
- hit_size_treshold: (int) Min number of hits for the sequences to be kept (default=10).
- codon: (flag) Indicate if codons should be provided in the output (default=False)
- b: (flag) Indicates if the step 1 should run in browser mode rather than with biopython (default=False).
- v: (flag) Indicates if the browser should be visible or be run in background (default=False).   
- c: (flag) Indicate if context should be taken into account for the mutational spectrum (default=False).
- h: (flag) Opens the Help section.

# Input/Output

## Step 1

This script has the particularity to have two possible variants. You can choose between using biopython to make BLAST queries or to use a browser automated approach. The [pros and cons](#automated-browser-vs-biopython) for both methods are listed at the end of the README if you are interested. The only difference between the two option is how to make the query. Once the file obtained, the two approaches use the same query parsing and processing steps.
 
### Input

The input sequences must be in the form of a control file. The control file have the following format:  

Species,Gene,Code,CDS,Sequence  
Query_1_Specie,Query_1_Gene,Query_1_Code,Query_1_CDS,Query_1_Sequence  
Query_2_Specie,Query_2_Gene,Query_2_Code,Query_2_CDS,Query_2_Sequence    

So the format is **comma delimited parameters** and **line separated queries**. Don't forget to put the first line with the name of the columns (comma separated).

### Columns Description:

- Species: (String) Name of the specie
- Gene: (String) Name of the gene
- Seq: (String) Sequence (either amino acids or nucleotide)
- Code: (int) Genetic code to use (1 is standard, 2 is vertebrate mithocondrial, etc. Check https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes#SG2 for the list of genetic codes id).
- CDS: (int) Indicate if the sequence is a coding sequence (in which case a 1 should be written and tblast will be used) or a non coding sequence (in which case a 0 should be written and blastn will be used).

### Input Example

Species,Gene,Code,CDS,Sequence   
Mus Musculus domesticus,NADH1,2,1,MPMANLLLLIVPILIAMAFLMLTERKILGYMQLRKGPNVVGPYGLLQPFADAMKLFTKEPLKPATSTITLYITAPTLALTIALLLWTPLPMPNPLVNLNLGLLFILATSSLAVYSILWSGWASNSNYALIGALRAVAQTISYEVTLAIILLSTLLMSGSFNLSTLITTQEHLWLLLPSWPLAMMWFISTLAETNRTPFDLAEGESELVSGFNIEYAAGPFALFFMAEYTNIIMMNTLTTTIFLGTTYDALSPELYTTYFVTKTLLLTSLFLWIRTAYPRFRYDQLMHLLWKNFLPLTLALLMWYVSMPITISSIPPQT   
Pan troglodytes,D-loop,2,0,tanccgntatgtatttcgtacattactgccagccaccatgaatattatatagtactataacacttaaccacctataacacataaaaacctacatccacattaaaactttaccccatgcttacaagcacgaacaataatcaacctccaactgtcaaacataacacacaactccaaagacactcctcccccaccccgatatcaacaaacctgacaatccttgatagtacatagtacatacagtcatataccgtacatagcacattacagtcaaatccattctcgcccccacggatgccccccctcnnnnagggg      
  

### Script Parameters

The script takes as parameters:  

- input_file: (String) Name of the control file containing the sequences.
- seq_out: (String) Name of the output file which will contain the sequences. It must be a .csv file (default=sequences.csv). 
- log_out: (String) Name of the output file which will contain the logs. It must be a .csv file (default=blast_logs.csv).
- hits_nb: (int) Max number of blast hits to return per query (default=100). It must be one of the following values: [10, 50, 100, 250, 500, 1000, 5000].
- hit_size_treshold: (int) Min number of hits for the sequences to be kept (default=10).
- h: (flag) Opens the Help section.
  
### Output  
  

The Step 1 outputs two files. One contains the sequences resulting from BLAST hits, the other is a log file that provide informations about initial queries. Both are pandas DataFrame saved as .csv files.

#### Sequences output
  
The file has five columns:

- Control_ID: (int) Indicate the id of the initial query. The id is the index in the log file (which is the same as in the Control file).  
- Seq_ID: (String) Indicate the id of the BLAST hit (in the blast sequence format). You can make a BLAST search based on this id to find extra information about the sequence.  
- Sequence: (String) Resulting sequence. Padding will be used so that each sequence has the same length.
- Start_Pos: (int) Position of the first nucleotide in the hit sequence, where the alignment start.
- Start_Pos: (int) Position of the last nucleotide in the hit sequence, where the alignment end.
  

#### Log output  
  
The file has seven columns, where the five first are identical to the Control file:

- Species: (String) Name of the specie.
- Gene: (String) Name of the gene.
- Seq: (String) Sequence (either amino acids or nucleotide).
- Code: (int) Genetic code to use for translation (1 is standard, 2 is vertebrate mithocondrial, etc. Check https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes#SG2 for the list of genetic codes id).
- CDS: (int) Indicate if the sequence is a coding sequence (in which case a 1 should be put and tblast will be used) or a non coding sequence (in which case a 0 should be put and blastn will be used).
- Align_Start: (int) Position of the first nucleotide in the query sequence, where the alignment start.
- Align_End" (int) Position of the last nucleotide in the query sequence, where the alignment end.
  


## Step 2
  

### Input

The Step 2 takes as mandatory input the two output files of the [Step 1](#step-1-1). Thus, look at the output of the previous step to see the input format for this one.  
  
### Script Parameters
  
- seq_input: (String) Name of the input file containing the sequences for each queries as specified in the output of Step 1. 
- log_input: (String) Name of the input file containing the log for each queries as specified in the output of Step 1. 
- output: (String) Name of the output file containing the deviations (default="deviations.json").
- neighbours_nb: (int) Number of desired neighbours (on the left and on the right) in the output (default : 2).
- codon: (flag) Indicate if codons should be provided in the output (default : False)
- h: (flag) Opens the Help section.

### Output

The Step 2 outputs a single file containing all the deviations (both observed and expected in the coding case). The file is a JSON with each index corresponding to an initial query and each value is a dictionnary with 2 to 3 fields depending if the original sequence was coding or not:  
- Consensus: (String) Consensus sequence for the given query.
- Observed_Deviations: (dict) Observed deviations represented as dictionnaries. 
- (Expected_Deviations: (dict) Expected deviations represented as dictionnaries. Present only in coding cases).
  
The deviations dictionnaries take the following form for the non-coding case:  
- Observed Deviations:  {Seq_ID,Codon_Pos,Nucleotide_Pos,From,To,MAF}.  

And the following one for the coding case:    
- Observed Deviations: {Seq_ID,Codon_Pos,Nucleotide_Pos,From,To,MAF,Amino_Acid_From,Amino_Acid_To,Synonymous}.
- Expected Deviations: {Codon_Pos,Nucleotide_Pos,From,To,Amino_Acid_From,Amino_Acid_To,Synonymous}.
  
The different fields are :
- Seq_ID: (int) Id of the sequence.
- Codon_Pos: (int) Codon position.
- Nucleotide_Pos: (int) Nucleotide position.
- From: (String) Original sub-sequence value in the consensus.
- To: (String) New sub-sequence value in the deviations.
- MAF: (float) Minor allele frequency.
- Amino_Acid_From: (String) Original amino acid in the consensus.
- Amino_Acid_To: (String) New amino acid in the deviation.
- Synonymous: (bool) Indicates if the mutation is synonymous or not.  


You can specify with the flag --codon if you want to add the codons (both From and To) to the deviations output.
  

## Step 3  
  

### Input

The Step 3 takes as mandatory input the output file of the [Step 2](#step-2-1). Thus, look at the output of the previous step to see the input format for this one.  

### Script Parameters:  
  
- deviations: (String) Name of  the input file containing the deviations as specified in the output of Step 2. 
- output: (String) Name of the output file. Be sure to add the .json extension at the end (default="mut_spec.json").
- c: Flag indicating if the context should be taken into account for the mutational spectrum generation (default=False).
- h: (flag) Opens the Help section. 
  
### Output  
  
The Step 3 outputs a single JSON file containg the mutational spectrums for every original queries. The output format depends depends if the context is specified or not. Recall that the context is simply the nucleotide on the left and the one on the right of the mutation. If the context is specified, then the file will be a dictionary whose keys are the query indices and values are all the potential SNP with context (so 192 entries per query). If context is not specified then it would be simply SNP (so 12 entries per query).

# Automated browser vs Biopython

## Browser
  
PROS:  
- Much faster (x15)
- More reliable with big input sequences (do not break because of CPU load)

CONS:
- If the NCBI BLAST web page changes, we need to change the code.
- Cannot specify genetic code (but it is a NCBI tblastn flaw, so not really a problem).
- Need to have chrome installed (selenium's installation is pretty straightforward), which can be weird to have on remote servers.
- Hardcoded waiting timing (e.g. to avoid clicking too fast), which can be sensitive to internet connection.
- Kind of a "workaround" when you have an API for this website. Thus less elegant.  
  
  
## Biopython  
  
PROS:
- Elegant.
- Easier in implementation in other machines.
- Can specify the genetic code.

CONS:
- Super slow (15min per queries of a 100 hits).
- Not robust with big input sequences (break because of CPU load). 
  
  
## Conclusion  
    
In my sense, the choice really depends of the user. I think that for something more "official", like the implementation of the script on remote servers, the API choice is better. But if it is for personnal use on a personnal machine, then the browser approach is more convenient as it gives results way faster, is more robust and is as easy of use.  
  
  
# Limitations

Unfortunately there exists a problem with BLAST queries over Biopython. The queuing time of BLAST queries has a high variability and it is not unlikely to wait 10-15 minutes **per** query. Thus if you have a large number of queries, you can expect to wait for a day before getting the output. This might not be a problem for you but if it is there are some workarounds:

- Make several smaller control files and launch the script multiple times in parallel (but be careful that BLAST does not IP ban you for spamming).
- If you have multiple queries for a single specie, use the single_specie parameter (!!! not implemented yet).
- Use the automated browser approach.

