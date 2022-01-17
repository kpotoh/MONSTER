import os

os.system("tblastn -db nt -query query_test.txt -remote -entrez_query \"Mus Musculus Domesticus[organism]\" -out OUTPUT_BLAST.xml -outfmt 5")

