import sys
import pandas as pd
import numpy as np
import math

######### Usage ######### 
## Convert nanopolish output to BED
# awk -v OFS="\t" '{ if (NR!=1) print $1,$3,$4,$5,$6,$2 }' test.methylation_calls.tsv > test.methylation_calls.bed

## Pad CGI bed to include shores
#bedtools slop -i grch38_cgi.bed -g grch38.chrom.sizes -b 2000 > grch38_cgi_shores.bed

## Restrict methylation calls to those falling outside of islands and shores
#bedtools intersect -a test.methylation_calls.bed -b grch38_cgi_shores.bed -v > test.methylation_calls.nocgi.bed

## Summarize methylation by read using this script
#python methylation_by_read.py test.methylation_calls.nocgi.bed test.methylation_by_read.tsv 2.5
########################### 

#infile = "test.methylation_calls.nocgi.bed"
#outfile = "test.methylation_by_read.tsv"
#log_lik_threshold = 2.5

infile = sys.argv[1]
outfile = sys.argv[2]
log_lik_threshold = float(sys.argv[3])

def iter_chunk_by_read_name(file):
    csv_reader = pd.read_csv(file, iterator=True, chunksize=1000000, sep='\t', header=None, names=["chromosome", "start", "end", "read_name", "log_lik_ratio", "strand"])
    last_read_in_chunk=pd.DataFrame()
    chunk_id = 0
    for chunk in csv_reader:
        chunk = pd.concat([last_read_in_chunk, chunk])
        methylation_calls = pd.DataFrame(chunk)
        # Set aside rows corresponding to the last read in the chunk so that they 
        # can be grouped with the next chunk
        last_read_name_in_chunk = methylation_calls.tail(1)["read_name"].iloc[0]
        last_read_in_chunk = methylation_calls.loc[methylation_calls['read_name'] == last_read_name_in_chunk]
        reads = methylation_calls.groupby(["read_name"])
        for read_name, cpgs in reads:
            if read_name != last_read_name_in_chunk:
                #print(chunk_id, read_name, len(cpgs.index))
                yield(cpgs)
        chunk_id = chunk_id + 1
    #print(last_read_name_in_chunk)
    yield(last_read_in_chunk)


of = open(outfile, 'w')

read_iter = iter_chunk_by_read_name(infile)

for read in read_iter:
    chromosome = read["chromosome"].iloc[0]
    read_name = read["read_name"].iloc[0]
    read_start_pos = min(read["start"])
    read_end_pos = max(read["end"])
    m = 0
    u = 0
    for index, locus in read.iterrows():
        #print(locus['read_name'], locus['log_lik_ratio'])
        if locus.log_lik_ratio > log_lik_threshold:
            m = m + 1
        if locus.log_lik_ratio < -log_lik_threshold:
            u = u + 1
    #print ('\t'.join([str(chromosome), str(read_start_pos), str(read_end_pos), read_name, str(m), str(u)]) + "\n")
    of.write('\t'.join([str(chromosome), str(read_start_pos), str(read_end_pos), read_name, str(m), str(u)]) + "\n")


