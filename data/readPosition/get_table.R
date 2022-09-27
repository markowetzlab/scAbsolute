# A script to take reference sizes and reference read start positions as input and convert to genomic co-ordinates

library(dplyr)
library(magrittr)
library(readr)
library(stringr)

args = commandArgs(trailingOnly = TRUE)

# debugging
#args = character(4) 
#args[1] = "~/Data/project1/20-align/prediction/UID-DLP-SA1044_SLX-A96139A_000603_R18-C44.bam.reference_sizes.txt"
#args[2] = "~/Data/project1/20-align/prediction/assembly.tsv"
#args[3] = "~/Data/project1/20-align/prediction/UID-DLP-SA1044_SLX-A96139A_000603_R18-C44.bam.deduplicated.read_start_sites.txt"
#args[4] = "out.position.tsv"

# get reference sizes
reference_sizes = readr::read_tsv(args[1], comment='',col_names=c('reference_name','reference_size'), col_types="cc")

# reformat
reference_sizes$reference_name = reference_sizes$reference_name %>% stringr::str_remove('SN:')
reference_sizes$reference_size = as.numeric(reference_sizes$reference_size %>% stringr::str_remove('LN:'))

# remove unwanted reference sequences
relevant_reference_sequences = readr::read_tsv(args[2], col_names=c('reference_name', "X2"), col_types="cd") # a single column file with column of relevant reference sequence identifiers
reference_sizes = reference_sizes %>% dplyr::filter(reference_name %in% relevant_reference_sequences$reference_name)

# calculate the cumulative sum of the reference sizes
reference_sizes$cumulative_sum = cumsum(reference_sizes$reference_size)

# transform from length to a cumulative start position
reference_sizes$reference_genomic_start = NA
reference_sizes$reference_genomic_start[1] = 0

for (i in 2:length(reference_sizes$cumulative_sum)) {
	reference_sizes$reference_genomic_start[i] = reference_sizes$cumulative_sum[i-1]
}
reference_sizes = reference_sizes %>% dplyr::select(-cumulative_sum)

# read in the read start positions
read_start_sites = readr::read_tsv(args[3], col_names=c('query_name','reference_name','query_start_position_chromosome', 'query_tlen'), col_types="ccdd")

# join the tables together
joined_table = dplyr::inner_join(read_start_sites, reference_sizes, by='reference_name')

# compute the query start position with respect to the genome
joined_table$query_start_position_genome = as.numeric(joined_table$reference_genomic_start) + as.numeric(joined_table$query_start_position_chromosome)

raw_read_table = joined_table %>% dplyr::select(reference_name, query_start_position_chromosome, query_tlen) %>% dplyr::rename("chromosome"="reference_name", "start"="query_start_position_chromosome", 'length'='query_tlen')

#write to disk
readr::write_tsv(raw_read_table, args[4], append=FALSE)
