# SIMILE: Discover shared genomic regions across a collection of metagenomic samples

Given a collection of metagenome-assembled-genomes (MAGs), extract the contigs that represent the genomic regions that are shared across a certain proportion ($\phi$) of the given MAGs i.e. metagenomic samples.

### Dependencies
1. Python 3.5 or later
2. C++11 or later

### Installation:
```bash
git clone https://github.com/KiranJavkar/SIMILE.git
cd PRAWNS
make
```

### Running the command:
```
python run_simile.py input.csv
```
Where ```input.csv``` comprises of 2 columns: sample_name, fasta_file_path
```
-bash-4.2$ python run_simile.py -h
usage: run_simile.py [-h] -i INPUT [-n [NCORES]] [-K [KMER_LEN]]
                     [-p [MIN_PERC]] [-f [SAMPLING_FREQUENCY]] [-b [BIN_SIZE]]
                     [-c [SAMPLED_KMER_CONTIG_PRESENCE]] [-m [MEM]]

SIMILE: Discover shared genomic regions across a collection of metagenomic
samples

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input csv file
  -n [NCORES], --ncores [NCORES]
                        Number of cores to be used (default: 8)
  -K [KMER_LEN], --kmer_len [KMER_LEN]
                        Length of kmers (default: 25)
  -p [MIN_PERC], --min_perc [MIN_PERC]
                        Minimum % of samples that should contain the shared
                        region (default: 5.0)
  -f [SAMPLING_FREQUENCY], --sampling_frequency [SAMPLING_FREQUENCY]
                        Proportion of k-mers to be sampled and selected for
                        comparisons (default 1.0)
  -b [BIN_SIZE], --bin_size [BIN_SIZE]
                        Size of bins for sampling the k-mers (default: 1000)
  -c [SAMPLED_KMER_CONTIG_PRESENCE], --sampled_kmer_contig_presence [SAMPLED_KMER_CONTIG_PRESENCE]
                        Minimum % of the sampled k-mers that should have
                        shared presence for the contig to be considered as a
                        shared contig (default: 40.0)
  -m [MEM], --mem [MEM]
                        Upper limit for RAM memory usage. Can be in
                        mb/MB/gb/GB/tb/TB (case insensitive), default unit is
                        MB. (default: 36000MB)
```


