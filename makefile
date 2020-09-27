DEST_DIR = ~/bin

CFLAGS =  -O3 -Wall -Wextra -std=c++11

ALL =   kmer_binned_sampling.o sampled_kmer_filtering.o kmer_bin_presence_generator.o shared_contig_extractor.o

all: $(ALL)

kmer_binned_sampling.o:
		g++ $(CFLAGS) -o kmer_binned_sampling.o kmer_binned_sampling.cpp

sampled_kmer_filtering.o:
		g++ $(CFLAGS) -o sampled_kmer_filtering.o sampled_kmer_filtering.cpp

kmer_bin_presence_generator.o:
		g++ $(CFLAGS) -o kmer_bin_presence_generator.o kmer_bin_presence_generator.cpp

shared_contig_extractor.o:
		g++ $(CFLAGS) -o shared_contig_extractor.o shared_contig_extractor.cpp

clean:
		rm -f $(ALL)

install:
		cp $(ALL) $(DEST_DIR)