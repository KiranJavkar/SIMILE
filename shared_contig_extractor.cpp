#include <iostream>
#include <cstdio>
#include <list>
#include <vector>
#include <map>
#include <string>
#include <cstring>
#include <cassert>
#include <fstream>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <algorithm>
#include <iterator>
#include <chrono>
#include <set>
#include <utility>
#include <math.h>

using namespace std;

typedef unsigned long long int ulli;


ulli strtoulli(string number){
    ulli result = 0;
    for(int pos = 0; pos< number.length(); pos++){
        result *=10;
        result += number[pos]-'0';
    }
    return result;
}


bool is_valid_split_string(string item){
    if(item.length()>1)
        return true;
    if(item.length()==0)
        return false;
    char c = item[0];
    switch(c){
        case ' ': case '\t': case '\n': case '\b': case '\r': case ',': case '.':
        return false;
        default:
        return true;
    }
}


template<typename Out>
void split(const string &s, char delim, Out result) {
    stringstream ss(s);
    string item;
    while (std::getline(ss, item, delim)) {
        if(is_valid_split_string(item))
            *(result++) = item;
    }
}


list<string> split(const string &s, char delim) {
    list<string> elems;
    split(s, delim, back_inserter(elems));
    return elems;
}


list<string> get_lines(const char* FILENAME){
    ifstream input_file(FILENAME);
    list<string> lines_list;
    string text;
    while (getline( input_file, text )){
        lines_list.push_back(text);
    }
    return lines_list;
}


list< ulli > get_contig_end_bin_list(string contig_len_filename){
    list<string> lines = get_lines(contig_len_filename.c_str());
    list<string> lsplit;
    list< ulli > end_bins_list;
    list<string> :: iterator line_it, it;

    for(line_it = lines.begin(); line_it != lines.end(); line_it++){
        lsplit = split(*line_it, '\t');
        it = lsplit.begin();
        it++; // end bin no is in the second column
        end_bins_list.push_back(strtoulli( *it ));
    }

    return end_bins_list;
}

list< bool > aggregate_contig_bin_presence( string assembly_map_dir, string assembly_idx_str, short prefix_count,
                                            float presence_threshold, unsigned long sampling_count, list< ulli > end_bins_list){
    cout<<"aggregate_contig_bin_presence started: "<<assembly_idx_str<<" "<<end_bins_list.size()<<" "<<presence_threshold<<"\n";
    // Load contig bin presence (for the current sample) for a prefix
    // Using the contig ends bins, update the filtered sampled kmers count representing that contig
    // Repeat for all prefixes

    // Use a threshold \gamma (%) that need to present for a contig to be deemed present
    // if a contig is comprised of b bins, each containing sampling_count kmers,
    //      the contig would be deemed to be present if it is represented by \gamma * sampled_count * b kmers

    ulli map_size, contig_no, current_bin_no;
    unsigned long bin_count;

    list<ulli> contig_kmer_count_list;
    for(contig_no=0; contig_no < end_bins_list.size(); contig_no++)
        contig_kmer_count_list.push_back(0);

    list<ulli>::iterator it_end_bin, it_contig_kmer_count;

    ifstream inFile;

    for(short prefix_idx = 0; prefix_idx < prefix_count; prefix_idx++){
        it_end_bin = end_bins_list.begin();
        it_contig_kmer_count = contig_kmer_count_list.begin();

        inFile.open(assembly_map_dir + assembly_idx_str + "_" + to_string(prefix_idx), ios::binary);
        inFile.read((char*) (&map_size), sizeof(map_size));

        while(map_size--){
            inFile.read((char*) (&current_bin_no), sizeof(current_bin_no));
            inFile.read((char*) (&bin_count), sizeof(bin_count));

            while(it_end_bin != end_bins_list.end() && *it_end_bin < current_bin_no){
                it_end_bin++;
                it_contig_kmer_count++;
            }
            if(it_end_bin != end_bins_list.end()){
                (*it_contig_kmer_count) += bin_count;
            }
            else
                cout<<"ALERT!!! ERROR IN PARSING CONTIG BIN MAPS!!! "<<assembly_idx_str<<" "<<prefix_idx<<"\n";
        }

        inFile.close();
        // remove((assembly_map_dir + assembly_idx_str + "_" + to_string(prefix_idx)).c_str());
    }


    // Determine contig shared presence

    it_end_bin = end_bins_list.begin();
    it_contig_kmer_count = contig_kmer_count_list.begin();

    ulli current_contig_bins, current_start_bin_no = 0, current_min_presence, selected_contigs = 0;
    bool presence;
    string out_string = "";
    contig_no = 0;

    list<bool> contig_presence_list;

    while(it_contig_kmer_count != contig_kmer_count_list.end()){
        current_contig_bins = *it_end_bin - current_start_bin_no + 1;
        current_min_presence = presence_threshold * sampling_count * current_contig_bins;
        if(current_min_presence < 3) // Ad hoc threshold used for boundary conditions
            current_min_presence = 3;

        presence = (*it_contig_kmer_count) >= current_min_presence;
        if(presence)
            selected_contigs++;
        contig_presence_list.push_back(presence);

        out_string += to_string(contig_no) + "\t" + to_string(current_contig_bins) + "\t" + to_string(*it_contig_kmer_count) + "\t";
        out_string += to_string(current_min_presence) + (presence?"\t1\n":"\t0\n");

        current_start_bin_no = *it_end_bin + 1;

        it_contig_kmer_count++;
        it_end_bin++;
        contig_no++;
    }

    ofstream outFile;
    outFile.open(assembly_map_dir + assembly_idx_str , ios::out);
    outFile << out_string;
    outFile.close();


    cout<<"aggregate_contig_bin_presence ended: "<<assembly_idx_str<<" "<<end_bins_list.size()<<" "<<presence_threshold<<" "<<selected_contigs<<"\n";

    return contig_presence_list;
}


void extract_shared_contig_fasta_sequences(string fasta_filename, list<bool> contig_presence_list, string shared_contigs_filename){
    cout<<"extract_shared_contig_fasta_sequences shared: "<<fasta_filename<<" "<<shared_contigs_filename<<"\n";

    list<bool>::iterator it_presence = contig_presence_list.begin();
    bool current_contig_presence, stop = false;

    ifstream input_file(fasta_filename);
    string line, out_string="";
    ofstream outFile;

    while (!stop && getline( input_file, line )){
        if(line[0] == '>'){
            if(out_string.length() > 0){
                outFile.open(shared_contigs_filename, ios::app);
                outFile << out_string;
                outFile.close();
                out_string = "";
            }

            if(it_presence == contig_presence_list.end()){
                stop=true;
                break;
            }

            current_contig_presence = *it_presence;

            if(*it_presence){
                // Shared contig header
                out_string = line + "\n";
            }
            it_presence++;
        }
        else{
            // Contig FASTA
            if(current_contig_presence)
                out_string += line + "\n";
        }
    }

    if(out_string.length() > 0){
        outFile.open(shared_contigs_filename, ios::app);
        outFile << out_string;
        outFile.close();
    }

    cout<<"extract_shared_contig_fasta_sequences ended: "<<fasta_filename<<" "<<shared_contigs_filename<<"\n";
}


int main(int argc, char** argv){
    cout<<"shared_contig_extractor.cpp "<<argv<<"\n";
    string assembly_idx_str = argv[1];
    short prefix_count = static_cast<short>(stoi(argv[2]));
    float presence_threshold = stof(argv[3]);
    presence_threshold /= 100.0;
    unsigned long sampling_count = stoul(argv[4]);
    string assembly_map_dir = argv[5];
    string contig_len_dir = argv[6];
    string fasta_filename = argv[7];
    string shared_contigs_filename = argv[8];

    list< ulli > end_bins_list = get_contig_end_bin_list(contig_len_dir + assembly_idx_str);

    list< bool > contig_presence_list  = aggregate_contig_bin_presence( assembly_map_dir, assembly_idx_str, prefix_count,
                                                                        presence_threshold, sampling_count, end_bins_list);

    extract_shared_contig_fasta_sequences(fasta_filename, contig_presence_list, shared_contigs_filename);

    return 0;
}