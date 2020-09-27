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

using Clock = std::chrono::steady_clock;
using std::chrono::time_point;
using std::chrono::duration_cast;
using std::chrono::milliseconds;
using namespace std;

typedef pair<string,short> pss;
typedef unsigned long long int ulli;
typedef pair<short,list<ulli> > pslistulli;
typedef pair<ulli,ulli> pulliulli;
typedef pair<unsigned long,ulli> pululli;
typedef pair<ulli,bool> pullib;
typedef pair<unsigned int,bool> puib;
typedef pair<int,bool> pib;
typedef tuple<ulli,int,bool> tup_uib;
typedef tuple<ulli,ulli,bool> tup_uub;
typedef pair<short,ulli> psulli;
typedef tuple<ulli,psulli,bool> tup_upsub;
typedef tuple<string,ulli,ulli,ulli> tup_suuu;
typedef tuple<unsigned int,ulli,ulli,ulli> tup_iuuu;
typedef pair<string, unsigned int> psui;
typedef tuple<unsigned int,bool,unsigned int,bool,float,float,short> tup_ibibffs;
typedef tuple<psulli,bool,unsigned int> tup_psubi;
typedef tuple<psulli,psulli,bool,bool> tup_psupsubb;
typedef tuple<ulli,tup_psupsubb> tup_ulli_tpsupsubb;
typedef tuple<tup_psubi,tup_psubi> tup_tpsubitpsubi;
typedef pair<psulli,vector<bool> > ppsuvecb;
typedef pair<tup_psupsubb,vector<bool> > ptuppsuvecb;
typedef pair<psulli,list<ulli> > ppsulistulli;
typedef pair<tup_psupsubb,list<ulli> > ptuppsulistulli;
typedef pair<psulli,ulli > ppsuulli;
typedef tuple<ulli,ulli,bool,bool> tup_uubb;
typedef tuple<unsigned long,unsigned long,bool,bool> tup_ululbb;
typedef tuple<ulli,ulli,int,bool> tup_uuib;
typedef tuple<ulli,bool,int> tup_ullibi;
typedef tuple<tup_ullibi,tup_ullibi> tup_tubitubi;
typedef tuple<unsigned long,bool,int> tup_ulbi;
typedef tuple<tup_ulbi,tup_ulbi> tup_tulbitulbi;
typedef pair<tup_ululbb,vector<bool> > ptupulvecb;
typedef pair<short,pullib> ps_pullib;
typedef pair<ulli,pullib> pullipullib;


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


list<ulli> get_current_prefix_filtered_suffixes(string filtered_kmer_dir, string prefix_idx_str){
    list<ulli> filtered_suffix_list;
    ulli filtered_kmer_count, filtered_suffix_idx;
    ifstream inFile_filtered;

    inFile_filtered.open(filtered_kmer_dir + prefix_idx_str, ios::binary);
    inFile_filtered.read((char*) (&filtered_kmer_count), sizeof(filtered_kmer_count));
    while(filtered_kmer_count--){
        inFile_filtered.read((char*) (&filtered_suffix_idx), sizeof(filtered_suffix_idx));
        filtered_suffix_list.push_back(filtered_suffix_idx);
    }

    cout<<"get_current_prefix_filtered_suffixes: "<<prefix_idx_str<<" "<<filtered_suffix_list.size()<<"\n";
    return filtered_suffix_list;
}


/*inline bool file_exists (string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}*/


void locate_and_save_filtered_kmer_bin_presence(string binned_kmer_dir, string assembly_idx_str, string prefix_idx_str,
                                                list<ulli>filtered_suffix_list, string contig_len_dir, string assembly_map_outdir){
    cout<<"locate_and_save_filtered_kmer_bin_presence started: "<<assembly_idx_str<<" "<<prefix_idx_str<<"\n";
    list<ulli>::iterator it_filtered_suffix;
    ulli binned_count, repeat_binned_count, check, repeat_check, assembly_suffix_idx, current_bin_no, repeat_count;
    ifstream inFile_binned_mer, inFile_binned_contig_bin, inFile_binned_mer_repeat, inFile_binned_contig_bin_repeat;

    map< ulli, unsigned long > contig_bin_representation_map;
    map< ulli, unsigned long >::iterator it_contig_bin_map;

    // Repeat kmers:
    list<ulli> repeat_kmer_suffixes_list, inner_bins_list;
    list< list<ulli> > repeat_kmer_bins_list;
    list< list<ulli> >::iterator it_repeat_bins_outerlist;
    list<ulli>::iterator it_repeat_suffix;

    inFile_binned_mer_repeat.open(binned_kmer_dir + assembly_idx_str + "_" + prefix_idx_str + "_mer_repeat", ios::binary);
    inFile_binned_contig_bin_repeat.open(binned_kmer_dir + assembly_idx_str + "_" + prefix_idx_str + "_bin_repeat", ios::binary);

    inFile_binned_mer_repeat.read((char*) (&repeat_binned_count), sizeof(repeat_binned_count));
    inFile_binned_contig_bin_repeat.read((char*) (&repeat_check), sizeof(repeat_check));

    cout<<"\t\t"<<repeat_binned_count<<" "<<repeat_check<<"\n";

    assert(repeat_binned_count==repeat_check);

    while(repeat_binned_count--){
        inFile_binned_mer_repeat.read((char*) (&assembly_suffix_idx), sizeof(assembly_suffix_idx));
        inFile_binned_contig_bin_repeat.read((char*) (&repeat_count), sizeof(repeat_count));

        repeat_kmer_suffixes_list.push_back(assembly_suffix_idx);

        while(repeat_count--){
            inFile_binned_contig_bin_repeat.read((char*) (&current_bin_no), sizeof(current_bin_no));
            inner_bins_list.push_back(current_bin_no);
        }
        repeat_kmer_bins_list.push_back(inner_bins_list);
        inner_bins_list.clear();
    }

    inFile_binned_mer_repeat.close();
    inFile_binned_contig_bin_repeat.close();



    list<ulli> end_bins_list = get_contig_end_bin_list(contig_len_dir + assembly_idx_str);
    list<ulli>::iterator it_end_bin = end_bins_list.begin();

    inFile_binned_mer.open(binned_kmer_dir + assembly_idx_str + "_" + prefix_idx_str + "_mer", ios::binary);
    inFile_binned_contig_bin.open(binned_kmer_dir + assembly_idx_str + "_" + prefix_idx_str + "_bin", ios::binary);

    inFile_binned_mer.read((char*) (&binned_count), sizeof(binned_count));
    inFile_binned_contig_bin.read((char*) (&check), sizeof(check));

    assert(binned_count==check);

    it_filtered_suffix = filtered_suffix_list.begin();
    it_repeat_suffix = repeat_kmer_suffixes_list.begin();
    it_repeat_bins_outerlist = repeat_kmer_bins_list.begin();

    while( (binned_count>0 || it_repeat_suffix != repeat_kmer_suffixes_list.end()) &&
            it_filtered_suffix != filtered_suffix_list.end() && it_end_bin != end_bins_list.end()){
        //TO DO: check if the filtered kmer is a repeat kmer: if yes, append the corresponding bins in a file:: use the bins from this file to update the maps later
        inFile_binned_mer.read((char*) (&assembly_suffix_idx), sizeof(assembly_suffix_idx));
        inFile_binned_contig_bin.read((char*) (&current_bin_no), sizeof(current_bin_no));
        binned_count--;

        while(it_end_bin != end_bins_list.end() && current_bin_no > *it_end_bin)
            it_end_bin++;

        if(it_end_bin == end_bins_list.end()){
            cout<<"ALERT!!! ERROR IN PARSING AND MAPPING CONTIG BINS: "<<current_bin_no<<" "<<assembly_idx_str<<" "<<prefix_idx_str<<"\n";
            continue;
        }

        while( it_filtered_suffix != filtered_suffix_list.end() && it_repeat_suffix != repeat_kmer_suffixes_list.end() &&
                (( binned_count>0 && assembly_suffix_idx > *it_repeat_suffix) || ( binned_count==0)) ){
            // Check if the repeat kmer is a filtered kmer
            if( *it_filtered_suffix < *it_repeat_suffix){
                it_filtered_suffix++;
                continue;
            }
            else{
                if( *it_filtered_suffix == *it_repeat_suffix){
                    // Add corresponding kmer bins into a list to be used to update the maps
                    inner_bins_list.insert(inner_bins_list.end(), (*it_repeat_bins_outerlist).begin(), (*it_repeat_bins_outerlist).end());
                    cout<<"\trepeat kmer: ("<<prefix_idx_str<<","<< *it_repeat_suffix <<"): "<<inner_bins_list.size()<<"\n";
                }
                // To be incremented even if the repeat kmer is not filtered selected
                it_repeat_suffix++;
                it_repeat_bins_outerlist++;
            }
        }

        while(it_filtered_suffix != filtered_suffix_list.end() && assembly_suffix_idx > *it_filtered_suffix)
            it_filtered_suffix++;

        if(it_filtered_suffix != filtered_suffix_list.end()){
            if(*it_filtered_suffix == assembly_suffix_idx){
                // current sampled k-mer is present in multiple samples
                
                it_contig_bin_map = contig_bin_representation_map.find( *it_end_bin );

                if(it_contig_bin_map == contig_bin_representation_map.end()){
                    contig_bin_representation_map[ *it_end_bin ] = 1;
                }
                else{
                    contig_bin_representation_map[ *it_end_bin ] += 1;
                }
            }
        }
    }

    inFile_binned_mer.close();
    inFile_binned_contig_bin.close();

    // Remaining repeat kmers
    while( it_filtered_suffix != filtered_suffix_list.end() && it_repeat_suffix != repeat_kmer_suffixes_list.end()){
        // Check if the repeat kmer is a filtered kmer
        if( *it_filtered_suffix < *it_repeat_suffix){
            it_filtered_suffix++;
            continue;
        }
        else{
            if( *it_filtered_suffix == *it_repeat_suffix){
                // Add corresponding kmer bins into a list to be used to update the maps
                inner_bins_list.insert(inner_bins_list.end(), (*it_repeat_bins_outerlist).begin(), (*it_repeat_bins_outerlist).end());
                cout<<"\trepeat kmer: ("<<prefix_idx_str<<","<< *it_repeat_suffix <<"): "<<inner_bins_list.size()<<"\n";
            }
            // To be incremented even if the repeat kmer is not filtered selected
            it_repeat_suffix++;
            it_repeat_bins_outerlist++;
            it_filtered_suffix++;
        }
    }

    cout<<"Unique kmers bin map size "<<assembly_idx_str + "_" + prefix_idx_str<<" : "<<contig_bin_representation_map.size()<<"\n";

    for(it_end_bin = inner_bins_list.begin(); it_end_bin != inner_bins_list.end(); it_end_bin++){
        it_contig_bin_map = contig_bin_representation_map.find( *it_end_bin );

        if(it_contig_bin_map == contig_bin_representation_map.end()){
            contig_bin_representation_map[ *it_end_bin ] = 1;
        }
        else{
            contig_bin_representation_map[ *it_end_bin ] += 1;
        }
    }

    cout<<"Total kmers bin map size "<<assembly_idx_str + "_" + prefix_idx_str<<" : "<<contig_bin_representation_map.size()<<"\n";


    // Save assembly prefix contig bin representation maps

    ofstream outFile;
    ulli map_size = contig_bin_representation_map.size();
    unsigned long bin_count;

    outFile.open(assembly_map_outdir + assembly_idx_str + "_" + prefix_idx_str, ios::binary);
    outFile.write((char*) (&map_size), sizeof(map_size));

    for(const auto &bincount_keyval_pair : contig_bin_representation_map){
        current_bin_no = bincount_keyval_pair.first;
        bin_count = bincount_keyval_pair.second;
        outFile.write((char*) (&current_bin_no), sizeof(current_bin_no));
        outFile.write((char*) (&bin_count), sizeof(bin_count));
    }

    outFile.close();

    remove((binned_kmer_dir + assembly_idx_str + "_" + prefix_idx_str + "_mer").c_str());
    remove((binned_kmer_dir + assembly_idx_str + "_" + prefix_idx_str + "_bin").c_str());
    remove((binned_kmer_dir + assembly_idx_str + "_" + prefix_idx_str + "_mer_repeat").c_str());
    remove((binned_kmer_dir + assembly_idx_str + "_" + prefix_idx_str + "_bin_repeat").c_str());

    cout<<"locate_and_save_filtered_kmer_bin_presence ended: "<<assembly_idx_str<<" "<<prefix_idx_str<<"\n";
}


void kmer_bin_presence_wrapper( string binned_kmer_dir, string filtered_kmer_dir, unsigned long assembly_count,
                                string prefix_idx_str, string contig_len_dir, string assembly_map_outdir){

    list<ulli> filtered_suffix_list = get_current_prefix_filtered_suffixes(filtered_kmer_dir, prefix_idx_str);

    for(unsigned long assembly_idx=0; assembly_idx < assembly_count; assembly_idx++)
        locate_and_save_filtered_kmer_bin_presence( binned_kmer_dir, to_string(assembly_idx), prefix_idx_str, filtered_suffix_list,
                                                    contig_len_dir, assembly_map_outdir);
}


int main(int argc, char** argv){
    cout<<"kmer_bin_presence_generator.cpp "<<argv<<"\n";
    unsigned long assembly_count = stoul(argv[1]);
    string prefix_idx_str = argv[2];
    string binned_kmer_dir = argv[3];
    string filtered_kmer_dir = argv[4];
    string assembly_map_outdir = argv[5];
    string contig_len_dir = argv[6];

    kmer_bin_presence_wrapper(  binned_kmer_dir, filtered_kmer_dir, assembly_count, prefix_idx_str,
                                contig_len_dir, assembly_map_outdir);

    return 0;
}