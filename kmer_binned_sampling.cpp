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


void get_prefix_dict(map<string, short> &prefix_dict){
    char nucleotides[] = "ACGT";
    short prefix_id_value = 0;
    string current_prefix = "AAAAA";
    for(short idx_0=0; idx_0<4; idx_0++){
        current_prefix[0] = nucleotides[idx_0];
        for(short idx_1=0; idx_1<4; idx_1++){
            current_prefix[1] = nucleotides[idx_1];
            for(short idx_2=0; idx_2<4; idx_2++){
                current_prefix[2] = nucleotides[idx_2];
                for(short idx_3=0; idx_3<4; idx_3++){
                    current_prefix[3] = nucleotides[idx_3];
                    for(short idx_4=0; idx_4<4; idx_4++){
                        current_prefix[4] = nucleotides[idx_4];
                        prefix_dict.insert(pss(current_prefix, prefix_id_value++));
                    }
                }
            }
        }
    }
}


ulli get_suffix_index(string suffix){
    ulli index = 0;
    for(string::iterator it=suffix.begin(); it!=suffix.end(); it++){
        // index*=4;
        index = index<<2;
        switch(*it){
            case 'A': case 'a': index += 0;
            break;
            case 'C': case 'c': index += 1;
            break;
            case 'G': case 'g': index += 2;
            break;
            case 'T': case 't': index += 3;
            break;
        }
    }
    return index;
}


string get_nucleotide_string(ulli value, short length){
    string result = string(length, 'A');
    int pos = length-1;
    while(value>0){
        char nuc;
        switch(value%4){
            case 0: nuc='A'; break;
            case 1: nuc='C'; break;
            case 2: nuc='G'; break;
            case 3: nuc='T'; break;
        }
        result[pos--] = nuc;
        value /= 4;
    }
    return result;
}


string get_reverse_kmer(string kmer){
    string rev_kmer = string(kmer.length(), 'A');
    unsigned short pos = kmer.length()-1;
    for(string::iterator it=kmer.begin(); it!=kmer.end(); it++){
        switch(*it){
            case 'A': rev_kmer[pos--] = 'T';
            break;
            case 'C': rev_kmer[pos--] = 'G';
            break;
            case 'G': rev_kmer[pos--] = 'C';
            break;
            case 'T': rev_kmer[pos--] = 'A';
            break;
        }
    }
    return rev_kmer;
}


list<string> get_fasta_filepaths(const char* FILENAME){
    ifstream input_file(FILENAME);
    string filepath;
    list<string> fasta_filepaths;
    while (getline( input_file, filepath )){
        fasta_filepaths.push_back(filepath);
    }
    return fasta_filepaths;
}


void init_kmer_map_group_lists(list< map< ulli, ulli > >& kmer_map_list, list< map< ulli, list<ulli> > >& repeat_kmer_map_list){
    // cout<<"init_kmer_map_group_lists\n";
    // No. of lists = total number of prefixes ==> 4^5 = 1024
    map< ulli, ulli > new_kmer_map;
    map< ulli, list<ulli> > new_repeat_kmer_map;
    for(short idx=0; idx<1024; idx++){
        kmer_map_list.push_back(new_kmer_map);
        repeat_kmer_map_list.push_back(new_repeat_kmer_map);
    }
    // cout<<"init_map_group_lists completed\n";
}


ps_pullib get_lex_min_kmer(string kmer, map<string, short> &prefix_dict, short suffix_len){
    string rev_kmer = get_reverse_kmer(kmer);

    short forward_prefix_idx = prefix_dict[kmer.substr(0, 5)];
    short reverse_prefix_idx = prefix_dict[rev_kmer.substr(0, 5)];

    if(forward_prefix_idx<reverse_prefix_idx)
        return ps_pullib( forward_prefix_idx, pullib( get_suffix_index(kmer.substr(5, suffix_len)), true ) );
    else if(forward_prefix_idx>reverse_prefix_idx)
        return ps_pullib( reverse_prefix_idx, pullib( get_suffix_index(rev_kmer.substr(5, suffix_len)), false ) );
    else{
        short forward_suffix_idx = get_suffix_index(kmer.substr(5, suffix_len));
        short reverse_suffix_idx = get_suffix_index(rev_kmer.substr(5, suffix_len));
        if(forward_suffix_idx<=reverse_suffix_idx)
            return ps_pullib( forward_prefix_idx, pullib(forward_suffix_idx, true) );
        else
            return ps_pullib( reverse_prefix_idx, pullib(reverse_suffix_idx, false) );
    }
}


// set<psulli> add_sampled_kmers(priority_queue< ppsuulli >& sampled_kmers_pq, list< map< ulli, ulli > >& kmer_map_list,
//                                 list< map< ulli, list<ulli> > >& repeat_kmer_map_list){
set<psulli> add_sampled_kmers(priority_queue< psulli >& sampled_kmers_pq, list< map< ulli, ulli > >& kmer_map_list,
                                list< map< ulli, list<ulli> > >& repeat_kmer_map_list, ulli current_bin_no){
    if(current_bin_no%100==0)
        cout<<"add_sampled_kmers started: "<<sampled_kmers_pq.size()<<" "<<current_bin_no<<"\n";

    set <psulli> current_remove_kmer_set;
    set <psulli>::iterator it_set;

    list< map< ulli, ulli > >::iterator map_list_it;
    map< ulli, ulli >::iterator it_map;

    list< map< ulli, list<ulli> > >::iterator repeat_map_list_it;
    // map< ulli, list<ulli> >::iterator it_repeat_map;

    // ppsuulli kmer_idx_bin_pair;
    psulli kmer_idx;


    while(!sampled_kmers_pq.empty()){

        // kmer_idx_bin_pair = sampled_kmers_pq.top();
        kmer_idx = sampled_kmers_pq.top();
        sampled_kmers_pq.pop();

        // kmer_idx = kmer_idx_bin_pair.first;

        map_list_it = kmer_map_list.begin();
        advance(map_list_it, kmer_idx.first );

        it_map = (*map_list_it).find( kmer_idx.second );

        if(it_map == (*map_list_it).end()){
            // New kmer
            // (*map_list_it)[ kmer_idx.second ] =  (kmer_idx_bin_pair.second);
            (*map_list_it)[ kmer_idx.second ] =  current_bin_no;
        }
        else{
            // Repeated kmer
            cout<<"\t\tRepeat kmer detected: ("<<kmer_idx.first<<","<<kmer_idx.second<<") in "<<current_bin_no<<"\t";
            it_set = current_remove_kmer_set.find(kmer_idx);

            repeat_map_list_it = repeat_kmer_map_list.begin();
            advance(repeat_map_list_it, kmer_idx.first );

            if(it_set == current_remove_kmer_set.end()){
                // first occurence of kmer repeat
                current_remove_kmer_set.insert(kmer_idx);

                cout<<"\tPrevious occurrence in :"<< it_map->second <<"\n";

                // list<ulli> repeat_kmer_bins;
                // repeat_kmer_bins.push_back( it_map->second ); // Previously identified bin
                // repeat_kmer_bins.push_back(current_bin_no);

                // (*repeat_map_list_it)[ kmer_idx.second ] = repeat_kmer_bins;
                (*repeat_map_list_it)[ kmer_idx.second ].push_back( it_map->second );
                (*repeat_map_list_it)[ kmer_idx.second ].push_back(current_bin_no);
            }
            else{
                cout<<"\tAppending "<<current_bin_no<<" to existing repeat kmer bins of ("<<kmer_idx.first<<","<<kmer_idx.second<<")\n";
                (*repeat_map_list_it)[ kmer_idx.second ].push_back(current_bin_no);
            }
        }
    }

    if(current_bin_no%100==0)
        cout<<"add_sampled_kmers ended: "<<sampled_kmers_pq.size()<<" "<<current_bin_no<<" "<<current_remove_kmer_set.size()<<"\n";
    return current_remove_kmer_set;
}


void generate_binned_sampled_kmers(string fasta_filename, list< map< ulli, ulli > >& kmer_map_list,
                                    list< map< ulli, list<ulli> > >& repeat_kmer_map_list, map<string, short> &prefix_dict, short kmer_len,
                                    unsigned long sampling_count, unsigned long bin_size, string contig_len_filename){
    cout<<"generate_binned_sampled_kmers started: "<<fasta_filename<<"\n";
    short suffix_len = kmer_len-5;
    ulli position = 1; // Following the convention from jellyfish, mummer etc.

    ifstream input_file(fasta_filename);
    if(!input_file){
        cout<<"ERROR!!! UNABLE TO ACCESS FASTA FILE!!! "<<fasta_filename<<"\n";
        return;
    }
    string line;
    bool list_reset = false; // Used for moving the kmer length window through vector
    ps_pullib kmer_idx_strand_pair;

    short prefix_idx;
    ulli suffix_idx, current_bin_no = 0;
    // bool strand;
    // priority_queue< ppsuulli > sampled_kmers_pq;
    priority_queue< psulli > sampled_kmers_pq;
    unsigned long bin_position = 0;
    psulli kmer_idx, current_top_kmer_idx;
    ppsuulli kmer_idx_bin_pair;

    ulli sampled_kmer_count = 0;


    list <char> kmer_list;
    set <psulli> remove_kmer_set, current_remove_kmer_set;
    // kmer_list.reserve(kmer_len);

    list< map< ulli, ulli > >::iterator map_list_it;
    map< ulli, ulli >::iterator it;

    ofstream outFile;
    outFile.open(contig_len_filename, ios::out);

    while (getline( input_file, line )){
        if(line[0]=='>'){

            if(position>1){
                // add the sampled kmers from previous contig bin in the binned maps

                sampled_kmer_count += sampled_kmers_pq.size();
                current_remove_kmer_set =  add_sampled_kmers(sampled_kmers_pq, kmer_map_list, repeat_kmer_map_list, current_bin_no);
                remove_kmer_set.insert(current_remove_kmer_set.begin(), current_remove_kmer_set.end());
                bin_position = 0;

                // Add the trailing nucleotide count (typically kmer_len-1) from previous contig
                position += kmer_list.size();
                // 1 less than the current position value marks the end of the previous contig
                outFile << current_bin_no++ << "\t" << (position-1) << "\n";
            }

            outFile<<line<<"\t";
            // cout<<line<<"\n";
            kmer_list.clear();
            // kmer_list.reserve(kmer_len);
            list_reset = true;
            continue;
        }

        // cout<<line<<"\n";
        // if(!kmer_list.empty()){
        //     string overflow(kmer_list.begin(), kmer_list.end());
        //     cout<<overflow<<"\n";
        // }

        if(line.length()<kmer_len){
            position += line.length();
            continue;
        }

        for(char const &c: line){
            kmer_list.push_back(c);

            if(kmer_list.size()==kmer_len){

                string kmer(kmer_list.begin(), kmer_list.end());
                // cout<<kmer<<"  ";
                kmer_idx_strand_pair = get_lex_min_kmer(kmer, prefix_dict, suffix_len);
                prefix_idx = kmer_idx_strand_pair.first;
                suffix_idx = (kmer_idx_strand_pair.second).first;
                // strand = (kmer_idx_strand_pair.second).second;

                kmer_idx = psulli(prefix_idx, suffix_idx);
                kmer_idx_bin_pair = ppsuulli(kmer_idx, current_bin_no);


                // Check if current kmer is to be chosen for sampling within the current bin

                if(sampled_kmers_pq.size() < sampling_count){
                    // sampled_kmers_pq.push(kmer_idx_bin_pair);
                    sampled_kmers_pq.push(kmer_idx);
                    // cout<<"\t("<<kmer_idx.first<<","<<kmer_idx.second<<") ";
                }
                else{
                    // current_top_kmer_idx = psulli( ((sampled_kmers_pq.top()).first).first, ((sampled_kmers_pq.top()).first).second);
                    current_top_kmer_idx = sampled_kmers_pq.top();

                    if(current_top_kmer_idx >= kmer_idx){
                        while(current_top_kmer_idx > kmer_idx && sampled_kmers_pq.size() >= sampling_count){
                            // cout<<"\tXXXX("<<current_top_kmer_idx.first<<","<<current_top_kmer_idx.second<<") ";
                            sampled_kmers_pq.pop();
                            // current_top_kmer_idx = psulli( ((sampled_kmers_pq.top()).first).first, ((sampled_kmers_pq.top()).first).second);
                            current_top_kmer_idx = sampled_kmers_pq.top();
                        }
                        // sampled_kmers_pq.push(kmer_idx_bin_pair);
                        sampled_kmers_pq.push(kmer_idx);
                        // cout<<"\t<--("<<kmer_idx.first<<","<<kmer_idx.second<<") ";
                    }
                }

                list_reset = false;
                // cout<<prefix_idx<<" "<<suffix_idx<<" "<<"  "<<position<<" "<<strand<<"\n";
            }
            if(!list_reset)
                if(!kmer_list.empty()){
                    kmer_list.pop_front();
                    position++;
                    bin_position++;

                    // If contig bin has been traversed, add the sampled kmers in the binned maps

                    if(bin_position == bin_size){

                        sampled_kmer_count += sampled_kmers_pq.size();
                        current_remove_kmer_set =  add_sampled_kmers(sampled_kmers_pq, kmer_map_list, repeat_kmer_map_list, current_bin_no);
                        remove_kmer_set.insert(current_remove_kmer_set.begin(), current_remove_kmer_set.end());
                        bin_position = 0;
                        current_bin_no++;
                    }
                }
        }
    }

    if(bin_position > 0){

        sampled_kmer_count += sampled_kmers_pq.size();
        current_remove_kmer_set =  add_sampled_kmers(sampled_kmers_pq, kmer_map_list, repeat_kmer_map_list, current_bin_no);
        remove_kmer_set.insert(current_remove_kmer_set.begin(), current_remove_kmer_set.end());
        bin_position = 0;
        current_bin_no++;
    }

    outFile << current_bin_no << "\t" << (position-1) << "\n"; // End position of last contig

    cout<<"Sampled kmers: "<<sampled_kmer_count<<"\tRepeated kmers: "<<remove_kmer_set.size()<<"\n";

    map_list_it = kmer_map_list.begin();
    short remove_prefix_idx = 0;
    for (auto const &kmer_idx : remove_kmer_set){
        prefix_idx = kmer_idx.first;
        suffix_idx = kmer_idx.second;
        advance(map_list_it, prefix_idx-remove_prefix_idx);
        remove_prefix_idx = prefix_idx;
        it = (*map_list_it).find(suffix_idx);
        // cout<< "Removed: " << get_nucleotide_string(prefix_idx, 5) << get_nucleotide_string(suffix_idx, suffix_len)<<"(";
        // cout<<(it->second).first<<","<<(it->second).second<<")"<<"\n";
        it = (*map_list_it).erase(it);
    }

    // cout<<"generate_binned_sampled_kmers ended: "<<position<<"\n";
}


void save_binned_positional_kmers(list< map< ulli,ulli > > kmer_map_list, list< map< ulli, list<ulli> > > repeat_kmer_map_list,
                                string out_dir, string assembly_idx_str){
    // cout<<"save_binned_positional_kmers started: "<<assembly_idx_str<<"\n";
    // Save the kmer suffixes and corresponding kmer positions in separate files
    // Allows direct parsing of the suffixes for kmer filtering
    short prefix_idx=0;
    ulli kmer_count = 0, repeat_kmer_count = 0, binned_count, suffix_idx, current_bin_no, repeat_count;
    list<ulli> repeat_kmer_bins;
    // pullib pos_strand_pair;
    ofstream outFile_suffix, outFile_pos;
    string filename_prefix = out_dir + assembly_idx_str + "_";

    for(list< map<ulli,ulli> >::iterator map_list_it = kmer_map_list.begin(); map_list_it != kmer_map_list.end(); map_list_it++){
        binned_count = (*map_list_it).size();
        kmer_count += binned_count;

        outFile_suffix.open(filename_prefix + to_string(prefix_idx) + "_mer", ios::binary);
        outFile_pos.open(filename_prefix + to_string(prefix_idx) + "_bin", ios::binary);

        outFile_suffix.write((char*) (&binned_count), sizeof(binned_count));
        outFile_pos.write((char*) (&binned_count), sizeof(binned_count));

        // Maps are sorted by default due to heap structure

        for(const auto &kmer_keyval_pair : (*map_list_it)){
            suffix_idx = kmer_keyval_pair.first;
            current_bin_no = kmer_keyval_pair.second;
            outFile_suffix.write((char*) (&suffix_idx), sizeof(suffix_idx));
            outFile_pos.write((char*) (&current_bin_no), sizeof(current_bin_no));
        }

        // for(map< ulli, pullib >::iterator it = (*map_list_it).begin(); it != (*map_list_it).end(); it++){
        //     ulli& suffix_idx = it->first;
        //     pullib& pos_strand_pair = it->second;
        //     outFile_suffix.write((char*) (&suffix_idx), sizeof(suffix_idx));
        //     outFile_pos.write((char*) (&pos_strand_pair), sizeof(pos_strand_pair));
        // }

        prefix_idx++;
        outFile_suffix.close();
        outFile_pos.close();
    }

    prefix_idx = 0;
    for(list< map< ulli, list<ulli> > >::iterator repeat_map_list_it = repeat_kmer_map_list.begin();
            repeat_map_list_it != repeat_kmer_map_list.end(); repeat_map_list_it++){

        binned_count = (*repeat_map_list_it).size();
        repeat_kmer_count += binned_count;

        outFile_suffix.open(filename_prefix + to_string(prefix_idx) + "_mer_repeat", ios::binary);
        outFile_pos.open(filename_prefix + to_string(prefix_idx) + "_bin_repeat", ios::binary);

        outFile_suffix.write((char*) (&binned_count), sizeof(binned_count));
        outFile_pos.write((char*) (&binned_count), sizeof(binned_count));

        // Maps are sorted by default due to heap structure

        for(const auto &kmer_keyval_pair : (*repeat_map_list_it)){
            suffix_idx = kmer_keyval_pair.first;
            repeat_kmer_bins = kmer_keyval_pair.second;
            repeat_count = repeat_kmer_bins.size();
            outFile_suffix.write((char*) (&suffix_idx), sizeof(suffix_idx));
            outFile_pos.write((char*) (&repeat_count), sizeof(repeat_count));
            for(list<ulli>::iterator it_bins = repeat_kmer_bins.begin(); it_bins != repeat_kmer_bins.end(); it_bins++){
                current_bin_no = *it_bins;
                outFile_pos.write((char*) (&current_bin_no), sizeof(current_bin_no));
            }
        }

        // for(map< ulli, pullib >::iterator it = (*repeat_map_list_it).begin(); it != (*repeat_map_list_it).end(); it++){
        //     ulli& suffix_idx = it->first;
        //     pullib& pos_strand_pair = it->second;
        //     outFile_suffix.write((char*) (&suffix_idx), sizeof(suffix_idx));
        //     outFile_pos.write((char*) (&pos_strand_pair), sizeof(pos_strand_pair));
        // }

        prefix_idx++;
        outFile_suffix.close();
        outFile_pos.close();
    }
    cout<<"save_binned_positional_kmers ended: "<<assembly_idx_str<<": "<<kmer_count<<" "<<repeat_kmer_count<<"\n";
}


int main(int argc, char** argv){
    cout<<"kmer_binned_sampling.cpp "<<argc<<"\n";
    // string input_file = argv[1];
    // short kmer_len = static_cast<short>(stoi(argv[2]));
    // unsigned long sampling_count = stoul(argv[3]);
    // unsigned long bin_size = stoul(argv[4]);
    // unsigned long start_assembly_idx = stoul(argv[5]);
    // string out_dir = argv[6];
    // string contig_len_dir = argv[7];

    string assembly_fasta = argv[1];
    short kmer_len = static_cast<short>(stoi(argv[2]));
    unsigned long sampling_count = stoul(argv[3]);
    unsigned long bin_size = stoul(argv[4]);
    string assembly_idx_str = argv[5];
    string out_dir = argv[6];
    string contig_len_dir = argv[7];

    map<string, short> prefix_dict;
    list< map< ulli, ulli > > kmer_map_list;
    list< map< ulli, list<ulli> > > repeat_kmer_map_list;

    time_point<Clock> start = Clock::now();
    time_point<Clock> end;

    // list<string> fasta_filepaths = get_fasta_filepaths(input_file.c_str());
    get_prefix_dict(prefix_dict);

    /*unsigned long assembly_idx = start_assembly_idx;
    string assembly_idx_str;

    for(list<string>::iterator it=fasta_filepaths.begin(); it!=fasta_filepaths.end(); it++){

        kmer_map_list.clear();
        repeat_kmer_map_list.clear();

        assembly_idx_str = to_string(assembly_idx++);
    
        init_kmer_map_group_lists(kmer_map_list, repeat_kmer_map_list);
    
        // generate_binned_positional_kmers((*it), kmer_map_list, prefix_dict, kmer_len, contig_len_dir+assembly_idx_str);
        generate_binned_sampled_kmers( (*it), kmer_map_list, repeat_kmer_map_list, prefix_dict, kmer_len, sampling_count,
                                        bin_size, contig_len_dir+assembly_idx_str);
    
        // save_binned_positional_kmers(kmer_map_list, out_dir, assembly_idx_str);
        save_binned_positional_kmers(kmer_map_list, repeat_kmer_map_list, out_dir, assembly_idx_str);
    }*/

    
    init_kmer_map_group_lists(kmer_map_list, repeat_kmer_map_list);

    generate_binned_sampled_kmers( assembly_fasta, kmer_map_list, repeat_kmer_map_list, prefix_dict, kmer_len, sampling_count,
                                        bin_size, contig_len_dir+assembly_idx_str);

    save_binned_positional_kmers(kmer_map_list, repeat_kmer_map_list, out_dir, assembly_idx_str);
      

    // string remove_input_file_cmd = "rm -rf " + input_file;
    // system(remove_input_file_cmd.c_str());

    // remove(input_file.c_str());

    end = Clock::now();
    milliseconds diff = duration_cast<milliseconds>(end - start);

    // cout<<"TIME taken to generate unique kmers from "<<input_file<<" : "<<diff.count() << "ms\n";
    cout<<"TIME taken to generate unique kmers from "<<assembly_fasta<<" : "<<diff.count() << "ms\n";

    return 0;
}