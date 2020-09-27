import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import Counter
import glob
import os, time
from collections import Counter
import math
import shlex, subprocess
from io import StringIO
from queue import Queue
from multiprocessing import Pool
import random
from scipy.spatial.distance import hamming
import pickle
import timeit
import time
from datetime import datetime
import argparse, sys
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA, TruncatedSVD
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import squareform, hamming
import psutil
import scipy
from scipy.sparse import csr_matrix, dok_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from joblib import Parallel, delayed, parallel_backend
import seaborn as sns


def read_file(filepath):
    f = open(filepath, 'r')
    lines = f.readlines()
    f.close()
    return lines


def write_file(filename, out_str):
    f = open(filename, 'w+')
    f.write(out_str)
    f.close()


def run_cpp_binaries(binary_file, *args):
    cmd = binary_file
    for arg in args:
        cmd += ' {}'.format(arg)
    print('run_cpp_binaries: ', cmd)
    try:
        p = subprocess.Popen(shlex.split(cmd))#, stdout=open(mum_results_file, 'w'))
        # print('waiting...', cmd)
        p.wait()
        print(cmd, ' output obtained')
    except Exception as e:
        print(cmd, e)


def create_dir(directory_name):
    cmd = "mkdir {}".format(directory_name)
    p = subprocess.Popen(shlex.split(cmd))
    p.wait()
    print(directory_name, " created")


def remove_file(filename):
    cmd = "rm -rf {}".format(filename)
    try:
        p = subprocess.Popen(shlex.split(cmd), shell=True)
        p.wait()
    except Exception as e:
        print("remove_file error", cmd, e)


def cpp_dir_input_setup(assemblies, fasta_file_list, ncores, contig_feature_partitions, min_presence_count):
    outdir = "SIMILE_results"
    if(os.path.exists(outdir)):
        outdir += "_{}".format(datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
    outdir += '/'
    create_dir(outdir)
    create_dir("{}binned_kmers/".format(outdir))
    create_dir("{}retained_binned_kmers_{}/".format(outdir, min_presence_count))
    create_dir("{}retained_binned_kmers_{}/assemblywise/".format(outdir, min_presence_count))
    create_dir("{}contig_lengths/".format(outdir))
    create_dir("{}shared_contigs/".format(outdir))


    assembly_count = len(assemblies)
    step = math.ceil(assembly_count/(float)(ncores))
    assembly_partition_pos_arr = np.arange(step, assembly_count, step)
    core_idx = 0
    partition_start = 0

    for partition_end in assembly_partition_pos_arr:
        opstr = '\n'.join(fasta_file_list[partition_start:partition_end])+'\n'
        write_file("{}input_binned_assemblies_{}.txt".format(outdir, core_idx), opstr)
        core_idx += 1
        partition_start = partition_end
    if(partition_start<assembly_count):
        opstr = '\n'.join(fasta_file_list[partition_start:assembly_count])+'\n'
        write_file("{}input_binned_assemblies_{}.txt".format(outdir, core_idx), opstr)

    opstr = '\n'.join(fasta_file_list) + '\n'
    write_file("{}all_assembly_filepaths.txt".format(outdir), opstr)


    contig_feature_step = math.ceil(assembly_count/(float)(contig_feature_partitions))
    contig_feature_partition_pos_list = [0]
    core_idx = 0
    partition_start = 0
    partition_end = contig_feature_step
    
    while(partition_end < assembly_count):
        contig_feature_partition_pos_list.append(partition_end)
        partition_start = partition_end
        partition_end += contig_feature_step
    contig_feature_partition_pos_list.append(assembly_count)

    return outdir, np.concatenate([[0],assembly_partition_pos_arr]), np.array(contig_feature_partition_pos_list)


def presence_vector_files(hash_string, hash_string_idx, in_dir, outdir='.'):
    cmd = "ls {}{}* > {}{}_group_files".format(in_dir, hash_string, outdir, hash_string_idx)
    print("presence_vector_files : ", cmd)
    try:
        p = subprocess.Popen(cmd, shell=True)
        p.wait()
        print("presence_vector_files created:", hash_string)
    except Exception as e:
        print("presence_vector_files :", cmd, e)


def combine_filemaps(in_dir, outdir='.'):
    cmd = "echo {}filemap_* | xargs cat > {}combined_map.csv".format(in_dir, outdir)
    print("combine_filemaps : ", cmd)
    try:
        p = subprocess.Popen(cmd, shell=True)
        p.wait()
        print("combine_filemaps combined")
    except Exception as e:
        print("combine_filemaps error :", cmd, e)


def get_intermediate_filemaps(filenames_str, outfilename):
    cmd = "cat {} > {}".format(filenames_str, outfilename)
    rm_cmd = "rm -f {}".format(filenames_str)
#     print("get_intermediate_filemaps : ", cmd)
    print("get_intermediate_filemaps : ", outfilename)
    try:
        p = subprocess.Popen(cmd, shell=True)
        p.wait()
        print("get_intermediate_filemaps combined")
        p = subprocess.Popen(rm_cmd, shell=True)
        p.wait()
    except Exception as e:
        print("get_intermediate_filemaps error :", cmd, e)


def combine_intermediate_filemaps(in_dir, outdir='.'):
    cmd = "echo {}intermediate_combined_* | xargs cat > {}final_combined_map.csv".format(in_dir, outdir)
    rm_cmd = "echo {}intermediate_combined_* | xargs rm -rf".format(in_dir)
#     print("combine_filemaps : ", cmd)
    print("combine_filemaps : ", outdir)
    try:
        p = subprocess.Popen(cmd, shell=True)
        p.wait()
        print("combine_intermediate_filemaps combined")
        p = subprocess.Popen(rm_cmd, shell=True)
        p.wait()
    except Exception as e:
        print("combine_intermediate_filemaps error :", cmd, e)


def get_hamming_distance(df, assembly_count, start_idx=1, end_idx=-1):
    nrows = df.shape[0]
    hamming_dist_mat = np.zeros((assembly_count, assembly_count))
    
    for i in range(start_idx-1, end_idx):
        for j in range(start_idx, assembly_count):
            pos = np.logical_xor(df[df.columns[i]], df[df.columns[j]])
            val = np.sum(df.blocks_count[pos])
            hamming_dist_mat[i][j] = val
            hamming_dist_mat[j][i] = val
    print("get_hamming_distance", start_idx, end_idx)
    # print(hamming_dist_mat)
    return hamming_dist_mat


def get_weighted_hamming_distance(df, assembly_count, start_idx=1, end_idx=-1):
    nrows = df.shape[0]
    hamming_dist_mat = np.zeros((assembly_count, assembly_count))
    
    for i in range(start_idx-1, end_idx):
        for j in range(start_idx, assembly_count):
            pos = np.logical_xor(df[df.columns[i]], df[df.columns[j]])
            val = np.sum(df.total_blocks_length[pos])
            hamming_dist_mat[i][j] = val
            hamming_dist_mat[j][i] = val
    print("get_weighted_hamming_distance", start_idx, end_idx)
    # print(hamming_dist_mat)
    return hamming_dist_mat


def generate_column_vector(filename, col_idx, outfilename):
    cmd = "cut -d, -f{} {} > {}".format(col_idx, filename, outfilename)
    print("generate_column_vector started:", cmd)
    try:
        p = subprocess.Popen(cmd, shell=True)
        p.wait()
        print("generate_column_vector ended:", col_idx)
    except Exception as e:
        print("generate_column_vector error:", cmd, e)


def generate_start_indices(filename, col_idx, outfilename):
    cmd = "cut -d, -f{} {} | awk '{{total += $0; $0 = total - $0}}1' > {}".format(col_idx, filename, outfilename)
    print("generate_start_indices started:", cmd)
    try:
        p = subprocess.Popen(cmd, shell=True)
        p.wait()
        print("generate_start_indices ended:", col_idx)
    except Exception as e:
        print("generate_start_indices error:", cmd, e)


def get_total_block_count(filename, col_idx):
    cmd = "cut -d, -f{} {} | awk '{{ sum += $1 }} END {{ print sum }}'".format(col_idx, filename)
    print("get_total_block_count: ", cmd)
    try:
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, error = p.communicate()
        print(out, error)
        return out.decode("utf-8").strip()
    except Exception as e:
        print("get_total_block_count ERROR: ", cmd, e)


def retain_presence_vectors(filename, col_idx_str, outfilename):
    cmd = "cut -d, -f{} --complement {} > {}".format(col_idx_str, filename, outfilename)
    print("retain_presence_vectors started:", cmd)
    try:
        p = subprocess.Popen(cmd, shell=True)
        p.wait()
        print("retain_presence_vectors ended:", col_idx_str)
    except Exception as e:
        print("retain_presence_vectors error:", cmd, e)


def get_pricipal_components_count(model, min_explained_variance=0.8):
    cumsum = 0
    for idx, val in enumerate(model.explained_variance_ratio_):
        cumsum += val
        if(cumsum>=min_explained_variance):
            print(idx, cumsum)
            return idx
            break


def get_kmeans_elbow_point(dim_reduced_data, linearity_limiting_threshold=1.1):
    previous_slope = np.inf
    previous_sum_of_sq = np.inf
    k = 0
    while(True):
        k += 1
        print(k)
        km = KMeans(n_clusters=k)
        km = km.fit(X)
        current_sum_of_sq = km.inertia_
        if(k==1):
            previous_sum_of_sq = current_sum_of_sq
            continue
        else:
            current_slope = previous_sum_of_sq - current_sum_of_sq
            if(previous_slope <= linearity_limiting_threshold*current_slope):
                return k-2
            previous_slope = current_slope
            previous_sum_of_sq = current_sum_of_sq


def kmeans_wrapper(X, dim_reduce_model, min_explained_variance=0.8, linearity_limiting_threshold=1.1):
    
    dim_reduce_model.fit(X)
    
    n_principal_components = get_pricipal_components_count(dim_reduce_model, min_explained_variance)
    X_dim_reduce = dim_reduce_model.transform(X)
    n_clusters = get_kmeans_elbow_point(X_dim_reduce[:,:n_principal_components], linearity_limiting_threshold)
    
    kmeans = KMeans(n_clusters=n_clusters)
    kmeans = kmeans.fit(X_dim_reduce[:,:n_principal_components])
    
    return kmeans, n_clusters


def generate_hamming_distance_dendrogram(hamming_distance_mat, labels, save_plot=True, reference='', outdir='',
                                         show_plot=False, figsize=(12,60), leaf_font_size=9, logtransform=False):
    if(logtransform):
        hamming_distance_mat = np.log(hamming_distance_mat+np.ones(hamming_distance_mat.shape))
    Z = linkage(hamming_distance_mat)
    fig, ax = plt.subplots(figsize=figsize)
    dendrogram(Z, labels=labels, orientation='left', ax=ax, leaf_font_size=leaf_font_size)
    colored_labels = ax.get_ymajorticklabels()
    for label in colored_labels:
        text = label.get_text()
        if(':R' in text):
            color = 'red'
        elif(':S' in text):
            color = 'blue'
        else:
            color = 'black'
        label.set_color(color)
    plt.tight_layout()
    if(show_plot):
        plt.show()
    if(save_plot):
        outfile = '{}hamming_distance_{}'.format(outdir, len(labels))
        if(len(reference)>0):
            outfile += '_{}'.format(reference)
        if(logtransform):
            outfile += '_logtransformed'
        outfile += '.png'
        fig.savefig(outfile, dpi=fig.dpi, bbox_inches="tight")


def str2bool(v):
    return v.lower() in ("y", "yes", "true", "t", "1")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='SIMILE: Discover shared genomic regions across a collection of '+
                                     'metagenomic samples')
    parser.add_argument('-i', '--input', required=True, help="Input csv file")
    parser.add_argument('-n', '--ncores', type=int, nargs='?', default=8,
                        help="Number of cores to be used (default: 8)")
    parser.add_argument('-K', '--kmer_len', type=int, nargs='?', default=25,
                        help="Length of kmers (default: 25)")
    parser.add_argument('-p', '--min_perc', type=float, nargs='?', default=5.0,
                        help="Minimum %% of samples that should contain the shared region (default: 5.0)")
    parser.add_argument('-f', '--sampling_frequency', type=float, nargs='?', default=1.0,
                        help="Proportion of k-mers to be sampled and selected for comparisons (default 1.0)")
    parser.add_argument('-b', '--bin_size', type=int, nargs='?', default=1000,
                        help="Size of bins for sampling the k-mers (default: 1000)")
    parser.add_argument('-c', '--sampled_kmer_contig_presence', type=float, nargs='?', default=40.0,
                        help="Minimum %% of the sampled k-mers that should have shared presence for the contig to be "+
                            "considered as a shared contig (default: 40.0)")
    parser.add_argument('-m', '--mem', nargs='?', default="36000MB", help="Upper limit for RAM memory usage.  " +
                        "Can be in mb/MB/gb/GB/tb/TB (case insensitive), default unit is MB. (default: 36000MB)")
    
    args = parser.parse_args()

    ncores = args.ncores #8 # 5

    input_pd = pd.read_csv(args.input, header=None)
    input_assemblies = input_pd[input_pd.columns[0]].values
    fasta_filepaths = input_pd[input_pd.columns[1]].values


    kmer_len = args.kmer_len
    prefix_len = 5 ## Fixed kmer prefix length
    min_presence_perc = args.min_perc
    prefix_count = pow(4, prefix_len)
    min_presence_count = math.ceil(len(input_assemblies)*min_presence_perc/100.0)
    if(min_presence_count < 2):
        min_presence_count = 2
    assembly_count = len(input_assemblies)
    # feature_partitions = max(ncores, math.ceil(assembly_count/10))

    mem = args.mem
    if(mem.isdigit()):
        available_memory = int(mem)
    else:
        assert(len(mem)>2)
        val = mem[:-2]
        units = mem[-2:]
        assert(val.isdigit())
        val = int(val)
        assert(units[0] in set(['m', 'M', 'g', 'G', 't', 'T']) and units[1] in set(['b', 'B']))
        # assert(units[0].isalpha() and units[1].isalpha())
        units = units.upper()
        if(units[0]=='M'):
            available_memory = val
        elif(units[0]=='G'):
            available_memory = val*1000
        else:
            available_memory = val*1000000

    max_assemblies_per_partition = available_memory/ncores - 2500
    if(max_assemblies_per_partition < 10):
        max_assemblies_per_partition = 10
#     else:
#         max_assemblies_per_partition = math.floor(math.sqrt(max_assemblies_per_partition/genome_len_mb))
#         max_assemblies_per_partition = max(10, max_assemblies_per_partition)

    feature_partitions = max(ncores, math.ceil(assembly_count/max_assemblies_per_partition))

    assemblies_per_partition = math.ceil(assembly_count/feature_partitions)
    sampling_frequency = args.sampling_frequency
    bin_size = args.bin_size
    sampling_count = max(math.ceil(sampling_frequency*bin_size/100.0), 2)
    contig_feature_partitions = max(feature_partitions, math.ceil(assembly_count/4.0))
    contig_presence_threshold = args.sampled_kmer_contig_presence
    seq2write_batch = 1000

    print(ncores, input_pd.shape, kmer_len, min_presence_perc, min_presence_count, assembly_count, sampling_frequency,
          bin_size, sampling_count, feature_partitions, contig_feature_partitions, contig_presence_threshold, available_memory)
            

    results_dir, binned_assembly_arr, contig_feature_partition_pos_arr = cpp_dir_input_setup(input_assemblies, fasta_filepaths, ncores,
                                                                                            feature_partitions, min_presence_count)
    print(binned_assembly_arr)
    print(contig_feature_partition_pos_arr)

    out_str = "#samples: {}\n".format(len(input_assemblies))
    out_str += "kmer_len: {}\n".format(kmer_len)
    out_str += "min_perc: {}\n".format(min_presence_perc)
    out_str += "#cores: {}\n".format(ncores)
    out_str += "min_presence_count: {}\n".format(min_presence_count)
    out_str += "available_memory: {}\n".format(available_memory)
    out_str += "sampling_frequency: {}\n".format(sampling_frequency)
    out_str += "bin_size: {}\n".format(bin_size)
    out_str += "contig_presence_threshold: {}\n".format(contig_presence_threshold)
    # out_str += "genome_len: {}\n".format(genome_len)

    write_file('{}params.txt'.format(results_dir), out_str)

    # start = time.time()
    # Parallel(n_jobs=ncores, prefer="threads")(
    #     delayed(run_cpp_binaries)(  "./kmer_binned_sampling.o", "{}input_binned_assemblies_{}.txt".format(results_dir, core_idx),
    #                                 kmer_len, sampling_count, bin_size, binned_assembly_arr[core_idx], "{}binned_kmers/".format(results_dir),
    #                                 "{}contig_lengths/".format(results_dir))
    #         for core_idx in range(len(binned_assembly_arr)))
    # end = time.time() # timeit.timeit()
    # print('TIME taken to sample and bin kmers from all files:', end-start)

    start = time.time()
    Parallel(n_jobs=ncores, prefer="threads")(
        delayed(run_cpp_binaries)(  "./kmer_binned_sampling.o", assembly_fasta, kmer_len, sampling_count, bin_size, assembly_idx,
                                    "{}binned_kmers/".format(results_dir), "{}contig_lengths/".format(results_dir))
            for assembly_idx, assembly_fasta in enumerate(fasta_filepaths))
    end = time.time() # timeit.timeit()
    print('TIME taken to sample and bin kmers from all files:', end-start)


    start = time.time()
    Parallel(n_jobs=ncores, prefer="threads")(
        delayed(run_cpp_binaries)(  "./sampled_kmer_filtering.o", assembly_count, prefix_idx, "{}binned_kmers/".format(results_dir),
                                    min_presence_count, "{}retained_binned_kmers_{}/".format(results_dir, min_presence_count))
            for prefix_idx in range(prefix_count))
    end = time.time() # timeit.timeit()
    print('TIME taken to get kmers present in at least {} perc of input files:'.format(min_presence_perc), end-start)


    start = time.time()
    Parallel(n_jobs=ncores, prefer="threads")(
        delayed(run_cpp_binaries)(  "./kmer_bin_presence_generator.o", assembly_count, prefix_idx, "{}binned_kmers/".format(results_dir),
                                    "{}retained_binned_kmers_{}/".format(results_dir, min_presence_count),
                                    "{}retained_binned_kmers_{}/assemblywise/".format(results_dir, min_presence_count),
                                    "{}contig_lengths/".format(results_dir))
            for prefix_idx in range(prefix_count))
    end = time.time() # timeit.timeit()
    print('TIME taken to get filtered prefix binned kmers\' assemblywise bin presence:', end-start)


    # contig_presence_threshold = 40.0
    start = time.time()
    Parallel(n_jobs=ncores, prefer="threads")(
        delayed(run_cpp_binaries)("./shared_contig_extractor.o", assembly_idx, prefix_count, contig_presence_threshold,
                                  sampling_count, "{}retained_binned_kmers_{}/assemblywise/".format(results_dir, min_presence_count),
                                  "{}contig_lengths/".format(results_dir), fasta_filepaths[assembly_idx],
                                  "{}shared_contigs/{}_shared.fasta".format(results_dir, input_assemblies[assembly_idx]))
            for assembly_idx in range(assembly_count))
    end = time.time() # timeit.timeit()
    print('TIME taken to fetch shared contigs:', end-start)

