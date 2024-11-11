from helpers import  read_data, rm_clusters, merge_clusters


# UPGMA algorithm implementation
def upgma(orig_dict, clusters_dict, newick_f_dict, otu_count):
    max_dist = float('inf')

    # otu_count is equivalent to the number of clusters to iterate over
    # corresponding to the individuals (OTUs) from the data set
    while otu_count > 1:
        # calculate distances between clusters
        min = max_dist
        min_i = 0
        min_j = 0
        clusters = list(clusters_dict.keys())
    
        # print(clusters)
        # print(clusters_dict)
        # print(orig_dict)
    
        for i in range(otu_count - 1):
            for j in range(i + 1, otu_count):
     
                print(f"Cluster {i}: {clusters[i]} Cluster {j}: {clusters[j]}")
                temp_dist = 0
                # determine distance between clusters as the average of all distances.
                for k in range(len(clusters[i])):
                    # print(len(clusters[j]))
                    # print(clusters[j])
                    for m in range(len(clusters[j])):
                        temp_dist += orig_dict[clusters_dict[clusters[i]][k]][clusters_dict[clusters[j]][m]]
                temp_dist /= (len(clusters[i]) * len(clusters[j]))
                if temp_dist < min:
                    min = temp_dist
                    min_i = i
                    min_j = j
                    
        # merge min_i and min_j into clusters
        cluster_i, cluster_j, merge = merge_clusters(min_i, min_j, clusters)

        merge_cluster_values(clusters_dict, min, cluster_i, cluster_j, merge)
            
        newick_f_dict[merge] = f"({newick_f_dict[cluster_i]}, {newick_f_dict[cluster_j]})"
        # print('The Newick Format dict is: ', newick_f_dict)

        # remove old clusters, decrement otu_count
        rm_clusters(clusters_dict, newick_f_dict, cluster_i, cluster_j)
        otu_count -= 1

    clusters = list(clusters_dict.keys()) 
    # print(clusters)
    # print(clusters[0])
    print("Newick Format: ", newick_f_dict[clusters[0]])


def merge_cluster_values(clusters_dict, min, cluster_i, cluster_j, merge):
    print(f"Merging {cluster_i}, {cluster_j} with distance {min}\n")
    i = 0
    for j in range(len(cluster_i)): # iterate over cluster i
        # assign the value from clusters_dict[cluster_i][j] to clusters_dict[merge][i]
        clusters_dict.setdefault(merge, {})[i] = clusters_dict[cluster_i][j]
        i += 1
    for j in range(len(cluster_j)): # iterate over cluster_j
            # similarly, assign value from clusters_dict[cluster_j][j] to clusters_dict[merge][i]
        clusters_dict[merge][i] = clusters_dict[cluster_j][j]
        i += 1


# Main
otu_count, codes, orig_dict, clusters_dict, newick_f_dict = read_data("upgma")
upgma(orig_dict, clusters_dict, newick_f_dict, otu_count)
