from helpers import rm_clusters, user_prompt



def read_data():
    orig_dict = {} 
    clusters_dict = {}
    newick_format = {} 
    codes_dict = {}
    while True:
        user_in = input(user_prompt)
        try:
            with open(user_in, 'r') as file:
                # operational taxonomic units count on first line
                otu_count = int(file.readline().strip())
                print(otu_count)

                # codes are on second line
                # a code is simply the name of an individiual
                codes = [code[0] for code in file.readline().strip().split()]
                # codes = file.readline().strip().split()
        
                # read in actual distance matrix from input file
                # populate 
                dist_mat = []
                for i in range(otu_count):
                    codes_dict[codes[i]] = i
                    newick_format[codes[i]] = codes[i]
                    clusters_dict[codes[i]] = {0: codes[i]}
                    dist_mat.append(list(map(float, file.readline().strip().split())))
                print("\n", " ".join(codes))
    

            for i in range(otu_count):
                print(f"{codes[i]} ", end="")
                for j in range(otu_count):
                    print(f"{dist_mat[i][j]} ", end="")
                    orig_dict.setdefault(codes[i], {})[codes[j]] = dist_mat[i][j]
                print()
            print()

            # bad code
            return otu_count, codes, orig_dict, clusters_dict, newick_format

        except FileNotFoundError:
            print(f"Error opening input file with name {user_in}, try again.")
            continue
        except Exception as e:
            print(f"Failure reading data from file: {e}")
            exit(1)


# UPGMA algorithm implementation
def upgma(orig_dict, clusters_dict, newick_format, otu_count):
    max_dist = float('inf')
    num_clusters = otu_count
    while num_clusters > 1:

        # calculate distances between clusters
        min = max_dist
        min_i = 0
        min_j = 0
        clusters = list(clusters_dict.keys())
    
        # print(clusters)
        # print(clusters_dict)
        # print(orig_dict)
    
        for i in range(num_clusters - 1):
            for j in range(i + 1, num_clusters):
     
                print(f"Cluster {i}: {clusters[i]} Cluster {j}: {clusters[j]}")
                temp_dist = 0
                # determine distance between clusters as the unweighted average of all distances.
                
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

        # merge min_i and min_j into clusters and add to hash tables
        cluster_i = clusters[min_i]
        cluster_j = clusters[min_j]
        merge = cluster_i + cluster_j


        merge_cluster_values(clusters_dict, min, cluster_i, cluster_j, merge)
            
        newick_format[merge] = f"({newick_format[cluster_i]}, {newick_format[cluster_j]})"
        # print('The Newick Format dict is: ', newick_format)

        # eliminate old clusters, decrement num_clusters
        rm_clusters(clusters_dict, newick_format, cluster_i, cluster_j)
        num_clusters -= 1

    clusters = list(clusters_dict.keys()) 
    # print(clusters)
    # print(clusters[0])
    print("Newick Format: ", newick_format[clusters[0]])


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
otu_count, codes, orig_dict, clusters_dict, newick_format = read_data()
upgma(orig_dict, clusters_dict, newick_format, otu_count)
