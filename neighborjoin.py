from helpers import read_data, rm_clusters, merge_clusters


def neighbor_joining(otu_count, codes, dist_mat):
    clusters_dict = init_clusters_dict(otu_count, codes, dist_mat)
    newick_f_dict = init_newick_f(codes) 
    max_dist = float('inf') # float('inf') enables a theoretical infinite upper bound

    # otu_count is equivalent to the number of clusters to iterate over
    # corresponding to each individual
    while otu_count > 2:
        # calculate r values
        clusters, r_values = calculate_r_values(clusters_dict, otu_count)

        # transformed distances calculation
        min_i, min_j = calculate_transformed_distances(clusters_dict, max_dist, otu_count, clusters, r_values)

        # merge the clusters with the min transformed distance
        cluster_i, cluster_j, merge = merge_clusters(min_i, min_j, clusters)

        # calculate branch lengths
        branch_1 = (clusters_dict[cluster_j][cluster_i] + r_values[min_i] - r_values[min_j]) / 2
        branch_2 = (clusters_dict[cluster_j][cluster_i] + r_values[min_j] - r_values[min_i]) / 2

        display_merge_info(cluster_i, cluster_j, branch_1, branch_2)

        clusters_dict[merge] = {}
        
        # Update clusters_dict with new distances
        for i in range(otu_count):
            if clusters[i] != cluster_i and clusters[i] != cluster_j:
                d1 = clusters_dict[clusters[i]][cluster_i]
                d2 = clusters_dict[clusters[i]][cluster_j]
                clusters_dict[merge][clusters[i]] = (d1 + d2 - clusters_dict[cluster_i][cluster_j]) / 2
                clusters_dict[clusters[i]][merge] = clusters_dict[merge][clusters[i]]

        # Update Newick format for the merged cluster
        newick_f_dict[merge] = f"({newick_f_dict[cluster_i]}, {newick_f_dict[cluster_j]})"

        # Remove the merged clusters from clusters_dict and newick_f_dict
        rm_clusters(clusters_dict, newick_f_dict, cluster_i, cluster_j)

        # Update list of clusters
        clusters = list(clusters_dict.keys())
        otu_count -= 1

    print("Distance btwn remaining clusters:")
    print(clusters_dict[clusters[0]][clusters[1]], clusters_dict[clusters[1]][clusters[0]])
    clusters= list(clusters_dict.keys())
    print_newick_f_dict(newick_f_dict, otu_count, clusters)


def init_newick_f(codes):
    newick_f_dict = {}
    for code in codes:
        newick_f_dict[code] = code
    return newick_f_dict


def init_clusters_dict(otu_count, codes, dist_mat):
    clusters_dict = {}
    for j, code in enumerate(codes):
        # codes are used as keys
        clusters_dict[code] = {}
        for i in range(otu_count):
            clusters_dict[code][codes[i]] = dist_mat[j][i]
    return clusters_dict


def print_newick_f_dict(newick_f_dict, otu_count, clusters):
    # maybe the worst code I have ever written
    print("\nNewick Format: (", end="")
    for j in range(otu_count):
        if j < otu_count - 1:
            print(f"{newick_f_dict[clusters[j]]}", sep="", end=", ")
        else:
            print(f"{newick_f_dict[clusters[j]]}", sep="", end="")
    print(')')


# prints out which clusters to merge and the distances
def display_merge_info(cluster_i, cluster_j, branch_1, branch_2):
    print(f"\nMerge {cluster_i}, {cluster_j}")
    print(f"Distance btwn {cluster_i} and ancestral node {branch_1}")
    print(f"Distance btwn {cluster_j} and ancestral node {branch_2}")


def calculate_transformed_distances(clusters_dict, max_dist, otu_count, clusters, r_values):
    min = max_dist
    min_i = 0
    min_j = 0

    transformed_distances = {}
    print("Transformed Distances:")
    for i in range(otu_count - 1):
        for j in range(i + 1, otu_count):
            if (i, j) not in transformed_distances:
                # calculate transformed distance
                transformed_distances[(i, j)] = clusters_dict[clusters[i]][clusters[j]] - r_values[i] - r_values[j]
            temp_td = transformed_distances[(i, j)]
            print(f"TransformedDist({clusters[i]}, {clusters[j]}) = {temp_td}", end=", ")

            if temp_td < min:
                min = temp_td
                min_i = i
                min_j = j

        print()
    return min_i, min_j


def calculate_r_values(clusters_dict, otu_count):
    clusters = list(clusters_dict.keys())
    r_values = []

    for i in range(otu_count):
        temp = 0
        for j in range(otu_count):
            temp += clusters_dict[clusters[i]].get(clusters[j], 0)
        r_values.append(temp / (otu_count - 2))

    print("\nR Values:")
    for i in range(otu_count):
        print(f"{clusters[i]} {r_values[i]}")
    print()
    return clusters,r_values

otu_count, codes, dist_mat = read_data("neighborjoin")


# Run neighbor-joining algorithm
neighbor_joining(otu_count, codes, dist_mat)