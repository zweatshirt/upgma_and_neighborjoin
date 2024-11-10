user_prompt = "Input distance matrix file name (make sure the file is in the project directory, or specify the path): "

# helper function to remove clusters from dictionaries
def rm_clusters(clusters, newick_format, cluster_i, cluster_j):
    clusters.pop(cluster_i, None)
    clusters.pop(cluster_j, None)
    newick_format.pop(cluster_i, None)
    newick_format.pop(cluster_j, None)