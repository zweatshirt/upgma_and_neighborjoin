from helpers import rm_clusters, user_prompt

def read_data():
    while True:
        try:# Get input file
            user_in = input(user_prompt)
            with open(user_in, 'r') as file:
                otu_count = int(file.readline().strip())

                # codes are on second line
                # a code is simply the name of an individiual
                codes = [code[0] for code in file.readline().strip().split()]

                # read in actual distance matrix from input file
                dist_mat = []
                for line in file:
                    dist_mat.append(list(map(float, line.strip().split())))
                
            return otu_count, codes, dist_mat
        except FileNotFoundError:
            print(f"Error opening input file with name {user_in}")
            continue
        except Exception as e:
            print(f"Failure reading data from file: {e}")
            exit(1)


def neighbor_joining(otu_count, codes, dist_mat):
    clusters_dict = {}
    for j, code in enumerate(codes):
        clusters_dict[code] = {}
        for i in range(otu_count):
            clusters_dict[code][codes[i]] = dist_mat[j][i]

    # Initialize the Newick format dictionary with each code as its own cluster
    newick_format = {}
    for code in codes:
        newick_format[code] = code

    max_dist = float('inf')
    num_clusters = otu_count  # initializing the number of clusters

    while num_clusters > 2:
        # Step 1: Calculate r values
        clusters = list(clusters_dict.keys())
        array_r_values = []

        for i in range(num_clusters):
            temp = 0
            for j in range(num_clusters):
                temp += clusters_dict[clusters[i]].get(clusters[j], 0)
            array_r_values.append(temp / (num_clusters - 2))  # with divide by 2 modification due to double counting

        # print("R Values")
        # for i in range(num_clusters):
        #     print(f"{clusters[i]} {array_r_values[i]}")

        # Step 2: Calculate transformed distances (TD)
        min = max_dist
        min_i = 0
        min_j = 0

        print("Transformed Distances Matrix:")
        for i in range(num_clusters - 1):
            for j in range(i + 1, num_clusters):
                temp_td = clusters_dict[clusters[i]][clusters[j]] - array_r_values[i] - array_r_values[j]
                print(f"TransformedDist({clusters[i]},{clusters[j]}) = {temp_td}", end=", ")

                if temp_td < min:
                    min = temp_td
                    min_i = i
                    min_j = j

            print()

        # Step 3: Merge the clusters with the min transformed distance
        cluster_i = clusters[min_i]
        cluster_j = clusters[min_j]
        merge = cluster_i + cluster_j

        # Calculate branch lengths
        branch_1 = (clusters_dict[cluster_j][cluster_i] + array_r_values[min_i] - array_r_values[min_j]) / 2
        branch_2 = (clusters_dict[cluster_j][cluster_i] + array_r_values[min_j] - array_r_values[min_i]) / 2


        print(f"\nMerge {cluster_i}, {cluster_j}")
        print(f"Distance btwn {cluster_i} and ancestral node {branch_1}")
        print(f"Distance btwn {cluster_j} and ancestral node {branch_2}")

        # Initialize the merged cluster in clusters_dict
        clusters_dict[merge] = {}
        
        # Update clusters_dict with new distances
        for i in range(num_clusters):
            if clusters[i] != cluster_i and clusters[i] != cluster_j:
                d1 = clusters_dict[clusters[i]][cluster_i]
                d2 = clusters_dict[clusters[i]][cluster_j]
                clusters_dict[merge][clusters[i]] = (d1 + d2 - clusters_dict[cluster_i][cluster_j]) / 2
                clusters_dict[clusters[i]][merge] = clusters_dict[merge][clusters[i]]

        # Update Newick format for the merged cluster
        newick_format[merge] = f"({newick_format[cluster_i]}, {newick_format[cluster_j]})"

        # Remove the merged clusters from clusters_dict and newick_format
        rm_clusters(clusters_dict, newick_format, cluster_i, cluster_j)

        # Update list of clusters
        clusters = list(clusters_dict.keys())
        num_clusters -= 1

    # Final results
    print("Distance btwn remaining clusters:")
    print(clusters_dict[clusters[0]][clusters[1]], clusters_dict[clusters[1]][clusters[0]])

    clusters= list(clusters_dict.keys())
 
    # maybe the worst code I have ever written
    print("\nNewick Format: (", end="")
    for j in range(num_clusters):
        print(f"{newick_format[clusters[j]]}", sep=" ", end="")
    print(')')

otu_count, codes, dist_mat = read_data()

# Print distance matrix
print("\n    ", " ".join(codes))
for i in range(otu_count):
    print(f"{codes[i]} ", " ".join(map(str, dist_mat[i])))

# Run neighbor-joining algorithm
neighbor_joining(otu_count, codes, dist_mat)