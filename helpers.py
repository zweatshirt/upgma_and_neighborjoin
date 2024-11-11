user_prompt = "Input distance matrix file name (make sure the file is in the project directory, or specify the path): "
UPGMA = "upgma"
NJ = "neighborjoin"

# helper function to remove clusters from dictionaries
def rm_clusters(clusters, newick_format, cluster_i, cluster_j):
    clusters.pop(cluster_i, None)
    clusters.pop(cluster_j, None)
    newick_format.pop(cluster_i, None)
    newick_format.pop(cluster_j, None)


def merge_clusters(min_i, min_j, clusters):
    cluster_i = clusters[min_i]
    cluster_j = clusters[min_j]
    merge = cluster_i + cluster_j
    return cluster_i,cluster_j,merge


# alg type must be either nj or upgma
def read_data(alg_type):
    if alg_type == UPGMA:
        return read_upgma_data()
    elif alg_type == NJ:
        return read_nj_data()


def read_upgma_data():
    orig_dict = {} 
    clusters_dict = {}
    newick_f_dict = {} 
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
                    newick_f_dict[codes[i]] = codes[i]
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
            return otu_count, codes, orig_dict, clusters_dict, newick_f_dict

        except FileNotFoundError:
            print(f"Error opening input file with name {user_in}, try again.")
            continue
        except Exception as e:
            print(f"Failure reading data from file: {e}")
            exit(1)


def read_nj_data():
    while True:
        try:# Get input file
            user_in = input(user_prompt)
            with open(user_in, 'r') as file:

                # OTUs are the individuals in the data set. 
                # OTU = Operational Taxonomic Unit
                otu_count = int(file.readline().strip())

                # codes are on second line
                # a code is simply the name of an individiual
                codes = [code[0] for code in file.readline().strip().split()]

                # read in actual distance matrix from input file
                dist_mat = []
                for line in file:
                    dist_mat.append(list(map(float, line.strip().split())))
            # Print distance matrix
            print_distance_matrix(otu_count, codes, dist_mat)
                
            return otu_count, codes, dist_mat
        except FileNotFoundError:
            print(f"Error opening input file with name {user_in}")
            continue
        except Exception as e:
            print(f"Failure reading data from file: {e}")
            exit(1)


def print_distance_matrix(otu_count, codes, dist_mat):
    print("\n    ", " ".join(codes))
    for i in range(otu_count):
        print(f"{codes[i]} ", " ".join(map(str, dist_mat[i])))

