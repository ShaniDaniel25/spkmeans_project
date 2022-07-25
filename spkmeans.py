import sys
import numpy as np
import spkmeans as spkmeans
import pandas as pd
EPS = 0.001
MAX_ITER = 300

def spkmeans_py(k, goal, vectors):
    k_zero = 1 if k == 0 else 0
    goals = {"wam", "ddg", "lnorm", "jacobi"}

    if goal == "spk":
        c_out = spkmeans.execute_spk(vectors, goal, k_zero)
        final_k = c_out[-1] if k_zero == 1 else k
        n = (len(c_out) - 1) / final_k if k_zero== 1 else len(c_out) / final_k
        T = np.ndarray(n, final_k)
        for i in range(n):
            for j in range(final_k):
                T[i][j] = c_out[final_k * i + j]
                
        indices = init_centroids(T, k)
        result_centroids = spkmeans.kmeans(final_k, MAX_ITER, EPS, T, indices, n)
        
        print(",".join(indices))
        for i in range(final_k):
            vec_buffer = []
            for j in range(final_k):
                to_append = ("%.4f" % result_centroids[final_k * i + j])
                vec_buffer.append(to_append)
            print(",".join(vec_buffer))

    elif goal in goals:
        c_out = spkmeans.execute_spk(vectors, goal, 0)
        n = len(vectors)
        mat = np.ndarray(n, n)
        for i in range(n):
            for j in range(n):
                mat[i][j] = c_out[n * 1 + j]
        
        if goal == 'jacobi':
            eigen_vals = []
            for i in range(n):
                eigen_vals.append(c_out[n*n + i])
            print(",".join(eigen_vals))
            
        for i in range(n):
            vec_buffer = []
            for j in range(n):
                to_append = ("%.4f" % mat[n * i + j])
                vec_buffer.append(to_append)
            print(",".join(vec_buffer))
        
    else:
        print("Invalid Input!")


def init_centroids(vectors, k):
    np.random.seed(0)
    N, M = vectors.shape
    assert 0 < k < N

    indices = [i for i in range(N)]
    ind_array = np.array(indices)
    result_indices = []
    centroids = []
    idx = np.random.choice(indices)
    result_indices.append(idx)
    centroids.append(vectors[idx])
    i = 1
    while i < k:
        probs = []
        D = []
        for vec in vectors:
            min_dist = np.sum((vec - centroids[0])**2)
            for j in range(1, i):
                dist = np.sum((vec - centroids[j])**2)
                if dist < min_dist:
                    min_dist = dist
            D.append(min_dist)

        for l in range(N):
            prob = D[l] / sum(D)
            probs.append(prob)

        i += 1
        idx = np.random.choice(ind_array, p=probs)
        result_indices.append(idx)
        centroids.append(vectors[idx])

    return result_indices

def init_vecs(file_name):
    vectors = pd.read_csv(file_name, header=None)
    vectors = vectors.to_numpy(dtype=float)
    vectors = vectors[vectors[:, 0].argsort()][:, 1:]
    return vectors


def check_file(file_name):
    name = file_name.split(".")
    if len(name) == 0 or (name[-1] != "txt" and name[-1] != "csv"):
        return False
    return True

def main(args):
    try:
        assert len(args) == 3 and args[0].isdigit() and check_file(args[2])
        k = int(args[0])
        goal = args[1]
        file_name = args[2]
        
        vectors = init_vecs(file_name)
        spkmeans_py(k, goal, vectors)
        
        
    except AssertionError:
        print("Invalid Input!")
        return
    
    except Exception:
        print("An Error Has Occurred")
        return

if __name__ == "__main__":
    main(sys.argv[1:])

        
        
    
    
    

