import sys
import numpy as np
import pandas as pd
import mykmeanssp as kmeansc


def k_means_pp(vectors, keys, k, max_iter, eps):
    np.random.seed(0)
    N, M = vectors.shape

    assert 0 < k < N and max_iter > 0 and eps >= 0
    ind_array = np.arange(N)
    result_indices = []
    centroids = []
    idx = np.random.choice(ind_array)
    result_indices.append(idx)
    centroids.append(vectors[idx])
    i = 1

    while i < k:
        D = []
        for vec in vectors:
            min_dist = np.min([np.sum((vec - x) ** 2) for x in centroids])
            D.append(min_dist)

        probs = [D[j] / sum(D) for j in range(N)]
        idx = np.random.choice(ind_array, p=probs)
        result_indices.append(idx)
        centroids.append(vectors[idx])
        i += 1

    vectors = vectors.flatten()
    centroids = kmeansc.fit(k, max_iter, eps, vectors, result_indices, M, N)

    result_keys = [str(keys[x]) for x in result_indices]
    print(",".join(result_keys))
    for j in range(k):
        vec_buffer = []
        for l in range(M):
            to_append = ("%.4f" % centroids[M * j + l])
            vec_buffer.append(to_append)
        print(",".join(vec_buffer))


def initialize_vecs(file_name_1, file_name_2):
    vectors_1 = pd.read_csv(file_name_1, header=None)
    vectors_2 = pd.read_csv(file_name_2, header=None)
    vectors = pd.merge(vectors_1, vectors_2, on=0)
    vectors = vectors.to_numpy(dtype=float)
    vectors = vectors[vectors[:, 0].argsort()]
    keys = np.array(vectors[:, 0], dtype=int)
    vectors = vectors[:, 1:]
    return vectors, keys


def check_file(file_name):
    name = file_name.split(".")
    if len(name) == 0 or (name[-1] != "txt" and name[-1] != "csv"):
        return False
    return True


def main(args):
    try:
        if len(args) == 4:
            k = int(args[0])
            max_iter = 300
            eps = float(args[1])
            file_name_1 = args[2]
            file_name_2 = args[3]

        elif len(args) == 5:
            k = int(args[0])
            max_iter = int(args[1])
            eps = float(args[2])
            file_name_1 = args[3]
            file_name_2 = args[4]

        else:
            print("Invalid Input!")
            return

        if not check_file(file_name_1) or not check_file(file_name_2):
            print("Invalid Input!")
            return

    except ValueError:
        print("Invalid Input!")
        return

    try:
        vectors, keys = initialize_vecs(file_name_1, file_name_2)

    except FileNotFoundError:
        print("Invalid Input!")
        return

    try:
        k_means_pp(vectors, keys, k, max_iter, eps)

    except AssertionError:
        print("Invalid Input!")
        return

    except Exception:
        print("An Error Has Occurred")


if __name__ == "__main__":
    main(sys.argv[1:])
