from os import listdir
from os.path import join


"""
Collects curves in files outputted by diffstat and compares them. The file structure is as follows:

    /1/
        curves_...
        curves_...
        ...
    /2/
        ...
    ...

The numbers 1, 2, ... represent the iteration, and the curves_... files contained are the output of diffstat for
various calls for that iteration. It is required that all curve files contain the same amount of curves.

On topic of comparing, specifically two things are done:

    1. For each iteration, it is printed how many curves are duplicated. This excludes any curves generated with
       coefficient and prime moduli range between 1 and 10. (Likely many curves will be duplicated here because of the
       small amount of possibilities.)
    2. For each iteration, it is printed how many (unique) curves it contains is shared with the (unique) curves of
       every other iteration. Specifically, if the set of curves of iteration i is A, and the set of curves of iteration
       j is B, then len(A intersect B) is printed.
"""


def _read_curves(file_name):
    """(str) -> list"""
    data = lambda line: tuple(int(i) for i in line.split())

    curves = []
    with open(file_name, "r") as file:
        for curve_data in file:
            e = data(curve_data)
            curves.append(e)

    return curves


def collect_curves(iterations_path, iterations_dirs):
    """(str, list) -> dict
    Returns a dict with keys as the iteration number and values as the recorded curves in a set.

    For each iteration i, all curves are collected from the recorded curves for that iteration into set's Si (so only
    unique ones are taken). The returned dict will have the key, value pairs i, Si. If a curve record for coefficients
    and prime moduli in the range 1 to 10 exists, it is skipped (see this module's description above).

    Assumes iterations_path has the directory structure explained in this module's description. It is not necessary that
    only curve records as outputted by diffstat are contained within the iteration folders (these will be skipped). It
    is assumed that all curve files in all iterations have the (globally) same amount of curves.
    """
    iterations = {int(i): set() for i in iterations_dirs}
    for iter in iterations:
        files = listdir(join(iterations_path, str(iter)))
        for file_name in files:
            # The endswith is required, at least on UNIX machines, to avoid .txt~ files (hidden files)?
            if file_name.startswith("curves") and ("curves_1_10_" not in file_name) and file_name.endswith(".txt"):
                curves = _read_curves(join(iterations_path, str(iter), file_name))
                iterations[iter] = iterations[iter].union(set(curves))

    return iterations


def pairwise_duplicates(iterations, iter_indices):
    """(dict, list) -> dict
    Returns the shared curves between each iteration in iterations.

    Assumes iter_indices is a sorted list of all the iteration numbers for iterations. Let i in iter_indices. For each
    j in iter_indices such that j > i, there is Sij equal to the intersection of the set of all curves in iteration i
    and the set of all curves in iteration j. The returned dict has key, value pairs i, Li, where Li is a dict with key,
    value pairs j, Sij.

    Assumes iterations is output of collect_curves().
    """
    similarities = {i: dict() for i in iter_indices}
    for i in iter_indices:
        # Iterations to compare i with
        compare = []
        for j in iter_indices:
            # Don't want to do (i, j) and (j, i), nor (i, i)
            # Can do this because iter_indices is sorted
            if j > i:
                compare.append(j)

        for j in compare:
            similarities[i].update({j: iterations[i].intersection(iterations[j])})

    return similarities


if __name__ == "__main__":
    ITERATIONS_PATH = join("E:/", "Lucas", "Dropbox", "school", "ellipcurves_matd92", "diffstat_iterations")
    # ITERATIONS_PATH = join("/home", "lucas", "Documents", "diffstat_iterations")
    ITERATION_DIRS = ["{}".format(i) for i in range(1, 11)]
    # The constant amount of curves in each curve file in each iteration
    NUM_CURVES_PER_RUN = 1000
    # The number of different curve files for each iteration, excluding for values 1--10
    NUM_RUNS = 5
    NUM_CURVES_PER_ITER = NUM_CURVES_PER_RUN * NUM_RUNS

    iterations = collect_curves(ITERATIONS_PATH, ITERATION_DIRS)

    iter_indices = list(iterations.keys())
    iter_indices.sort()

    percent_per_iter = lambda n: 100 * (n / NUM_CURVES_PER_ITER)

    print("\ncurves unique in each iteration:")
    for i in iter_indices:
        num_unique = len(iterations[i])
        print("iteration {}: {} / {}% unique".format(i, num_unique, percent_per_iter(num_unique)))

    min_unique = min([len(s) for s in iterations.values()])
    print("\nminimum unique curves in each iteration: {} / {}%".format(min_unique, percent_per_iter(min_unique)))

    # Pairwise comparison between iterations
    similarities = pairwise_duplicates(iterations, iter_indices)

    print("\npairwise duplicates between iterations:")
    num_shared_ = []
    for i in similarities:
        print("iteration {}:".format(i))
        for j in similarities[i]:
            num_shared = len(similarities[i][j])
            print("\tshared with iteration {}: {} / {}%".format(j, num_shared, percent_per_iter(num_shared)))
            num_shared_.append(num_shared)

    max_shared = max(num_shared_)
    print("\nmax pairwise duplicate between iterations: {} (of {}) / {}%\n".format(max_shared, NUM_CURVES_PER_ITER,
                                                                                   percent_per_iter(max_shared)))
