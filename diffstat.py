from time import time
from datetime import datetime
from random import seed, randint
import argparse

from sympy import sieve, legendre_symbol, randprime


def _write_curve(file, e):
    """(TextIOWrapper, tuple) -> NoneType
    Writes the curve e onto file.

    Note that there is no newline at the end. This is left to _gen_x_coords(*e) when the last x--coordinate of e is
    written in by _write_x_coord().
    """
    file.write("{} {} {}".format(*e))


def _write_x_coord(file, x):
    """(TextIOWrapper, int) -> NoneType
    Writes x onto file.

    Note that there is no newline at the end. This is left to _gen_x_coords(*e) when the last x--coordinate of e is
    written.
    """
    file.write(" {}".format(x))


def _gen_x_coords(b, c, p, file=None):
    """(int, int, int, TextIOWrapper|NoneType) -> generator
    Returns a generator of the unique x-coordinates of the points on the curve given by (b, c, p).

    Assumes 0 <= b, c < p, and that p is an odd prime. The generator is sorted in ascending order. The generated
    x-coordinates are unique in that if two points on the curve share an x-coordinate, then that x-coodinate is
    generated only once. If there are no x-coordinates, then an empty generator is returned--p > 3 guarantees
    a nonempty generator (follows from Hasse's Theorem).

    If file is given, then the generator also writes the found x--coordinates to it. In this case it is assumed that
    _write_curve(file, (b, c, p)) has already been called. Also a newline will be written after the last x--coordinate
    has been written.

    See gen_rand_curves() for more information on the topic of curves and points.
    """
    # Use Legendre symbols (can do this since p is an odd prime)
    # The right-hand-side of curve given by (b, c, p): y^2 = x^3 + bx + c mod p
    y_sqrd = lambda x: (pow(x, 3) + b * x + c) % p
    x = 0
    while x < p:
        _y_sqrd = y_sqrd(x)
        if legendre_symbol(_y_sqrd, p) != -1:
            yield x
            if file:
                _write_x_coord(file, x)
        x += 1

    if file:
        file.write("\n")


def _avg(iter):
    """(iterable) -> float
    Returns the average of the elements in iterable.

    Assumes iter contains numbers, and that len(iter) > 0.
    """
    return sum(iter) / len(iter)


def _is_singular(b, c, p):
    """(int, int, int) -> bool
    Returns True iff 4b^3 + 27c^2 == 0 mod p.
    """
    return ((4 * pow(b, 3) + 27 * pow(c, 2)) % p) == 0


def _rand_prime(lower_bd, upper_bd):
    """(int, int) -> int
    Returns a random prime p such that lower_bd <= p < 2 * upper_bd.

    Assumes lower_bd > 1, and that lower_bd <= upper_bd.
    """
    x = randint(lower_bd, upper_bd)
    # This is guaranteed to return a prime because the range is of the form (y, 2 * y) for y > 1--see Bertrand's
    # postulate.
    p = randprime(x, 2 * x)
    # Since lower_bd <= x <= upper_bd, and since x <= p < 2 * x, therefore lower_bd <= p < 2 * upper_bd
    return p


# TODO: Look for references to this function generating curves with p satisfying prime_min <= p <= prime_max, because
# TODO: now it is ... p < prime_max.
def gen_rand_curves(coeff_min, coeff_max, prime_min, prime_max, num_curves):
    """(int, int, int, int, int) -> generator
    Returns a generator of random tuples (b, c, p) defining elliptic curves within these parameters.

    Specifically the curves are cubics E: y^2 = x^3 + bx + c mod p, where x, y, b, c are in Z_p, the integers mod p,
    where p > 3 (so as to avoid fields characteristic 2 and 3 (i.e. Z_2 and Z_3), and empty curves (follows from
    Hasse's Theorem)), and 4b^3 + 27c^2 != 0 mod p (so that E is non-singular). It is assumed that prime_min > 3. It is
    not required that prime_min and prime_max be primes. However prime_max must be divisible by 2. Each tuple (b, c, p)
    satisfies, in addition to the above:

    - coeff_min <= b, c <= coeff_max
    - prime_min <= p < prime_max
    - num_curves of these tuples are generated
    """
    seed()

    # This needs to be an integer for _rand_prime() to work
    randprime_bound = prime_max / 2

    i = 0
    while i < num_curves:
        b = randint(coeff_min, coeff_max)
        c = randint(coeff_min, coeff_max)
        p = _rand_prime(prime_min, randprime_bound)

        # Ensure curve is non-singular
        # TODO: It appears that only curves of the form E: y^2 = x^3 - 3x +/- 2 mod p for any prime p are singular.
        # TODO: If this is the case, then probably can change random number generating to avoid these and then remove
        # TODO: the singular checks, probably making this much faster.
        if not _is_singular(b, c, p):
            yield b, c, p
            i += 1
        else:
            # Try again. (It appears the probability of this happening is really low so this shouldn't be a
            # computational issue in the way of often getting stuck finding a non-singular curve.)
            continue


def max_distance(x_coords):
    """(generator) -> int
    Returns the maximum distance between the x-coordinates in x_coords.

    Assumes x_coords is the return value of _gen_x_coords(b, c, p) for some b, c, p satisfying that 0 <= b, c < p, and
    that p > 3 is prime. In letting E be the curve given by (b, c, p), if P = (finite points on E), then
    max{|P[i].x - P[i + 1].x|} is returned. For example, if P == ((0, 5), (1, 2), (4, 7)), then 3 is returned. The
    requirement that p > 3 ensures that there is at least one finite point on E.

    Performance note: Iterating through x_coords may have the generator write its contents to a file.

    See gen_rand_curves() for more information on the topic of curves and points.
    """
    max_diff = 0
    # Since p > 3, therefore x_coords has at least one point
    x_left = next(x_coords)
    for x_right in x_coords:
        new_diff = abs(x_left - x_right)
        max_diff = max(new_diff, max_diff)
        x_left = x_right

    return max_diff


def arg_parser():
    """(NoneType) -> argparse.ArgumentParser"""
    parser = argparse.ArgumentParser("Point density statistics for elliptic curves.")

    parser.add_argument("coeff_min", type=int,
                        help="Minimum value for coefficients defining generated elliptic curves.")
    parser.add_argument("coeff_max", type=int,
                        help="Maximum value for coefficients defining generated elliptic curves.")
    parser.add_argument("prime_min", type=int,
                        help="Minimum value for prime moduli defining generated elliptic curves.")
    parser.add_argument("prime_max", type=int,
                        help="Maximum value for prime moduli defining generated elliptic curves.")
    parser.add_argument("num_curves", type=int, help="Number of curves to be generated.")

    return parser


# Deprecated
def max_index(iter):
    """(iterable) -> tuple
    Returns (max(iter), index) where iter[index]==max(iter).

    Assumes there is an ordering on iter, and that len(iter)>=1.
    """
    length = len(iter)
    if length == 1:
        return iter[0], 0

    i = 0
    max_elem = iter[0], 0
    while i < length:
        elem = iter[i]
        if elem > max_elem[0]:
            max_elem = elem, i
        i += 1

    return max_elem


if __name__ == "__main__":
    # Collect arguments
    parser = arg_parser()

    ARGS = parser.parse_args()

    COEFF_MIN = ARGS.coeff_min
    COEFF_MAX = ARGS.coeff_max
    PRIME_MIN = ARGS.prime_min
    PRIME_MAX = ARGS.prime_max
    NUM_CURVES = ARGS.num_curves

    # File to write results and statistics in
    TIME_ID = str(datetime.now()).replace(":", "-").split()
    RESULTS_FILE_NAME = "results_{}_{}_{}_{}_{}_{}_{}.txt".format(COEFF_MIN, COEFF_MAX, PRIME_MIN, PRIME_MAX,
                                                                  NUM_CURVES, *TIME_ID)
    # File to write generated curves in
    CURVES_FILE_NAME = "curves_{}_{}_{}_{}_{}_{}_{}.txt".format(COEFF_MIN, COEFF_MAX, PRIME_MIN, PRIME_MAX, NUM_CURVES,
                                                                *TIME_ID)

    # Curves are written progressively to file as they are generated to avoid a large memory footprint
    curves_file = open(CURVES_FILE_NAME, "w")

    # Generate random curves with these specifications
    # Count time elapsed for upcoming calculations
    START_TIME = time()

    # Calculate maximum distances for all curves
    # Keep track of these distances
    max_distance_list = []
    curves = gen_rand_curves(COEFF_MIN, COEFF_MAX, PRIME_MIN, PRIME_MAX, NUM_CURVES)
    # Empty line to pad upcoming "% done" updates for each curve
    print()
    for e in curves:
        _write_curve(curves_file, e)
        b, c, p = e
        x_coords = _gen_x_coords(b, c, p, curves_file)
        max_distance_list.append(max_distance(x_coords))
        print("{}% done".format((len(max_distance_list) / NUM_CURVES) * 100))

    curves_file.close()

    # Results and statistics
    TOTAL_TIME = time() - START_TIME

    TIME_TAKEN_SEC = TOTAL_TIME
    TIME_TAKEN_MIN = TIME_TAKEN_SEC / 60
    TIME_TAKEN_HOUR = TIME_TAKEN_MIN / 60
    TIME_TAKEN_PER_CURVE_SEC = TOTAL_TIME / NUM_CURVES
    TIME_TAKEN_PER_CURVE_MIN = TIME_TAKEN_PER_CURVE_SEC / 60

    STATS_MSG = """\nTotal time taken: {} hours / {} minutes / {} seconds
                \nAverage time taken per curve: {} minutes / {} seconds
                \nMaximum of the maximum distances: {}
                \nAverage of the maximum distances: {}""".format(TIME_TAKEN_HOUR, TIME_TAKEN_MIN, TIME_TAKEN_SEC,
                                                                 TIME_TAKEN_PER_CURVE_MIN, TIME_TAKEN_PER_CURVE_SEC,
                                                                 max(max_distance_list), _avg(max_distance_list))

    print(STATS_MSG)

    with open(RESULTS_FILE_NAME, "w") as file:
        file.write("{} {} {} {} {}\n".format(COEFF_MIN, COEFF_MAX, PRIME_MIN, PRIME_MAX, NUM_CURVES))
        file.write(STATS_MSG)
