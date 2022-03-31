from diffstat import _is_singular

from ellipcurve import Curve


"""
Generates all elliptic curves mod a given prime, prints the number of curve pairs that have equal j-invariants and
are equal (as sets), and these latter curves themselves (as tuples (b, c, p)).
"""


def _gen_curves(p):
    """(int) -> generator
    Returns a generator of all elliptic curves (b, c, p) mod p, for p prime.

    Assumes p is prime. Hence the returned list contains curves with b-coefficient i, c-coefficient j, for all
    i, j mod p, assuming that such a curve exists.
    """
    for i in range(p):
        for j in range(p):
            e = i, j, p
            if not _is_singular(*e):
                yield e


def j_invariant(b, c, p):
    """(int, int, int) -> int
    Returns the j-invariant of the curve given by (b, c, p).

    Assumes p is prime and that (b, c, p) denotes an elliptic curve (this guarantees that division by zero does not
    occur).
    """
    return (1728 * ((4 * pow(b, 3)) / (4 * pow(b, 3) + 27 * pow(c, 2)))) % p


def gen_eq_jvar(p):
    """(int) -> generator
    Returns a generator of all curve pairs (e1, e2) mod p that have equal j-invariants.

    Assumes p is prime. Specifically the curves are taken from _gen_curves(p).
    """
    # Skip curves so it is not the case that both (e1, e2) and (e2, e1) are yielded
    skip = 0
    for e1 in _gen_curves(p):
        i = skip
        for e2 in _gen_curves(p):
            if i > 0:
                # Skip, already did (e2, e1) (otherwise would have returned (e1, e2)
                i -= 1
                continue
            elif e1 == e2:
                # Don't want to consider pairs (e, e)--trivially they have equal j-invariant
                continue
            else:
                if j_invariant(*e1) == j_invariant(*e2):
                    yield (e1, e2)

        skip += 1


def gen_eq_curves(curve_pairs):
    """(generator) -> generator
    Returns a generator of all curve pairs (e1, e2) in curve_pairs that are equal (set-wise).

    For curves with large primes this will take awhile.
    """
    for e1, e2 in curve_pairs:
        curve1 = Curve(*e1)
        curve2 = Curve(*e2)
        # For large prime this will take awhile
        curve1.gen_pts()
        curve2.gen_pts()

        if curve1.pts == curve2.pts:
            yield (e1, e2)


if __name__ == "__main__":
    PRIME = 37

    curves_eq_jvar = gen_eq_jvar(PRIME)

    # Find the curves that have equal j-invariant AND are equal to each other and print them
    # Also find number of these curves, and print those
    num_eq_and_jvar = 0
    # Padding
    print()
    for e1, e2 in gen_eq_curves(curves_eq_jvar):
        print(e1, e2)
        num_eq_and_jvar += 1

    print("\nnum equal curves e1, e2 mod {} with equal j-invariant: {}".format(PRIME, num_eq_and_jvar))
