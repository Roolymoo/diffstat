from math import factorial, floor, log

from random import randint, choice

from sympy import Number, sqrt_mod, totient, legendre_symbol, isprime

import matplotlib.pyplot as plt


def _compute_prod(p, n, comp):
    """(Point, int, dict) -> Point
    If n is in comp, then return comp[n] == (2^n) * p. Otherwise computes and returns 2^n * p.

    Assumes n >= 0.
    """
    if n in comp:
        return comp[n]
    else:
        # n - 1 cannot be < 0, because n would have to be 0 at some point and the above code would run instead
        prev_prod = _compute_prod(p, n - 1, comp)
        return prev_prod + prev_prod


# TODO: change name to _compute_slope_data?
def _compute_slope(p1, p2):
    """(Point, Point) -> tuple
    (slope of p1 and p2, True) is returned if such slope is defined, otherwise (factoring data, False) is returned.

    Assumes p1, p2 share a common Curve. If the denominator in computing the slope has a multiplicative inverse mod the
    common modulus N, then (slope, True) is returned. If the denominator equals 0, then (0, False) is returned--this
    indicates a point at infinity. Otherwise (gcd(denominator, N), False) is returned. This gcd is greater than 1 and
    so can be used for factoring N.
    """
    gcd = Number.gcd

    # Assumed that p1, p2 share common Curve
    curve = p1.curve
    mod = curve.mod
    b = curve.b

    # Distinct point case
    if p1 != p2:
        numerator = (p2.y - p1.y) % mod
        denominator = (p2.x - p1.x) % mod
    # Equal point case
    else:
        numerator = (3 * pow(p1.x, 2) + b) % mod
        denominator = (2 * p1.y) % mod

    # Inverse of denominator
    if denominator == 0:
        # Special case of no multiplicative inverse
        # TODO: change to 0, True--really the slope is defined, it is an infinite slope
        return 0, False
    # Since denominator is not 0, computing gcd is meaningful now
    div = gcd(denominator, mod)
    if div != 1:
        # No multiplicative inverse mod mod
        return div, False

    # Otherwise can compute the multiplicative inverse of the denominator (using Euler's theorem)
    inv_denominator = pow(denominator, totient(mod) - 1, mod)

    slope = (numerator * inv_denominator) % mod
    return slope, True


def _is_congr(x, y, mod):
    """(int, int, int) -> bool
    Returns true iff x is congruent to y mod mod.

    Assumes mod is a positive integer.
    """
    return (x % mod) == (y % mod)


def _mul(p, n):
    """(Point, int) -> Point
    Returns p scaled by n, i.e. n * p.

    Assumes n >= 0.
    """
    # Points at infinity are the additive identity--see Point
    if not p.is_finite():
        return p

    # TODO: find out if 0 * infinity = 0
    if n == 0:
        return Infinity()
    elif n == 1:
        return p

    # Keep track of largest computation (largest n) 2^n * p to reuse for (possible) next larger computation
    # comp == {power n: 2^n * p}
    comp = {0: p}

    # Compute product as sum of form 2^a p + 2^b p + ... + 2^z p in binary (i.e. n in binary)
    # Start with point at infinite, it being the additive identity like 0 is usually used as
    sum = Infinity()
    remainder = n % 2
    quotient = int(n / 2)
    power = 0
    # Loop guaranteed to execute at least once because n >= 2
    while quotient != 0:
        # If remainder == 0, don't add anything, only increases power with another 2
        if remainder != 0:
            # Add (2^power) * p
            q = _compute_prod(p, power, comp)
            comp = {power: q}
            sum += q

        remainder = quotient % 2
        quotient = int(quotient / 2)
        power += 1

    # Collect the 2's
    return sum + _compute_prod(p, power, comp)


def _add(p1, p2, slope_data):
    """(Point, Point, tuple) -> Point|NoneType
    Returns the Point that is the sum of p1 and p2 if it is defined, otherwise None is returned.

    Assumes p1, p2 share a common Curve. Assumes slope_data == _compute_slope(p1, p2). If slope_data == (int, True),
    i.e. the slope is defined, then the Point that is their sum is returned. If (0, False), then an instance of
    Infinity is returned. Otherwise the sum is not defined in which case None is returned.
    """
    slope, is_valid = slope_data
    if not is_valid:
        if slope == 0:
            # Denominator 0
            return Infinity()
        else:
            # Slope is not defined
            return None

    # Otherwise slope is finite
    # Assumes that p1 and p2 share common Curve
    mod = p1.curve.mod

    sum_x = (pow(slope, 2) - p1.x - p2.x) % mod
    sum_y = (slope * (p1.x - sum_x) - p1.y) % mod

    return Point(sum_x, sum_y, p1.curve)


def _get_rand_pt(e):
    """(Curve) -> Point|NoneType
    Attempts to find a random point on the curve. If successful, returns it, otherwise None is returned.

    Priority is put on finite points such that only an infinite point will be returned if that is the only point on the
    curve. If the curve has NOT had its points generated, a point will be found without generating all points.
    Furthermore no points found to lie on the curve are added to the curves points.
    """
    pts_length = len(e.pts)
    if pts_length > 0:
        # The curve already has had points generated
        if pts_length > 1:
            # Priority is put on finite points
            # The point at infinity always is at e.pts[0]
            return choice(e.pts[1:])
        else:
            # Only point is the one at infinity
            return e.pts[0]

    # The curve has not had its points generated: find a point without much work
    x_coord = randint(0, e.mod - 1)
    limit = 1000  # This limit does fail for certain curves, but the probability is unknown
    i = 0
    while i < limit:
        # TODO: the following overlaps with Curve.gen_pts()
        # See if x_coord + i makes any point on the curve
        f = lambda x: (pow(x, 3) + e.b * x + e.c) % e.mod  # The curve
        squares = sqrt_mod(f(x_coord + i), e.mod, True)  # List of all y such that y^2 % e.mod == f(x_coord + i)
        if not squares:
            # Didn't work, move on and try others
            i += 1
        else:
            # Found some points, arbitrarily take the first one
            return Point(x_coord, squares[0], e)

    # Passed limit and found nothing, give up
    return None


def _factor(n):
    """(int) -> tuple|NoneType
    Attempts to factor n into two factors. If successful, returns them as a tuple, otherwise None is returned.

    Assumes n is a positive integer."""
    coeff_lower_bdd, coeff_upper_bdd = 0, 100
    num_curves = 5
    fact = 5

    # Create specified amount of random Curve's mod n
    curves = []
    for i in range(num_curves):
        b = randint(coeff_lower_bdd, coeff_upper_bdd)
        c = randint(coeff_lower_bdd, coeff_upper_bdd)
        curves.append(Curve(b, c, n))

    # Cycle through curves trying to find a factor of n
    for e in curves:
        # Pick a random (finite) point
        p = _get_rand_pt(e)
        if p is None:
            # Failed finding a random point, try next curve
            continue

        # Start scaling p to find a factor of n
        sum = p
        for i in range(2, factorial(fact)):
            slope, is_valid = _compute_slope(sum, p)
            # TODO: some redundent code here with _add
            if not is_valid:
                if slope == 0:
                    # Then sum == infinite, hence this curve attempt failed, try next curve
                    break
                else:
                    # Found a factor
                    return slope, int(n / slope)
            else:
                # Sum is finite, continue
                sum = _add(sum, p, (slope, is_valid))

    # Failed trying to find a factor
    return None


def _num_pts_leg(curve, start=0):
    """(Curve) -> int
    Returns the number of points on curve, including infinite, using Legendre symbols.

    An optional argument start specifies to only count the points with x-coordinate at least start. See Theorem 4.14,
    p. 105 in Washington for the theorem statement and proof.
    """
    b = curve.b
    c = curve.c
    mod = curve.mod

    # The right-hand-side of curve: y^2 = x^3 + bx + c
    y_sqrd = lambda x: (pow(x, 3) + b * x + c) % mod

    x = start
    sum = 0
    while x < mod:
        _y_sqrd = y_sqrd(x)
        sum += legendre_symbol(_y_sqrd, mod)
        x += 1

    # The original algorithm adds mod, but that would be too much because we use [start, mod) which might not be
    # [0, mod) so (mod - start) is used instead
    return sum + (mod - start) + 1


# Old num_pts method trying to be clever, didnt work bc assumed finite field with 2^q elements could be represented as
# Z_p for some p...
#
#
# length = len(self.pts)
# if length != 0:
#     return length
#
# # Use other elementary methods
# p = self.mod
#
# # Use base method to calculate points up to a certain point
# base, power = _log_approx(p)
# approx_num_pts = _num_pts_base(self, base, power)
#
# # Use Legendre method to calculate the rest
# # Have calculated up to pow(base, power), exclusively
# return approx_num_pts + _num_pts_leg(self, start=pow(base, power))
def _num_pts_base(curve, base, power):
    """(Curve, int, int) -> int
    Returns the number of points on curve using the base method.

    Assumes curve.mod is a power of an nonnegative integer. See p. 102--104 in Washington for the algorithm.
    """
    # Get the number of points for the base
    base_curve = Curve(curve.b, curve.c, base)
    base_num_pts = _num_pts_leg(base_curve)

    # Start recurrence relation
    a = base + 1 - base_num_pts

    s_before = 2
    s_current = a
    i = 2
    while i <= power:
        s_current = a * s_current - base * s_before
        i += 1

    return pow(base, power) + 1 - s_current


# TODO: allow to go > p, because not doing this removes possible better n
def _log_approx(p):
    """(int) -> tuple
    Returns n in {2, 3, 5, 7} such that p - pow(n, floor(log_n (p))) is minimal WRT these choices of n.

    Also returns Assumes p is prime. The idea is that pow(n, floor(log_n (p))) best approximates p while still be < p.
    """
    primes = {3, 5, 7}
    int_log = lambda p, n: floor(log(p, n))
    diff = lambda p, n, _int_log: p - pow(n, _int_log)

    _int_log = int_log(p, 2)
    best_approx = diff(p, 2, _int_log), 2, _int_log
    for n in primes:
        _int_log = int_log(p, n)
        new_diff = diff(p, 2, _int_log)
        if new_diff > best_approx[0]:
            best_approx = new_diff, n, _int_log

    return best_approx[1:]


def _gen_x_coords(pts):
    """(list) -> list
    Returns a list of the unique x-coordinates from the points in pts.

    Assumes pts is a list of Point's. The returned list preserves the ordering of pts."""
    x_coords = []
    # [1:] to pass Infinity
    for p in pts[1:]:
        # Avoid duplicates--generally, each x-coordinate has two y-coordinates
        if p.x not in x_coords:
            x_coords.append(p.x)

    return x_coords


class Curve:
    """Equation of the form y^2 = x^3 + bx + c (mod n) where n is a positive integer, known as the Weierstrass equation
     for elliptic curves, over the integers mod n. The curve is considered in the standard 2-dimensional Euclidean
     space. Points can be tested whether they satisfy the equation. All points that satisfy the equation can be
     generated, and the equation can be graphed.

     Public instance variables:
     b -- x-term coefficient
     c -- constant coefficient
     mod -- modulus
     pts -- instances of Point that satisfy the equation

     Public methods:
     gen_pts() -- calculates all points that satisfy the equation and assigns them as a list to pts. The point at
                  infinity is included at index 0
     graph() -- plots all points that satisfy the equation using matplotlib and displays the (interactive) plot
     """

    def __init__(self, b, c, mod):
        """(Curve, int, int, int) -> NoneType
        x-term coefficient b, constant term coefficient c, curve modulus mod.
        """
        self.b = b
        self.c = c
        self.mod = mod
        self.pts = []  # List of points on the curve

    def __repr__(self):
        """(Curve) -> str
        Returns str of form, E : y^2 = x^3 + bx + c mod n.

        The bx term or c term is eliminated if b or c is 0, respectively. Furthermore  if b==1, then the 1 is
        eliminated (i.e. x instead of 1x).
        """
        s = "E : y^2 = x^3 + {}x + {} mod {}".format(self.b, self.c, self.mod)
        if self.b == 0:
            s = s.replace(" + {}x".format(self.b), "")
        if self.c == 0:
            s = s.replace(" + {}".format(self.c), "")
        if self.b == 1:
            s = s.replace(" + {}x".format(self.b), " + x")

        return s

    def __str__(self):
        """(Curve) -> str"""
        return self.__repr__()

    def __eq__(self, e2):
        """(Curve, Curve) -> bool
        Equality iff coefficients and the moduli match.
        """
        return self.b == e2.b and self.c == e2.c and self.mod == e2.mod

    def __contains__(self, p):
        """(Curve, Point) -> bool
        Returns True iff p satisfies the equation defined by self, regardless of whether it is stored in self.pts.
        """
        # Always contains a point at infinity
        if not p.is_finite():
            return True
        else:
            x, y = p.x, p.y
            b, c, mod = self.b, self.c, self.mod
            return _is_congr(pow(y, 2), pow(x, 2) + b * x + c, mod)

    def gen_pts(self):
        """(Curve) -> NoneType
        Computes all points (x,y) satisfying the equation defined by self and appends them to self.pts as Points.

        First point appended is the one at infinity.
        """
        mod = self.mod
        f = lambda x: (pow(x, 3) + self.b * x + self.c) % mod  # y^2 = f(x), equation defined by self

        self.pts.append(Infinity())

        # The x-coordinates of possible points are limited to 0, ..., mod - 1
        # Only such x-coordinates x such that the congruence equation y^2 = f(x) has solutions y are points
        x = 0
        while x != mod:
            squares = sqrt_mod(f(x), mod, True)  # List of all y such that y^2 % mod == f(x)
            for sq in squares:
                self.pts.append(Point(x, sq, self))

            x += 1

    def graph(self):
        """(Curve) -> NoneType
        Plots self.pts, except the point at infinity, using matplotlib and displays the resulting plot.

        If self.pts is None, then self.gen_pts() is called. Otherwise it is assumed that all desired points are in
        self.pts.
        """
        if not self.pts:
            self.gen_pts()

        # Only work with self.pts[1:] because by construction the point at infinity is always at the front
        x_coords = [p.x for p in self.pts[1:]]
        y_coords = [p.y for p in self.pts[1:]]
        x_axis_padding = 1
        y_axis_padding = 1
        x_range = [-x_axis_padding, self.mod]
        y_range = [-y_axis_padding, max(y_coords) + y_axis_padding]
        x_label = "x-axis"
        y_label = "y-axis"
        title = self.__str__()

        # Start graphing
        plt.title(title)
        # Axes labels
        plt.xlabel(x_label)
        plt.ylabel(y_label)

        # Axes
        # x-axis
        plt.plot(x_range, [0, 0], "k-")
        # y-axis
        plt.plot([0, 0], y_range, "k-")
        # Range for axes
        plt.axis(x_range + y_range)

        # Grid
        plt.grid()

        # Points
        plt.plot(x_coords, y_coords, "ko")

        plt.show()

    def min_distance(self):
        """(Curve) -> int
        Returns the minimum distance between the x-coordinates of the points on self.

        If the points on self have not been generated, gen_pts() will be called and afterward self.pts will be
        cleared.
        """
        # Check if pts have already been generated, and temporarily generate them if they haven't
        pts_already_gen = True
        if len(self.pts) == 0:
            pts_already_gen = False
            self.gen_pts()

        # Set comprehension first to lose duplicates, then list so can index later, [1:] to pass Infinity
        x_coords = list({p.x for p in self.pts[1:]})

        # By construction x_coords is already sort in ascending order. This means that the minimum difference can be
        # calculated pairwise for len(x_cords) - 1 total comparisons
        if len(self.pts) == 1:
            # No differences can be calculated
            return 0

        min_diff = abs(x_coords[0] - x_coords[1])
        if min_diff == 1:
            # Can't get smaller than 1
            pass
        else:
            # Try to find a smaller difference
            i = 1
            while i <= len(x_coords) - 2:
                new_diff = abs(x_coords[i] - x_coords[i + 1])
                min_diff = min(new_diff, min_diff)
                if min_diff == 1:
                    break
                i += 1

        # Clear self.pts if they were not generated before this function
        if not pts_already_gen:
            self.pts = []

        return min_diff

    def max_distance(self):
        """(Curve) -> int
        Returns the maximum distance from the x-coordinate of a point to the that of the point to the right of it.

        Specifically, if P = {finite points on self}, then max{|P.x - (right point of P).x|} is returned. For example,
        if self.pts == [inf, (0, 5), (1, 2), (4, 7)], then 3 is returned. If the points on self have not been generated,
        then all x-coordinates are generated via Legendre symbols if self.mod is an odd prime (faster), or
        self.gen_pts() is called (slower).
        """
        # TODO: Maybe _gen_x_coords -> _gen_x_coords_from_pts, and the below code for generating the x-coordinates can
        # TODO: become _gen_x_coords.
        # TODO: Make more efficient by making a helper function version that assumes that p is an odd prime and the
        # TODO: points have not already been generated.
        # TODO: Repeating code with min_distance.
        # Generate x-coordinates
        x_coords = []
        delete_pts = False
        if len(self.pts) != 0:
            # pts already generated, use those
            x_coords = _gen_x_coords(self.pts)
        elif isprime(self.mod) and (self.mod != 2):
            # Generate x-coordinates using Legendre symbols
            # The right-hand-side of curve: y^2 = x^3 + bx + c
            y_sqrd = lambda x: (pow(x, 3) + self.b * x + self.c) % self.mod
            x = 0
            while x < self.mod:
                _y_sqrd = y_sqrd(x)
                if legendre_symbol(_y_sqrd, self.mod) != -1:
                    x_coords.append(x)
                x += 1
        else:
            # Generate all points
            self.gen_pts()
            delete_pts = True
            x_coords = _gen_x_coords(self.pts)

        # By construction x_coords is already sorted in ascending order
        if len(x_coords) <= 1:
            # No differences can be calculated
            return 0

        max_diff = abs(x_coords[0] - x_coords[1])
        # Try to find a larger difference
        i = 1
        while i <= len(x_coords) - 2:
            new_diff = abs(x_coords[i] - x_coords[i + 1])
            max_diff = max(new_diff, max_diff)
            i += 1

        # Clear self.pts if they were generated for temporary use
        if delete_pts:
            self.pts = []

        return max_diff

    def num_pts(self):
        """(Curve) -> int
        Returns the number of points on self, counting infinite.

        Uses len() if self.pts != [], else Legendre symbols are used if self is mod a prime, or if the mod is composite
        the points are temporarily generated and their length is taken, and the points are then forgotten.
        """
        length = len(self.pts)
        if length != 0:
            return length

        # Use other elementary methods
        if isprime(self.mod):
            return _num_pts_leg(self)
        else:
            self.gen_pts()
            length = len(self.pts)
            self.pts = []
            return length


class Point:
    """A point (x,y) satisfying the equation defined by an instance of Curve. Addition p1 + p2 and multiplication by a
    positive integer n * p == p * n are defined. Addition satisfies abelian group axioms where points at infinity
    (instances of Infinity) are the identity. Hence any finite point P added with one at infinity is P again, and sums
    of points at infinity is the point at infinity again.

    Public instance variables:
    x -- first coordinate
    y -- second and final coordinate
    curve -- Curve instance (x,y) satisfies defined equation of

    Public methods:
    is_finite() -- returns True iff self is not the subclass Infinity
    """

    def __init__(self, x, y, e):
        """(Point, int, int, Curve) -> NoneType
        Interprets x, y mod e.mod.
        """
        mod = e.mod

        self.x = x % mod
        self.y = y % mod
        self.curve = e

    def __repr__(self):
        """(Point) -> str
        Returns string '(self.x, self.y)'.
        """
        return "({}, {})".format(self.x, self.y)

    def __str__(self):
        """(Point) -> str"""
        return self.__repr__()

    def __eq__(self, p2):
        """(Point, Point) -> bool
        Returns True iff self and p2 are equal points on the same curve.

        If p2 is an instance of Infinity, or if p2 lies on a different Curve, then False is returned. Otherwise self
        and p2 lie on the same Curve, in which case return True iff the coordinates are equal (recall that the
        coordinates are already modulated, and that since the Curve's are the same they therefore share the same
        modulus).
        """
        if not p2.is_finite() or (self.curve != p2.curve):
            return False
        else:
            return self.x == p2.x and self.y == p2.y

    def __add__(self, p2):
        """(Point, Point) -> Point|NoneType
        Returns the Point that is the sum of self and p2 if it is defined, otherwise None.

        If self and p2 lie on different Curve's, then None is returned. Otherwise the sum is attempted to be computed.
        If the sum is not defined, then None is returned. Otherwise sum is returned.
        """
        # Infinite case
        if not p2.is_finite():
            return p2 + self

        # Finite case, so check curves match first
        if self.curve != p2.curve:
            return None

        # Curves match, try to compute sum
        slope_data = _compute_slope(self, p2)
        return _add(self, p2, slope_data)

    def __mul__(self, n):
        """(Point, int) -> Point
        Return self scaled by n.

        Assumes N >= 0.
        """
        return _mul(self, n)

    def __rmul__(self, n):
        """(Point, int) -> Point"""
        return self.__mul__(n)

    def is_finite(self):
        """(Point) -> bool
        Returns true iff self is not infinite.
        """
        return not isinstance(self, Infinity)


# Maybe treat Infinity as a number instead of a point, and consider infinity = (infinity, infinity) like in the book
# Also it may be confusing that Infinity is a subclass of Point, e.g. needing to add in (finite) when talking about
# Points often.
class Infinity(Point):
    """The point at infinity added to the standard 2-dimensional Euclidean space, rigorously defined via 2-dimensional,
    real projective geometry. Here it is uniquely characterized as a Point that lies on all Curve's and is the additive
    identity for all Point's on a Curve--see Point.

    Public instance variables:
    None.

    Public methods:
    None.
    """

    def __init__(self):
        """(Infinity) -> NoneType"""
        # No instance variables
        ...

    def __repr__(self):
        """(Infinity) -> str
        Returns the string 'inf'."""
        return "inf"

    def __eq__(self, p2):
        """(Infinity, Point) -> bool
        Returns True iff p2 is Infinity."""
        return not p2.is_finite()

    def __add__(self, p2):
        """(Infinity, Point) -> Point
        Returns p2 if p2 is a (finite) Point, else returns self."""
        if p2.is_finite():
            return p2
        else:
            return self


def read_curves(file_name):
    """(str) -> list
    Returns a list of Curve's generated from the contents file_name.

    Assumes file_name corresponds to an existing file, and that it includes its extension. If file is empty (of Curve
    data), then an empty list is returned. Assumes proper formatting:

    <b_value> <c_value> <mod_value>
    <b_value> <c_value> <mod_value>
    ...
    \n

    """
    curve_list = []
    with open(file_name, "r") as file:
        for curve_data in file:
            if curve_data == "\n":
                # EOF
                break

            b, c, mod = [int(data) for data in curve_data.split()]
            e = Curve(b, c, mod)
            curve_list.append(e)

    return curve_list


# TODO: fix issue with newlines (seems to have to do with append behaviour, maybe append adds a newline when there is
# TODO: an existing newline, because it doesn't seem to add a newline if there isn't one already there)
def write_curves(curve_list, file_name):
    """(list, str) -> NoneType
    Writes the Curve's in curve_list to file_name.

    The formatting and data written is specified in read_curves(). If file_name exists, then the Curve data is appended
    to the existing data, where it is assumed that file_name contains Curve data only, as per the proper formatting.
    Otherwise file_name is created. No newlines are put before the data, but after each curve data written a newline
    is inserted (thus the empty line in an existing file is preserved)."""
    with open(file_name, "a") as file:
        for e in curve_list:
            file.write("{} {} {}\n".format(e.b, e.c, e.mod))

        # EOF indicator line
        file.write("\n")


def factor(n):
    """(int) -> dict
    Returns the prime factorization of n.

    Assumes n is a positive integer. The prime factorization is a dict with keys as prime factors and their values as
    their multiplicities. If no prime factors can be found, the empty dict is returned.
    """
    old_factorization = {i: 1 for i in _factor(n)}
    while True:
        new_factorization = old_factorization.copy()
        for k in new_factorization:
            more_factors = _factor(k)
            if more_factors is None:
                # Failed or no more factors
                continue
            else:
                # Remove old factor
                multiplicity = new_factorization.pop(k)
                for l in more_factors:
                    if l in new_factorization:
                        new_factorization[l] += multiplicity
                    else:
                        new_factorization.update({l: multiplicity})

        if new_factorization == old_factorization:
            # Failed to find more factors
            return new_factorization
        else:
            # Found some new factors, try to find more
            old_factorization = new_factorization.copy()


def min_mult(p):
    """(Point|Infinity) -> int
    Returns the smallest nonnegative integer n such that n * p == infinity."""
    if not p.is_finite():
        return 0

    # incrementally check until infinity is found
    n = 1
    while (n * p).is_finite():
        n += 1

    return n
