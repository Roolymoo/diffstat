import unittest

from ellipcurve import Curve, Point, Infinity, min_mult, _is_congr, _compute_prod, _compute_slope


# TODO: examples not using primes


class TestComputeProd(unittest.TestCase):

    def setUp(self):
        self.e1 = Curve(4, 4, 5)
        self.p1 = Point(1, 2, self.e1)

    def test_compute_prod_zero(self):
        n = 0
        comp = {0: self.p1}

        product = self.p1

        self.assertEqual(_compute_prod(self.p1, n, comp), product, "Should return P when N == 0")

    def test_compute_prod_one(self):
        n = 1
        comp = {0: self.p1}

        product = Point(2, 0, self.e1)  # == P + P, verified by hand

        self.assertEqual(_compute_prod(self.p1, n, comp), product, "Should return P + P when n == 1")

    def test_compute_prod_small1(self):
        n = 2
        comp = {0: self.p1}

        product = Infinity()  # == 2^2 * P, verified by hand

        self.assertEqual(_compute_prod(self.p1, n, comp), product)

    def test_compute_prod_small2(self):
        n = 2
        comp = {2: Infinity()}

        product = Infinity()  # Verified by hand

        self.assertEqual(_compute_prod(self.p1, n, comp), product, "Should return comp[2] when comp = {2: ...}")


class TestComputeSlope(unittest.TestCase):

    def setUp(self):
        self.e1 = Curve(4, 4, 5)
        self.e2 = Curve(4, 4, 2773)  # Washington & Trappe, p. 353
        self.p1 = Point(1, 2, self.e1)
        self.p2 = Point(4, 3, self.e1)
        self.p3 = Point(1, 3, self.e2)
        self.p4 = Point(1, 3, self.e1)
        self.p5 = Point(1771, 705, self.e2)
        self.inf = Infinity()

    def test_finite_distinct_points_prime_mod(self):
        """Computing slope of distinct points where the slope is finite."""
        ret = 2, True
        self.assertEqual(ret, _compute_slope(self.p1, self.p2))

    def test_finite_same_points_composite_mod(self):
        """Computing slope of equal points where the slope is finite."""
        ret = 2312, True
        self.assertEqual(ret, _compute_slope(self.p3, self.p3))

    def test_zero_denominator(self):
        ret = 0, False
        self.assertEqual(ret, _compute_slope(self.p1, self.p4),
                         "Should return (0, False) when the denominator in computing the slope is congruent to 0")

    def test_no_mult_inverse_composite_mod(self):
        ret = 59, False
        self.assertEqual(ret, _compute_slope(self.p5, self.p3),
                         "Should return (the gcd, False) when denominator is nonzero and has no multiplicative "
                         "inverse")


class TestIsCongr(unittest.TestCase):

    def setUp(self):
        self.mod = 2
        self.n1 = 5
        self.n2 = 21
        self.n3 = -5
        self.n4 = 100
        self.n5 = 0

    def test_basic1(self):
        self.assertTrue(_is_congr(self.n1, self.n2, self.mod))

    def test_basic2(self):
        self.assertFalse(_is_congr(self.n1, self.n4, self.mod))

    def test_negative(self):
        self.assertTrue(_is_congr(self.n1, self.n3, self.mod))

    def test_zero(self):
        self.assertTrue(_is_congr(self.n4, self.n5, self.mod))


class TestMul(unittest.TestCase):

    def setUp(self):
        self.e1 = Curve(4, 4, 5)
        self.p1 = Point(1, 2, self.e1)
        self.inf = Infinity()

    def test_zero_LHS(self):
        n = 0

        product = self.inf

        self.assertEqual(n * self.p1, product, "0 * P == inf")

    def test_zero_RHS(self):
        n = 0

        product = self.inf

        self.assertEqual(self.p1 * n, product, "P * 0 == inf")

    def test_one_LHS(self):
        n = 1

        product = self.p1

        self.assertEqual(n * self.p1, product, "1 * P == P")

    def test_one_RHS(self):
        n = 1

        product = self.p1

        self.assertEqual(self.p1 * n, product, "P * 1 == P")

    def test_two_LHS(self):
        n = 2

        product = self.p1 + self.p1

        self.assertEqual(n * self.p1, product, "2 * P == P + P")

    def test_two_RHS(self):
        n = 2

        product = self.p1 + self.p1

        self.assertEqual(self.p1 * n, product, "P * 2 == P + P")

    def test_small_LHS(self):
        n = 5

        product = self.p1  # Verified by hand

        self.assertEqual(n * self.p1, product)

    def test_small_RHS(self):
        n = 5

        product = self.p1  # Verified by hand

        self.assertEqual(self.p1 * n, product)


class TestCurve(unittest.TestCase):

    def setUp(self):
        self.e1 = Curve(4, 4, 5)
        self.e2 = Curve(4, 4, 7)
        self.e3 = Curve(4, 3, 5)
        self.e4 = Curve(3, 4, 5)
        self.inf = Infinity()

    def test_repr(self):
        s = "E : y^2 = x^3 + 4x + 4 mod 5"
        self.assertEqual(repr(self.e1), s)

    def test_repr_zero_b_coeff(self):
        e = Curve(0, 4, 5)
        s = "E : y^2 = x^3 + 4 mod 5"
        self.assertEqual(repr(e), s)

    def test_repr_zero_c_coeff(self):
        e = Curve(4, 0, 5)
        s = "E : y^2 = x^3 + 4x mod 5"
        self.assertEqual(repr(e), s)

    def test_repr_zero_b_c_coeff(self):
        e = Curve(0, 0, 5)
        s = "E : y^2 = x^3 mod 5"
        self.assertEqual(repr(e), s)

    def test_eq_distinct_b_coeff(self):
        self.assertNotEqual(self.e1, self.e4)

    def test_eq_distinct_c_coeff(self):
        self.assertNotEqual(self.e1, self.e3)

    def test_eq_distinct_prime(self):
        self.assertNotEqual(self.e1, self.e2)

    def test_eq_same(self):
        self.assertEqual(self.e1, self.e1)

    def test_contains_inf(self):
        self.assertIn(self.inf, self.e1, "Curves should contain the point at infinity")

    def test_contains_point_true(self):
        p = Point(1, 2, self.e1)
        self.assertIn(p, self.e1)

    def test_contains_point_false(self):
        p = Point(1, 4, self.e1)
        self.assertNotIn(p, self.e1)

    def test_gen_pts_small(self):
        e = self.e1

        pts = [Infinity(), Point(0, 2, e), Point(0, 3, e), Point(1, 2, e),
               Point(1, 3, e), Point(2, 0, e), Point(4, 2, e), Point(4, 3, e)]

        e.gen_pts()

        self.assertEqual(e.pts, pts)


class TestInfinity(unittest.TestCase):

    def setUp(self):
        self.inf = Infinity()
        self.inf1 = Infinity()
        self.inf2 = Infinity()
        self.e = Curve(4, 4, 5)
        self.p = Point(0, 2, self.e)  # Verified by hand

    def test_subclass(self):
        self.assertIsInstance(self.inf, Point, "Instances of Infinity should be instances of Point")

    def test_repr(self):
        s = "inf"
        self.assertEqual(repr(self.inf1), s)

    def test_eq_inf_diff_instance(self):
        self.assertEqual(self.inf1, self.inf2, "Different instances of Infinity should be equal")

    def test_eq_inf_same_instance(self):
        self.assertEqual(self.inf1, self.inf1, "Same instances of Infinity should be equal")

    def test_eq_pts(self):
        self.assertNotEqual(self.inf1, self.p, "Infinity never equals Points")

    def test_add_same_instance(self):
        self.assertEqual(self.inf, self.inf1 + self.inf1, "Sum of same instance of Infinity should be Infinity")

    def test_add_diff_instance(self):
        self.assertEqual(self.inf, self.inf1 + self.inf2, "Sum of different instances of Infinity should be Infinity")

    def test_add_inf_and_point(self):
        self.assertEqual(self.p, self.inf1 + self.p, "Sum of Infinity and a Point (in this order) should be the Point")

    def test_add_point_and_inf(self):
        self.assertEqual(self.p, self.p + self.inf1, "Sum of a Point and Infinity (in this order) should be the Point")


class TestMinMult(unittest.TestCase):

    def setUp(self):
        self.e1 = Curve(4, 4, 5)
        self.p1 = Point(1, 2, self.e1)
        self.inf = Infinity()

    def test_inf(self):
        self.assertEqual(min_mult(self.inf), 0, "min_mult(infinity) == 0")

    def test_small(self):
        self.assertEqual(min_mult(self.p1), 4)  # Verified by hand


class TestPoint(unittest.TestCase):

    def setUp(self):
        self.e1 = Curve(4, 4, 5)
        self.e2 = Curve(4, 4, 7)
        self.p1 = Point(1, 2, self.e1)
        self.p2 = Point(4, 3, self.e1)
        self.p3 = Point(0, 2, self.e2)  # Verified by hand
        self.p4 = Point(1, 3, self.e1)
        self.inf = Infinity()

    def test_init_not_modulated(self):
        x, y = 26, -8
        p = Point(x, y, self.e1)  # self.p1 written funny
        self.assertEqual(p.x, self.p1.x, "Points constructed from x-coordinates not yet modulated should be modulated")
        self.assertEqual(p.y, self.p1.y, "Points constructed from y-coordinates not yet modulated should be modulated")

    def test_init_modulated(self):
        x, y = 1, 2
        self.assertEqual(self.p1.x, x, "Points constructed from modulated x-coordinates should have same x-coordinates")
        self.assertEqual(self.p1.y, y, "Points constructed from modulated y-coordinates should have same y-coordinates")

    def test_eq_inf(self):
        self.assertNotEqual(self.p1, self.inf, "Points should never equal Infinity")

    def test_eq_diff_curve(self):
        self.assertNotEqual(self.p1, self.p3, "Points are different curves should not be equal")

    def test_eq_diff_pts_same_curve(self):
        self.assertNotEqual(self.p1, self.p2, "Different points on the same curve should not be equal")

    def test_eq_same_pts_same_curve(self):
        self.assertEqual(self.p1, self.p1, "Same points on same curve should be equal")

    def test_repr(self):
        s = "(0, 2)"
        self.assertEqual(repr(self.p3), s)

    def test_add_finite_distinct(self):
        sum = Point(4, 2, self.e1)  # Verified by hand

        self.assertEqual(self.p1 + self.p2, sum)

    def test_add_finite_equal(self):
        sum = Point(2, 0, self.e1)  # Verified by hand

        self.assertEqual(self.p1 + self.p1, sum)

    def test_add_infinite_and_finite_RHS(self):
        sum = self.p1

        self.assertEqual(self.p1 + self.inf, sum)

    def test_add_infinite_and_finite_LHS(self):
        sum = self.p1

        self.assertEqual(self.inf + self.p1, sum)

    def test_add_infinite_slope(self):
        # Sufficient for points to share x-coordinate
        self.assertEqual(self.p1 + self.p4, self.inf, "Addition of points that make an infinite slope "
                                                      "(zero denominator) should create sum of Infinity")

    def test_is_finite_true(self):
        self.assertTrue(self.p1.is_finite(), "Points not instances of Infinity should be finite")

    def test_is_finite_false(self):
        self.assertFalse(self.inf.is_finite(), "Instances of Infinity should not be finite")


if __name__ == "__main__":
    unittest.main()
