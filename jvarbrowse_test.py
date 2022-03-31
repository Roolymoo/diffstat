import unittest

from diffstat import _is_singular

from jvarbrowse import _gen_curves, j_invariant, gen_eq_jvar


class TestGenCurves(unittest.TestCase):
    def setUp(self):
        self.p = 5
        self.curves = _gen_curves(self.p)

    def test_nonsingular(self):
        for e in self.curves:
            self.assertFalse(_is_singular(*e))

    def test_coefficients(self):
        for e in self.curves:
            b, c, p = e
            self.assertLessEqual(0, b, "0 should be lower bound")
            self.assertLessEqual(0, c, "0 should be lower bound")
            self.assertLess(b, p, "p should be strict upper bound")
            self.assertLess(c, p, "p should be strict upper bound")

    def test_format(self):
        for e in self.curves:
            self.assertIsInstance(e, tuple)
            self.assertEqual(len(e), 3)
            self.assertEqual(e[-1], self.p)

    def test_len_upper_bd(self):
        # There are (self.p)^2 possible permutations of b, p
        self.assertLessEqual(len(list(self.curves)), pow(self.p, 2), "(self.p)^2 should be an upper bound")


class TestJInvariant(unittest.TestCase):
    def setUp(self):
        self.curve1 = 0, 5, 7
        self.curve2 = 13, 0, 17

    def test_basic(self):
        # Elliptic Curves Numbery Theory and Cryptography, p. 46
        self.assertEqual(j_invariant(*self.curve1), 0)
        self.assertEqual(j_invariant(*self.curve2), 11)  # 1728 % 17 == 11


class TestGenEqJvar(unittest.TestCase):
    def setUp(self):
        self.p = 5

    def test_eq_jvar(self):
        # Dependency j_invariant
        for e1, e2 in gen_eq_jvar(self.p):
            self.assertEqual(j_invariant(*e1), j_invariant(*e2))

    def test_format(self):
        for elem in gen_eq_jvar(self.p):
            self.assertIsInstance(elem, tuple)
            self.assertEqual(len(elem), 2)

            for e in elem:
                self.assertIsInstance(e, tuple)
                self.assertEqual(len(e), 3)
                self.assertEqual(e[-1], self.p)

    def test_repeats(self):
        seen = []
        for e1, e2 in gen_eq_jvar(self.p):
            self.assertNotIn((e1, e2), seen)
            self.assertNotIn((e2, e1), seen)
            self.assertNotEqual(e1, e2)
            seen.append((e1, e2))


if __name__ == "__main__":
    unittest.main()
