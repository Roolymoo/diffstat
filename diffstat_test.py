import unittest
import os
import os.path
from random import randint

from sympy import isprime

from diffstat import _is_singular, _avg, _gen_x_coords, _write_curve, _write_x_coord, _rand_prime, gen_rand_curves, \
    max_distance, arg_parser


def _read_curves(file_name):
    """(str) -> tuple
    Returns (generating args, curve_list), read from file_name.

    Assumes file_name contains its extension, that it exists, and is formatted as follows:

    <coeff_min> <coeff_max> <prime_min> <prime_max> <num_curves>
    <b1> <c1> <p1>
    <b2> <c2> <p2>
    ...

    Generating args is a tuple of the form (coeff_min, coeff_max, prime_min, prime_max, num_curves).
    """
    data = lambda line: tuple(int(i) for i in line.split())
    curve_list = []
    with open(file_name, "r") as file:
        args = data(file.readline())
        for curve_data in file:
            e = data(curve_data)
            curve_list.append(e)

    return args, curve_list


def _write_rand_curves(file_name, args):
    """(str, tuple) -> NoneType
    Writes to file_name the curves given by gen_rand_curves(*args).

    Assumes file_name contains its file extension. Assumes the curves in curves are of the form (b, c, p). file_name
    does not need to exist, and it will be overwritten if it does exist. See _read_curves() for the output format.
    """
    with open(file_name, "w") as file:
        file.write("{} {} {} {} {}".format(*args))

        curves = gen_rand_curves(*args)
        for e in curves:
            file.write("\n{} {} {}".format(*e))


class TestReadCurves(unittest.TestCase):
    def setUp(self):
        self.file1 = "test_read_curves_basic.txt"

    def test_basic(self):
        args = 1, 1000, 4, 1000, 10
        curve_list = [(861, 533, 433), (596, 512, 269), (136, 576, 151), (594, 75, 179), (171, 272, 853),
                      (707, 951, 73), (845, 983, 691), (177, 872, 389), (166, 732, 67), (509, 539, 947)]

        self.assertEqual((args, curve_list), _read_curves(self.file1))


class TestWriteRandCurves(unittest.TestCase):
    def setUp(self):
        self.args1 = 1, 100, 4, 100, 5
        self.file1 = "test_write_rand_curves_temp{}.txt".format(randint(100000, 200000))

    def tearDown(self):
        os.remove(self.file1)

    def test_format(self):
        data = lambda s: tuple(int(i) for i in s.split())

        _write_rand_curves(self.file1, self.args1)
        with open(self.file1, "r") as file:
            args = data(file.readline())
            self.assertEqual(args, self.args1)

            for curve_data in file:
                curve = data(curve_data)
                self.assertEqual(len(curve), 3)


class TestWriteCurve(unittest.TestCase):
    def setUp(self):
        self.curves1 = [(1, 2, 3), (4, 5, 6), (7, 8, 9)]
        self.file1 = "test_write_curve1_temp{}.txt".format(randint(100000, 200000))

    def tearDown(self):
        os.remove(self.file1)

    def test_basic(self):
        data = lambda s: tuple(int(i) for i in s.split())

        # Write the curves with newlines separating them
        with open(self.file1, "w") as file:
            for e in self.curves1:
                _write_curve(file, e)
                file.write("\n")

        # Test to see if written correctly
        curves = []
        with open(self.file1, "r") as file:
            for line in file:
                curves.append(data(line))

        self.assertEqual(curves, self.curves1)


class TestWriteXCoord(unittest.TestCase):
    def setUp(self):
        self.x_coords1 = (1, 2, 3, 4, 5, 6)
        self.file1 = "test_write_x_coords1_temp{}.txt".format(randint(100000, 200000))

    def tearDown(self):
        os.remove(self.file1)

    def test_basic(self):
        # Note: Dependency on _write_curve()
        # Write the x--coordinates with a dummy curve already written
        with open(self.file1, "w") as file:
            _write_curve(file, (0, 0, 0))
            for x in self.x_coords1:
                _write_x_coord(file, x)

        # Read the results
        with open(self.file1, "r") as file:
            results = file.readline()
            line_after_results = file.readline()

        self.assertEqual(results, "0 0 0 1 2 3 4 5 6")
        self.assertEqual(line_after_results, "")


class TestIsSingular(unittest.TestCase):
    def setUp(self):
        # Non-singular
        self.e1 = (861, 533, 433)
        # Singular
        # Any curve E: y^2 = x^3 - 3x + 2 mod p works
        self.e2 = (8, 2, 11)

    def test_is_nonsingular(self):
        self.assertFalse(_is_singular(*self.e1))

    def test_is_singular(self):
        self.assertTrue(_is_singular(*self.e2))


class TestAvg(unittest.TestCase):
    def setUp(self):
        self.pts1 = [1, 2, 3, 4, 5]
        self.pts2 = [1, 2, 3, 4, 5, 6]

    def test_int(self):
        self.assertEqual(_avg(self.pts1), 3.0)

    def test_nonint(self):
        self.assertEqual(_avg(self.pts2), 3.5)

    def test_zero1(self):
        self.assertEqual(_avg([0, 0, 0]), 0.0)

    def test_zero2(self):
        self.assertEqual(_avg([-5, 5]), 0.0)


class TestGenXCoords(unittest.TestCase):
    def setUp(self):
        self.e1 = 2, 4, 5
        self.e2 = 861, 533, 433
        self.e3 = 4, 4, 5
        # From Trappe and Washington pp. 348--349
        self.e1_x_coords = [0, 2, 4]
        # From Trappe and Washington pp. 352--353
        self.e3_x_coords = [0, 1, 2, 4]
        self.file1 = "test_gen_x_coords_temp{}.txt".format(randint(100000, 200000))

    def test_unique(self):
        x_coords = list(_gen_x_coords(*self.e2))
        self.assertCountEqual(x_coords, set(x_coords))

    def test_basic1(self):
        self.assertEqual(list(_gen_x_coords(*self.e1)), self.e1_x_coords)

    def test_basic2(self):
        self.assertEqual(list(_gen_x_coords(*self.e3)), self.e3_x_coords)

    def test_single_curve_file(self):
        # Test _gen_x_coords() writing x--coordinates for one curve
        # Note: Dependency on _write_curve(), _write_x_coord()

        # Have _gen_x_coords write the x--coordinates
        with open(self.file1, "w") as file:
            _write_curve(file, self.e1)
            b, c, p = self.e1
            x_coords = _gen_x_coords(b, c, p, file)
            # Iterate through x_coords to have the x--coordinates written
            for x in x_coords:
                pass

        # Read the results
        with open(self.file1, "r") as file:
            results = file.readline()
            line_after_results = file.readline()

        # Cleanup
        os.remove(self.file1)

        self.assertEqual(results, "2 4 5 0 2 4\n")
        self.assertEqual(line_after_results, "")

    def test_multilpe_curve_file(self):
        # Test _gen_x_coords() writing x--coordinates for multiple curves curve
        # Note: Dependency on _write_curve(), _write_x_coord()

        # Have _gen_x_coords write the x--coordinates
        with open(self.file1, "w") as file:
            for e in self.e1, self.e3:
                _write_curve(file, e)
                b, c, p = e
                x_coords = _gen_x_coords(b, c, p, file)
                # Iterate through x_coords to have the x--coordinates written
                for x in x_coords:
                    pass

        # Read the results
        with open(self.file1, "r") as file:
            results_e1 = file.readline()
            results_e3 = file.readline()
            line_after_results = file.readline()

        # Cleanup
        os.remove(self.file1)

        self.assertEqual(results_e1, "2 4 5 0 2 4\n")
        self.assertEqual(results_e3, "4 4 5 0 1 2 4\n")
        self.assertEqual(line_after_results, "")


class TestRandPrime(unittest.TestCase):
    def setUp(self):
        # (lower bd, upper bd, num primes)
        self.ranges = (4, 5, 10), (4, 500, 100), (5, 500000, 1000)

    def test_is_prime(self):
        for lower_bd, upper_bd, num_primes in self.ranges:
            for i in range(num_primes):
                p = _rand_prime(lower_bd, upper_bd)
                self.assertTrue(isprime(p),
                                "args={}, {}, {}, prime={}".format(lower_bd, upper_bd, num_primes, p))

    def test_is_prime_equal_bds(self):
        bds = ((100, 100), (158, 158), (167, 167), (10000, 10000))
        for bd_pair in bds:
            p = _rand_prime(*bd_pair)
            self.assertTrue(isprime(p), "bd pair: ({}, {})".format(*bd_pair))

    def test_range(self):
        # Assert that the generated random primes are within the specified range
        for lower_bd, upper_bd, num_primes in self.ranges:
            for i in range(num_primes):
                p = _rand_prime(lower_bd, upper_bd)
                self.assertLessEqual(lower_bd, p)
                self.assertLess(p, 2 * upper_bd)

    def test_repeats(self):
        # Ensure primes are not repeated often (for large enough ranges)
        lower_bd, upper_bd, num_primes = self.ranges[-1]
        # Do a few runs to try to ensure accuracy (e.g. no flukes)
        for i in range(5):
            unique_primes = set(_rand_prime(lower_bd, upper_bd) for i in range(num_primes))
            self.assertGreaterEqual(len(unique_primes), num_primes * 0.75, "> 25% repeats")

    def test_random(self):
        # Ensure randomness (e.g. one batch gathered won't be too similar to another batch)
        lower_bd, upper_bd, num_primes = self.ranges[-1]
        # Do a few runs to try to ensure accuracy (e.g. no flukes)
        for i in range(5):
            primes1 = set(_rand_prime(lower_bd, upper_bd) for i in range(num_primes))
            primes2 = set(_rand_prime(lower_bd, upper_bd) for i in range(num_primes))
            common_primes = primes1 & primes2
            self.assertLess(len(common_primes), num_primes * 0.25, ">= 25% primes shared")


# Test classes for gen_rand_curves are split up due to computational intensity of the function and setUp() methods
class TestGenRandCurvesRandom(unittest.TestCase):
    def setUp(self):
        self.args = 1, 1000000, 4, 1000000, 1000

    def test_random(self):
        # Ensure generated batches of curves are sufficiently different
        num_curves = self.args[-1]
        for i in range(5):
            batch1 = set(gen_rand_curves(*self.args))
            batch2 = set(gen_rand_curves(*self.args))
            common_curves = batch1 & batch2
            self.assertLess(len(common_curves), num_curves * 0.25)


# Test classes for gen_rand_curves are split up due to computational intensity of the function and setUp() methods
class TestGenRandCurvesBasic(unittest.TestCase):
    def setUp(self):
        args1, curve_list1 = _read_curves("test_rand_curves_basic1_TEMP.txt")
        args2, curve_list2 = _read_curves("test_rand_curves_basic2_TEMP.txt")
        args3, curve_list3 = _read_curves("test_rand_curves_basic3_TEMP.txt")

        self.args_curves = (args1, curve_list1), (args2, curve_list2), (args3, curve_list3)

    def test_num_curves(self):
        for args, curve_list in self.args_curves:
            num_curves = args[-1]
            self.assertEqual(num_curves, len(curve_list))

    def test_ranges(self):
        for args, curve_list in self.args_curves:
            coeff_min, coeff_max, prime_min, prime_max = args[:-1]

            for e in curve_list:
                b, c, p = e
                self.assertGreaterEqual(b, coeff_min, "b < coeff_min")
                self.assertGreaterEqual(c, coeff_min, "c < coeff_min")
                self.assertLessEqual(b, coeff_max, "b > coeff_max")
                self.assertLessEqual(c, coeff_max, "c > coeff_max")
                self.assertGreaterEqual(p, prime_min, "p < prime_min")
                self.assertLess(p, prime_max, "c >= prime_max")

    def test_primes(self):
        for args, curve_list in self.args_curves:
            for e in curve_list:
                p = e[-1]
                self.assertTrue(isprime(p))

    def test_nonsingular(self):
        # Note: Dependency on _is_singular()
        for args, curve_list in self.args_curves:
            for e in curve_list:
                self.assertFalse(_is_singular(*e))


# TODO: Another (larger) basic test
class TestMaxDistance(unittest.TestCase):
    # Note: Dependency on _gen_x_coords()

    def setUp(self):
        self.e1 = 2, 4, 5
        self.e2 = 4, 4, 5
        self.e3 = 1, 1, 3
        self.file1 = "test_max_distance_temp{}.txt".format(randint(100000, 200000))

    def test_two_pts(self):
        x_coords = _gen_x_coords(*self.e3)
        self.assertEqual(max_distance(x_coords), 1)

    def test_basic1(self):
        x_coords = _gen_x_coords(*self.e1)
        self.assertEqual(max_distance(x_coords), 2)

    def test_basic2(self):
        x_coords = _gen_x_coords(*self.e2)
        self.assertEqual(max_distance(x_coords), 2)


class TestArgParser(unittest.TestCase):
    def setUp(self):
        self.args = ["1", "1000", "4", "1000", "100"]

    def test_basic(self):
        parser = arg_parser()
        parsed_args = parser.parse_args(self.args)
        self.assertEqual([int(i) for i in self.args], [parsed_args.coeff_min, parsed_args.coeff_max,
                                                       parsed_args.prime_min, parsed_args.prime_max,
                                                       parsed_args.num_curves])


if __name__ == "__main__":
    # Write output of gen_rand_curves beforehand for testing use (computationally intensive to have setUp generate
    # this for each testing function, so just do it once now)
    ARGS_BASIC1 = 1, 1000, 4, 1000, 10
    ARGS_BASIC2 = 1, 100000, 4, 100000, 1000
    ARGS_BASIC3 = 1, 1000000, 4, 1000000, 1000
    BASIC_NAME1 = "test_rand_curves_basic1_TEMP.txt"
    BASIC_NAME2 = "test_rand_curves_basic2_TEMP.txt"
    BASIC_NAME3 = "test_rand_curves_basic3_TEMP.txt"
    _write_rand_curves(BASIC_NAME1, ARGS_BASIC1)
    _write_rand_curves(BASIC_NAME2, ARGS_BASIC2)
    _write_rand_curves(BASIC_NAME3, ARGS_BASIC3)

    # Don't exit so can clean up after
    unittest.main(exit=False)

    # Clean up
    os.remove(BASIC_NAME1)
    os.remove(BASIC_NAME2)
    os.remove(BASIC_NAME3)
