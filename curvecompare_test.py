import unittest
from os.path import join

from curvecompare import _read_curves, collect_curves, pairwise_duplicates


class TestReadCurves(unittest.TestCase):
    def setUp(self):
        # This file does have a duplicated curve
        self.file_name = join("test_data", "test_iterations", "1", "curves_1_10_4_10_10_2015-07-12_15-54-36.471535.txt")
        self.curves = _read_curves(self.file_name)

    def test_all_curves_read(self):
        expected_curves = [(2, 6, 7), (2, 10, 7), (6, 4, 7), (4, 5, 5), (6, 4, 7), (3, 10, 7), (6, 10, 7), (1, 6, 5),
                           (9, 2, 5), (5, 7, 7)]
        # Order that the curves are read should be the same
        self.assertListEqual(expected_curves, self.curves)


class TestCollectCurves(unittest.TestCase):
    def setUp(self):
        # At least one of the test iteration curve files has a duplicated curve
        self.iterations_path = join("test_data", "test_iterations")
        self.iterations_dirs = ["1", "2", "3"]
        self.collected_curves = collect_curves(self.iterations_path, self.iterations_dirs)

    def test_all_dirs_read(self):
        # Dependency curvecompare._read_curves
        dirs = [int(i) for i in self.iterations_dirs]
        dirs.sort()

        dirs_read = list(self.collected_curves.keys())
        dirs_read.sort()

        self.assertListEqual(dirs, dirs_read)

    def test_all_curves_read(self):
        # Dependency curvecompare._read_curves
        # This also tests that the curves_1_10_ files are not read, as they should not be
        self.maxDiff = None
        expected_curves = {}
        for i in self.iterations_dirs:
            expected_curves.update({int(i): set(_read_curves(join(self.iterations_path, i, "all_curves.txt")))})

        self.assertDictEqual(expected_curves, self.collected_curves)


class TestPairwiseDuplicates(unittest.TestCase):
    def setUp(self):
        self.iterations_path = join("test_data", "test_iterations")
        self.iterations_dirs = ["1", "2", "3"]
        self.iterations = collect_curves(self.iterations_path, self.iterations_dirs)
        self.iter_indices = [1, 2, 3]
        self.similarities = pairwise_duplicates(self.iterations, self.iter_indices)

    def test_keys(self):
        keys = set(self.similarities.keys())
        self.assertSetEqual(set(self.similarities.keys()), {1, 2, 3})

    def test_values_keys(self):
        # Assert all pairwise combinations are done
        self.assertSetEqual(set(self.similarities[1].keys()), {2, 3})
        self.assertSetEqual(set(self.similarities[2].keys()), {3})
        self.assertSetEqual(set(self.similarities[3].keys()), set())

    def test_duplicates(self):
        # Assert the results are correct
        self.assertSetEqual(self.similarities[1][2], {(80, 92, 11)})
        self.assertSetEqual(self.similarities[1][3], set())
        self.assertSetEqual(self.similarities[2][3], set())


if __name__ == "__main__":
    unittest.main()
