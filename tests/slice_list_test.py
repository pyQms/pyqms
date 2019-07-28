#!/usr/bin/env python
# encoding: utf-8
"""

Test _slice_list

#NOTE bisect.bisect uses first the first value of a tuple for sorting and then
the second if the first value is completely identical!

"""
import pyqms
import unittest

SPECTRUM = [
    (1000, 1337),
    (1001, 1337),
    (1002, 1337),
    (1003, 1337),
    (1004, 1337),
    (1005, 1337),
]

TESTS = [
    {"input": [(1002, 0.001), (1003, 0.001)], "output": 5},
    {"input": [(1002, 0.001), (1003, 1337.001)], "output": 6},
    {"input": [(999, 0.001), (1001, 0.001)], "output": 3},
    {"input": [(999, 0.001), (999, 0.001)], "output": 2},
    {"input": [(1006, 0.001), (1006, 0.001)], "output": 2},
]

MOLECULES = ["TEST"]
CHARGES = [2]
lib = pyqms.IsotopologueLibrary(molecules=MOLECULES, charges=CHARGES)


def extend_kb_with_fixed_labels_test():
    for test_dict in TESTS:
        yield checker_function, test_dict


def checker_function(test_dict):
    # with a higher second tuple value, bisect selects the next position...
    assert len(lib._slice_list(SPECTRUM, test_dict["input"])) == test_dict["output"]
    return


class TestResults(unittest.TestCase):
    def setUp(self):
        pass

    def crash_test(self):
        """
        Testing to slice with an tuple which can not be ordered
        # __Modified CF__ ... tuple can be slice now
        """
        pass
        # with self.assertRaises(SystemExit) as system_exit_check:
        #     lib._slice_list(
        #         SPECTRUM,
        #         [ (1002, 0.001), (1003, 0.001)]
        #     )
        # self.assertEqual(
        #     system_exit_check.exception.code,
        #     1
        # )


if __name__ == "__main__":
    slice_test()
