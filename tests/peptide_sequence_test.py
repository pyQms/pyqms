import pyqms
from pathlib import Path

tests_path = Path(__file__).parent
app_path = tests_path.parent

unimod_file_list = [
    tests_path.joinpath("data", "test_usermod.xml"),
    app_path.joinpath("pyqms", "kb", "ext", "unimod.xml"),
]

TESTS = [
    {
        "input": "TESTPEPTIDE",
        "output": {
            "cc": "C(50)H(79)N(11)O(24)",
            "abun": [52441, 31423, 11803, 3313, 763, 134, 23, 2, 0, 0, 0],
            "tmzs": {609769, 609770, 609771, 609772, 609773, 609774, 609775},
        },
    },
    {
        "input": "TESTPEPTIDE#Oxidation:1",
        "output": {
            "cc": "C(50)H(79)N(11)O(25)",
            "abun": [52314, 31367, 11894, 3374, 787, 140, 24, 3, 0, 0, 0],
            "tmzs": {617767, 617768, 617769, 617770, 617771, 617772, 617773},
        },
    },
    {
        "input": "TESTPEPTIDE#SILAC TMT:0",
        "output": {
            "cc": "C(52)H(99)13C(10)15N(1)N(12)O(26)",
            "abun": [12, 2247, 49894, 30980, 12102, 3539, 849, 156, 28, 3, 0, 0, 0],
            "tmzs": {726859, 726860, 726861, 726862, 726863, 726864, 726865, 726866},
        },
    },
]


def peptide_sequence_test():
    for test in TESTS:
        yield sequence_to_isotopes, test["input"], test["output"]


def sequence_to_isotopes(sequence, results):
    lib = pyqms.IsotopologueLibrary(
        molecules=[sequence], charges=[2], unimod_files=unimod_file_list
    )
    assert results["cc"] in lib
    iso_data = lib[results["cc"]]["env"][(("N", "0.000"),)]
    assert iso_data["abun"] == results["abun"]
    # test first non-None set of tmzs
    for i in range(len(iso_data[2]["tmzs"])):
        if iso_data[2]["tmzs"][i] is not None:
            break
    assert iso_data[2]["tmzs"][i] == results["tmzs"]
