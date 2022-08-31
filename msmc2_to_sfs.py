# from msmc2 to unfolded sfs
import argparse
import sys
import pytest

def parse_args(args):
    # Instantiate the parser

    parser = argparse.ArgumentParser(
        description="""
        Convert an MSMC2 file to an unfolded SFS 
        using alleles from least one outgroup.
        """
    )
    parser.add_argument(
        "--msmc2_file",
        type=str,
        required=True,
        help="""
        Path to the msmc2 file
        """,
    )
    parser.add_argument(
        "--allele_total",
        type=int,
        required=True,
        help="""
        The total allele count to be expected in the final column. 
        """,
    )
    parser.add_argument(
        "-g",
        "--outgroup_index",
        type=int,
        nargs="*",
        required=True,
        help="""
        The zero indexed position(s) of allele(s) to be treated as the outgroup.
        The two alleles must match for the locus to be counted in the SFS.
        Positions must be space separated. 
        """,
    )
    return parser.parse_args()


class SFS:
    def __init__(self, allele_total, outgroup_index):
        self.allele_total = allele_total
        self.outgroup_index = outgroup_index
        self.sfs = [0] * (1 + allele_total - len(outgroup_index))

    def process_line(self, line):
        chrom, pos, invariant, alleles, *extra = line.split()
        if len(extra) > 0:
            raise ValueError(f"chromosome {chrom} position {pos} has too many columns.")
        if len(alleles) != self.allele_total:
            raise ValueError(
                f"chromosome {chrom} position {pos} has "
                f"found {len(alleles)} alleles expected {self.allele_total}."
            )
        anc_allele = get_anc(alleles, self.outgroup_index)
        bi_alleles = get_alleles(alleles)
        if anc_allele and bi_alleles:
            derived = list(bi_alleles - anc_allele)[0]
            self.sfs[alleles.count(derived)] += 1
            self.sfs[0] += int(invariant) - 1


# helper functions
def get_alleles(alleles):
    if len(set(alleles)) == 2:
        return set(alleles)
    return None


def get_anc(alleles, anc_indx):
    anc_set = set([alleles[i] for i in anc_indx])
    if len(anc_set) == 1:
        return anc_set
    else:
        return None


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    sfs_data = SFS(args.allele_total, args.outgroup_index)
    with open(args.msmc2_file, "r") as f:
        for line in f:
            sfs_data.process_line(line.strip())
    print(" ".join([str(i) for i in sfs_data.sfs]))

##############
# UNIT TESTS #
##############
class TestClass:
    def test_get_alleles1(self):
        alleles = "AAAAAA"
        assert get_alleles(alleles) == None

    def test_get_alleles2(self):
        alleles = "AAAATAA"
        assert get_alleles(alleles) == {"A", "T"}

    def test_get_anc1(self):
        alleles = "TATATAA"
        anc_indx = [0, 2, 4]
        assert get_anc(alleles, anc_indx) == {"T"}

    def test_get_anc2(self):
        alleles = "TATAAAA"
        anc_indx = [0, 2, 4]
        assert get_anc(alleles, anc_indx) == None

    def test_process_line(self):
        test_data = SFS(5, [0, 1])
        line = "8       264084  63      TTAAA"
        test_data.process_line(line.strip())
        line = "8       264084  63      TTATA"
        test_data.process_line(line.strip())
        assert test_data.sfs == [124, 0, 1, 1]

    def test_wronginput(self):
        with pytest.raises(
            ValueError, match="chromosome 1 position 8 has too many columns."
        ):
            test_data = SFS(5, [0, 1])
            line = "1       8       264084  63      TAAA"
            test_data.process_line(line.strip())

    def test_wrongsize(self):
        with pytest.raises(
            ValueError,
            match="chromosome 8 position 264084 has found 4 alleles expected 5.",
        ):
            test_data = SFS(5, [0, 1])
            line = "8       264084  63      TAAA"
            test_data.process_line(line.strip())
