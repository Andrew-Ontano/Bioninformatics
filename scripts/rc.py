#!/usr/bin/env python3

import sys
from Bio.Seq import Seq

try:
    sequence = Seq(sys.argv[1])
    print(sequence.reverse_complement())
except ValueError:
    print('-')