from collections import deque
from typing import List
import numpy as np

StringList = List[str]


def get_kmers(inpt_seq: str,
              k: int = 3,
              looping: bool = False) -> StringList:
    '''Breaks down the provided inut sequence in the appropriate length
    kmers for ingestion by a CA rule of a specific order (where order = k)
    eg. input = '101', k=3, looping = False
        kmers -> ['101']

    eg. input = '101', k=3, looping = True
        kmers -> ['110', '101', '011']

    eg. input = '10110', k = 3, looping = False
        kmers -> ['101', '011', '110']

    eg  input = '10110', k = 3, looping = False
        kmers -> ['010', '101', '011', '110', '101']

    Parameters:
    -----------
    inpt_seq: str
        sequence of input states to be processed
    k: int
        length of kmer to be produced
    looping: Bool
        loop sequence or not
        note: if looping: Count(kmers) = length of inpt_seq
              if not looping: Count(kmers) = (length of input_seq - k + 1)

    Returns:
    --------
    kmers: list of strings
        list of kmers as described above
    '''
    # trivial cases
    if k > len(inpt_seq):
        return None

    if looping:
        inpt_length = len(inpt_seq)
        inpt = deque(inpt_seq)
        shift = int(k / 2) if k % 2 else int(np.floor((k - 1) / 2))
        inpt.rotate(shift)

        kmers = []
        for _ in range(inpt_length):
            kmers.append(list(inpt)[:k])
            inpt.rotate(-1)

    else:
        kmers = zip(*[inpt_seq[idx:] for idx in range(k)])

    kmers = [''.join(kmer) for kmer in kmers]
    return kmers
