from itertools import product
from collections import deque
import numpy as np


def get_lookup_indicies(inpt_states, order=3):
    '''Generates lookup indicies for the set of CA rules defined by the
    input_states & specific rule order.

    Parameters:
    -----------
    inpt_states: str
        set of states (in string form) that form the lookup indicies for
        any of the CA rules
        eg '01'
    order: int
        length of the permutations of the inpt_states
        eg for rule_order = 3, a possible index is '101'
        note: order is equivalent to the k-mer length of equal to order

    Returns:
    --------
    indicies_map: dict
        map of the possible input indicies to their index in the rule
        note: order of indices preserves the order of the input_states
        eg.
            - states = '01', order = 2 -> ['00', '01', '10', '11']
            - states = '10', order = 2 -> ['11', '10', '01', '00']

    Throws:
    -------

    '''
    kmers = product(inpt_states, repeat=order)
    indicies_map = {''.join(kmer) for kmer in kmers}
    return indicies_map


def get_binary_rule_rep(rule_number, rule_order):
    ''' Generates the correctly justified binary representation of a CA rule
    based on the rule order.
    note: This is only intended for binary states.
          A generalization will be made for more possible states
          in future notebooks.

    Parameters:
    -----------
    rule_number: int
        The number of the rule to be converted into binary
    order: int
        length of the permutations of the inpt_states
        eg for rule_order = 3, a possible index is '101'
        note: order is equivalent to the k-mer length of equal to order

    Returns:
    --------
    binary_rep: str
        binary representation of the rule_number

    Throws:
    -------
    AssertionError: if the requested rulenumber is larger than the possible
                    number of CA rules given the rule order
                    (assumbes binary states by default)
    '''
    base = 2
    rule_length = base ** rule_order
    rule_count = base ** rule_length

    assert rule_number < rule_count, f'For rule order = {rule_order} & \
    states = [0, 1], rule_number must be [0,{rule_count})'

    rule_rep = bin(rule_number)[2:]
    rule_rep = rule_rep.rjust(rule_length, '0')
    return rule_rep


def get_kmers(inpt_seq, k=3, looping=False):
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
    if looping:
        inpt_length = len(inpt_seq)
        inpt = deque(inpt_seq)
        shift = int(k / 2) if k % 2 else int(np.floor((k - 1) / 2))
        inpt.rotate(shift)

        kmers = []
        for _ in range(inpt_length):
            kmers.append(list(inpt)[k:])
            inpt.rotate(-1)

    else:
        kmers = zip(*[inpt_seq[idx:] for idx in range(k)])

    kmers = [''.join(kmer) for kmer in kmers]
    return kmers


def process_inpt(inpt_seq,
                 rule_number,
                 input_states='01',
                 order=3,
                 looping=True):
    '''
    Process the input sequence through the indicated CA rule and
    produces the corrsponding output sequence:
    notes:
        - if looping is selected, the output sequence will be the same
        length as the inputsequence and adequately alligned
        ! input states should not be changed, they will be used
        in a next version.
    Parameters:
    -----------
    inpt_seq: str
        sequence of input states to be processed
    rule_number: int
        The number of the rule to be used to convert the inpt_seq
    inpt_states: str
        set of states (in string form) that form the inpt_seq & CA Rule
    order: int
        length of the permutations of the inpt_states
        eg for rule_order = 3, a possible index is '101'
        note: order is equivalent to the k-mer length of equal to order
    looping: Bool
        loop sequence or not

    Returns:
    --------
    output_seq: str
        sequence of output states after the inpt_seq has been processed
        through the CA Rule
    '''

    # get rule representation (sequence of binary states)
    rule_seq = get_binary_rule_rep(rule_number=rule_number,
                                   rule_order=order)

    # break input sequence into k-mers
    inpt_kmers = get_kmers(inpt_seq=inpt_seq,
                           k=order,
                           looping=looping)

    # translate inpt_kmers into lookup_indicies
    lookup_indices = get_lookup_indicies(input_states='01', order=order)
    output_indicies = [lookup_indices[kmer] for kmer in inpt_kmers]
    output_seq = ''.join([rule_seq[idx] for idx in output_indicies])

    return output_seq
