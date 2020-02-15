from itertools import product


def get_lookup_indicies(states: str,
                        order: int = 3) -> dict:
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
    kmers = product(states, repeat=order)
    indicies_map = {''.join(kmer): idx for idx, kmer in enumerate(kmers)}
    return indicies_map
