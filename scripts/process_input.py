from int_to_based_string import int_to_based_string
from get_kmers import get_kmers
from get_lookup_indicies import get_lookup_indicies


def process_input(inpt_seq: str,
                  rule_number: int,
                  states: str,
                  order: int = 3,
                  looping: bool = True) -> str:
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

    # get fixed CA rule indicies mapping
    lookup_indices = get_lookup_indicies(states=states, order=order)

    # get rule representation: rule sequence length = length of lookup_indicies
    rule_seq = int_to_based_string(number=rule_number,
                                   states=states,
                                   length=len(lookup_indices))

    # break input sequence into k-mers
    inpt_kmers = get_kmers(inpt_seq=inpt_seq,
                           k=order,
                           looping=looping)

    # translate inpt_kmers into lookup_indicies
    output_indicies = [lookup_indices[kmer] for kmer in inpt_kmers]
    output_seq = ''.join([rule_seq[idx] for idx in output_indicies])

    return output_seq
