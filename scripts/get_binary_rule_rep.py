def get_binary_rule_rep(rule_number: int,
                        rule_order: int) -> str:
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
