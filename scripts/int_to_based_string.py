def int_to_based_string(number: int,
                        states: str,
                        length: int = 0) -> str:
    """Generates the quivalent string representation of an integer in the
    provided base states

    Parameters:
    -----------
    number: int
        integer value to be represented
    states: string
        ordered states to use for representation
    length: int
        required length of representation

    Returns:
    --------
    sequence: string
        based & justified representation of input number
    Raises:
    -------
    * None
    """

    base = len(states)
    fill_state = states[0]
    (div, mod) = divmod(number, base)
    if div > 0:
        sequence = int_to_based_string(div, states=states) + states[mod]
        return sequence.rjust(length, fill_state)
    return states[mod] if not length else states[mod].rjust(length, fill_state)
