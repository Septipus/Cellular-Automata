def based_string_to_int(sequence: str,
                        states: str) -> int:
    """Generates integer value of provided sequence and states

    Parameters:
    -----------
    sequence: string
        representation of input number
    states: string
        ordered states used to make sequence

    Returns:
    --------
    number: int
        integer value of sequence

    Raises:
    -------
    * None
    """

    base = len(states)
    digits = [states.index(val) for val in sequence[::-1]]
    number = 0
    for idx, digit in enumerate(digits):
        number += digit * (base ** idx)
    return number
