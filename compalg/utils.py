"""
Utils file.
"""


def generate_list(generator, first_element, limit):
    element = first_element
    result = []

    for _ in range(limit):
        result.append(element)
        element = generator(element)

    return result


def filter_with_mask(mask, elements):
    mask_bits = list(bin(mask)[2:])
    return map(lambda p: p[1], filter(lambda p: p[0] == '1', zip(mask_bits, elements)))
