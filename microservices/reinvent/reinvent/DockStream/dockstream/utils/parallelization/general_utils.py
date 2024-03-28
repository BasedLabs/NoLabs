import math


def split_into_sublists(input_list, partitions=None, slice_size=None):
    if partitions is None and slice_size is None:
        raise ValueError("Either specify partitions or slice size.")

    return_list = []
    start_indices = []    # store the global index of the first ligand

    len_input = len(input_list)

    if partitions is not None:
        chunk_size = int(math.ceil(len_input / partitions))
    else:
        chunk_size = slice_size

    for i in range(0, len_input, chunk_size):
        start_indices.append(i)
        return_list.append(input_list[i:i + chunk_size])
    return start_indices, return_list


def get_progress_bar_string(done, total, prefix="", suffix="", decimals=1, length=100, fill='█'):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (done / float(total)))
    filledLength = int(length * done // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    return f"{prefix}|{bar}| {percent}% {suffix}"
