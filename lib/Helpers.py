from typing import List


def rev_check(collection, val, min_value):
    collection_copy = collection.copy()
    for i in range(1, len(collection_copy) + 1):
        if len(collection_copy) == 1:
            if collection_copy[0] - val >= min_value:
                collection_copy[0] -= val
                return collection_copy
            return []

        if collection_copy[-i] - val >= 0:
            collection_copy[-i] -= val
        else:
            return []

    return collection_copy


def update_array(arr1, arr2):
    for index, item in enumerate(arr2):
        arr1[index] = arr2[index]


def separate_array_values(arr: List[int], min_value: int, max_value: int):
    for i in range(1, len(arr)):
        diff = arr[i] - arr[i - 1]
        a = min_value - diff
        half = a / 2
        m = False
        if diff < min_value:
            if arr[i] + half > max_value:
                half += (arr[i] + half) - max_value
                m = True
            col = rev_check(arr[0:i], half, min_value)
            if len(col) == 0:
                arr[i] = arr[i - 1] + min_value
            else:
                update_array(arr, col)
                if m:
                    arr[i] = max_value
                else:
                    arr[i] = arr[i] + half
    return arr


def get_center_value(collection: List[List[str]]) -> List[int]:
    output = []

    for item in collection:
        center_value = int(item[0]) + (int(item[1]) - int(item[0]))
        output.append(center_value)

    return output


def get_fasta_length(file_name: str):
    length = 0

    with open(file_name, 'r') as file:
        for line in file:
            if line.startswith('>'):
                continue

            line = line.rstrip()

            length += len(line)

    return length


def change_str_to_int(collection: List[List[str]], indexes_to_change: List[int]):
    output = []

    for line_splitted in collection:
        line_mod = []
        for index, item in enumerate(line_splitted):
            if index in indexes_to_change:
                line_mod.append(int(line_splitted[index]))
            else:
                line_mod.append(line_splitted[index])
        output.append(line_mod)

    return output

# _min_distance = 10
# _max_distance = 100
# _data = [2,5,8]
#
#
# ou = separate_array_values(_data, _min_distance, _max_distance)
# print(ou)
