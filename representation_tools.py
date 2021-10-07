from tabulate import tabulate


def format_numpy_array_digits(numpy_array,
                              decimal_spaces,
                              should_near_zero_be_adjusted_to_low_value=True,
                              do_the_round=False):
    the_divisor = 10 ** decimal_spaces
    the_low_limit = 10 ** -decimal_spaces
    the_return = numpy_array * the_divisor // 1 / the_divisor
    if do_the_round:
        the_divisor = the_divisor * 10
        the_return = ((numpy_array * the_divisor % 10) >= 5) \
            * the_low_limit \
            + the_return
    if should_near_zero_be_adjusted_to_low_value:
        the_low_limit_compare = the_low_limit if not do_the_round else the_low_limit / 2
        the_return = (numpy_array < the_low_limit_compare) * (numpy_array > 0) \
            * the_low_limit \
            + the_return
    return the_return


def print_data_in_table(data, top_down_headers=None, left_right_headers=None, dimensions=1):
    if left_right_headers is None:
        left_right_headers = ["data"]
    if dimensions < 1 or dimensions > 2:
        print(f'{dimensions} dimensions are not supported in print_data_in_table.')
        return
    if dimensions == 1:
        if top_down_headers is None:
            print(tabulate(data, left_right_headers))
        else:
            new_data = []
            for index in range(len(data)):
                new_data.append([top_down_headers[index], data[index]])
            print(tabulate(new_data, ["name"] + left_right_headers))
        return
    elif dimensions == 2:
        # we substract one because top_down_headers will be added
        if len(data) % (len(left_right_headers) - 1) > 0:
            print(f'Given data cannot be converted to {len(left_right_headers)} column table.')
            return
        elif len(left_right_headers) == 1:
            print(f'1 column header can not form 2 dimensional table.')
            return
        if top_down_headers is None:
            print(f'top down headers are must for 2 dimensional table. Else use 1 dimensional')
            return
        new_data = []
        the_index = 0
        the_index_of_header = 0
        the_length = len(left_right_headers) - 1
        while the_index + the_length < len(data)+1:
            new_data.append(
                [top_down_headers[the_index_of_header]]
                + data[the_index: the_index + the_length])
            the_index += the_length
            the_index_of_header += 1
        print(tabulate(new_data, left_right_headers))
        return


def print_in_phylip(phylip_array):
    the_row_count = phylip_array[0]
    print(the_row_count)
    the_start = True
    for the_row in phylip_array:
        if the_start:
            the_start = False
            continue
        print(' '.join(the_row))
