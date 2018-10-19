'''
utils.py

Module with various useful functions

---
Created: Vadim Pribylov, 2018-10-19
'''


def numlist(itemlist):
    '''
    Returns numbered list
    '''
    i = 0
    try:
        for item in itemlist:
            print(str(i) + '\t' + item)
    #        print(str_format.format(i, item)
            i += 1
    except TypeError:
        raise TypeError(
                'expected itemlist to be list-type, '
                'got {}'.format(type(itemlist))
            )


def tableconvert(sep='\t|', *args):
    '''
    Returns text table
    '''
    table = []
    for item_row in args:
        row = ''
        for item_col in item_row:
            if type(item_col) == str:
                row = row + item_col + '{}'.format(sep)
            else:
                row = row + '{0:3.2f}{1}'.format(item_col, sep)
        table.append(row)
    return table


def mdtableconvert(header, rownames, *args):
    '''
    Returns markdown text table

    Args:
        header:list-type - table header
        rownames:list-type - names of the row to place in the first column
    '''
    table = []
    row = '| '
    for item in header:
        row += str(item) + ' | '
    table.append(row)
    row = '| '
    for item in header:
        row += ' --- |'
    table.append(row)
    for item_row, item_rowname in zip(args, rownames):
        row = '| '
        row += str(item_rowname) + ' | '
        for item_col in item_row:
            row += '{:3.4f}'.format(item_col)
            row += ' | '
        table.append(row)
    return table


def print_csv(*args):
    '''
    Returns csv data from input
    '''
    table = []
    for item_row in args:
        row = ''
        for item_col in item_row:
            if isinstance(item_col, str):
                pass
            else:
                item_col = str(item_col)
            row = row + item_col + ','
        table.append(row[:-1])  # truncate trailing separator
    return table


def loctime2int(hours, minutes, seconds, mod):
    '''
    Converts time to second count from midnight

    mod affects the second division

    reverse function is `int2loctime`
    '''
    if minutes > 60 or seconds > 60 or hours > 24:
        raise AttributeError('incorrect time format')
    if minutes < 0 or seconds < 0 or hours < 0:
        raise AttributeError('negative time error')
    return mod * 3600 * hours + mod * 60 * minutes + mod * seconds


def int2loctime(time, mod):
    '''
    Converts second count from midnight into time
    returns (hours, minutes, seconds)

    mod affects the second division

    reverse function is `loctime2int`
    '''
    hours = time // (3600*mod)
    minutes = (time - hours*3600*mod) // (60*mod)
    seconds = (time - hours*3600*mod - minutes*60*mod) / mod
    return (hours, minutes, seconds)
