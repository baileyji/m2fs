def typeClearance(type):
    type=type.lower()
    if type=='c':
        return 25.0
    elif type=='g':
        return 33.0
    elif type in ['o','s','a']:
        return 6.5
    else:
        raise ValueError

def min_space(kind, other):
    """ types are C, S, O, A, G
        returns the minimum allowable seperation in asec
        """
    return typeClearance(kind)+typeClearance(other)