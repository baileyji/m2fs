from jbastro.misc import derangify

def get_seqnos(listfile):
    ret = set()
    with open(listfile, 'r') as lf:
        for l in lf:
            if l[0] in '1234567890':
                range = l.split()[0]
                ret.update(map(str, derangify(range)))
            elif len(l) > 1 and l[0] in 'RBrb' and l[1] in '1234567890':
                range = l[1:].split()[0]
                ret.update(map(lambda x: l[0].lower() + str(x), derangify(range)))
    return list(ret)


def get_seqnos_advanced(listfile):
    ret = set()
    datestr = ''
    crmap = {}
    with open(listfile, 'r') as lf:
        for l in lf:
            images = []
            docr = 'docr' in l.lower().split()
            if l.startswith('#'):
                newdate = parse_date_comment(l)
                if newdate is not None:
                    datestr = newdate.strftime('ut%Y%m%d')
            elif l[0] in '1234567890':
                nums = derangify(l.split()[0])
                images = [(datestr, '{}{:04}'.format(color, x)) for color in 'rb' for x in nums]
            elif len(l) > 1 and l[0] in 'RBrb' and l[1] in '1234567890':
                nums = derangify(l[1:].split()[0])
                color = l[0].lower()
                images = [(datestr, '{}{:04}'.format(color, x)) for x in nums]
            for i in images:
                crmap[i] = docr
            ret.update(images)
    return [(i[0], i[1], crmap[i]) for i in ret]
