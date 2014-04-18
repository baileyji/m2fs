from plate import Plate
from jbastro.astroLib import sexconvert

ROTATOR_SETTING='-7.18'

def target(name, ra, de, ep=2000.0):
    """name, 'hh mm ss.s', 'dd.mm.ss.s', ep=2000.0"""
    return {'name':name,'ra':ra,'de':de,'epoch':ep}





def generate_summary_file(platefiles):
    tlist=[]
    platerec='{name:<10} {ns:>2} {offset:>5}\n'
    setuprec=('     {name:<11} {ra:<12} {de:<12} '
             '{epoch:<11} {sidereal_time:<11} {airmass:<11} {n:<11}\n')

    file_lines=[]

    for f in platefiles:
        try:
            p=Plate(f)
            if p.file_version=='0.1':
                raise Exception("Can't process v0.1"+f)
        except Exception, e:
            print 'Platefile Error: {}'.format(e)
            continue

        file_lines.append(platerec.format(name=p.name, ns=p.n_setups,
                                          offset=p.standard_offset))

        file_lines.append(setuprec.format(name='Name', ra='RA',
                                          de='DE',epoch='Epoch',
                                          sidereal_time='ST',
                                          airmass='Air', n='N'))
        for sname, setup in p.setups.iteritems():
            file_lines.append(setuprec.format(**setup.attrib))
            attrib=setup.attrib.copy()
            attrib['name']=p.name+' '+attrib['name']
            tlist.append(attrib)
    return file_lines, tlist



def write_summary_file(sfile, platefiles):
    file_lines, tlist=generate_summary_file(platefiles)
    with open(sfile,'w') as fp:
        for l in file_lines:
            fp.write(l)

def generate_tlist_file(recs):
    file_lines=[]
    """rec iterable of dicts with 'name' 'ra', 'de', & 'epoch'"""
    obsfmt=('{n:<3} {id:<30} {ra:<12} {de:<12} {eq:<11} {pmRA:<11} {pmDE:<11} '
    '{irot:<11} {rotmode:<11} {gra1:<11} {gde1:<11} {geq1:<11} '
    '{gra2:<11} {gde2:<11} {geq2:<11}')
    header=obsfmt.format(n='#',
                        id='ID',
                        ra='RA',
                        de='DE',
                        eq='Eq',
                        pmRA='pmRA',
                        pmDE='pmDE',
                        irot='Rot',
                        rotmode='Mode',
                        gra1='GRA1',
                        gde1='GDE1',
                        gra2='GRA2',
                        gde2='GDE2',
                        geq2='GEQ2',
                        geq1='GEQ1')

    file_lines.append(header+'\n')
    obsfmt=('{n:<3} {id:<30} {ra:<12} {de:<12} {eq:<11} {pmRA:<11.2f} '
    '{pmDE:<11.2f} {irot:<11} {rotmode:<11} {gra1:<11} {gde1:<11} {geq1:<11} '
    '{gra2:<11} {gde2:<11} {geq2:<11}')
    for i,r in enumerate(recs):
        s=obsfmt.format(n=i+1,
                        id=r['name'].replace(' ', '_'),
                        ra=sexconvert(r['ra'],ra=True,dtype=str),
                        de=sexconvert(r['de'],dtype=str),
                        eq=r['epoch'],
                        pmRA=r.get('pmra', 0.0),
                        pmDE=r.get('pmdec', 0.0),
                        irot=ROTATOR_SETTING,
                        rotmode='EQU',
                        gra1=sexconvert(0,dtype=str),
                        gde1=sexconvert(0,dtype=str),
                        gra2=sexconvert(0,dtype=str),
                        gde2=sexconvert(0,dtype=str),
                        geq2=0,
                        geq1=0)
        file_lines.append(s+'\n')
    return file_lines


def write_target_list(tfile, recs):
    """rec iterable of dicts with 'name' 'ra', 'de', & 'epoch'"""
    file_lines=generate_tlist_file(recs)
    with open(tfile,'w') as fp:
        for l in file_lines:
            fp.write(l)

if __name__ == '__main__':
    import glob, sys
    fname=sys.argv[1]
    sfile=fname+'_summ.txt'
    tfile=fname+'_tlist.txt'
    files = glob.glob('*.plate')
    
    trec=write_summary_file(sfile, files)
    write_target_list(tfile, trec+extra)
