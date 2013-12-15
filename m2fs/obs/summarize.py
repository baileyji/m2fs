from astropy.io import fits
from flask import Flask, render_template


app = Flask(__name__)
records=[]
keys=[]
default_keys=['FILENAME','OBJECT', 'EXPTIME', 'BINNING', 'SLIDE',
          'HI-AZIM', 'HI-ELEV','LO-ELEV', 'FOCUS', 'FILTER']
keys=[]
@app.route('/')
def summary_table():
    global keys
    return render_template('obs_log.html', column_names=keys, data=records)


def run(files, userkeys=None, prepend=[], append=[], filterfunc=lambda x:True):
    """filterfunc gets a fits header and returns true or false"""
    global keys, records
    if not userkeys:
        keys=default_keys
    else:
        keys=userkeys
    keys=prepend+keys+append
    print "Summarizing {} files.".format(len(files))
    for f in files:
        header=fits.getheader(f)
        if filterfunc(header):
            records.append([header[k] for k in keys])
    #app.debug=True
    app.run()
