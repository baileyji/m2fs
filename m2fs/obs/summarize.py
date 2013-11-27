from astropy.io import fits
from flask import Flask, render_template

app = Flask(__name__)

default_keys=['FILENAME','OBJECT', 'EXPTIME', 'BINNING', 'SLIDE',
          'HI-AZIM', 'HI-ELEV', 'FOCUS', 'FILTER']
keys=[]
@app.route('/')
def summary_table():
    global keys
    return render_template('obs_log.html', column_names=keys, data=records)


def make_table(files, userkeys=None):
    global keys, records
    if not userkeys:
        keys=default_keys
    else:
        keys=userkeys
    records=[[fits.open(f)[0].header[k] for k in keys] for f in files]
    app.debug=True
    app.run()