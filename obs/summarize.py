from astropy.io import fits
from flask import Flask, render_template

app = Flask(__name__)

default_keys=['OBJECT', 'EXPTIME', 'BINNING', 'SLIDE',
          'HI-AZIM', 'HI-ELEV', 'FOCUS', 'FILTER']

@app.route('/')
def summary_table():
    return render_template('obs_log.html', column_names=keys, data=records)


def make_table(files, keys=None):
    if not keys:
        keys=default_keys
    records=[[fits.open(f)[0].header[k] for k in keys] for f in files]
    app.debug=True
    app.run()