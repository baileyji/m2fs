"""This scratchpad assumes the use of m2fs_reduce through stacking, any dark subtraction or
additional bias subtraction already done

What remains is to
- identify spectral traces for extraction: mark_orders()
- associate traces with fibers
- scattered light: estimate_background_scatter()- scattered light: estimate_background_scatter()
- run the REDUCE slit curvature finding stuff, sorting out which files this is run on its TBD

- flat fielding
- wavelength calibration
- science spectrum extraction
- continuum normilization

"""


# this is through CR and stacking, now to trace≤ wavecal, extract etc with pyreduce

from m2fs.seqno import *
from pyreduce.reduce import *
import pyreduce
from pyreduce.instruments.mspec_ifum import MSPEC_IFUM

list_file = '/Users/one/no_backup/lcoreduction/ifum/list_sept'
raw_dir = '/Users/one/no_backup/lcoreduction/ifumraw/'
out_dir = '/Users/one/no_backup/lcoreduction/ifum/ifum_sept'

dry_run = False  # args.dry_run
overwrite = False
force_gzip = True
do_cosmic_stack = False

dset = M2FSIFUMListfileDataset(list_file, raw_dir, out_dir)
msf = dset.msf
seq_strs = dset.seq_strs



instrument = MSPEC_IFUM()
configuration = pyreduce.configuration.get_configuration_for_instrument(instrument, plot=1)

# Ansgar's steps args: instrument, mode, target, night, output_dir, order_range, **config


# orders
config = configuration['orders']

order_img = dset['ut20220901/r1449'].r[0]
order_img = fits.open(order_img.merged_file())  # a twilight or LED flat
ohead = order_img['science'].header

# pyreduce was doing this for all the images under the hood
order_img = np.ma.masked_array(order_img['science'].data, mask=order_img['mask'].data)

mapping = dict(opower='degree', sigma='split_sigma')
config.update({k: config[v] for k, v in mapping.items()})

# Ansgar's code bounces back and forth between orders, column_range and orders=orders, column_range
orders = mark_orders(order_img, **config)  # orders=orders, column_range

# scatter
config = configuration['scatter']
scatter_img = msf.fits  # an image we want to correct for scattered light (likely a science or sky)
shead = scatter_img['bias_corrected'].header
scatter_img = np.ma.masked_array(scatter_img['bias_corrected'].data, mask=scatter_img['mask'].data)
scatter = estimate_background_scatter(scatter_img, orders[0], column_range=orders[1], **config)

# curvature
config = configuration['curvature']
order_range = None  # this isn't a config parameter in pyreduce but rather a user param

orig = msf.fits  # an arc image we want to use to estimate the slit curvature
head = orig['bias_corrected'].header
orig = np.ma.masked_array(orig['bias_corrected'].data, mask=orig['mask'].data)

extraction_kwargs = dict(extraction_width=config["extraction_width"], sigma_cutoff=config["extraction_cutoff"],
                         collapse_function=config["collapse_function"],
                         plot=config['plot'], plot_title=config['plot_title'], order_range=order_range,
                         extraction_type=config['extraction_method'])

extracted, _, _, _ = extract(orig, orders[0],
                             gain=head["e_gain"], readnoise=head["e_readn"], dark=head["e_drk"],
                             column_range=orders[1], **extraction_kwargs)

mapping = dict(fit_degree='degree', curvature_mode='dimensionality', sigma_cutoff='curvature_cutoff',
               extraction_type='extraction_method')
config.update({k: config[v] for k, v in mapping.items()})
module = CurvatureModule(orders[0], column_range=orders[1], **config)
curvature = module.execute(extracted, orig)  # curvature = tilt, shear

# normalize flat field
config = configuration['norm_flat']
flat = msf.fits  # an image we want to use to estimate the blaze with, will twilight work, quartz, probably not led
fhead = orig['bias_corrected'].header
flat = np.ma.masked_array(flat['bias_corrected'].data, mask=flat['mask'].data)

order_range = None
config['order_range'] = order_range
mapping = dict(lambda_sf="smooth_slitfunction", lambda_sp="smooth_spectrum", osample="oversampling",
               # sigma_cutoff="extraction_cutoff"
               extraction_type='extraction_method')
config.update({k: config[v] for k, v in mapping.items()})

# if threshold is smaller than 1, assume percentage value is given
threshold = config['threshold']
threshold = np.percentile(flat, threshold * 100) if threshold <= 1 else threshold

norm, _, blaze, _ = extract(
    flat,
    orders[0],
    gain=fhead["e_gain"],
    readnoise=fhead["e_readn"],
    dark=fhead["e_drk"],
    column_range=orders[1],
    scatter=scatter,
    threshold=threshold,
    tilt=curvature[0],
    shear=curvature[1],
    **extraction_kwargs,
)

blaze = np.ma.filled(blaze, 0)
norm = np.nan_to_num(norm, nan=1)
norm_flat = norm, blaze

# wavecal_master
config = configuration['wavecal_master']
order_range = None  # this isn't a config parameter in pyreduce but rather a user param

orig = msf.fits  # an arc image we want to use for the wavelength solution
thead = orig['bias_corrected'].header
orig = np.ma.masked_array(orig['bias_corrected'].data, mask=orig['mask'].data)
orig /= norm_flat[0]

extraction_kwargs = dict(extraction_width=config["extraction_width"], sigma_cutoff=config["extraction_cutoff"],
                         collapse_function=config["collapse_function"],
                         plot=config['plot'], plot_title=config['plot_title'], order_range=order_range,
                         extraction_type=config['extraction_method'])

thar, _, _, _ = extract(orig, orders[0],
                        gain=head["e_gain"], readnoise=head["e_readn"], dark=head["e_drk"],
                        column_range=orders[1], tilt=curvature[0], shear=curvature[1], **extraction_kwargs)
wavecal_master = thar, thead

# wavecal_init
config = configuration['wavecal_init']
thar, thead = wavecal_master
# Get the initial wavelength guess from the instrument
wave_range = instrument.get_wavelength_range(thead, None)

module = WavelengthCalibrationInitializeModule(*config)
linelist = module.execute(thar, wave_range)
wavecal_init = linelist

# wavecal
config = configuration['wavecal']

thar, thead = wavecal_master
linelist = wavecal_init

module = WavelengthCalibrationModule(**config)
wave, coef = module.execute(thar, linelist)

# science
config = configuration['science']

img = msf.fits  # an image we want to extract
head = img['bias_corrected'].header
img = np.ma.masked_array(img['bias_corrected'].data, mask=img['mask'].data)
img /= norm_flat[0]

mapping = dict(lambda_sf="smooth_slitfunction", lambda_sp="smooth_spectrum", osample="oversampling",
               sigma_cutoff="extraction_cutoff", extraction_type='extraction_method')
config.update({k: config[v] for k, v in mapping.items()})
config['order_range'] = order_range

data, unc, blaze, cr = extract(img, orders[0], column_range=orders[1],
                               gain=head["e_gain"],
                               readnoise=head["e_readn"],
                               dark=head["e_drk"],
                               tilt=curvature[0], shear=curvature[1],
                               **config)

spec, sigma, _, cr = data, unc, blaze, cr

# continuum normalization
config = configuration['continuum']
#TBD
