__author__ = 'omniscope'
import numpy as np
import optparse, sys, os

script_path = os.path.dirname(os.path.realpath(__file__))
labels = ['Synchrotron', 'CMB', 'HI', 'Dust1', 'Dust2', 'Free-Free']
n_comp = len(labels)
kB = 1.38065e-23
C = 2.99792e8
h = 6.62607e-34
T = 2.725
hoverk = h / kB

def K_CMB2MJysr(K_CMB, nu):#in Kelvin and Hz
    B_nu = 2 * (h * nu)* (nu / C)**2 / (np.exp(hoverk * nu / T) - 1)
    conversion_factor = (B_nu * C / nu / T)**2 / 2 * np.exp(hoverk * nu / T) / kB
    return  K_CMB * conversion_factor * 1e20#1e-26 for Jy and 1e6 for MJy

def K_RJ2MJysr(K_RJ, nu):#in Kelvin and Hz
    conversion_factor = 2 * (nu / C)**2 * kB
    return  K_RJ * conversion_factor * 1e20#1e-26 for Jy and 1e6 for MJy


o = optparse.OptionParser()
o.add_option('-f', '--frequency', action='store', type='float', default=100., help='Frequency of the map in GHz, default 100 GHz.')
o.add_option('-r', '--resolution', action='store', type='float', default=0, help='Required resolution in arcminutes. The output resolution will be either 5 degrees or 0.8 degrees below 10 GHz, either 5 degrees or 0.4 degrees above 10 GHz.')
o.add_option('-u', '--unit', action='store', default='MJysr', help='Output unit, default MJysr. Other options include TCMB for CMB temperatures in Kelvin or TRJ for Rayleigh-Jeans temperatures in Kelvin.')
o.add_option('-o', '--outputpath', action='store', default=None, help='Path to store the output map, including file name.')
o.add_option('--ring', action='store_true', default=False, help='Output HEALPIX RING format. Default is NEST.')

opts, args = o.parse_args(sys.argv[1:])
freq = opts.frequency
resolution = opts.resolution
unit = opts.unit
oppath = opts.outputpath
convert_ring = opts.ring

#checking inputs
if oppath is None:
    oppath = script_path + '/output/gsm2016_%.3eghz_%s_healpy%s.txt'%(freq, unit, ['NEST', 'RING'][int(convert_ring)])
elif os.path.isfile(oppath):
    print "PATH ERROR: %s already exists."%oppath
if unit not in ['MJysr', 'TCMB', 'TRJ']:
    print "UNIT ERROR: %s not supported. Only MJysr, TCMB, TRJ are allowed."%unit
if resolution == 0:
    resolution = 300

if resolution < 300:
    nside = 1024
    if freq < 10:
        op_resolution = 56
    else:
        op_resolution = 24
else:
    nside = 64
    op_resolution = 300


print '###Input Parameters###'
print 'Frequency: %.6f GHz.'%freq
print 'Requested Resolution: %.3f arcmin.'%resolution
print 'Unit: %s.'%unit
print 'Output Resolution: %.3f arcmin.'%op_resolution
print 'Output HEALPIX Format: nside = %i, NEST.'%nside
print 'Output Path: ' + oppath
print '######################'
sys.stdout.flush()

print 'Reading GSM Data...',
sys.stdout.flush()
if resolution < 300:
    map_ni = np.array([np.fromfile(script_path + '/data/highres_%s_map.bin'%lb, dtype='float32') for lb in labels])
else:
    map_ni = np.loadtxt(script_path + '/data/lowres_maps.txt')
spec_nf = np.loadtxt(script_path + '/data/spectra.txt')
nfreq = spec_nf.shape[1]

print 'Computing New GSM...',
sys.stdout.flush()

left_index = -1
for i in range(nfreq - 1):
    if freq >= spec_nf[0, i] and freq <= spec_nf[0, i + 1]:
        left_index = i
        break
if left_index < 0:
    print "FREQUENCY ERROR: %.2e GHz is outside supported frequency range of %.2e GHz to %.2e GHz."%(freq, spec_nf[0, 0], spec_nf[0, -1])

interp_spec_nf = np.copy(spec_nf)
interp_spec_nf[0:2] = np.log10(interp_spec_nf[0:2])
x1 = interp_spec_nf[0, left_index]
x2 = interp_spec_nf[0, left_index + 1]
y1 = interp_spec_nf[1:, left_index]
y2 = interp_spec_nf[1:, left_index + 1]
x = np.log10(freq)
interpolated_vals = (x * (y2 - y1) + x2 * y1 - x1 * y2) / (x2 - x1)
result = np.sum(10.**interpolated_vals[0] * (interpolated_vals[1:, None] * map_ni), axis=0)
print 'Done.'

print 'Unit Conversion at %.2f GHz:'%freq,
sys.stdout.flush()
if unit == 'TCMB':
    conversion = 1. / K_CMB2MJysr(1., 1e9 * freq)
elif unit == 'TRJ':
    conversion = 1. / K_RJ2MJysr(1., 1e9 * freq)
else:
    conversion = 1.
result *= conversion
print '1 MJysr == %.2e %s'%(conversion, unit)
sys.stdout.flush()

if convert_ring:
    print 'Converting to HEALPIX RING...',
    sys.stdout.flush()
    try:
        import healpy as hp
    except:
        print "Healpy package not found. Cannot convert to HEALPIX RING format."
    result = hp.reorder(result, n2r=True)
    print 'done.'
    sys.stdout.flush()

print 'Outputting Result to %s...'%oppath,
sys.stdout.flush()
if oppath[-4:] == '.npz':
    np.savez(oppath, result.astype('float32'))
elif oppath[-4:] == '.bin':
    result.astype('float32').tofile(oppath)
else:
    print '(Saving as a text file may be very slow. For faster speed, try .bin for single precision binary or .npz for compressed numpy format.)...',
    sys.stdout.flush()
    np.savetxt(oppath, result, fmt='%.3e')

print 'All done.'
sys.stdout.flush()
