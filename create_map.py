import numpy as np
import optparse, sys, os

__author__ = 'omniscope'

o = optparse.OptionParser()
o.add_option('-f', '--frequency', action='store', type='float', default=100., help='Frequency of the map in GHz, default 100 GHz.')
o.add_option('-r', '--resolution', action='store', type='float', default=0, help='Required resolution in arcminutes. The output resolution will be either 5 degrees or 0.8 degrees below 10 GHz, either 5 degrees or 0.4 degrees above 10 GHz.')
o.add_option('-u', '--unit', action='store', default='MJy/sr', help='Output unit, default MJy/sr. Other options include TCMB for CMB temperatures in Kelvin or TRJ for Rayleigh-Jeans temperatures in Kelvin.')
o.add_option('-o', '--outputpath', action='store', default=None, help='Path to store the output map.')

opts, args = o.parse_args(sys.argv[1:])
freq = opts.frequency
resolution = opts.resolution
unit = opts.unit
oppath = opts.outputpath

#checking inputs
if oppath is None:
    oppath = os.path.realpath(__file__) + '/gsm2016_%.3eghz_%s_healpynest.txt'%(freq, unit)
elif os.path.isfile(oppath):
    print "PATH ERROR: %s already exists."%oppath
if resolution == 0:
    resolution = 300

print '###Input Parameters###'
print 'Frequency: %.6f GHz'%freq
print 'Resolution: %.3f arcmin'%resolution
print 'Output Unit: ' + unit
print 'Output Path: ' + oppath
sys.stdout.flush()

print '###Reading GSM Data###'
sys.stdout.flush()

print '###Computing New GSM###'
sys.stdout.flush()
result = 0.

print '###Unit Conversion###'
sys.stdout.flush()

print '###Relative Component RMS###'
sys.stdout.flush()


print '###Outputting Result###'
sys.stdout.flush()
np.savetxt(oppath, result)