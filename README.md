# HETobs
This package is intended to suppliment [Panacea](https://github.com/grzeimann/Panacea).  Currently, the package has tools for aquisition camera (ACAM) reductions which include overscan subtraction, trimming, bias subtraction, sky background modeling and removal, source detection, and astrometry.  The last task is accomplished using an archival catalog match to sources extracted from the image starting from a relatively good guess.  The stability and robustness of the algorithm is not yet tested.

## Getting Started
```
cd WHERE_YOU_WANT_FOLDER
git clone https://github.com/grzeimann/HETobs.git
cd HETobs
```
Since hetobs.py, which is the functioning program, is in test mode, to utilize this package merely edit lines 35-37:
```
rootdir = '/Users/gregz/cure/lrs2_raw/20170721'
observation = ('/Users/gregz/cure/lrs2_raw/20170721/lrs2/lrs20000004'
               '/exp01/lrs2/20170721T045123.5_066LL_sci.fits')
```
The rootdir variable should point to a date folder for one night's worth of HET data, and the observation should point to an example frame on that night for which you want to test the program.  Then simply run it:
```
python hetobs.py
```          
The output is ...

### Prerequisites

The prerequisites are almost all bundled within [anaconda](https://anaconda.org/anaconda/python).  Specifically:

```
numpy
matplotlib
astropy
photutils
httplib
urllib
```

The program doesn't require pyhetdex yet, but it may soon.  Follow the link for instructions if you want to install pyhetdex: http://www.mpe.mpg.de/~montefra/documentation/pyhetdex/0.12.0/


## Authors

* **Greg Zeimann**

## License

No details yet.

