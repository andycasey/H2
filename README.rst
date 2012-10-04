This repository contains much of the data and code that was used for
the paper by Crighton et al., "A high molecular fraction in a
sub-damped absorber at z=0.56". There's no need to have git installed
to download the repository, click on 'Downloads' at the top right of
this page and you can download the whole repository as a .zip or
.tar.gz file.

Below is a short explanation of the files included.

If you need any help understanding what was done for the paper, or
getting the code to run, please contact me (neilcrighton .at. gmail
.dot. com). If you find any of the data or code in the repository
useful, please let me know by email. Finally, if you use any of the
code or data for a publication, please cite the paper
(http://arxiv.org/abs/1210.0905).


**data/**

This directory contains data and catalogues used in the paper.

- C.txt
 
  Combined COS spectrum, with columns wavelength (Angstroms), flux,
  error (1 sigma), continuum. The raw exposures and calibration
  files are available from the HST archive.
 
- q0107c_HIRES.txt
 
  Combined HIRES spectrum, with columns wavelength (Angstroms),
  flux, error (1 sigma), continuum. The raw exposures will be
  available in the Keck archive
  (http://www2.keck.hawaii.edu/koa/public/koa.php) in late 2012/
  early 2013. Email me if you need them before this.
 
- qsoc_G190H.txt
 
  Combined FOS spectrum, with columns wavelength (Angstroms), Flux,
  error (1 sigma), continuum. The raw exposures and calibration
  files are available from the HST archive.

- K.fits.fz

  K-band image of Q 0107-0232 compressed with fpack
  (http://heasarc.nasa.gov/fitsio/fpack/). The individual raw
  exposures are available in the ESO archive.

- f_H2_DLA.tbl
 
  Table in the excellent IPAC format
  (http://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html)
  of molecular fraction measurements in high redshift DLAs. The
  references from which these are taken are given in the table header.

- Richter03_t2.tbl

  Table 2 from Richter et al. 2003, ApJ, 586, 230.

- Savage77_tab1.tbl

  Table 1 from Savage et al. 1977, ApJ, 216, 291.

- Welty12_LMC.tbl and Welty12_SMC.tbl

  The LMC and SMC measurements from Tables 1 and 2 in Welty et
  al. 2012, ApJ, 745, 173

- observed_logN/
  
  Text files with the column densities of all the transitions
  measured for the sub-DLA, for all components summed and for each
  component individually. See the README in this directory.

**cloudy/**

Example cloudy input files for each of the scenarios described in
Section 4.5 of the paper.

**scripts/**

Python scripts used to perform many of the calculations in the
paper.

**modules/**

This contains most of the modules needed by the various python
scripts. You will also need the following packages installed:

- Numpy and Scipy,    http://numpy.scipy.org/
- Matplotlib,         http://matplotlib.org/
- Pyfits,             http://www.stsci.edu/institute/software_hardware/pyfits
- Atpy,               http://atpy.github.com/

The Enthought Python distribution
(http://www.enthought.com/products/epd_free.php) is a quick way to
install most of these.


**mcmc/**

The samples from the monte carlo markov chains generated with Emcee,
and the python module used to define the likelihood.
