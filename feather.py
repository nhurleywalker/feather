#!/usr/bin/env python

from scipy import fftpack
import astropy.io.fits as fits
import numpy as np
#import pylab as py
import os
import sys
# Regridding
import montage_wrapper as montage

import argparse

__author__ = "Natasha Hurley-Walker"
__date__ = "09/03/2018"


# Definitions
sig2fwhm = (2.*np.sqrt(2.*np.log(2.)))

def gaussian2d(x, y, mux, muy, sigmax, sigmay, theta):
   a = np.cos(theta)**2 / (2*sigmax**2) + np.sin(theta)**2 / (2*sigmay**2)
   b = -np.sin(2*theta) / (4*sigmax**2) + np.sin(2*theta) / (4*sigmay**2)
   c = np.sin(theta)**2 / (2*sigmax**2) + np.cos(theta)**2 / (2*sigmay**2)
   g = np.exp(-(a*(x-mux)**2 + 2*b*(x-mux)*(y-muy) + c*(y-muy)**2))
   return g

def fft(data):
    data_fft=fftpack.fftshift(
        fftpack.fft2(
            fftpack.ifftshift(data)
            )
        )
    return data_fft

def ifft(data):
    data_fft=np.real(
        fftpack.fftshift(
            fftpack.ifft2(
                fftpack.ifftshift(data)
                )
            )
        )
    return data_fft


def exportfits(data, header, outfile):
    print "writing to "+outfile
    # Convert back to float32 to make the files smaller
    new = fits.PrimaryHDU(data.astype(np.float32) ,header=header) #create new hdu
    newlist = fits.HDUList([new]) #create new hdulist
    newlist.writeto(outfile, overwrite = True)

def feather(lowres, highres, exportpsf, sdfactor, regrid):
    # Import high-res
    hdu_highres = fits.open(highres)
    # Single -> Double, remove redundant axes
    hr = np.squeeze(np.squeeze(hdu_highres[0].data.astype(float)))
    # Replace NaNs with zeros
    replace = np.isnan(hr)
    hr[replace]=0.0

    try:
        bmaj_highres = hdu_highres[0].header["BMAJ"]
    except KeyError:
        print "No valid beam in header of "+highres
        sys.exit(1)
    bmin_highres = hdu_highres[0].header["BMIN"]
    bpa_highres = hdu_highres[0].header["BPA"]

# Montage copies ALL the fits keys across, including the beam values! So we need to replace those with the original values
# Test whether the file has a beam
    original = fits.open(lowres)
    try:
        bmaj_lowres = original[0].header["BMAJ"]
    except KeyError:
       print "No valid beam in header of "+lowres
       sys.exit(1)

    # Regrid low-res
    if regrid:
        lowres_rg = lowres.replace(".fits", "_montaged.fits")
        if not os.path.exists(lowres_rg):
            montage.mGetHdr(highres,"temp.txt")
            montage.reproject(lowres,lowres_rg,header="temp.txt",exact_size=True)
        else:
            print "Will not overwrite existing regridded image "+lowres_rg
    else:
        print "Not regridding; expecting image co-ordinates to match exactly."
        lowres_rg = lowres
# TODO: add a test for image co-ordinate match rather than letting it get to the FT then fail

    hdu_lowres = fits.open(lowres_rg)
    newhdr = hdu_lowres[0].header
    for fitskey in ["BMAJ", "BMIN", "BPA"]:
        newhdr[fitskey] = original[0].header[fitskey]
#        print fitskey, original[0].header[fitskey]
    try:
        naxis4 = newhdr["NAXIS4"]
    except:
        naxis4 = None
    if naxis4:
        newhdr["NAXIS"] = 4

    hdu_lowres.writeto(lowres_rg, overwrite = True)

    # Import regridded low-res
    hdu_lowres = fits.open(lowres_rg)
    # Single -> Double, remove redundant axes
    lr = np.squeeze(np.squeeze(hdu_lowres[0].data.astype(float)))
    # Replace NaNs with zeros
    replace = np.isnan(lr)
    lr[replace]=0.0

    bmaj_lowres = hdu_lowres[0].header["BMAJ"]
    bmin_lowres = hdu_lowres[0].header["BMIN"]
    bpa_lowres = hdu_lowres[0].header["BPA"]


# TODO: enable diagnostic plots
# if plots
#    py.figure(1)
#    py.clf()
#    py.imshow(np.log10(highres))
#    py.savefig("highres.png")

    hr_fft = fft(hr)
    lr_fft = fft(lr)

    # According to https://casa.nrao.edu/docs/taskref/feather-task.html
    # Scale the low-resolution image by the ratio of the volumes of the two clean beams
    # (high-res / low-res)

    ratio = (sdfactor*bmaj_highres*bmin_highres) / (bmaj_lowres*bmin_lowres)
    print ratio

    #Add to this, the uv-grid of the high-resolution image, scaled by  
    #                   (1-wt) where "wt" is the Fourier transform of the "clean beam"
    #                   defined in the low-resolution image.  

    # Make a model image of low-resolution psf
    xmax = hdu_lowres[0].header["NAXIS1"]
    ymax = hdu_lowres[0].header["NAXIS2"]
    try:
        pix2deg = hdu_lowres[0].header["CDELT2"]
    except KeyError:
        pix2deg = hdu_lowres[0].header["CD2_2"]

    x, y = np.meshgrid(np.linspace(0,xmax,xmax), np.linspace(0,ymax,ymax))

    sigmax = bmaj_lowres / (pix2deg * sig2fwhm)
    sigmay = bmin_lowres / (pix2deg * sig2fwhm)
    mux = xmax/2 + 0.5
    muy = ymax/2 + 0.5

    g = gaussian2d(x, y, mux, muy, sigmax, sigmay, np.deg2rad(bpa_lowres))
    g_fft = fft(g)

    if exportpsf:
        exportfits(g, hdu_highres[0].header, lowres.replace(".fits","_psf.fits"))
        exportfits(np.real(g_fft), hdu_highres[0].header, lowres.replace(".fits","_psf_fft_real.fits"))
        try:
            exportfits(np.imag(g_fft), hdu_highres[0].header, lowres.replace(".fits","_psf_fft_imag.fits"))
        except:
# I get some weird errors when I try to export an incredibly tiny imaginary part, but this works:
            exportfits(np.zeros(g_fft.shape), hdu_highres[0].header, lowres.replace(".fits","_psf_fft_imag.fits"))

    # Add together
    comb_fft = ratio * lr_fft + (1 - (g_fft/np.nanmax(g_fft))) * hr_fft

    # Inverse FFT
    comb = ifft(comb_fft)
    exportfits(comb, hdu_highres[0].header, highres.replace(".fits","+")+lowres)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("files to read")
    group1.add_argument('--lowres', dest='lowres', default=None, \
                        help="Low-resolution file to read")
    group1.add_argument('--highres', dest='highres', default=None, \
                        help="High-resolution file to read")
    group2 = parser.add_argument_group("functions to perform")
    group2.add_argument('--regrid', dest='regrid', default=False, action='store_true', \
                        help="Perform regridding with montage to match the low-resolution image to the high-resolution image; \
                              if this is not selected then the image co-ordinates must already match exactly. (default = False)")
    group2.add_argument('--psf', dest='psf', default=False, action='store_true', \
                        help="Save the intermediate PSF FFT information (default = False)")
    group2.add_argument('--plot', dest='makeplots', default=False, action='store_true', \
                        help="Make png plots of the FITS images and FFTs (Currently disabled)")
    group3 = parser.add_argument_group("fiddle factors")
    group3.add_argument('--sdfactor', dest='sdfactor', default=1.0, type=float, \
                        help="Float by which to multiply the low-resolution data (default = 1.0)")

    options = parser.parse_args()

    if len(sys.argv) <= 2:
        parser.print_help()
        sys.exit()

    # Perform the feathering
    feather(options.lowres, options.highres, options.psf, options.sdfactor, options.regrid)
