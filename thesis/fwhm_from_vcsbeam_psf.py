
import numpy as np
import argparse
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from vcstools.analyse_psf import read_vcsbeam_psf

# np.set_printoptions(threshold=np.inf)

def plot_vcsbeam_psf(psf_file, output_name="vcsbeam_psf.png", normalise=False, vmin=None):
    """Plot the PSF data file output from vcsbeam's mwa_tied_array_beam_psf

    Parameters
    ----------
    psf_file : `str`
        The PSF data file output from vcsbeam's mwa_tied_array_beam_psf.
    output_name : `str`, optional
        Output plot name. |br| Default: vcsbeam_psf.png.
    normalise : `boolean`, optional
        Normalise the PSF. |br| Default: `False`.
    vmin : `float`, optional
        Minimum value of plot. |br| Default: `None`.
    """
    ra, dec, power = read_vcsbeam_psf(psf_file)
    # print(ra)
    # print(dec)
    # print(power)

    if normalise:
        power = power / np.amax(power)

    # fig, ax = plt.subplots()
    # im = ax.pcolormesh(ra, dec, power, norm=colors.LogNorm(), cmap='plasma', vmin=vmin)

    # ax.set_xlabel(r"Right Ascension (hours)")
    # ax.set_ylabel(r"Declination ($^{\circ}$)")
    # fig.colorbar(im,#spacing='uniform', shrink = 0.65, #ticks=[2., 10., 20., 30., 40., 50.],
    #              label="Normalised array factor power")
    # fig.savefig(output_name)

    return ra, dec, power

if __name__ == '__main__':
    # Option parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("--psf", type=str, help="PSF file from VCSBeam's mwa_tied_array_beam_psf", default=None, required=True)
    parser.add_argument("--output", type=str, help="Output basename. Default: vcsbeam_psf", default="vcsbeam_psf")
    args = parser.parse_args()

    ra, dec, power = plot_vcsbeam_psf(args.psf, output_name=f"{args.output}.png", normalise=True, vmin=0.5)
    #convert ra from hours to deg
    ra *= 15
    # print(ra)
    # print(dec)
    # print(power)

    ra_middle, dec_middle = ra.shape
    ra_middle  = ra_middle // 2
    dec_middle = dec_middle // 2
    # print(power[dec_middle, ra_middle])


    # dec line (ra same so it is vertical lines)
    # print(dec[dec_middle,:])
    # print(power[dec_middle,:])
    spline = UnivariateSpline(np.array(dec[dec_middle,:]), np.array(power[dec_middle,:]) - 0.5, s=0)
    enter_beam, exit_beam = spline.roots()
    dec_fwhm = exit_beam - enter_beam
    # print(f"FHWM (dec): {dec_fwhm:.3f} deg")


    # ra line (dec same so it is horizontal lines)
    # print(ra[:,ra_middle])
    # print(power[:,ra_middle])
    spline = UnivariateSpline(np.array(ra[:,ra_middle]), np.array(power[:,ra_middle]) - 0.5, s=0)
    enter_beam, exit_beam = spline.roots()
    ra_fwhm = exit_beam - enter_beam
    # print(f"FHWM (ra):  {ra_fwhm:.3f} deg")

    fwhm = (ra_fwhm + dec_fwhm) / 2
    print(f"FWHM:       {fwhm*60:.1f} arcmin")

    fig, ax = plt.subplots()
    ax.plot(dec[dec_middle,:] - dec[dec_middle, ra_middle], power[dec_middle,:], label="dec")
    ax.plot(ra[:,ra_middle]   - ra[dec_middle, ra_middle],  power[:,ra_middle],  label="ra")
    fig.savefig("FWHM_test.png")

