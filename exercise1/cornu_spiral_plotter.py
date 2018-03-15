#!/usr/local/bin/python3
# c-basic-offset: 4; tab-width: 4; indent-tabs-mode: nil
# vi: set shiftwidth=4 tabstop=4 expandtab
# :indentSize=4:tabSize=4:noTabs=true

from multiprocessing.pool import Pool
import scipy.integrate as integrate
from math import atan2
from numpy import pi, cos, sin, sqrt, abs, sign, linspace
from os import mkdir
from os.path import exists
import matplotlib.pyplot as plt


def s_integrand(arg: float):
    """Shorthand function."""
    return sin((pi*arg**2)/2)


def c_integrand(arg: float):
    """Shorthand function."""
    return cos((pi*arg**2)/2)


def integral(f: callable, upper_limit: float,
             lower_limit: float = 0) -> tuple():
    """Wrapper function. Pre-applies the integral to shorten the code"""
    return integrate.quad(lambda arg: f(arg), lower_limit, upper_limit)


def fresnel_c(u: float, lower_limit: float = 0.0) -> tuple():
    """Shorthand function."""
    return integral(c_integrand, u, lower_limit=lower_limit)


def fresnel_s(u: float, lower_limit: float = 0.0) -> tuple():
    """Shorthand function. """
    return integral(s_integrand, u, lower_limit=lower_limit)


def map_to_array_mp(f: callable, args: iter) -> iter:
    """A multiprocess functional mapping a function to an iterable object

    Parameters
    ----------
    f: function
    args: iterable

    Returns
    -------
    list or list of tuples, depending on the output of function.
    """
    with Pool() as p:
        xs = p.map(f, args)
    if isinstance(xs[0], tuple):
        x1, x2 = [x[0] for x in xs], [x[1] for x in xs]
        return x1, x2
    else:
        return xs


def plot_cornu_spiral(rng: float = 20, rate: int = 2**12) -> None:
    """Produce a Cornu spiral with given range and resolution

    Parameters
    ----------
    rng: float64
    Range of the integration variable to cover
    rate: int
    Number of sampling points on the given range.
    """
    samples = linspace(-rng, rng, rate)
    plt.figure()
    y, _ = map_to_array_mp(fresnel_s, samples)
    x, _ = map_to_array_mp(fresnel_c, samples)
    plt.plot(x, y)
    plt.xlabel('$C(u)$')
    plt.ylabel('$S(u)$')
    plt.title("The Cornu spiral\n "
              "a geometric representation of the Fresnel integrals")
    ticks_u = list(sign(x)*sqrt(abs(x)) for x in range(-5, 5))
    ticks_y, _ = map_to_array_mp(fresnel_s, ticks_u)
    ticks_x, _ = map_to_array_mp(fresnel_c, ticks_u)
    plt.plot(ticks_x, ticks_y, 'k+')
    if not exists('figures'):
        mkdir('figures')
    plt.savefig('figures/Cornu-spiral-plot.pdf')
    plt.show()


class DiffractionPattern:
    """A class to encapsulate the notion of a (Fresnel) diffraction pattern"""

    def __init__(self, screen_distance: float, wavelength: float = 1.0,
                 slit_width: float = 10.0) -> None:
        self.slit_width = slit_width
        self.wavelength = wavelength
        self.screen_distance = screen_distance
        self.scale = sqrt(2/(self.wavelength*self.screen_distance))

    def amplitude_at(self, x: float) -> float:
        """Amplitude on screen at position x"""
        re = self.real(x)
        im = self.imaginary(x)
        return sqrt(re**2 + im**2)

    def imaginary(self, x: float) -> float:
        """Imaginary part of the Fresnel integral at x"""
        high = (x + self.slit_width/2)*self.scale
        low = (x - self.slit_width/2)*self.scale
        im = fresnel_s(high, low)[0]
        return im*self.scale

    def real(self, x: float) -> float:
        """Real part of the Fresnel integral at x"""
        high = x + self.slit_width/2
        low = x - self.slit_width/2
        re = fresnel_c(high, low)[0]
        return re*self.scale

    def phase_at(self, x: float) -> float:
        """Complex phase of the diffraction pattern at position x"""
        re = self.real(x)
        im = self.imaginary(x)
        return atan2(im, re)


def plot_pattern(distance: float, axarr) -> None:
    """Plot a Fresnel diffraction pattern given distance to aperture."""
    data = linspace(-20, 20, 5000)
    p = DiffractionPattern(distance)
    amplitude = map_to_array_mp(p.amplitude_at, data)
    phase = map_to_array_mp(p.phase_at, data)
    axarr[0].plot(data, amplitude, label=distance)
    axarr[1].plot(data, phase, label=distance)


def plot_diffraction_patterns(distances: iter = None) -> None:
    """Plot Fresnel diffraction patterns, at different distances from
    the aperture.

    Parameters
    ----------
    distances: list
    list of distances for which to make plots.
    """
    if distances is None:
        distances = [30, 50, 100]
    fig, axes = plt.subplots(2, sharex=True)
    for d in distances:
        plot_pattern(d, axes)
    plt.legend(loc='best')
    plt.xlabel('Position on screen / cm')
    axes[0].set_ylabel(r'Intensity / W $m^{-2}$')
    axes[1].set_ylabel(r'Phase / radians')
    axes[0].set_title('Distribution of intensity in the diffraction pattern(s)')
    axes[1].set_title('Phase as a function of position')
    plt.savefig('figures/Fresnel-intensity.pdf')
    plt.show()


if __name__ == '__main__':
    plot_cornu_spiral()
    plot_diffraction_patterns()
