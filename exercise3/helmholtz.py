#!/usr/local/bin/python3
# c-basic-offset: 4; tab-width: 4; indent-tabs-mode: nil
# vi: set shiftwidth=4 tabstop=4 expandtab
# :indentSize=4:tabSize=4:noTabs=true

from functools import partial
from genericpath import exists
from multiprocessing.pool import Pool
import matplotlib.pyplot as plt
from numpy import array, cross, dot, sqrt, pi, cos, sin, linspace, zeros, \
    shape, meshgrid, vstack
from os import mkdir


def save_figure(file_name):
    """Routine to quickly save figures"""
    if not exists('figures/'):
        mkdir('figures/')
    plt.savefig('figures/%s.pdf' % file_name)


class Wire:
    """Abstract wire class """

    def __init__(self, segments=None):
        self.segments = segments


class Segment:
    """Segment. Usually something for which the field can be directly
    evaluated. """

    def __init__(self, r_0: array, dl: array, current):
        self.r_0 = r_0
        self.dl = dl
        self.current = current


class StraightWire(Wire):
    """CLass for Straight wire. """

    def __init__(self, current, dl, r_0):
        super().__init__([Segment(current=current, dl=dl, r_0=r_0)])


class CircularWire(Wire):
    """Class for Circular wire. Represents segments using straight wire
    objects. """

    def __init__(self, current=1, radius=1,
                 centre_location: array = array([0, 0, 0]),
                 resolution: int = 2**6):
        self.current = current
        self.radius = radius
        self.centre_location = centre_location
        self.resolution = resolution
        super().__init__(list(self.populate_segments(resolution)))

    def populate_segments(self, resolution):
        """
        Splits the circle into 2*resolution StraightWire objects of equal
        length positioned on the vertices of a regular 2*resolution-gon.

        Since we want to have imposable symmetries, the segments are
        produced not in succession but rather with the opposite first.
        This can cause underflow problems in rare cases.
        """
        theta = pi/resolution
        normal = array([0, self.radius, 0])
        chord = array([-2*self.radius*sin(theta/2), 0, 0], )

        cf = cos(theta)
        sf = sin(theta)
        rotation_matrix = array([[cf, -sf, 0], [sf, cf, 0], [0, 0, 1]], )
        for i in range(0, resolution):
            yield StraightWire(current=self.current, dl=chord, r_0=normal)
            yield StraightWire(current=self.current, dl=-chord,
                               r_0=-normal)
            chord = dot(rotation_matrix, chord)
            normal = dot(rotation_matrix, normal)


def biot_savart(at, wire: StraightWire):
    """Evaluate field of straight wire segment"""
    seg = wire.segments[0]
    dr = at - seg.r_0
    mod_r = sqrt(dot(dr, dr))
    if mod_r == 0:
        return array([0, 0, 0])
    return seg.current*cross(seg.dl, dr)/(4*pi*mod_r**3)


def field(wire, position):
    """Evaluate the field of a wire position position(s)."""
    if shape(position) == (3,):
        f = partial(biot_savart, position)
        s = sum(map(f, wire.segments))
        return s
    else:
        fld = partial(field, wire)
        with Pool() as p:
            ans = p.map(fld, position)
        return array(ans)


def superimpose(wires, at):
    """Evaluate the field of Wires at given position(s)"""
    f = partial(lambda x, y: field(y, x), at)
    return sum(map(f, wires))


def gen_z_spaced(low=0., high=2., number_of_samples=50):
    """Helper function to generate """
    m = zeros(shape=(number_of_samples, 3))
    number_of_samples = linspace(low, high, number_of_samples)
    m[:, 2] += number_of_samples
    return m


def single_coil_on_axis(resolutions=None):
    """Investigate the agreement of theoretical prediction of on-axis
    field with numerical result"""
    if resolutions is None:
        resolutions = [2**r for r in [4, 6, 7, 10]]
    ts = [CircularWire(resolution=r) for r in resolutions]
    for t in ts:
        m = gen_z_spaced(number_of_samples=50)
        data = field(t, m)
        plt.xlabel('z / m')
        plt.ylabel(r'$B_z$ / T')
        plt.title(
            'On-axis value of magnetic field.\n Comparison of theory '
            'with simulated data')
        plt.plot(m[:, 2], data[:, 2], 'k+', label='data', )
        theory = 1/(2*(1 + m[:, 2]**2)**(3/2))
        plt.plot(m[:, 2], theory, 'r-', label='theory')
        plt.legend(loc='best')
        save_figure('single_coil_on_axis-values-' + str(t.resolution))
        plt.show()
        plt.title('On axis value of magnetic field.\n Residuals in '
                  'proportion to theoretical value')
        plt.xlabel('z / m')
        plt.ylabel(r'$\frac{\Delta B_z} {B_z}$ / arb. units')
        residuals = (theory - data[:, 2])/theory
        plt.plot(m[:, 2], residuals)
        save_figure('single_coil_on_axis-residues' + str(t.resolution))
        plt.show()


def generate_yz_space(n=20, low: float = -2, high: float = 2) -> array:
    """Helper function. Generates thw section through the y-z plane"""
    x = zeros(1)
    y = linspace(low, high, n)
    z = linspace(low, high, n)
    grid_x, grid_y, grid_z = meshgrid(x, y, z)
    grid_x, grid_y, grid_z = grid_x.reshape(n**2, ), grid_y.reshape(
        n**2, ), grid_z.reshape(n**2, )
    return vstack((grid_x, grid_y, grid_z)).T


def longest(vectors: array) -> array:
    """Finds the longest vector in an array"""
    moduli = [sqrt(dot(v, v)) for v in vectors]
    return max(moduli)


def yz_coil(m, coils=None, low: float = -2., high: float = 2.,
            reference_point: array = None) -> float:
    """Plots the magnetic field in the y-z plane, for a given set of
    coils"""
    if coils is None:
        coils = [CircularWire()]
    fig = plt.figure()
    ax = fig.gca()
    ax.set_xlabel('y / m')
    ax.set_ylabel('z / m')
    plt.title('Plot of the magnetic field.')
    args = generate_yz_space(n=m, low=low, high=high)
    res = superimpose(coils, args)
    result_grid_y, result_grid_z = res[:, 1:3].T.reshape(2, m, m)
    arg_grid_y, arg_grid_z = args[:, 1:3].T.reshape(2, m, m)
    ax.quiver(arg_grid_y, arg_grid_z, result_grid_y, result_grid_z)
    save_figure(str(len(coils)) + '_coils_yz_section')
    plt.show()
    if reference_point is not None:
        deltas = res - superimpose(coils, reference_point)
        return longest(deltas)


def helmholtz_coils(d: float = 1/2) -> None:
    """Investigate the behaviour of helmholtz coils"""
    wires = [CircularWire(centre_location=array([0, 0, z])) for z in
             [-d, d]]
    samples = gen_z_spaced(low=-1, high=1, number_of_samples=40)
    data = superimpose(wires, samples)[:, 2]  # z component
    plt.plot(samples[:, 2], data)
    act = (4/5)**(3/2)
    # print(act - data[10]) # -1.2 e-13
    plt.title('Value of field for Helmholtz coils')
    plt.ylabel(r'$B_z$ / T')
    plt.xlabel('z / m')
    plt.plot([-0.5, 0.5], [act, act], 'k+')
    save_figure('helmholtz_coils_on_axis')
    plt.show()
    print(yz_coil(25, wires, -.05, .05, reference_point=array([0, 0, 0])))
    # 0.00419 T


def many_coils_on_axis(number: int = 3, d: float = 5) -> None:
    """Generate a plot of field on axis of a many-coil system"""
    wires = [CircularWire(centre_location=array([0, 0, z])) for z in
             linspace(-d, d, number)]
    samples = gen_z_spaced(low=-1, high=1, number_of_samples=25)
    data = superimpose(wires, samples)[:, 2]  # z component
    plt.plot(samples[:, 2], data, label=number)


def investigate_many_coils(numbers: list = None) -> None:
    if numbers is None:
        numbers = [3, 5, 12]
    for n in numbers:
        many_coils_on_axis(n)
    plt.legend(loc='best')
    plt.title('Field generated by many coils')
    plt.ylabel(r'$B_z$ / T')
    plt.xlabel('z / m')
    save_figure('multiple_coils_on_axis')
    plt.show()


if __name__ == '__main__':
    single_coil_on_axis()
    yz_coil(25)
    helmholtz_coils()
    investigate_many_coils()
