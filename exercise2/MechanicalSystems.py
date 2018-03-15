#!/usr/local/bin/python3
# c-basic-offset: 4; tab-width: 4; indent-tabs-mode: nil
# vi: set shiftwidth=4 tabstop=4 expandtab
# :indentSize=4:tabSize=4:noTabs=true

from numpy import linspace, sin, cos, pi, zeros, fmod
from scipy.integrate import odeint
from matplotlib import pyplot as plt


class MechanicalSystem:
    """
    General mechanical 1D system.

    A class definition of a general dynamical system, whose motion is
    defined by a single second order ODE.

    It is assumed that the ODe has form
    $$
    y'' = a_1(y') + a_0(y) + b(t)
    $$

    Attributes
    ----------
    a1 : Function
    Coefficient of first derivative (damping term)

    a0 : function
    Coefficient of second derivative (restoring force)

    b :  function
    Coefficient of free term (driving force)

    ans :  function, optional
    Analytical solution if one is known

    """

    y_name = 'y'
    ydot_name = r'$\dot{' + y_name + '} $'
    total_energy_stored = None
    initial_energy = 0
    y_unit = 'arb. units'

    def __init__(self, a1: callable, a0: callable, b: callable, y0: iter,
                 ans: callable = None):
        self.a1 = a1
        self.a0 = a0
        self.b = b
        self.y0 = y0
        self.ans = ans
        self.rate = 500
        self.samples = None
        self.simulation_data = None

    def to_coupled_linear(self, y: tuple, t: float) -> list:
        r"""
        Parse second order ODE to two first order coupled.

        One of those ODE's is just the definition of the new variable
        as the derivative.

        Parameters
        ----------
        y : tuple
        has form (\theta, \dot{\theta})

        t : numerical value
        placeholder variable for re-evaluating the functions

        Returns
        -------
        out : list
        [y[1], y[1] in terms of ODE]

        """
        return [y[1], self.a1(y[1]) + self.a0(y[0]) + self.b(t)]

    def simulate(self, rate=500, duration=10.0):
        """
        Perform simulation with given parameters.

        Parameters
        ----------
        rate : int, optional
        How many samples per unit time to simulate

        duration : int, optional
        How long to run simulation for in seconds

        """
        self.rate = rate
        self.samples = linspace(0, duration, rate*duration)
        self.simulation_data = odeint(self.to_coupled_linear, self.y0,
                                      self.samples)

    def draw(self, time_shown, show_analytical=False, painter=plt,
             fmt='k', label=None):
        if label is None:
            label = self.y_name + ' - numerical'
        """Plot the state of the current simulation. Self-explanatory"""
        if self.simulation_data is None:
            print('Please run the simulation first. ')
        else:
            beg = int(time_shown[0]*self.rate)
            fin = int(time_shown[1]*self.rate)
            xs = self.samples[beg:fin]
            ys = self.simulation_data[beg:fin, 0]
            painter.plot(xs, ys, fmt, label=label)
            if show_analytical and self.ans is not None:
                yas = self.ans(xs)
                painter.plot(xs, yas, 'b-',
                             label=self.y_name + ' - analytical')

    def draw_energies(self, show_analytical=False, painter=plt):
        if self.simulation_data is None:
            print('Please run the simulation first')
        else:
            if self.total_energy_stored is not None:
                xs = self.samples[:]
                ys = self.simulation_data[:, 0]
                yds = self.simulation_data[:, 1]
                energy = zeros(len(xs))
                painter.plot(xs, 1 -
                             self.total_energy_stored(ys, yds) /
                             self.initial_energy,
                             'b', label=self.y_name + ' - numerical')
                if show_analytical:
                    painter.plot(xs, energy, 'k-',
                                 label=self.y_name + ' - analytical')

    def draw_phase_space_plot(self, painter=plt, label=''):
        if self.simulation_data is None:
            print('Please run the simulation first')
        else:
            painter.plot(overwrap(self.simulation_data[:,0]),
                         self.simulation_data[
                :,1], label=label)


def overwrap(data):
    return fmod(data + pi, 2*pi) - pi


class Pendulum(MechanicalSystem):
    """
    A type of Mechanical system.

    Parameters
    ----------
    w_0 : float
    Natural frequency
    w_d : float
    Driving frequency
    q : float
    Dissipation Constant
    f : float
    Driving Force coefficient

    """
    y_unit = 'rad.'
    y_name = r'$\theta$'
    ydot_name = r'$\dot{\theta}$'

    def __init__(self, w_0=1, w_d=2/3, q=0, f=0, y0=None):
        """Initialise pendulum class."""
        if y0 is None:
            y0 = [0.01, 0]
        MechanicalSystem.__init__(self, lambda y_dot: -q*y_dot,
                                  lambda y: -(w_0**2)*sin(y),
                                  lambda t: f*sin(w_d*t), y0,
                                  ans=lambda t: y0[0]*cos(w_0*t) + (
                                          y0[1]/w_0)*sin(w_0*t))
        # Since we haven't made the small angle approximation,
        # need proper potential.
        self.total_energy_stored = lambda y, yd: (yd**2)/2 + (
                w_0**2*(1 - cos(y)))
        self.initial_energy = 1/2*(y0[1])**2 + w_0**2*(1 - cos(y0[0]))