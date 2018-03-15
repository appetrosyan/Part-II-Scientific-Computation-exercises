#!/usr/local/bin/python3
# c-basic-offset: 4; tab-width: 4; indent-tabs-mode: nil
# vi: set shiftwidth=4 tabstop=4 expandtab
# :indentSize=4:tabSize=4:noTabs=true

from genericpath import exists
from os import mkdir
from MechanicalSystems import Pendulum
import matplotlib.pyplot as plt
from numpy import pi, mean, diff, ravel, nonzero, logical_and, linspace


def compare_with_theory(cycles=None, p=None, check=True, sampling_rate=500,
                        bounds=None):
    """A Helper function to generate the plots comparing numerical
    integration to theory"""
    if cycles is None:
        cycles = [10, 100, 1000]
    if p is None:
        p = Pendulum(1, 2/3, 0, 0, y0=[0.01, 0])
    if bounds is None:
        bounds = [[t, t + 10] for t in [2*pi*x for x in cycles]]

    p.simulate(rate=sampling_rate, duration=(cycles[-1] + 10.0)*2*pi)
    fig, axes = plt.subplots(len(cycles))
    try:
        _ = iter(axes)
    except TypeError:
        axes = [axes]
    plt.suptitle('Trace of the numerical solution')
    plt.xlabel('t/s')
    plt.ylabel('deflection / radians')
    for (b, ax) in zip(bounds, axes):
        p.draw(b, painter=ax, show_analytical=check, fmt='r')
    plt.legend(loc='best')
    save_figure('comparison_with_analytical-deflection')
    plt.show()

    plt.figure()
    p.draw_energies(show_analytical=check)
    plt.xlabel('time/s')
    plt.ylabel('fraction of energy lost (dimensionless)')
    plt.title('Total energy evolution, for small initial deflection ')
    save_figure('comparison_with_analytical-energy')
    plt.show()


def save_figure(file_name):
    if not exists('figures/'):
        mkdir('figures/')
    plt.savefig('figures/%s.pdf' % file_name)


# Core Task 1.2: Finding the period of oscillations and plotting
def freq_from_crossings(data, sampling_rate=500):
    """
    Estimate frequency by counting zero crossings.

    Returns the most likely frequency in Hertz

    """
    indices, = nonzero(ravel(logical_and(data[1:] >= 0, data[:-1] < 0)))
    crossings = [i - data[i]/(data[i + 1] - data[i]) for i in indices]
    return sampling_rate/mean(diff(crossings))


def investigate_period_amplitude(at=pi/2):
    """"Plot period's dependence on initial conditions"""
    domain = linspace(0.0001, pi, 50)
    period_from_position, period_from_velocity = [], []
    sampling_rate = 100
    for i in domain:
        p = Pendulum(1, 2/3, 0, 0, [i, 0])
        p.simulate(rate=sampling_rate, duration=100*2*pi)
        position_data = p.simulation_data[:, 0]
        velocity_data = p.simulation_data[:, 1]
        period_from_position.append(
            1/freq_from_crossings(position_data, sampling_rate))
        period_from_velocity.append(
            1/freq_from_crossings(velocity_data, sampling_rate))
    fig, ax = plt.subplots()
    plt.plot(domain, period_from_position[:], 'b.', label=r'from $\theta$')
    plt.plot(domain, period_from_velocity[:], 'k+', label=r'from $\dot{'
                                                          r'\theta}$')
    plt.xlabel(r'$\theta_0$ / radians')
    plt.xticks([pi/8, pi/4, pi/2, pi])
    ax.set_xticklabels(
        [r'$\frac{\pi}{8}$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$',
         r'$\pi$'])
    plt.ylabel('Period / s')

    p = Pendulum(1, 2/3, 0, 0, [at, 0])
    p.simulate(rate=500, duration=1000)
    ys = p.simulation_data[:, 1]  # Use velocities, because of zero offset
    period = 1/freq_from_crossings(ys, 500)
    plt.legend(loc='best')
    plt.plot(at, period, 'ko')
    plt.annotate(r'Period at $\pi/2$ = ' + "{:.3f}".format(period),
                 xy=(at, period), xytext=(0, 30),
                 arrowprops=dict(facecolor='black', shrink=0.05))
    plt.title("Dependence of period on initial conditions. ")
    save_figure('period_vs_amplitude')
    plt.show()
    return period


def investigate_damping(qs=None, sampling_rate=500):
    """Investigate the period's dependence on damping"""
    if qs is None:
        qs = [1, 5, 10]
    fig, axes = plt.subplots(len(qs), 2, sharey='col', sharex='row')
    plt.subplots_adjust(wspace=0.5, hspace=0.05)
    plt.suptitle('Comparison of different damping regimes')
    for q, ax0, ax1 in zip(qs, axes[:, 0], axes[:, 1]):
        pend = Pendulum(q=q)
        pend.simulate(rate=sampling_rate)
        length = max(qs)*1.71
        pend.draw([0, length], painter=ax0)
        ax0.set_title('q = ' + str(q), y=0.18, x=1.05)
        pend.draw_energies(painter=ax1)
    last = len(qs) - 1
    axes[int(last/2), 0].set_ylabel(r'$\theta$ / radians', y=0.5)
    axes[last, 0].set_xlabel('t/s')
    axes[int(last/2), 1].set_ylabel(
        r'fraction of energy lost (dimensionless)', y=0.5)
    axes[last, 1].set_xlabel('t/s')
    save_figure('damping')
    plt.show()


def investigate_driving(fs=None, rate=500, filename='driving'):
    """Investigate the effect of different driving amplitudes"""
    data = []
    if fs is None:
        fs = [0.5, 1.2, 1.44, 1.465]
    fig, axes = plt.subplots(len(fs), 2, sharex='col')
    plt.subplots_adjust(wspace=1, hspace=0.05)
    plt.suptitle('Comparison of different driving regimes')
    for f, ax0, ax1 in zip(fs, axes[:, 0], axes[:, 1]):
        p = Pendulum(f=f)
        time = 100
        p.simulate(rate=rate, duration=time)
        sim_data = p.simulation_data[:, 1]
        data.append(1/freq_from_crossings(sim_data))
        ax0.set_title('f = ' + str(f), y=0.18, x=1.05)
        p.draw(time_shown=[0, time], painter=ax0)
        p.draw_energies(painter=ax1)
    last = len(fs) - 1
    axes[int(last/2) + 1, 0].set_ylabel(r'$\theta$ / radians')
    axes[int(last/2) + 1, 1].set_ylabel(
        r'fraction of energy lost (dimensionless)')
    axes[last, 1].set_xlabel('t/s')
    axes[last, 0].set_xlabel('t/s')
    save_figure(filename + '-deflections_and_energies')
    plt.show()

    plt.figure()
    plt.plot(fs, data, 'k+')
    plt.title("Dependence of period on Driving amplitude. ")
    plt.xlabel('Driving amplitude / (rad $s^{-2}$)')
    plt.ylabel('Period/s')
    save_figure(filename + '-period_vs_amplitude')
    plt.show()


def investigate_sensitivity():
    p1 = Pendulum(1, 2/3, 0.5, 1.2, [0.2, 0])
    p2 = Pendulum(1, 2/3, 0.5, 1.2, [0.20001, 0])
    p1.simulate(500, 1000)
    p2.simulate(500, 1000)
    p1.draw([80, 100], fmt='r', label='0.2')
    p2.draw([80, 100], fmt='b', label=r'0.2 + $\epsilon$')
    plt.title('Comparison of solutions with different initial conditions')
    plt.xlabel('t/s')
    plt.ylabel('deflection/radians')
    plt.legend(loc='best')
    save_figure('sensitivity')
    plt.show()

    plt.figure()
    p1.draw_phase_space_plot()
    p2.draw_phase_space_plot()
    plt.legend(loc='best')
    save_figure('sensitivity-phase')
    plt.show()


def investigate_chaos():
    p = Pendulum()
    p.simulate(duration=100)
    p.draw_phase_space_plot(label='No damping')
    for q in [.2, 1, 1.5, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2]:
        p = Pendulum(q=q)
        p.simulate(duration=100)
        p.draw_phase_space_plot(label='q = ' + str(q))
    plt.legend(loc='best')
    plt.xlabel(Pendulum.y_name)
    plt.ylabel(Pendulum.ydot_name)
    save_figure('phase_space-damping')
    plt.show()
    for f in [0, 0.02, 0.01, 0.03, 0.05]:
        p = Pendulum(f=f)
        p.simulate(duration=100)
        p.draw_phase_space_plot(label='f =' + str(f))
    plt.xlabel(Pendulum.y_name)
    plt.ylabel(Pendulum.ydot_name)
    plt.legend(loc='best')
    save_figure('phase_space-light_driving')
    plt.show()
    for f in [0, 0.5]:
        p = Pendulum(f=f)
        p.simulate(duration=100)
        p.draw_phase_space_plot(label='f =' + str(f))
    plt.xlabel(Pendulum.y_name)
    plt.ylabel(Pendulum.ydot_name)
    plt.legend(loc='best')
    save_figure('phase_space-heavy_driving')
    plt.show()


if __name__ == '__main__':
    compare_with_theory()
    investigate_period_amplitude()
    investigate_damping()
    investigate_driving()
    investigate_driving([0.01, 0.02, 0.05, 0.1], filename='weak_driving')
    investigate_sensitivity()
    investigate_chaos()
