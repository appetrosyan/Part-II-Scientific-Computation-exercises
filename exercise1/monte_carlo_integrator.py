#!/usr/local/bin/python3
# c-basic-offset: 4; tab-width: 4; indent-tabs-mode: nil
# vi: set shiftwidth=4 tabstop=4 expandtab
# :indentSize=4:tabSize=4:noTabs=true
from functools import partial
from multiprocessing.pool import Pool
import matplotlib.pyplot as plt
from numpy import pi, sum, sin, exp, log, vectorize, average, std, fromfunction,\
        polyfit, poly1d, sqrt
from numpy.random import uniform as uniform
from os.path import exists
from os import mkdir


# Core Task 1: Monte Carlo integration

S = pi/8
D = 8


def f_example(arg: float) -> float:
    """Definition of integrand, given in exercise."""
    return (10**6)*sin(sum(arg, axis=1))


def expect(fun: callable, args, norm: float = 1):
    """Compute expectation value.

    Parameters
    ----------
    fun: callable
    of which expectation value is computed
    args: iterable
    of sample points.
    norm: float, optional
    argument to normalise the output with.

    """
    return sum(fun(args))*norm


def monte_carlo_integrate(number_of_samples: int,
                          integrand: callable = f_example,
                          dimensionality_of_space: int = D,
                          box_side_length: int = S, stats=False):
    """
    Integrate a dimensionality_of_space-dimensional function on a box of size s.

    Uses monte_carlo integration techniques.
    Parameters
    ----------
    number_of_samples : int
    Number of dimensionality_of_space-dimensional sample vectors to use.
    integrand : callable, optional
    function to be integrated.
    dimensionality_of_space : int, optional
    dimensionality of the integral.
    box_side_length : float, optional
    size of integration box.
    stats
    Whether to calculate the error as well as the value.
    Returns
    -------
    (integral, error): tuple(float, float)
    given stats=True
    or
    integral: float

    """
    samples = uniform(low=0.0, high=box_side_length,
                      size=(number_of_samples + 1, dimensionality_of_space))
    volume = box_side_length**dimensionality_of_space
    reciprocal = 1.0/number_of_samples
    mean = expect(integrand, samples, reciprocal)
    msq = expect(lambda arg: integrand(arg)**2, samples, reciprocal)
    integral = mean*volume
    error = volume*sqrt((msq - mean)*reciprocal)
    if stats:
        return integral, error
    else:
        return integral


def find_best_value_mp(samples_per_iteration, num_of_iterations=25,
                       stats=False):
    """
    Iterate over monte-carlo integrations to compute an estimate of error.
    Parameters
    ----------
    samples_per_iteration: int
    num_of_iterations: int
    stats: bool

    Returns
    -------

    """

    data = [samples_per_iteration for _ in range(num_of_iterations)]
    with Pool() as p:
        xs = p.map(partial(monte_carlo_integrate, stats=stats), data)
    if stats:
        integrals = [x[0] for x in xs]
        theoretical_errors, norm = [x[1] for x in xs], sqrt(1)/num_of_iterations
        return average(integrals), std(integrals), average(
            theoretical_errors)*norm
    else:
        return average(xs), std(xs)


# --Tabulating the results for plotting
def tabulate_estimates(samples: iter, n: int, with_stats: bool = False):
    """Helper function to tabulate results many calculations."""
    est = vectorize(lambda arg: find_best_value_mp(arg, n, stats=with_stats))
    return est(samples)


def plot_integral_value(number_of_iterations, table_of_estimates) -> None:
    """Plot the convergence of the integral as a function of number of
    iterations. """
    plt.figure()
    plt.xlabel('Number of sample points')
    plt.ylabel('Integral/arb. units')
    plt.title('Monte-carlo integration.\n Value estimates and convergence to '
              'known analytical result')
    plt.yscale('log')
    plt.xscale('log')
    plt.errorbar(number_of_iterations, table_of_estimates[0],
                 yerr=table_of_estimates[1], fmt='+', color='k', capsize=6,
                 label='Errors from standard deviation of estimates')
    if len(table_of_estimates) > 2:
        plt.errorbar(number_of_iterations, table_of_estimates[0],
                     yerr=table_of_estimates[2], fmt='+', color='0.5',
                     capsize=3, label='Errors from numerical error estimate')

    plt.axhline(y=537.1873411, color='0.3', linestyle='-',
                label='Analytical solution')
    plt.legend(loc='best')
    if not exists('figures'):
        mkdir('figures')
    plt.savefig('figures/Mc-values.pdf')
    plt.show()


def log_log_error_plot(number_of_iterations, table_of_estimates, labels=None):
    """Produce a log-log plot of the obtained numerical integration errors
    as a function of the number of iterations. Includes linear regression
    analysis, and comparison to higher order estimates given from Monte-Carlo
    integration theory."""
    if labels is None:
        labels = ['Monte Carlo integration - errors',
                  'Number of sample points / arb units',
                  'Error in integral / arb units']
    plt.figure()
    plt.title(labels[0])
    plt.xlabel(labels[1])
    plt.ylabel(labels[2])
    plt.yscale('log')
    plt.xscale('log')

    plt.plot(number_of_iterations, table_of_estimates[1], 'k+',
             label='Error estimates by standard deviation of integrals')
    label, y = linear_fit_on_log_log(number_of_iterations,
                                     table_of_estimates[1])
    plt.plot(number_of_iterations, y, 'k--', label=r'Best fit to $\sigma$')
    plt.text(1.25, 1.25, label)

    if len(table_of_estimates) > 2:
        plt.plot(number_of_iterations, table_of_estimates[2], 'bx',
                 label='Error estimates by numerical error estimate')
        label, y = linear_fit_on_log_log(number_of_iterations,
                                         table_of_estimates[2])
        plt.plot(number_of_iterations, y, 'b-.', label=r'Best fit to $\sigma$')
        plt.text(2.25, 2.25, label, color='b')
    plt.legend(loc='best')
    if not exists('figures'):
        mkdir('figures')
    plt.savefig('figures/Mc-error-bar-plots.pdf')
    plt.show()


def linear_fit_on_log_log(x_vector, y_vector):
    """Helper function to rescale a linear regression in logarithmic space
    back to linear space, to give an accurate plot. """
    x = log(x_vector)
    y = log(y_vector)
    coefficients = polyfit(x, y, 1)
    label = '$\sigma =' + '%.3f'%exp(coefficients[
                                         1]) + '\cdot ' + "n_{it}^{" + f'{coefficients[0]:.3f}' + "}" + "$"
    fit = poly1d(coefficients)
    y = exp(fit(x))
    return label, y


if __name__ == '__main__':
    sample_sets = fromfunction(lambda i: 2**i, (13,), dtype=int)
    table = tabulate_estimates(sample_sets, 40, with_stats=True)
    log_log_error_plot(sample_sets, table)
    plot_integral_value(sample_sets, table)