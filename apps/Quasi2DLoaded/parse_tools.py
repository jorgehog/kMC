__author__ = 'jorgehog'

from numpy import *

def interpolate(value_array, time_array, time, start=0):

    if start >= len(time_array):
        return start, 0, 0

    i = start
    while time_array[i] < time:
        i += 1

        if i == len(time_array):
            return i, 0, 0

    if i == 0:

        if time != time_array[0]:
            return 0, 0, 0

        return i, 1, value_array[0]

    incline = (value_array[i] - value_array[i-1])/(time_array[i] - time_array[i-1])

    return i, 1, value_array[i-1] + incline*(time - time_array[i-1])


def time_derivative(time, data, pad=False):

    N = len(time)

    if pad:
        derivative = zeros(N)
    else:
        derivative = zeros(N-2)

    for i in range(1, N-1):
        derivative[i] = (data[i+1] - data[i-1])/(time[i+1]-time[i-1])

    if pad:
        derivative[0] = derivative[1]
        derivative[-1] = derivative[-2]

    return derivative


def combine_results(time_array, value_array, do_roughen=True):

    times = set()
    for time in time_array:
        time = asarray(time)
        times = times.union(set(time))

    times = array(sorted(list(times)))

    if do_roughen:
        times = roughen_single(times)

    all_f_t = zeros([len(value_array), len(times)])
    scales = zeros(len(times))
    starts = zeros(len(value_array))

    length = 40
    for ti, time in enumerate(times):
        for i in range(len(value_array)):
            loc, scale, all_f_t[i][ti] = interpolate(value_array[i], time_array[i], time, starts[i])

            scales[ti] += scale
            starts[i] = loc

        if scales[ti] == 0:
            print "ERROR IN SCALING", ti, time
            sys.exit(1)

        s = (ti*length)/(len(times)-1)
        print "\rcombining: [%s%s]" % ("-"*s, " "*(length-s)),
        sys.stdout.flush()

    return times, all_f_t.sum(0)/scales


def roughen_single(_array):

    idx = [0]
    k = 1
    N = len(_array)

    while round(k) < N:
        idx.append(int(round(k)))
        k += 0.001*k*scalar_log(1 + k)

    idx = unique(asarray(idx, dtype=int))

    return asarray(_array)[idx]


def roughen(*args):

    n_args = len(args)

    new = []
    for _array in args:
        new.append(roughen_single(_array))

    return new
