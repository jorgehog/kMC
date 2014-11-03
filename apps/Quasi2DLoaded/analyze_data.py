from distutils.fancy_getopt import neg_alias_re
import h5py
from hgext.mq import clone
import sys
import re
from gst._gst import TIME_ARGS
from matplotlib.artist import allow_rasterization
from matplotlib.pylab import *
from math import log as scalar_log
from numpy import where, asarray, random, log, empty



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


def get_w_beta(rms, times, do_plot=False):

    t, f_t = combine_results(times, rms, do_roughen=True)

    from scipy import polyfit

    M = 10000000000000000
    M2 = M

    logt = log(t)
    logf = log(f_t)
    N = len(logt)

    I = None
    I2 = None
    R = None
    R2 = None
    for i in range(N/10, N):
        r, residuals, c, c, c = polyfit(logt[:i], logf[:i], 1, full=True)

        if residuals[0] < M:
            M = residuals[0]
            R = r
            I = i

        i2 = N - i
        r2, residuals2, c, c, c = polyfit(logt[i2:-1], logf[i2:-1], 1, full=True)

        if residuals2[0] < M2:
            M2 = residuals2[0]
            R2 = r2
            I2 = i2

    if do_plot:
        cla()
        loglog(t, f_t)

        figure()
        plot(logt, logf)
        plot(logt[:I], R[1] + logt[:I]*R[0])
        plot(logt[I2:], R2[1] + logt[I2:]*R2[0])

        print "w inf = %g, beta = %g" % (R[0], R2[1])

        show()

    return R[0], R2[1]


def main():

    f = h5py.File(sys.argv[1], 'r')

    all_data = {}

    c0 = 0.4
    em0 = 1
    t0 = 1
    n0 = 1
    l0 = 128

    total = len(f.keys())
    k = 0
    for l, run in f.items():

        print "l=%4s %2d / %2d Parsing %d items ..." % (l, k + 1, total, len(run.items()))

        for potential, data in run.items():

            if not "concentration" in potential:
                potential += "_concentration_0.4"

            eb, t, em, c = [float(re.findall("%s\_(\d+\.?\d*)" % ID, potential)[0]) for ID in ["Eb",
                                                                                               "beta",
                                                                                               "EsMax",
                                                                                               "concentration"]]

            n = re.findall("\_n_(\d+)", potential)

            if n:
                n = int(n[0])
            else:
                n = ""

            #
            # # print c,c0, em, em0*eb, t, t0, l, l0
            # if c != c0 or em != em0*eb or t != t0:
            #     # if c != c0:
            #     #     print "c"
            #     # if em != em0*eb:
            #     #     print "em"
            #     # if t != t0:
            #     #     print "t"
            #     # if l != l0:
            #     #     print "l"
            #     continue

            if eb not in all_data.keys():
                all_data[eb] = {}

            if t not in all_data[eb].keys():
                all_data[eb][t] = {}
            if em not in all_data[eb][t].keys():
                all_data[eb][t][em] = {}
            if c not in all_data[eb][t][em].keys():
                all_data[eb][t][em][c] = {}
            if l not in all_data[eb][t][em][c].keys():
                all_data[eb][t][em][c][l] = {}

            all_data[eb][t][em][c][l][n] = {}

            for i, name in enumerate(data["ignisEventDescriptions"][0]):

                if name.split("@")[0] == "SurfaceSize":
                    all_data[eb][t][em][c][l][n][str(name).split("@")[0]] = data["ignisData"][i, :]

            all_data[eb][t][em][c][l][n]["name"] = potential

            # all_data[eb][t][em][c][l][n]["AutoCorr"] = array(data["AutoCorr"])

            all_data[eb][t][em][c][l][n]["AutoCorr"] = array([0])

            if "eqConc" in data.attrs.keys():
                all_data[eb][t][em][c][l][n]["eqConc"] = data.attrs["eqConc"]


        k += 1

    f.close()

    all_times = []
    all_lengths = []
    all_rms = []

    all_temps = []
    all_c_eq = []
    all_em_maxes = []
    all_eb = []

    all_beta_eb = []
    all_cumz = []
    all_sizes = []
    # cum_mat = zeros(shape=(5,5))

    i = 0
    for eb, init in sorted(all_data.items(), key=lambda x: x[0]):
        j = 0

        for t, rest in sorted(init.items(), key=lambda x: x[0]):
            for em, rest2 in sorted(rest.items(), key=lambda x: x[0]):
                for c, rest3 in sorted(rest2.items(), key=lambda x: x[0]):

                    name = rest3.values()[0].values()[0]["name"].split("_n_")[0]

                    for length, all_runs in sorted(rest3.items(), key=lambda x: x[0]):

                        # time_array = [data["Time"] for data in all_runs.values()]
                        # rms_array = [data["heightRMS"] for data in all_runs.values()]
                        # h_array = [data["height"] for data in all_runs.values()]
                        # acf = [data["AutoCorr"] for data in all_runs.values()]
                        # cum = [data["cumulant"] for data in all_runs.values()]
                        sizes = [data["SurfaceSize"] for data in all_runs.values()]

                        try:
                            c_eq = sum([data["eqConc"] for data in all_runs.values()])/len(all_runs.values())
                        except KeyError:
                            c_eq = 0

                        print eb, t

                        # cumz = sum([c[where(c != 0)].mean() for c in cum])/len(cum)
                        surfacesize = sum(s[-1] for s in sizes)/len(sizes)
                        # cum_mat[i, j] = cumz

                        all_beta_eb.append(t*eb)
                        # all_cumz.append(cumz)
                        all_sizes.append(surfacesize)

                        all_c_eq.append(c_eq)
                        continue

                        all_em_maxes.append(em)
                        all_temps.append(t)
                        all_eb.append(eb)

                        toss = []
                        for i, t in enumerate(time_array):
                            if where(t < 0)[0].size != 0:
                                print "Negative times in run ", i
                                toss.append(i)

                        for i in toss[::-1]:
                            time_array.pop(i)
                            rms_array.pop(i)

                        t, f_t = combine_results(roughen(*time_array), roughen(*rms_array), do_roughen=True)

                        figure(0)
                        hold("on")
                        loglog(t, f_t, label=length)

                        t, f_t = combine_results(roughen(*time_array), [h_i/t_i for h_i, t_i in zip(h_array, time_array)])

                        K = 10
                        L = len(t)/K

                        figure(1)
                        hold("on")
                        plot(t[L:], f_t[L:], label=length)

                        if int(length) == 512:

                            a = []

                            figure(2)
                            hold("on")

                            for t_i, h_i in zip(time_array, h_array):
                                v = h_i/t_i

                                l_i = len(h_i)/K
                                plot(t_i[l_i:], v[l_i:])
                                a.append(v[l_i:].mean())

                            figure(3)
                            hist(a, max([len(time_array)/5, 20]))

                        D = 20
                        acf_tot = zeros(D)
                        for acf_i in acf:
                            acf_tot += acf_i[:D]
                        acf_tot /= len(acf)

                        figure(4)
                        hold("on")
                        plot(acf_tot, label=length)

                    figure(0)
                    title(name)

                    legend(loc=2)
                    xlabel("t [s * R_0]")
                    ylabel("<RMS(h)>")
                    draw()
                    savefig(name + "_logT_logRMS.png")

                    figure(1)
                    legend(loc=2)
                    xlabel("t [s * R_0]")
                    ylabel("<v(t)>")
                    xscale('log')
                    draw()
                    savefig(name + "_logT_v.png")

                    figure(2)
                    draw()
                    savefig("test.png")

                    figure(3)
                    draw()
                    savefig("v0dist.png")

                    figure(4)
                    legend()
                    xlabel("dx")
                    ylabel("acf(dx)")
                    draw()
                    savefig("acf.png")

                    # show()


                    # winf, beta = get_w_beta(all_rms, all_times)
            j += 1
        i += 1

    close("all")
    figure()
    plot(all_beta_eb, all_sizes, 'k*')
    xlabel("beta*eb")
    ylabel("<s>")
    savefig("sizes_betaeb.png")

    save("all_beta_eb", all_beta_eb)
    save("all_sizes", all_sizes)

    # figure()
    # imshow(cum_mat, extent=(0.1, 0.5, 0.1, 0.9))
    # xlabel("beta")
    # ylabel("eb")
    # colorbar()
    # savefig("kumulant_beta_eb.png")

    # figure()
    # plot(all_beta_eb, all_cumz, 'k*')
    # xlabel("beta*eb")
    # ylabel("k")
    # savefig("kumulant_betaeb.png")

    save("all_eqc", all_c_eq)

    figure()
    plot(all_beta_eb, all_c_eq, 'kx')
    xlabel("Eb/kT")
    ylabel("C_0/l_0")
    savefig("eqc_betaeb.png")
    show()

    return

    # if len(all_temps) == len(all_c_eq) and all_temps:
    #     all_temps_inv = [1.0/t for t in all_temps]
    #     plot(all_temps, all_c_eq, 'kx')
    #     xlabel("kT/E_0")
    #     ylabel("C_0/l_0")
    # elif len(all_em_maxes) == len(all_c_eq) and all_em_maxes:
    # plot(all_em_maxes, all_c_eq, 'kx')
    # xlabel("Em/E_0")
    # ylabel("C_0/l_0")
    plot(all_eb, all_c_eq, 'kx')
    xlabel("Eb/E0")
    ylabel("C_0/l_0")
    show()

    #
    # close("all")
    # for t, rest in all_data.items():
    #     for em, rest2 in rest.items():
    #         for c, lengths in rest2.items():
    #             figure()
    #
    #             name = lengths.values()[0]["name"]
    #
    #             hold("on")
    #             for length, data in lengths.items():
    #                 plot(data["Time"], data["height"]/data["Time"], label=length)
    #
    #             xscale("log")
    #             title(name)
    #
    #             legend(loc=2)
    #             xlabel("t [s * R_0]")
    #             ylabel("<v(t)>")
    #             draw()
    #
    #
    # show()




if __name__ == "__main__":
    main()