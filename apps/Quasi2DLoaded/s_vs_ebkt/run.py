from sys import argv, path
from os.path import join, split, exists
from os import (system,
                chdir,
                getcwd,
                remove)
import time

path.append(join(getcwd(), "..", ".."))

from run_app_paramloop import ParameterSet, ParameterSetController, quick_replace


def run_kmc(proc, combination, this_dir, path, app, dump_file_name):
    time.sleep(proc/10.)

    print "Running ",
    for value in combination:
        print "%.3f" % value,
    print

    chdir(path)
    system("./%s %d >> %s_%d.txt" % (app, proc, dump_file_name, proc))
    chdir(this_dir)


def main():
    this_dir = getcwd()

    path, app = split(argv[1])

    dump_file_name = "/tmp/kmc_dump"

    if exists(dump_file_name):
        remove(dump_file_name)

    cfg = join(path, "infiles", app + ".cfg")


    if len(argv) > 2:
        outpath = argv[2]
        quick_replace(cfg, "path", outpath)

    controller = ParameterSetController()

    E0_values = ParameterSet(cfg, "Ew0dL\s*\=\s*(.*)\;")
    E0_values.initialize_set_incr(0.05, 0.251, 0.025)

    alpha_values = ParameterSet(cfg, "alpha\s*=\s*(.*)\;")
    alpha_values.initialize_set_incr(0.2, 2.01, 0.2)

    controller.register_parameter_set(E0_values)
    controller.register_parameter_set(alpha_values)

    controller.run(run_kmc, this_dir, path, app, dump_file_name, n_procs=8)


if __name__ == "__main__":
    main()