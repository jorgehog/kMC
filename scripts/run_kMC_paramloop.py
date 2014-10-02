from sys import argv
from os.path import join, split
from os import (system,
                chdir,
                getcwd,
                remove)

from run_app_paramloop import *


def run_kmc(controller, this_dir, path, app):

    print "Running ", controller
    chdir(path)
    system("./%s >> /tmp/kmc_dump.txt" % app)
    chdir(this_dir)


def main():
    this_dir = getcwd()
    path, app = split(argv[1])

    remove("/tmp/kmc_dump.txt")

    cfg = join(path, "infiles", app + ".cfg")
    repeat_count = 10

    controller = ParameterSetController()

    lengths = ParameterSet(128, 512, 2, cfg, "BoxSize\s*=\s*\[(\d+), .*\]", lambda p, i: p*i)
    temperatures = ParameterSet(0.1, 5.0, 0.05, cfg, "beta\s*=\s*(.*)\;")
    es_maxes = ParameterSet(2.0, 20.0, 1.0, cfg, "EsMax\s*\=\s*(.*)\;")
    concentrations = ParameterSet(0.4, 0.4, 0.1, cfg, "SaturationLevel\s*=\s*(.*)\;")
    run_id = ParameterSet(1, repeat_count, 1, cfg, "runID\s*=\s*(.*)\;")
    eb_values = ParameterSet(0.01, 1, 0.05, cfg, "Eb\s*=\s*(.*)\;")


    #controller.register_parameter_set(lengths)
    #controller.register_parameter_set(temperatures)
    # controller.register_parameter_set(es_maxes)
    #controller.register_parameter_set(concentrations)
    controller.register_parameter_set(eb_values)
    controller.register_parameter_set(run_id)

    controller.run(run_kmc, controller, this_dir, path, app)


if __name__ == "__main__":
    main()