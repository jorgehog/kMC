from sys import argv
from os.path import join, split
from os import (system,
                chdir,
                getcwd)


from run_app_paramloop import *

def run_kmc(controller, this_dir, path, app):

    print "Running ", controller
    chdir(path)
    system("./%s > /tmp/kmc_dump.txt" % app)
    chdir(this_dir)

def main():
    this_dir = getcwd()
    path, app = split(argv[1])

    cfg = join(path, "infiles", app + ".cfg")
    nLoops = 300

    controller = ParameterSetController()

    lengths = ParameterSet(128, 1024, 2, cfg, "BoxSize\s*=\s*\[(\d+), .*\]", lambda p, i: p*i)
    temperatures = ParameterSet(0.1, 1.0, 0.1, cfg, "beta\s*=\s*(.*)\;")
    es_maxes = ParameterSet(10.0, 10.0, 2.0, cfg, "EsMax\s*\=\s*(.*)\;", lambda p, i: p*i)
    concentrations = ParameterSet(0.4, 0.4, 0.1, cfg, "SaturationLevel\s*=\s*(.*)\;")
    IDS = ParameterSet(1, nLoops, 1, cfg, "runID\s*=\s*(.*)\;")

    controller.register_parameter_set(lengths)
    controller.register_parameter_set(temperatures)
    controller.register_parameter_set(es_maxes)
    controller.register_parameter_set(concentrations)
    controller.register_parameter_set(IDS)

    controller.run(run_kmc, controller, this_dir, path, app)


if __name__ == "__main__":
    main()