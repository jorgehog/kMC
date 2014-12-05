from distutils.command.config import dump_file
from sys import argv, path
from os.path import join, split, exists
from os import (system,
                chdir,
                getcwd,
                remove)

from run_app_paramloop import *


def run_kmc(controller, this_dir, path, app, dump_file_name, **kwargs):

    print "Running ", controller
    chdir(path)
    system("./%s >> %s" % (app, dump_file_name))
    chdir(this_dir)


def main():
    this_dir = getcwd()

    path, app = split(argv[1])

    dump_file_name = "/tmp/kmc_dump.txt"

    if exists(dump_file_name):
        remove(dump_file_name)

    cfg = join(path, "infiles", app + ".cfg")

    controller = ParameterSetController()

    no_eqVal = ParameterSet(0, 0, 1, cfg, "concEquil\s*\=\s*(.*)\;")
    temperatures = ParameterSet(0.1, 1., 0.1, cfg, "beta\s*=\s*(.*)\;")
    concentrations = ParameterSet(0.1, 0.9, 0.01, cfg, "SaturationLevel\s*=\s*(.*)\;")
    eb_values = ParameterSet(0.1, 1., 0.1, cfg, "Eb\s*=\s*(.*)\;")

    controller.register_parameter_set(no_eqVal)
    controller.register_parameter_set(temperatures)
    controller.register_parameter_set(concentrations)
    controller.register_parameter_set(eb_values)

    controller.run(run_kmc, controller, this_dir, path, app, dump_file_name)

    print "equilibriating"

    eqController = ParameterSetController()

    eqVal = ParameterSet(1, 1, 1, cfg, "concEquil\s*\=\s*(.*)\;")
    reset = ParameterSet(0, 0, 1, cfg, "reset\s*\=\s*(.*)\;")

    eqController.register_parameter_set(eqVal)
    eqController.register_parameter_set(reset)

    eqController.register_parameter_set(temperatures)
    eqController.register_parameter_set(eb_values)

    eqController.run(run_kmc, eqController, this_dir, path, app, dump_file_name, ask=False)

if __name__ == "__main__":
    main()