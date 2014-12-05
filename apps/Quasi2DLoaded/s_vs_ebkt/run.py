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

    eqVal = ParameterSet(1, 1, 1, cfg, "concEquil\s*\=\s*(.*)\;")
    concMult = ParameterSet(0.5, 0.5, 1, cfg, "concMult\s*\=\s*(.*)\;")
    h0_values = ParameterSet(2, 32, 2, cfg, "h0\s*\=\s*(.*)\;", lambda x, i: i*x)
    eb_values = ParameterSet(0.2, 5, 0.1, cfg, "Eb\s*=\s*(.*)\;")

    #controller.register_parameter_set(eqVal)
   # controller.register_parameter_set(concMult)
    #controller.register_parameter_set(h0_values)
    controller.register_parameter_set(eb_values)

    controller.run(run_kmc, controller, this_dir, path, app, dump_file_name)

if __name__ == "__main__":
    main()