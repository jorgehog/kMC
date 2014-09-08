from sys import argv
from os.path import join
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
    app = "Quasi2DLoaded"
    path = "../../build-kMC-Desktop_Qt_5_3_GCC_64bit-Release/apps/Quasi2DLoaded"
    cfg = join(path, "infiles", app + ".cfg")

    controller = ParameterSetController()

    lengths = ParameterSet(100, 150, 50, cfg, "BoxSize\s*=\s*\[(\d+), .*\]")
    temperatures = ParameterSet(0.5, 1.0, 0.5, cfg, "beta\s*=\s*(.*)\;")

    controller.register_parameter_set(lengths)
    controller.register_parameter_set(temperatures)

    controller.run(run_kmc, controller, this_dir, path, app)


if __name__ == "__main__":
    main()