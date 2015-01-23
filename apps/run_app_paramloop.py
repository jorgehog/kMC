from itertools import product, combinations
from re import sub, findall
import threading
import time

class Worker(threading.Thread):

    def __init__(self, controller, exec_program_rule, set, proc, *args, **kwargs):

        self.controller = controller
        self.exec_program_rule = exec_program_rule
        self.set = set
        self.proc = proc

        self.args = args
        self.kwargs = kwargs

        super(Worker, self).__init__()

    def write_combination(self, combination):

        for parameter_set, value in zip(self.controller.parameter_sets,
                                        combination):
            parameter_set.write_value(value, self.proc)

    def run(self):
        for combination in self.set:
            self.write_combination(combination)

            self.exec_program_rule(self.proc, combination, *self.args, **self.kwargs)

class ParameterSet:

    def __init__(self,
                 config_filename,
                 variable_pattern_in_config,
                 regex_flags=[]):

        self.config_filename = config_filename
        self.variable_pattern_in_config = variable_pattern_in_config

        self.regex_flags = regex_flags

        self.current_index = 0

        self.child = None
        self.parent = None

    def initialize_set(self, set):

        self.set = set

    def initialize_set_incr(self,
                            start,
                            stop,
                            increment,
                            increment_rule=lambda previous, increment: previous + increment):

        set = []

        value = start

        while (value <= stop):

            set.append(value)

            value = increment_rule(value, increment)

        self.initialize_set(set)

    @staticmethod
    def invert_regex(pattern):

        inverted_pattern = ""

        if pattern[0] != "(":
            inverted_pattern += "("

        for i in range(0, len(pattern)):

            if pattern[i] == "(":
                inverted_pattern += ")"

            elif pattern[i] == ")":
                inverted_pattern += "("

            else:
                inverted_pattern += pattern[i]

        if pattern[-1] != ")":
            inverted_pattern += ")"

        return inverted_pattern.replace("()", "")

    def get_config_name(self, ending):
        if "." in self.config_filename:

            split = self.config_filename.split(".")
            pre = split[:-1]
            suf = split[-1]

            return ".".join(pre) + ending + "." + suf

        else:

            return self.config_filename + ending

    def write_value(self, value, proc):

        ending = "_%d" % proc
        config_out = self.get_config_name(ending)

        with open(self.config_filename, 'r') as config_file:

            prev_config_raw_text = config_file.read()

        new_config_raw_text = sub(self.invert_regex(self.variable_pattern_in_config),
                                  r"%s%s%s" % (r"\g<1>", value, r"\g<2>"),
                                  prev_config_raw_text,
                                  *self.regex_flags)

        with open(config_out, 'w') as config_file:

            config_file.write(new_config_raw_text)


class ParameterSetController:

    def __init__(self):
        self.parameter_sets = []

    def register_parameter_set(self, parameter_set):

        if len(self.parameter_sets) != 0:
            parameter_set.parent = self.parameter_sets[-1]
            self.parameter_sets[-1].child = parameter_set

        self.parameter_sets.append(parameter_set)

    def prepare(self):

        for parameter_set in self.parameter_sets:
            parameter_set.prepare()

    def run(self, execute_program_rule, *args, **kwargs):

        n = 1
        for parameter_set in self.parameter_sets:
            n *= len(parameter_set.set)

        skip = False
        if "ask" in kwargs:
            if kwargs["ask"] is False:
                skip = True
            kwargs.pop("ask")

        n_procs = 1
        if "n_procs" in kwargs:

            n_procs = kwargs.pop("n_procs")

            if n_procs <= 1:
                n_procs = 1


        if skip:
            ans = ""
        else:
            ans = raw_input("Press either: \n\tenter to run %d simulations."
                            "\n\tt [time per sim in seconds] for a time estimate. "
                            "\n\tq to end.\n" % n)

        if ans.startswith("t "):
            time = float(ans.lstrip("t "))
            t = n*time/n_procs

            unit = ""
            if t > 3600:
                t /= 3600.
                unit = "hours"

            elif t > 60:
                t /= 60.
                unit = "minutes"

            else:
                unit = "seconds"

            print "Expected time: %g %s" % (t, unit)

            return self.run(execute_program_rule, *args, **kwargs)

        elif ans == "q":
            print "exiting.."
            return

        combinations = list(product(*[pset.set for pset in self.parameter_sets]))

        n_per_proc = len(combinations)/n_procs

        remainder = len(combinations) - n_per_proc*n_procs

        all_threads = []

        prev = 0
        for proc in range(n_procs):

            N = n_per_proc

            if proc < remainder:
                N += 1

            chunk = combinations[prev:prev+N]

            prev += N

            thread = Worker(self, execute_program_rule, chunk, proc, *args, **kwargs)
            all_threads.append(thread)

            thread.start()

def quick_replace(cfg, name, value):

    cfg_str = ""
    with open(cfg, 'r') as f_read:
        cfg_str = f_read.read()

    repl = sub("(%s\s*=\s*[\"\']?).*?([\"\']?;)" % name, "\g<1>%s\g<2>" % str(value), cfg_str)

    with open(cfg, 'w') as f_write:
        f_write.write(repl)



def exec_test_function(proc, combination):

    time.sleep(proc)

    testfile_name = "/tmp/test_paramloop%d.cfg" % proc

    with open(testfile_name, 'r') as testfile:
        testfile_raw_text = testfile.read()

    results = findall(r"a\s*\=\s*\[\s*(\d+)\s*\,\s*(.*)\s*\,\s*(\d+)\s*\]", testfile_raw_text)[0]

    for i, result in enumerate(results):
        if float(result) == combination[i]:
            print "combination", combination, "success."
        else:
            print "combination", combination, "failed."


def testbed():

    testfile_name = "/tmp/test_paramloop.cfg"

    with open(testfile_name, 'w') as testfile:

        testfile.write("a=[0,1,2]\n")

    set1 = ParameterSet(testfile_name, r"a\=\[(\d+)\,.*\,\d+\]")
    set2 = ParameterSet(testfile_name, r"a\=\[\d+\,(.*)\,\d+\]")
    set3 = ParameterSet(testfile_name, r"a\=\[\d+\,.*\,(\d+)\]")

    set1.initialize_set_incr( 0,   3,   1)
    set2.initialize_set([-2, -0.5, 1])
    set3.initialize_set_incr( 0, 100,  30, )


    controller = ParameterSetController()

    controller.register_parameter_set(set1)
    controller.register_parameter_set(set2)
    controller.register_parameter_set(set3)

    controller.run(exec_test_function, ask=False, n_procs=15)


if __name__ == "__main__":
    testbed()


