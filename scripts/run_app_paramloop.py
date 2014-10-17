from hgext.mq import prev
from numpy.dual import inv
from re import sub, findall


class ParameterSet:

    def __init__(self,
                 start,
                 stop,
                 increment,
                 config_filename,
                 variable_pattern_in_config,
                 increment_rule=lambda previous, increment: previous + increment,
                 regex_flags=[]):

        self.start = start
        self.stop = stop
        self.increment = increment

        self.config_filename = config_filename
        self.variable_pattern_in_config = variable_pattern_in_config

        self.increment_rule = increment_rule
        self.regex_flags = regex_flags

        self.parent = None
        self.current_value = None

    def __len__(self):
        n = 0

        value = self.start

        while value <= self.stop:
            value = self.increment_rule(value, self.increment)
            n += 1

        return n

    def reset(self):

        self.write_value(self.start)

    def finished(self):

        return self.current_value > self.stop

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

    def write_value(self, value, prev_config_raw_text=None):

        self.current_value = value

        if not prev_config_raw_text:

            with open(self.config_filename, 'r') as config_file:

                prev_config_raw_text = config_file.read()

        new_config_raw_text = sub(self.invert_regex(self.variable_pattern_in_config),
                                  r"%s%s%s" % (r"\g<1>", value, r"\g<2>"),
                                  prev_config_raw_text,
                                  *self.regex_flags)

        with open(self.config_filename, 'w') as config_file:

            config_file.write(new_config_raw_text)

    def increment_value(self):

        with open(self.config_filename, 'r') as config_file:

            prev_config_raw_text = config_file.read()

        #try...except...
        prev_value_string = findall(self.variable_pattern_in_config, prev_config_raw_text)[0]

        self.write_value(self.increment_rule(self.current_value, self.increment), prev_config_raw_text)

    def iterate(self):

        self.increment_value()

        if self.parent:
            if self.finished():
                self.reset()
                self.parent.iterate()


class ParameterSetController:

    def __init__(self):
        self.parameter_sets = []

    def register_parameter_set(self, parameter_set):

        if len(self.parameter_sets) != 0:
            parameter_set.parent = self.parameter_sets[-1]

        self.parameter_sets.append(parameter_set)

    def prepare(self):

        for parameter_set in self.parameter_sets:
            parameter_set.prepare()

    def run(self, execute_program_rule, *args, **kwargs):

        n = 1

        for parameter_set in self.parameter_sets:
            n *= len(parameter_set)
            parameter_set.reset()

        ans = raw_input("Press either: \n\tenter to run %d simulations."
                        "\n\tt [time per sim in seconds] for a time estimate. "
                        "\n\tq to end.\n" % n)

        if ans.startswith("t "):
            time = float(ans.lstrip("t "))
            t = n

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

        while not self.parameter_sets[0].finished():

            execute_program_rule(*args, **kwargs)

            self.parameter_sets[-1].iterate()

    def __str__(self):

        s = ""

        for parameter_set in self.parameter_sets:

            s += " " + str(parameter_set.current_value)

        return s


def exec_test_function():

    testfile_name = "/tmp/test_paramloop.cfg"

    with open(testfile_name, 'r') as testfile:

        testfile_raw_text = testfile.read()

    print testfile_raw_text.strip("\n"), " - ",

    results = findall(r"a\s*\=\s*\[\s*(\d+)\s*\,\s*(.*)\s*\,\s*(\d+)\s*\]", testfile_raw_text)[0]

    for result in results:
        print result,
    print


def testbed():

    testfile_name = "/tmp/test_paramloop.cfg"

    with open(testfile_name, 'w') as testfile:

        testfile.write("a=[0,1,2]\n")

    set1 = ParameterSet( 0,   3,   1, testfile_name, r"a\=\[(\d+)\,.*\,\d+\]")
    set2 = ParameterSet(-2,   2, 1.5, testfile_name, r"a\=\[\d+\,(.*)\,\d+\]")
    set3 = ParameterSet( 0, 100,  30, testfile_name, r"a\=\[\d+\,.*\,(\d+)\]")

    controller = ParameterSetController()

    controller.register_parameter_set(set1)
    controller.register_parameter_set(set2)
    controller.register_parameter_set(set3)

    controller.run(exec_test_function)


if __name__ == "__main__":
    testbed()


