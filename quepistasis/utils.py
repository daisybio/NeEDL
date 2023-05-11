import datetime
import json
import numpy as np


def get_ising_cost(h, J, variable_assignments):
    """

    :param h:
    :param J:
    :param variable_assignments:
    :return:
    """
    cost = 0
    for x, hcoeff in zip(variable_assignments, h):
        cost += x * hcoeff
    for ((s,e), jcoeff) in J.items():
        cost += variable_assignments[s] * variable_assignments[e] * jcoeff
    return cost


class MyJsonEncoder(json.JSONEncoder):
    """The class is used to transform NumPy objects and datetime objects into a format that can be fed to json"""

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, datetime.timedelta):
            return str(obj)
        return super(MyJsonEncoder, self).default(obj)


class Trace:
    """Trace the execution of the quantum annealing"""

    def __init__(self, path, print=False) -> None:
        """Constructor.
        :param path: relative or absolute path to save the files, including the first part of the file name.
            Passing './dwave' results in the creation of two files dwave_<timestamp>.json and dwave_<timestamp>.txt
            in the directory of the current execution.
        :param print: true if the tool is allowed to print on stdout.
        """
        self.info = {}
        self.info_text = ""
        self.print = print
        self.path = path

    def add(self, section, subsection, key, value):
        """Add new information
        :param section: primary tag of the information (str)
        :param subsection: secondary tag of the information, optional (str or None)
        :param key: name of the information (str)
        :param value: content of the information (obj)
        :return None
        """
        if type(value) == np.ndarray:
            value = value.tolist()

        # save for json
        if section not in self.info:
            self.info[section] = {}
        if subsection is not None:
            if subsection not in self.info[section]:
                self.info[section][subsection] = {}
            self.info[section][subsection][key] = value
        else:
            self.info[section][key] = value

        # save for text
        if subsection is not None:
            this_text = f"{section:40s} :: {subsection:30s} :: {key:40s} = {value}\n"
        else:
            this_text = f"{section:40s} :: {'':30s} :: {key:40s} = {value}\n"
        self.info_text += this_text
        if self.print:
            print(this_text, end='', flush=True)

    def save(self):
        """Save to json and txt"""
        json.dump(self.info, open(f"{self.path}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S_%f')}.json", "w"),
                  cls=MyJsonEncoder)
        open(f"{self.path}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S_%f')}.txt", "w").writelines(self.info_text)
