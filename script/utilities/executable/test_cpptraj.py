from unittest import TestCase
from script.utilities.executable.cpptraj import Cpptraj


class TestCpptraj(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestCpptraj, self).__init__(*args, **kwargs)
        self.cpptraj = Cpptraj()
