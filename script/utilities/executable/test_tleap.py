from unittest import TestCase
from script.utilities.executable.tleap import TLeap


class TestTLeap(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestTLeap, self).__init__(*args, **kwargs)
        self.tleap = TLeap()
