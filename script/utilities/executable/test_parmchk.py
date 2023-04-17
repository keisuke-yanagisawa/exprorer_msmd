from unittest import TestCase
from utilities.executable.parmchk import Parmchk


class TestParmchk(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestParmchk, self).__init__(*args, **kwargs)
        self.parmchk = Parmchk()
