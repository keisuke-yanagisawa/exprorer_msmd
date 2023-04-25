import subprocess
from unittest import TestCase
from script.utilities.executable.execute import Command


class TestExecute(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestExecute, self).__init__(*args, **kwargs)

    def test_execute(self):
        command = Command(comm="echo Hello World!")
        self.assertEqual(str(command), "echo Hello World!")
        ret = command.run()
        self.assertEqual(ret, "Hello World!\n")

    def test_execute_error(self):
        command = Command(comm="echo Hello World!; exit 1")
        self.assertEqual(str(command), "echo Hello World!; exit 1")
        with self.assertRaises(subprocess.CalledProcessError):
            command.run()
