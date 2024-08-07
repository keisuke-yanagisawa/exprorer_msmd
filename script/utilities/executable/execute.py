import subprocess

from ..logger import logger


class Command(object):
    def __init__(self, comm):
        self.comm = comm

    def __str__(self):
        return self.comm

    def run(self):
        try:
            res = subprocess.run(self.comm, shell=True, stdout=subprocess.PIPE, check=True)
        except Exception as e:
            logger.error(e)
            raise e
        return res.stdout.decode("utf-8")
