import subprocess

class Command(object):
  def __init__(self, comm):
    self.comm = comm

  def __str__(self):
    return self.comm
  
  def run(self):
    try:
      res = subprocess.run(self.comm, shell=True, stdout=subprocess.PIPE)
    except Exception as e:
      print(e)
    return res.stdout.decode("utf-8")