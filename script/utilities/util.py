import random
import string
from os.path import expanduser, expandvars


def expandpath(path):
    path = expanduser(path)
    return expandvars(path)


def randomname(n):
    randlst = [random.choice(string.ascii_letters + string.digits) for i in range(n)]
    return "".join(randlst)


def dat_dump(dat):
    for section in dat.sections():
        print("[%s]" % section)
        for option in dat[section]:
            print("  %s : %s" % (option, dat[section][option]))
