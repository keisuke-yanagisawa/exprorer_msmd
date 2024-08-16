import collections.abc
import os
import random
import string
from pathlib import Path
from typing import Union

from .logger import logger


def update_dict(d: dict, u: dict):
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            d[k] = update_dict(d.get(k, {}), v)  # type: ignore
        else:
            d[k] = v
    return d


def getabsolutepath(path: Path) -> Path:
    """
    Get absolute path from relative path
    """
    path = expandpath(path)
    if not path.is_absolute():  # relative path
        path = os.getcwd() / path
    return path


def expandpath(path: Path) -> Path:
    """
    Expand ~ and $HOME and other environment variables
    """
    path = path.expanduser()
    return Path(os.path.expandvars(path))


def randomname(n: int):
    """
    Generate random string of length n
    """
    randlst = [random.choice(string.ascii_letters + string.digits) for i in range(n)]
    return "".join(randlst)


def dat_dump(dat):
    for section in dat.sections():
        logger.debug("[%s]" % section)
        for option in dat[section]:
            logger.debug("  %s : %s" % (option, dat[section][option]))


# Notation 1-2 => 1,2  5-9:2 => 5,7,9
def expand_index(ind_info: Union[str, int]):
    # single element will be parsed as integer in yaml.safe_load
    # Thus this process is needed
    if isinstance(ind_info, int):
        return [ind_info]

    ret = []
    elems = ind_info.split(",")
    for elem in elems:
        if elem.find("-") == -1:
            elem = int(elem)
            ret.append(elem)
        elif elem.find(":") == -1:
            st, ed = elem.split("-")
            st, ed = int(st), int(ed)
            ret.extend(list(range(st, ed + 1)))
        else:
            window, offset = elem.split(":")
            st, ed = window.split("-")
            st, ed, offset = int(st), int(ed), int(offset)
            ret.extend(list(range(st, ed + 1, offset)))
    return ret
