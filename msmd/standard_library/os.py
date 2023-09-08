import os


def getabsolutepath(path: str) -> str:
    """
    Get absolute path from relative path
    """
    path = __expandpath(path)
    if not path.startswith("/"):  # relative path
        path = os.getcwd() + "/" + path
    return path


def __expandpath(path: str) -> str:
    """
    Expand ~ and $HOME and other environment variables
    """
    path = os.path.expanduser(path)
    return os.path.expandvars(path)
