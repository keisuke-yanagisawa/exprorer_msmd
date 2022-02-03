# https://qiita.com/yopya/items/63155923602bf97dec53
from logging import Formatter, handlers, StreamHandler, getLogger, DEBUG

class Logger:
    def __init__(self, name=__name__, path=""):
        self.logger = getLogger(name)
        self.logger.setLevel(DEBUG)
        formatter = Formatter("[%(asctime)s] [%(process)d] [%(name)s] [%(levelname)s] %(message)s")

        # stdout
        handler = StreamHandler()
        handler.setLevel(DEBUG)
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)

        # file
        if(path != ""):
            handler = handlers.RotatingFileHandler(filename = path,
                                                   maxBytes = 1048576,
                                                   backupCount = 3)
            handler.setLevel(DEBUG)
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)

    def debug(self, msg):
        self.logger.debug(msg)

    def info(self, msg):
        self.logger.info(msg)

    def warn(self, msg):
        self.logger.warning(msg)

    def error(self, msg):
        self.logger.error(msg)

    def critical(self, msg):
        self.logger.critical(msg)