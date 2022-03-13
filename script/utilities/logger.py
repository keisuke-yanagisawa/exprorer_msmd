# https://qiita.com/yopya/items/63155923602bf97dec53
from logging import Formatter, handlers, StreamHandler, getLogger
from logging import DEBUG, ERROR, WARN, INFO, CRITICAL

__DEFAULT_LOG_LEVEL__ = INFO

class Logger:
    def __init__(self, path=""):
        self.logger = getLogger()
        self.logger.setLevel(__DEFAULT_LOG_LEVEL__)
        formatter = Formatter("[%(asctime)s] [%(levelname)s] %(message)s")

        # stdout
        handler = StreamHandler()
        handler.setLevel(__DEFAULT_LOG_LEVEL__)
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)

        # file
        if(path != ""):
            handler = handlers.RotatingFileHandler(filename = path,
                                                   maxBytes = 1048576,
                                                   backupCount = 3)
            handler.setLevel(__DEFAULT_LOG_LEVEL__)
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

    def setLevel(self, level):
        if level.lower() == "debug":
            self.logger.setLevel(DEBUG)
        elif level.lower() == "info":
            self.logger.setLevel(INFO)
        elif level.lower() == "warn" \
            or level.lower() == "warning":
            self.logger.setLevel(WARN)
        elif level.lower() == "error":
            self.logger.setLevel(ERROR)
        elif level.lower() == "critical":
            self.logger.setLevel(CRITICAL)
        else:
            self.warn(f"level '{level}' is invalid.")
            self.warn(f"logger level has not been changed.")

logger = Logger()