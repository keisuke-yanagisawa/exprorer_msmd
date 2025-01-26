# https://qiita.com/yopya/items/63155923602bf97dec53
from logging import CRITICAL, DEBUG, ERROR, INFO, WARN, Formatter, StreamHandler, getLogger, handlers

__DEFAULT_LOG_LEVEL__ = WARN


class Logger:
    def __init__(self, path=""):
        self.logger = getLogger()
        self.logger.setLevel(__DEFAULT_LOG_LEVEL__)
        formatter = Formatter("[%(asctime)s] [%(levelname)s] %(message)s")

        # stdout
        self.handler = StreamHandler()
        self.handler.setLevel(__DEFAULT_LOG_LEVEL__)
        self.handler.setFormatter(formatter)
        self.logger.addHandler(self.handler)

        # file
        if path != "":
            self.handler = handlers.RotatingFileHandler(filename=path, maxBytes=1048576, backupCount=3)
            self.handler.setLevel(__DEFAULT_LOG_LEVEL__)
            self.handler.setFormatter(formatter)
            self.logger.addHandler(self.handler)

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
        if not isinstance(level, str):
            raise TypeError(f"ログレベルは文字列で指定してください: {level}")

        level_lower = level.lower().strip()
        valid_levels = {
            "debug": DEBUG,
            "info": INFO,
            "warn": WARN,
            "warning": WARN,
            "error": ERROR,
            "critical": CRITICAL
        }

        if level_lower not in valid_levels:
            raise ValueError(
                f"無効なログレベルです: {level}\n"
                f"有効なログレベル: {', '.join(sorted(set(valid_levels.keys())))}"
            )

        log_level = valid_levels[level_lower]
        self.logger.setLevel(log_level)
        self.handler.setLevel(log_level)
        if log_level == DEBUG:
            self.debug("logging level: DEBUG")


logger = Logger()
