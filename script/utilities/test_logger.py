import os
import tempfile
import re
import pytest
from logging import DEBUG, INFO, WARN, ERROR, CRITICAL, getLogger, StreamHandler, handlers
from .logger import Logger, __DEFAULT_LOG_LEVEL__


@pytest.fixture(autouse=True)
def clean_handlers():
    """
    Fixture to clean up logger handlers before and after each test
    Automatically applied to all tests with autouse=True
    """
    logger = getLogger()
    # Setup: Clear existing handlers
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    yield
    
    # Cleanup: Clear handlers after test
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)


@pytest.fixture
def temp_dir():
    """Fixture providing a temporary directory for testing"""
    with tempfile.TemporaryDirectory() as temp_dir:
        yield temp_dir


@pytest.fixture
def log_file(temp_dir):
    """Fixture providing a log file path for testing"""
    return os.path.join(temp_dir, "test.log")


@pytest.fixture
def logger():
    """Fixture providing a logger with standard output only"""
    return Logger()


@pytest.fixture
def file_logger(log_file):
    """Fixture providing a logger with file output"""
    return Logger(log_file)


class TestLoggerInitialization:
    """Logger initialization test class"""

    def test_default_initialization(self, logger):
        """Default initialization check"""
        assert logger.logger.getEffectiveLevel() == __DEFAULT_LOG_LEVEL__
        assert isinstance(logger.handler, StreamHandler)
        # Check for standard output handler
        assert any(isinstance(h, StreamHandler) and not isinstance(h, handlers.RotatingFileHandler) 
                  for h in logger.logger.handlers)

    def test_file_initialization(self, file_logger, log_file):
        """File output logger initialization check"""
        handlers_list = file_logger.logger.handlers
        # Check for standard output handler
        assert any(isinstance(h, StreamHandler) and not isinstance(h, handlers.RotatingFileHandler) 
                  for h in handlers_list)
        # Check for file output handler
        assert any(isinstance(h, handlers.RotatingFileHandler) for h in handlers_list)
        assert os.path.exists(log_file)


class TestLoggerLevels:
    """Log level test class"""

    @pytest.mark.parametrize("level_name,level_const", [
        ("debug", DEBUG),
        ("info", INFO),
        ("warn", WARN),
        ("warning", WARN),
        ("error", ERROR),
        ("critical", CRITICAL),
        # Case insensitive
        ("DEBUG", DEBUG),
        ("INFO", INFO),
        ("WARN", WARN),
        ("ERROR", ERROR),
        ("CRITICAL", CRITICAL),
        # Allow surrounding spaces
        (" info ", INFO),
        ("  debug  ", DEBUG),
        ("\tdebug\n", DEBUG)  # Allow tabs and newlines
    ])
    def test_valid_log_levels(self, logger, level_name, level_const):
        """Valid log level setting test"""
        logger.setLevel(level_name)
        assert logger.logger.getEffectiveLevel() == level_const
        assert logger.handler.level == level_const

    @pytest.mark.parametrize("invalid_level,expected_error", [
        ("", ValueError),
        (None, TypeError),
        ("INVALID", ValueError),
        (123, TypeError),
        ("debugger", ValueError),  # Invalid string
        ("inf o", ValueError),  # Space within word
        ("debug.info", ValueError),  # Contains invalid character
        ("DEBUG_LEVEL", ValueError)  # Contains invalid character
    ])
    def test_invalid_log_levels(self, logger, invalid_level, expected_error):
        """Invalid log level setting test"""
        with pytest.raises(expected_error) as exc_info:
            logger.setLevel(invalid_level)
        
        if expected_error == ValueError:
            # Check that error message includes list of valid log levels
            assert "Valid log levels:" in str(exc_info.value)
        elif expected_error == TypeError:
            assert "Log level must be specified as a string" in str(exc_info.value)


class TestLoggerOutput:
    """Log output test class"""

    @pytest.mark.parametrize("log_method,expected_level", [
        ("debug", "DEBUG"),
        ("info", "INFO"),
        ("warn", "WARNING"),
        ("error", "ERROR"),
        ("critical", "CRITICAL")
    ])
    def test_log_methods(self, file_logger, log_file, log_method, expected_level):
        """Test output for each log method"""
        file_logger.setLevel("debug")
        test_message = f"Test {log_method} message"
        
        method = getattr(file_logger, log_method)
        method(test_message)
        
        with open(log_file, 'r') as f:
            content = f.read()
            assert test_message in content
            assert f"[{expected_level}]" in content

    def test_log_format(self, file_logger, log_file):
        """Log format verification test"""
        test_message = "Test format message"
        file_logger.setLevel("info")  # Set to INFO level
        file_logger.info(test_message)

        with open(log_file, 'r') as f:
            content = f.read().strip()
            # Modified timestamp pattern (includes comma)
            timestamp_pattern = r"\[\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3}\]"
            level_pattern = r"\[INFO\]"
            pattern = f"{timestamp_pattern} {level_pattern} {test_message}"
            match = re.search(pattern, content)
            assert match is not None, \
                f"Log format does not match.\nExpected: {pattern}\nActual: {content}"

    def test_file_rotation(self, file_logger, log_file):
        """File rotation verification test"""
        file_logger.setLevel("debug")
        large_msg = "x" * 1000
        rotation_count = 1100  # Equivalent to 1.1MB

        # Check for log rotation occurrence
        for i in range(rotation_count):
            file_logger.debug(f"{large_msg} - {i}")

        # Check backup files
        backup_files = [f"{log_file}.{i}" for i in range(1, 4)]
        assert any(os.path.exists(bf) for bf in backup_files), "Backup files were not created"
        
        # Check main log file size
        assert os.path.getsize(log_file) <= 1048576, "Main log file exceeds maximum size"


class TestLoggerHandlers:
    """Logger handler test class"""

    def test_handler_synchronization(self, logger):
        """Handler level synchronization test"""
        test_levels = ["debug", "info", "warn", "error", "critical"]
        
        for level in test_levels:
            logger.setLevel(level)
            # Check only logger.handler level
            assert logger.handler.level == logger.logger.level, \
                f"Handler is not synchronized at level {level}"

    def test_handler_formatter(self, logger):
        """Handler formatter verification test"""
        expected_format = "[%(asctime)s] [%(levelname)s] %(message)s"
        assert logger.handler.formatter._fmt == expected_format, \
            f"Format does not match. Expected: {expected_format}, Actual: {logger.handler.formatter._fmt}"