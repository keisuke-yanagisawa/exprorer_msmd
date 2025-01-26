import os
import tempfile
import re
import pytest
from logging import DEBUG, INFO, WARN, ERROR, CRITICAL, getLogger, StreamHandler, handlers
from .logger import Logger, __DEFAULT_LOG_LEVEL__


@pytest.fixture(autouse=True)
def clean_handlers():
    """
    各テストの前後でロガーのハンドラーをクリーンアップするfixture
    autouse=Trueで全テストに自動適用
    """
    logger = getLogger()
    # セットアップ: 既存のハンドラーをクリア
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    yield
    
    # クリーンアップ: テスト後にハンドラーをクリア
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)


@pytest.fixture
def temp_dir():
    """テスト用の一時ディレクトリを提供するfixture"""
    with tempfile.TemporaryDirectory() as temp_dir:
        yield temp_dir


@pytest.fixture
def log_file(temp_dir):
    """テスト用のログファイルパスを提供するfixture"""
    return os.path.join(temp_dir, "test.log")


@pytest.fixture
def logger():
    """標準出力のみのロガーを提供するfixture"""
    return Logger()


@pytest.fixture
def file_logger(log_file):
    """ファイル出力付きのロガーを提供するfixture"""
    return Logger(log_file)


class TestLoggerInitialization:
    """ロガーの初期化に関するテストクラス"""

    def test_default_initialization(self, logger):
        """デフォルト初期化の確認"""
        assert logger.logger.getEffectiveLevel() == __DEFAULT_LOG_LEVEL__
        assert isinstance(logger.handler, StreamHandler)
        # 標準出力用ハンドラーの存在確認
        assert any(isinstance(h, StreamHandler) and not isinstance(h, handlers.RotatingFileHandler) 
                  for h in logger.logger.handlers)

    def test_file_initialization(self, file_logger, log_file):
        """ファイル出力付きロガーの初期化確認"""
        handlers_list = file_logger.logger.handlers
        # 標準出力用ハンドラーの存在確認
        assert any(isinstance(h, StreamHandler) and not isinstance(h, handlers.RotatingFileHandler) 
                  for h in handlers_list)
        # ファイル出力用ハンドラーの存在確認
        assert any(isinstance(h, handlers.RotatingFileHandler) for h in handlers_list)
        assert os.path.exists(log_file)


class TestLoggerLevels:
    """ログレベルに関するテストクラス"""

    @pytest.mark.parametrize("level_name,level_const", [
        ("debug", DEBUG),
        ("info", INFO),
        ("warn", WARN),
        ("warning", WARN),
        ("error", ERROR),
        ("critical", CRITICAL),
        # 大文字小文字は区別しない
        ("DEBUG", DEBUG),
        ("INFO", INFO),
        ("WARN", WARN),
        ("ERROR", ERROR),
        ("CRITICAL", CRITICAL),
        # 前後のスペースは許容
        (" info ", INFO),
        ("  debug  ", DEBUG),
        ("\tdebug\n", DEBUG)  # タブや改行も許容
    ])
    def test_valid_log_levels(self, logger, level_name, level_const):
        """有効なログレベルの設定テスト"""
        logger.setLevel(level_name)
        assert logger.logger.getEffectiveLevel() == level_const
        assert logger.handler.level == level_const

    @pytest.mark.parametrize("invalid_level,expected_error", [
        ("", ValueError),
        (None, TypeError),
        ("INVALID", ValueError),
        (123, TypeError),
        ("debugger", ValueError),  # 無効な文字列
        ("inf o", ValueError),  # 単語内のスペース
        ("debug.info", ValueError),  # 無効な文字を含む
        ("DEBUG_LEVEL", ValueError)  # 無効な文字を含む
    ])
    def test_invalid_log_levels(self, logger, invalid_level, expected_error):
        """無効なログレベルの設定テスト"""
        with pytest.raises(expected_error) as exc_info:
            logger.setLevel(invalid_level)
        
        if expected_error == ValueError:
            # エラーメッセージに有効なログレベルの一覧が含まれていることを確認
            assert "有効なログレベル:" in str(exc_info.value)
        elif expected_error == TypeError:
            assert "ログレベルは文字列で指定してください" in str(exc_info.value)


class TestLoggerOutput:
    """ログ出力に関するテストクラス"""

    @pytest.mark.parametrize("log_method,expected_level", [
        ("debug", "DEBUG"),
        ("info", "INFO"),
        ("warn", "WARNING"),
        ("error", "ERROR"),
        ("critical", "CRITICAL")
    ])
    def test_log_methods(self, file_logger, log_file, log_method, expected_level):
        """各ログメソッドの出力テスト"""
        file_logger.setLevel("debug")
        test_message = f"Test {log_method} message"
        
        method = getattr(file_logger, log_method)
        method(test_message)
        
        with open(log_file, 'r') as f:
            content = f.read()
            assert test_message in content
            assert f"[{expected_level}]" in content

    def test_log_format(self, file_logger, log_file):
        """ログフォーマットの確認テスト"""
        test_message = "Test format message"
        file_logger.setLevel("info")  # INFOレベルに設定
        file_logger.info(test_message)

        with open(log_file, 'r') as f:
            content = f.read().strip()
            # タイムスタンプパターンを修正（カンマを含む）
            timestamp_pattern = r"\[\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3}\]"
            level_pattern = r"\[INFO\]"
            pattern = f"{timestamp_pattern} {level_pattern} {test_message}"
            match = re.search(pattern, content)
            assert match is not None, \
                f"ログフォーマットが一致しません。\n期待: {pattern}\n実際: {content}"

    def test_file_rotation(self, file_logger, log_file):
        """ファイルローテーションの確認テスト"""
        file_logger.setLevel("debug")
        large_msg = "x" * 1000
        rotation_count = 1100  # 1.1MB相当

        # ログローテーションの発生確認
        for i in range(rotation_count):
            file_logger.debug(f"{large_msg} - {i}")

        # バックアップファイルの確認
        backup_files = [f"{log_file}.{i}" for i in range(1, 4)]
        assert any(os.path.exists(bf) for bf in backup_files), "バックアップファイルが作成されていません"
        
        # メインログファイルのサイズ確認
        assert os.path.getsize(log_file) <= 1048576, "メインログファイルが最大サイズを超えています"


class TestLoggerHandlers:
    """ロガーのハンドラーに関するテストクラス"""

    def test_handler_synchronization(self, logger):
        """ハンドラーのレベル同期確認テスト"""
        test_levels = ["debug", "info", "warn", "error", "critical"]
        
        for level in test_levels:
            logger.setLevel(level)
            # logger.handlerのレベルのみを確認
            assert logger.handler.level == logger.logger.level, \
                f"レベル {level} でハンドラーが同期されていません"

    def test_handler_formatter(self, logger):
        """ハンドラーのフォーマッター確認テスト"""
        expected_format = "[%(asctime)s] [%(levelname)s] %(message)s"
        assert logger.handler.formatter._fmt == expected_format, \
            f"フォーマットが一致しません。期待: {expected_format}, 実際: {logger.handler.formatter._fmt}"