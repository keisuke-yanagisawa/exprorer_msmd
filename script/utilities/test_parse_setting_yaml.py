from unittest import TestCase
from script.utilities.util import parse_yaml, expand_index


class TestParseNormalSettingYaml(TestCase):

    def __init__(self, *args, **kwargs):
        super(TestParseNormalSettingYaml, self).__init__(*args, **kwargs)
        self.yaml_path = 'script/utilities/test_data/normal_setting.yaml'

    def test_parse_yaml(self):
        """check normal yaml file can be parsed"""
        parse_yaml(self.yaml_path)

    def test_parse_settings(self):
        """check some fields in normal yaml file"""
        setting = parse_yaml(self.yaml_path)
        self.assertEqual(setting['general']['name'], 'normal_protocol')
        self.assertEqual(len(setting['exprorer_msmd']['sequence']), 12)

    def test_default_values(self):
        setting = parse_yaml(self.yaml_path)
        self.assertEqual(setting["general"]["multiprocessing"], True)


class TestParseSettingYamlErrorCases(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestParseSettingYamlErrorCases, self).__init__(*args, **kwargs)

    # this test case is not needed since the parser itself works well
    # the error will be raised when the setting is used
    #
    # def test_empty_yaml(self):
    #     with self.assertRaises(ValueError):
    #         parse_yaml('script/utilities/test_data/empty.yaml')

    def test_yaml_file_does_not_exist(self):
        with self.assertRaises(FileNotFoundError):
            parse_yaml('script/utilities/test_data/does_not_exist.yaml')

    def test_yaml_file_is_directory(self):
        with self.assertRaises(IsADirectoryError):
            parse_yaml('script/utilities/test_data')

    def test_yaml_file_is_not_yaml(self):
        with self.assertRaises(ValueError):
            parse_yaml('script/utilities/test_data/normal_setting.json')


class TestExpandIndex(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestExpandIndex, self).__init__(*args, **kwargs)

    def test_single_element(self):
        self.assertEqual(expand_index('1'), [1])

    def test_range_elements(self):
        self.assertEqual(expand_index('1-10'), list(range(1, 11)))
        self.assertEqual(expand_index('1-10:2'), list(range(1, 11, 2)))
        self.assertEqual(expand_index('100-110'), list(range(100, 111)))

    def test_commaseparated_elements(self):
        self.assertEqual(expand_index('1,2,3'), [1, 2, 3])
        self.assertEqual(expand_index("1,5-10"), [1, 5, 6, 7, 8, 9, 10])
