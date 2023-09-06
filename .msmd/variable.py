import abc


class VariableInterface(abc.ABC):
    @staticmethod
    @abc.abstractmethod
    def _validation(var):
        pass

    @abc.abstractmethod
    def get(self):
        pass
