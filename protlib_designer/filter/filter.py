from abc import ABC, abstractmethod


class Filter(ABC):
    @abstractmethod
    def filter(self):
        pass

    @abstractmethod
    def __str__(self):
        pass
