from abc import ABC, abstractmethod


class Generator(ABC):
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def update_generator_before_generation(self):
        pass

    @abstractmethod
    def generate_one_solution(self, **kwargs) -> dict:
        pass

    @abstractmethod
    def update_generator(self):
        pass

    @abstractmethod
    def __str__(self):
        return "Generator"
