from protlib_designer.filter.filter import Filter


class NoFilter(Filter):
    def __init__(self):
        super().__init__()

    def filter(self, solution):
        return True

    def __str__(self):
        return "NoFilter()"
