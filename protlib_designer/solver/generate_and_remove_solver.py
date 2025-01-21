import time

from protlib_designer import logger


class GenerateAndRemoveSolver:
    def __init__(
        self, generator, filter, length_of_library, maximum_number_of_iterations
    ):
        self.generator = generator
        self.filter = filter
        self.length_of_library = length_of_library
        self.maximum_number_of_iterations = maximum_number_of_iterations
        self.list_of_solution_dicts = []
        self.wallclock_time_start = None

    def run(self):
        """Generate and remove solutions until the library is full or the maximum number of iterations is reached."""
        number_of_solutions = 0
        iteration = 0

        self.wallclock_time_start = time.time()
        while (
            number_of_solutions < self.length_of_library
            and iteration < self.maximum_number_of_iterations
        ):
            self.generator.update_generator_before_generation(iteration=iteration)
            solution_dict = self.generator.generate_one_solution(iteration=iteration)
            try:
                solution = solution_dict.get("solution")
            except AttributeError:
                logger.info(f"No solution found for iteration {iteration}. Exiting.")
                break
            if self.filter.filter(solution):
                self.list_of_solution_dicts.append(solution_dict)
                number_of_solutions += 1
            self.generator.update_generator(iteration=iteration)
            iteration += 1

    def __str__(self):
        return f"GenerateAndRemoveSolver(generator={self.generator}, filter={self.filter}, length_of_library={self.length_of_library}, maximum_number_of_iterations={self.maximum_number_of_iterations})"
