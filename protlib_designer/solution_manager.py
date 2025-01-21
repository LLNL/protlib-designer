import time
from pathlib import Path

import pandas as pd

from protlib_designer import logger


class SolutionManager:
    def __init__(self, solver):
        self.solver = solver
        self.solutions_df = None

    def process_solutions(self):
        solutions = self.solver.list_of_solution_dicts
        self.solutions_df = pd.DataFrame(solutions)
        if not self.solutions_df.empty:
            logger.info("All Solutions:")
            for idx, row in self.solutions_df.iterrows():
                logger.info(
                    f"Solution {idx}: {row['solution']}, Objective: {row['objective_value']}"
                )

    def output_results(self):
        config = self.solver.generator.config
        output_folder = config.get("output_folder", ".")
        output_path = Path(output_folder)
        if not output_path.exists():
            output_path.mkdir(parents=True, exist_ok=True)
        # Rename `solution` to `Mutation` for consistency with other outputs.
        self.solutions_df.rename(columns={"solution": "Mutation"}, inplace=True)
        self.solutions_df.to_csv(output_path / "solutions.csv", index=False)
        logger.info(f"Solutions saved to {output_path / 'solutions.csv'}")

        cpu_times = self.solver.generator.cpu_times
        if cpu_times and config.get("debug", 0) > 0:
            cpu_times_df = pd.DataFrame(
                {"Iteration": range(len(cpu_times)), "CPU Time": cpu_times}
            )
            cpu_times_df.to_csv(output_path / "cpu_times.csv", index=False)
            logger.info("CPU times saved")

        total_cpu_time = sum(cpu_times)
        total_wallclock_time = time.time() - self.solver.wallclock_time_start
        logger.info(f"Total CPU Time: {total_cpu_time}")
        logger.info(f"Total Wallclock Time: {total_wallclock_time}")
