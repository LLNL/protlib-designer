import platform
from multiprocessing import cpu_count
from pathlib import Path

import click

from protlib_designer import logger
# from protlib_designer.config import ProtlibDesignerConfig
from protlib_designer.dataloader import DataLoader
from protlib_designer.filter.no_filter import NoFilter
from protlib_designer.generator.ilp_generator import ILPGenerator
from protlib_designer.solution_manager import SolutionManager
from protlib_designer.solver.generate_and_remove_solver import \
    GenerateAndRemoveSolver
from protlib_designer.utils import format_and_validate_parameters, write_config

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument("data", nargs=1, type=click.Path(exists=True), required=True)
@click.argument("nb_iterations", default=10, nargs=1, type=int)
@click.option("--min_mut", default=1, nargs=1, type=int)
@click.option("--max_mut", default=4, nargs=1, type=int)
@click.option(
    "--output-folder", default="lp_solution", nargs=1, type=click.Path(exists=False)
)
@click.option("--forbidden-aa", type=str)
@click.option("--max-arom-per-seq", type=int)
@click.option("--dissimilarity-tolerance", default=0.0, type=float)
@click.option("--interleave-mutant-order", default=False, type=bool)
@click.option("--force-mutant-order-balance", default=False, type=bool)
@click.option("--schedule", default=0, type=int)
@click.option("--schedule-param", type=str)
@click.option("--objective-constraints", type=str)
@click.option("--objective-constraints-param", type=str)
@click.option("--weighted-multi-objective", default=False, type=bool)
@click.option("--debug", default=0, type=int)
@click.option("--data-normalization", default=False, type=bool)
def run_protlib_designer(
    data: str,
    nb_iterations: int,
    min_mut: int,
    max_mut: int,
    output_folder: str,
    forbidden_aa: str,
    max_arom_per_seq: int,
    dissimilarity_tolerance: float,
    interleave_mutant_order: bool,
    force_mutant_order_balance: bool,
    schedule: int,
    schedule_param: str,
    objective_constraints: str,
    objective_constraints_param: str,
    weighted_multi_objective: bool,
    debug: int,
    data_normalization: bool,
):
    """Integer Linear Programming-based Optimization for Protein Library Design.

    \b
    Parameters
    ----------
    data : click.Path
        Path to the sum of single point mutations file
    nb_iterations : int
        Number of iterations
    min_mut : int
        Minimum number of mutations
    max_mut : int
        Maximum number of mutations
    output_folder : click.Path
        Output directory prefix
    forbidden_aa : str
        String of forbidden amino acid mutations at any position seperated by a comma.
        Pass as a comma seperated list (e.g C,K)
    max_arom_per_seq : int
        Maximum number of aromatic residues per sequence.
    dissimilarity_tolerance : float
        dissimilarity_tolerance parameter for diversity
    interleave_mutant_order : bool
        Make the method find first 1 mutation, then 2 mutations , ..., then max_mut mutations,
        then 1 mutation, then 2 mutations, ...
    force_mutant_order_balance : bool
        Force balance of number of mutations if interleave_mutant_order is True
    schedule : int
        Schedule for diversity 0, 1, or 2.
        0 : No schedule
        1 : Remove the commonest mutation every p0 iterations and remove the commonest position every p1 iterations
        2 : Remove the mutation if it appears more than p0 times and remove the position if it appears more than p1
        times
    schedule_param : str
        The interpretation of this parameter depends on the value of the schedule parameter:
        schedule=0 : No parameters
        schedule=1 : 'p0,p1' where p0 = Number of iterations to remove the commonest mutation and p1 = Number of iterations
        to remove the commonest position
        schedule=2 : 'p0,p1' where p0 = Number of occurrences of mutation to remove it and p1 = Number of occurrences of
        position to remove it
    objective_constraints : str
        Objective constraint. This is a list of objectives to constraint. It can be a list of names, e.g., "ddg1_stability,ddg2_binding".
        It can also be a list of indices, e.g., "1,2"
    objective_constraints_param : str
        Objective constraint parameters. This is a list of parameters for the objective constraints
    weighted_multi_objective : bool
        Use a weighted multi-objective formulation.
        If False, then reduce the multi-objective matrix using SVD and use the rank-1 approximation
        as the objective matrix
    debug : int
        The code will print debug information based on the value of this parameter.
        0 : No debug
        > 0 : Information about the ILP problem constraints is printed. The CPU time for each iteration is saved.
        > 1 : The ILP problems are saved to disk
        > 2 : The trace of the ILP solver is printed
    data_normalization : bool
        Normalize the data to be between 0 and 1
    """

    logger.info(f"System Type: {platform.platform()}")
    logger.info(f"Processors (logical cores): {cpu_count()}")
    logger.info(f"Python Version: {platform.python_version()}")

    # config = (
    #         ProtlibDesignerConfig.read_config(config) if config else ProtlibDesignerConfig.get_default_config()
    # )

    (
        _,
        nb_iterations,
        forbidden_aa,
        max_arom_per_seq,
        schedule,
        schedule_param,
        objective_constraints,
        objective_constraints_param,
    ) = format_and_validate_parameters(
        data,
        nb_iterations,
        forbidden_aa,
        max_arom_per_seq,
        schedule,
        schedule_param,
        objective_constraints,
        objective_constraints_param,
    )

    config = {
        "output_folder": output_folder,
        "data": data,
        "min_mut": min_mut,
        "max_mut": max_mut,
        "nb_iterations": nb_iterations,
        "forbidden_aa": forbidden_aa,
        "max_arom_per_seq": max_arom_per_seq,
        "dissimilarity_tolerance": dissimilarity_tolerance,
        "interleave_mutant_order": interleave_mutant_order,
        "force_mutant_order_balance": force_mutant_order_balance,
        "schedule": schedule,
        "schedule_param": schedule_param,
        "objective_constraints": objective_constraints,
        "objective_constraints_param": objective_constraints_param,
        "weighted_multi_objective": weighted_multi_objective,
        "debug": debug,
        "data_normalization": data_normalization,
    }

    data_loader = DataLoader(data)
    data_loader.load_data()
    config = data_loader.update_config_with_data(config)

    output_path = Path(output_folder)
    if not output_path.exists():
        output_path.mkdir(parents=True, exist_ok=True)
    logger.info(f"Created directory {output_folder}")
    write_config(config, output_path)

    # Create the ILP generator.
    ilp_generator = ILPGenerator(data_loader, config)

    # Create filter.
    no_filter = NoFilter()

    # Create the solver.
    generate_and_remove_solver = GenerateAndRemoveSolver(
        ilp_generator, no_filter, length_of_library=10, maximum_number_of_iterations=50
    )

    # Run the solver.
    generate_and_remove_solver.run()

    # Process the solutions.
    solution_manager = SolutionManager(generate_and_remove_solver)
    solution_manager.process_solutions()
    solution_manager.output_results()


if __name__ == "__main__":
    run_protlib_designer()
