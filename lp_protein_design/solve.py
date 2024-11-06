# Built in libraries
import json
import platform
import string
import sys
import time
from multiprocessing import cpu_count
from pathlib import Path

# External Dependencies
import click
import numpy as np
import pandas as pd  # CSV reading
import pulp
from numpy.linalg import matrix_rank, svd

# Internal Dependencies
from lp_protein_design import logger

# TODO : Move amino acids list to utils/mutation.py
amino_acids = [
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
]

aromatic_amino_acids = [
    "F",
    "Y",
    "W",
]

BANNER = "########################################"

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
def ilp(
    data,
    nb_iterations,
    min_mut,
    max_mut,
    output_folder,
    forbidden_aa,
    max_arom_per_seq,
    dissimilarity_tolerance,
    interleave_mutant_order,
    force_mutant_order_balance,
    schedule,
    schedule_param,
    objective_constraints,
    objective_constraints_param,
    weighted_multi_objective,
    debug,
    data_normalization,
):
    """Integer Linear Programming-based Optimization for Antibody Design

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

    (
        data_df,
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

    targets = data_df.columns[1:].values.tolist()
    logger.info(f"Targets: {targets}")
    logger.info(f"Number of targets: {len(targets)}")

    positions, wildtype_position_amino = extract_positions_and_wildtype_amino_from_data(
        data_df
    )

    logger.info(f"Detected positions: {positions}")
    logger.info(f"Number of unique detected positions: {len(positions)}")
    logger.info(f"Detected wild type amino acid: {wildtype_position_amino}")

    # Check that max_mut is less than the number of positions
    if max_mut > len(positions) and interleave_mutant_order:
        logger.warning(
            f"Max number of mutations ({max_mut}) is greater than the number of positions ({len(positions)}). \
Setting max_mut to {len(positions)}."
        )
        max_mut = len(positions)

    given_path = Path(output_folder)
    if not given_path.exists():
        given_path.mkdir(parents=True, exist_ok=True)
    logger.info(f"Created directory {output_folder}")

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
    write_config(config, given_path)

    problem = pulp.LpProblem("GUIDE_Antibody_Optimization", pulp.LpMinimize)
    logger.info("Linear programming problem initialized")
    solver_msg = debug > 2

    # TODO : Add a solver parameter
    solver_list = pulp.listSolvers(onlyAvailable=True)
    logger.info(f"Available solvers: {solver_list}")
    # PULP_CBC_CMD : This solver uses a precompiled version of COIN-OR Branch and Cut solver (CBC) provided with the package.
    # https://coin-or.github.io/Cbc/intro.html
    solver = pulp.PULP_CBC_CMD(msg=solver_msg, keepFiles=False)

    # Building of the symbolic variable vector X and core objective matrix A.
    # Zero-pad the data with missing position-amino acid pairs.
    # The final length of the data is len(positions) * (len(amino_acids) - 1).
    # The added position-amino acid pairs are later constrained to be zero.
    # Loop over the positions and amino acids and check if the row exists.
    x_vars = []
    x_vars_dict = {}
    aromatic_vars = []  # Save the aromatic variables for later constraint
    zero_enforced_mutations = []
    data_df_padded = []  # The padded dataframe
    for position in positions:
        wt = wildtype_position_amino[position]
        for aa in amino_acids:
            if aa == wt:
                continue

            mutation_name = f"{wt}{position}{aa}"

            # Create the x variable
            x_var = pulp.LpVariable(f"X_{mutation_name}", cat="Binary")
            x_vars.append(x_var)
            x_vars_dict[mutation_name] = x_var

            # Save the aromatic variables
            if aa in aromatic_amino_acids:
                aromatic_vars.append(x_var)

            # [CONSTRAINT]: Force forbidden amino acids to be zero
            if aa in forbidden_aa:
                problem += x_var == 0, f"forbidden_{mutation_name}"
                if debug > 0:
                    logger.debug(
                        f"Adding constraint (forbidden amino acid): {x_var} == 0"
                    )

            # Check if row exists in the input dataframe
            if mutation_name in data_df["MutationHL"].values:
                # Extract the row from the dataframe in a dictionary format
                row = data_df[data_df["MutationHL"] == mutation_name].to_dict(
                    "records"
                )[0]
                # Append the row to the padded dataframe
                data_df_padded.append(row)
            else:  # The row does not exist in the input dataframe
                # Add 0-vector row for the new mutation
                new_row = {"MutationHL": mutation_name}
                # Save the position and aa to add X_pos_a = 0 constraint later in the script.
                zero_enforced_mutations.append((wt, position, aa))
                for curr_target in targets:
                    new_row[curr_target] = 0
                # Append the row to the padded dataframe
                data_df_padded.append(new_row)

                # [CONSTRAINT]: Force the missing position-amino acid pairs to be zero
                problem += x_var == 0, f"missing_{mutation_name}"
                if debug > 0:
                    logger.debug(
                        f"Adding constraint (mutation missing from data): {x_var} == 0"
                    )

    # Convert the padded dataframe to a pandas dataframe
    data_df = pd.DataFrame(data_df_padded)
    del data_df_padded

    logger.info(f"Number of x-variables: {len(x_vars)}")

    # Check the dimensions of the dataframe after adding the missing position-amino acid pairs.
    if data_df.shape[0] != len(positions) * (len(amino_acids) - 1):
        logger.error(
            f"Error adding missing position-amino acid pairs. Expected {len(positions) * (len(amino_acids) - 1)} rows. \
            Got {data_df.shape[0]} rows."
        )
        exit()

    # Check that data_df["MutationHL"].values is equivalent (ordered in the same way) as x_vars
    for index, x_var in enumerate(x_vars):
        mutation_name = x_var.getName().split("_")[1]
        if mutation_name != data_df["MutationHL"].values[index]:
            logger.error(
                f"Error adding missing position-amino acid pairs. Expected {mutation_name}. \
                Got {data_df['MutationHL'].values[index]}"
            )
            exit()

    data_df.to_csv(Path(given_path) / "padded_data.csv", index=False)

    objective_constraint_row_indices = []
    if len(objective_constraints) > 0:
        if objective_constraints[0].isdigit():
            objective_constraint_row_indices = [
                int(val) for val in objective_constraints
            ]
        else:
            cols = list(data_df.columns)
            for constraint_name in objective_constraints:
                # Subtract 1 because index 0 is the string "Mutation" not a target and will not be in the A matrix.
                try:
                    row_number = cols.index(constraint_name)
                    objective_constraint_row_indices.append(row_number)
                except ValueError:
                    logger.error(
                        f"Constraint '{constraint_name}' was not found in the input file."
                    )
                    exit()

    A_original = data_df[data_df.columns[1:]].values.transpose()

    if objective_constraint_row_indices is not None:
        constraint_data_df = data_df[data_df.columns[objective_constraint_row_indices]]
        data_df = data_df.drop(
            data_df.columns[objective_constraint_row_indices], axis=1
        )

    # Tranpose the data to get the matrix A
    A_auxiliar = data_df[data_df.columns[1:]].values.transpose()
    if data_normalization:
        min_vals = np.min(A_auxiliar, axis=1)
        max_vals = np.max(A_auxiliar, axis=1)
        # normalize the data to be between 0 and 1
        A = (A_auxiliar - min_vals[:, np.newaxis]) / (max_vals - min_vals)[
            :, np.newaxis
        ]
        # project to be between -1 and 1
        A = 2 * A - 1
        logger.info("Data was normalized to be between -1 and 1")
    else:
        A = A_auxiliar
    del A_auxiliar

    logger.info(f"Dimensions of raw data after transposing: {A.shape}")

    single_objective = len(targets) == 1
    if single_objective:
        assert A.shape[0] == 1
        if weighted_multi_objective:
            logger.warning(
                "Only one target was found. Setting weighted_multi_objective = False."
            )
            weighted_multi_objective = False

    # Get the low rank decomposition
    rank_of_matrix = matrix_rank(A)
    logger.info(f"The rank of the matrix is {rank_of_matrix}")

    # If we do not use the weighted multi-objective formulation, then we need to compute the SVD
    # to get the rank-1 approximation of the objective matrix.
    if not single_objective and not weighted_multi_objective:
        logger.info(
            "We have a multi-objective problem and not using the weighted multi-objective formulation.\n\
Computing the SVD to get the rank-1 approximation of the objective matrix."
        )
        try:
            u, sigma, vt = svd(A, full_matrices=False)
        except Exception:
            logger.error("Unable to compute the SVD")
            sys.exit(1)
        rr = np.mean(u[:, 0]) * sigma[0] * vt[0, :]

    # We have 20 amino acids minus one for the wild type
    NUM_AVAILABLE_AA = 19

    # [CONSTRAINT]: The number of mutations per position is at most 1
    position_to_x_vars_dict = {
        parse_mutation(x_vars[i].getName().split("_")[1])[1]: x_vars[
            i : i + NUM_AVAILABLE_AA
        ]
        for i in range(0, len(x_vars), NUM_AVAILABLE_AA)
    }
    for position, x_vars_in_pos in position_to_x_vars_dict.items():
        assert len(x_vars_in_pos) == 19
        problem += pulp.lpSum(x_vars_in_pos) <= 1, f"tot_mut_at_{position}"
        if debug > 0:
            logger.debug(
                f"Adding constraint (number of mutations per position <= 1): {pulp.lpSum(x_vars_in_pos)} <= 1"
            )

    # [CONSTRAINT]: The number of total mutations is at least min_mut and at most max_mut
    problem += pulp.lpSum(x_vars) >= min_mut, "min_mut"
    problem += pulp.lpSum(x_vars) <= max_mut, "max_mut"
    if debug > 0:
        logger.debug(
            f"Adding constraint (number of total mutations >= {min_mut}): {pulp.lpSum(x_vars)} >= {min_mut}"
        )
        logger.debug(
            f"Adding constraint (number of total mutations <= {max_mut}): {pulp.lpSum(x_vars)} <= {max_mut}"
        )

    # [CONSTRAINT]: Aromatic var count is less than user set value max_arom_per_seq
    if max_arom_per_seq is not None:
        problem += (
            pulp.lpSum(aromatic_vars) <= max_arom_per_seq,
            "max_arom_per_seq",
        )
        if debug > 0:
            logger.debug(
                f"Adding constraint (aromatic var count <= {max_arom_per_seq}): {pulp.lpSum(aromatic_vars)}\
<= {max_arom_per_seq}"
            )

    # [CONSTRAINT]: The objective constraints
    if objective_constraint_row_indices is not None:
        for index, constraint in enumerate(objective_constraint_row_indices):
            A_column = constraint_data_df.iloc[:, index].values
            problem += (
                pulp.lpSum(A_column * x_vars) <= objective_constraints_param[index],
                f"objective_constraint_{constraint}",
            )
            if debug > 0:
                logger.debug(
                    f"Adding constraint (objective constraint {constraint}): {pulp.lpSum(A_column * x_vars)}\
<= {objective_constraints_param[index]}"
                )

    # Static objective problem
    if single_objective:
        problem += pulp.lpSum(np.dot(A, x_vars))
    elif not weighted_multi_objective:
        problem += pulp.lpSum(rr * x_vars)

    # Iterative loop for diversity inducing optimization
    many_hot_encoded_solutions = []
    list_of_solution_dicts = []
    cpu_times = []
    positions_in_solution_counts = {}
    mutations_in_solution_counts = {}
    wallclock_time_start = time.time()
    logger.info(f"{BANNER}\nStarting {nb_iterations} iterations\n{BANNER}")
    for iteration in range(nb_iterations):

        logger.info(f"Iteration {iteration}")

        if interleave_mutant_order:
            curr_max_mut = min_mut + (iteration % (max_mut - min_mut + 1))
            problem.constraints["max_mut"].changeRHS(curr_max_mut)
            if debug > 0:
                logger.debug(
                    f"Changing constraint max_mut constraint to {curr_max_mut}"
                )
            if force_mutant_order_balance:
                problem.constraints["min_mut"].changeRHS(curr_max_mut)
                if debug > 0:
                    logger.debug(
                        f"Changing constraint min_mut constraint to {curr_max_mut}"
                    )

        if weighted_multi_objective:
            weights = np.random.random(size=A.shape[0])
            # normalize to add up to 1
            weights = weights / np.sum(weights)
            dict_weights = dict(zip(targets, weights))
            if debug > 0:
                logger.debug(f"Target weights: {dict_weights}")
            problem += pulp.lpSum(weights * np.dot(A, x_vars))

        if debug > 1:
            if iteration == 0:
                dir_problems = Path(given_path) / "ilp_problems"
                # Create the directory to save the problem files to
                if not dir_problems.exists():
                    dir_problems.mkdir(parents=True, exist_ok=True)
            # Save the problem to a file
            problem.writeLP(str(dir_problems / f"problem_{iteration}.lp"))

        # Solve the problem
        cpu_time_start = time.time()
        # return the number of function evaluations
        status = problem.solve(solver)

        if status != 1:
            logger.error(
                f"Error solving problem on iteration: {iteration}. Error Status: {pulp.LpStatus[status]}"
            )
            break
        cpu_time = time.time() - cpu_time_start
        cpu_times.append(cpu_time)

        # Post processing the solution
        many_hot_encoded_solution = [x.value() for x in x_vars]
        x_mutations_in_solution = [x for x in x_vars if x.value() == 1]
        parsed_mutations_in_solution = [
            parse_mutation(x.getName().split("_")[1]) for x in x_mutations_in_solution
        ]
        formatted_mutations_in_solution = ",".join(
            [f"{x[0]}{x[1]}{x[2]}" for x in parsed_mutations_in_solution]
        )
        # Check if the solution is valid
        for mutation in parsed_mutations_in_solution:
            if mutation in zero_enforced_mutations:
                logger.error(
                    f"Solution {formatted_mutations_in_solution} is not valid. Skipping."
                )
                continue

        # Logging
        many_hot_encoded_solutions.append(many_hot_encoded_solution)
        objective_values_all = np.dot(A, many_hot_encoded_solution)
        solution_dict = {
            "solution": formatted_mutations_in_solution,
            "objective_value": problem.objective.value(),
            "reason_selected": f"ilp_solve_iteration_{iteration}",
        }
        if weighted_multi_objective:
            solution_dict["weights"] = weights
        if data_normalization or objective_constraint_row_indices is not None:
            objective_values_all_original = np.dot(
                A_original, many_hot_encoded_solution
            )
        for target, obj_value in zip(targets, objective_values_all):
            solution_dict[target] = obj_value
        if data_normalization or objective_constraint_row_indices is not None:
            for target, obj_value in zip(targets, objective_values_all_original):
                solution_dict[f"{target}_original"] = obj_value
        list_of_solution_dicts.append(solution_dict)
        logger.info(
            f"Runtime: {cpu_time}, Solution: {x_mutations_in_solution}, Objective: {problem.objective.value()}"
        )

        ##################################################################################################
        # Prepare constraints for next iteration
        for wt, pos, mutation_aa in parsed_mutations_in_solution:
            if pos in positions_in_solution_counts:
                positions_in_solution_counts[pos] += 1
            else:
                positions_in_solution_counts[pos] = 1

            mutation = f"{wt}{pos}{mutation_aa}"
            if mutation in mutations_in_solution_counts:
                mutations_in_solution_counts[mutation] += 1
            else:
                mutations_in_solution_counts[mutation] = 1

        # [CONSTRAINT]: Remove from search space the ball of radius dissimilarity_tolerance around the solution
        # (x_vars - previous_sol)^2 >= 1 + dissimilarity_tolerance
        # x_vars^2 - 2 * x_vars.previous_sol + previous_sol^2 >= 1 + dissimilarity_tolerance
        # Sum(x_vars) - 2 * x_vars.previous_sol  >= 1 + dissimilarity_tolerance - Sum(previous_sol) (binary variables)
        dissimilarity_rhs = (
            1 + dissimilarity_tolerance - pulp.lpSum(many_hot_encoded_solution)
        )
        dissimilarity_lhs = pulp.lpSum(x_vars) - 2 * np.dot(
            x_vars, many_hot_encoded_solution
        )
        problem += (
            dissimilarity_lhs >= dissimilarity_rhs,
            f"dissimilarity_from_iter_{iteration}",
        )
        if debug > 0:
            logger.debug(
                f"Adding constraint (dissimilarity): {dissimilarity_lhs} >= {dissimilarity_rhs}"
            )

        if schedule == 0:
            continue

        if schedule == 1:

            if iteration % schedule_param[0] == 0:
                sorted_mutation_counts = sorted(
                    mutations_in_solution_counts.items(),
                    key=lambda x: x[1],
                    reverse=True,
                )
                mutation_to_eliminate = sorted_mutation_counts[0][0]

                # [CONSTRAINT]: Remove the commonest mutation every p0 iterations
                problem += (
                    x_vars_dict[mutation_to_eliminate] == 0,
                    f"rm_{mutation_to_eliminate}_from_iter_{iteration}",
                )
                if debug > 0:
                    logger.debug(
                        f"Adding constraint (eliminate mutation): {x_vars_dict[mutation_to_eliminate]} == 0"
                    )
                mutations_in_solution_counts[mutation_to_eliminate] = -1

            if iteration % schedule_param[1] == 0:
                sorted_position_counts = sorted(
                    positions_in_solution_counts.items(),
                    key=lambda x: x[1],
                    reverse=True,
                )
                position_to_eliminate = sorted_position_counts[0][0]

                # [CONSTRAINT]: Remove the commonest position every p1 iterations
                problem += (
                    pulp.lpSum(position_to_x_vars_dict[position_to_eliminate]) == 0,
                    f"rm_\
{position_to_eliminate}_from_iter_{iteration}",
                )
                if debug > 0:
                    logger.debug(
                        f"Adding constraint (eliminate position): \
{pulp.lpSum(position_to_x_vars_dict[position_to_eliminate])} == 0"
                    )
                positions_in_solution_counts[position_to_eliminate] = -1

        elif schedule == 2:

            for mutation, count in mutations_in_solution_counts.items():
                if count > schedule_param[0]:
                    # [CONSTRAINT]: Remove the mutation if it appears more than p0 times
                    problem += (
                        x_vars_dict[mutation] == 0,
                        f"rm_{mutation}_from_iter_{iteration}",
                    )
                    if debug > 0:
                        logger.debug(
                            f"Adding constraint (eliminate mutation): {x_vars_dict[mutation]} == 0"
                        )
                    mutations_in_solution_counts[mutation] = -1

            for position, count in positions_in_solution_counts.items():
                if count > schedule_param[1]:
                    # [CONSTRAINT]: Remove the position if it appears more than p1 times
                    problem += (
                        pulp.lpSum(position_to_x_vars_dict[position]) == 0,
                        f"rm_\
{position}_from_iter_{iteration}",
                    )
                    if debug > 0:
                        logger.debug(
                            f"Adding constraint (eliminate position): \
{pulp.lpSum(position_to_x_vars_dict[position])} == 0"
                        )
                    positions_in_solution_counts[position] = -1

    # Post processing
    if list_of_solution_dicts:
        df_solutions = pd.DataFrame(list_of_solution_dicts)

        df_solutions.rename(columns={"solution": "MutationHL"}).to_csv(
            Path(given_path) / "solutions.csv", index=False
        )

        sol_list_sorted = zip(
            df_solutions["solution"].values, df_solutions["objective_value"].values
        )
        logger.info("All Solutions: ")
        for solution, obj_value in sol_list_sorted:
            solution_str = f"{solution} \t Objective Value: {round(obj_value, 3)}"
            logger.info(f"\t {solution_str}")
    else:
        logger.warn("No feasible solutions found.")

    if cpu_times:
        if debug > 0:
            logger.info("CPU times:")
            for index, rt in enumerate(cpu_times):
                logger.info(f"\tIteration {index}: {rt} seconds.")
            if debug > 1:
                # create a dataframe [Iteration, CPU time]
                cpu_times_df = pd.DataFrame(
                    {"Iteration": list(range(nb_iterations)), "CPU Time": cpu_times}
                )
                cpu_times_df.to_csv(Path(given_path) / "cpu_times.csv", index=False)
        logger.info("CPU time stats (per iteration):")
        logger.info(f"\tMin: {round(np.min(cpu_times), 2)} Seconds")
        logger.info(f"\tMean: {round(np.mean(cpu_times), 2)} Seconds")
        logger.info(f"\tMax: {round(np.max(cpu_times), 2)} Seconds")
        logger.info(f"CPU time (seconds): {str(np.sum(cpu_times))}")
        logger.info(
            f"Wallclock time (seconds): {str(time.time() - wallclock_time_start)} "
        )


def format_and_validate_parameters(
    data,
    nb_iterations,
    forbidden_aa,
    max_arom_per_seq,
    schedule,
    schedule_param,
    objective_constraints,
    objective_constraints_param,
):
    """Format the parameters."""
    if nb_iterations > 5000:
        logger.warning("Maximum number of iterations is 5000. Setting to 5000.")
        nb_iterations = 5000

    if forbidden_aa is not None:
        forbidden_aa = [x.strip() for x in forbidden_aa.split(",")]
        logger.info(f"Forbidden amino acids: {forbidden_aa}")
    else:
        forbidden_aa = []

    if max_arom_per_seq is not None and max_arom_per_seq < 0:
        logger.error("max_arom_per_seq must be a positive integer.")
        exit()

    if schedule_param is not None:
        schedule_param = [int(val) for val in schedule_param.split(",")]
    else:
        schedule_param = []
    validate_schedule_param(schedule, schedule_param)

    if objective_constraints is not None:
        objective_constraints = objective_constraints.split(",")
    else:
        objective_constraints = []

    if objective_constraints_param is not None:
        objective_constraints_param = [
            float(val) for val in objective_constraints_param.split(",")
        ]
    else:
        objective_constraints_param = []
    validate_objective_constraints(objective_constraints, objective_constraints_param)

    data_df = pd.read_csv(data)
    validate_data(data_df)

    return (
        data_df,
        nb_iterations,
        forbidden_aa,
        max_arom_per_seq,
        schedule,
        schedule_param,
        objective_constraints,
        objective_constraints_param,
    )


def validate_data(df: pd.DataFrame):
    """Validate the data file. The data file must have the following columns:
    MutationHL, Target1, Target2, ..., TargetN

    Parameters
    ----------
    df : pd.DataFrame
        The dataframe containing the data.
    """
    if "MutationHL" not in df.columns:
        logger.error("Data file must have a MutationHL column.")
        sys.exit(2)

    if len(df.columns) < 2:
        logger.error(
            "Data file must have at minimum the MutationHL column and at least one Objective/Target."
        )
        sys.exit(3)


def validate_schedule_param(schedule: int, schedule_param: list):
    """Validate the scheduling parameters.

    Parameters
    ----------
    schedule : int
        The schedule to use.
    schedule_param : list
        The scheduling parameters.
    """
    if schedule == 0 and schedule_param:
        logger.error(
            "Scheduling = 0 needs no parameters. Please use 'schedule-param = none'."
        )
        exit()

    if schedule == 1 and len(schedule_param) != 2:
        logger.error(
            "Scheduling = 1 need 2 parameters in the form 'p0,p1'. p0 : Number of iterations to kill the commonest \
                mutation. p1 : Number of iterations to kill the commonest positions."
        )
        exit()

    if schedule == 2 and len(schedule_param) != 2:
        logger.error(
            "Scheduling = 2 need 2 parameters in the form 'p0,p1'. p0 : Number of occurrences of mutation to kill it. \
                p1 : Number of occurrences of position to kill it"
        )
        exit()


def validate_objective_constraints(
    objectives_constraints: list, objectives_constraints_param: list
):
    """Validate the objective constraints.

    Parameters
    ----------
    objectives_constraints : list
        The list of objective constraints.
    objectives_constraints_param : list
        The list of parameters for the objective constraints.
    """
    if len(objectives_constraints) != len(objectives_constraints_param):
        logger.error(
            "Number of objective constraint parameters does not match number of objective constraints"
        )
        exit()


def extract_positions_and_wildtype_amino_from_data(df: pd.DataFrame):
    """Extract the positions at which mutations occur
    and the wild type amino acid at that position.

    Parameters
    ----------
    df : pd.DataFrame
        The dataframe containing the data.
    """
    mutation_full = df["MutationHL"].values.tolist()
    positions = []  # Positions that have mutations
    wildtype_position_amino = {}  # Position to wild type amino acid mapping
    for mutation in mutation_full:

        wildtype_amino, position, _ = parse_mutation(mutation)

        positions.append(position)
        if (
            position in wildtype_position_amino
            and wildtype_position_amino[position] != wildtype_amino
        ):
            logger.error(
                f"Conflicting information: Wild type amino at position {position} \
                said to be {wildtype_position_amino[position]} and {wildtype_amino}"
            )
            exit()

        # Save the wild type amino acid at this position
        wildtype_position_amino[position] = wildtype_amino

    # Get distinct positions
    positions = list(set(positions))

    # Order the positions in ascending order
    # Consider positions like H28 < H100A
    # positions = sorted(positions)
    positions_df = pd.DataFrame.from_dict(
        {
            i: {
                "chain": list(position)[0],
                "pos": int(position[1:].rstrip(string.ascii_uppercase)),
                "pos_extra": position[1:].lstrip("0123456789"),
            }
            for i, position in enumerate(positions)
        },
        orient="index",
    )

    positions_df = positions_df.sort_values(
        by=["chain", "pos", "pos_extra"],
        ascending=[True, True, True],
    )

    # Get the order by merging the strings
    positions = [
        f"{row['chain']}{row['pos']}{row['pos_extra']}"
        for _, row in positions_df.iterrows()
    ]

    return positions, wildtype_position_amino


# Mikel Jan 25, 2024. Copied from campaign_tools/utils/mutation.py
def parse_mutation(mutation: str):
    """Parse wild-type, position, mutated amino acid information
    from compact mutation notation.

    Parameters
    ----------
    mutation : str

    Returns
    -------
    Tuple
        (WT, pos, mutation) as str values.

    """
    return mutation[0], mutation[1:-1], mutation[-1]


def extract_mutation_key(mutationbreak):
    chars = list(mutationbreak)  # Convert string to list of chars
    return f"{chars[1]}_{chars[2]}"


def write_config(config: dict, output_dir: Path):
    """Write the parameters used to execute the program to disk.

    Parameters
    ----------
    config : dict
        The parameters used to execute the program.
    output_dir : str
        The output directory to write the config file to.
    """
    filepath = output_dir / "config.json"

    with open(filepath, "w") as fp:
        json.dump(config, fp, indent=4)


if __name__ == "__main__":
    ilp()
