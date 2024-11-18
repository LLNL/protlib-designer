import json
import string
import sys
from pathlib import Path

import pandas as pd

from lp_protein_design import logger

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
