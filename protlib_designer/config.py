import pandas as pd
from pydantic import BaseModel
from yaml import Loader, load

from protlib_designer import logger
from protlib_designer.utils import (
    validate_data,
    validate_objective_constraints,
    validate_schedule_param,
)


class ProtlibDesignerConfig(BaseModel):
    data: str
    nb_iterations: int
    min_mut: int
    max_mut: int
    output_folder: str
    forbidden_aa: str
    max_arom_per_seq: int
    dissimilarity_tolerance: float
    interleave_mutant_order: bool
    force_mutant_order_balance: bool
    schedule: int
    schedule_param: str
    objective_constraints: str
    objective_constraints_param: str
    weighted_multi_objective: bool
    debug: int
    data_normalization: bool

    @staticmethod
    def get_default_config():
        raise NotImplementedError()

    @staticmethod
    def parse_config(config):
        """Format the input config into a ProtlibDesignerConfig object."""
        pass

    @staticmethod
    def read_config(config):
        if isinstance(config, str):
            with open(config, encoding="utf-8") as f:
                config = load(f, Loader=Loader)
        elif isinstance(config, dict):
            config = config
        else:
            raise TypeError(
                f"config parameter should be str or dict type, not {type(config)}"
            )

        return ProtlibDesignerConfig(**config)


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
    """Format and validate the parameters.

    Parameters
    ----------
    data : str
        The path to the data file.
    nb_iterations : int
        The number of iterations to run.
    forbidden_aa : str
        The forbidden amino acids.
    max_arom_per_seq : int
        The maximum number of aromatic amino acids per sequence.
    schedule : int
        The schedule to use.
    schedule_param : str
        The scheduling parameters.
    objective_constraints : str
        The objective constraints.
    objective_constraints_param : str
        The objective constraints parameters.
    """
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
