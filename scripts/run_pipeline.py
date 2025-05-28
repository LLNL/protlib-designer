#!/usr/bin/env python
import click
import pandas as pd
from pathlib import Path
from functools import reduce

from protlib_designer import logger
from protlib_designer.dataloader import DataLoader
from protlib_designer.filter.no_filter import NoFilter
from protlib_designer.generator.ilp_generator import ILPGenerator
from protlib_designer.scorer.plm_scorer import PLMScorer
from protlib_designer.scorer.ifold_scorer import IFOLDScorer
from protlib_designer.solution_manager import SolutionManager
from protlib_designer.solver.generate_and_remove_solver import GenerateAndRemoveSolver
from protlib_designer.utils import (
    format_and_validate_protlib_designer_parameters,
    extract_sequence_from_pdb,
    write_config,
    cif_to_pdb,
)

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


def format_and_validate_pipeline_parameters(
    sequence: str = None,
    pdb_path: str = None,
    positions: list = None,
    plm_model_names: list = None,
    plm_model_paths: list = None,
    ifold_model_name: str = None,
    ifold_model_path: str = None,
):
    """Validate the parameters for the Protlib Designer pipeline.
    Parameters
    ----------
    sequence : str, optional
        Protein sequence for PLM scoring.
    pdb_path : str, optional
        Path to PDB file for IFold scoring.
    positions : list, optional
        List of positions to score.
    plm_models : list, optional
        List of PLM model names.
    plm_model_paths : list, optional
        List of PLM model paths.
    ifold_model_name : str, optional
        IFold model name.
    ifold_model_path : str, optional
        IFold model path.
    """

    if not sequence and not pdb_path:
        raise ValueError("You must provide either a sequence or a PDB file path.")

    if pdb_path and not Path(pdb_path).exists():
        raise ValueError(f"The provided PDB path does not exist: {pdb_path}")

    if not positions or len(positions[0]) <= 1:
        raise ValueError(
            "You must provide positions correctly formatted {Wildtype}{Chain}{Position} or {*}{Chain}{*}."
        )

    if not sequence:
        if not pdb_path or not Path(pdb_path).exists():
            raise ValueError(
                f"There is no way to extract the sequence. You must provide a sequence directly or extract it from a PDB file. Provided sequence: {sequence}, PDB path: {pdb_path}."
            )

        pdb_chain = positions[0][1]
        logger.warning(
            f"No sequence provided. Extracting the chain ID {pdb_chain} from the first position {positions[0]}."
        )
        sequence = extract_sequence_from_pdb(pdb_path, pdb_chain)
        logger.warning(
            f"No sequence provided. Extracting sequence: {sequence} from PDB file {pdb_path} and chain {pdb_chain}."
        )

    if (
        len(positions) == 1
        and positions[0].startswith("*")
    ):
        logger.warning(
            f"A placeholder position was provided {positions}, generating positions for the entire sequence."
        )
        # if the last position is *, generate positions for the entire sequence
        if positions[0].endswith("*"):
            positions = [
                f"{sequence[i]}{positions[0][1]}{i+1}" for i in range(len(sequence))
            ]
        elif positions[0][2] == "{" and positions[0][-1] == "}":
            range_str = positions[0][3:-1]
            start, end = map(int, range_str.split("-"))
            positions = [
                f"{sequence[i]}{positions[0][1]}{i+1}" for i in range(start - 1, end)
            ]
        else:
            raise ValueError(
                f"Invalid position format: {positions[0]}. Expected format: {{*}}{{Chain}}{{*}} or {{*}}{{Chain}}{{start-end}}."
            )
        logger.warning(f"Inferred positions: {positions}")

    return sequence, positions


def append_plm_scores_to_dataframes(
    dataframes: list,
    sequence: str,
    positions: list,
    plm_model_names: list = None,
    plm_model_paths: list = None,
    chain_type: str = "heavy",
    score_type: str = "minus_llr",
    mask: bool = True,
    mapping: str = None,
):
    """
    Append PLM scores to the provided dataframes.
    """
    for model_name in plm_model_names:
        plm_scorer = PLMScorer(
            model_name=model_name,
            model_path=None,
            score_type=score_type,
            mask=True,
            mapping=None,
        )
        df = plm_scorer.get_scores(sequence, list(positions), chain_type)
        dataframes.append(df)

    for model_path in plm_model_paths:
        plm_scorer = PLMScorer(
            model_name=model_path,
            model_path=model_path,
            score_type=score_type,
            mask=mask,
            mapping=mapping,
        )
        df = plm_scorer.get_scores(sequence, list(positions), chain_type)
        dataframes.append(df)


def append_ifold_scores_to_dataframes(
    dataframes: list,
    pdb_path: str,
    positions: list,
    ifold_model_name: str = None,
    ifold_model_path: str = None,
    seed: int = None,
    score_type: str = "minus_llr",
):
    """
    Append IFold scores to the provided dataframes.
    """

    ifold_scorer = IFOLDScorer(
        seed=seed,
        model_name=ifold_model_name,
        model_path=ifold_model_path,
        score_type=score_type,
    )
    df = ifold_scorer.get_scores(pdb_path, list(positions))
    dataframes.append(df)


def combine_dataframes(dataframes: list):
    """
    Combine multiple dataframes on the 'Mutation' column.
    """

    for df in dataframes:
        if 'Mutation' not in df.columns:
            logger.error(
                "Data file must have at minimum the Mutation column and at least one Objective/Target."
            )
            return
        if len(df) != len(dataframes[0]):
            logger.error("Data files must have the same number of rows.")
            return

    return reduce(
        lambda left, right: pd.merge(left, right, on='Mutation', how='left'),
        dataframes,
    )


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument("positions", type=str, nargs=-1)
@click.option("--sequence", type=str, help="Protein sequence for PLM scoring")
@click.option("--pdb-path", type=str, help="Path to PDB file for IFold scoring")
@click.option("--plm-model-names", type=str, multiple=True, help="PLM model names")
@click.option("--plm-model-paths", type=str, multiple=True, help="PLM model paths")
@click.option("--ifold-model-name", type=str, default=None, help="IFold model name")
@click.option("--ifold-model-path", type=str, default=None, help="IFold model path")
@click.option(
    "--score-type",
    type=click.Choice(["minus_ll", "minus_llr"]),
    default="minus_llr",
    help="Score type for scoring",
)
@click.option(
    "--intermediate-output",
    type=str,
    default="combined_scores.csv",
    help="Path to save intermediate scores",
)
@click.option(
    "--plm-chain-type", type=str, default="heavy", help="Chain type for PLM scoring"
)
@click.option(
    "--plm-mask/--no-plm-mask", default=True, help="Whether to mask wildtype amino acid"
)
@click.option(
    "--plm-mapping", type=str, default=None, help="Mapping file for PLM scoring"
)
@click.option(
    "--ifold-seed", type=int, default=None, help="Random seed for IFold scoring"
)
@click.option(
    "--nb-iterations",
    default=10,
    type=int,
    help="Number of iterations for the designer",
)
@click.option("--min-mut", default=1, type=int, help="Minimum number of mutations")
@click.option("--max-mut", default=4, type=int, help="Maximum number of mutations")
@click.option(
    "--output-folder",
    default="lp_solution",
    type=click.Path(exists=False),
    help="Output folder for the designer results",
)
@click.option(
    "--forbidden-aa", type=str, help="Comma-separated list of forbidden amino acids"
)
@click.option(
    "--max-arom-per-seq",
    type=int,
    help="Maximum number of aromatic residues per sequence",
)
@click.option(
    "--dissimilarity-tolerance",
    default=0.0,
    type=float,
    help="Dissimilarity tolerance for the designer",
)
@click.option(
    "--interleave-mutant-order",
    default=False,
    type=bool,
    help="Interleave mutant order in the designer",
)
@click.option(
    "--force-mutant-order-balance",
    default=False,
    type=bool,
    help="Force balance in mutant order in the designer",
)
@click.option(
    "--schedule",
    default=0,
    type=int,
    help="Schedule type for the designer (0: no schedule, 1: remove commonest mutation/position every p0/p1 iterations, 2: remove mutation/position if it appears more than p0/p1 times)",
)
@click.option(
    "--schedule-param", type=str, help="Parameters for the schedule (e.g., 'p0,p1')"
)
@click.option(
    "--objective-constraints",
    type=str,
    help="Objective constraints for the designer (e.g., 'constraint1,constraint2')",
)
@click.option(
    "--objective-constraints-param",
    type=str,
    help="Parameters for the objective constraints (e.g., 'param1,param2')",
)
@click.option(
    "--weighted-multi-objective",
    default=True,
    type=bool,
    help="Use weighted multi-objective optimization in the designer",
)
@click.option("--debug", default=0, type=int, help="Debug level")
@click.option(
    "--data-normalization",
    default=False,
    type=bool,
    help="Normalize data before running the designer",
)
def run_pipeline(
    positions,
    sequence,
    pdb_path,
    plm_model_names,
    plm_model_paths,
    ifold_model_name,
    ifold_model_path,
    score_type,
    intermediate_output,
    plm_chain_type,
    plm_mask,
    plm_mapping,
    ifold_seed,
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
    """
    Run the complete Protlib Designer pipeline:\n
    1. Generate mutation scores using PLM and/or IFold
    2. Run the designer to create optimal protein libraries

    You must specify either a sequence (for PLM) or a PDB path (for IFold) or both.
    """

    logger.info("Starting the Protlib Designer pipeline...")

    sequence, positions = format_and_validate_pipeline_parameters(
        sequence,
        pdb_path,
        positions,
        plm_model_names,
        plm_model_paths,
        ifold_model_name,
        ifold_model_path,
    )

    dataframes = []

    # Step 1: Run PLM scoring if a plm model is provided
    if plm_model_names or plm_model_paths:
        logger.info("Running PLM Scorer...")
        logger.info(f"Sequence: {sequence}")
        logger.info(f"Positions: {positions}")

        append_plm_scores_to_dataframes(
            dataframes,
            sequence,
            positions,
            plm_model_names=plm_model_names,
            plm_model_paths=plm_model_paths,
            chain_type=plm_chain_type,
            score_type=score_type,
            mask=plm_mask,
            mapping=plm_mapping,
        )

        logger.info(f"PLM scoring completed with {len(dataframes)} models")

    # Step 2: Run IFold scoring if PDB path is provided
    if pdb_path:
        logger.info("Running IFOLD Scorer...")
        logger.info(f"PDB Path: {pdb_path}")
        logger.info(f"Positions: {positions}")

        # If the PDB path is provided, convert it to a PDB file if it's in CIF format
        if pdb_path.endswith(".cif"):
            logger.info("Converting CIF to PDB...")
            pdb_path = cif_to_pdb(pdb_path)
            logger.info(f"Converted CIF to PDB: {pdb_path}")

        append_ifold_scores_to_dataframes(
            dataframes,
            pdb_path,
            positions,
            ifold_model_name=ifold_model_name,
            ifold_model_path=ifold_model_path,
            seed=ifold_seed,
            score_type=score_type,
        )
        logger.info(f"IFold scoring completed with {len(dataframes)} models")

    if not dataframes:
        logger.error("No scores were generated. Check your input parameters.")
        return

    # Merge the dataframes into one, over the Mutation column
    logger.info("Combining scores from different models...")
    final_df = combine_dataframes(dataframes)
    logger.info("Scores combined successfully")

    # Save the combined scores
    final_df.to_csv(intermediate_output, index=False)
    logger.info(f"Combined scores saved to {intermediate_output}")

    # Step 3: Run the protlib designer
    logger.info("Running Protlib Designer...")

    # Format the input and validate parameters
    config, _ = format_and_validate_protlib_designer_parameters(
        output_folder,
        intermediate_output,
        min_mut,
        max_mut,
        nb_iterations,
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
    )

    # Load the data
    data_loader = DataLoader(intermediate_output)
    data_loader.load_data()
    config = data_loader.update_config_with_data(config)

    # Create the output directory
    output_path = Path(output_folder)
    if not output_path.exists():
        output_path.mkdir(parents=True, exist_ok=True)
    logger.info(f"Created directory {output_folder}")
    write_config(config, output_path)

    # Create the ILP generator
    ilp_generator = ILPGenerator(data_loader, config)

    # Create filter
    no_filter = NoFilter()

    # Create the solver
    generate_and_remove_solver = GenerateAndRemoveSolver(
        ilp_generator,
        no_filter,
        length_of_library=nb_iterations,
        maximum_number_of_iterations=2 * nb_iterations,
    )

    # Run the solver
    generate_and_remove_solver.run()

    # Process the solutions
    solution_manager = SolutionManager(generate_and_remove_solver)
    solution_manager.process_solutions()
    solution_manager.output_results()

    logger.info("Pipeline completed successfully!")


if __name__ == "__main__":
    run_pipeline()
