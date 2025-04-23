import click
import pandas as pd

from protlib_designer import logger
from protlib_designer.scorer.ifold_scorer import IFOLDScorer

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument("pdb_path", type=str, required=True)
@click.argument("positions", type=str, nargs=-1, required=True)
@click.option("--seed", type=int, default=None, help="Random seed for reproducibility.")
@click.option(
    "--model-name",
    "model_name",
    type=str,
    default=None,
    help="Model version names to use.",
)
@click.option(
    "--model-path",
    "model_path",
    type=str,
    default=None,
    help="Model weight paths corresponding to each model name.",
)
@click.option(
    "--score-type",
    type=click.Choice(["minus_ll", "minus_llr"]),
    default="minus_llr",
    help="Type of score to compute.",
)
@click.option(
    "--output-file",
    type=str,
    default="ifold_scores.csv",
    help="Output CSV file for combined scores.",
)
def run_ifold_scorer(
    pdb_path, positions, seed, model_name, model_path, score_type, output_file
):
    """
    Compute in silico mutagenesis scores using ProteinMPNN via IFOLDScorer.

    \b
    Parameters
    ----------
    pdb_path : str
        Path to the PDB file.
    positions : list[str]
        Positions to mutate, format WTCHAINPDBINDEX (e.g., EH1).
    seed : int
        Random seed for reproducibility.
    model_names : tuple[str]
        List of model version names.
    model_paths : tuple[str]
        List of model weight paths; must correspond one-to-one with model names.
    score_type : str
        Score type to compute: minus_ll or minus_llr.
    output_file : str
        Path to save the combined CSV of scores.
    """

    logger.info("Running IFOLD Scorer...")
    logger.info(f"PDB Path: {pdb_path}")
    logger.info(f"Positions: {positions}")

    scorer = IFOLDScorer(
        seed=seed, model_name=model_name, model_path=model_path, score_type=score_type
    )
    df = scorer.get_scores(pdb_path, list(positions))
    dataframes = [df]
    if not dataframes:
        logger.error("No dataframes to combine.")
        return

    # Merge all dataframes on the "Mutation" column
    combined_df = dataframes[0]
    for df in dataframes[1:]:
        combined_df = pd.merge(combined_df, df, on="Mutation")

    combined_df.to_csv(output_file, index=False)
    logger.info(f"Combined scores saved to {output_file}")

    logger.info("IFOLD Scorer run completed.")


if __name__ == "__main__":
    run_ifold_scorer()
