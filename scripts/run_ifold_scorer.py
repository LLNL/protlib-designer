import click
import pandas as pd

from protlib_designer import logger
from protlib_designer.scorer.ifold_scorer import IFOLDScorer

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('pdb_path', type=str, required=True)
@click.argument('positions', type=str, nargs=-1, required=True)
@click.option('--seed', type=int, default=None, help="Random seed for reproducibility.")
@click.option(
    '--model-name',
    'model_names',
    type=str,
    multiple=True,
    required=True,
    help="Model version names to use.",
)
@click.option(
    '--model-path',
    'model_paths',
    type=str,
    multiple=True,
    required=True,
    help="Model weight paths corresponding to each model name.",
)
@click.option(
    '--score-type',
    type=click.Choice(['minus_ll', 'minus_llr']),
    default='minus_llr',
)
@click.option(
    '--output-file',
    type=str,
    default='ifold_scores.csv',
    help="Output CSV file for combined scores.",
)
def run_ifold_scorer(
    pdb_path, positions, seed, model_names, model_paths, score_type, output_file
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
    if len(model_names) != len(model_paths):
        logger.error(
            "The number of --model-name entries must match the number of --model-path entries."
        )
        return

    dataframes = []
    for name, path in zip(model_names, model_paths):
        scorer = IFOLDScorer(
            seed=seed, model_name=name, model_path=path, score_type=score_type
        )
        df = scorer.get_scores(pdb_path, list(positions))
        dataframes.append(df)

    if not dataframes:
        logger.error("No dataframes to combine.")
        return

    # Merge all dataframes on the "Mutation" column
    combined_df = dataframes[0]
    for df in dataframes[1:]:
        combined_df = pd.merge(combined_df, df, on="Mutation")

    combined_df.to_csv(output_file, index=False)
    logger.info(f"Combined scores saved to {output_file}")


if __name__ == "__main__":
    run_ifold_scorer()
