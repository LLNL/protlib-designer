import click
import pandas as pd

from protlib_designer import logger
from protlib_designer.scorer.plm_scorer import PLMScorer

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('sequence', type=str, required=True)
@click.argument('positions', type=str, required=True, nargs=-1)
@click.option('--models', 'model_names', type=str, multiple=True, required=True)
@click.option('--chain-type', type=str, required=False, default="heavy")
@click.option('--model-paths', 'model_paths', type=str, multiple=True)
@click.option(
    '--score-type',
    type=click.Choice(['minus_ll', 'minus_llr']),
    default='minus_llr',
)
@click.option('--mask/--no-mask', default=True)
@click.option('--mapping', type=str, default=None)
@click.option('--output-file', type=str, default='plm_scores.csv')
def run_plm_scorer(
    sequence,
    positions,
    model_names,
    chain_type,
    model_paths,
    score_type,
    mask,
    mapping,
    output_file,
):
    """
    Compute in silico mutagenesis scores using Protein Language Models (PLM).

    \b
    Parameters
    ----------
    sequence : str
        The protein sequence.
    positions : str
        The positions to mutate.
    model_names : str
        The model names.
    chain_type : str
        The chain type.
    model_paths : str
        The model paths.
    score_type : str
        The score type.
    mask : bool
        Whether to mask the wild type amino acid.
    mapping : str
        The mapping.
    """

    logger.info("Running PLM Scorer...")
    logger.info(f"Sequence: {sequence}")
    logger.info(f"Positions: {positions}")

    dataframes = []

    for model_name in model_names:
        plm_scorer = PLMScorer(
            model_name=model_name,
            model_path=None,
            score_type=score_type,
            mask=mask,
            mapping=mapping,
        )
        df = plm_scorer.get_scores(sequence, list(positions), chain_type)
        dataframes.append(df)

    for model_path in model_paths:
        plm_scorer = PLMScorer(
            model_name=model_path,
            model_path=model_path,
            score_type=score_type,
            mask=mask,
            mapping=mapping,
        )
        df = plm_scorer.get_scores(sequence, list(positions), chain_type)
        dataframes.append(df)

    if not dataframes:
        logger.error("No dataframes to combine.")
        return

    # Merge the dataframes over the column "Mutation"
    combined_df = None
    for i, df in enumerate(dataframes):
        combined_df = df if i == 0 else pd.merge(combined_df, df, on="Mutation")

    combined_df.to_csv(output_file, index=False)
    logger.info(f"Combined scores saved to {output_file}")

    logger.info("PLM Scorer run completed.")


if __name__ == "__main__":
    run_plm_scorer()
