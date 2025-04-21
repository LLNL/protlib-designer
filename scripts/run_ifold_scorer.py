import click
import pandas as pd

from protlib_designer import logger
from protlib_designer.scorer.ifold_scorer import IFOLDScorer

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('pdb_file', type=str, required=True)
@click.argument('positions', type=str, required=True, nargs=-1)
@click.option('--models', 'model_names', type=str, multiple=True, required=False, default=["v_48_002"])
@click.option("--model-paths", "model_paths", type=str, default="/usr/workspace/vaccines/proteinmpnn_weights/vanilla_model_weights/")
@click.option(
    '--score-type',
    type=click.Choice(['minus_ll', 'minus_llr']),
    default='minus_llr',
)
@click.option('--output-file', type=str, default='plm_scores.csv')
def run_ifold_scorer(
    pdb_file,
    positions,
    model_names,
    model_paths,
    score_type,
    output_file
):

    dataframes = []

    for model_name in model_names:
        ifold_scorer = IFOLDScorer(
            model_name=model_name,
            model_path=model_paths,
            score_type=score_type
        )
        df = ifold_scorer.get_scores(pdb_file, list(positions))
        dataframes.append(df)

    if not dataframes:
        logger.error("No dataframes to combine.")

    # Merge the dataframes over the column "Mutation"
    combined_df = None
    for i, df in enumerate(dataframes):
        combined_df = df if i == 0 else pd.merge(combined_df, df, on="Mutation")

    combined_df.to_csv(output_file, index=False)
    logger.info(f"Combined scores saved to {output_file}")


if __name__ == "__main__":
    run_ifold_scorer()
