import string

# import Dict
from typing import Any, Dict

import pandas as pd

from protlib_designer import logger
from protlib_designer.utils import parse_mutation


def extract_positions_and_wildtype_amino_from_data(df: pd.DataFrame):
    """Extract the positions at which mutations occur
    and the wild type amino acid at that position.

    Parameters
    ----------
    df : pd.DataFrame
        The dataframe containing the data.
    """
    mutation_full = df["Mutation"].values.tolist()
    positions = []  # Positions that have mutations.
    wildtype_position_amino = {}  # Position to wild type amino acid mapping.
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

        # Save the wild type amino acid at this position.
        wildtype_position_amino[position] = wildtype_amino

    # Get distinct positions.
    positions = list(set(positions))

    # Order the positions in ascending order.
    # Consider positions like H28 < H100A.
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

    # Get the order by merging the strings.
    positions = [
        f"{row['chain']}{row['pos']}{row['pos_extra']}"
        for _, row in positions_df.iterrows()
    ]

    return positions, wildtype_position_amino


class DataLoader:
    def __init__(self, data_path):
        self.data_path = data_path
        self.data_df = None
        self.targets = None
        self.positions = None
        self.wildtype_position_amino = None

    def load_data(self):
        self.data_df = pd.read_csv(self.data_path)
        self.targets = self.data_df.columns[1:].values.tolist()
        (
            self.positions,
            self.wildtype_position_amino,
        ) = extract_positions_and_wildtype_amino_from_data(self.data_df)
        logger.info(f"Targets: {self.targets}")
        logger.info(f"Number of targets: {len(self.targets)}")
        logger.info(f"Detected positions: {self.positions}")
        logger.info(f"Number of unique detected positions: {len(self.positions)}")
        logger.info(f"Detected wild type amino acid: {self.wildtype_position_amino}")

    def update_config_with_data(self, config: Dict[str, Any]):
        # Check that max_mut is less than the number of positions.
        if (
            config["max_mut"] > len(self.positions)
            and config["interleave_mutant_order"]
        ):
            logger.warning(
                f"Max number of mutations ({config['max_mut']}) is greater than the number of positions ({len(self.positions)}). \
Setting max_mut to {len(self.positions)}."
            )
            config["max_mut"] = len(self.positions)

        return config
