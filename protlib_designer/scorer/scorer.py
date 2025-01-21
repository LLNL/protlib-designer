from abc import ABC, abstractmethod
import numpy as np
from numpy.typing import NDArray
from typing import List


def score_function(
    position_logps: NDArray, wt_logps: NDArray, score_type: str = "llr"
) -> NDArray:
    """Calculate the score for the positions."""
    if score_type == "ll":
        position_score = position_logps
    elif score_type == "minus_ll":
        position_score = -position_logps
    elif score_type == "llr":
        position_score = position_logps - wt_logps
    elif score_type == "minus_llr":
        position_score = -position_logps + wt_logps
    elif score_type == "probs":
        position_score = np.exp(position_logps)
    return position_score


def from_user_input_to_scorer_input(positions: List[str], mapping: dict) -> List[str]:
    """Use the mapping to convert the positions to the numbering used by the scorer.

    Parameters
    ----------
    positions : list
        List of positions.
    mapping : dict
        mapping : dict
            A dictionary to map the positions from the user input to the numbering used by the scorer. Example:
            {'H': {'H1': 1, 'H5': 2, 'H6': 3}, 'L': {'L1': 1, 'L2': 2, 'L3': 3}}
    """
    for i, position in enumerate(positions):
        wt, chain, pos = position[0], position[1], position[2:]
        mapped_pos = mapping[chain][str(chain) + str(pos)]
        positions[i] = wt + chain + str(mapped_pos)
    return positions


def from_scorer_output_to_user_output(mutation: str, mapping: dict) -> str:
    """Use the mapping to convert the positions from the numbering used by the scorer to the user input.

    Parameters
    ----------
    mutations : list
        List of mutations.
    mapping : dict
        Mapping of chains and positions.
    """

    wt, chain, pos, aa = mutation[0], mutation[1], mutation[2:-1], mutation[-1]
    chain_pos = list(mapping[chain].keys())[
        list(mapping[chain].values()).index(int(pos))
    ]
    return wt + chain_pos + aa


class Scorer(ABC):
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def load_model(self):
        """Load the model."""
        pass

    @abstractmethod
    def prepare_input(self, **kwargs):
        """Prepare the input for the model."""
        pass

    @abstractmethod
    def forward_pass(self):
        """Perform the forward pass."""
        pass

    @abstractmethod
    def get_scores(self):
        """Get the scores."""
        pass

    @abstractmethod
    def __str__(self):
        return "Generator"
