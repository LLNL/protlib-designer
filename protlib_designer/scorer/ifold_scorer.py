import torch
import pandas as pd
from typing import List

from transformers import AutoModelForMaskedLM, AutoTokenizer

from protlib_designer.utils import amino_acids
from protlib_designer.scorer.scorer import (
    score_function,
    from_user_input_to_scorer_input,
    from_scorer_output_to_user_output,
    Scorer,
)
from protlib_designer.scorer.pmpnn.protein import Protein, from_pdb_dir
from protlib_designer.scorer.pmpnn.runner import ProteinMPNNRunner
from protlib.designer.scorer.pmpnn.utils import parse_PDB
from protlib_designer.scorer.pmpnn.extra_utils import (
    ALPHABET,
    assigned_fixed_chain,
    make_fixed_positions_dict,
)


class IFOLDScorer(Scorer):
    def __init__(
        self,
        seed: int | None = None,
        model_name: str = "v_48_002",
        model_path: str = "/usr/workspace/vaccines/proteinmpnn_weights/vanilla_model_weights/",
        score_type: str = "minus_llr",
    ):
        self.seed = seed
        self.model_name = model_name
        self.model_path = model_path
        self.score_type = score_type

        self.model = self.load_model()

    def load_model(self):
        model = ProteinMPNNRunner(
            seed=self.seed,
            model_name=self.model_name,
            model_weights_path=self.model_path,
        )
        return model

    def prepare_input(self, pdb_path: str, positions: List[str]):
        """Prepare the input for the model.

        Parameters
        ----------

        pdb_path : str
            Path to PDB structure file to read from.
        positions : list
            Positions on the sequence to be used to generate the score.
            Positions must be in the following format: {WT}{CHAIN}{PDBINDEX}.
            Note: PDBINDEX is 1-indexed, that is, the first position is 1. For example, the first positions in
            the list of positions are [EH1, VH2, QH3, ...].        
        chain_type : str
            Required parameter which specifies the chain type being analyzed.

        Returns
        -------

        proteins : list
            This returns a list of `Protein` inputs for ProteinMPNN.
        chains : list
            List of chains used, should all be same chain.
        locs : list
            List of positions to score
        wildtype_dict : dict
            Dictionary mapping {position: wildtype residue}
        """
        chains = [position[1] for position in positions]
        if len(set(chains)) > 1:
            raise ValueError(
                "All positions must have the same chain letter. Please provide positions with the same chain type."
            )
        locs = [int(position[2:]) for position in positions]
        wildtype_dict = {int(position[2:]): position[0] for position in positions}
        pdb_dict = parse_PDB(pdb_path)
        proteins = []
        for i, chain in enumerate(chains):
            fixed_chains = assigned_fixed_chain(pdb_dict, design_chain_list=[chain])
            fixed_positions = make_fixed_positions_dict(pdb_dict, locs[i], [chain], specify_non_fixed=True)
            protein = Protein.from_pdb(pdb_path, chain_id_dict=fixed_chains, fixed_positions_dict=fixed_positions)
            proteins.append(protein)
        return proteins, chains, locs, wildtype_dict

    def forward_pass(self, proteins: List[Protein]):
        log_prob_list = []
        for protein in proteins:
            log_probs, S, mask, design_mask, chain_order = self.model.get_probabilities(protein, conditional_probs=True)
            mask = mask.bool().squeeze()
            log_probs = log_probs.contiguous().view(-1, log_probs.size(-1))
            log_probs = log_probs[mask].cpu().double()
            log_prob_list.append(log_probs)
        return log_prob_list

    def get_scores(self, pdb_path: str, positions: list[str]):
        aas = list(ALPHABET)
        proteins, chains, locs, wildtype_dict = self.prepare_input(pdb_path, positions)
        logps = self.forward_pass(proteins)
        mutation2score = {}
        for bi, posn in enumerate(locs):
            wildtype_aa = wildtype_dict[posn]
            seq_index = posn - 1 # need to do 0 based index
            position_logps = logps[bi][seq_index].numpy()
            wt_aa_id = aas.index(wildtype_aa)
            wt_logps = logps[bi][wt_aa_id].numpy()
            position_scores = list(
                score_function(position_logps, wt_logps, score_type=self.score_type)
            )
            for i, amino_acid in enumerate(aas):
                mutation = f"{wildtype_aa}{chains[0]}{posn}{amino_acid}"
                mutation2score[mutation] = position_scores[i]

        return pd.DataFrame(
            mutation2score.items(),
            columns=["Mutation", f"ProteinMPNN_{self.score_type}"]
        )

    def __str__(self):
        return super().__str__() + ": IFOLD Scorer"
