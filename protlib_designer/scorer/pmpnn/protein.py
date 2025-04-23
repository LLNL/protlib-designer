from dataclasses import dataclass

import torch
from protlib_designer.scorer.pmpnn.utils import tied_featurize, parse_PDB


@dataclass
class Protein:
    structure_feats: torch.Tensor
    sequence_feats: torch.Tensor
    mask: torch.Tensor
    lengths: torch.Tensor
    chain_M: torch.Tensor
    chain_encoding_all: torch.Tensor
    chain_list_list: list[list]
    visible_list_list: list[list]
    masked_list_list: list[list]
    masked_chain_length_list_list: list[list]
    chain_M_pos: torch.Tensor
    omit_AA_mask: torch.Tensor
    residue_idx: torch.Tensor
    dihedral_mask: torch.Tensor
    tied_pos_list_of_lists_list: list[list[list]]
    pssm_coef: torch.Tensor
    pssm_bias: torch.Tensor
    pssm_log_odds_all: torch.Tensor
    bias_by_res_all: torch.Tensor
    tied_beta: torch.Tensor
    device: torch.device | str
    name: str | None = None

    @classmethod
    def get_features(
        cls,
        protein_data: dict,
        device: torch.device | str = "cuda",
        chain_id_dict: dict | None = None,
        fixed_positions_dict: dict | None = None,
        omit_AA_dict: dict | None = None,
        tied_positions_dict: dict | None = None,
        pssm_dict: dict | None = None,
        bias_by_res_dict: dict | None = None,
        ca_only: dict | None = None,
    ):
        if type(device) is str and "cuda" in device.lower():
            if torch.cuda.is_available():
                device = torch.device(device)
            else:
                device = torch.device("cpu")

        if type(protein_data) is not list:
            protein_data = [protein_data]
        outputs = tied_featurize(
            protein_data,
            device,
            chain_id_dict,
            fixed_position_dict=fixed_positions_dict,
            omit_AA_dict=omit_AA_dict,
            tied_positions_dict=tied_positions_dict,
            pssm_dict=pssm_dict,
            bias_by_res_dict=bias_by_res_dict,
            ca_only=ca_only,
        )
        return cls(*outputs, device=device, name=protein_data[0]["name"])

    @classmethod
    def from_pdb(
        cls,
        pdb_path: str,
        device: torch.device | str = "cuda",
        chain_id_dict: dict | None = None,
        fixed_positions_dict: dict | None = None,
        omit_AA_dict: dict | None = None,
        tied_positions_dict: dict | None = None,
        pssm_dict: dict | None = None,
        bias_by_res_dict: dict | None = None,
        ca_only: dict | None = None,
    ):
        pdb_dict_list = parse_PDB(pdb_path)
        name = pdb_dict_list[0]["name"]

        if type(device) is str and "cuda" in device.lower():
            if torch.cuda.is_available():
                device = torch.device(device)
            else:
                device = torch.device("cpu")

        outputs = tied_featurize(
            pdb_dict_list,
            device,
            chain_id_dict,
            fixed_position_dict=fixed_positions_dict,
            omit_AA_dict=omit_AA_dict,
            tied_positions_dict=tied_positions_dict,
            pssm_dict=pssm_dict,
            bias_by_res_dict=bias_by_res_dict,
            ca_only=ca_only,
        )
        return cls(*outputs, device=device, name=name)
