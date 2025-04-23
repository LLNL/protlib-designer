import numpy as np


ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"


def assigned_fixed_chain(pdb_dict: dict, design_chain_list: list[str] = None) -> dict:
    all_chain_list = [
        item[-1:] for item in list(pdb_dict) if item[:9] == "seq_chain"
    ]  # ['A','B', 'C',...]

    if not design_chain_list:
        design_chain_list = ["A"]  # manually specify, e.g.

    fixed_chain_list = [
        letter for letter in all_chain_list if letter not in design_chain_list
    ]  # fix/do not redesign these chains

    return {pdb_dict["name"]: (design_chain_list, fixed_chain_list)}


def make_fixed_positions_dict(
    pdb_dict, position_list, design_chain_list, specify_non_fixed=True
):
    all_chain_list = [item[-1:] for item in list(pdb_dict) if item[:9] == "seq_chain"]
    fixed_position_dict = {}
    if not specify_non_fixed:
        for i, chain in enumerate(design_chain_list):
            fixed_position_dict[chain] = position_list[i]
        for chain in all_chain_list:
            if chain not in design_chain_list:
                fixed_position_dict[chain] = []
    else:
        assert len(position_list) == len(design_chain_list)
        for chain in all_chain_list:
            seq_length = len(pdb_dict[f"seq_chain_{chain}"])
            all_residue_list = (np.arange(seq_length) + 1).tolist()
            if chain not in design_chain_list:
                fixed_position_dict[chain] = all_residue_list
            else:
                idx = np.argwhere(np.array(design_chain_list) == chain)[0][0]
                fixed_position_dict[chain] = list(
                    set(all_residue_list) - set(position_list[idx])
                )
    return {pdb_dict["name"]: fixed_position_dict}