import json
import glob
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


# Not used in the code
def parse_multiple_chains(
    input_path=None, input_list=None, save_path=None, ca_only=False
):
    alpha_1 = list("ARNDCQEGHILKMFPSTWYV-")
    states = len(alpha_1)
    alpha_3 = [
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "GAP",
    ]

    aa_1_N = {a: n for n, a in enumerate(alpha_1)}
    aa_3_N = {a: n for n, a in enumerate(alpha_3)}
    aa_N_1 = dict(enumerate(alpha_1))
    aa_1_3 = dict(zip(alpha_1, alpha_3))
    aa_3_1 = {b: a for a, b in zip(alpha_1, alpha_3)}

    def AA_to_N(x):
        # ["ARND"] -> [[0,1,2,3]]
        x = np.array(x)
        if x.ndim == 0:
            x = x[None]
        return [[aa_1_N.get(a, states - 1) for a in y] for y in x]

    def N_to_AA(x):
        # [[0,1,2,3]] -> ["ARND"]
        x = np.array(x)
        if x.ndim == 1:
            x = x[None]
        return ["".join([aa_N_1.get(a, "-") for a in y]) for y in x]

    def parse_PDB_biounits(x, atoms=["N", "CA", "C"], chain=None):
        """
        input:  x = PDB filename
                atoms = atoms to extract (optional)
        output: (length, atoms, coords=(x,y,z)), sequence
        """
        xyz, seq, min_resn, max_resn = {}, {}, 1e6, -1e6
        for line in open(x, "rb"):
            line = line.decode("utf-8", "ignore").rstrip()

            if line[:6] == "HETATM" and line[17 : 17 + 3] == "MSE":
                line = line.replace("HETATM", "ATOM  ")
                line = line.replace("MSE", "MET")

            if line[:4] == "ATOM":
                ch = line[21:22]
                if ch == chain or chain is None:
                    atom = line[12 : 12 + 4].strip()
                    resi = line[17 : 17 + 3]
                    resn = line[22 : 22 + 5].strip()
                    x, y, z = [float(line[i : (i + 8)]) for i in [30, 38, 46]]

                    if resn[-1].isalpha():
                        resa, resn = resn[-1], int(resn[:-1]) - 1
                    else:
                        resa, resn = "", int(resn) - 1
                    #         resn = int(resn)
                    if resn < min_resn:
                        min_resn = resn
                    if resn > max_resn:
                        max_resn = resn
                    if resn not in xyz:
                        xyz[resn] = {}
                    if resa not in xyz[resn]:
                        xyz[resn][resa] = {}
                    if resn not in seq:
                        seq[resn] = {}
                    if resa not in seq[resn]:
                        seq[resn][resa] = resi

                    if atom not in xyz[resn][resa]:
                        xyz[resn][resa][atom] = np.array([x, y, z])

        # convert to numpy arrays, fill in missing values
        seq_, xyz_ = [], []
        try:
            for resn in range(min_resn, max_resn + 1):
                if resn in seq:
                    for k in sorted(seq[resn]):
                        seq_.append(aa_3_N.get(seq[resn][k], 20))
                else:
                    seq_.append(20)
                if resn in xyz:
                    for k in sorted(xyz[resn]):
                        for atom in atoms:
                            if atom in xyz[resn][k]:
                                xyz_.append(xyz[resn][k][atom])
                            else:
                                xyz_.append(np.full(3, np.nan))
                else:
                    for atom in atoms:
                        xyz_.append(np.full(3, np.nan))
            return np.array(xyz_).reshape(-1, len(atoms), 3), N_to_AA(np.array(seq_))
        except TypeError:
            return "no_chain", "no_chain"

    pdb_dict_list = []
    c = 0

    init_alphabet = [
        "A",
        "B",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "J",
        "K",
        "L",
        "M",
        "N",
        "O",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "U",
        "V",
        "W",
        "X",
        "Y",
        "Z",
        "a",
        "b",
        "c",
        "d",
        "e",
        "f",
        "g",
        "h",
        "i",
        "j",
        "k",
        "l",
        "m",
        "n",
        "o",
        "p",
        "q",
        "r",
        "s",
        "t",
        "u",
        "v",
        "w",
        "x",
        "y",
        "z",
    ]
    extra_alphabet = [str(item) for item in list(np.arange(300))]
    chain_alphabet = init_alphabet + extra_alphabet

    if input_path is not None:
        biounit_names = glob.glob(input_path + "/*.pdb")
    else:
        biounit_names = input_list

    for biounit in biounit_names:
        my_dict = {}
        s = 0
        concat_seq = ""
        concat_N = []
        concat_CA = []
        concat_C = []
        concat_O = []
        concat_mask = []
        coords_dict = {}
        for letter in chain_alphabet:
            sidechain_atoms = ["CA"] if ca_only else ["N", "CA", "C", "O"]
            xyz, seq = parse_PDB_biounits(biounit, atoms=sidechain_atoms, chain=letter)
            if type(xyz) != str:
                concat_seq += seq[0]
                my_dict["seq_chain_" + letter] = seq[0]
                coords_dict_chain = {}
                if ca_only:
                    coords_dict_chain["CA_chain_" + letter] = xyz.tolist()
                else:
                    coords_dict_chain["N_chain_" + letter] = xyz[:, 0, :].tolist()
                    coords_dict_chain["CA_chain_" + letter] = xyz[:, 1, :].tolist()
                    coords_dict_chain["C_chain_" + letter] = xyz[:, 2, :].tolist()
                    coords_dict_chain["O_chain_" + letter] = xyz[:, 3, :].tolist()
                my_dict["coords_chain_" + letter] = coords_dict_chain
                s += 1
        fi = biounit.rfind("/")
        my_dict["name"] = biounit[(fi + 1) : -4]
        my_dict["num_of_chains"] = s
        my_dict["seq"] = concat_seq
        if s < len(chain_alphabet):
            pdb_dict_list.append(my_dict)
            c += 1

    if save_path:
        with open(save_path, "w") as f:
            for entry in pdb_dict_list:
                f.write(json.dumps(entry) + "\n")

    return pdb_dict_list
