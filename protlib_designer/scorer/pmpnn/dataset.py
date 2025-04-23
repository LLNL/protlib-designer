import json
import time


class StructureDataset:
    def __init__(
        self,
        jsonl_file,
        verbose=True,
        truncate=None,
        max_length=100,
        alphabet="ACDEFGHIKLMNPQRSTVWYX-",
    ):
        alphabet_set = set(list(alphabet))
        discard_count = {"bad_chars": 0, "too_long": 0, "bad_seq_length": 0}

        with open(jsonl_file) as f:
            self.data = []

            lines = f.readlines()
            start = time.time()
            for i, line in enumerate(lines):
                entry = json.loads(line)
                seq = entry["seq"]
                name = entry["name"]

                if bad_chars := set(list(seq)).difference(alphabet_set):
                    if verbose:
                        print(name, bad_chars, entry["seq"])
                    discard_count["bad_chars"] += 1

                elif len(entry["seq"]) <= max_length:
                    self.data.append(entry)
                else:
                    discard_count["too_long"] += 1
                # Truncate early
                if truncate is not None and len(self.data) == truncate:
                    return

                if verbose and (i + 1) % 1000 == 0:
                    elapsed = time.time() - start
                    print(
                        "{} entries ({} loaded) in {:.1f} s".format(
                            len(self.data), i + 1, elapsed
                        )
                    )
            if verbose:
                print("discarded", discard_count)

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        return self.data[idx]


class StructureDatasetPDB:
    def __init__(
        self,
        pdb_dict_list,
        verbose=True,
        truncate=None,
        max_length=100,
        alphabet="ACDEFGHIKLMNPQRSTVWYX-",
    ):
        alphabet_set = set(list(alphabet))
        discard_count = {"bad_chars": 0, "too_long": 0, "bad_seq_length": 0}

        self.data = []

        start = time.time()
        for i, entry in enumerate(pdb_dict_list):
            seq = entry["seq"]
            # name = entry["name"]

            bad_chars = set(list(seq)).difference(alphabet_set)
            if bad_chars:
                discard_count["bad_chars"] += 1

            elif len(entry["seq"]) <= max_length:
                self.data.append(entry)
            else:
                discard_count["too_long"] += 1
            # Truncate early
            if truncate is not None and len(self.data) == truncate:
                return

            if verbose and (i + 1) % 1000 == 0:
                elapsed = time.time() - start
                print(
                    "{} entries ({} loaded) in {:.1f} s".format(
                        len(self.data), i + 1, elapsed
                    )
                )

            # print('Discarded', discard_count)

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        return self.data[idx]
