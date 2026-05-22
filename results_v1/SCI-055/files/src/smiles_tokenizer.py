"""SMILES Tokenizer for seq2seq retrosynthesis model."""

import re
from typing import List, Dict, Optional


SMILES_PATTERN = re.compile(
    r"(\[[^\]]+\]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\|\/|:|~|@|\?|>>?|\*|\$|%[0-9]{2}|[0-9])"
)


class SMILESTokenizer:
    """Tokenizer that splits SMILES strings into chemically meaningful tokens."""

    PAD = "<PAD>"
    SOS = "<SOS>"
    EOS = "<EOS>"
    UNK = "<UNK>"

    def __init__(self):
        self.token2idx: Dict[str, int] = {}
        self.idx2token: Dict[int, str] = {}
        self._special_tokens = [self.PAD, self.SOS, self.EOS, self.UNK]
        for i, tok in enumerate(self._special_tokens):
            self.token2idx[tok] = i
            self.idx2token[i] = tok

    def tokenize(self, smiles: str) -> List[str]:
        tokens = SMILES_PATTERN.findall(smiles)
        return tokens

    def build_vocab(self, smiles_list: List[str]):
        all_tokens = set()
        for smi in smiles_list:
            all_tokens.update(self.tokenize(smi))
        for tok in sorted(all_tokens):
            if tok not in self.token2idx:
                idx = len(self.token2idx)
                self.token2idx[tok] = idx
                self.idx2token[idx] = tok

    def encode(self, smiles: str, max_len: Optional[int] = None) -> List[int]:
        tokens = self.tokenize(smiles)
        ids = [self.token2idx.get(self.SOS)]
        for t in tokens:
            ids.append(self.token2idx.get(t, self.token2idx[self.UNK]))
        ids.append(self.token2idx.get(self.EOS))
        if max_len:
            ids = ids[:max_len]
            ids += [self.token2idx[self.PAD]] * (max_len - len(ids))
        return ids

    def decode(self, ids: List[int]) -> str:
        tokens = []
        for i in ids:
            tok = self.idx2token.get(i, self.UNK)
            if tok == self.EOS:
                break
            if tok not in (self.PAD, self.SOS):
                tokens.append(tok)
        return "".join(tokens)

    @property
    def vocab_size(self) -> int:
        return len(self.token2idx)
