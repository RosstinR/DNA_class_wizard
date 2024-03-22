#!/usr/bin/env python
# coding: utf-8

## Helpful Variables

standard_code = {
    "UUU": "F",
    "UUC": "F",
    "UUA": "L",
    "UUG": "L",
    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",
    "UAU": "Y",
    "UAC": "Y",
    "UAA": "*",
    "UAG": "*",
    "UGA": "*",
    "UGU": "C",
    "UGC": "C",
    "UGG": "W",
    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",
    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAU": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AUU": "I",
    "AUC": "I",
    "AUA": "I",
    "AUG": "M",
    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAU": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGU": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",
    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAU": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

aa_mol_weights = {
    "A": 89.09,
    "C": 121.15,
    "D": 133.1,
    "E": 147.13,
    "F": 165.19,
    "G": 75.07,
    "H": 155.16,
    "I": 131.17,
    "K": 146.19,
    "L": 131.17,
    "M": 149.21,
    "N": 132.12,
    "P": 115.13,
    "Q": 146.15,
    "R": 174.2,
    "S": 105.09,
    "T": 119.12,
    "V": 117.15,
    "W": 204.23,
    "X": 0,
    "Y": 181.19,
}


## class seq


class seq:

    def __init__(self, name, organism, sequence, type):
        self.name = name
        self.organism = organism
        self.sequence = sequence
        self.type = type

    def info(self):
        print(self.name)
        print(self.organism)
        print(self.sequence)
        print(self.type)

    def length(self):
        print(len(self.sequence))

    def fasta_out(self):
        f = open("{}.fa".format(self.name), "w")
        f.write(
            ">"
            + self.name
            + "_"
            + self.organism
            + "_"
            + self.type
            + "\n"
            + self.sequence
        )
        f.close()


## class protein


class protein(seq):
    def __init__(self, name, organism, sequence, type, size):
        self.size = size
        super().__init__(name, organism, sequence, type)

    def info(self):
        print(self.name)
        print(self.organism)
        print(self.sequence)
        print(self.type)
        print(self.size)

    def mol_weight(self):
        mol_weight = sum(aa_mol_weights[i] for i in self.sequence)
        print(mol_weight)

    def fasta_out(self):
        f = open("{}.fa".format(self.name), "w")
        f.write(
            ">"
            + self.name
            + "_"
            + self.organism
            + "_"
            + self.type
            + "\n"
            + self.sequence
            + "\n"
            + self.size
        )
        f.close()


## class nucleotide


class nucleotide(seq):
    def __init__(self, name, organism, sequence, type):
        super().__init__(name, organism, sequence, type)

    def gc_content(self):
        gc_content = (
            (self.sequence.count("C") + self.sequence.count("G"))
            / len(self.sequence)
            * 100
        )
        print(gc_content)


## class DNA


class DNA(nucleotide):
    def __init__(self, name, organism, sequence, type):
        super().__init__(name, organism, sequence, type)

    def transcribe(self):
        transcribe = self.sequence.replace("T", "U")
        print(transcribe)

    def reverse_complement(self):
        new_letters = {"A": "T", "C": "G", "G": "C", "T": "A"}
        reverse_complement = "".join(
            new_letters.get(i, i) for i in (self.sequence[::-1])
        )
        print(reverse_complement)

    def six_frames(self):
        frames = []
        new_letters = {"A": "T", "C": "G", "G": "C", "T": "A"}
        reverse_complement_frames = "".join(
            new_letters.get(i, i) for i in (self.sequence[::-1])
        )
        for n in range(3):
            frames.append(self.sequence[n:])
        for n in range(3):
            frames.append(reverse_complement_frames[n:])
        print(frames)


## class RNA


class RNA(nucleotide):
    def __init__(self, name, organism, sequence, type):
        super().__init__(name, organism, sequence, type)

    def start(self):
        start = self.sequence.find("AUG")
        print(start)

    def translate(self):
        start = self.sequence.find("AUG")
        translate = "".join(
            [
                standard_code.get(self.sequence[start : start + 3])
                for start in range(start, len(self.sequence) - 2, 3)
            ]
        )
        print(translate)


## test

uidA_DNA = DNA(
    name="uidA",
    organism="bacteria",
    type="DNA",
    sequence="CGCATGTTACGTCCTGTAGAAACCCCAACCCGTGAAATCAAAAAA",
)

uidA_DNA.fasta_out()
uidA_DNA.six_frames()
uidA_DNA.reverse_complement()
uidA_DNA.transcribe()


uidA_RNA = RNA(
    name="uidA_RNA",
    organism="bacteria",
    type="RNA",
    sequence="CGCAUGUUACGUCCUGUAGAAACCCCAACCCGUGAAAUCAAAAAA",
)

uidA_RNA.fasta_out()
uidA_RNA.translate()


uidA_protein = protein(
    name="uidA_protein",
    organism="bacteria",
    type="protein",
    size="27",
    sequence="MLRPVETPTREIKK",
)

uidA_protein.fasta_out()
uidA_protein.mol_weight()
