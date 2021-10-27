# Scoring of multiple sequence alignment

The alignment scoring is calculated with the help of [**pyMSA**](https://pypi.org/project/pyMSA/) package, which contains a few internal distance matrices (*PAM250*, *Blosum62*). Also **pyMSA** allows to create custom distance matrix with the help of **_FileMatrix_** class, by simply specifying path to the substitution matrix (*ROCS* in our case).

However, such matrix will not be useful for scoring of the MSA, which include non-natural amino acids (AAs). Algorithm of non-natural AAs encoding is a little bit tricky. It does not take non-natural monomers directly from *ROCS* or *monomers_map.txt* files. Current algorithm randomly assigns non-natural AAs to the unicode characters in range of the first 256 symbols. This results in some inconsistency and the need to create distance matrix from custom substitution matrix.

## Explanation of inconsistency

*ROCS* file contains square matrix of single elements with scores for their alignment. Each symbol, encodes a certain non-natural monomer. These can be decoded with the help of *monomers_map.txt* file.
Below is example of contents of *monomers_map.txt* file:

	symbol	Unicode	polymertype	letter3	name							monomertype
	...
	Trp1Me	0100	PEPTIDE		WA7	(2~{S})-2-amino-3-(1-methylindol-3-yl)propanoic acid	Backbone
	ClAc	0178	PEPTIDE		GX7	2-chloroacetic acid					Backbone
	...

However, this API dynamically allocates non-natural AAs to different unicode symbols.
Below is example of contents of custom substitution matrix for one type of polymer (e.g. *PEPTIDE1*):

	...
	0x1b 0x0f -1   # 4Pal 	x ClAc
	0x41 0x0f  5   # A 	x ClAc
	0x01 0x0f  9   # Ac 	x ClAc
	...

And here is example of contents of custom substitution matrix for another type of polymer (e.g. *PEPTIDE2*):

	...
	0x41 0x01  5   # A 	x ClAc
	0x03 0x01  0   # Ahp 	x ClAc
	0x05 0x01 -3   # Bip 	x ClAc
	...

Look at *"ClAc"* monomer. It has **_"0178"_** code in *monomers_map.txt*, but **_"0x0f"_** and **_"0x01"_** codes in custom substitution matrix. So, non-natural AAs are encoded with unicode symbols from custom substitution matrix, and trying to score the alignment with the help of *ROCS* substitution matrix will not succeed.

## Custom class for the generation of distance matrix

ScoringUtils.py contains new class **_CustomMatrix_**, which allows creation of distance matrix from the custom substitution matrix. This is also much faster, as *ROCS* matrix can be as big, as *1500x1500* (dimensions increase with the addition of new monomers), while custom substitution matrix has maximum dimensions of *256x256*.