# Coding Challenge - Finding two instances of cas9 from gbff files

### Intro:
This was a coding challenge for a SWE position for a Series A funded biotech company (as of December 2020), founded by a Professor and graduates from MIT. Given two GenBank Flat Files ([`gbff`](https://www.ncbi.nlm.nih.gov/datasets/docs/about-ncbi-gbff/)) that contain nucleotide sequences and its metadata, can we predict two instances of cas9? 


### Descriptions:
The developed code is designed to 
1. compute set of features useful for cas9 prediction for each coding sequence (CDS) 
     * Number of nucleotides to the nearest repeated sequence. If not observed, displays `infinity`
     * Protein size. Simplified to the length of the amino acids
     * Useful metadata. Can be obtained by using `SeqIO.parse(gbff_directory, format='genbank')`
2. Querry (1) into a DataFrame and output a csv (row=CDS, column=features)
3. For each CDS near the repeated region, find whether the candidate amino acid sequence is similar to cas9 statistically using NCBI's [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome). 


### Requirements
1. macOS Mojave 10.14.6
2. Python       3.7.6 
3. BioPython    1.78 (`!pip3 install biopython`)
4. Pandas       1.1.5


### How to run
1. Download the source code, `solution.py`
2. Download the data, `GCF_001239625.1_7068_7_24_genomic.gbff` and `GCF_002014815.1_ASM201481v1_genomic.gbff`
3. Open terminal, navigate to the folder that has 1. and 2. above
4. Use appropriate arguments for the following results:


### Input arguments
```
-h --help:                      optional, help command
-b --blast:                     optional, run ncbi's BLAST database. Used to check the filter for discovering new filter (Cas9) sequences
                                          takes foerever to run
-g --gbff <directory>:          required, enter gbff file directory for analysis
-c --check_common_features:     optional, used to check the most common features of the CDS.
-a --add_features               optional, used to add features to add to the output csv. If entering multiple features, please use the format "-a feature1 -a feature2 ...".
-k --keyword:                   optional, used for determining a set of criteria
```

### Usage
1. Given a gbff file, compute the most common coding sequence (CDS) features<br /> 
    `!python3 solution.py -g <directory_to_gbff_file> -c`
2. Given a gbff file and feature(s), output a csv with rows=CDS and column=features. Default features will be read for each CDS with the argument `-a default`. Additional features are welcome with the arguments `-a <feature1> -a <feature2>`. Duplicate features are omitted.<br /> 
    `!python3 solution.py -g <directory_to_gbff_file> -a default`<br /> 
    or
      `!python3 solution.py -g <directory_to_gbff_file> -a <feature1> -a <feature2> ... -a <featureN>`<br /> 
      or
    `!python3 solution.py -g <directory_to_gbff_file> -a default -a <feature1> -a <feature2> ... -a <featureN>`<br /> 
3. Given a gbff file and keyword, output csv that has a non-case sensitive match with the keyword (in the given assignment, it was cas9) in each feature. keyword must be cas\d, ex. "cas9". Output's a filtered csv that only contains rows with a match <br /> 
    `!python3 solution.py -g <directory_to_gbff_file> -k <keyword>`
4. Given a gbff file and keyword, output a textfile sending each CDS's translation (Amino acid sequence) to NCBI's BLAST database. Only CDS near the repeated regions are sent to the database (since CRIPSR sequences repeats). Output's a textfile with 
    a. `Sequences producing significant alignments`
    b. `Score (Bits)`
    c. `E-value`
    d. `Max (percent) identical` <br /> 

    `!python3 solution.py -g <directory_to_gbff_file> -k <keyword> -b`


### Output Example
```
(base) Gs-MacBook-Pro: ~george$ python3 solution.py -g /Users/george/Desktop/GCF_002014815.1_ASM201481v1_genomic.gbff -k cas9 -b -a default
----------------------------------------------------------------------------------------------------
Input Arguments: 

GBFF:                  /Users/george/Desktop/GCF_002014815.1_ASM201481v1_genomic.gbff
CHECK_COMMON_FEATURES: False
ADDITIONAL FEATURE(S): ['default']
KEYWORD:               cas9
RUN BLAST              True
----------------------------------------------------------------------------------------------------
Obtaining CSV for the features  1.default_values  

default_values are: 1.organism  2.record_id  3.feature_location  4.locus_tag  5.protein_id  6.gene  7.protein_size  8.distance_to_repeat_region  9.translation  

Organizing to DataFrame and searching for repeat region (for CRISPR identifier..)

Found repeat region! Used for identifying Cas candidate sequences


Found repeat region! Used for identifying Cas candidate sequences

Appeneded a total of 29 records, 3624 features
Saving the DataFrame to "output/CDS_features/CDS_features_default_GCF_002014815.1_ASM201481v1_genomic_gbff_22122020_003537.csv"
----------------------------------------------------------------------------------------------------
Searching for keyword cas9 in DataFrame..

cas9 is present at row index 714. Protein ID is WP_010922251.1 

Saving DataFrame with the filter/keyword: cas9
Filter "cas9" found from feature "gene"
Saving filtered DataFrame with keyword cas9 to "output/filter_cas9/cas9_22122020_003537.csv"
----------------------------------------------------------------------------------------------------
Searching for potential cas9 candidate sequences via BLAST

Only consider CDS near repeated region if any..
Sending the Amino acid sequence to BLAST and retrieving the results...
Near repeated region sequence found!
Blasting protein_id: "WP_002989532.1", Amino acid: "MAFGENGPRKKTTFEKVTMFVVILMVLVTVGGLIASALSVLM"
This may take a while
```
    


```python

```
