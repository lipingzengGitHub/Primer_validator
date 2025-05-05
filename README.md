## WGS Primer Validation Tool

A flexible and advanced primer validation pipeline designed for validating primer pairs for RSV Whole Genome Sequencing (WGS). This tool integrates primer property calculations, in silico PCR with adjustable mismatch tolerance, off-target checking using BLAST, and automatic scoring of primer pairs to assist in selecting the best candidates.

###  Features
- Batch validation of primer pairs (F and R)
- In silico PCR simulation with user-defined mismatch tolerance
- Automatic BLAST database creation and off-target checking
- Comprehensive scoring system based on Tm difference, off-target results, and amplifiability
- Easy-to-read Excel report generation with scores and rankings

###  Requirements
- Python 3.8+
- Biopython
- Pandas
- BLAST+ (makeblastdb, blastn)

Install dependencies:

```bash
pip install biopython pandas
sudo apt-get install ncbi-blast+  # for Ubuntu/Debian
```

###  Input Files
- **Primer CSV file** (comma separated F,R pairs):

```
ATGCAGCTGATCGATCGTAG,CGATGCTAGCTACGATCGA
GCTAGCTAGCTGACTGCTA,TACGATCGATCGGCTAGCA
```

- **Reference genome FASTA file**

###  Usage

```bash
python primer_validator.py --primer_file primers.csv --reference RSV.fasta --max_mismatch 2
```

**Parameters:**
- `--primer_file`: Path to CSV file with primer pairs
- `--reference`: Path to reference genome FASTA file
- `--max_mismatch`: Maximum mismatches allowed in in silico PCR (default: 2)

If you run the program without parameters:

```bash
python primer_validator.py
```

The terminal will automatically print the usage help message with available arguments and instructions.

###  Output

`RSV_WGS_primer_validation_full.xlsx` (Excel file)

| Forward Primer | Reverse Primer | Fwd_Tm | Rev_Tm | Amplifiable | Off Target? | Score |
|----------------|----------------|--------|--------|-------------|-------------|-------|
| ATGCAGCT...    | CGATGCTAG...   | 60.2   | 59.8   | True        | No          | 98    |

- Ranked by Score, best primer pairs first

###  Notes
- The tool will automatically create BLAST index files if not present
- Fuzzy matching supports degenerate bases (N -> ATCG)
- Best suited for viral and pathogen WGS primer screening workflows

###  Application Scenario: Amplicon-based Whole Genome Sequencing (WGS)
This tool is particularly designed for **amplicon-based WGS**, where primers are used to amplify tiled regions across the viral genome (such as RSV) before Illumina sequencing. Compared to shotgun WGS, amplicon-based WGS is more efficient for smaller genomes, low input RNA/DNA samples, and provides uniform coverage essential for pathogen surveillance and variant detection.

Typical Applications:
- RSV, Influenza and SARS-CoV-2 genome sequencing
- Viral genomic surveillance
- Targeted sequencing workflows
