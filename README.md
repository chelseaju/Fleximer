# Fleximer
Fleximer: Accurate Quantification of RNA-Seq via Variable Length k-mers

## Download
```
git clone https://github.com/chelseaju/Fleximer.git
```

## Installation
Fleximer is implemented in C++11, and requires compiler g++ 4.7 and above. Please run the following commands to compile the executables:
```
cd Fleximer
make
```
Binary files are also provided:
```
sigmer_generation, sigmer_selection, sigmer_count, EM
```

## Data
Fleximer contains two stages: *sig-mer* preparation and RNA-Seq quantification. In *sig-mer* preparation, Fleximer takes a predefined partition of transcriptome. The partition divide transcripts into a set of non-overlapping clusters in a specialized fasta format. A regular fasta format contains a header line of transcript ID follow by one sequence:
```
>seq1
ACGTACGTA
```
A specialized fasta format contains a header line of a cluster ID and a set of member transcript IDs, follow by member transcript sequence. Each member is separated by '|':
```
>clusterX|seq1|seq2
ACGTACGTA|ACGAACGTA
```
The partition can be based on sequence similarity, or functional annotation. By default, we partition the transcriptome based on sequence similarity as described in [RNA-Skim](https://github.com/chelseaju/RNASkim). A human transcriptome partition is provided in /Data/reference/human_60mers.fa.

In RNA-Seq quantification stage, Fleximer takes the RNA-Seq data as input. Several examples simulated from [polyester](https://github.com/alyssafrazee/polyester) are provided in Data/test.

## *Sig-mer* Preparation
Fleximer relies on the structure and properties of suffix tree to discover all the possible *sig-mers*. It then constructs *sig-mers* overlap graphs to guide the selection. sigmer_generation takes the specialized fasta file of transcriptome partition as input, and writes all *sig-mers* to an intermediate file. sigmer_selection takes the output from sigmer_generation and a desired read length to produces a list of representative *sig-mers* for the quantification stage.

```
./sigmer_generation Data/reference/human_60mers.fa Data/sigmers.txt
./sigmer_selection Data/reference/human_60mers.fa Data/sigmers.txt 75 Data/selected_sigmer 20
```

## RNA-Seq Quantification
Fleximer determines the reads origin through the representative *sig-mers*. It uses the Aho-Corasick algorithm to determine the *sig-mers* occurrences for each read. For reads that may come from more than one origin, Fleximer uses the EM algorithm to resolve the ambiguity.
```
./sigmer_count Data/reference/human_60mers.fa Data/selected_sigmer.txt Data/test/len75_sample_46/eClass_count.txt Data/test_len75_sample_46/sample_01.fq
./EM Data/reference/human_60mers.fa Data/test/len75_sample_46/eClass_count.txt 75 Data/test/len75_sample_46/abundance.txt 
```

The final transcript abundance for each transcript is reported in read count, RPKM, and TPM.
