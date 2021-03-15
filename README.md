# Vibrio-campbellii-DS40M4-bulk-gene-deletion


## Usage

The script is run on combined fasta and GFF files and as output will create an excel file.

```bash
cat CP030788.1.fa CP030789.1.fa CP030790.1.fa > combined.fa
cat CP030788.1.gff3 CP030789.1.gff3 CP030790.1.gff3 > combined.gff3
./make-candidates.pl combined.fa combined.gff3
```

The script depends on the following perl modules:

* BioPerl
* Array::Circular
* Excel::Writer::XLSX
