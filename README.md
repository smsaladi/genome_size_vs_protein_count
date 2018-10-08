# Genome size vs. gene count

A Shiny app that plots an up-to-date version of the total number of annotated genes
in genomes submitted to GenBank as a function of genome size (based on data provided by
[NCBI genome](https://www.ncbi.nlm.nih.gov/genome) reports). The app stays up-to-date by
retrieving data from the [NCBI FTP server](https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS)
 upon loading.

Image on [Wikipedia](https://en.wikipedia.org/wiki/Genome_size#/media/File:Genome_size_vs_number_of_genes.svg)
is a frozen version of one created by this app (using the Download SVG button).

URL: [https://saladi.shinyapps.io/genome_size_vs_gene_count/](https://saladi.shinyapps.io/genome_size_vs_gene_count/)

## Notes

Commits to branches are automatically deployed to shiny via a
[Travis CI](https://travis-ci.org/smsaladi/genome_size_vs_gene_count) build.

