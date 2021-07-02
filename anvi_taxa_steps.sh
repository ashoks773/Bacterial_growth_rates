#!/bin/bash

anvi-run-hmms -c contigs.db -I Bacteria_71
anvi-run-scg-taxonomy -c contigs.db

anvi-estimate-scg-taxonomy -c contigs.db --metagenome-mode
anvi-estimate-scg-taxonomy -c contigs.db -p SAMPLES-MERGED/PROFILE.db --metagenome-mode --compute-scg-coverages 

anvi-estimate-scg-taxonomy -c contigs.db -p SAMPLES-MERGED/PROFILE.db -C metabat
anvi-estimate-genome-completeness -c contigs.db -p SAMPLES-MERGED/PROFILE.db -C metabat
