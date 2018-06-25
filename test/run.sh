# example for paired tumor and normal bams 
../msisensor msi -d example.microsate.sites -n example.normal.bam -t example.tumor.bam -e example.bed -o example.paired.output -l 1 -q 1 -b 2

# example for tumor only bams
../msisensor msi -d example.microsate.sites -t example.tumor.bam -e example.bed -o example.tumor.output -l 1 -q 1 -b 2

