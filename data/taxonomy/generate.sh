#!/bin/bash

wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xvf  taxdump.tar.gz names.dmp nodes.dmp
rm taxdump.tar.gz

