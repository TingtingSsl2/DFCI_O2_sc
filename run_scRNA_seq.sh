#!/bin/sh

sbatch -p short -t 6:00:00 --mem=64G -e %j.err -o %j.out /home/tiz228/scripts/scRNA_seq.sh
