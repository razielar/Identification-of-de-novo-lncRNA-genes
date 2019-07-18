#!/bin/bash

mkdir qsub_logs

echo "host: $HOSTNAME" > snakemake.log; #rsync -aP --exclude '.git' /users/rg/jlagarde/projects/software_utilities/LR-Seq/ ./; 
snakemake -p --reason --latency-wait 100 -j 4500 --jobname {rulename}_{jobid}.sh --cluster-config cluster_config.json --max-jobs-per-second 10 --drmaa " -V -q {cluster.queue} -l virtual_free={cluster.virtual_free} -l h={cluster.node} -l h_rt={cluster.h_rt}  -o {cluster.out} -e {cluster.err} {cluster.threads}" --rerun-incomplete  2>> snakemake.log;
tail -n2 snakemake.log | mail -s "snakemake run test finished" raziel.amador@crg.es

mv snakemake.log logs/; mv qsub_logs/ logs/

#Generate DAG: 

snakemake -p --forceall --rulegraph > test_dag.dot
cat test_dag.dot | dot -Tsvg > test_dag.svg
rm test_dag.dot








