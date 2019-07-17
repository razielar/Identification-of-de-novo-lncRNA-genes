#!/bin/bash

mv qsub_logs qsub_logs_bkp; mkdir qsub_logs; rm -rf qsub_logs_bkp/ &
#--restart-times 2
# do not use --keep-going, it makes snakemake not delete output files from failed jobs!!
echo "host: $HOSTNAME" >  snakemake.log; rsync -aP --exclude '.git' /users/rg/jlagarde/projects/software_utilities/LR-Seq/ ./; 
snakemake -p --reason --latency-wait 100 --configfile config_cls3.json -s master.smk 
	-j 4500 --jobname {rulename}.{jobid}._{config["PROJECT_NAME"]}_.sh --cluster-config cluster_config.json 
	--max-jobs-per-second 10 --drmaa 
	" -V -q {cluster.queue} -l virtual_free={cluster.virtual_free} -l h={cluster.node} -l h_rt={cluster.h_rt}  -o {cluster.out} 
	-e {cluster.err} {cluster.threads} -P {cluster.project}" --rerun-incomplete  2>> snakemake.log; 
	tail -n2 snakemake.log |mail -s "snakemake run cls3 finished" julien.lagarde@crg.es

#generate DAG
# snakemake -p --configfile config_cls3.json -s master.smk  --forceall --rulegraph > cls3_dag.dot

# cat cls3_dag.dot | dot -Tsvg > cls3_dag.svg
