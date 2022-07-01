# sRNadatabace Data Collection Pipeline

Project dedicated to running the data collection for the sRNA database.

Make sure your .config/snakemake/slurm_all/ folder is setup correctly

Recommended execution:
```console
user@main:~$ snakemake -p -j 999 --profile slurm_all  --latency-wait 60  --use-conda   --group-components high_IO_mini_job=30
```
