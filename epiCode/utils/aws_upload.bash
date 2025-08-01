# AWS

# backup Skynet

aws s3 sync /home/michalula/code  s3://streets-lab-data/michalula/Skynet/code

aws s3 sync /home/michalula/data  s3://streets-lab-data/michalula/Skynet/data

aws s3 sync /home/michalula/miniconda3  s3://streets-lab-data/michalula/Skynet/miniconda3



aws s3 sync /var/lib/minknow/data/MR_NT_Day28_07212025/no_sample_id/20250721_2353_MN36507_FBD30760_a93f7f57/ s3://streets-lab-data/michalula/nanoporedata/0721202_Day28_Tcells_NonTreated/run1_same_flowcell/

aws s3 sync /var/lib/minknow/data/MR_rerun_same_flowcell_NT_Day28_07222025/no_sample_id/20250722_1343_MN36507_FBD30760_cf22e917  s3://streets-lab-data/michalula/nanoporedata/0721202_Day28_Tcells_NonTreated/run2_same_flowcell/

aws s3 sync /var/lib/minknow/data/queued_reads/MR_T_cells_CRoff_Day28_07212025/  s3://streets-lab-data/michalula/nanoporedata/0721202_Day28_Tcells_CRISPRoff/run1_same_flowcell/

aws s3 sync /var/lib/minknow/data/MR_runLonger_Tcells_CRoff_Day28_07222025/no_sample_id/20250722_1518_MN31715_FBD30913_b890222c/  s3://streets-lab-data/michalula/nanoporedata/0721202_Day28_Tcells_CRISPRoff/run2_same_flowcell/


# aws s3 cp s3://streets-lab-data/michalula/nanoporedata/20241226_nCATS_K562_ZFPOFFpostSort_HIGH/no_sample_id/20241227_0140_MN36507_FAW26101_7cfd685d/ ./20241226_nCATS_K562_ZFPOFFpostSort_HIGH
