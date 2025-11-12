#!/bin/bash
#SBATCH --job-name=diego
#SBATCH --account=biology
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4           # four BPP runs
#SBATCH --cpus-per-task=4             # each uses 6 threads internally
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=24:00:00
#SBATCH --chdir=/gscratch/leache/adam/ana/diego/ancestor_down
#SBATCH --output=multi-serial.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=leache@uw.edu
#SBATCH --export=ALL

/gscratch/leache/adam/ana/bpp/bpp --cfile=/gscratch/leache/adam/ana/diego/ancestor_down/bpp1.ctl --no-pin &
/gscratch/leache/adam/ana/bpp/bpp --cfile=/gscratch/leache/adam/ana/diego/ancestor_down/bpp2.ctl --no-pin &
/gscratch/leache/adam/ana/bpp/bpp --cfile=/gscratch/leache/adam/ana/diego/ancestor_down/bpp3.ctl --no-pin &
/gscratch/leache/adam/ana/bpp/bpp --cfile=/gscratch/leache/adam/ana/diego/ancestor_down/bpp4.ctl --no-pin &
wait
