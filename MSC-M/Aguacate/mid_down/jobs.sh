#!/bin/bash
#SBATCH --job-name=aguacate
#SBATCH --account=biology
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4           # four BPP runs
#SBATCH --cpus-per-task=3             # each uses 6 threads internally
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=48:00:00
#SBATCH --chdir=/gscratch/leache/adam/ana/aguacate/mid_down/
#SBATCH --output=multi-serial.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=leache@uw.edu
#SBATCH --export=ALL

/gscratch/leache/adam/ana/bpp/bpp --cfile=/gscratch/leache/adam/ana/aguacate/mid_down/bpp1.ctl --no-pin &
/gscratch/leache/adam/ana/bpp/bpp --cfile=/gscratch/leache/adam/ana/aguacate/mid_down/bpp2.ctl --no-pin &
/gscratch/leache/adam/ana/bpp/bpp --cfile=/gscratch/leache/adam/ana/aguacate/mid_down/bpp3.ctl --no-pin &
/gscratch/leache/adam/ana/bpp/bpp --cfile=/gscratch/leache/adam/ana/aguacate/mid_down/bpp4.ctl --no-pin &

wait
