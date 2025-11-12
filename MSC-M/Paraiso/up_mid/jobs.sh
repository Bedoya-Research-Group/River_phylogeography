#!/bin/bash
#SBATCH --job-name=paraiso
#SBATCH --account=leache
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4           # four BPP runs
#SBATCH --cpus-per-task=3             # each uses 6 threads internally
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=48:00:00
#SBATCH --chdir=/gscratch/leache/adam/ana/paraiso/up_mid
#SBATCH --output=multi-serial.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=leache@uw.edu
#SBATCH --export=ALL

/gscratch/leache/adam/ana/bpp/bpp --cfile=/gscratch/leache/adam/ana/paraiso/up_mid/bpp1.ctl --no-pin &
/gscratch/leache/adam/ana/bpp/bpp --cfile=/gscratch/leache/adam/ana/paraiso/up_mid/bpp2.ctl --no-pin &
/gscratch/leache/adam/ana/bpp/bpp --cfile=/gscratch/leache/adam/ana/paraiso/up_mid/bpp3.ctl --no-pin &
/gscratch/leache/adam/ana/bpp/bpp --cfile=/gscratch/leache/adam/ana/paraiso/up_mid/bpp4.ctl --no-pin &

wait
