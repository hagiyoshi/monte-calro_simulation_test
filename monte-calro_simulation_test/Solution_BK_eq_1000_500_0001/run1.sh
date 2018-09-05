#!/bin/bash
#--- パラメータ設定部
#SBATCH -J BK
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH -p M

#--- ジョブスクリプト本体
#export OMP_NUM_THREADS=4
pwd
#cd /xc/home/yoshikazu.hagiwara/OLD/non_local/bms/bms2/gg-H
#touch numb1.txt
srun ./RK_BK_eq_rel_simp
