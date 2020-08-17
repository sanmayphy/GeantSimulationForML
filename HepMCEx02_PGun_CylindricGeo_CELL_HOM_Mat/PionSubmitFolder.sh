#Submit N jobs for a given foder
# for S queue use 3000 jobs
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1


Njobs=300
TargetFolder=$1

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 DIRECTORY" >&2
  exit 1
fi

#EOfolder=/storage/agrp/antonc/LOGFILES
for i in `seq $Njobs`; do
  run_number=$i
  echo qsub -e $TargetFolder/$run_number -o $TargetFolder/$run_number -N N${run_number} -q N -l nodes=1:ppn=2 -l mem=2gb  -v run_number=$run_number,DIR=$TargetFolder Pionrun.sh
  mkdir $TargetFolder/$run_number
  qsub -e $TargetFolder/$run_number -o $TargetFolder/$run_number -N N${run_number} -q N -l nodes=1:ppn=2 -l mem=2gb  -v run_number=$run_number,DIR=$TargetFolder Pionrun.sh
done


