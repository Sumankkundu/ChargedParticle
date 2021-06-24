#!/bin/bash
# Random Seed
#define the directory
rm -rf incard_CP5
mkdir incard_CP5

present_dir=$PWD
input_dir="$PWD/incard_CP5"

#Name of the parameter and define directory
ParaName="CP5_Default"
rm -rf $present_dir/output_$ParaName
mkdir $present_dir/output_$ParaName
output_dir=$present_dir/output_$ParaName
mkdir $output_dir/logs
log_dir=$output_dir/logs

#declear -a para_vals
Ranseed=(10 1000000 2000000 3000000 40000000 50000000 60000000 7000000) #8000000 90000000 100000000 110000000)  #

npar=${#Ranseed[@]}  #number of value to check
#para_var=(1)   # Name the run
para_var=(1 2 3 4 5 6 7 8) #9 10 11 12)   # Name the run

#create input card for pythia.
cd $input_dir
rm *cmnd
cd $log_dir
rm *log

cd $present_dir
rm HardQCD
make HardQCD

rm -rf Run_$ParaName
mkdir Run_$ParaName
cd Run_$ParaName
#loop over the variable
for((i =1 ;i <${npar}+1 ; i++))
do
var=${para_var[$i-1]}
val=${Ranseed[$i-1]}
cat>$input_dir/incard_$var.cmnd<<%
!inputcard_$var
Main:numberOfEvents= 1000000
#-----------------Random Seed
Random:setSeed = on
Random:seed = $val
#-----------------CP5 Tune
Tune:pp =14
Tune:ee = 7

#PartonLevel:MPI = on
#PartonLevel:ISR = on
#PartonLevel:FSR = off
%

cp $present_dir/HardQCD $present_dir/Run_$ParaName/pyrun_${ParaName}_$var   #copy the exqutive
#create the run script
cat>run_${ParaName}_$var.sh<<%
#!/bin/bash
$present_dir/Run_$ParaName/./pyrun_${ParaName}_$var $input_dir/incard_$var.cmnd $output_dir/pp_13_${ParaName}_$var.root>>$log_dir/run_$ParaName$var.log
%

pwd
chmod u+x run_${ParaName}_$var.sh

./run_${ParaName}_$var.sh &

done   #End of loop on the parameter value

cp -r $input_dir $output_dir
cp $present_dir/complie.sh $output_dir
echo complete all
