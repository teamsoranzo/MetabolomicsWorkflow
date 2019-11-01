path=/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetalworker_inverse_norm

for i in {1..22} "X"
do
echo "Met #lines" > ${path}/checkfilelinesChr${i}.txt
wc -l   ${path}/Chr${i}/listMet*/*score.txt | awk '{print $1}' | sort | uniq -c | tail -n+2 >> ${path}/checkfilelinesChr${i}.txt
echo "Chr${i} Done"
done

