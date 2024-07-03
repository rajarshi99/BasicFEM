seq 2 7 | while read num
do
		echo $num
		python fem.py $num | tee -a log_runs.txt
done
