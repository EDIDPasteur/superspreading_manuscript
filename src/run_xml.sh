#!/usr/bin/bash

#SBATCH -q long
#SBATCH -p long
#SBATCH --ntasks 1
#SBATCH -N 1
#SBATCH --job-name="Run xml"
#SBATCH --output %A_%a_BEAST_cs_sim.out
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --time 7-00:00:00

module load beagle-lib/ beast/v2.7.3
module load R

export _JAVA_OPTIONS="-Xmx100G"

FILES=($1/*.xml)
INPUTFILE=${FILES[$SLURM_ARRAY_TASK_ID]}


get_min_ess() {
	log_file="$1"
	if [[ ! -f "$log_file" ]]; then
		echo 0
		return
	fi

	ess_values=$(loganalyser -t "$log_file" 2>/dev/null \
        | awk 'NR > 1 {print $9}' \
        | grep -Eo '^[0-9]+(\.[0-9]+)?$')

	if [[ -z "$ess_values" ]]; then
		echo 0
		return
	fi
	min_ess=$(echo "$ess_values" | sort -n | head -n 1)
	echo "$min_ess"
}


echo "Starting BEAST for $INPUTFILE"
base="${INPUTFILE%.xml}"
log_name="${base}.log"

beast -threads 2 -overwrite -working $INPUTFILE & pid=$!

while true; do
	sleep 900
	ess=$(get_min_ess "$log_name")
	echo "Current lowest ESS for $INPUTFILE: $ess"

	if (( $(echo "$ess >= 200" | bc -l) )); then
		echo "All ESS >= 200. Stopping BEAST for $INPUTFILE"
		kill "$pid"
		wait "$pid" 2>/dev/null
		break
	fi
done


