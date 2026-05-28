
readarray -t lines < "biosample_ids.txt"

for line in "${lines[@]}"; do
    echo "$line"    
    SRR=$(elink -db biosample -id $line -target sra | efetch -format docsum | xtract -pattern Runs -element Run@acc | tr '\n' ',')
    echo "$line: $SRR" >> sample_to_run_mapping.txt
    sleep 1 # Respect NCBI servers to prevent IP blocks
done

sed 's/.*: //g' sample_to_run_mapping.txt | tr ',' '\n' | tr '\t' '\n' | grep . > runs.txt

readarray -t runs < "runs.txt"

for run in "${runs[@]}"; do
    echo "$run"
    esearch -db sra -query $run | efetch -format runinfo >> runinfo.csv
    sleep 1
done
