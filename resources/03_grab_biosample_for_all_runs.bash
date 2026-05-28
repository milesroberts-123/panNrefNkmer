readarray -t runs < "runs_srr_only.txt"

sub=("${runs[@]:11515}")

for run in "${sub[@]}"; do
    echo "$run"
    [[ "$run" == "NA" ]] && continue
    esearch -db sra -query $run | efetch -format runinfo >> runinfo_all_runs_3.csv
    sleep 1
done
