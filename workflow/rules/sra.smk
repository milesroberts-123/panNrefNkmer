rule sra:
	output:
		temp("raw_reads/{ID}_1.fastq.gz"),
		temp("raw_reads/{ID}_2.fastq.gz")
	log:
		"logs/sra/{ID}.log"
	conda:
		"../envs/sra.yaml"
	shell:
		"""
		# download data
		fastq-dump --gzip --split-e {wildcards.ID} &> {log}
		
		# move data to folder
		mv {wildcards.ID}_1.fastq.gz raw_reads/
		mv {wildcards.ID}_2.fastq.gz raw_reads/
		"""
