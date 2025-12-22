# Simulating a bidirectional gene flow scenario and testing the bidirectional admixture proportion estimation for "Notes on f4-ratio estimation", Kalle Leppälä.
# This log records how the tree sequences created in bidirectional.slim are recapitulated with coalescent simulations,
# then equipped with random mutations, and combined into binary PLINK files of prefix "merged".
# To make pyslim, tskit and msprime to work in my Windows machine,
# before running the python code I type the following two lines in Windows Command Prompt (you have to find your own way):
# C:\Users\kalle\Miniconda3\Scripts\activate.bat
# conda activate

import pyslim, tskit, msprime, io

subset_inds = list(range(0, 10)) + list(range(10000, 10010)) + list(range(20000, 20010)) + list(range(30000, 30010)) + list(range(40000, 40010))

for i in range(0, 10):
    ts = tskit.load(f"chromosome_{i + 1}.trees")
    recap = pyslim.recapitate(ts, ancestral_Ne = 10000, recombination_rate = 1e-8, random_seed = 1)
    mutated = msprime.sim_mutations(recap, rate = 1e-8, random_seed = 1, keep = True)
    subset_nodes = [n for ind in subset_inds for n in mutated.individual(ind).nodes]
    subset_ts = mutated.simplify(subset_nodes)
    # While writing the VCF files I give custom names to variants, so that I can later merge the files without collisions. 
    buffer = io.StringIO()
    subset_ts.write_vcf(buffer, contig_id = f"chr{i + 1}")
    buffer.seek(0)
    with open(f"chr{i + 1}.vcf", "w") as out_vcf:
        for line in buffer:
            if line.startswith("#"):
                out_vcf.write(line)
            else:
                fields = line.strip().split("\t")
                pos = fields[1]
                fields[2] = f"chr{i + 1}_{pos}"
                out_vcf.write("\t".join(fields) + "\n")

# The file merge_list.txt contains nine rows, from "chr2" to "chr10" (without the quotation marks).
# In Windows Command Prompt I use PLINK to merge the ten data sets into one:
# for /L %i in (1, 1, 10) do plink --vcf chr%i.vcf --double-id --make-bed --out chr%i
# plink --bfile chr1 --merge-list merge_list.txt --make-bed --out merged
# This analysis is concluded in the file bidirectional.R (using admixtools to calculate the f-statistics necessary for admixture proportion estimation).