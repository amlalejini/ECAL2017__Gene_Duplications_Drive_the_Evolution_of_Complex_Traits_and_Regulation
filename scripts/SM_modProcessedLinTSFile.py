"""
This script takes information from fdom_lineage .csv data file and appends it appropriately to
fdom .csv data file.

For now:
 * Add m50_phenotype_score
 * Add m100_phenotype_score
 * Add m1000_phenotype_score
"""
import os, csv
from sets import Set

treatment_map = {
    "Slip_DUP": "Slip-duplicate",
    "Slip_SCRAM": "Slip-scramble",
    "Slip_RAND": "Slip-random",
    "Slip_NOP": "Slip-NOP",
    "Slip_SCATTER": "Slip-scatter",
    "HIGHMUT": "High mutation rate",
    "Slip_NONE": "Baseline"
}


def main():
    # Some settings
    exp_base_dir = "/Users/amlalejini/DataPlayground/slip_muts/iter_2"
    processed_dir = os.path.join(exp_base_dir, "extra_data_processed")
    lin_data_fname = "extra_runs_lineage_score_ts.csv"
    content = ""
    # Check uniqueness of each line.
    unique = Set()
    with open(os.path.join(processed_dir, lin_data_fname), "r") as fp:
        lines = fp.read().split("\n")
        content += lines[0] + ",fancy_name\n"
        for i in range(1, len(lines)):
            line = lines[i]
            time = line.split(",")[0]
            rep = line.split(",")[-2]
            unique_name = time + "___" + rep
            fancy_name = treatment_map[line.split(",")[1].split("__")[1]]
            if unique_name in unique:
                print "Oh no! %s" % unique_name
                continue
            unique.add(unique_name)
            content += "%s,%s\n" % (line, fancy_name)
    with open("MODIFIED__" + lin_data_fname, "w") as fp:
        fp.write(content)

if __name__ == "__main__":
    main()
