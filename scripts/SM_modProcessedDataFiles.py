"""
This script takes information from fdom_lineage .csv data file and appends it appropriately to
fdom .csv data file.

For now:
 * Add m50_phenotype_score
 * Add m100_phenotype_score
 * Add m1000_phenotype_score
"""
import os, csv

collect = [-1000]
question_info = {
    "Q1": {
        "final_update": 200000,
        "grab_scores_at": [200000 + val for val in collect],
    },
    "Q2": {
        "final_update": 200000,
        "grab_scores_at": [200000 + val for val in collect]
    },
    "Q3": {
        "final_update": 400000,
        "grab_scores_at": [400000 + val for val in collect]
    }
}

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
    fdom_data_fname = "extra_runs__slipmuts_iter2_fdom.csv"
    relevant_lin_data = {}

    # 1) Extract info we want from lineage time series.
    with open(os.path.join(processed_dir, lin_data_fname), "r") as fp:
        reader = csv.reader(fp, delimiter = ",")
        header = reader.next()
        col_lookup = {header[i]: i for i in range(0, len(header))}
        for row in reader:
            u = int(row[col_lookup["update"]])
            q = row[col_lookup["treatment"]][:2]
            rep = row[col_lookup["rep"]]
            if not rep in relevant_lin_data: relevant_lin_data[rep] = {}
            if u in question_info[q]["grab_scores_at"]:
                relevant_lin_data[rep][u] = row[col_lookup["score"]]

    content = ""
    with open(os.path.join(processed_dir, fdom_data_fname), "r") as fp:
        reader = csv.reader(fp, delimiter = ",")
        header = reader.next()
        col_lookup = {header[i]: i for i in range(0, len(header))}
        # New header
        content += ",".join(header) + "," + ",".join(["%d_phenotype_score" % val for val in collect]) + ",min_genome_size,fancy_name" + "\n"
        for row in reader:
            rep = row[col_lookup["rep_id"]]
            q = row[col_lookup["treatment"]][:2]
            treat = row[col_lookup["treatment"]]
            op = treat.split("__")[1]
            mingen = treat.split("__")[-1].split("_")[-1]
            # Get rep
            content += ",".join(map(str, row)) + "," + ",".join([str(relevant_lin_data[rep][question_info[q]["final_update"] + val]) for val in collect]) + "," + mingen + "," + treatment_map[op] + "\n"

    with open("MODIFIED__" + fdom_data_fname, "w") as fp:
        fp.write(content)

if __name__ == "__main__":
    main()
