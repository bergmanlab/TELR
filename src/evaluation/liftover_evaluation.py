#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import json
import glob
import numpy as np
import matplotlib.pyplot as plt

def evaluate_liftover(json_file):

    with open(json_file) as data:
        data = json.load(data)
    for item in data:
        if item["num_hits"] == 1:
            print("yes")
        break

    df = pd.read_json(json_file, orient='columns')

    df.hist(column='num_hits', bins = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
    df['num_hits'].value_counts()

    # generate summary
    nonref_total = 0
    ref_total = 0
    unlift_total = 0
    nonref_comments = dict()
    ref_comments = dict()
    unlift_comments = dict()
    for item in data:
        if item["num_hits"] == 1:
            nonref_total = nonref_total + 1
        else:
            infos = item["report"]
            if len(infos) > 1:
                print(item["ID"])
            for info in infos:
                if info["type"] == "non-reference":
                    if "comment" in info.keys():
                        if info["comment"] not in nonref_comments.keys():
                            nonref_comments[info["comment"]] = 1
                        else:
                            nonref_comments[info["comment"]] = nonref_comments[info["comment"]] + 1
                if info["type"] == "reference":
                    ref_total = ref_total + 1
                    if "comment" in info.keys():
                        if info["comment"] not in ref_comments.keys():
                            ref_comments[info["comment"]] = 1
                        else:
                            ref_comments[info["comment"]] = ref_comments[info["comment"]] + 1
                if info["type"] == "unlifted":
                    unlift_total = unlift_total + 1
                    if "comment" in info.keys():
                        if info["comment"] not in unlift_comments.keys():
                            unlift_comments[info["comment"]] = 1
                        else:
                            unlift_comments[info["comment"]] = unlift_comments[info["comment"]] + 1
    summary_data = {
        "non-reference": {"total": nonref_total, "comments": nonref_comments},
        "reference": {"total": ref_total, "comments": ref_comments},
        "unlifted": {"total": unlift_total, "comments": unlift_comments}
    }

    summary_report = "/scratch/sh60271/s2rplus/liftover_0627/out/liftover_summary.json"
    with open(summary_report, "w") as output:
        json.dump(summary_data, output, indent=4, sort_keys=False)

    df_new = df.loc[df['num_hits'] == 0]
    result = df_new.to_dict(orient="records")
    json_out = "/scratch/sh60271/s2rplus/liftover_0622/out/num_hit_zero.json"
    with open(json_out, "w") as output:
        json.dump(result, output, indent=4, sort_keys=False)


if __name__ == '__main__':
    evaluate_liftover()