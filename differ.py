# Author Daniel Hakim
# 10/28/2021

import pandas as pd

ms = pd.read_csv("MS assembled to NC_008261 Annotations_genes.csv")
hhc = pd.read_csv("HHC assembled to NC_008261 Annotations_genes.csv")

print(ms)
print(hhc)

# Appears the two csv files are sorted in descending order by column Minimum, the position at which the snp starts
# We can therefore iterate both lists at once in descending order

ms_index = 0
hhc_index = 0

ms_unmatched = 0
hhc_unmatched = 0
matches = 0

all_matches = []

# While there are remaining elements in both lists
while ms_index < ms.shape[0] and hhc_index < hhc.shape[0]:
    index_sum = ms_index + hhc_index
    if index_sum % 1000 == 0:
        print("Progress: " + str(int((index_sum / (ms.shape[0] + hhc.shape[0])) * 100)) + "%")
    # Get the current snp from each file
    ms_snp = ms.iloc[ms_index]
    hhc_snp = hhc.iloc[hhc_index]

    # if the files have the same minimum column, they might be the same, check them
    if ms_snp["Minimum"] == hhc_snp["Minimum"]:

        # While the files are sorted by Minimum, there can be multiple snps occurring at the same position
        # Find all rows in both files that have the same minimum, then we'll pairwise compare them
        ms_i = ms_index
        while ms.iloc[ms_i]["Minimum"] == ms_snp["Minimum"]:
            ms_i += 1

        hhc_i = hhc_index
        while hhc.iloc[hhc_i]["Minimum"] == hhc_snp["Minimum"]:
            hhc_i += 1

        if ms_i - 1 != ms_index or hhc_i - 1 != hhc_index:
            print("Comparing Ranges", ms_index, "-", ms_i, " vs ", hhc_index, "-", hhc_i)

        left_match_set = [0] * (ms_i - ms_index)
        right_match_set = [0] * (hhc_i - hhc_index)

        for left_index in range(ms_index, ms_i):
            for right_index in range(hhc_index, hhc_i):
                left = ms.iloc[left_index]
                right = hhc.iloc[right_index]
                if left['Maximum'] == right['Maximum'] and left['Change'] == right['Change']:
                    all_matches.append((left_index, right_index))
                    left_match_set[left_index - ms_index] = 1
                    right_match_set[right_index - hhc_index] = 1
                    matches += 1  # Hmm, is this a useful metric?

        ms_unmatched += len(left_match_set) - sum(left_match_set)
        hhc_unmatched += len(right_match_set) - sum(right_match_set)
        ms_index = ms_i
        hhc_index = hhc_i
    elif ms_snp["Minimum"] > hhc_snp["Minimum"]:
        # Since both files are descending, this ms_snp cannot be matched in the hhc file,
        # skip it and move on to the next one
        ms_index = ms_index + 1
        ms_unmatched += 1
    elif ms_snp["Minimum"] < hhc_snp["Minimum"]:
        # Since both files are descending, this hhc_snp cannot be matched in the ms file,
        # skip it and move on to the next one
        hhc_index = hhc_index + 1
        hhc_unmatched += 1
    else:
        # I don't think we can get here without parsing errors
        print("Huh?")
        print(ms_snp["Minimum"])
        print(hhc_snp["Minimum"])
        raise Exception("What is in this column??")

while ms_index < ms.shape[0]:
    # Leftover rows from the MS file
    ms_unmatched += 1
    ms_index += 1

while hhc_index < hhc.shape[0]:
    # Leftover rows from the hhc file
    hhc_unmatched += 1
    hhc_index += 1

print("MATCHES:", matches)
print("HHC_ONLY:", hhc_unmatched)
print("MS_ONLY:", ms_unmatched)

print(all_matches[:10])
ms_indexes = [x[0] for x in all_matches]
ms_same = ms.iloc[ms_indexes]

hhc_indexes = [x[1] for x in all_matches]
hhc_same = hhc.iloc[hhc_indexes]

ms_same.to_csv("ms_shared.csv")
hhc_same.to_csv("hhc_shared.csv")
with open("matches.txt", "w") as f:
    for m in all_matches:
        f.write(str(m[0] + 2) + ", " + str(m[1] + 2) + "\n")

ms_diff_set = [True] * ms.shape[0]
hhc_diff_set = [True] * hhc.shape[0]
for m in all_matches:
    ms_diff_set[m[0]] = False
    hhc_diff_set[m[1]] = False

ms_diff = ms[ms_diff_set]
hhc_diff = hhc[hhc_diff_set]

ms_diff.to_csv("ms_diff.csv")
hhc_diff.to_csv("hhc_diff.csv")

