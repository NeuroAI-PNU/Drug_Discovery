import numpy as np
import pandas as pd



if __name__=="__main__":
    for i in range(40):
        data = open(f"./docking_data/docking_score_based_tanimoto_similarity/docking_score_{i}.txt", "r", encoding="UTF8")
        lines = data.readlines()
        score = []
        for x in range(len(lines)):
            if "mode" in lines[x]:
                if lines[x+4][0] == "2":
                    score.append(lines[x+4][8:12])
                elif lines[x+4][0] != "2":
                    score.append(lines[x+3][8:12])

        if i==0:
            merged_data = score
        else:
            merged_data = merged_data + score

        print(len(score))
        merged_docking_data = pd.DataFrame(merged_data)
        
    merged_docking_data.to_csv("./docking_data/parsed_docking_data_tanimoto_similarity/merged_docking_data_based_tanimoto_similarity.csv", index=False)