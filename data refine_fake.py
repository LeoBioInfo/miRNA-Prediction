from tqdm import tqdm  ### 오래걸리는 반복문의 진행시간을 알려줌

import csv
from seqfold import fold
from Bio import SeqIO  # biopython 모듈 불러오기
import RNA

result = []

for seq in tqdm(list(SeqIO.parse("pseudo_hairpin_human.fasta", "fasta"))):
    seq = str(seq.seq)
    seq = seq.replace("\\", "").replace("}", "")
    result.append(["SQ", seq])
    # vienna fold
    (ss, mfe) = RNA.fold(seq)  # Vienna fold 실행
    result.append(["VF", ss])  # vienna fold정보
    result.append(["FE", mfe])  # free energy

    # seq fold
    folding = fold(seq)
    for line in folding:
        line = str(line)
        line = line.strip().split("   ")
        line.insert(0, "FD")
        line.append("\n")
        if line[0] == "FD" and (len(line[1]) >= 4):  ## miRNA 3자리인 경우 스플릿 제대로 되지 않아 수정
            miRNAsq = line[1].split("  ")
            line[1] = miRNAsq[0]
            line.insert(2, miRNAsq[1])
            result.append(line)
        else:
            result.append(line)

    result.append(["//", "\n"])

with open('fold_pseudo_miRNA.csv', "w") as f:
    wr = csv.writer(f)
    for line in result:
        wr.writerow(line)
