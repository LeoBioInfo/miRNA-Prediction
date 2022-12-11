import csv
from Bio.SeqUtils import GC  ## GC 계산 위해

onedata = []
dataline = [0, 0, 0, 0, 0, 0, 0, 0, 0]
fulldata = [["ISmiRNA", "MFEI", "NMFE", "sqlen", "GC", "hairpinMFE", "interior_loopMFE", "STACKMFE", "BULGEMFE"],
            ["d", "c", "c", "c", "c", "c", "c", "c", "c"], ["class", "", "", "", ""]]
with open("fold_miRNA.csv", "r") as f:  # fold 된 miRNA 데이터 불러오기
    data = csv.reader(f)
    for row in data:
        if row[0] == "FD" or row[0] == "SQ" or row[0] == "VF" or row[0] == "FE":
            onedata.append(row)
        elif row[0] == "//":  # 하나의 데이터 모음 끝남
            onedata.append(row)
            IRcount = 10
            BUcount = 10
            STcount = 1
            for line in onedata:
                if line[0] == "SQ":
                    fullseq = line[1]
                    sqlen = len(fullseq)  ## sq 길이
                    GCper = GC(fullseq)
                elif line[0] == "VF":
                    continue
                elif line[0] == "FE":
                    MFE = float(line[1])

                elif line[0] == "FD":
                    if ("INTERIOR_LOOP" in line[3]):
                        IRcount += 1
                    elif ("STACK" in line[3]):
                        STcount += 1
                    elif ("BULGE" in line[3]):
                        BUcount += 1

                    elif ("HAIRPIN" in line[3]):

                        hairpinlen = int(line[2])-int(line[1])
            # 순서 : "ISmiRNA", "MFEI", "NMFE", "sqlen", "GC", "hairpinMFE", "interior_loopMFE", "STACKMFE", "BULGEMFE"
            dataline[0] = "T"  # 맨앞:인덱스0 : miRNA 인지 아닌지 구별 : T이 miRNA :: ISmiRNA
            dataline[1] = MFE / GCper  # 인덱스 1: MFE/GC
            dataline[2] = MFE / sqlen  # 인덱스 2 : 표준화 MFE = MFE/sqlen
            dataline[3] = sqlen  # 인덱스 3 : 시퀀스 길이
            dataline[4] = GCper  # 인덱스 4 : GC 비율
            dataline[5] = MFE / hairpinlen  # 헤어핀 시퀀스 길이로 정규화
            dataline[6] = MFE / IRcount  # interior_loop 개수
            dataline[7] = MFE / STcount  # STACK 개수
            dataline[8] = MFE / BUcount  # BULGE 개수
            fulldata.append(dataline)
            onedata = []
            dataline = [0, 0, 0, 0, 0, 0, 0, 0, 0]
# 뒤에 fake data 붙이기
with open("fold_pseudo_miRNA.csv", "r") as f:
    data = csv.reader(f)
    for row in data:
        if row[0] == "FD" or row[0] == "SQ" or row[0] == "VF" or row[0] == "FE":
            onedata.append(row)
        elif row[0] == "//":  # 하나의 데이터 모음 끝남
            onedata.append(row)
            IRcount = 10
            BUcount = 10
            STcount = 1

            for line in onedata:
                if line[0] == "SQ":
                    fullseq = line[1]
                    sqlen = len(fullseq)  ## sq 길이
                    GCper = GC(fullseq)

                    dataline[2] = sqlen  # 인덱스 2: 시퀀스 길이
                elif line[0] == "VF":
                    continue
                elif line[0] == "FE":
                    MFE = float(line[1])

                elif line[0] == "FD":
                    if ("INTERIOR_LOOP" in line[3]):
                        IRcount += 1
                    elif ("STACK" in line[3]):
                        STcount += 1
                    elif ("BULGE" in line[3]):
                        BUcount += 1

                    elif ("HAIRPIN" in line[3]):

                        hairpinloc = int(line[1])  # 인덱스 4번 : 헤어핀 sq 위치
            # 순서 : "ISmiRNA", "MFEI", "NMFE", "sqlen", "GC", "hairpinMFE", "interior_loopMFE", "STACKMFE", "BULGEMFE"
            dataline[0] = "F"  # 맨앞:인덱스0 : miRNA 인지 아닌지 구별 : T이 miRNA :: ISmiRNA
            dataline[1] = MFE / GCper  # 인덱스 1: MFE/GC
            dataline[2] = MFE / sqlen  # 인덱스 2 : 표준화 MFE = MFE/sqlen
            dataline[3] = sqlen  # 인덱스 3 : 시퀀스 길이
            dataline[4] = GCper  # 인덱스 4 : GC 비율
            dataline[5] = MFE / hairpinlen  # 헤어핀 시퀀스 길이로 정규화
            dataline[6] = MFE / IRcount  # interior_loop 개수
            dataline[7] = MFE / STcount  # STACK 개수
            dataline[8] = MFE / BUcount  # BULGE 개수
            fulldata.append(dataline)
            onedata = []
            dataline = [0, 0, 0, 0, 0, 0, 0, 0, 0]

with open("is pre-miRNA.tab", "w", newline="") as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(fulldata)
