# 처음 miRNA base 와 fake miRNA 파일을 가져와서 사용할수 있는 csv 파일 형태로 정제
# fold 정보 추가

## 제 컴퓨터에는 따로 가상환경에 해당 패키지를 설치하여 사용하였습니다, 코드를 실행하여 확인해야 할 경우 아래 패키지등을 설치해야 합니다.
#pip install tqdm # for 문의 진행사항 확인 해주는 패키지
# pip install seqfold # RNA fold 정보의 hairpin위치와 붙은 염기 등등을 알려주는 seqfold모듈 설치
# pip install ViennaRNA # RNA fold 정보의 일렬 배열과 free energy 등을 알려주는 vienna fold 설치


from tqdm import tqdm  # 오래걸리는 for 문의 진행상황 확인을 위해
from seqfold import fold  # seqfold
import RNA  # vienna RNA
import csv  # csv파일로 최종 정제

# miRNA 데이터에서 인간의 데이터만을 남기기
oneRNA = []
hsamiRNA = []
with open('miRNA.dat', "r") as f:  # 파일 불러오기
    for i in tqdm(range(13000000)):
        # 각 줄을 불러옴 , 데이터는 약 1200만줄, 대량의 데이터라 readlines 로 읽어오는게 어려워 한줄씩 읽어오기 위해 for문을 이렇게 구성하였음
        line = f.readline()
        line = line.split("   ")  # 띄어쓰기 3 칸을 기준으로 좌측 헤더와 데이터를 분리
        if line[0] == 'XX\n':  # 필요없는 줄 제거
            continue
        oneRNA.append(line)  # 하나의 서열에 대한 데이터를 oneRNA 리스트에 담음
        if line[0] == "//\n":  ## 하나의 RNA 데이터가 끝남
            if ("Homo sapiens" in oneRNA[2][1]):  # 그 RNA 데이터가 인간 것인지 확인
                hsamiRNA.append(oneRNA)  # 인간 것이라면 새로운 리스트에 담음
                oneRNA = []  # 리스트 초기화
            else:
                oneRNA = []  # 사람의 데이터가 아니라면 그냥 리스트 초기화

# hsamiRNA 안에는 3중 리스트의 형태로 사람의 데이터만이 모여있음

# 사람 miRNA 의 folding data 만들어 csv파일로 저장하기
miRNA = []
oneRNA = []
sq = ""
for folder in tqdm(hsamiRNA):
    for line in folder:
        if (line[0] == "SQ"):  # 시퀸스 헤더 남김
            seq = []            #시퀸스 리스트를 만들어 시퀀스 정보를 모으려고 함
            seq.append("SQ")    # 헤더 추가
        elif (line[0] == ""):  # 실제 시퀀스정보를 시퀀스 리스트에 추가
            sq += line[1]
            sq = sq.strip().replace(" ", "").upper()  # seq 띄어쓰기 없애기, 대문자로 정렬
        elif (line[0] == '//\n'):  # 하나의 데이터가 끝남을 확인

            seq.append(sq)  # 시퀀스를 위 시퀀스 라인에 추가함
            oneRNA.append(seq)  # 시퀸스를 하나의 RNA 정보 리스트에 추가
            # 폴딩 시작 : seq fold / ViennaRNA 2가지 폴딩 패키지에서 나오는 정보를 이용할것임
            # vienna fold
            (ss, mfe) = RNA.fold(sq)  # Vienna fold 실행하여 fold 정보 저장
            oneRNA.append(["VF", ss])  # vienna fold정보 : )(. 3가지 형태로 시퀀스를 나타냄
            oneRNA.append(["FE", mfe])  # free energy : minimum Free Energy 를 반환함
            # seq fold
            folding = fold(sq)  # seq fold 실행

            for line in folding:  # 폴딩 데이터는 리스트 안의 스트럭처 형태로 반환됨
                line = str(line)  # 각 줄의 폴딩 데이터
                line = line.strip() #폴딩 데이터의 좌우 공백 제거
                line = line.split("   ") #폴딩 데이터를 각 폴딩의 위치 시퀀스 / 모양으로 분할함
                if len(line[0]) >= 4:  # 3자리수 이상일경우 split이 정상적으로 이루어지지 않음
                    miRNAsq = line[0].split("  ") # 3자리수 이상일때는 2개의 공백으로 나누기
                    line[0] = miRNAsq[0]          # 위 나누어진 데이터를 일반적인 2자리수
                    line.insert(1, miRNAsq[1])    # 데이터와 같은 리스트 인덱스 위치로 만들어줌
                line.insert(0, "FD") # 좌측에 헤더 삽입
                line.append("\n")
                oneRNA.append(line)
            sq = "" # 시퀀스 초기화
            seq.append("\n") # 마지막에 줄나누기 삽입
            oneRNA.append(["//", "\n"]) #하나의 데이터가 끝났음을 알리는 // 삽입
            miRNA.append(oneRNA)
            oneRNA = []
            seq = []
        else:
            continue

with open("fold_miRNA.csv", "w", newline="") as f:
    wr = csv.writer(f)
    for row in miRNA:
        for line in row:
            wr.writerow(line)

# miRNA 에는 전체 시퀀스, folding정보가 들어있음 ,
