import pandas as pd
import os

# 1. 기존 통합 메타데이터 로드
meta_path = "version-2_integrated/integrated_metadata.tsv"
df = pd.read_csv(meta_path, sep="	")

# 2. ID 변환 규칙 정의
# OSD-772 샘플들은 뒤에 시퀀싱 정보가 붙어야 함.
# 외부 데이터(Terrestrial_SRR...)는 그대로임.

def fix_id(sample_id):
    if sample_id.startswith("Terrestrial_"):
        return sample_id
    
    # OSD-772 샘플들에 대한 매핑 (파일명 리스트 osd_r1_list.txt 활용)
    return sample_id_map.get(sample_id, sample_id)

# 파일명에서 정확한 ID 매핑 생성
sample_id_map = {}
with open("osd_r1_list.txt", "r") as f:
    for line in f:
        path = line.strip()
        filename = os.path.basename(path)
        # GAmplicon_ 뒤부터 _R1_ 전까지가 QIIME2 인식 ID임
        # 예: GLDS-675_GAmplicon_Q1-1-Swab1_S1_L001_R1_raw.fastq.gz -> Q1-1-Swab1_S1_L001
        parts = filename.split("GAmplicon_")[1].split("_R1_")[0]
        original_id = parts.split("_S")[0] # Q1-1-Swab1
        sample_id_map[original_id] = parts

df['sample-id'] = df['sample-id'].apply(fix_id)

# 3. 저장
df.to_csv(meta_path, sep="	", index=False)
print(f"✅ 메타데이터 Sample ID 동기화 완료: {meta_path}")
