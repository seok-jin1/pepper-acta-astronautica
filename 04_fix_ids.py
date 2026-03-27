import pandas as pd

meta_path = "version-2_integrated/integrated_metadata.tsv"
df = pd.read_csv(meta_path, sep="	")

# 에러 메시지에서 확인된 정확한 ID 목록
missing_map = {
    'Q1_Seed1': 'Q1_Seed1_S8_L001',
    'Q1_Seed2': 'Q1_Seed2_S9_L001',
    'Q1_Seed3': 'Q1_Seed3_S10_L001',
    'Q2_Seed1': 'Q2_Seed1_S12_L001',
    'Q3_Seed1': 'Q3_Seed1_S16_L001',
    'Q3_Seed2': 'Q3_Seed2_S17_L001',
    'Q3_Seed3': 'Q3_Seed3_S18_L001',
    'Q4_Seed1': 'Q4_Seed1_S22_L001',
    'Q4_Seed2': 'Q4_Seed2_S23_L001',
    'Q4_Seed3': 'Q4_Seed3_S24_L001'
}

def final_fix(sid):
    return missing_map.get(sid, sid)

df['sample-id'] = df['sample-id'].apply(final_fix)
df.to_csv(meta_path, sep="	", index=False)
print("✅ 최종 Sample ID 보정 완료!")
