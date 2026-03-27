import pandas as pd

# 1. OSD-772 메타데이터 로드
osd_meta = pd.read_csv("osd772_metadata.tsv", sep="	")
# 'spaceflight' 컬럼 기반으로 'Study_Group' 정의
def define_group(row):
    if row['spaceflight'] == 'Ground_Control':
        return 'Ground_Seed'
    else:
        return 'Space_Flight'

osd_meta['Study_Group'] = osd_meta.apply(define_group, axis=1)
osd_meta['Environment'] = osd_meta.apply(lambda x: 'Space_Hydroponic' if x['Study_Group'] == 'Space_Flight' else 'Terrestrial_Seed', axis=1)
osd_meta['Variety'] = 'NuMex_Trinidad_SC'

# 2. 외부 데이터(Terrestrial) 메타데이터 생성
ext_mapping = []
with open("external_mapping.txt", "r") as f:
    for line in f:
        run_acc, alias = line.strip().split("	")
        variety = alias.split("_")[0] # XY6, XY21 등 추출
        ext_mapping.append({
            'sample-id': f"Terrestrial_{run_acc}",
            'spaceflight': 'Ground_Control',
            'sample_type': 'Rhizosphere',
            'location': 'Earth_Organic_Farm',
            'sanitization': 'Not Applicable',
            'description': f'Terrestrial Soil from {variety}',
            'Study_Group': 'Terrestrial_Soil',
            'Environment': 'Terrestrial_Field',
            'Variety': variety
        })
ext_meta = pd.DataFrame(ext_mapping)

# 3. 통합
integrated_meta = pd.concat([osd_meta, ext_meta], ignore_index=True)

# 4. 저장
output_path = "version-2_integrated/integrated_metadata.tsv"
integrated_meta.to_csv(output_path, sep="	", index=False)

print(f"✅ 통합 메타데이터 생성 완료: {output_path}")
print(f"📊 총 샘플 수: {len(integrated_meta)}")
print(integrated_meta['Study_Group'].value_counts())
