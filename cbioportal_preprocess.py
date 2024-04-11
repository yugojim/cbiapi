# 匯入套件
import json
import pandas as pd
import os
import csv
import codecs
from seqslab import hive
import requests
import sys
import numpy as np
import re

# 開始進行資料清理
print('--------------------- Start cBio Portal Preprocess. Clean and Generate Data for uploading to cBio Portal ---------------------')

# 確認目前 python 環境是否正確
print(f'Resent working environment is: {sys.executable}\n')

# 定義會使用到的函數
# 定義 json 轉成 dataframe 的函數
def json_to_dataframe (json_file_path):
    # 建立一個空列表來儲存JSON資料
    json_data_list = []
    # 打開 JSON 檔案並逐行讀取其中的內容
    with open(json_file_path, 'r', encoding='utf-8') as json_file:
        for line in json_file:
            json_data = json.loads(line)
            json_data_list.append(json_data)
    return pd.DataFrame(json_data_list)

# 對 dataframe 做基本分析
def analyze_dataframe(dataframe, dataframe_name=None):
    # 所有欄位名稱
    column_names = dataframe.columns.tolist()
    # 資料筆數
    data_count = len(dataframe)
    # 每個欄位的資料類型
    data_types = dataframe.dtypes
    # 是否有缺失值
    has_missing_values = dataframe.isnull().any()
    # 每個欄位的缺失值數量
    missing_value_counts = dataframe.isnull().sum()
    # 總共有多少筆缺失值
    total_missing_values = dataframe.isnull().sum().sum()

    if dataframe_name:
        print(f"{dataframe_name} 所有欄位名稱:")
    else:
        print("所有欄位名稱:")
    print(column_names)
    print("----------------------------------------")
    
    if dataframe_name:
        print(f"{dataframe_name} 資料筆數:")
    else:
        print("資料筆數:")
    print(data_count)
    print("----------------------------------------")
    
    if dataframe_name:
        print(f"{dataframe_name} 每個欄位的資料類型:")
    else:
        print("每個欄位的資料類型:")
    print(data_types)
    print("----------------------------------------")
    
    if dataframe_name:
        print(f"{dataframe_name} 是否有缺失值:")
    else:
        print("是否有缺失值:")
    print(has_missing_values)
    print("----------------------------------------")
    
    if dataframe_name:
        print(f"{dataframe_name} 每個欄位的缺失值數量:")
    else:
        print("每個欄位的缺失值數量:")
    print(missing_value_counts)
    print("----------------------------------------")
    
    if dataframe_name:
        print(f"{dataframe_name} 總共有多少筆缺失值:")
    else:
        print("總共有多少筆缺失值:")
    print(total_missing_values)

# 定義函數以處理 'MAF' 列，如果 list 中沒有值，則把 'MAF' 的欄位改成 None
def replace_empty_list_with_none(maf_list):
    if not maf_list:
        return None
    return maf_list

# 定義針對一些特殊情況做處理的函數
def preprocess_tumor_sample_barcode(df):
    # 替換空格和特定字符串
    df['Tumor_Sample_Barcode'] = df['Tumor_Sample_Barcode'].str.replace(' ', '_')
    df['Tumor_Sample_Barcode'] = df['Tumor_Sample_Barcode'].str.replace('修正報告', '_revised')
    df['Tumor_Sample_Barcode'] = df['Tumor_Sample_Barcode'].str.replace('拷貝', '_copy')
    df['Tumor_Sample_Barcode'] = df['Tumor_Sample_Barcode'].str.replace('初步', '_preliminary')
    df['Tumor_Sample_Barcode'] = df['Tumor_Sample_Barcode'].str.replace('_v1_04d005dc-90d7-47ee-8252-0d9e0d45ecb7_', '_')
    # 替換指定字串
    replacement_string = 'TBB'
    df['Tumor_Sample_Barcode'] = df['Tumor_Sample_Barcode'].str.replace('Exclude_variant_in_Taiwan_BioBank_with_over_1%_allele_frequency', replacement_string)
    df['Tumor_Sample_Barcode'] = df['Tumor_Sample_Barcode'].str.replace('Exclude_variant_in_Taiwan_BioBank_with_over_1%_allele_frequency.', replacement_string)
    df['Tumor_Sample_Barcode'] = df['Tumor_Sample_Barcode'].str.replace('Exclude_variant_in_Taiwan_BioBank_with_over_1__allele_frequency', replacement_string)
    df['Tumor_Sample_Barcode'] = df['Tumor_Sample_Barcode'].str.replace('Exclude_variant_in_Taiwan_BioBank_with_over_1__allele_frequency.', replacement_string)
    df['Tumor_Sample_Barcode'] = df['Tumor_Sample_Barcode'].str.replace('Exclude_variant_in_Taiwan_BioBank_with_over_1__AF', replacement_string)
    df['Tumor_Sample_Barcode'] = df['Tumor_Sample_Barcode'].str.replace('.', '')
    # 使用str.slice截斷字串
    max_length = 40
    df['Tumor_Sample_Barcode'] = df['Tumor_Sample_Barcode'].str.slice(0, max_length)
    return df

# 定義針對一些特殊情況做處理的函數
def preprocess_patient_id(df):
    # 替換空格和特定字符串
    df['PATIENT_ID'] = df['PATIENT_ID'].str.replace(' ', '_')
    df['PATIENT_ID'] = df['PATIENT_ID'].str.replace('修正報告', '_revised')
    df['PATIENT_ID'] = df['PATIENT_ID'].str.replace('拷貝', '_copy')
    df['PATIENT_ID'] = df['PATIENT_ID'].str.replace('初步', '_preliminary')
    df['PATIENT_ID'] = df['PATIENT_ID'].str.replace('_v1_04d005dc-90d7-47ee-8252-0d9e0d45ecb7_', '_')
    # 替換指定字串
    replacement_string = 'TBB'
    df['PATIENT_ID'] = df['PATIENT_ID'].str.replace('Exclude_variant_in_Taiwan_BioBank_with_over_1%_allele_frequency', replacement_string)
    df['PATIENT_ID'] = df['PATIENT_ID'].str.replace('Exclude_variant_in_Taiwan_BioBank_with_over_1%_allele_frequency.', replacement_string)
    df['PATIENT_ID'] = df['PATIENT_ID'].str.replace('Exclude_variant_in_Taiwan_BioBank_with_over_1__allele_frequency', replacement_string)
    df['PATIENT_ID'] = df['PATIENT_ID'].str.replace('Exclude_variant_in_Taiwan_BioBank_with_over_1__allele_frequency.', replacement_string)
    df['PATIENT_ID'] = df['PATIENT_ID'].str.replace('Exclude_variant_in_Taiwan_BioBank_with_over_1__AF', replacement_string)
    df['PATIENT_ID'] = df['PATIENT_ID'].str.replace('.', '')
    # 使用str.slice截斷字串
    max_length = 40
    df['PATIENT_ID'] = df['PATIENT_ID'].str.slice(0, max_length)
    return df

# 定義處理 'HGVSp_Short' 這個欄位的函數
def clean_HGVSp_Short(df):
    # 將加號 '+' 替換為 ' plus '
    df['HGVSp_Short'] = df['HGVSp_Short'].str.replace('+', ' plus ')
    # 把欄位裡面的值都轉換成字串
    df['HGVSp_Short'] = df['HGVSp_Short'].astype(str)
    all_strings = empty_df['HGVSp_Short'].apply(lambda x: isinstance(x, str)).all()
    print("All values in 'HGVSp_Short' are converted to strings:", all_strings)
    # 清洗'HGVSp_Short'欄位的資料，只保留左右小括弧中間的字串
    df['HGVSp_Short'] = df['HGVSp_Short'].apply(lambda x: re.search(r'\((.*?)\)', x).group(1) if re.search(r'\((.*?)\)', x) else x)
    # 清洗'HGVSp_Short'欄位的資料，只保留左右中括弧中間的字串，中括弧裡面若有分號 ';' 則用空格取代
    df['HGVSp_Short'] = df['HGVSp_Short'].apply(lambda x: re.search(r'\[(.*?)\]', x).group(1) if re.search(r'\[(.*?)\]', x) else x)
    # 將分號 ';' 替換為空格 ' '
    df['HGVSp_Short'] = df['HGVSp_Short'].str.replace(';', ' ')
    return df

# 定義函數來匯出報告資料是否有缺漏，特別是針對 maf 中的 mutations 資料的檢查結果
def export_dataframe_to_files(index_list, patient_id_list, filename_prefix, output_folder, file_formats=['csv', 'xlsx']):
    # 建立 DataFrame 用來儲存資料
    df_export = pd.DataFrame({'Index': index_list, 'Patient_id': patient_id_list})
    
    for file_format in file_formats:
        # 匯出文件
        if file_format == 'csv':
            file_extension = 'csv'
        elif file_format == 'xlsx':
            file_extension = 'xlsx'
        
        file_name = f"{filename_prefix}_in_maf.{file_extension}"
        file_path = f"{output_folder}/{file_name}"
        
        if file_format == 'csv':
            df_export.to_csv(file_path, index=False)
        elif file_format == 'xlsx':
            df_export.to_excel(file_path, index=False)

# 定義一個函數來針對 cancer type 執行提取操作，cancer type 放在第一個逗點之後
def extract_cancer_type(x):
    if x is not None:
        parts = x.split(',')
        if len(parts) > 1:
            return parts[1].strip()
    return 'NA'

# 定義一個函數來自動精簡化 'REPORT_TEST_ASSAY' 欄位的名稱
def simplify_report_names(df, assay_list):
    for index, value in enumerate(df["REPORT_TEST_ASSAY"]):
        lower_value = value.lower()
        matching_reports = [report for report in assay_list if report.lower() in lower_value]
        if matching_reports:
            df.at[index, "REPORT_TEST_ASSAY"] = matching_reports[0]
    return df

# 從 Seqslab 撈取需要用到的資料
print('------------------------------------ Start to Download Data from Seqslab ------------------------------------')
# 指定要下载的文件的URL和文件名
file_info = [
    {'url': 'https://vghtpe685c2storage.blob.core.windows.net/seqslab/drs/usr_admin_vghtpe/fTJGZDsbPR74YvfR4ZJMs7/vghtpe_report_1_maf_icd10_icd_o.json?st=2023-10-30T01%3A31%3A47Z&se=2023-11-30T00%3A00%3A00Z&sp=rle&spr=https&sv=2020-06-12&sr=b&sig=43w8YchneqGXnDvAd4AqvsPruT5u7yL6Gi1uXvRLjhY%3D',
     'save_name': 'vghtpe_report_1_maf_icd10_icd_o.json'},
    {'url': 'https://vghtpe685c2storage.blob.core.windows.net/seqslab/drs/usr_admin_vghtpe/fTJGZDsbPR74YvfR4ZJMs7/vghtpe_report_2_maf_icd10_icd_o.json?st=2023-10-30T01%3A31%3A47Z&se=2023-11-30T00%3A00%3A00Z&sp=rle&spr=https&sv=2020-06-12&sr=b&sig=BLjP4LBmf6pN0cWIiCETD4WxDhxQr0M1jFKVJr4sdZE%3D',
     'save_name': 'vghtpe_report_2_maf_icd10_icd_o.json'}
]

# 指定文件的保存路径
save_path = 'path_to_save_file'  # 替換為實際要保存文件的路徑(從 seqslab 下載下來的 json 檔存放的路徑)

# 循環迭代 url 列表以下載 .json 文件
for file in file_info:
    file_url = file['url']
    save_name = file['save_name']
    save_path_full = os.path.join(save_path, save_name)

    # 取得要儲存文件的目錄
    save_dir = os.path.dirname(save_path_full)

    # 檢查該目錄是否存在，若不存在則創建該目錄
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    # 發送 GET 請求並下載文件
    response = requests.get(file_url)

    # 檢查響應狀態碼
    if response.status_code == 200:
        with open(save_path_full, 'wb') as file:
            file.write(response.content)
        print(f'document is successfully download to : {save_path_full}')
        # print('\n') # 空行
    else:
        print(f'failed to download: {response.status_code}')
        # print('\n') # 空行

# 要讀取的JSON文件路徑列表
json_file_paths = [
    'path_to_save_file/vghtpe_report_1_maf_icd10_icd_o.json',
    'path_to_save_file/vghtpe_report_2_maf_icd10_icd_o.json'
]

# 使用字典來儲存 dataframes
dataframes = {}

for i, json_file_path in enumerate(json_file_paths):
    dataframe_name = f'df_report_{i+1}'  # 使用i來生成dataframe名稱
    dataframes[dataframe_name] = json_to_dataframe(json_file_path)

# 存取 dataframes
df_report_1 = dataframes['df_report_1']
df_report_2 = dataframes['df_report_2']

# 合併兩個 DataFrame，將缺失值填充為 'NA'
dataframes = {'df_report_1': df_report_1, 'df_report_2': df_report_2}
concat_df = pd.concat(dataframes, axis=0, join='outer', ignore_index=True, sort=False)
# # 確保缺失的值為 'NA'
# concat_df.fillna('NA', inplace=True)

# 指定最後要拿來做資料處理的 dataframe
df_report = concat_df # 從這裡去修改要用哪一張表來製作上傳 cBio Portal 的資料

# 印出兩筆資料(兩堆報告)的欄位
print('-----------------------------------------------------------------------------------------')
print(f'df_report_1 columns: {df_report_1.columns}')
print(f'df_report_1 shape: {df_report_1.shape}')
print('-----------------------------------------------------------------------------------------')
print(f'df_report_2 columns: {df_report_2.columns}')
print(f'df_report_2 shape: {df_report_2.shape}')
print('-----------------------------------------------------------------------------------------')
print(f'concat_df columns: {concat_df.columns}')
print(f'concat_df shape: {concat_df.shape}')
print('-----------------------------------------------------------------------------------------')
print('\n') # 空行

# 檢查資料總表的狀況
analyze_dataframe(df_report, dataframe_name='df_report(從 Seqslab 抓取到的資料總表)')

# 檢查是否有重複的 label
duplicates = df_report[df_report.label.duplicated()]
if duplicates.empty:
    print("No duplicate label in df_report.\n")
else:
    print("There are duplicate labels in df_report：")
    print(duplicates)
    print('\n') # 空行

print('------------------------------------ Finish Downloading Data from Seqslab ------------------------------------')
print('\n') # 空行
print('\n') # 空行

# 讀取 config.json 裡面定義的變數(參數)設定，包含資料讀取與儲存路徑
# 讀取JSON文件
with open('config.json', 'r') as config_file:
    config = json.load(config_file)

# 定義資料讀取與儲存路徑
output_folder = config["output_folder"]
output_folder_case = config["output_folder_case"]
check_data_folder = config["check_data_folder"]

# 定義 meta 檔裡面需要的變數 (study)
type_of_cancer = config["type_of_cancer"]
cancer_study_identifier = config["cancer_study_identifier"]
name = config["name"]
description = config["description"]
add_global_case_list = config["add_global_case_list"]
groups = config["groups"]

# 定義 meta 檔裡面需要的變數 (mutations)
genetic_alteration_type_mu = config["genetic_alteration_type_mu"]
stable_id_mu = config["stable_id_mu"]
datatype_mu = config["datatype_mu"]
show_profile_in_analysis_tab = config["show_profile_in_analysis_tab_mu"]
profile_name_mu = config["profile_name_mu"]
profile_description = config["profile_description_mu"]
data_filename_mu = config["data_filename_mu"]

# 定義 meta 檔裡面需要的變數 (patient)
genetic_alteration_type_pa = config["genetic_alteration_type_pa"]
datatype_pa = config["datatype_pa"]
data_filename_pa = config["data_filename_pa"]

# 定義 meta 檔裡面需要的變數 (sample)
genetic_alteration_type_sa = config["genetic_alteration_type_sa"]
datatype_sa = config["datatype_sa"]
data_filename_sa = config["data_filename_sa"]

# 定義 meta 檔裡面需要的變數 (cna)
genetic_alteration_type_cna = config["genetic_alteration_type_cna"]
datatype_cna = config["datatype_cna"]
stable_id_cna = config["stable_id_cna"]
show_profile_in_analysis_tab = config["show_profile_in_analysis_tab_cna"]
profile_name_cna = config["profile_name_cna"]
profile_description_cna = config["profile_description_cna"]
data_filename_cna = config["data_filename_cna"]

# 定義 meta 檔裡面需要的變數 (sv)
genetic_alteration_type_sv = config["genetic_alteration_type_sv"]
datatype_sv = config["datatype_sv"]
stable_id_sv = config["stable_id_sv"]
show_profile_in_analysis_tab = config["show_profile_in_analysis_tab_sv"]
profile_name_sv = config["profile_name_sv"]
profile_description_sv = config["profile_description_sv"]
data_filename_sv = config["data_filename_sv"]

# 定義 meta 檔裡面需要的變數 (cancer type)
genetic_alteration_type_ctype = config["genetic_alteration_type_ctype"]
datatype_ctype = config["datatype_ctype"]
data_filename_ctype = config["data_filename_ctype"]

# 如果資料夾(路徑)不存在，則創建資料夾
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
# 如果資料夾(路徑)不存在，則創建資料夾
if not os.path.exists(output_folder_case):
    os.makedirs(output_folder_case)
# 如果資料夾(路徑)不存在，則創建資料夾
if not os.path.exists(check_data_folder):
    os.makedirs(check_data_folder)

# 自動生成 meta 檔
print('------------------------- Start to Check Meta Data is Successfully Generated or Not -------------------------')
# 指定输出文件路徑(meta_study)
output_file_meta_study = "meta_study.txt"
# 拼接完整的路徑名稱
output_path = os.path.join(output_folder, output_file_meta_study)
# 寫到 meta
with open(output_path, "w", encoding="utf-8") as file:
    file.write(f"type_of_cancer: {type_of_cancer}\n")
    file.write(f"cancer_study_identifier: {cancer_study_identifier}\n")
    file.write(f"name: {name}\n")
    file.write(f"description: {description}\n")
    file.write(f"add_global_case_list: {add_global_case_list}\n")
    file.write(f"groups: {groups}\n")
print(f"meta study is saved to '{output_path}'")
# print('\n') # 空行

# 指定输出文件路徑(meta_mutations_extended)
output_file_mutations_extended = "meta_mutations_extended.txt"
# 拼接完整的路徑名稱
output_path = os.path.join(output_folder, output_file_mutations_extended)
# 寫到 meta
with open(output_path, "w", encoding="utf-8") as file:
    file.write(f"cancer_study_identifier: {cancer_study_identifier}\n")
    file.write(f"genetic_alteration_type: {genetic_alteration_type_mu}\n")
    file.write(f"stable_id: {stable_id_mu}\n")
    file.write(f"datatype: {datatype_mu}\n")
    file.write(f"show_profile_in_analysis_tab: {show_profile_in_analysis_tab}\n")
    file.write(f"profile_description: {profile_description}\n")
    file.write(f"profile_name: {profile_name_mu}\n")
    file.write(f"data_filename: {data_filename_mu}\n")
print(f"meta mutations extended is saved to '{output_path}'")
# print('\n') # 空行

# 指定输出文件路徑(meta_clinical_patient)
output_file_clinical_patient = "meta_clinical_patient.txt"
# 拼接完整的路徑名稱
output_path = os.path.join(output_folder, output_file_clinical_patient)
# 寫到 meta
with open(output_path, "w", encoding="utf-8") as file:
    file.write(f"cancer_study_identifier: {cancer_study_identifier}\n")
    file.write(f"genetic_alteration_type: {genetic_alteration_type_pa}\n")
    file.write(f"datatype: {datatype_pa}\n")
    file.write(f"data_filename: {data_filename_pa}\n")
print(f"meta clinical patient is saved to '{output_path}'")
# print('\n') # 空行

# 指定输出文件路徑(meta_clinical_sample)
output_file_clinical_sample = "meta_clinical_sample.txt"
# 拼接完整的路徑名稱
output_path = os.path.join(output_folder, output_file_clinical_sample)
# 寫到 meta
with open(output_path, "w", encoding="utf-8") as file:
    file.write(f"cancer_study_identifier: {cancer_study_identifier}\n")
    file.write(f"genetic_alteration_type: {genetic_alteration_type_sa}\n")
    file.write(f"datatype: {datatype_sa}\n")
    file.write(f"data_filename: {data_filename_sa}\n")
print(f"meta clinical sample is saved to '{output_path}'")
# print('\n') # 空行

# 指定输出文件路徑(meta_cna)
output_file_cna = "meta_cna.txt"
# 拼接完整的路徑名稱
output_path = os.path.join(output_folder, output_file_cna)
# 寫到 meta
with open(output_path, "w", encoding="utf-8") as file:
    file.write(f"cancer_study_identifier: {cancer_study_identifier}\n")
    file.write(f"genetic_alteration_type: {genetic_alteration_type_cna}\n")
    file.write(f"datatype: {datatype_cna}\n")
    file.write(f"stable_id: {stable_id_cna}\n")
    file.write(f"show_profile_in_analysis_tab: {show_profile_in_analysis_tab}\n")
    file.write(f"profile_name: {profile_name_cna}\n")
    file.write(f"profile_description: {profile_description_cna}\n")
    file.write(f"data_filename: {data_filename_cna}\n")
print(f"meta cna is saved to '{output_path}'")
# print('\n') # 空行

# 指定输出文件路徑(meta_sv)
output_file_sv = "meta_sv.txt"
# 拼接完整的路徑名稱
output_path = os.path.join(output_folder, output_file_sv)
# 寫到 meta
with open(output_path, "w", encoding="utf-8") as file:
    file.write(f"cancer_study_identifier: {cancer_study_identifier}\n")
    file.write(f"genetic_alteration_type: {genetic_alteration_type_sv}\n")
    file.write(f"datatype: {datatype_sv}\n")
    file.write(f"stable_id: {stable_id_sv}\n")
    file.write(f"show_profile_in_analysis_tab: {show_profile_in_analysis_tab}\n")
    file.write(f"profile_name: {profile_name_sv}\n")
    file.write(f"profile_description: {profile_description_sv}\n")
    file.write(f"data_filename: {data_filename_sv}\n")
print(f"meta sv is saved to '{output_path}'")
# print('\n') # 空行

# 指定输出文件路徑(meta_cancer_type)
output_file_ctype = "meta_cancer_type.txt"
# 拼接完整的路徑名稱
output_path = os.path.join(output_folder, output_file_ctype)
# 寫到 meta
with open(output_path, "w", encoding="utf-8") as file:
    file.write(f"genetic_alteration_type: {genetic_alteration_type_ctype}\n")
    file.write(f"datatype: {datatype_ctype}\n")
    file.write(f"data_filename: {data_filename_ctype}\n")
print(f"meta cancer type is saved to '{output_path}'")
# print('\n') # 空行

# 製做 cancer type 資料 (靠自己手寫定義)
# 定義 cancer type 檔裡面需要的變數
type_of_cancer = "vghtpe" # 共用(必須跟 cancer_type.txt 裡面的第一個欄位相同)
type_of_cancer_full_name = "VGH-TPE Studies" # 定義的癌症名稱全名
show_color = "yellowgreen" # 要顯現出來的顏色
type_of_cancer_parent = "Other" # 癌症種類的母類別(有限制，要去原始程式碼改)

# 指定输出文件路徑(meta_cancer_type)
output_file_data_ctype = "cancer_type.txt"
# 拼接完整的路徑名稱
output_path = os.path.join(output_folder, output_file_data_ctype)
# 將文字寫入 .txt 文件檔中
with open(output_path, "w", encoding="utf-8") as file:
    file.write(f"{type_of_cancer}	")
    file.write(f"{type_of_cancer_full_name}	")
    file.write(f"{show_color}	")
    file.write(f"{type_of_cancer_parent}\n")
print(f"data cancer type is saved to '{output_path}'")
# 結束製作 meata data 的步驟
print('----------------------------------------- Finish Checking Meta Data -----------------------------------------')
print('\n') # 空行
print('\n') # 空行


# 製做 maf 資料 (資料來源: seqslab 的 label, maf)
print('---------------------- Start to Make Mutation Data(MAF) by Using Label and MAF Columns ----------------------')
# 建立 dataFrame
df = df_report[['label', 'maf']].copy()
df = df.rename(columns={'label': 'PATIENT_ID', 'maf': 'MAF'})
# 使用 apply 方法在 'MAF' 列上應用函數
df['MAF'] = df['MAF'].apply(replace_empty_list_with_none)
# 定義針對一些特殊情況做處理的函數，刪除冗言贅詞
df = preprocess_patient_id(df) # 使用前面定義的函數
# 查看 dataFrame 的形狀
rows, columns = df.shape
print(f'此次匯入的 maf 原始資料筆數: {rows}')
print(f'此次匯入的 maf 原始資料欄位數: {columns}')
print('\n') # 空行
# 確認基因變異資料
df.to_csv(check_data_folder + '/maf_from_seqslab.csv', index=False)
df.to_excel(check_data_folder + '/maf_from_seqslab.xlsx', index=False)

# 查看資料基本特性
analyze_dataframe(df, dataframe_name='df(load label and maf)')
print('\n') # 空行

# 建立一個用製做 maf (包含 mutation_extended 和 mutation_mskcc) 的 dataframe
# 新增 dataframe 的欄位
keys = ['Chromosome', 'End_Position', 'HGVSp_Short', 'Hugo_Symbol', 'Mutation_Status', 'NCBI_Build', 'RefSeq', 'Reference_Allele', 'Start_Position', 'Strand', 'Tumor_Sample_Barcode', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'Variant_Classification', 'Variant_Type', 'Protein_Change']
# 使用 kyes 來建立一個 dataFrame
empty_df = pd.DataFrame(columns=keys)
# 檢查、印出空白的 dataFrame
print(f'empty_df columns: {empty_df}')
# print('\n') # 空行

# 假设你有一個 dataFrame df 包含 'MAF' 列
# 找出 'MAF' 列中等於 NaN 的值
nan_indices = df['MAF'].isna()
# 統計 'MAF' 列中等於 NaN 的值的數量
nan_count = nan_indices.sum()
# 將 'MAF' 列中等於 NaN 的值替換成 None
df.loc[nan_indices, 'MAF'] = None
# 印出 NaN 值的數量
print(f'df 中 NaN 值的數量: {nan_count}')
# print('\n') # 空行

# 把對應的值填入 empty_df 來產出 maf
print('---------------- Mapping MAF Data, Please wait for few Seconds to few Minutes -----------------')
# 建立一個 List 來儲存沒有 mutation 資料的病例
no_mutations_patient = []
no_mutations_patient_origin = []
no_mutations_id = []
no_mutations_data = []
no_tumor_sample_barcode_patient_origin = set()
no_tumor_sample_barcode_id = set()
data_less_than_10_columns_patient_origin = set()
data_less_than_10_columns_id = set()
no_ncbi_build_patient_origin = set()
no_ncbi_build_id = set()
no_chromosome_patient_origin = set()
no_chromosome_id = set()

# 初始化 total_mutations
total_mutations = 0

# 有幾筆資料
for i in range(rows):
    # 檢查 df.MAF[i] 是否為 None
    if df.MAF[i] is not None and df.MAF[i] != 'NA': # 之前是寫 if df.MAF[i] is not None:
        # 使用 json.dumps() 將列表轉換為 JSON 格式的字符串
        df.MAF[i] = json.dumps(df.MAF[i])
        df.MAF[i] = df.MAF[i].replace("'", '"')
        # df.MAF[i] = df.MAF[i].replace(" ", "")

        # 使用 json.loads() 將 JSON 字串解析為 Python list
        data_list = json.loads(df.MAF[i])
        # 單筆資料裡面有多少次基因變異
        mutation_numbers = len(data_list)
        # print(f'第 {i+1} 筆資料有 {mutation_numbers} 變異點')
        # print(f'第 {i+1} 筆資料的 PATIENT_ID: {df.PATIENT_ID[i]}')
        # print(f'第 {i+1} 筆資料的 PATIENT_ID: {empty_df.Tumor_Sample_Barcode[i]}')
        for j in range(mutation_numbers):
            # 檢查 data_list[j] 這個字典裡有沒有 'Tumor_Sample_Barcode' 這個索引，如果沒有，用 PATIENT_ID 代替，並把原始 label 記起來
            if 'Tumor_Sample_Barcode' not in data_list[j]:
                data_list[j]['Tumor_Sample_Barcode'] = df.PATIENT_ID[i]
                # 將沒有 'Tumor_Sample_Barcode' 的資料添加到集合，以便後續檢查資料，確保單一不重複
                no_tumor_sample_barcode_patient_origin.add(df_report['label'].iloc[i])
                no_tumor_sample_barcode_id.add(i)
            # 檢查 data_list[j] 這個字典裡有沒有 'NCBI_Build' 這個索引，如果沒有，用 GRCH37 代替，並把原始 label 記起來
            if 'NCBI_Build' not in data_list[j]:
                data_list[j]['NCBI_Build'] = 'GRCH37'
                # 將沒有 'NCBI_Build' 的資料添加到集合，以便後續檢查資料，確保單一不重複
                no_ncbi_build_patient_origin.add(df_report['label'].iloc[i])
                no_ncbi_build_id.add(i)
            # 檢查 data_list[j] 這個字典裡有沒有 'Chromosome' 這個索引，如果沒有，把原始 label 記起來
            if 'Chromosome' not in data_list[j]:
                # data_list[j]['Chromosome'] = '7'
                # 將沒有 'Chromosome' 的資料添加到集合，以便後續檢查資料，確保單一不重複
                no_chromosome_patient_origin.add(df_report['label'].iloc[i])
                no_chromosome_id.add(i)
            if len(data_list[j]) < 10:
                data_less_than_10_columns_patient_origin.add(df_report['label'].iloc[i])
                data_less_than_10_columns_id.add(i)
            # 將字典裡面每個索引的值填入 DataFrame 中對應的位置
            for key, value in data_list[j].items():
                if key in empty_df.columns:
                    empty_df.at[j+total_mutations, key] = value
                else:
                    empty_df.at[j+total_mutations, key] = 'NA'
        # 現在累積多少筆基因變異的資料            
        total_mutations += mutation_numbers
    else:
        # 處理 df.MAF[i] 為 None 的情况
        # print(f"第 {i+1} 筆病例的 df.MAF[{i+1}] 是 None type，表示這筆病例沒有基因變異的資料")
        no_mutations_patient.append(df.PATIENT_ID[i]) # 儲存沒有基因變異資料的 PATIENT_ID(精簡)
        no_mutations_patient_origin.append(df_report['label'].iloc[i]) # 儲存沒有基因變異資料的 PATIENT_ID(原始)
        no_mutations_id.append(i) # 儲存沒有基因變異資料的 PATIENT_ID 的 index
        no_mutations_data.append(df_report['maf'].iloc[i]) # 儲存沒有基因變異資料的 PATIENT_ID 的 欄位資料，確認真的沒有資料

    # print(f'現在已經儲存 {total_mutations} 筆變異資料')
    # print('----------------------------------------------------------------------------')

print(f'全部總共有 {total_mutations} 筆變異資料')
print('\n') # 空行

# 將記錄沒有 tumor sample barcode, ncbi_build, chromosome 以及欄位數小於 10 的資料轉成列表的型態
no_tumor_sample_barcode_patient_origin = list(no_tumor_sample_barcode_patient_origin)
no_tumor_sample_barcode_id = list(no_tumor_sample_barcode_id)
no_ncbi_build_patient_origin = list(no_ncbi_build_patient_origin)
no_ncbi_build_id = list(no_ncbi_build_id)
no_chromosome_patient_origin = list(no_chromosome_patient_origin)
no_chromosome_id = list(no_chromosome_id)
data_less_than_10_columns_patient_origin = list(data_less_than_10_columns_patient_origin)
data_less_than_10_columns_id = list(data_less_than_10_columns_id)

# 建立 DataFrame 用來儲存沒有基因變異資料的病例
df_check_no_mutations = pd.DataFrame({'Index': no_mutations_id, 'Patient_id': no_mutations_patient_origin, 'Data': no_mutations_data})
# 匯出 CSV 文件
df_check_no_mutations.to_csv(check_data_folder + '/no_mutations_in_maf.csv', index=False)
# 匯出 Excel 文件
df_check_no_mutations.to_excel(check_data_folder + '/no_mutations_in_maf.xlsx', index=False)

# 使用函數匯出 check data 的檢查結果
export_dataframe_to_files(no_tumor_sample_barcode_id, no_tumor_sample_barcode_patient_origin, 'no_tumor_sample_barcode', check_data_folder)
export_dataframe_to_files(no_ncbi_build_id, no_ncbi_build_patient_origin, 'no_ncbi_build', check_data_folder)
export_dataframe_to_files(no_chromosome_id, no_chromosome_patient_origin, 'no_chromosome', check_data_folder)
export_dataframe_to_files(data_less_than_10_columns_id, data_less_than_10_columns_patient_origin, 'data_less_than_10_columns', check_data_folder)

print('-----------------------------------------------------------------------------------------------')
# 印出沒有基因變異資料的病人代碼(病例)
print(f'沒有基因變異資料的病例總共有: {len(no_mutations_patient)} 筆')
print(f'沒有基因變異資料的病例: {no_mutations_patient}')
print('-----------------------------------------------------------------------------------------------')
# 印出沒有 tumor sample barcode 的病人代碼(病例)
print(f'沒有 tumor sample barcode 的病例總共有: {len(no_tumor_sample_barcode_patient_origin)} 筆')
print(f'沒有 tumor sample barcode 的病例: {no_tumor_sample_barcode_patient_origin}')
print('-----------------------------------------------------------------------------------------------')
# 印出沒有 ncbi build id 的病人代碼(病例)
print(f'沒有 tumor sample barcode 的病例總共有: {len(no_ncbi_build_patient_origin)} 筆')
print(f'沒有 tumor sample barcode 的病例: {no_ncbi_build_patient_origin}')
print('-----------------------------------------------------------------------------------------------')
# 印出沒有 chromosome 的病人代碼(病例)
print(f'沒有 tumor sample barcode 的病例總共有: {len(no_chromosome_patient_origin)} 筆')
print(f'沒有 tumor sample barcode 的病例: {no_chromosome_patient_origin}')
print('-----------------------------------------------------------------------------------------------')
# 印出資料欄位少於 10 個的病人代碼(病例)
print(f'資料欄位少於 10 個的病例總共有: {len(data_less_than_10_columns_patient_origin)} 筆')
print(f'資料欄位少於 10 個的病例: {data_less_than_10_columns_patient_origin}')
print('-----------------------------------------------------------------------------------------------')
print('\n') # 空行

# 使用 `reorder_levels` 方法重新排列 DataFrame 列
order = ['Hugo_Symbol', 'NCBI_Build', 'Chromosome', 'Start_Position', 'End_Position', 'Strand', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode', 'Mutation_Status', 'HGVSp_Short', 'RefSeq', 'Protein_Change']
empty_df = empty_df[order]
# 使用 Pandas 的 .loc 方法來選擇並且替換符合條件的值
empty_df.loc[empty_df['NCBI_Build'] == 'GRCH37', 'NCBI_Build'] = 'GRCh37'
# 使用函數處理 DataFrame，處理 'HGVSp_Short' 欄位裡面儲存的資料
empty_df = clean_HGVSp_Short(empty_df)
# 確認最終 empty_df 的形狀
empty_df_rows, empty_df_columns = empty_df.shape
# print(f'mutation data 的資料筆數: {empty_df_rows}')
# print(f'mutation data 的資料欄位數: {empty_df_columns}')
# print('\n') # 空行
# 定義針對一些特殊情況做處理的函數，刪除冗言贅詞
empty_df = preprocess_tumor_sample_barcode(empty_df) # 使用前面定義的函數

# 把表格中填入 'na' 的地方改成 'NA'
empty_df = empty_df.replace('na', 'NA')
# 刪除沒有 'Tumor_Sample_Barcode' 的資料，並用 'Na' 填補遺失值
empty_df = empty_df[empty_df['Tumor_Sample_Barcode'].notna()]
empty_df = empty_df.fillna('NA')
# 將 dataframe 儲存成用空白分割的 .txt 文件
final_df = empty_df
print(f'shape of final_df mapping from maf: {final_df.shape}')
# 匯出 CSV 文件
final_df.to_csv(check_data_folder + '/final_df_mutation_data.csv', index=True)
# 匯出 Excel 文件
final_df.to_excel(check_data_folder + '/final_df_mutation_data.xlsx', index=True)
# 指定要储存的文件名和資料夾
output_file_1 = 'data_mutations_extended.txt'
output_file_2 = 'data_mutations_mskcc.txt'
# 拼接完整的路徑名稱
output_path_1 = os.path.join(output_folder, output_file_1)
output_path_2 = os.path.join(output_folder, output_file_2)
# 儲存資料
final_df.to_csv(output_path_1, sep='	', index=False, encoding='utf-8')
final_df.to_csv(output_path_2, sep='	', index=False, encoding='utf-8')
print(f'mutation data is saved to {output_path_1}')
print('---------------------- Finish Making Mutation Data(MAF) and Saving to the output path -----------------------')
print('\n') # 空行
print('\n') # 空行


# 製做 patient 和 sample 資料 (資料來源: seqslab 的 report, cnacer_type)
print('-------------------- Start to Make Patient Data by Using Report and Cancer Type Columns ---------------------')
# 先處理 report 資料
# 建立 dataFrame
df_rep = df_report[['report']].copy()

# 解析 JSON 資料並建立 dataFrame
data = []
for i in range(rows):
    data.append(df_rep.report[i])

# 建立 dataFrame
df_rep = pd.DataFrame([d['ner'] for d in data])
# 變更欄位的名稱
df_rep = df_rep.rename(columns={'REPORT_DIAGNOSIS': 'DIAGNOSIS', 'REPORT_PATIENT_ID': 'PATIENT_ID', 'REPORT_PATIENT_NAME': 'PATIENT_NAME'})
# 用 'NA' 填補遺失值
df_rep = df_rep.fillna('NA')
# 定義針對一些特殊情況做處理的函數，刪除冗言贅詞
df_rep = preprocess_patient_id(df_rep) # 使用前面定義的函數
# 查看 dataFrame 的形狀
rows, columns = df_rep.shape
print(f'report 原始的資料筆數: {rows}')
print(f'report 原始的資料欄位數: {columns}')
print('\n') # 空行

# 查看資料基本特性
analyze_dataframe(df_rep, dataframe_name='df_rep(load report)') # 使用前面定義的函數
print('\n') # 空行

print('-----------------------------------------------------------------------------------------------')
# 檢查 'REPORT_TEST_ASSAY' 是否都有值
missing_values = df_rep['REPORT_TEST_ASSAY'].isnull().sum()
if missing_values == 0:
    print("All 'REPORT_TEST_ASSAY' values are present.")
else:
    print(f"There are {missing_values} missing 'REPORT_TEST_ASSAY' values.")
print('-----------------------------------------------------------------------------------------------')
# 統計 'REPORT_TEST_ASSAY' 不同值的出現次數
assay_counts = df_rep['REPORT_TEST_ASSAY'].value_counts()
print("Counts of different 'REPORT_TEST_ASSAY' values:")
print(assay_counts)
print('-----------------------------------------------------------------------------------------------')

# 統計沒有基因變異資料的 'REPORT_TEST_ASSAY' 不同值的出現次數
# 篩選出沒有基因變異資料的資料
no_mutations_df = df_rep[df_rep['PATIENT_ID'].isin(no_mutations_patient)]
# 統計 'REPORT_TEST_ASSAY' 不同值的出現次數
assay_counts = no_mutations_df['REPORT_TEST_ASSAY'].value_counts()
# print("Counts of different 'REPORT_TEST_ASSAY' values for patients with no mutations:")
# print(assay_counts)
# print('\n') # 空行

# 印出沒有基因變異資料的病人代碼(病例)
no_name = []
for i in range(rows):
    if df_rep.PATIENT_NAME[i] == 'NA':
        # print(df.PATIENT_ID[i])
        no_name.append(df_rep.PATIENT_ID[i])
    else:
        pass

print(f'沒有病患名字的病例總共有: {len(no_name)} 筆')
print(f'沒有病患名字的病例: {no_name}')
print('-----------------------------------------------------------------------------------------------')

#先把原本的 df_rep 儲存起來備用
df_rep_origin = df_rep
for patient_id in no_mutations_patient:
    # 找到沒有基因變異資料的病例號索引
    indices_to_remove = df_rep[df_rep['PATIENT_ID'] == patient_id].index
    # 從 dataFrame 中刪除這些沒有基因變異資料的資料
    df_rep = df_rep.drop(indices_to_remove)
df_rep = df_rep.reset_index(drop=True)

# 檢查刪除沒有基因變異資料後 dataFrame 的形狀
rows, columns = df_rep.shape
# print(f'report 刪除沒有基因變異的資料筆數: {rows}')
# print(f'report 刪除沒有基因變異的資料資料欄位數: {columns}')
# print('\n') # 空行

# 檢查 'REPORT_ID' 是否有重複的資料
duplicates_report_id = df_rep['REPORT_ID'].duplicated()
# 檢查 'PATIENT_ID' 是否有重複的資料
duplicates_patient_id = df_rep['PATIENT_ID'].duplicated()
# 印出檢查的結果
# print(f'重複的 REPORT_ID: {df_rep[duplicates_report_id]}')
# print(f'重複的 PATIENT_ID: {df_rep[duplicates_patient_id]}')
# print('\n') # 空行

# 負責做基因檢測報告的公司或實驗室名稱列表，紀錄於 'REPORT_TEST_ASSAY' 欄位中，此列表用來做報告名稱的精簡化
assay_list = ["ACTG", "ACTOnco+", "FoundationOne", "BRCA", "Focus", "Archer", "Guardant360", "Myeloid", "Tumor Mutation Load"]
# 修改 FoundationOneLiquidDx 的值，變成 FoundationOne LiquidDx
df_rep.loc[df_rep['REPORT_TEST_ASSAY'] == "FoundationOneLiquidDx", 'REPORT_TEST_ASSAY'] = "FoundationOne Liquid CDx"
df_rep_origin.loc[df_rep_origin['REPORT_TEST_ASSAY'] == "FoundationOneLiquidDx", 'REPORT_TEST_ASSAY'] = "FoundationOne Liquid CDx"
# df_rep = simplify_report_names(df_rep, assay_list) # 根據前面定義的函數
# df_rep_origin = simplify_report_names(df_rep_origin, assay_list) # 根據前面定義的函數
# # 印出 REPORT_TEST_ASSAY 欄位裡面的值精簡化完成
# print('values in the column "REPORT_TEST_ASSAY" are renamed based on assay list')
# print('\n') # 空行

# 再處理 cancer_type 資料
print('------------------------ Start to Make Sample Data by Using Label and ICD_O Columns -------------------------')
# 建立 dataFrame
df_cancer = df_report[['label', 'icd_o']].copy()
df_cancer = df_cancer.rename(columns={'label': 'PATIENT_ID', 'icd_o': 'CANCER_TYPE'})
# 定義針對一些特殊情況做處理的函數，刪除冗言贅詞
df_cancer = preprocess_patient_id(df_cancer) # 使用前面定義的函數
# 將 'N/A' 和 'NaN' 替換成 None
df_cancer['CANCER_TYPE'] = df_cancer['CANCER_TYPE'].replace(['N/A', 'NaN'], None)
# 使用 apply 方法應用提取操作
df_cancer['CANCER_TYPE'] = df_cancer['CANCER_TYPE'].apply(extract_cancer_type)
# 查看資料基本特性
analyze_dataframe(df_cancer, dataframe_name='df_cancer(load label and icd_o to extract cancer type)') # 使用前面定義的函數

# 將 df 和 df_cancer 兩個 dataframe 用病人代碼(病例)作為索引合併起來
merged_df_origin = df_rep_origin.merge(df_cancer, on='PATIENT_ID', how='inner') # 保留 df_rep_orgin 和 df_cancer 都有的 'PATIENT_ID'
merged_df = df_rep.merge(df_cancer, on='PATIENT_ID', how='inner') # 保留 df_rep 和 df_cancer 都有的 'PATIENT_ID'

# 開始製做 patient 和 sample 資料，從 merge_df 這張大表中取值
# 製做 patient 資料，從 merge_df 中提取對應的欄位值
clinical_patient_df = merged_df[['PATIENT_ID', 'DIAGNOSIS']].copy()
clinical_patient_df_origin = merged_df_origin[['PATIENT_ID', 'DIAGNOSIS']].copy()
clinical_patient_df['STAGE'] = 'NA'
clinical_patient_df_origin['STAGE'] = 'NA'
# clinical_patient_df['PATIENT_DISPLAY_NAME'] = df['PATIENT_NAME']
clinical_patient_df['AGE'] = 'NA' 
clinical_patient_df_origin['AGE'] = 'NA' 
clinical_patient_df['SEX'] = 'NA' 
clinical_patient_df_origin['SEX'] = 'NA' 
clinical_patient_df['ETHNICITY'] = 'NA'
clinical_patient_df_origin['ETHNICITY'] = 'NA'

# 重新排列 dataFrame 的欄位順序
order = ['PATIENT_ID', 'DIAGNOSIS', 'STAGE', 'AGE', 'SEX', 'ETHNICITY']
clinical_patient_df = clinical_patient_df[order]
clinical_patient_df_origin = clinical_patient_df_origin[order]
# 用 'NA' 填補空白值
clinical_patient_df.fillna('NA', inplace=True)
clinical_patient_df_origin.fillna('NA', inplace=True)
# print(f'clinical_patient_df info: {clinical_patient_df.info()}')
print(f'clinical_patient_df_origin (exist no mutations data) info: {clinical_patient_df_origin.info()}')
print('\n') # 空行

# 製做 sample 資料，從 merge_df 中提取對應的欄位值
clinical_sample_df = merged_df[['PATIENT_ID', 'CANCER_TYPE', 'REPORT_TEST_ASSAY']].copy()
clinical_sample_df_origin = merged_df_origin[['PATIENT_ID', 'CANCER_TYPE', 'REPORT_TEST_ASSAY']].copy()
clinical_sample_df['SAMPLE_ID'] = merged_df['PATIENT_ID'] # df['REPORT_ID']
clinical_sample_df_origin['SAMPLE_ID'] = merged_df_origin['PATIENT_ID'] # df['REPORT_ID']
# clinical_sample_df['CANCER_TYPE'] = merged_df['CANCER_TYPE'] 
clinical_sample_df['CANCER_TYPE_DETAILED'] = 'NA' 
clinical_sample_df_origin['CANCER_TYPE_DETAILED'] = 'NA' 
# clinical_sample_df['SOMATIC_STATUS'] = 'NA' 
clinical_sample_df['REPORT_TEST_ASSAY'] = merged_df['REPORT_TEST_ASSAY']
clinical_sample_df_origin['REPORT_TEST_ASSAY'] = merged_df_origin['REPORT_TEST_ASSAY']

# 重新排列 DataFrame 列
order = ['PATIENT_ID', 'SAMPLE_ID', 'CANCER_TYPE', 'CANCER_TYPE_DETAILED', 'REPORT_TEST_ASSAY']
clinical_sample_df = clinical_sample_df[order]
clinical_sample_df_origin = clinical_sample_df_origin[order]
# 用 'NA' 填補空白值
clinical_sample_df.fillna('NA', inplace=True)
clinical_sample_df_origin.fillna('NA', inplace=True)
# print(f'clinical_sample_df info: {clinical_sample_df.info()}')
print(f'clinical_sample_df (exist no mutations data) info: {clinical_sample_df_origin.info()}')
print('\n') # 空行

# 指定要储存的文件名和資料夾
# output_folder = 'report_1'
output_file_3 = 'data_clinical_patient.txt'
output_file_4 = 'data_clinical_sample.txt'
# 拼接完整的路徑名稱
output_path_3 = os.path.join(output_folder, output_file_3)
output_path_4 = os.path.join(output_folder, output_file_4)

# 儲存資料
sep = '	'
with open(output_path_3, 'w', encoding="utf-8", newline='\r\n') as f:
    # 添加四行以#開頭的註解
    f.write(f'#Patient Identifier{sep}Diagnosis{sep}Overall Survival Status{sep}Age{sep}Sex{sep}ETHICITY\n')
    f.write(f'#Identifier to uniquely specify a patient.{sep}Overall survival in months since initial diagonosis.{sep}Overall patient survival status.{sep}Age at which a condition or disease was first diagnosed.{sep}Sex.{sep}Ethicity.\n')
    f.write(f'#STRING{sep}STRING{sep}NUMBER{sep}NUMBER{sep}STRING{sep}STRING\n')
    f.write(f'#1{sep}1{sep}1{sep}1{sep}1{sep}1\n')
    # 删除額外的換行符號
    clinical_patient_csv = clinical_patient_df_origin.to_csv(sep='	', index=False, encoding='utf-8')
    clinical_patient_csv = clinical_patient_csv.replace('\r\n', '\n')
    f.write(clinical_patient_csv)

# 儲存資料
with open(output_path_4, 'w', encoding="utf-8", newline='\r\n') as f:
    # 添加四行以#開頭的註解
    f.write(f'#Patient Identifier{sep}Sample Identifier{sep}Sample Cancer Type{sep}Cancer Type Detailed{sep}Report Test Assay\n')
    f.write(f'#Identifier to uniquely specify a patient.{sep}A unique sample identifier.{sep}Cancer Type.{sep}Cancer Type Detailed.{sep}Report Assay.\n')
    f.write(f'#STRING{sep}STRING{sep}STRING{sep}STRING{sep}STRING\n')
    f.write(f'#1{sep}1{sep}1{sep}1{sep}1\n')
    # 删除額外的換行符號
    clinical_sample_csv = clinical_sample_df_origin.to_csv(sep='	', index=False, encoding='utf-8')
    clinical_sample_csv = clinical_sample_csv.replace('\r\n', '\n')
    f.write(clinical_sample_csv)

print(f'patient data is saved to {output_path_3}')
print(f'sample data is saved to {output_path_4}')
print('----------------- Finish Making Patient and Sample Data(MAF) and Saving to the output path ------------------')
print('\n') # 空行
print('\n') # 空行

# 製作 cases sequenced (放在 case lists 資料夾)
print('----------------- Start to Make Cases Sequenced by Using Sample Data Generated Previously -------------------')
# 從 clinical_sample_df 這個 dataframe 中提取 'SAMPLE_ID' 的值
case_list_sequenced_df = clinical_sample_df[['SAMPLE_ID']].copy()
# 取得 'SAMPLE_ID' 列的值
sample_ids = case_list_sequenced_df['SAMPLE_ID']
# 定義一個變數來儲存 cases sequenced
text_data_sequenced = '	'.join(str(sample_id) for sample_id in sample_ids)

# 製做 cna (資料來源: seqslab 的 label, cna, cnv)
print('------ Start to Make CNA Data by Using CNA, CNV Columns and Mutation Data(MAF) Generated Previously ---------')
# 需要用到前面製做完成的 final_df dataframe
cna_df = final_df[['Hugo_Symbol', 'Tumor_Sample_Barcode']].copy()
cna_df = cna_df[cna_df['Tumor_Sample_Barcode'] != 'Na']
# cna_df = preprocess_tumor_sample_barcode(cna_df)

# 建立一個全零矩陣，索引為'Hugo_Symbol'，欄位為'Tumor_Sample_Barcode'
matrix = pd.DataFrame(np.zeros((len(cna_df['Hugo_Symbol'].unique()), len(cna_df['Tumor_Sample_Barcode'].unique())), dtype=int),
                     index=cna_df['Hugo_Symbol'].unique(),
                     columns=cna_df['Tumor_Sample_Barcode'].unique())

# 重命名索引為欄位名稱
matrix.index.name = 'Hugo_Symbol'

# 檢查 'REPORT_ID' 是否有重複的資料
duplicates_report_id = df_rep['REPORT_ID'].duplicated()
# 檢查 'PATIENT_ID' 是否有重複的資料
duplicates_patient_id = df_rep['PATIENT_ID'].duplicated()
# 印出檢查的結果
# print(f'重複的 REPORT_ID: {df_rep[duplicates_report_id]}')
# print(f'重複的 PATIENT_ID: {df_rep[duplicates_patient_id]}')
# print('\n') # 空行

# 建立 dataFrame
# 建立 dataFrame，只包含 'label' 列
df = df_report[['label']].copy()
df = df.rename(columns={'label': 'PATIENT_ID'})
if 'cna' not in df_report: # 如果 'cna' 列不存在，將其加到 df 並填 None
    df['CNA'] = None
else: # 如果 'cna' 列存在，將其值赋给 df
    df['CNA'] = df_report['cna']
if 'cnv' not in df_report: # 如果 'cnv' 列不存在，將其加到 df 並填 None
    df['CNV'] = None
else: # 如果 'cnv' 列存在，將其值赋给 df
    df['CNV'] = df_report['cnv']
df = df.rename(columns={'label': 'PATIENT_ID', 'cna': 'CNA', 'cnv': 'CNV'})
# 定義針對一些特殊情況做處理的函數，刪除冗言贅詞
df = preprocess_patient_id(df)
# 用迴圈移除掉沒有基因變異資料的病例號索引
for patient_id in no_mutations_patient:
    # 找到沒有基因變異資料的病例號索引
    indices_to_remove = df[df['PATIENT_ID'] == patient_id].index
    # 從 dataFrame 中刪除這些沒有基因變異資料的資料
    df = df.drop(indices_to_remove)
df = df.reset_index(drop=True)

# 建立一個新的欄位 'Final_Value'，根據條件選取值
df['Final_Value'] = np.where(df['CNA'].notna(), df['CNA'], np.where(df['CNV'].notna(), df['CNV'], None))
# 刪除 'Final_Value' 欄位中沒有值的資料
clean_df = df.dropna(subset=['Final_Value'])
clean_df = clean_df.reset_index(drop=True)

# 最後 CNA 的資料數
rows, columns = clean_df.shape
print(f'cna data 的資料筆數: {rows}')
print(f'cna data 的資料欄位數: {columns}')

print('---------------- Mapping CNA Data, Please wait for few Seconds to few Minutes -----------------')
for i in range(rows):
    # 使用 json.dumps() 將列表轉換為 JSON 格式的字符串
    clean_df.Final_Value[i] = json.dumps(clean_df.Final_Value[i])
    clean_df.Final_Value[i] = clean_df.Final_Value[i].replace("'", '"')
    # clean_df.Final_Value[i] = clean_df.Final_Value[i].replace(" ", "")
    data_list = json.loads(clean_df.Final_Value[i])
    gene_number = len(data_list)
    # print(f'第 {i+1} 筆資料 {clean_df.PATIENT_ID[i]} 出現 {gene_number} 個基因')
    for j in range(gene_number):
        data = data_list[j]
        # 處理 CNA 的 COPY_NUMBER
        copy_number = data['ner']['COPY_NUMBER']
        matches = [float(match.group()) for match in re.finditer(r'-?\d+\.\d+|-?\d+', copy_number)]
        if matches:
            copy_number = round(matches[0])
            if copy_number == 0:
                copy_number = -2
            elif copy_number == 1:
                copy_number = -1
            elif copy_number == 2:
                copy_number = 0
            elif 2 < copy_number <= 3.5:
                copy_number = 1
            else:
                copy_number = 2
        
        gene = data['ner']['GENE']
        # variant_form = data['ner']['VARIANT_FORM'] # 不是每一筆資料都有，所以不要亂存

        # 儲存 CNA 資料到 matrix 裡面
        matrix[clean_df.PATIENT_ID[i]][gene] = copy_number

        # print('----------------------------------')
        # print(f'第 {i+1} 筆資料的 第 {j+1} 個基因:')
        # print("COPY_NUMBER:", copy_number)
        # print("GENE:", gene)
        # print("VARIANT_FORM:", variant_form) # 不是每一筆資料都有，所以不要亂印
    # print('-----------------------------------------------------')

# print('\n') # 空行

# 儲存 cna 矩陣資料的.txt文件
# 指定要储存的文件名和資料夾
output_file_cna = 'data_cna.txt'
# 拼接完整的路徑名稱
output_file_path = os.path.join(output_folder, output_file_cna)
matrix.to_csv(output_file_path, sep='	', encoding='utf-8')
print(f'cna data is saved to {output_file_path}')
print('--------------------------- Finish Making CNA Data and Saving to the output path ----------------------------')
print('\n') # 空行
print('\n') # 空行

# 製做 cases cna (放在 case lists 資料夾)
# 取得 'SAMPLE_ID' 列的值
sample_ids = case_list_sequenced_df['SAMPLE_ID']
# 定義一個變數來儲存 cases cna
text_data_cna = '	'.join(str(sample_id) for sample_id in sample_ids)

# 製做 sv (資料來源: seqslab 的 label, fusion)
print('--------------------------- Start to Make SV Data by Using Label, Fusion Columns ----------------------------')
# 建立 dataFrame
sv_df = df_report[['label', 'fusion']].copy()
sv_df = sv_df.rename(columns={'label': 'PATIENT_ID', 'fusion': 'FUSION'})
# 定義針對一些特殊情況做處理的函數，刪除冗言贅詞
sv_df = preprocess_patient_id(sv_df)

# 檢查 dataframe 的形狀和缺漏值情況
print(f'sv_df shape (origin): {sv_df.shape}')
# print(f'sv_df info: {sv_df.info()}')
# print('\n') # 空行

# 保留 'FUSION' 欄位不為空的資料
sv_df = sv_df[sv_df['FUSION'].notna()]
# 重新设置索引
sv_df = sv_df.reset_index(drop=True)
# 再次檢查 dataframe 的形狀和缺漏值情況
print(f'sv_df shape (after remove "na"): {sv_df.shape}')

# 把沒有基因變異資料的的病例號刪除，重新獲取 dataframe 的資料筆數與欄位數目
for patient_id in no_mutations_patient:
    # 找到沒有基因變異資料的病例號索引
    indices_to_remove = sv_df[sv_df['PATIENT_ID'] == patient_id].index
    # 從 dataFrame 中刪除這些沒有基因變異資料的資料
    sv_df = sv_df.drop(indices_to_remove)
sv_df = sv_df.reset_index(drop=True)
print(f'sv_df shape (after drop no fusion data): {sv_df.shape}')
print('\n') # 空行
rows, columns = sv_df.shape

# 新增 dataframe 的欄位
keys = ['Sample_Id', 'SV_Status', 'Site1_Hugo_Symbol', 'Site1_Region', 'Site1_Region_Number', 'Site2_Hugo_Symbol', 'Site2_Region', 'Site2_Region_Number']
# 使用 kyes 來建立一個 dataFrame
save_sv_df = pd.DataFrame(columns=keys)
# 印出空白的 dataFrame
print(f'save_sv_df columns: {save_sv_df}')
print('\n') # 空行

# 定義正則表達式來儲存 SV 資料
print('----------------- Mapping SV Data, Please wait for few Seconds to few Minutes -----------------')
pattern = r"([a-z]+) (\d+) of the ([A-Z0-9]+) gene and ([a-z]+) (\d+) of the ([A-Z0-9]+) gene"
count = 1

for i in range(rows):
    # 使用 json.dumps() 將列表轉換為 JSON 格式的字符串
    sv_df.FUSION[i] = json.dumps(sv_df.FUSION[i])
    sv_df.FUSION[i] = sv_df.FUSION[i].replace("'", '"')
    # sv_df.FUSION[i] = sv_df.FUSION[i].replace(" ", "")
    data_list = json.loads(sv_df.FUSION[i])
    # print(f'第 {i+1} 筆資料裡面有 {len(data_list)} 組成對的基因資料')
    for j in range(len(data_list)):
        data = data_list[j]
        # 處理 NER 的 FUSION_BREAKPOINT 
        text_data = data['ner']['FUSION_BREAKPOINT']
    
        # 建立一個列表來儲存結果
        matches = []
        
        match = re.search(pattern, text_data)
        if match:
            sv = {
                'site_1_region': match.group(1),
                'site_1_region_number': int(match.group(2)),
                'site_1_hugo_symbol': match.group(3),
                'site_2_region': match.group(4),
                'site_2_region_number': int(match.group(5)),
                'site_2_hugo_symbol': match.group(6)
            }
            matches.append(sv)

        # 印出結果
        for idx, match in enumerate(matches, start=1):
            # print('----------------------------------')
            # 儲存到 save_sv_df
            save_sv_df.at[count, 'Sample_Id'] = sv_df.PATIENT_ID[i]
            save_sv_df.at[count, 'SV_Status'] = 'Somatic'
            save_sv_df.at[count, 'Site1_Region'] = match['site_1_region']
            save_sv_df.at[count, 'Site1_Region_Number'] = match["site_1_region_number"]
            save_sv_df.at[count, 'Site1_Hugo_Symbol'] = match["site_1_hugo_symbol"]
            save_sv_df.at[count, 'Site2_Region'] = match["site_2_region"]
            save_sv_df.at[count, 'Site2_Region_Number'] = match["site_2_region_number"]
            save_sv_df.at[count, 'Site2_Hugo_Symbol'] = match["site_2_hugo_symbol"]

            # 印出數值
            # print(f"Match {j+1}:")
            # print(f'Site1_Region: {match["site_1_region"]}')
            # print(f'Site1_Region_Number: {match["site_1_region_number"]}')
            # print(f'Site1_Hugo_Symbol: {match["site_1_hugo_symbol"]}')
            # print(f'Site2_Region: {match["site_2_region"]}')
            # print(f'Site2_Region_Number: {match["site_2_region_number"]}')
            # print(f'Site2_Hugo_Symbol: {match["site_2_hugo_symbol"]}')
        count += 1
    # print('-----------------------------------------------------')

# print('\n') # 空行

# 填補缺漏值
save_sv_df.fillna('NA', inplace=True)
# 指定要储存的文件名和資料夾
output_file_sv = 'data_sv.txt'
# 拼接完整的路徑名稱
output_file_path = os.path.join(output_folder, output_file_sv)
save_sv_df.to_csv(output_file_path, sep='	', encoding='utf-8', index=False)
print(f'sv data is saved to {output_file_path}')
print('---------------------------- Finish Making SV Data and Saving to the output path ----------------------------')
print('\n') # 空行
print('\n') # 空行

# 製做 cases sv (放在 case lists 資料夾)
# 取得 'SAMPLE_ID' 列的值
sample_ids = save_sv_df['Sample_Id']
# 定義一個變數來儲存 cases sv
text_data_sv = '	'.join(str(sample_id) for sample_id in sample_ids)

# 製做各種不同資料的 case list (包含 cases sequenced, cases cna 和 cases sv)
# 定義用來製做 cases sequenced 檔案的變數
stable_id_se = cancer_study_identifier + "_sequenced"
case_list_name_se = "Samples with mutation data"
case_list_description_se = "Samples with mutation data"
case_list_category_se = "all_cases_with_mutation_data"
case_list_ids_se = text_data_sequenced

# 定義用來製做 cases cna 檔案的變數
stable_id_cna = cancer_study_identifier + "_cna"
case_list_name_cna = "Samples with cna data"
case_list_description_cna = "Samples with cna data"
case_list_category_cna = "all_cases_with_cna_data"
case_list_ids_cna = text_data_cna

# 定義用來製做 cases sv 檔案的變數
stable_id_sv = cancer_study_identifier + "_sv"
case_list_name_sv = "Samples with sv data"
case_list_description_sv = "Samples with sv data"
case_list_category_sv = "all_cases_with_sv_data"
case_list_ids_sv = text_data_sv

# 指定输出文件路徑(cases_sequenced)
output_file_cases_sequenced = "cases_sequenced.txt"
# 拼接完整的路徑名稱
output_path = os.path.join(output_folder_case, output_file_cases_sequenced)
# 寫到 case sequenced
with open(output_path, "w", encoding="utf-8") as file:
    file.write(f"cancer_study_identifier: {cancer_study_identifier}\n")
    file.write(f"stable_id: {stable_id_se}\n")
    file.write(f"case_list_name: {case_list_name_se}\n")
    file.write(f"case_list_description: {case_list_description_se}\n")
    file.write(f"case_list_category: {case_list_category_se}\n")
    file.write(f"case_list_ids:	{case_list_ids_se}\n")
print(f"cases sequenced is saved to '{output_path}'")
# print('\n') # 空行

# 指定输出文件路徑(cases_cna)
output_file_cases_cna = "cases_cna.txt"
# 拼接完整的路徑名稱
output_path = os.path.join(output_folder_case, output_file_cases_cna)
# 寫到 case cna
with open(output_path, "w", encoding="utf-8") as file:
    file.write(f"cancer_study_identifier: {cancer_study_identifier}\n")
    file.write(f"stable_id: {stable_id_cna}\n")
    file.write(f"case_list_name: {case_list_name_cna}\n")
    file.write(f"case_list_description: {case_list_description_cna}\n")
    file.write(f"case_list_category: {case_list_category_cna}\n")
    file.write(f"case_list_ids:	{case_list_ids_cna}\n")
print(f"cases cna is saved to '{output_path}'")
# print('\n') # 空行

# 指定输出文件路徑(cases_sv)
output_file_cases_sv = "cases_sv.txt"
# 拼接完整的路徑名稱
output_path = os.path.join(output_folder_case, output_file_cases_sv)
# 寫到 case sv
with open(output_path, "w", encoding="utf-8") as file:
    file.write(f"cancer_study_identifier: {cancer_study_identifier}\n")
    file.write(f"stable_id: {stable_id_sv}\n")
    file.write(f"case_list_name: {case_list_name_sv}\n")
    file.write(f"case_list_description: {case_list_description_sv}\n")
    file.write(f"case_list_category: {case_list_category_sv}\n")
    file.write(f"case_list_ids:	{case_list_ids_sv}\n")
print(f"cases sv is saved to '{output_path}'\n")
print('----------------- Finish Making Cases Sequenced, CNA, and SV, and Saving to the output path -----------------')
print('\n') # 空行
print('\n') # 空行

print('Data for uploading to cBio Portal is successfully generated.')
print('Please check the "data_upload" folder in the "cbioportal-docker-compose/study" directory.')
print('-------------------- Finish cBio Portal Preprocess. Data is Generated and Uploaded to cBio Portal Process --------------------')