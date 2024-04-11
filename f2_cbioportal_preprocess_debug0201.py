# ======================================== 匯入套件 ========================================
# 匯入套件
import json
import pandas as pd
import os
import csv
import codecs
import requests
import sys
import numpy as np
import re
import argparse
from io import StringIO
from io import BytesIO
from datetime import datetime
import pytz

# ======================================== 定義函數 ========================================
# 解析 Json 的檔案名稱來找出 Path Number 和 MP Number
def parse_file_name(file_name):
    # 定義檔案名稱的正則表達式
    pattern = re.compile(r'^(.*?)_\((.*?)\)\.json$')
    
    # 使用不同的 pattern 來解析
    match = pattern.match(file_name)
    
    # 定義檔案名稱的正則表達式
    # pattern_1 = re.compile(r'(.*?)\.json')
    # pattern_2 = re.compile(r'(.*?)_(.*?)\.json')
    # pattern_3 = re.compile(r'(.*?)_(.*?)_(.*?)\.json')

    # 使用不同的 pattern 來解析
    # match_1 = pattern_1.match(file_name)
    # match_2 = pattern_2.match(file_name)
    # match_3 = pattern_3.match(file_name)

    # 選擇有效的匹配
    # match = match_3 or match_2 or match_1

    if match:
        # 根據匹配的 pattern，取得對應的捕獲組
        path_number = match.group(1)
        mp_number = match.group(2)
        # path_number = match.group(1)
        # mp_number = match.group(2) if match_2 or match_3 else None
        return path_number.replace(" ", ""), mp_number.replace("(", "").replace(")", "").replace(" ", "")
    else:
        return None, None

# 讀取 Json 檔
def read_json_files(directory_path):
    # 用來存放所有 JSON 檔案資料的列表
    all_data = []
    filenames = []
    path_numbers = []
    mp_numbers = []
    problematic_files = []  # 用來記錄有問題的(無法讀取)的 JSON 檔

    # 列出目錄下的所有檔案
    for filename in os.listdir(directory_path):
        # 檢查是否為 .json 檔案
        if filename.endswith(".json"):
            file_path = os.path.join(directory_path, filename)
            
            # 讀取 JSON 檔案內容
            with open(file_path) as f:
                try:
                    data = json.load(f)
                    # 如果 JSON 文件內容是空白的，把有問題的文件名記錄下來
                    if not data:
                        problematic_files.append(filename)
                        continue

                    # 印出檔案名稱和類型
                    #print(f"File: {filename}, Type: {type(data)}")

                    # 解析 Json 的檔案名稱來找出 Path Number 和 MP Number
                    path_number, mp_number = parse_file_name(filename)
                    filenames.append(filename)
                    path_numbers.append(path_number)
                    mp_numbers.append(mp_number)

                    # 將資料添加到列表中
                    all_data.append(data)

                except json.JSONDecodeError as e:
                    # 把無法解析的 JSON 文件名記錄下來
                    problematic_files.append(filename)
                    #print(f"Error decoding JSON file {filename}: {e}")

    # 返回所有 JSON 檔案資料的列表
    return all_data, filenames, path_numbers, mp_numbers, problematic_files

# 拼湊出 Patient Id
def find_patient_id(path_number, mp_number):
    if path_number is not None and mp_number is not None:
        patient_id = path_number + '_' + mp_number
        return patient_id.replace("(", "").replace(")", "")
    elif path_number is not None and mp_number is None:
        patient_id = path_number
        return patient_id.replace("(", "").replace(")", "")
    elif path_number is None and mp_number is not None:
        patient_id = path_number
        return patient_id.replace("(", "").replace(")", "")
    else:
        #print("Both path_number and mp_number are none.")
        return patient_id.replace("(", "").replace(")", "")

# 移除文字中的底線
def replace_underscore_with_space(input_text):
    if isinstance(input_text, int):
        return input_text
    elif isinstance(input_text, float):
        return input_text
    else:
        input_text = input_text.replace('_', ' ').replace('\r\n', '').replace('\n', '')
        # 使用原始字符串方式
        formatted_text = r"{}".format(input_text)
        return formatted_text

# 提取 position 欄位中的 chromosome、start_position 和 end_position
def extract_chromosome_start_and_end_position(position_str):
    # 先檢查 position_str 不等於 None
    if position_str is None:
        #print("position_str is None.")
        return None, None, None

    # 提取 position 欄位中的 chromosome、start_position 和 end_position 的正則表達式
    pattern = re.compile(r'(?P<chromosome>chr[XY\d]+):(?P<start_position>\d+)(?:_(?P<end_position>\d+))?')
    match = pattern.match(position_str)
    if match:
        chromosome = match.group('chromosome')
        start_position = match.group('start_position')
        end_position = match.group('end_position')

        # 如果 end_position 不存在，則設為和 start_position 相同的值(需要跟醫生確認)
        end_position = end_position if end_position is not None else start_position
        return chromosome, start_position, end_position
    else:
        return None, None, None

# 提取 cds_effect 欄位中的 reference_allele、tumor_seq_allele1 和 tumor_seq_allele2
def extract_ref_allele_tumor_seq_allele1_and_allele2(cds_effect_str):
    # 先檢查確認 cds_effect_str 是否為 None
    if cds_effect_str is None:
        #print("cds_effect_str is None.")
        return None, None, None

    match = re.search(r'([AUCGTaucgt])>([AUCGTaucgt])', cds_effect_str)
    if match:
        reference_allele, tumor_seq_allele2 = match.groups()
        tumor_seq_allele1 = reference_allele
        return reference_allele, tumor_seq_allele1, tumor_seq_allele2

    return None, None, None

# 找 gender
def find_gender(data):
    # 在 final_report 中查找
    if "ResultsReport" in data and "ResultsPayload" in data["ResultsReport"]:
        final_report_gender = data["ResultsReport"]["ResultsPayload"].get("FinalReport", {}).get("PMI", {}).get("Gender")
        if final_report_gender:
            if final_report_gender.lower() == 'f':
                final_report_gender = 'Female'
                return final_report_gender
            elif final_report_gender.lower() == 'm':
                final_report_gender = 'Male'
                return final_report_gender
            else:
                return final_report_gender
    
    # # 在 variant_report 中查找
    # variant_report_gender = data["ResultsReport"]["ResultsPayload"].get("variant_report", {}).get("gender")
    # if variant_report_gender:
    #     if variant_report_gender.lower() == 'f':
    #         variant_report_gender = 'Female'
    #         return variant_report_gender
    #     elif variant_report_gender.lower() == 'm':
    #         variant_report_gender = 'Male'
    #         return variant_report_gender
    #     else:
    #         return variant_report_gender
    
    # 如果都找不到，返回 None 或者其他適當的預設值
    return None

# 找 report_test_assay
def find_test_assay(data):
    # 在 final_report 中查找
    if "ResultsReport" in data and "ResultsPayload" in data["ResultsReport"]:
        final_report_test_assay = data["ResultsReport"]["ResultsPayload"].get("FinalReport", {}).get("Sample", {}).get("TestType")
        if final_report_test_assay:
            return replace_underscore_with_space(final_report_test_assay)
    
    # 在 variant_report 中查找
    variant_report_test_assay = data["ResultsReport"]["ResultsPayload"].get("variant_report", {}).get("test_type")
    if variant_report_test_assay:
        return replace_underscore_with_space(variant_report_test_assay)
    
    # 如果都找不到，返回 None 或者其他適當的預設值
    return None

# 找 diagnosis
def find_diagnosis(data):
    # # 在 final_report 中查找
    # if "ResultsReport" in data and "ResultsPayload" in data["ResultsReport"]:
    #     variant_report_pathology_diagnosis = data["ResultsReport"]["ResultsPayload"].get("variant_report", {}).get("pathology_diagnosis")
    #     if variant_report_pathology_diagnosis:
    #         return replace_underscore_with_space(variant_report_pathology_diagnosis)
            
    # 在 variant_report 中查找
    if "ResultsReport" in data and "ResultsPayload" in data["ResultsReport"]:
        final_report_diagnosis = data["ResultsReport"]["ResultsPayload"].get("FinalReport", {}).get("PMI", {}).get("SubmittedDiagnosis") # 主要
        if final_report_diagnosis:
            return replace_underscore_with_space(final_report_diagnosis)
        
    # 如果都找不到，返回 None 或者其他適當的預設值
    return None

# 找 cancer_type
def find_cancer_type(data):
    # # 在 final_report 中查找
    # if "ResultsReport" in data and "ResultsPayload" in data["ResultsReport"]:
    #     final_report_cancer_type = data["ResultsReport"]["ResultsPayload"].get("FinalReport", {}).get("PMI", {}).get("Disease")
    #     if final_report_cancer_type:
    #         return replace_underscore_with_space(final_report_cancer_type)
    
    # # 在 variant_report 中查找
    # variant_report_cancer_type = data["ResultsReport"]["ResultsPayload"].get("variant_report", {}).get("disease")
    # if variant_report_cancer_type:
    #     return replace_underscore_with_space(variant_report_cancer_type)
    
    # 在 variant_report 中查找
    if "ResultsReport" in data and "ResultsPayload" in data["ResultsReport"]:
        final_report_cancer_type = data["ResultsReport"]["ResultsPayload"].get("FinalReport", {}).get("PMI", {}).get("TumorType")
        if final_report_cancer_type:
            return replace_underscore_with_space(final_report_cancer_type)
    
    # 如果都找不到，返回 None 或者其他適當的預設值
    return None

# # 找 cancer_type_detailed (暫時應該是用不到了)
# def find_cancer_type_detailed(data):
#     # 在 final_report 中查找
#     if "ResultsReport" in data and "ResultsPayload" in data["ResultsReport"]:
#         final_report_cancer_type_detailed = data["ResultsReport"]["ResultsPayload"].get("FinalReport", {}).get("PMI", {}).get("DiseaseOntology")
#         if final_report_cancer_type_detailed:
#             return replace_underscore_with_space(final_report_cancer_type_detailed)

#     # 在 variant_report 中查找
#     variant_report_cancer_type_detailed = data["ResultsReport"]["ResultsPayload"].get("variant_report", {}).get("disease_ontology")
#     if variant_report_cancer_type_detailed:
#         return replace_underscore_with_space(variant_report_cancer_type_detailed)
    
#     # 如果都找不到，返回 None 或者其他適當的預設值
#     return None

# 自動跟常見癌症列表比較，如果傳入的 cancer type 符合其中一個癌症類別，則回傳該類別
def find_cancer_category(cancer_type, cancer_type_dict):
    # Standardize the input cancer_type by converting to lowercase and replacing underscores and hyphens
    # standardized_cancer_type = cancer_type.lower().replace('_', ' ').replace('-', ' ')
    standardized_cancer_type = cancer_type.replace('_', ' ').replace('-', ' ')

    # 和常見癌症類別做比較(不需要此步驟，直接從 ngs.csv 拿 cancer type 這個欄位即可)
    # for category, cancer_list in cancer_type_dict.items():
    #     # Standardize each item in the cancer_list for comparison
    #     standardized_cancer_list = [item.lower().replace('_', ' ').replace('-', ' ') for item in cancer_list]

    #     if standardized_cancer_type in standardized_cancer_list:
    #         return category
    #     else:
    #         return cancer_type
    # return cancer_type

    return standardized_cancer_type

# 找 percent_tumor_nuclei
def find_percent_tumor_nuclei(data):
    # # 在 final_report 中查找
    # if "ResultsReport" in data and "ResultsPayload" in data["ResultsReport"]:
    #     final_report_percent_tumor_nuclei = data["ResultsReport"]["ResultsPayload"].get("FinalReport", {}).get("PMI", {}).get("PercentTumorNuclei")
    #     if final_report_percent_tumor_nuclei is not None and final_report_percent_tumor_nuclei != "-":
    #         return replace_underscore_with_space(final_report_percent_tumor_nuclei)
    
    # # 在 variant_report 中查找
    # variant_report_percent_tumor_nuclei = data["ResultsReport"]["ResultsPayload"].get("variant_report", {}).get("percent_tumor_nuclei")
    # if variant_report_percent_tumor_nuclei is not None and final_report_percent_tumor_nuclei != "-":
    #     return replace_underscore_with_space(variant_report_percent_tumor_nuclei)
    
    # 在 final_report 中查找
    if "ResultsReport" in data and "ResultsPayload" in data["ResultsReport"]:
        final_report_percent_tumor_nuclei = data["ResultsReport"]["ResultsPayload"].get("FinalReport", {}).get("Sample", {}).get("TumorPurity")
        if final_report_percent_tumor_nuclei is not None and final_report_percent_tumor_nuclei != "-":
            return replace_underscore_with_space(final_report_percent_tumor_nuclei)
    
    # # 在 variant_report 中查找
    # variant_report_percent_tumor_nuclei = data["ResultsReport"]["ResultsPayload"].get("variant_report", {}).get("tumor_purity")
    # if variant_report_percent_tumor_nuclei is not None and final_report_percent_tumor_nuclei != "-":
    #     return replace_underscore_with_space(variant_report_percent_tumor_nuclei)
    
    # 如果都找不到，返回 None 或者其他適當的預設值
    return None

# 儲存 sample dataframe
def save_sample_data(data_sample_df, output_folder, output_file_data_clinical_sample):
    # # 將缺失值填充為 'NA'，並在原地進行修改
    # data_sample_df.fillna('NA', inplace=True)
    # 空白值、遺失值填上 'NA'
    data_sample_df = data_sample_df.fillna('NA')
    data_sample_df = data_sample_df.replace('N/A', 'NA')
    # 拼接完整的路徑名稱
    output_file_path = os.path.join(output_folder, output_file_data_clinical_sample)

    # 儲存資料
    sep = '	'
    with open(output_file_path, 'w', encoding="utf-8", newline='\r\n') as f:
        # 添加四行以#開頭的註解
        f.write(f'#Patient Identifier{sep}Sample Identifier{sep}Path Number{sep}MP Number{sep}Cell Block Number{sep}Tumor Type{sep}Tumor Purity{sep}Report Assay{sep}Report Date\n')
        f.write(f'#Identifier to uniquely specify a patient.{sep}A unique sample identifier.{sep}A patient identifier.{sep}A report identifier.{sep}A cell block number.{sep}Tumor type.{sep}Proportion of tumor cells within a tissue sample.{sep}Report Assay.{sep}Report Date.\n')
        f.write(f'#STRING{sep}STRING{sep}STRING{sep}STRING{sep}STRING{sep}STRING{sep}NUMBER{sep}STRING{sep}STRING\n')
        f.write(f'#1{sep}1{sep}1{sep}1{sep}1{sep}1{sep}1{sep}1{sep}1\n')
        # 儲存 .txt 文件
        data_sample_df = data_sample_df.to_csv(sep='	', index=False, encoding='utf-8')
        # 删除額外的換行符號
        data_sample_df = data_sample_df.replace('\r\n', '\n')
        f.write(data_sample_df)
    #print(f'data sample is saved to {output_file_path}')

# 儲存 sample dataframe
def save_patient_data(data_patient_df, output_folder, output_file_data_clinical_patient):
    # # 空白值、遺失值填上 'NA'
    # data_patient_df = data_patient_df.fillna('NA', inplace=True)
    # 空白值、遺失值填上 'NA'
    data_patient_df = data_patient_df.fillna('NA')
    data_patient_df = data_patient_df.replace('N/A', 'NA')
    # 拼接完整的路徑名稱
    output_file_path = os.path.join(output_folder, output_file_data_clinical_patient)

    # 儲存資料
    sep = '	'
    with open(output_file_path, 'w', encoding="utf-8", newline='\r\n') as f:
        # 添加四行以#開頭的註解
        f.write(f'#Patient Identifier{sep}Diagnosis{sep}Overall Survival Status{sep}Age{sep}Gender{sep}ETHICITY\n')
        f.write(f'#Identifier to uniquely specify a patient.{sep}Type of disease determined through pathology examination.{sep}Overall patient survival status.{sep}Age at which a condition or disease was first diagnosed.{sep}The gender of the patient.{sep}Ethicity.\n')
        f.write(f'#STRING{sep}STRING{sep}NUMBER{sep}NUMBER{sep}STRING{sep}STRING\n')
        f.write(f'#1{sep}1{sep}1{sep}1{sep}1{sep}1\n')
        # 儲存 .txt 文件
        data_patient_df = data_patient_df.to_csv(sep='	', index=False, encoding='utf-8')
        # 删除額外的換行符號
        data_patient_df = data_patient_df.replace('\r\n', '\n')
        f.write(data_patient_df)
    #print(f'data patient is saved to {output_file_path}')

# 正規化 percent_reads 的格式，應該要取到小數點後兩位，並包含百分比符號
def format_percent_reads(percent_reads):
    # Ensure percent_reads is a string
    percent_reads = str(percent_reads)

    # Define the regular expression pattern
    pattern = re.compile(r'^\d+\.\d{2}%$')

    # Check if the input matches the pattern
    if pattern.match(percent_reads):
        # Input already matches the format, no adjustment needed
        return percent_reads
    else:
        # Try to extract the numerical part and format it with two decimal places and a percentage sign
        try:
            numerical_part = float(re.search(r'\d+\.\d+', percent_reads).group())
            adjusted_percent_reads = f"{numerical_part:.2f}%"
            return adjusted_percent_reads
        except (AttributeError, ValueError):
            # If extraction or conversion fails, return None (invalid format)
            return percent_reads

# 找 mutations 資料
def extract_mutations(data):
    # 檢查 short_variant_list 是否存在
    if "ResultsReport" not in data or "ResultsPayload" not in data["ResultsReport"] or "variant_report" not in data["ResultsReport"]["ResultsPayload"]:
        return None, None, None, None, None, None, None, None, None, None, None, None, None

    # 從 variant_report 中取得 short_variants 字典，如果不存在則設為空字典
    if data["ResultsReport"]["ResultsPayload"]["variant_report"].get("short_variants", {}) == None:
        return None, None, None, None, None, None, None, None, None, None, None, None, None

    short_variant_list = data["ResultsReport"]["ResultsPayload"]["variant_report"].get("short_variants", {}).get("short_variant", [])

    # 檢查 short_variant_list 是否只有一筆資料，這種情況回傳的格式會是一個字典而不是列表
    if type(short_variant_list) == dict:
        #print("This short_variant_list only have one variant.")
        new_short_variant_list = []
        new_short_variant_list.append(short_variant_list)
        short_variant_list = new_short_variant_list

    # 檢查 short_variant_list 是否有值
    if not short_variant_list:
        return None, None, None, None, None, None, None, None, None, None, None, None, None

    # 建立列表來儲存提取到的變數值
    hugo_symbols = []
    chromosomes = []
    start_positions = []
    end_positions = []
    strands = []
    variant_classifications = []
    reference_alleles = []
    tumor_seq_allele1s = []
    tumor_seq_allele2s = []
    hgvsp_shorts = []
    refseqs = []
    allele_fractions = []
    percent_readss = []

    for dic in short_variant_list:
        hugo_symbol = dic.get("gene")
        chromosome, start_position, end_position = extract_chromosome_start_and_end_position(dic.get("position", ""))
        strand = dic.get("strand")
        variant_classification = dic.get("functional_effect")
        reference_allele, tumor_seq_allele1, tumor_seq_allele2 = extract_ref_allele_tumor_seq_allele1_and_allele2(dic.get("cds_effect", ""))
        hgvsp_short = dic.get("protein_effect")
        # 去掉文字中的換行符號
        if hgvsp_short is not None:
            hgvsp_short = hgvsp_short.replace('\r\n', '').replace('\n', '') 
        refseq = dic.get("transcript")
        # 去掉文字中的換行符號
        if refseq is not None:
            refseq = refseq.replace('\r\n', '').replace('\n', '')
        allele_fraction = dic.get("allele_fraction")
        percent_reads = dic.get("percent_reads")
        # 修正 percent_reads 的值
        if percent_reads == '_' or percent_reads == '-':
            percent_reads = 'NA'
        # 去掉文字中的換行符號
        if percent_reads is not None:
            percent_reads = percent_reads.replace('\r\n', '').replace('\n', '')

        # 提取特定變數，並將其儲存在列表中
        hugo_symbols.append(hugo_symbol)
        chromosomes.append(chromosome)
        start_positions.append(start_position)
        end_positions.append(end_position)
        strands.append(strand)
        variant_classifications.append(variant_classification)
        reference_alleles.append(reference_allele)
        tumor_seq_allele1s.append(tumor_seq_allele1)
        tumor_seq_allele2s.append(tumor_seq_allele2)
        hgvsp_shorts.append(hgvsp_short)
        refseqs.append(refseq)
        allele_fractions.append(allele_fraction)
        percent_readss.append(format_percent_reads(percent_reads))

    return hugo_symbols, chromosomes, start_positions, end_positions, strands, variant_classifications, reference_alleles, tumor_seq_allele1s, tumor_seq_allele2s, hgvsp_shorts, refseqs, allele_fractions, percent_readss

# 印出解析出來的 mutations
def print_extracted_mutations(hugo_symbols, chromosomes, start_positions, end_positions, strands, variant_classifications, reference_alleles, tumor_seq_allele1s, tumor_seq_allele2s, hgvsp_shorts, refseqs, allele_fractions):
    for i in range(len(hugo_symbols)):
        print(f"Index {i}:")
        #print("Hugo Symbol:", hugo_symbols[i])
        #print("Chromosome:", chromosomes[i])
        #print("Start Position:", start_positions[i])
        #print("End Position:", end_positions[i])
        #print("Strand:", strands[i])
        #print("Variant Classification:", variant_classifications[i])
        #print("Reference Allele:", reference_alleles[i])
        #print("Tumor Seq Allele1:", tumor_seq_allele1s[i])
        #print("Tumor Seq Allele2:", tumor_seq_allele2s[i])
        #print("HGVSp Short:", hgvsp_shorts[i])
        #print("RefSeq:", refseqs[i])
        #print("Allele Fraction:", allele_fractions[i])
        #print("----------------------------------------")

# 處理 Dataframe 中 'HGVSp_Short' 這個欄位的函數
def clean_HGVSp_Short(df):
    # 將加號 '+' 替換為 ' plus '
    df['HGVSp_Short'] = df['HGVSp_Short'].str.replace('+', ' plus ')
    # 把欄位裡面的值都轉換成字串
    df['HGVSp_Short'] = df['HGVSp_Short'].astype(str)
    all_strings = df['HGVSp_Short'].apply(lambda x: isinstance(x, str)).all()
    #print("All values in 'HGVSp_Short' are converted to strings:", all_strings)
    # 清洗'HGVSp_Short'欄位的資料，只保留左右小括弧中間的字串
    df['HGVSp_Short'] = df['HGVSp_Short'].apply(lambda x: re.search(r'\((.*?)\)', x).group(1) if re.search(r'\((.*?)\)', x) else x)
    # 清洗'HGVSp_Short'欄位的資料，只保留左右中括弧中間的字串，中括弧裡面若有分號 ';' 則用空格取代
    df['HGVSp_Short'] = df['HGVSp_Short'].apply(lambda x: re.search(r'\[(.*?)\]', x).group(1) if re.search(r'\[(.*?)\]', x) else x)
    # 將分號 ';' 替換為空格 ' '
    df['HGVSp_Short'] = df['HGVSp_Short'].str.replace(';', ' ')
    # 把 'Protein_Change' 的值等同於 'HGVSp_Short' 的值
    df['Protein_Change'] = df['HGVSp_Short']
    return df

# 修改 mutations Dataframe 中有規則可循的欄位
def post_process_mutations(df):
    df['NCBI_Build'] = 'GRCh37'
    df['Strand'] = '+' 
    df['Mutation_Status'] = 'Somatic' 
    df['Variant_Type'] = 'INS'
    return df

# 把 Dataframe 中的 NaN 改成 None
def replace_nan_with_none(df):
    # 使用 replace 函數將 NaN 替換為 None
    df.replace({pd.NA: None, pd.NaT: None}, inplace=True)
    return df

# 儲存 mutations data dataframe
def save_mutations_data(data_mutations_df, output_folder, output_file_data_mutations_extended):
    # 空白值、遺失值填上 'NA'
    data_mutations_df = data_mutations_df.fillna('NA')
    data_mutations_df = data_mutations_df.replace('N/A', 'NA')
    # 拼接完整的路徑名稱
    output_file_path = os.path.join(output_folder, output_file_data_mutations_extended)
    # 儲存 .txt 文件
    data_mutations_df.to_csv(output_file_path, sep='	', encoding='utf-8', index=False)
    #print(f'data mutations is saved to {output_file_path}')

# 找 cna 資料
def extract_cna(data):
    # 檢查 copy_number_alteration_list 是否存在
    if "ResultsReport" not in data or "ResultsPayload" not in data["ResultsReport"] or "variant_report" not in data["ResultsReport"]["ResultsPayload"]:
        return None, None, None, None, None

    # 從 variant_report 中取得 copy_number_alterations 字典，如果不存在則設為空字典
    if data["ResultsReport"]["ResultsPayload"]["variant_report"].get("copy_number_alterations", {}) == None:
        return None, None, None, None, None

    copy_number_alteration_list = data["ResultsReport"]["ResultsPayload"]["variant_report"].get("copy_number_alterations", {}).get("copy_number_alteration", [])

    # 檢查 copy_number_alteration_list 是否只有一筆資料，這種情況回傳的格式會是一個字典而不是列表
    if type(copy_number_alteration_list) == dict:
        #print("This copy_number_alteration_list only have one copy_number.")
        new_copy_number_alteration_list = []
        new_copy_number_alteration_list.append(copy_number_alteration_list)
        copy_number_alteration_list = new_copy_number_alteration_list
    
    # 檢查 copy_number_alteration_list 是否有值
    if not copy_number_alteration_list:
        return None, None, None, None, None

    # 建立列表來儲存提取到的變數值
    copy_numbers = []
    genes = []
    chromosomes = []
    start_positions = []
    end_positions = []

    for dic in copy_number_alteration_list:
        copy_number = dic.get("copy_number")
        gene = dic.get("gene")
        chromosome, start_position, end_position = extract_chromosome_start_and_end_position(dic.get("position", ""))

        # 提取特定變數，並將其儲存在列表中
        copy_numbers.append(copy_number)
        genes.append(gene)
        chromosomes.append(chromosome)
        start_positions.append(start_position)
        end_positions.append(end_position)

    return copy_numbers, genes ,chromosomes, start_positions, end_positions

# 印出解析出來的 cna
def print_extracted_cna(copy_numbers, genes, chromosomes, start_positions, end_positions):
    for i in range(len(copy_numbers)):
        print(f"Index {i}:")
        #print("Copy Number:", copy_numbers[i])
        #print("Gene:", genes[i])
        #print("Chromosome:", chromosomes[i])
        #print("Start Position:", start_positions[i])
        #print("End Position:", end_positions[i])
        #print("----------------------------------------")

# 儲存 cna data (matrix_cna 矩陣) dataframe
def save_cna_data(matrix_cna, output_folder, output_file_data_cna):
    # 拼接完整的路徑名稱
    output_file_path = os.path.join(output_folder, output_file_data_cna)
    # 儲存 .txt 文件
    matrix_cna.to_csv(output_file_path, sep='	', encoding='utf-8', index=True)
    #print(f'data cna is saved to {output_file_path}')

# 對 cna 資料中的 copy_number 做數值上的轉換
def convert_copy_number(copy_number):
    try:
        # 將 copy_number 轉換成浮點數
        copy_number_float = float(copy_number)

        # 根據條件設定 save_copy_number
        if copy_number_float > 4:
            save_copy_number = 2
        elif 2 <= copy_number_float <= 4:
            save_copy_number = 1
        elif copy_number_float == 2 or copy_number_float is None:
            save_copy_number = 0
        elif copy_number_float == 1:
            save_copy_number = -1
        elif copy_number_float == 0:
            save_copy_number = -2
        else:
            # 如果不符合上述條件，設定為預設值
            save_copy_number = 'NA'  # 可以修改為其他適當的預設值
    except (ValueError, TypeError):
        # 轉換失敗時的處理
        save_copy_number = 'NA'  # 可以修改為其他適當的預設值
        #print('Failed to convert copy number.')
        #print(f'Type of compy number: {type(copy_number)}')
    
    return save_copy_number

# 印出經過轉換的 copy_number，儲存成 save_copy_number
def print_save_copy_number(data_cna_df):
    for i in range(len(data_cna_df)):
        print('index:', i)
        #print('patient_id:', data_cna_df.PATIENT_ID[i])
        #print('gene:', data_cna_df.Gene[i])
        #print('copy_number:', data_cna_df.Copy_Number[i])
        #print('save_copy_number:', convert_copy_number(data_cna_df.Copy_Number[i]))
        #print('----------------------------------')

# 提取 description 欄位中的 site1_region_number 和 site2_region_number
def extract_site1_region_number_and_site2_region_number(description):
    # 若檢查到 description 為 None 則都回傳 None
    if description is None:
        #print("Description is None.")
        return None, None, None, None
    
    # 提取 position 欄位中的 chromosome、start_position 和 end_position 的正則表達式
    pattern = re.compile(r'([A-Za-z0-9]+)_([A-Za-z0-9]+)\.([A-Za-z])(\d+)([A-Za-z])(\d+)')
    matches = pattern.findall(description)
    if matches:
        for match in matches:
            gene1 = match[0]
            site1_region_number = match[3]
            gene2 = match[1]
            site2_region_number = match[5]
            return gene1, site1_region_number, gene2, site2_region_number
    else:
        #print("No site region number found.")
        return None, None, None, None

# 找 sv 資料
def extract_sv(data):
    # 檢查 structural_variant_list 是否存在
    if "ResultsReport" not in data or "ResultsPayload" not in data["ResultsReport"] or "variant_report" not in data["ResultsReport"]["ResultsPayload"]:
        return None, None, None, None, None, None, None, None, None

    # 從 variant_report 中取得 rearrangements 字典，如果不存在則設為空字典
    if data["ResultsReport"]["ResultsPayload"]["variant_report"].get("rearrangements", {}) == None:
        return None, None, None, None, None, None, None, None, None

    structural_variant_list = data["ResultsReport"]["ResultsPayload"]["variant_report"].get("rearrangements", {}).get("rearrangement", [])

    # 檢查 structural_variant_list 是否只有一筆資料，這種情況回傳的格式會是一個字典而不是列表
    if type(structural_variant_list) == dict:
        #print("This structural_variant_list only have one structural variant gene.")
        new_structural_variant_list = []
        new_structural_variant_list.append(structural_variant_list)
        structural_variant_list = new_structural_variant_list
    
    # 檢查 copy_number_alteration_list 是否有值
    if not structural_variant_list:
        return None, None, None, None, None, None, None, None, None

    # 建立列表來儲存提取到的變數值
    descriptions = []
    gene1s = []
    site1_region_numbers = []
    gene2s = []
    site2_region_numbers = []
    other_genes = []
    targeted_genes = []
    pos1s = []
    pos2s = []

    for dic in structural_variant_list:
        description = dic.get("description")
        gene1, site1_region_number, gene2, site2_region_number = extract_site1_region_number_and_site2_region_number(dic.get("description", ""))
        # 把 'N/A' 替換成 'NA'
        if gene1 == 'N/A':
            gene1 = gene1.replace('N/A', 'NA')
        if site1_region_number == 'N/A':
            site1_region_number = site1_region_number.replace('N/A', 'NA')
 	          # site1_region_number = '2'
        if gene2 == 'N/A':
            gene2 = gene2.replace('N/A', 'NA')
        if site2_region_number == 'N/A':
            site2_region_number = site2_region_number.replace('N/A', 'NA')
	          # site2_region_number = '4' 
        other_gene = dic.get("other_gene")
        if other_gene == 'N/A':
            other_gene = other_gene.replace('N/A', 'NA')
        targeted_gene = dic.get("targeted_gene")
        if targeted_gene == 'N/A':
            targeted_gene = targeted_gene.replace('N/A', 'NA')
        pos1 = dic.get("pos1")
        pos2 = dic.get("pos2")

        # 提取特定變數，並將其儲存在列表中
        descriptions.append(description)
        gene1s.append(gene1)
        site1_region_numbers.append(site1_region_number)
        gene2s.append(gene2)
        site2_region_numbers.append(site2_region_number)
        other_genes.append(other_gene)
        targeted_genes.append(targeted_gene)
        pos1s.append(pos1)
        pos2s.append(pos2)
    #return 'NA', 'NA', 'NA', 'NA', 'NA', other_genes ,targeted_genes, 'NA', 'NA'
    return descriptions, gene1s, site1_region_numbers, gene2s, site2_region_numbers, other_genes ,targeted_genes, pos1s, pos2s
    #return descriptions, gene1s, site1_region_numbers, gene2s, site2_region_numbers, other_genes ,targeted_genes, pos1s, pos2s

# 修改 sv Dataframe 中有規則可循的欄位
def post_process_sv(df):
    df['SV_Status'] = 'Somatic' 
    df['Site1_Region'] = 'exon' 
    df['Site2_Region'] = 'exon' 
    return df

# 儲存  data sv 的 dataframe
def save_sv_data(data_sv_df, output_folder, output_file_data_sv):
    # 空白值、遺失值填上 'NA'
    data_sv_df = data_sv_df.fillna('NA')
    data_sv_df = data_sv_df.replace('N/A', 'NA')
    data_sv_df = data_sv_df.drop_duplicates()
    # 拼接完整的路徑名稱
    output_file_path = os.path.join(output_folder, output_file_data_sv)
    # 儲存 .txt 文件
    data_sv_df.to_csv(output_file_path, sep='	', encoding='utf-8', index=False)
    #print(f'data sv is saved to {output_file_path}')

# # 定義一個函數來處理 Tumor_Purity 欄位的值(從 ngs.csv)
# def format_tumor_purity(value):
#     # 檢查值是否為 NaN
#     if pd.isna(value):
#         return value
#     try:
#         # 嘗試將值轉換為浮點數
#         float_value = float(value)
#         # 使用 f-string 格式化帶有百分比符號的字串，四捨五入到小數點第一位
#         formatted_value = f"{float_value:.1f} %"
#         return formatted_value
#     except (ValueError, TypeError):
#         # 如果轉換失敗或值為非數字，則返回原始值
#         return value

# # 定義一個函數來處理 Tumor_Purity 欄位的值(從 json)
# def format_tumor_purity(value):
#     # 檢查值是否為 NaN
#     if pd.isna(value):
#         return value
#     try:
#         # 嘗試將值轉換為浮點數
#         float_value = float(value)
#         # 如果是整數，直接返回整數
#         if float_value.is_integer():
#             return int(float_value)
#         # 否則轉換為整數並四捨五入
#         integer_value = round(float_value)
#         return integer_value
#     except (ValueError, TypeError):
#         # 如果轉換失敗或值為非數字，則返回原始值
#         return value

# 定義一個函數來處理 Tumor_Purity 欄位的值(從 json)
def format_tumor_purity(value):
    # 檢查值是否為 NaN
    if pd.isna(value):
        return value
    try:
        # 先檢查是否為數字
        if isinstance(value, (int, float)):
            # 如果是數字，直接返回整數或四捨五入後的整數
            if float(value).is_integer():
                return int(value)
            else:
                return round(float(value))
        else:
            # 如果是字符串，移除百分比符號並將值轉換為浮點數
            value = value.replace('%', '')
            float_value = float(value)
            # float_value = float(value.rstrip('%'))
            # 如果是整數，直接返回整數
            if float_value.is_integer():
                return int(float_value)
            # 否則轉換為整數並四捨五入
            integer_value = round(float_value)
            return integer_value
    except (ValueError, TypeError):
        # 如果轉換失敗或值為非數字，則返回原始值
        return 'NA'

# 定義一個函數，將底線後面的部分用小括號刮起來
def add_parentheses(input_string):
    # 使用正則表達式找到底線後的部分
    pattern = re.compile(r'_(\w+)')
    match = pattern.search(input_string)
    
    # 如果找到匹配，則加上小括號
    if match:
        replacement = f"_{match.group(1)}"
        result = input_string.replace(match.group(), replacement)
        result = result.replace('_', '_(') + ')'
        return result
    else:
        return input_string

# ======================================== 讀取 config.json 中定義的變數 ========================================
# 讀取JSON文件
with open('config.json', 'r', encoding='utf-8') as config_file:
    config = json.load(config_file)

# 定義資料讀取與儲存路徑
output_folder = config["output_folder"]
output_folder_case = config["output_folder_case"]
# check_data_folder = config["check_data_folder"]

# 定義所有的 meta 檔的儲存名稱
output_file_meta_study = config["output_file_meta_study"]
output_file_meta_mutations_extended = config["output_file_meta_mutations_extended"]
output_file_meta_clinical_patient = config["output_file_meta_clinical_patient"]
output_file_meta_clinical_sample = config["output_file_meta_clinical_sample"]
output_file_meta_cna = config["output_file_meta_cna"]
output_file_meta_sv = config["output_file_meta_sv"]
output_file_meta_ctype = config["output_file_meta_ctype"]

# 定義所有的 data 檔的儲存名稱
output_file_data_mutations_extended = config["output_file_data_mutations_extended"]
output_file_data_clinical_patient = config["output_file_data_clinical_patient"]
output_file_data_clinical_sample = config["output_file_data_clinical_sample"]
output_file_data_cna = config["output_file_data_cna"]
output_file_data_sv = config["output_file_data_sv"]
output_file_data_ctype = config["output_file_data_ctype"]

# 定義所有的 case lists 檔的儲存名稱
output_file_cases_sequenced = config["output_file_cases_sequenced"]
output_file_cases_cna = config["output_file_cases_cna"]
output_file_cases_sv = config["output_file_cases_sv"]

# 定義所有的 case lists 檔裡面需要用到的變數 (sequenced, cna, sv)
case_list_name_se =  config["case_list_name_se"]
case_list_description_se = config["case_list_description_se"]
case_list_category_se = config["case_list_category_se"]
case_list_name_cna = config["case_list_name_cna"]
case_list_description_cna = config["case_list_description_cna"]
case_list_category_cna = config["case_list_category_cna"]
case_list_name_sv = config["case_list_name_sv"]
case_list_description_sv = config["case_list_description_sv"]
case_list_category_sv = config["case_list_category_sv"]

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

# 定義 cancer type 檔裡面需要的變數
type_of_cancer_full_name = config["type_of_cancer_full_name"] # 定義的癌症名稱全名
show_color = config["show_color"] # 要顯現出來的顏色
type_of_cancer_parent = config["type_of_cancer_parent"] # 癌症種類的母類別(有限制，要去原始程式碼改)

# 指定目錄路徑，該路徑放置所有的病人病歷 json 檔
directory_path = config["directory_path"]

# 讀取常見癌症種類字典
common_cancer_type_dict = config["common_cancer_type"]

# 取得 ngs csv url 不變的部分
base_ngs_url = config["base_ngs_url"]
base_ngs_directory = config["base_ngs_directory"]

# 取得 json 的 url 以及儲存 json 的 url
base_download_json_url = config["base_download_json_url"]
save_json_directory =  config["save_json_directory"]
save_wrong_json_directory = config["save_wrong_json_directory"]

# 如果資料夾(路徑)不存在，則創建資料夾
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
# 如果資料夾(路徑)不存在，則創建資料夾
if not os.path.exists(output_folder_case):
    os.makedirs(output_folder_case)
# # 如果資料夾(路徑)不存在，則創建資料夾
# if not os.path.exists(check_data_folder):
#     os.makedirs(check_data_folder)
# 如果資料夾(路徑)不存在，則創建資料夾
if not os.path.exists(save_json_directory):
    os.makedirs(save_json_directory)
# 如果資料夾(路徑)不存在，則創建資料夾
if not os.path.exists(save_wrong_json_directory):
    os.makedirs(save_wrong_json_directory)

# ======================================== 自動生成 meata 檔 ========================================
# 拼接完整的路徑名稱(meta_study)
output_path = os.path.join(output_folder, output_file_meta_study)
# 寫到 meta
with open(output_path, "w", encoding="utf-8") as file:
    file.write(f"type_of_cancer: {type_of_cancer}\n")
    file.write(f"cancer_study_identifier: {cancer_study_identifier}\n")
    file.write(f"name: {name}\n")
    file.write(f"description: {description}\n")
    file.write(f"add_global_case_list: {add_global_case_list}\n")
    file.write(f"groups: {groups}\n")
#print(f"meta study is saved to '{output_path}'")

# 拼接完整的路徑名稱(meta_mutations_extended)
output_path = os.path.join(output_folder, output_file_meta_mutations_extended)
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
#print(f"meta mutations extended is saved to '{output_path}'")

# 拼接完整的路徑名稱(meta_clinical_patient)
output_path = os.path.join(output_folder, output_file_meta_clinical_patient)
# 寫到 meta
with open(output_path, "w", encoding="utf-8") as file:
    file.write(f"cancer_study_identifier: {cancer_study_identifier}\n")
    file.write(f"genetic_alteration_type: {genetic_alteration_type_pa}\n")
    file.write(f"datatype: {datatype_pa}\n")
    file.write(f"data_filename: {data_filename_pa}\n")
#print(f"meta clinical patient is saved to '{output_path}'")

# 拼接完整的路徑名稱(meta_clinical_sample)
output_path = os.path.join(output_folder, output_file_meta_clinical_sample)
# 寫到 meta
with open(output_path, "w", encoding="utf-8") as file:
    file.write(f"cancer_study_identifier: {cancer_study_identifier}\n")
    file.write(f"genetic_alteration_type: {genetic_alteration_type_sa}\n")
    file.write(f"datatype: {datatype_sa}\n")
    file.write(f"data_filename: {data_filename_sa}\n")
#print(f"meta clinical sample is saved to '{output_path}'")

# 拼接完整的路徑名稱(meta_cna)
output_path = os.path.join(output_folder, output_file_meta_cna)
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
#print(f"meta cna is saved to '{output_path}'")

# 拼接完整的路徑名稱(meta_sv)
output_path = os.path.join(output_folder, output_file_meta_sv)
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
#print(f"meta sv is saved to '{output_path}'")

# 拼接完整的路徑名稱(meta_cancer_type)
output_path = os.path.join(output_folder, output_file_meta_ctype)
# 寫到 meta
with open(output_path, "w", encoding="utf-8") as file:
    file.write(f"genetic_alteration_type: {genetic_alteration_type_ctype}\n")
    file.write(f"datatype: {datatype_ctype}\n")
    file.write(f"data_filename: {data_filename_ctype}\n")
#print(f"meta cancer type is saved to '{output_path}'")

# 拼接完整的路徑名稱(meta_cancer_type)
output_path = os.path.join(output_folder, output_file_data_ctype)
# 將文字寫入 .txt 文件檔中
with open(output_path, "w", encoding="utf-8") as file:
    file.write(f"{type_of_cancer}	")
    file.write(f"{type_of_cancer_full_name}	")
    file.write(f"{show_color}	")
    file.write(f"{type_of_cancer_parent}\n")
#print(f"data cancer type is saved to '{output_path}'")
#print("\n")

# ======================================== 測試是否可以讀取到從指令下的參數 (以下尚未開發完成) ========================================
# #print("############################# Start getting the NGS csv from the URL #############################")
# # 創建 ArgumentParser 對象
# parser = argparse.ArgumentParser(description='My Script Description')

# try:
#     # 添加命令行參數
#     parser.add_argument('Arg', type=str, help='This argument, please enter the file name of NGS csv, which is used to map and generate txt files for uploading to cbio portal.')

#     # 解析命令行參數
#     args = parser.parse_args()

#     # 訪問參數值
#     #print(f"Your csv file name is: {args.Arg}")
#     #print("\n")

# except argparse.ArgumentError as e:
#     #print(f"Error: {e}")
#     parser.#print_help()

# # 開始抓 NGS csv 檔案
# #print("Start getting the NGS csv in save_csv.")
# # 添加命令行參數
# parser.add_argument('Arg', type=str, help='This argument, please enter the file name of NGS csv, which is used to map and generate txt files for uploading to cbio portal.')
# # 解析命令行參數
# args = parser.parse_args()
# # 訪問參數值
# csv_file_name = args.Arg
# #print(f"Your csv file name is: {csv_file_name}")
# # 組合成完整的 URL
# csv_directory = f"{base_ngs_directory}/{csv_file_name}"
# #print(f'This is the directory used to get NGS csv: {csv_directory}')

# # 組合成完整的 URL ############################################################################
# csv_directory = base_ngs_url
# #print(f'This is the directory used to get NGS csv: {csv_directory}')

# # 使用 requests 下載 CSV 檔案
# response = requests.get(csv_directory)

# # 檢查是否成功取得檔案
# if response.status_code == 200:
#     #print("Check Point -0: Successfully get the NGS csv file.")
#     # 根據文件副檔名格式不同選擇讀取的方式，讀取成 Dataframe
#     if csv_directory.endswith('.csv'):
#         ngs_df = pd.read_csv(BytesIO(response.content))
#     elif csv_directory.endswith('.xlsx'):
#         ngs_df = pd.read_excel(BytesIO(response.content))
#     else:
#         #print("Unsupported file type.")
#         ngs_df = pd.DataFrame()

#     # 印出 DataFrame
#     #print(ngs_df, "\n")

#     # 獲取當前 UTC 時間
#     utc_now = datetime.utcnow()

#     # 設置時區為台灣
#     taipei_timezone = pytz.timezone('Asia/Taipei')
#     taipei_now = utc_now.replace(tzinfo=pytz.utc).astimezone(taipei_timezone)

#     # 將日期和時間轉換為字串格式
#     timestamp_str = taipei_now.strftime("%Y-%m-%d_%H%M%S")

#     # 檢查目錄是否存在
#     if not os.path.exists(base_ngs_directory):
#         # 如果不存在，創建目錄
#         os.makedirs(base_ngs_directory)

#     # 將 DataFrame 儲存成 CSV 文件
#     ngs_df.to_csv(f'{base_ngs_directory}/{timestamp_str}.csv', index=False, encoding='utf-8')
#     #print(f"Saved DataFrame to {base_ngs_directory}")

# else:
#     #print(f"Failed to fetch the CSV file. Status code: {response.status_code}", "\n")

# # # 讀取 NGS CSV 檔案到 DataFrame
# # ngs_df = pd.read_csv(csv_directory)

# # 顯示 DataFrame
# #print('Check Point -1: Check ngs.csv is successfully download and converted into dataframe.')
# #print(f'original ngs df columns: {ngs_df.columns}')
# #print(f'original ngs df shape: {ngs_df.shape}')
# #print(f'original ngs df: {ngs_df}', "\n")

# # 選擇要保留的欄位並複製到新的 DataFrame
# new_ngs_df = ngs_df[["分生號碼", "報告號碼", "病患姓名", "檢體別", "蠟塊號", "Tumor purity %", "Tumor type", "Diagnosis", "檢測項目", "臨床主治醫師", "報告日期"]].copy()

# # 修改欄位名稱
# new_ngs_df = new_ngs_df.rename(columns={
#     "分生號碼": "MP_Number",
#     "報告號碼": "Path_Number",
#     "病患姓名": "PATIENT_NAME",
#     "檢體別": "SPECIMEN",
#     "蠟塊號": "Cell_Block_Number",
#     "Tumor purity %": "Tumor_Purity",
#     "Tumor type": "CANCER_TYPE",
#     "Diagnosis": "DIAGNOSIS",
#     "檢測項目": "REPORT_ASSAY",
#     "臨床主治醫師": "Attending_Physician",
#     "報告日期": "DATE"
# })

# # 新增一個叫 PATIENT_ID 的欄位
# new_ngs_df["PATIENT_ID"] = new_ngs_df["Path_Number"] + "_" + "(" + new_ngs_df["MP_Number"] + ")"
# # new_ngs_df["PATIENT_ID"] = new_ngs_df["Path_Number"] + "_" + new_ngs_df["MP_Number"]

# # 重新排序欄位順序
# new_ngs_df = new_ngs_df[["PATIENT_ID", "Path_Number", "MP_Number", "SPECIMEN", "Cell_Block_Number", "Tumor_Purity", "CANCER_TYPE", "DIAGNOSIS", "REPORT_ASSAY", "Attending_Physician", "DATE"]]
# new_ngs_df["Tumor_Purity"] = new_ngs_df["Tumor_Purity"].apply(format_tumor_purity)

# # 顯示 DataFrame
# #print('Check Point -2: Check the original dataframe is converted into new dataframe by removing and renaming some columns.')
# #print(f'new ngs df columns: {new_ngs_df.columns}')
# #print(f'new ngs df shape: {new_ngs_df.shape}')
# #print(f'new ngs df df: {new_ngs_df}', "\n")

# # 假設 new_ngs_df 是包含 "PATIENT_ID" 的 DataFrame
# for patient_id in new_ngs_df["PATIENT_ID"]:
#     # 組合下載的 URL
#     download_url = f"{base_download_json_url}/{patient_id}.json"
#     #print(download_url)
#     # 組合保存的檔案路徑
#     save_path = os.path.join(save_json_directory, f"{patient_id}.json")
#     #print(save_path)
#     # 檢查檔案是否已存在
#     if os.path.exists(save_path):
#         #print(f"JSON for patient_id {patient_id} already exists. Skipping download.")
#     else:
#         #print(f"You need to go to {download_url} and download this json file to continue the following work.")
#         # 下載 JSON 檔案
#         response = requests.get(download_url)
#         # 檢查下載是否成功 (HTTP 狀態碼 200 表示成功)
#         if response.status_code == 200:
#             # 寫入 JSON 檔案
#             with open(save_path, 'w', encoding='utf-8') as file:
#                 # file.write(response.content)
#                 file.write(response.text)
#             #print(f"Downloaded and saved: {save_path}")
#         else:
#             #print(f"Failed to download for patient_id: {patient_id}")
#
# # 使用正規表達式移除左右兩邊的括號(為了最後和 data_sample_df 參照，所以調整 'PATIENT_ID' 的格式)
# new_ngs_df["PATIENT_ID"] = new_ngs_df["PATIENT_ID"].apply(lambda x: re.sub(r'[()]', '', str(x)))
# # 到這裡結束 ############################################################################

# # 發送 HTTP 請求取得 JSON 檔案列表
# response = requests.get(base_download_json_url)

# # 檢查 HTTP 狀態碼是否為 200 OK
# if response.status_code == 200:
#     try:
#         # 使用 try-except 捕捉 JSONDecodeError
#         json_files = response.json()
#     except json.JSONDecodeError:
#         #print("無法解析 JSON 格式的內容。")
#         json_files = []
# else:
#     #print(f"HTTP 請求失敗，狀態碼：{response.status_code}")
#     json_files = []

# # 初始化計數變數
# download_count = 0

# # 迭代下載每個 JSON 檔案
# for json_file in json_files:
#     file_url = base_download_json_url + json_file
#     file_path = os.path.join(save_json_directory, json_file)

#     # 下載 JSON 檔案
#     response = requests.get(file_url)
    
#     # 檢查 HTTP 狀態碼是否為 200 OK
#     if response.status_code == 200:
#         try:
#             # 使用 try-except 捕捉 JSONDecodeError
#             with open(file_path, 'w') as file:
#                 json.dump(response.json(), file)
#             download_count += 1
#             #print(f"下載並儲存 {json_file} 完成。")
#         except json.JSONDecodeError:
#             #print(f"無法解析 {json_file} 的 JSON 格式內容。")
#     else:
#         #print(f"HTTP 請求失敗，狀態碼：{response.status_code}，無法下載 {json_file}。")

# # 顯示下載總數
# #print(f"總共下載了 {download_count} 個 JSON 檔案。")

#print("############################# Start getting report json files in directory 'json' #############################")
# # 取得目錄下所有檔案
json_files = [f for f in os.listdir(base_download_json_url) if f.endswith('.json')]
#print(f'There are {len(json_files)} in this directory.')
save_json_directory = base_download_json_url
#print("############################# Finish getting report json files in directory 'json' ############################")

# #print("############################# Finish getting the NGS csv from the URL ############################")
#print("\n")
# ======================================== 建立空白的 dataframe 來儲存 parsing 出來的資訊 ========================================
# 用列表建立 Dataframe 表格中所需要的欄位
data_sample_list = ["PATIENT_ID", "SAMPLE_ID", "PATH_NUMBER", "MP_NUMBER", "CELL_BLOCK_NUMBER", "TUMOR_TYPE", "TUMOR_PURITY", "REPORT_ASSAY", "REPORT_DATE"]
data_patient_list = ["PATIENT_ID", "DIAGNOSIS", "STAGE", "AGE", "GENDER", "ETHNICITY"]
data_mutations_list = ["Hugo_Symbol", "NCBI_Build", "Chromosome", "Start_Position", "End_Position", "Strand", "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Allele_Fraction", "Tumor_Sample_Barcode", "Mutation_Status", "HGVSp_Short", "RefSeq", "Protein_Change", "Percent_Reads"]
data_cna_list = ["PATIENT_ID", "Copy_Number", "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position"]
data_sv_list = ["Sample_Id", "SV_Status", "Site1_Hugo_Symbol", "Site1_Region", "Site1_Region_Number", "Site2_Hugo_Symbol", "Site2_Region", "Site2_Region_Number"]
# 建立空白的 DataFrame
data_sample_df = pd.DataFrame(columns=data_sample_list)
data_patient_df = pd.DataFrame(columns=data_patient_list)
data_mutations_df = pd.DataFrame(columns=data_mutations_list)
data_cna_df = pd.DataFrame(columns=data_cna_list)
data_sv_df = pd.DataFrame(columns=data_sv_list)

# ======================================== 讀取病人病例報告的 json 檔 ========================================
# 呼叫函數並取得所有 JSON 檔案的資料
# all_json_data, filenames, path_numbers, mp_numbers = read_json_files(directory_path)
all_json_data, filenames, path_numbers, mp_numbers, problematic_files = read_json_files(save_json_directory)

# 印出 path number 跟 mp number
# #print('All json reports:', all_json_data)
#print('Check Point - 2: Check the json reports can be successfully download.')
# #print('Filename of json reports:', filenames)
# #print('Path number of json reports:', path_numbers)
# #print('MP number of json reports:', mp_numbers)
#print('These are json reports which can not be read:', problematic_files)
#print("\n")

# ======================================== 開始進行資料 parsing ========================================
# 印出此行開始將病歷報告的 json 檔轉換至可以上傳到 cbio portal 的資料格式
#print("################ Start generating data and meta file for uploading to cbio portal ################")
for i in range(len(all_json_data)):
    data = all_json_data[i]
    path_number = path_numbers[i]
    mp_number = mp_numbers[i]
    # 印出每一個 json report 的內容
    # #print(f'now processed data: {data}')
    #print(f'now processed data path number: {path_number}')
    #print(f'now processed data mp number: {mp_number}')
    #print("---------------------------------------------------------------------------------")
    # ================================= 開始做 patient 和 sample 資料 =================================
    # 用函數去找指定的值
    patient_id = find_patient_id(path_number, mp_number)
    sample_id =  mp_number
    cancer_type = find_cancer_type(data)
    # cancer_type_detailed = find_cancer_type_detailed(data)
    diagnosis = find_diagnosis(data)
    report_test_assay = find_test_assay(data)
    percent_tumor_nuclei = find_percent_tumor_nuclei(data)
    percent_tumor_nuclei = format_tumor_purity(percent_tumor_nuclei) # 確保回傳的是 int
    # # 判斷是否是數字
    # if percent_tumor_nuclei.isdigit():
    #     # 如果是數字，直接轉換為整數
    #     percent_tumor_nuclei = int(percent_tumor_nuclei)
    #     # #print("是數字:", percent_tumor_nuclei)
    # else:
    #     # 如果不是數字，去掉百分比符號再轉換為整數
    #     percent_tumor_nuclei = int(percent_tumor_nuclei.rstrip('%'))
    #     # #print("不是數字，去掉百分比符號後的數字:", percent_tumor_nuclei)
    gender = find_gender(data)
    # 檢查找到的值
    #print(f'patient_id: {patient_id}')
    #print(f'sample_id: {sample_id}')
    #print(f'cancer_type: {cancer_type}')
    # #print(f'cancer_type_detailed: {cancer_type_detailed}')
    #print(f'diagnosis: {diagnosis}')
    #print(f'report_assay: {report_test_assay}')
    #print(f'tumor_purity: {percent_tumor_nuclei}')
    #print(f'gender: {gender}')
    
    # 建立字典，將變數的值對應到 DataFrame 的列名
    new_sample_data = {
        "PATIENT_ID": find_patient_id(path_number, mp_number),
        "SAMPLE_ID": find_patient_id(path_number, mp_number),
        "PATH_NUMBER": path_number,
        "MP_NUMBER": mp_number,
        "TUMOR_TYPE": find_cancer_type(data), # find_cancer_category((lambda data: find_diagnosis(data) if (result := find_cancer_type(data)) is None else result)(data), common_cancer_type_dict)
        "TUMOR_PURITY": percent_tumor_nuclei, # find_percent_tumor_nuclei(data)
        "REPORT_ASSAY": find_test_assay(data)
    }

    new_patient_data = {
        "PATIENT_ID": find_patient_id(path_number, mp_number),
        "DIAGNOSIS": find_diagnosis(data),
        "STAGE": None,
        "AGE": None,
        "GENDER": find_gender(data),
        "ETHNICITY": None
    }

    # 將字典轉換為 DataFrame
    new_data_sample_df = pd.DataFrame([new_sample_data])
    new_data_patient_df = pd.DataFrame([new_patient_data])
    # 將新的數據逐筆添加到原有的 DataFrame
    data_sample_df = pd.concat([data_sample_df, new_data_sample_df], ignore_index=True)
    data_patient_df = pd.concat([data_patient_df, new_data_patient_df], ignore_index=True)
    #print("---------------------- Finish parsing patient and sample ------------------------")

    # ================================= 開始做 mutations 資料 =================================
    # 使用函數找 mutations
    hugo_symbols, chromosomes, start_positions, end_positions, strands, variant_classifications, reference_alleles, tumor_seq_allele1s, tumor_seq_allele2s, hgvsp_shorts, refseqs, allele_fractions, percent_readss = extract_mutations(data)

    if hugo_symbols is not None:
        # 檢查列表長度是否一致
        list_lengths = [len(lst) for lst in [hugo_symbols, chromosomes, start_positions, end_positions, strands, variant_classifications, reference_alleles, tumor_seq_allele1s, tumor_seq_allele2s, hgvsp_shorts, refseqs, allele_fractions, percent_readss]]
        
        if all(length == list_lengths[0] for length in list_lengths):
            # 確認新資料的數量
            rows = len(hugo_symbols)
            #print(f'There are {rows} gene mutations data.')
            
            # 創建一個字典，鍵是 DataFrame 的欄位名稱，值是新資料
            new_mutations_data = {
                "Hugo_Symbol": hugo_symbols,
                "Chromosome": chromosomes,
                "Start_Position": start_positions,
                "End_Position": end_positions,
                "Strand": strands,
                "Variant_Classification": variant_classifications,
                "Reference_Allele": reference_alleles,
                "Tumor_Seq_Allele1": tumor_seq_allele1s,
                "Tumor_Seq_Allele2": tumor_seq_allele2s,
                "HGVSp_Short": hgvsp_shorts,
                "Protein_Change": hgvsp_shorts,
                "RefSeq": refseqs,
                "Allele_Fraction": allele_fractions,
                "Tumor_Sample_Barcode": find_patient_id(path_number, mp_number),
                "Percent_Reads": percent_readss
            }
            
            # 將新資料串聯到現有 DataFrame
            data_mutations_df = pd.concat([data_mutations_df, pd.DataFrame(new_mutations_data)], ignore_index=True)
            # 使用函數處理 'HGVSp_Short' 欄位中的值
            data_mutations_df = clean_HGVSp_Short(data_mutations_df)
            # 把 NaN 改成 None
            data_mutations_df = replace_nan_with_none(data_mutations_df)
            # 修改 mutations Dataframe 中有規則可循的欄位
            data_mutations_df= post_process_mutations(data_mutations_df)
            # 印出分隔列代表製作 mutations 階段結束
            #print("---------------------- Finish parsing mutations ---------------------------------")
        
        else:
            print("Warning: The lengths of the lists for mutations are not consistent.")
            # 印出分隔列代表製作 mutations 階段結束
            #print("---------------------- Finish parsing mutations ---------------------------------")
    else:
        #print("Warning: This data has no mutations.")
        # 印出分隔列代表製作 mutations 階段結束
        #print("---------------------- Finish parsing mutations ---------------------------------")
        pass

    # ================================= 開始做 cna 資料 =================================
    # 使用函數找 cna
    copy_numbers, genes, chromosomes, start_positions, end_positions = extract_cna(data)

    if copy_numbers is not None:
        # 檢查列表長度是否一致
        list_lengths = [len(lst) for lst in [copy_numbers ,genes, chromosomes, start_positions, end_positions]]

        if all(length == list_lengths[0] for length in list_lengths):
            # 確認新資料的數量
            rows = len(copy_numbers)
            #print(f'There are {rows} gene cna data.')
            
            # 創建一個字典，鍵是 DataFrame 的欄位名稱，值是新資料
            new_cna_data = {
                "Copy_Number": copy_numbers,
                "Hugo_Symbol": genes,
                "Chromosome": chromosomes,
                "Start_Position": start_positions,
                "End_Position": end_positions,
                "PATIENT_ID": find_patient_id(path_number, mp_number)
            }
            
            # 將新資料串聯到現有 DataFrame
            data_cna_df = pd.concat([data_cna_df, pd.DataFrame(new_cna_data)], ignore_index=True)
            # 把 NaN 改成 None
            data_cna_df = replace_nan_with_none(data_cna_df)
            # 印出分隔列代表製作 cna 階段結束
            #print("---------------------- Finish parsing cna ---------------------------------------")
        
        else:
            print("Warning: The lengths of the lists for cna are not consistent.")
            # 印出分隔列代表製作 cna 階段結束
            #print("---------------------- Finish parsing cna ---------------------------------------")
    else:
        #print("Warning: This data has no copy number variants.")
        # 印出分隔列代表製作 cna 階段結束
        #print("---------------------- Finish parsing cna ---------------------------------------")
        pass

    # ================================= 開始做 cv 資料 ==================================
    # 使用函數找 sv
    descriptions, gene1s, site1_region_numbers, gene2s, site2_region_numbers, other_genes ,targeted_genes, pos1s, pos2s = extract_sv(data)

    if descriptions is not None:
        # 檢查列表長度是否一致
        list_lengths = [len(lst) for lst in [descriptions, other_genes ,targeted_genes, pos1s, pos2s]]
        
        if all(length == list_lengths[0] for length in list_lengths):
            # 確認新資料的數量
            rows = len(descriptions)
            #print(f'There are {rows} gene sv data.')
            
            # 創建一個字典，鍵是 DataFrame 的欄位名稱，值是新資料
            # Sample_Id	SV_Status	Site1_Hugo_Symbol	Site1_Region	Site1_Region_Number	Site2_Hugo_Symbol	Site2_Region	Site2_Region_Number
            new_sv_data = {
                "Site1_Hugo_Symbol": other_genes,
                "Site1_Region_Number": site1_region_numbers,
                "Site2_Hugo_Symbol": targeted_genes,
                "Site2_Region_Number": site2_region_numbers,
                "Sample_Id": find_patient_id(path_number, mp_number)
            }
            
            # 將新資料串聯到現有 DataFrame
            data_sv_df = pd.concat([data_sv_df, pd.DataFrame(new_sv_data)], ignore_index=True)
            # 把 NaN 改成 None
            data_sv_df = replace_nan_with_none(data_sv_df)
            # 修改 sv Dataframe 中有規則可循的欄位
            data_sv_df = post_process_sv(data_sv_df)
            # 印出分隔列代表製作 sv 階段結束
            #print("---------------------- Finish parsing sv ----------------------------------------")
        
        else:
            print("Warning: The lengths of the lists for sv are not consistent.")
            # 印出分隔列代表製作 sv 階段結束
            #print("---------------------- Finish parsing sv ----------------------------------------")
    else:
        #print("Warning: This data has no structural variants.")
        # 印出分隔列代表製作 sv 階段結束
        #print("---------------------- Finish parsing sv ----------------------------------------")
        pass
        
# ================================= 開始做 cna 矩陣 =================================
# 建立一個全零矩陣，索引為'Hugo_Symbol'，欄位為 'PATIENT_ID'
matrix_cna = pd.DataFrame(np.zeros((len(data_cna_df['Hugo_Symbol'].unique()), len(data_cna_df['PATIENT_ID'].unique())), dtype=int),
                     index=data_cna_df['Hugo_Symbol'].unique(),
                     columns=data_cna_df['PATIENT_ID'].unique())

# 重命名索引為欄位名稱
matrix_cna.index.name = 'Hugo_Symbol'

for j in range(len(data_cna_df)):
    patient_id = data_cna_df.PATIENT_ID[j]
    copy_number = data_cna_df.Copy_Number[j]
    save_copy_number = convert_copy_number(copy_number)
    gene = data_cna_df.Hugo_Symbol[j]

    # 將 save_copy_number 值填入 matrix 中對應的位置
    matrix_cna.loc[gene, patient_id] = save_copy_number

# 印出分隔列代表製作 cna 矩陣階段結束
#print("---------------------- Finish generating cna matrix -----------------------------")

# ================================= 儲存 sample 資料 ====================================
# # 創建一個字典，以'PATIENT_ID'為鍵，'DATE'為值
# date_sample_mapping_date = new_ngs_df.set_index('PATIENT_ID')['DATE'].to_dict()
# date_sample_mapping_path_number = new_ngs_df.set_index('PATIENT_ID')['Path_Number'].to_dict()
# date_sample_mapping_mp_number = new_ngs_df.set_index('PATIENT_ID')['MP_Number'].to_dict()
# date_sample_mapping_cell_block_number = new_ngs_df.set_index('PATIENT_ID')['Cell_Block_Number'].to_dict()
# # 使用 map 函數將對應的日期填入 'REPORT_DATE'
# data_sample_df['REPORT_DATE'] = data_sample_df['PATIENT_ID'].map(date_sample_mapping_date)
# data_sample_df['PATH_NUMBER'] = data_sample_df['PATIENT_ID'].map(date_sample_mapping_path_number)
# data_sample_df['MP_NUMBER'] = data_sample_df['PATIENT_ID'].map(date_sample_mapping_mp_number)
# data_sample_df['CELL_BLOCK_NUMBER'] = data_sample_df['PATIENT_ID'].map(date_sample_mapping_cell_block_number)
# # 確認 data_sample_df 的值都有正確
save_sample_data(data_sample_df, output_folder, output_file_data_clinical_sample)
#print(f"原始 data_sample_df 的形狀: {data_sample_df.shape}")
# ================================= 儲存 patient 資料 ===================================
save_patient_data(data_patient_df, output_folder, output_file_data_clinical_patient)
#print(f"原始 data_patient_df 的形狀: {data_patient_df.shape}")
# ================================= 儲存 mutations 資料 =================================
#print(f"原始 data_mutations_df 的形狀: {data_mutations_df.shape}")

# 挑出 "Hugo_Symbol" 欄位值為 'NA'、缺失值或是 None 的行
na_hugo_symbol_rows = data_mutations_df[data_mutations_df['Hugo_Symbol'].isin(['NA', np.nan, None])]
# 取得相應的 "Tumor_Sample_Barcode"
na_hugo_symbol_tumor_sample_barcodes = na_hugo_symbol_rows['Tumor_Sample_Barcode'].tolist()
# #print(f"這些是 data_mutations_df 中 'Hugo_Symbol' 欄位值為 'NA'、缺失值或是 None 的行，對應的 Tumor_Sample_Barcode: {na_hugo_symbol_tumor_sample_barcodes}")
#print(f"這些 'Hugo_Symbol' 欄位值為 'NA'、缺失值或是 None 的行總共有幾行: {len(na_hugo_symbol_tumor_sample_barcodes)}")
# 用集合的方式去除掉重複的值
na_hugo_symbol_tumor_sample_barcodes = list(set(na_hugo_symbol_tumor_sample_barcodes))
#print(f"這些是 data_mutations_df 中 'Hugo_Symbol' 欄位值為 'NA'、缺失值或是 None 的行，對應的 Tumor_Sample_Barcode(不重複): {na_hugo_symbol_tumor_sample_barcodes}")
#print(f"這些 'Hugo_Symbol' 欄位值為 'NA'、缺失值或是 None 的資料，總共有幾筆: {len(na_hugo_symbol_tumor_sample_barcodes)}")
# 刪除所有 "Hugo_Symbol" 欄位值為 'NA'、缺失值或是 None 的資料(以 Hugo_Symbol 為單位)
data_mutations_df = data_mutations_df[~data_mutations_df['Tumor_Sample_Barcode'].isin(na_hugo_symbol_tumor_sample_barcodes)]
#print(f"刪除有問題的 Hugo_Symbol 之後 data_mutations_df 的形狀: {data_mutations_df.shape}")

# # 檢查 data_mutations_df 裡面每一行是不是都有 18 個值
# missing_tumor_sample_barcodes = data_mutations_df[data_mutations_df.apply(lambda row: len(row) < 18, axis=1)]['Tumor_Sample_Barcode'].tolist()
# #print(f"這些是 data_mutations_df 中有問題的行，數值不足 18 個，與表格欄位不相符: {missing_tumor_sample_barcodes}")
# #print(f"這些是有問題，數值不足 18 個，的行，總共有幾行: {len(missing_tumor_sample_barcodes)}")
# # 用集合的方式去除掉重複的值
# missing_tumor_sample_barcodes = list(set(missing_tumor_sample_barcodes))
# #print(f"這些是 data_mutations_df 中有問題的資料，數值不足 18 個，與表格欄位不相符: {missing_tumor_sample_barcodes}")
# #print(f"這些是有問題，數值不足 18 個，的資料，總共有幾筆: {len(missing_tumor_sample_barcodes)}")
# # 刪除所有數值不足 18 個的行
# data_mutations_df = data_mutations_df[~data_mutations_df['Tumor_Sample_Barcode'].isin(missing_tumor_sample_barcodes)]
# #print(f"刪除有問題的 Tumor_Sample_Barcode 之後 data_mutations_df 的形狀: {data_mutations_df.shape}")

# # 確認 data_mutations_df 裡面的值是否正確，若有不正確，就整筆資料刪掉(相同的Tumor_Sample_Barcode)
# na_tumor_sample_barcodes = data_mutations_df[data_mutations_df.apply(lambda row: any(val in ['NA', np.nan, None] for val in row), axis=1)]['Tumor_Sample_Barcode'].tolist()
# # #print(f"這些是 data_mutations_df 中有問題的行，有值等於'NA'、缺失值或是 None，照理說基因資料每一個欄位都要有值: {na_tumor_sample_barcodes}")
# #print(f"這些是有問題的行，總共有幾行: {len(na_tumor_sample_barcodes)}")
# # 用集合的方式去除掉重複的值
# na_tumor_sample_barcodes = list(set(na_tumor_sample_barcodes))
# #print(f"這些是 data_mutations_df 中有問題的資料，有值等於'NA'，照理說基因資料每一個欄位都要有值: {na_tumor_sample_barcodes}")
# #print(f"這些是有問題的資料，總共有幾筆: {len(na_tumor_sample_barcodes)}")
# # 刪除所有含有 'NA'、缺失值或是 None 值的 tumor_sample_barcode 相同的行
# data_mutations_df = data_mutations_df[~data_mutations_df['Tumor_Sample_Barcode'].isin(na_tumor_sample_barcodes)]
# #print(f"刪除有問題的 tumor_sample_barcode 之後 data_mutations_df 的形狀: {data_mutations_df.shape}")

save_mutations_data(data_mutations_df, output_folder, output_file_data_mutations_extended)
# ================================= 儲存 cna 資料 =======================================
save_cna_data(matrix_cna, output_folder, output_file_data_cna)
# ================================= 儲存 sv 資料 ========================================
#print(f"原始 data_sv_df 的形狀: {data_sv_df.shape}")

# # 挑出 "Tumor_Seq_Allele2 " 欄位值為 'NA'、缺失值或是 None 的行
# na_tumor_seq_alle2_rows = data_sv_df[data_sv_df['Tumor_Seq_Allele2 '].isin(['NA', np.nan, None])]
# # 取得相應的 "Sample_Id"
# na_tumor_seq_alle2_rows = na_tumor_seq_alle2_rows['Sample_Id'].tolist()
# #print(f"這些是 data_sv_df 中 'Tumor_Seq_Allele2 ' 欄位值為 'NA'、缺失值或是 None 的行，對應的 Sample_Id: {na_tumor_seq_alle2_rows}")
# #print(f"這些 'Tumor_Seq_Allele2 ' 欄位值為 'NA'、缺失值或是 None 的行總共有幾行: {len(na_tumor_seq_alle2_rows)}")
# # 用集合的方式去除掉重複的值
# na_tumor_seq_alle2_rows = list(set(na_tumor_seq_alle2_rows))
# #print(f"這些是 data_sv_df 中 'Tumor_Seq_Allele2 ' 欄位值為 'NA'、缺失值或是 None 的行，對應的 Sample_Id(不重複): {na_tumor_seq_alle2_rows}")
# #print(f"這些 'Tumor_Seq_Allele2 ' 欄位值為 'NA'、缺失值或是 None 的資料，總共有幾筆: {len(na_tumor_seq_alle2_rows)}")
# # 刪除所有 "Tumor_Seq_Allele2 " 欄位值為 'NA'、缺失值或是 None 的資料(以 Tumor_Seq_Allele2  為單位)
# data_sv_df = data_sv_df[~data_sv_df['Sample_Id'].isin(na_tumor_seq_alle2_rows)]
# #print(f"刪除有問題的 Tumor_Seq_Allele2  之後 data_sv_df 的形狀: {data_sv_df.shape}")

# # 檢查 data_sv_df 裡面每一行是不是都有 8 個值
# missing_sample_ids = data_sv_df[data_sv_df.apply(lambda row: len(row) < 8, axis=1)]['Sample_Id'].tolist()
# #print(f"這些是 data_sv_df 中有問題的行，數值不足 8 個，與表格欄位不相符: {missing_sample_ids}")
# #print(f"這些是有問題的行，總共有幾行: {len(missing_sample_ids)}")
# # 用集合的方式去除掉重複的值
# missing_sample_ids = list(set(missing_sample_ids))
# #print(f"這些是 data_sv_df 中有問題的資料，數值不足 8 個，與表格欄位不相符: {missing_sample_ids}")
# #print(f"這些是有問題的資料，總共有幾筆: {len(missing_sample_ids)}")
# # 刪除所有數值不足 8 個的行
# data_sv_df = data_sv_df[~data_sv_df['Sample_Id'].isin(missing_sample_ids)]
# #print(f"刪除有問題的 sample_id 之後 data_sv_df 的形狀: {data_sv_df.shape}")

# 確認 data_sv_df 裡面的值是否正確，若有不正確，就整筆資料刪掉(相同的Sample_Id)
na_sample_id = data_sv_df[data_sv_df.apply(lambda row: any(val in ['NA', np.nan, None] for val in row), axis=1)]['Sample_Id'].tolist()
# #print(f"這些是 data_sv_df 有問題的行，有值等於'NA'、缺失值或是 None，照理說基因資料每一個欄位都要有值: {na_sample_id}")
#print(f"這些是有問題的行，總共有幾行: {len(na_sample_id)}")
# 用集合的方式去除掉重複的值
na_sample_id = list(set(na_sample_id))
#print(f"這些是 data_sv_df 有問題的資料，有值等於'NA'、缺失值或是 None，照理說基因資料每一個欄位都要有值: {na_sample_id}")
#print(f"這些是有問題資料，總共有幾筆: {len(na_sample_id)}")
# 刪除所有含有 'NA'、缺失值或是 None 值的 Sample_Id 相同的行
data_sv_df = data_sv_df[~data_sv_df['Sample_Id'].isin(na_sample_id)]
#print(f"刪除有問題的 sample_id 之後 data_sv_df 的形狀: {data_sv_df.shape}")

save_sv_data(data_sv_df, output_folder, output_file_data_sv)
#print("---------------------- Finish saving data files ---------------------------------")
# 印出此行代表所有上傳到 cbio portal 的資料前處理步驟結束
#print("################ Finish generating data and meta file for uploading to cbio portal ###############")
#print("\n")

#print("################ Start saving wrong mutations and sv data ################")
# 把錯誤的資料儲存下來，包含錯誤的 mutations 和 sv
# na_sample_id = add_parentheses(na_sample_id)
na_sample_id = [re.sub(r'_(\w+)', r'_(\1)', s) for s in na_sample_id]
# na_hugo_symbol_tumor_sample_barcodes = add_parentheses(na_hugo_symbol_tumor_sample_barcodes)
na_hugo_symbol_tumor_sample_barcodes = [re.sub(r'_(\w+)', r'_(\1)', s) for s in na_hugo_symbol_tumor_sample_barcodes]

# 建立字典
data_to_save = {
    "wrong_sv": na_sample_id,
    "wrong_mutations": na_hugo_symbol_tumor_sample_barcodes
}

# 獲取當前 UTC 時間
utc_now = datetime.utcnow()

# 設置時區為台灣
taipei_timezone = pytz.timezone('Asia/Taipei')
taipei_now = utc_now.replace(tzinfo=pytz.utc).astimezone(taipei_timezone)

# 將日期和時間轉換為字串格式
timestamp_str = taipei_now.strftime("%Y-%m-%d_%H%M%S")
# timestamp_str = taipei_now.strftime("%Y-%m-%d_%H%M%S%f") # 精確到毫秒

save_wrong_json_file = f'{save_wrong_json_directory}/{timestamp_str}_wrong.json'

# 寫入 JSON 檔案
with open(save_wrong_json_file, 'w', encoding='UTF-8') as json_file:
    json.dump(data_to_save, json_file, indent=4)

#print(f"已將 mutation 或是 sv 有誤的 json 報告儲存到: {save_wrong_json_file}")
#print("################ Finish saving wrong mutations and sv data ###############")
#print("\n")

# ======================================== 自動生成 case lists 檔 ========================================
#print("################ Start generating case lists for uploading to cbio portal ################")
# 從 data_sample_df 這個 dataframe 中提取 'SAMPLE_ID' 的值
case_list_sequenced_df = data_sample_df[['SAMPLE_ID']].copy()
# 取得 'SAMPLE_ID' 列的值
sample_ids = case_list_sequenced_df['SAMPLE_ID']
# 定義一個變數來儲存 cases sequenced
text_data_sequenced = '	'.join(str(sample_id) for sample_id in sample_ids)

# 從 data_cna_df 這個 dataframe 中提取 'PATIENT_ID' 的值
sample_ids = data_cna_df['PATIENT_ID'].copy().unique()
# 定義一個變數來儲存 cases cna
text_data_cna = '	'.join(str(sample_id) for sample_id in sample_ids)

# 從 data_sv_df 這個 dataframe 中提取 'Sample_Id' 的值
sample_ids = data_sv_df['Sample_Id'].copy()
# 定義一個變數來儲存 cases sv
text_data_sv = '	'.join(str(sample_id) for sample_id in sample_ids)

# 定義用來製做 cases sequenced 檔案的變數
stable_id_se = cancer_study_identifier + "_sequenced"
case_list_ids_se = text_data_sequenced
#print(f'Length of case list sequenced: {len(case_list_ids_se)}')
# 定義用來製做 cases cna 檔案的變數
stable_id_cna = cancer_study_identifier + "_cna"
case_list_ids_cna = text_data_cna
#print(f'Length of case list cna: {len(case_list_ids_cna)}')
# 定義用來製做 cases sv 檔案的變數
stable_id_sv = cancer_study_identifier + "_sv"
case_list_ids_sv = text_data_sv
#print(f'Length of case list sv: {len(case_list_ids_sv)}')

if len(case_list_ids_se) != 0:
    # 拼接完整的路徑名稱(cases_sequenced)
    output_path = os.path.join(output_folder_case, output_file_cases_sequenced)
    # 寫到 case sequenced
    with open(output_path, "w", encoding="utf-8") as file:
        file.write(f"cancer_study_identifier: {cancer_study_identifier}\n")
        file.write(f"stable_id: {stable_id_se}\n")
        file.write(f"case_list_name: {case_list_name_se}\n")
        file.write(f"case_list_description: {case_list_description_se}\n")
        file.write(f"case_list_category: {case_list_category_se}\n")
        file.write(f"case_list_ids:	{case_list_ids_se}\n")
    #print(f"cases sequenced is saved to '{output_path}'")

if len(case_list_ids_cna) != 0:
    # 拼接完整的路徑名稱(cases_cna)
    output_path = os.path.join(output_folder_case, output_file_cases_cna)
    # 寫到 case cna
    with open(output_path, "w", encoding="utf-8") as file:
        file.write(f"cancer_study_identifier: {cancer_study_identifier}\n")
        file.write(f"stable_id: {stable_id_cna}\n")
        file.write(f"case_list_name: {case_list_name_cna}\n")
        file.write(f"case_list_description: {case_list_description_cna}\n")
        file.write(f"case_list_category: {case_list_category_cna}\n")
        file.write(f"case_list_ids:	{case_list_ids_cna}\n")
    #print(f"cases cna is saved to '{output_path}'")

if len(case_list_ids_sv) != 0:
    # 拼接完整的路徑名稱(cases_sv)
    output_path = os.path.join(output_folder_case, output_file_cases_sv)
    # 寫到 case sv
    with open(output_path, "w", encoding="utf-8") as file:
        file.write(f"cancer_study_identifier: {cancer_study_identifier}\n")
        file.write(f"stable_id: {stable_id_sv}\n")
        file.write(f"case_list_name: {case_list_name_sv}\n")
        file.write(f"case_list_description: {case_list_description_sv}\n")
        file.write(f"case_list_category: {case_list_category_sv}\n")
        file.write(f"case_list_ids:	{case_list_ids_sv}\n")
    #print(f"cases sv is saved to '{output_path}'")
#print("################ Finish generating case lists for uploading to cbio portal ###############")
#print("\n")

print("##################################################################################################")
print("################ Finish .txt files generating process for uploading to cbioportal ################")
print("##################################################################################################")
print("\n")
