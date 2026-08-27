import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import openpyxl
import pickle
import os
import gc
from sklearn.linear_model import ElasticNet
from sklearn.cross_decomposition import PLSRegression
from sklearn.linear_model import BayesianRidge
from sklearn.linear_model import LinearRegression
from sklearn.neighbors import KNeighborsRegressor
from sklearn.svm import SVR
from sklearn.linear_model import Ridge
from sklearn.linear_model import Lasso
from sklearn.neural_network import MLPRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.tree import ExtraTreeRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import AdaBoostRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import BaggingRegressor
from sklearn.model_selection import train_test_split
from scipy.stats import pearsonr
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import r2_score, mean_squared_error
from skimage.metrics import normalized_root_mse as compare_nrmse
from scipy.stats import spearmanr
import csv
import os
import sys
import json
import time

def validate_sequence(sequence):
    """验证和清理蛋白质序列"""
    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
    cleaned = ''.join(aa for aa in sequence.upper() if aa in valid_aa)
    return cleaned

def main(input_sequence=None):
    print("创建必要的输出目录...")
    os.makedirs('./output_encoding', exist_ok=True)
    os.makedirs('./predict_data', exist_ok=True)
    os.makedirs('./predict_encoding', exist_ok=True)
    os.makedirs('./data', exist_ok=True)
    print("目录创建完成")

    # 如果有传入序列，使用传入的序列
    if input_sequence:
        print(f"原始传入序列: {input_sequence}")
        
        # 检测是否是未替换的n8n变量
        if "{{" in input_sequence and "}}" in input_sequence:
            print("检测到未替换的n8n变量，使用默认序列")
            # 使用默认序列
            default_sequence = "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
            cleaned_sequence = default_sequence
        else:
            # 清理和验证序列
            cleaned_sequence = validate_sequence(input_sequence)
        
        print(f"最终使用序列: {cleaned_sequence}")
        
        with open('./data/wt.txt', 'w') as f:
            f.write(cleaned_sequence)

if __name__ == "__main__":
    start_time = time.time()
    print(f"脚本开始执行: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    
    # 支持命令行参数
    if len(sys.argv) > 1:
        # 直接获取序列参数
        sequence = sys.argv[1]
        main(sequence)
    else:
        main()

    # 检查必要文件
    required_files = [
        './data/protein_encoding.xlsx',
        './data/aaindex id.csv', 
        './data/input_data.csv',
        './data/wt.txt'
    ]
    
    missing_files = []
    for file_path in required_files:
        if not os.path.exists(file_path):
            missing_files.append(file_path)
            print(f"错误: 找不到必要文件 {file_path}")
    
    if missing_files:
        print(f"缺少文件: {missing_files}")
        error_result = {
            "status": "error",
            "error": f"缺少必要文件: {', '.join(missing_files)}",
            "missing_files": missing_files
        }
        print(json.dumps(error_result))
        exit(1)

    print("所有必要文件都存在")
    gc.collect()

    try:
        # 读取Excel文件数据
        wb = openpyxl.load_workbook('./data/protein_encoding.xlsx')
        sheet = wb['蛋白编码']

        cells = sheet['B2':'V21']
        all_aaindex = []

        for r in cells:
            index_list = []
            for c in r:
                index_list.append(c.value)
            aa_list = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','X']
            aaindex = dict(zip(aa_list, index_list))
            all_aaindex.append(aaindex)

        N = 1024

        # 读取aaindex id
        file3 = open('./data/aaindex id.csv', 'r')
        aaindex_id = []
        line3 = file3.readline()
        aaindex_id.append(line3.strip())

        while line3 != '':
            line3 = file3.readline()
            aaindex_id.append(line3.strip())
        aaindex_id.pop()
        print("aaindex_id:", aaindex_id)

        # 读取输入数据
        data = pd.read_csv('./data/input_data.csv', sep=',')
        all_seqs = data['Sequence']

        for aaindex,encoding_name in zip(all_aaindex,aaindex_id):
            print("current aaindex is: ", aaindex)
            print("current encoding name is: ", encoding_name)
            print("*"*30)
            ls=[]
            for k,v in all_seqs.items():
                numeric_list = []
                seq = v

                for aa in seq:
                    index = aaindex[aa]
                    numeric_list.append(index)
                average = sum(numeric_list)/len(numeric_list)
                numeric_list[:] = [i - average for i in numeric_list]

                zero_padding_list = [0] * (N-len(numeric_list))
                numeric_list.extend(zero_padding_list)

                x = np.arange(len(numeric_list))
                half_x = x[range(int(N / 2))]
                fft = np.fft.fft(numeric_list)
                abs_fft = np.abs(fft)
                abs_fft = abs_fft/1024
                half_y = abs_fft[range(int(N/2))]
                ls.append(half_y)
            
            # 确保A在循环外部被定义
            if ls:  # 如果ls不为空
                B = np.array(ls)
                A = pd.DataFrame(B)
            else:  # 如果ls为空，创建空的DataFrame
                A = pd.DataFrame()
                print("警告: ls列表为空，创建空的DataFrame")
            
            print(A.shape)
            A.to_csv(f'./output_encoding/{encoding_name}_result_of_encoding.csv')
            gc.collect()

        filter = [".csv"]
        FitnessBest = []
        x_type = []

        def all_path(dirname):
          for maindir, subdir, file_name_list in os.walk(dirname):
            for filename in file_name_list:
              x_type.append(filename)
              apath = os.path.join(maindir, filename)
              FitnessBest.append(apath)
          return FitnessBest

        print(all_path(r"output_encoding"))

        R2_list = []
        p=0

        for i in FitnessBest:
          filepath = i

          # load data
          data_x = pd.read_csv(filepath)
          data = pd.read_csv('./data/input_data.csv')
          data_y = data['dep']

          # split data
          x_train, x_test, y_train, y_test = train_test_split( data_x, data_y, test_size=0.20, random_state=15)
          gc.collect()
          # model fit ： You can choose different algorithms for the model.
          model = GradientBoostingRegressor()
          model.fit(x_train, y_train)
          y_pred = model.predict(x_test)

          # calculate R2, rmse, mae, pccs
          R2 = round(r2_score(y_test, y_pred),4)
          mse = mean_squared_error(y_test, y_pred)
          rmse = np.sqrt(mse)
          mae = mean_absolute_error(y_test, y_pred)
          y_test2 = np.array(y_test.squeeze(axis=None))
          y_test3 = y_test2.astype(float)
          nrmse = compare_nrmse(y_test3, y_pred)
          pccs = pearsonr(y_test3,y_pred)
          correlation, p_value = spearmanr(y_test,y_pred)
          R2_list.append(R2)

          print(p)
          print('Test RMSE: %.4f , NRMSE: %.4f , R2: %.4f ,R:%s, P:%.4f' % (rmse, nrmse, R2, pccs, correlation))
          p=p+1

        print('*' * 30)
        print('R2:',max(R2_list))
        print('encoding_now',x_type[R2_list.index(max(R2_list))])

        filepath = os.path.join("./output_encoding", x_type[R2_list.index(max(R2_list))])

        # load data
        data_x = pd.read_csv(filepath)
        data = pd.read_csv('./data/input_data.csv')
        data_y = data['dep']

        # split data
        x_train, x_test, y_train, y_test = train_test_split(data_x, data_y, test_size=0.20, random_state=15)

        model_name = 'your_model'
        file_model = model_name + '.pickle'

        model = GradientBoostingRegressor()
        model.fit(x_train, y_train)

        with open(file_model, 'wb') as fw:
            pickle.dump(model, fw)

        y_pred = model.predict(x_test)
        gc.collect()
        print(y_test)
        print(y_pred)

        R2 = round(r2_score(y_test, y_pred),4)
        mse = mean_squared_error(y_test, y_pred)
        rmse = np.sqrt(mse)
        mae = mean_absolute_error(y_test, y_pred)
        y_test2 = np.array(y_test.squeeze(axis=None))
        y_test3 = y_test2.astype(float)
        nrmse = compare_nrmse(y_test3, y_pred)
        pccs = pearsonr(y_test3,y_pred)
        correlation, p_value = spearmanr(y_test,y_pred)

        print('Test RMSE: %.4f , NRMSE: %.4f , R2: %.4f ,R:%s, P:%.4f' % (rmse, nrmse, R2, pccs, correlation))
        print('*' * 30)

        with open('./best_encoding_result.txt', 'w') as f:
            max_r2 = max(R2_list)
            best_encoding = x_type[R2_list.index(max_r2)]
            encoding_name = best_encoding.split('_')[0]
            f.write(f'{encoding_name}')

        # 生成突变序列
        with open("./data/wt.txt", "r") as f:
            wt_seq = f.read().strip()

        seq = list(wt_seq)
        amino_acids = ['G', 'A', 'V', 'L', 'I', 'P', 'F', 'Y', 'W', 'S',
                       'T', 'C', 'M', 'N', 'Q', 'D', 'E', 'K', 'R', 'H']

        # 生成突变体文件（确保无空行）
        with open("mutations.txt", 'w') as fw:
            for i in range(len(seq)):
                original = seq[i]
                for aa in amino_acids:
                    if aa == original: continue
                    seq[i] = aa
                    fw.write(f"{''.join(seq)}\n")
                seq[i] = original

        # 转换CSV（过滤空行）
        with open("mutations.txt", 'r') as txt_file, \
                open("your_data.csv", 'w', newline='') as csv_file:

            writer = csv.writer(csv_file)
            writer.writerow(["ID", "Sequence"])

            line_counter = 0
            for idx, line in enumerate(txt_file, start=1):
                clean_line = line.strip()
                if clean_line:
                    line_counter += 1
                    writer.writerow([line_counter, clean_line])

        txt_path = "best_encoding_result.txt"
        csv_dir = "./predict_data"
        csv_name = "aaindex id.csv"
        csv_path = os.path.join(csv_dir, csv_name)
        with open(txt_path, "r", encoding="utf-8") as f:
            txt_content = f.read().strip()

        # 3. 将内容写入CSV的A1单元格
        with open(csv_path, "w", encoding="utf-8", newline="") as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow([txt_content])

        wb  =openpyxl.load_workbook('./predict_data/protein_encoding.xlsx')
        sheet=wb['蛋白编码']

        cells=sheet['B2':'V567']
        all_aaindex=[]

        for r in cells:
            index_list = []
            for c in r:
                index_list.append(c.value)
            aa_list=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','X']
            aaindex=dict(zip(aa_list,index_list))
            all_aaindex.append(aaindex)

        N = 1024

        file3=open('./predict_data/aaindex id.csv', 'r')
        aaindex_id=[]
        line3=file3.readline()
        aaindex_id.append(line3.strip())

        while line3 !='':
            line3=file3.readline()
            aaindex_id.append(line3.strip())
        aaindex_id.pop()
        print(aaindex_id)

        data = pd.read_csv('your_data.csv', sep=',')
        all_seqs = data['Sequence']

        for aaindex,predict_encoding_name in zip(all_aaindex,aaindex_id):
            print("current aaindex is: ", aaindex)
            print("current encoding name is: ", predict_encoding_name)
            print("*"*30)
            ls=[]
            for k,v in all_seqs.items():
                numeric_list = []
                seq = v

                for aa in seq:
                    index = aaindex[aa]
                    numeric_list.append(index)
                average = sum(numeric_list)/len(numeric_list)
                numeric_list[:] = [i - average for i in numeric_list]

                zero_padding_list = [0] * (N-len(numeric_list))
                numeric_list.extend(zero_padding_list)

                x = np.arange(len(numeric_list))
                half_x = x[range(int(N / 2))]
                fft = np.fft.fft(numeric_list)
                abs_fft = np.abs(fft)
                abs_fft = abs_fft/1024
                half_y = abs_fft[range(int(N/2))]
                ls.append(half_y)
                B = np.array(ls)
                A = pd.DataFrame(B)
            print(A.shape)
            A.to_csv(f'./predict_encoding/{predict_encoding_name}_result_of_encoding.csv')
            gc.collect()

        # 待预测的酶序列编码后csv文件 —— 显式使用训练阶段选出的最优编码，避免依赖循环尾值的隐含假设
        best_encoding_name = ''
        if os.path.exists('best_encoding_result.txt'):
            with open('best_encoding_result.txt', 'r', encoding='utf-8') as f:
                best_encoding_name = f.read().strip()
        if not best_encoding_name:
            best_encoding_name = predict_encoding_name  # 回退：未记录最优编码时使用循环最后一个编码
        predict_filepath = f'./predict_encoding/{best_encoding_name}_result_of_encoding.csv'
        if not os.path.exists(predict_filepath):
            # 回退：最优编码文件不存在时，尝试循环最后一个编码
            predict_filepath = f'./predict_encoding/{predict_encoding_name}_result_of_encoding.csv'
        data_x = pd.read_csv(predict_filepath)

        #读取编码前的酶序列
        data = pd.read_csv(r'.\your_data.csv', sep=',')
        data_sequnece = data['Sequence']

        model_name = 'your_model'
        file_model = model_name + '.pickle'
        outputpath = 'result.csv'

        # 加载 model
        with open(file_model, "rb") as f:
            model = pickle.load(f)

        gc.collect()
        #预测
        y_pred = model.predict(data_x)
        print(y_pred)
        gc.collect()

        RESULT = pd.DataFrame()
        df_y_test = data_sequnece.to_frame(name='Sequnece')
        RESULT['Sequence'] = df_y_test.reset_index(drop=True)
        df_y_pred = pd.DataFrame(data=y_pred[0:], columns=['yred'])
        RESULT['Y_PRED'] = df_y_pred
        RESULT.to_csv(outputpath, sep=',', index=False, header=True)

        result_df = pd.read_csv('result.csv')
        top_10_sequences = result_df.sort_values(by='Y_PRED', ascending=False).head(10)
        top_10_sequences.to_csv('top_10_sequences.csv', index=False)
        gc.collect()

        # 读取现有结果文件
        with open('best_encoding_result.txt', 'r') as f:
            best_encoding = f.read().strip()
        
        # 读取top序列结果
        top_df = pd.read_csv('top_10_sequences.csv')
        
        # 读取野生型序列
        with open("./data/wt.txt", "r") as f:
            wt_sequence = f.read().strip()

        # 构建突变体分析结果
        top_mutants = []
        for i, row in top_df.head(10).iterrows():
            mutant_seq = row['Sequence']
            score = row['Y_PRED']
            
            # 找出突变位置
            mutations = []
            for pos, (wt_aa, mut_aa) in enumerate(zip(wt_sequence, mutant_seq)):
                if wt_aa != mut_aa:
                    mutations.append(f"{wt_aa}{pos+1}{mut_aa}")
            
            top_mutants.append({
                "rank": i + 1,
                "sequence": mutant_seq,
                "score": float(score),
                "mutations": mutations,
                "mutation_count": len(mutations)
            })

        # 构建最终结果
        result = {
            "status": "success",
            "best_encoding": best_encoding,
            "top_sequences": "突变体分析完成",
            "message": "蛋白质突变体优化预测完成",
            
            "mutant_analysis": {
                "wild_type_sequence": wt_sequence,
                "wild_type_length": len(wt_sequence),
                "total_mutants_analyzed": len(top_df),
                "top_mutants": top_mutants,
                "best_mutant_score": top_mutants[0]["score"] if top_mutants else 0,
                "best_mutant_mutations": top_mutants[0]["mutations"] if top_mutants else []
            }
        }
        
        print(json.dumps(result))
        
    except (KeyError, ValueError, FileNotFoundError, OSError) as e:
        import traceback
        error_result = {
            "status": "error",
            "error": str(e),
            "error_type": type(e).__name__,
            "traceback": traceback.format_exc()
        }
        print(json.dumps(error_result, ensure_ascii=False))
    except Exception as e:
        import traceback
        error_result = {
            "status": "error",
            "error": f"未预期的错误: {str(e)}",
            "error_type": type(e).__name__,
            "traceback": traceback.format_exc()
        }
        print(json.dumps(error_result, ensure_ascii=False))
    
    execution_time = time.time() - start_time
    print(f"脚本执行完成，总耗时: {execution_time:.2f}秒")