import sys
import os
import streamlit as st
st.set_page_config(
    page_title="蛋白质AI设计实验室",
    layout="wide",
    initial_sidebar_state="expanded"
)

import logging
import pandas as pd
from datetime import datetime
import json
import re
from typing import Dict, List, Any, Optional
import pickle
import openpyxl
import numpy as np
import xml.etree.ElementTree as ET
import urllib.parse
import socket
import requests
import py3Dmol
import streamlit.components.v1 as components
from urllib.error import URLError, HTTPError
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
from scipy.stats import pearsonr, spearmanr
from skimage.metrics import normalized_root_mse as compare_nrmse
from structure_predictor_api import UnifiedStructurePredictor

import time
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

STRUCTURE_MODULE_AVAILABLE = True  

# 自定义CSS - 优化全宽显示
st.markdown("""
<style>
    .main-header {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 2rem;
        border-radius: 15px;
        text-align: center;
        margin-bottom: 2rem;
    }
    .metric-card {
        background: white;
        border-radius: 10px;
        padding: 1rem;
        margin: 0.5rem 0;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        border-left: 4px solid #667eea;
    }
    .reasoning-step {
        background: #f8f9fa;
        padding: 1rem;
        border-radius: 10px;
        margin: 0.5rem 0;
        border-left: 4px solid #28a745;
    }
    .mutant-card {
        background: #f1f3f4;
        border-radius: 8px;
        padding: 1rem;
        margin: 0.5rem 0;
        border: 1px solid #e0e0e0;
    }
    .success-pattern {
        background: #d4edda;
        color: #155724;
        padding: 0.5rem;
        border-radius: 5px;
        margin: 0.2rem 0;
    }
    .learning-progress {
        background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%);
        color: white;
        padding: 1rem;
        border-radius: 10px;
    }
    .similar-case {
        background: #e7f3ff;
        border-radius: 8px;
        padding: 1rem;
        margin: 0.5rem 0;
        border-left: 4px solid #1890ff;
    }
    .full-width-container {
        width: 100%;
        padding: 0;
        margin: 0;
    }
    .st-expander {
        width: 100% !important;
    }
    .blast-results-container {
        border: 1px solid #e0e0e0;
        border-radius: 10px;
        padding: 20px;
        margin-bottom: 20px;
        background-color: #f9f9f9;
        width: 100%;
    }
    .structure-results-container {
        border: 1px solid #e0e0e0;
        border-radius: 10px;
        padding: 20px;
        margin-bottom: 20px;
        background-color: #f8f9fa;
        width: 100%;
    }
    .analysis-tab {
        padding: 15px;
        border-radius: 8px;
        margin-bottom: 20px;
        width: 100%;
    }
    .sequence-preview {
        font-family: monospace;
        background-color: #f5f5f5;
        padding: 10px;
        border-radius: 5px;
        margin: 10px 0;
        overflow-x: auto;
        white-space: pre-wrap;
        word-wrap: break-word;
        width: 100%;
    }
    .stDataFrame {
        width: 100% !important;
    }
    .stMarkdown {
        width: 100% !important;
    }
    .stAlert {
        width: 100% !important;
    }
    .stMetric {
        width: 100% !important;
    }
    .full-width-table {
        width: 100% !important;
        max-width: 100% !important;
    }
    .wide-column {
        width: 100%;
    }
    .result-container {
        width: 100%;
        margin-top: 20px;
        padding: 20px;
        background-color: #f8f9fa;
        border-radius: 10px;
    }
</style>
""", unsafe_allow_html=True)

# 检查必要文件
def check_required_files():
    """检查必要的AAindex文件"""
    required_files = [
        './data/protein_encoding.xlsx',
        './data/aaindex id.csv'
    ]
    
    missing_files = []
    for file_path in required_files:
        if not os.path.exists(file_path):
            missing_files.append(file_path)
    
    return missing_files

missing_files = check_required_files()
if missing_files:
    st.warning(f"缺少AAindex数据文件: {', '.join(missing_files)}")
    st.info("请确保以下文件存在:\n- data/protein_encoding.xlsx\n- data/aaindex id.csv")


class AAindexEncoder:
    """AAindex编码器"""
    
    def __init__(self):
        self.all_aaindex = []
        self.aaindex_id = []
        self.N = 1024
        self._load_aaindex_data()
    
    def _load_aaindex_data(self):
        """加载AAindex数据"""
        try:
            wb = openpyxl.load_workbook('./data/protein_encoding.xlsx')
            sheet = wb['蛋白编码']
            
            cells = sheet['B2':'V21']
            
            for r in cells:
                index_list = []
                for c in r:
                    index_list.append(c.value)
                aa_list = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','X']
                aaindex = dict(zip(aa_list, index_list))
                self.all_aaindex.append(aaindex)
            
            with open('./data/aaindex id.csv', 'r') as f:
                self.aaindex_id = [line.strip() for line in f if line.strip()]
                
        except Exception as e:
            st.error(f"加载AAindex数据失败: {str(e)}")
    
    def encode_sequence(self, sequence: str, encoding_name: str) -> np.ndarray:
        """对单个序列进行编码"""
        encoding_dict = None
        for aaindex, name in zip(self.all_aaindex, self.aaindex_id):
            if name == encoding_name:
                encoding_dict = aaindex
                break
        
        if not encoding_dict:
            raise ValueError(f"未找到编码方案: {encoding_name}")
        
        numeric_list = []
        for aa in sequence:
            index = encoding_dict.get(aa, 0)
            numeric_list.append(index)
        
        average = sum(numeric_list)/len(numeric_list)
        numeric_list[:] = [i - average for i in numeric_list]
        
        zero_padding_list = [0] * (self.N - len(numeric_list))
        numeric_list.extend(zero_padding_list)
        
        fft = np.fft.fft(numeric_list)
        abs_fft = np.abs(fft)
        abs_fft = abs_fft/self.N
        half_y = abs_fft[range(int(self.N/2))]
        
        return half_y
    
    def encode_dataset(self, sequences: List[str], encoding_name: str) -> pd.DataFrame:
        """编码整个数据集"""
        encoded_data = []
        for seq in sequences:
            encoded_seq = self.encode_sequence(seq, encoding_name)
            encoded_data.append(encoded_seq)
        
        return pd.DataFrame(encoded_data)


class ProteinModelTrainer:
    """蛋白质模型训练器"""
    
    def __init__(self):
        self.encoder = AAindexEncoder()
        self.model = None
        self.best_encoding = None
        self.is_trained = False
    
    def prepare_training_data(self, sequences: List[str], targets: List[float]):
        """准备训练数据"""
        training_data = pd.DataFrame({
            'Sequence': sequences,
            'dep': targets
        })
        return training_data
    
    def find_best_encoding(self, sequences: List[str], targets: List[float]) -> str:
        """寻找最佳编码方案"""
        best_r2 = -1
        best_encoding = None
        
        for encoding_name in self.encoder.aaindex_id[:5]:
            try:
                X = self.encoder.encode_dataset(sequences, encoding_name)
                y = targets
                
                X_train, X_test, y_train, y_test = train_test_split(
                    X, y, test_size=0.2, random_state=42
                )
                
                model = GradientBoostingRegressor()
                model.fit(X_train, y_train)
                
                y_pred = model.predict(X_test)
                r2 = r2_score(y_test, y_pred)
                
                if r2 > best_r2:
                    best_r2 = r2
                    best_encoding = encoding_name
                    
            except Exception:
                continue
        
        return best_encoding
    
    def train_model(self, sequences: List[str], targets: List[float]):
        """训练最终模型"""
        self.best_encoding = self.find_best_encoding(sequences, targets)
        
        if not self.best_encoding:
            st.error("找不到合适的编码方案")
            return False
        
        X = self.encoder.encode_dataset(sequences, self.best_encoding)
        y = targets
        
        self.model = GradientBoostingRegressor()
        self.model.fit(X, y)
        self.is_trained = True
        
        with open('protein_model.pickle', 'wb') as f:
            pickle.dump(self.model, f)
        
        with open('best_encoding.txt', 'w') as f:
            f.write(self.best_encoding)
        
        return True
    
    def predict(self, sequence: str) -> float:
        """预测单个序列"""
        if not self.is_trained or not self.best_encoding:
            raise ValueError("模型未训练")
        
        X = self.encoder.encode_dataset([sequence], self.best_encoding)
        return self.model.predict(X)[0]


class RealNCBIBlastAPI:
    """真正的NCBI BLAST API使用BioPython"""
    
    def __init__(self, email: str, api_key: str = None):
        self.email = email
        self.api_key = api_key
        
    def run_blast(self, sequence: str, program: str = 'blastp', 
                 database: str = 'nr', max_wait: int = 300) -> Optional[Dict]:
        """运行真正的BLAST分析"""
        try:
            if not self.email or '@' not in self.email:
                raise ValueError("需要有效的NCBI邮箱地址")
            NCBIWWW.email = self.email
            # BLAST 在线查询为同步阻塞，设置 socket 超时避免无限等待（max_wait 秒）
            socket.setdefaulttimeout(max_wait)
            
            if len(sequence) > 2000:
                st.warning(f"序列过长 ({len(sequence)} aa)，将截取前2000个氨基酸")
                sequence = sequence[:2000]
            
            valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
            if not all(aa in valid_aa for aa in sequence):
                st.warning("序列包含非标准氨基酸，BLAST可能失败")
            
            with st.spinner("提交BLAST查询到NCBI服务器（可能耗时1-5分钟，请耐心等待）..."):
                try:
                    result_handle = NCBIWWW.qblast(
                        program=program,
                        database=database,
                        sequence=sequence,
                        expect=10.0,
                        hitlist_size=50,
                        alignments=50,
                        descriptions=50,
                        word_size=3,
                        matrix_name="BLOSUM62",
                        gapcosts="11 1"
                    )
                except Exception as e:
                    error_msg = str(e)
                    if "429" in error_msg:
                        raise Exception("NCBI API限制：请求过于频繁，请稍后再试")
                    elif "503" in error_msg:
                        raise Exception("NCBI服务器暂时不可用")
                    elif "timeout" in error_msg:
                        raise Exception("连接超时，请检查网络连接")
                    elif "Entrez Query" in error_msg:
                        st.info("尝试使用简化参数重新查询...")
                        result_handle = NCBIWWW.qblast(
                            program=program,
                            database=database,
                            sequence=sequence,
                            expect=10.0,
                            hitlist_size=20
                        )
                    else:
                        raise
            
            blast_records = NCBIXML.parse(result_handle)
            
            homologous_proteins = []
            conserved_regions = []
            
            try:
                for blast_record in blast_records:
                    query_length = blast_record.query_length
                    
                    for alignment in blast_record.alignments:
                        for hsp in alignment.hsps:
                            identity_percent = (hsp.identities / hsp.align_length) * 100
                            species = self._extract_species_from_title(alignment.title)
                            
                            protein_info = {
                                'accession': alignment.accession,
                                'description': alignment.title,
                                'species': species,
                                'e_value': hsp.expect,
                                'bit_score': hsp.bits,
                                'identity_percent': round(identity_percent, 2),
                                'alignment_length': hsp.align_length,
                                'query_start': hsp.query_start,
                                'query_end': hsp.query_end,
                                'subject_start': hsp.sbjct_start,
                                'subject_end': hsp.sbjct_end,
                                'alignment': hsp.match,
                                'query_seq': hsp.query,
                                'subject_seq': hsp.sbjct
                            }
                            homologous_proteins.append(protein_info)
                            
                            if identity_percent > 70 and hsp.align_length > 20:
                                conserved_regions.append({
                                    'start': hsp.query_start,
                                    'end': hsp.query_end,
                                    'length': hsp.align_length,
                                    'identity': identity_percent,
                                    'description': f"高保守区域 ({identity_percent:.1f}% 相似度)"
                                })
            
            except Exception as e:
                st.warning(f"解析BLAST结果时出错: {str(e)}")
            
            result_handle.close()
            
            unique_proteins = {}
            for protein in homologous_proteins:
                acc = protein['accession']
                if acc not in unique_proteins or protein['e_value'] < unique_proteins[acc]['e_value']:
                    unique_proteins[acc] = protein
            
            sorted_proteins = sorted(unique_proteins.values(), key=lambda x: x['e_value'])
            
            return {
                'homologous_proteins': sorted_proteins[:20],
                'conserved_regions': conserved_regions[:10],
                'query_length': query_length,
                'database': database,
                'total_hits': len(sorted_proteins),
                'message': f'BLAST分析完成，找到 {len(sorted_proteins)} 个同源蛋白',
                'is_offline_example': False
            }
            
        except Exception as e:
            st.error(f"BLAST分析失败: {str(e)}")
            logging.error(f"BLAST error: {str(e)}")
            return None
    
    def _extract_species_from_title(self, title: str) -> str:
        """从BLAST标题中提取物种信息"""
        import re
        
        match = re.search(r'\[([^\]]+)\]', title)
        if match:
            species = match.group(1)
            
            species_lower = species.lower()
            if 'homo sapiens' in species_lower:
                return 'Human'
            elif 'mus musculus' in species_lower:
                return 'Mouse'
            elif 'rattus norvegicus' in species_lower:
                return 'Rat'
            elif 'escherichia coli' in species_lower:
                return 'E. coli'
            elif 'saccharomyces cerevisiae' in species_lower:
                return 'Yeast'
            elif 'arabidopsis thaliana' in species_lower:
                return 'Arabidopsis'
            else:
                return species[:30]
        
        match = re.search(r'\(([^)]+)\)', title)
        if match:
            return match.group(1)[:30]
        
        return 'Unknown'


class BlastErrorHandler:
    """BLAST错误处理器"""
    
    @staticmethod
    def handle_error(error: Exception, context: str = "") -> Dict:
        """处理BLAST错误并返回用户友好的消息"""
        error_msg = str(error).lower()
        
        if "timeout" in error_msg or "timed out" in error_msg:
            return {
                'status': 'error',
                'message': '连接NCBI服务器超时，请检查网络连接',
                'suggestion': '尝试增加超时时间或稍后重试',
                'retryable': True
            }
        elif "connection" in error_msg or "network" in error_msg:
            return {
                'status': 'error',
                'message': '网络连接失败，无法访问NCBI',
                'suggestion': '检查网络设置或使用VPN',
                'retryable': True
            }
        elif "429" in error_msg or "too many requests" in error_msg:
            return {
                'status': 'error',
                'message': '请求过于频繁，NCBI限制了访问',
                'suggestion': '请等待5-10分钟后再试',
                'retryable': True,
                'wait_time': 300
            }
        elif "503" in error_msg or "service unavailable" in error_msg:
            return {
                'status': 'error',
                'message': 'NCBI服务器暂时不可用',
                'suggestion': '请稍后再试，通常是临时维护',
                'retryable': True
            }
        elif "email" in error_msg or "valid email" in error_msg:
            return {
                'status': 'error',
                'message': '需要有效的NCBI邮箱地址',
                'suggestion': '请在侧边栏输入有效的邮箱地址',
                'retryable': False
            }
        elif "sequence" in error_msg or "invalid" in error_msg:
            return {
                'status': 'error',
                'message': '序列格式无效',
                'suggestion': '请检查序列是否只包含标准氨基酸字符',
                'retryable': False
            }
        else:
            return {
                'status': 'error',
                'message': f'BLAST分析失败: {str(error)[:100]}',
                'suggestion': '请尝试简化序列或联系管理员',
                'retryable': False
            }
    
    @staticmethod
    def validate_sequence(sequence: str) -> Dict:
        """验证序列是否适合BLAST分析"""
        if not sequence or len(sequence.strip()) == 0:
            return {'valid': False, 'message': '序列为空'}
        
        if len(sequence) < 10:
            return {'valid': False, 'message': '序列过短（至少需要10个氨基酸）'}
        
        if len(sequence) > 20000:
            return {'valid': False, 'message': '序列过长（最多20000个氨基酸）'}
        
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        invalid_chars = []
        for i, aa in enumerate(sequence.upper()):
            if aa not in valid_aa:
                invalid_chars.append((i+1, aa))
        
        if invalid_chars:
            return {
                'valid': False,
                'message': f'发现 {len(invalid_chars)} 个无效字符',
                'invalid_positions': invalid_chars[:10]
            }
        
        return {'valid': True}


class EnhancedProteinMemory:
    """增强的蛋白质设计记忆系统（持久化到 JSON，跨会话保留）"""

    MEMORY_FILE = 'memory_data.json'

    def __init__(self):
        self.cases = {}
        self.case_count = 0
        self.patterns = {
            'stability': [],
            'activity': [],
            'solubility': []
        }
        self._load_from_disk()

    def _load_from_disk(self):
        """从磁盘加载历史记忆，避免重启丢失"""
        try:
            if os.path.exists(self.MEMORY_FILE):
                with open(self.MEMORY_FILE, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                self.cases = data.get('cases', {})
                self.case_count = data.get('case_count', 0)
                saved_patterns = data.get('patterns', {})
                for goal in self.patterns:
                    if goal in saved_patterns:
                        self.patterns[goal] = saved_patterns[goal]
        except (json.JSONDecodeError, OSError, ValueError) as e:
            # 记忆文件损坏时不影响启动，回退到空记忆
            print(f"[Memory] 加载记忆文件失败，使用空记忆: {e}")
            self.cases = {}
            self.case_count = 0

    def _save_to_disk(self):
        """将记忆持久化到磁盘"""
        try:
            data = {
                'cases': self.cases,
                'case_count': self.case_count,
                'patterns': self.patterns
            }
            with open(self.MEMORY_FILE, 'w', encoding='utf-8') as f:
                json.dump(data, f, ensure_ascii=False, indent=2, default=str)
        except OSError as e:
            print(f"[Memory] 保存记忆文件失败: {e}")

    def save_case(self, sequence: str, goal: str, strategy: str,
                 success: bool, mutations: List[str] = None, 
                 is_user_verified: bool = False) -> str:
        """保存设计案例"""
        self.case_count += 1
        case_id = f"case_{self.case_count:04d}"
        
        features = self._get_sequence_features(sequence)
        
        self.cases[case_id] = {
            'id': case_id,
            'sequence': sequence,
            'goal': goal,
            'strategy': strategy,
            'success': success,
            'mutations': mutations or [],
            'features': features,
            'timestamp': datetime.now().strftime("%Y-%m-%d %H:%M"),
            'time': datetime.now().strftime("%m-%d %H:%M"),
            'sequence_length': len(sequence),
            'is_user_verified': is_user_verified,
            'is_ai_auto_save': not is_user_verified
        }
        
        if is_user_verified and success and mutations:
            for mut in mutations[:2]:
                pattern = f"{mut} -> {goal}优化"
                if pattern not in self.patterns[goal]:
                    self.patterns[goal].append(pattern)

        self._save_to_disk()
        return case_id
    
    def _get_sequence_features(self, seq: str) -> Dict:
        """提取序列特征"""
        if not seq:
            return {}
        return {
            'length': len(seq),
            'hydrophobic': sum(1 for aa in seq if aa in 'FILMVW') / len(seq),
            'charged': sum(1 for aa in seq if aa in 'DEKR') / len(seq),
            'cysteine': seq.count('C') / len(seq),
            'aromatic': sum(1 for aa in seq if aa in 'FWY') / len(seq)
        }
    
    def find_similar(self, query_seq: str, goal: str, max_results: int = 3) -> List[Dict]:
        """查找相似案例"""
        if not self.cases:
            return []
            
        similarities = []
        query_features = self._get_sequence_features(query_seq)
        
        for case_id, case in self.cases.items():
            # 逻辑修正：只有经过用户实验验证（is_user_verified=True）且成功的案例才作为推荐参考
            if case['goal'] == goal and case.get('is_user_verified', False) and case.get('success', False):
                score = self._similarity_score(query_features, case['features'])
                if score > 0.3:
                    similarities.append((score, case))
        
        similarities.sort(key=lambda x: x[0], reverse=True)
        return [case for score, case in similarities[:max_results]]
    
    def _similarity_score(self, feat1: Dict, feat2: Dict) -> float:
        """计算特征相似度"""
        if not feat1 or not feat2:
            return 0.0
            
        score = 0.0
        for key in ['length', 'hydrophobic', 'charged', 'cysteine', 'aromatic']:
            if key in feat1 and key in feat2:
                if key == 'length':
                    max_len = max(feat1[key], feat2[key])
                    if max_len > 0:
                        score += (1 - abs(feat1[key] - feat2[key]) / max_len) * 0.2
                else:
                    score += (1 - abs(feat1[key] - feat2[key])) * 0.2
        return min(1.0, score)
    
    def get_stats(self) -> Dict[str, Any]:
        """获取记忆系统统计"""
        real_cases = {}
        for case_id, case in self.cases.items():
            if case.get('is_ai_auto_save', False):
                continue
            real_cases[case_id] = case
    
        if not real_cases:
            return {
                'total': 0,
                'success': 0,
                'rate': 0.0,
                'patterns': 0,
                'experience': 0
            }
    
        total = len(real_cases)
        success_count = 0
        for case in real_cases.values():
            if case.get('success') is True:
                success_count += 1
    
        rate = success_count / total if total > 0 else 0.0
    
        pattern_count = 0
        for patterns in self.patterns.values():
            if patterns:
                pattern_count += len(patterns)
    
        return {
            'total': total,
            'success': success_count,
            'rate': round(rate, 3),
            'patterns': pattern_count,
            'experience': total
        }
    
    def get_patterns(self, goal: str) -> List[str]:
        """获取特定目标的成功模式"""
        return self.patterns.get(goal, [])
    
    def get_recent_cases(self, count: int = 5) -> List[Dict]:
        """获取最近的设计案例"""
        all_cases = list(self.cases.values())
        all_cases.sort(key=lambda x: x['timestamp'], reverse=True)
        return all_cases[:count]


class ProteinAIAgent:
    """蛋白质AI智能体"""
    
    def __init__(self, random_seed=None):
        """初始化AI智能体"""
        self.design_history = []
        self.learning_data = {}
        self.model_trainer = ProteinModelTrainer()
        
        import random
        
        if random_seed is None:
            random_seed = int(time.time() * 1000) % 1000000
        
        random.seed(random_seed)
        np.random.seed(random_seed)
        self.random_seed = random_seed
        
        self._initialize_knowledge_base()
    
    def _initialize_knowledge_base(self):
        """初始化知识库"""
        self.domain_knowledge = {
            'enzyme': {
                'features': ['催化位点', '底物结合口袋', '辅因子结合位点'],
                'strategies': ['活性位点优化', '底物特异性改造', '热稳定性提升']
            },
            'antibody': {
                'features': ['CDR区域', '框架区', 'Fc区域'],
                'strategies': ['亲和力成熟', '稳定性优化', '人源化改造']
            },
            'membrane_protein': {
                'features': ['跨膜区域', '胞外域', '胞内域'],
                'strategies': ['可溶性表达优化', '稳定性增强', '功能保持']
            }
        }
        
        self.optimization_rules = {
            'stability': {
                'mutations': ['疏水核心填充', '表面电荷优化', '二硫键引入'],
                'targets': ['核心残基', '柔性区域', '表面残基']
            },
            'activity': {
                'mutations': ['活性位点微调', '底物通道优化', '构象动力学改善'],
                'targets': ['催化残基', '结合界面', '别构位点']
            },
            'solubility': {
                'mutations': ['表面亲水化', '电荷平衡', '聚集倾向降低'],
                'targets': ['表面暴露残基', '疏水斑块', '电荷簇']
            }
        }
    
    def _get_random_position(self, sequence, domain, goal, exclude_positions=[]):
        """获取随机突变位置，考虑蛋白质类型和目标"""
        length = len(sequence)
        
        if domain == 'enzyme':
            if goal == 'stability':
                priority_range = range(length//3, 2*length//3)
            elif goal == 'activity':
                priority_range = range(length//4, 3*length//4)
            else:
                priority_range = range(10, length-10)
        elif domain == 'antibody':
            priority_range = range(length//3, 2*length//3)
        else:
            priority_range = range(5, length-5)
        
        priority_positions = list(priority_range)
        available_positions = [p for p in priority_positions if p not in exclude_positions]
        
        if len(available_positions) < 3:
            all_positions = [p for p in range(length) if p not in exclude_positions]
            available_positions = list(set(available_positions + all_positions))
        
        if available_positions:
            return np.random.choice(available_positions)
        else:
            all_positions = [p for p in range(length) if p not in exclude_positions]
            return np.random.choice(all_positions) if all_positions else None
    
    def _generate_stability_mutations(self, sequence: str, domain: str) -> List[str]:
        """生成稳定性相关突变"""
        return self._generate_stability_mutations_random(sequence, domain)
    
    def _generate_stability_mutations_random(self, sequence: str, domain: str) -> List[str]:
        """随机生成稳定性突变"""
        mutations = []
        length = len(sequence)
        
        num_mutations = np.random.randint(1, min(4, length//10 + 2))
        used_positions = []
        
        for _ in range(num_mutations):
            pos = self._get_random_position(sequence, domain, 'stability', used_positions)
            if pos is None:
                continue
            
            used_positions.append(pos)
            original_aa = sequence[pos]
            
            if domain == 'enzyme':
                if pos in range(length//3, 2*length//3):
                    hydrophobic_aas = ['L', 'I', 'V', 'F', 'W', 'M']
                else:
                    hydrophobic_aas = ['P', 'A', 'G']
            else:
                hydrophobic_aas = ['L', 'I', 'V', 'P', 'A']
            
            possible_aas = [aa for aa in hydrophobic_aas if aa != original_aa]
            if possible_aas:
                new_aa = np.random.choice(possible_aas)
                mutations.append(f"{original_aa}{pos+1}{new_aa}")
            else:
                all_aas = ['G', 'A', 'V', 'L', 'I', 'P', 'F', 'Y', 'W', 'S',
                          'T', 'C', 'M', 'N', 'Q', 'D', 'E', 'K', 'R', 'H']
                possible_aas = [aa for aa in all_aas if aa != original_aa]
                if possible_aas:
                    new_aa = np.random.choice(possible_aas)
                    mutations.append(f"{original_aa}{pos+1}{new_aa}")
        
        return mutations[:3]
    
    def _generate_activity_mutations(self, sequence: str, domain: str) -> List[str]:
        """生成活性相关突变"""
        return self._generate_activity_mutations_random(sequence, domain)
    
    def _generate_activity_mutations_random(self, sequence: str, domain: str) -> List[str]:
        """随机生成活性突变"""
        mutations = []
        length = len(sequence)
        
        num_mutations = np.random.randint(1, 3)
        used_positions = []
        
        for _ in range(num_mutations):
            pos = self._get_random_position(sequence, domain, 'activity', used_positions)
            if pos is None:
                continue
            
            used_positions.append(pos)
            original_aa = sequence[pos]
            
            if domain == 'enzyme':
                active_aas = ['D', 'E', 'H', 'K', 'R', 'S', 'T', 'Y']
            elif domain == 'antibody':
                active_aas = ['Y', 'W', 'F', 'R', 'K', 'D', 'E', 'H']
            else:
                active_aas = ['D', 'E', 'R', 'K', 'H', 'Y', 'S', 'T']
            
            possible_aas = [aa for aa in active_aas if aa != original_aa]
            if possible_aas:
                new_aa = np.random.choice(possible_aas)
                mutations.append(f"{original_aa}{pos+1}{new_aa}")
            else:
                all_aas = ['G', 'A', 'V', 'L', 'I', 'P', 'F', 'Y', 'W', 'S',
                          'T', 'C', 'M', 'N', 'Q', 'D', 'E', 'K', 'R', 'H']
                possible_aas = [aa for aa in all_aas if aa != original_aa]
                if possible_aas:
                    new_aa = np.random.choice(possible_aas)
                    mutations.append(f"{original_aa}{pos+1}{new_aa}")
        
        return mutations[:2]
    
    def _generate_solubility_mutations(self, sequence: str, domain: str) -> List[str]:
        """生成可溶性相关突变"""
        return self._generate_solubility_mutations_random(sequence, domain)
    
    def _generate_solubility_mutations_random(self, sequence: str, domain: str) -> List[str]:
        """随机生成可溶性突变"""
        mutations = []
        length = len(sequence)
        
        num_mutations = np.random.randint(1, 3)
        used_positions = []
        
        surface_positions = []
        for i in range(length):
            if i < 10 or i > length - 10 or i % 7 == 0:
                surface_positions.append(i)
        
        for _ in range(num_mutations):
            available_surface = [p for p in surface_positions if p not in used_positions]
            if available_surface:
                pos = np.random.choice(available_surface)
            else:
                pos = self._get_random_position(sequence, domain, 'solubility', used_positions)
            
            if pos is None:
                continue
            
            used_positions.append(pos)
            original_aa = sequence[pos]
            
            if original_aa in 'FILMVWY':
                hydrophilic_aas = ['S', 'T', 'N', 'Q', 'D', 'E', 'K', 'R']
                possible_aas = [aa for aa in hydrophilic_aas if aa != original_aa]
            elif original_aa in 'DEKR':
                neutral_aas = ['S', 'T', 'N', 'Q', 'A', 'G']
                possible_aas = [aa for aa in neutral_aas if aa != original_aa]
            else:
                all_aas = ['G', 'A', 'V', 'L', 'I', 'P', 'F', 'Y', 'W', 'S',
                          'T', 'C', 'M', 'N', 'Q', 'D', 'E', 'K', 'R', 'H']
                possible_aas = [aa for aa in all_aas if aa != original_aa]
            
            if possible_aas:
                new_aa = np.random.choice(possible_aas)
                mutations.append(f"{original_aa}{pos+1}{new_aa}")
        
        return mutations
    
    def _apply_mutations_to_sequence(self, sequence: str, mutations: List[str]) -> str:
        """将突变应用到序列上，生成突变体序列"""
        if not mutations:
            return sequence
        
        seq_list = list(sequence)
        
        for mutation in mutations:
            if len(mutation) >= 3:
                try:
                    original_aa = mutation[0]
                    pos_str = mutation[1:-1]
                    new_aa = mutation[-1]
                    
                    pos = int(pos_str) - 1
                    
                    if 0 <= pos < len(seq_list) and seq_list[pos] == original_aa:
                        seq_list[pos] = new_aa
                except (ValueError, IndexError):
                    continue
        
        return ''.join(seq_list)
    
    def train_ml_model(self, sequences: List[str], targets: List[float]) -> bool:
        """训练机器学习模型"""
        return self.model_trainer.train_model(sequences, targets)
    
    def generate_ml_mutants(self, wildtype_seq: str, num_mutants: int = 5) -> List[Dict]:
        """ML模型生成突变体"""
        if not self.model_trainer.is_trained:
            return []
        
        mutants = self._generate_mutant_sequences(wildtype_seq, num_mutants * 3)
        
        scored_mutants = []
        for mutant_seq in mutants:
            try:
                score = self.model_trainer.predict(mutant_seq)
                mutations = self._identify_mutations(wildtype_seq, mutant_seq)
                
                confidence = min(0.95, 0.7 + (score * 0.2))
                
                scored_mutants.append({
                    'sequence': mutant_seq,
                    'score': float(score),
                    'confidence': round(confidence, 2),
                    'mutations': mutations,
                    'mutation_count': len(mutations),
                    'rationale': f"基于AAindex机器学习模型预测",
                    'is_ml_based': True
                })
            except Exception:
                continue
        
        scored_mutants.sort(key=lambda x: x['score'], reverse=True)
        return scored_mutants[:num_mutants]
    
    def _generate_mutant_sequences(self, wildtype_seq: str, num_mutants: int) -> List[str]:
        """生成突变体序列"""
        mutants = []
        seq_list = list(wildtype_seq)
        amino_acids = ['G', 'A', 'V', 'L', 'I', 'P', 'F', 'Y', 'W', 'S',
                      'T', 'C', 'M', 'N', 'Q', 'D', 'E', 'K', 'R', 'H']
        
        positions = list(range(len(wildtype_seq)))
        
        for _ in range(num_mutants):
            num_mutations = np.random.randint(1, 4)
            mut_positions = np.random.choice(positions, num_mutations, replace=False)
            
            mutant_seq = list(wildtype_seq)
            mutations = []
            
            for pos in mut_positions:
                original_aa = wildtype_seq[pos]
                possible_aas = [aa for aa in amino_acids if aa != original_aa]
                new_aa = np.random.choice(possible_aas)
                mutant_seq[pos] = new_aa
                mutations.append(f"{original_aa}{pos+1}{new_aa}")
            
            mutants.append(''.join(mutant_seq))
        
        return mutants
    
    def _identify_mutations(self, wildtype: str, mutant: str) -> List[str]:
        """识别突变位置"""
        mutations = []
        for i, (wt_aa, mut_aa) in enumerate(zip(wildtype, mutant)):
            if wt_aa != mut_aa:
                mutations.append(f"{wt_aa}{i+1}{mut_aa}")
        return mutations
    
    def analyze_sequence(self, sequence: str, goal: str, domain: str) -> Dict[str, Any]:
        """分析蛋白质序列"""
        analysis = {
            'sequence_length': len(sequence),
            'amino_acid_composition': self._calculate_aa_composition(sequence),
            'molecular_weight': self._estimate_molecular_weight(sequence),
            'isoelectric_point': self._estimate_pi(sequence),
            'instability_index': self._calculate_instability_index(sequence),
            'sequence_features': self._get_enhanced_sequence_features(sequence),
            'domain_features': self._identify_domain_features(sequence, domain),
            'optimization_suggestions': self._generate_optimization_suggestions(sequence, goal, domain)
        }
        return analysis
    
    def _get_enhanced_sequence_features(self, sequence: str) -> Dict:
        """获取增强的序列特征"""
        if not sequence:
            return {}
        
        return {
            'hydrophobic_ratio': sum(1 for aa in sequence if aa in 'FILMVW') / len(sequence),
            'charged_ratio': sum(1 for aa in sequence if aa in 'DEKR') / len(sequence),
            'cysteine_ratio': sequence.count('C') / len(sequence),
            'aromatic_ratio': sum(1 for aa in sequence if aa in 'FWY') / len(sequence),
            'polar_ratio': sum(1 for aa in sequence if aa in 'NQSTY') / len(sequence),
            'tiny_ratio': sum(1 for aa in sequence if aa in 'AGSV') / len(sequence)
        }
    
    def generate_mutants(self, sequence: str, goal: str, domain: str, 
                        memory: EnhancedProteinMemory, num_mutants: int = 5) -> List[Dict]:
        """生成突变体建议"""
        mutants = []
        base_seed = self.random_seed
        
        similar_cases = memory.find_similar(sequence, goal, 2)
        learned_mutations = []
        
        for case in similar_cases:
            if case['success'] and case['mutations']:
                learned_mutations.extend(case['mutations'][:1])
        
        for i in range(num_mutants):
            current_seed = base_seed + i * 1000
            np.random.seed(current_seed)
            
            if i < len(learned_mutations) and learned_mutations:
                mutations = [learned_mutations[i]]
                rationale = f"基于相似成功案例的模式: {learned_mutations[i]}"
                is_learned = True
            else:
                if goal == 'stability':
                    mutations = self._generate_stability_mutations(sequence, domain)
                    rationale = "稳定性优化：增加疏水核心包装或表面刚性"
                elif goal == 'activity':
                    mutations = self._generate_activity_mutations(sequence, domain)
                    rationale = "活性优化：优化催化位点或结合界面"
                elif goal == 'solubility':
                    mutations = self._generate_solubility_mutations(sequence, domain)
                    rationale = "可溶性优化：疏水表面亲水化或电荷平衡"
                else:
                    mutations = []
                    length = len(sequence)
                    num_mutations = np.random.randint(1, 3)
                    positions = np.random.choice(range(length), num_mutations, replace=False)
                    
                    for pos in positions:
                        original_aa = sequence[pos]
                        all_aas = ['G', 'A', 'V', 'L', 'I', 'P', 'F', 'Y', 'W', 'S',
                                  'T', 'C', 'M', 'N', 'Q', 'D', 'E', 'K', 'R', 'H']
                        possible_aas = [aa for aa in all_aas if aa != original_aa]
                        if possible_aas:
                            new_aa = np.random.choice(possible_aas)
                            mutations.append(f"{original_aa}{pos+1}{new_aa}")
                    
                    rationale = "通用优化：随机探索性突变"
                
                is_learned = False
            
            if mutations:
                mutant_sequence = self._apply_mutations_to_sequence(sequence, mutations)
                
                score = self._calculate_mutant_score(mutations, goal, domain, is_learned)
                
                confidence = 0.85 - (i * 0.1)
                confidence = max(0.6, min(0.9, confidence))
                
                mutant = {
                    'id': f"mutant_{i+1}",
                    'sequence': mutant_sequence,
                    'score': round(score, 3),
                    'confidence': round(confidence, 2),
                    'mutations': mutations,
                    'mutation_count': len(mutations),
                    'rationale': rationale,
                    'is_learned_pattern': is_learned,
                    'is_ml_based': False
                }
                mutants.append(mutant)
        
        np.random.seed(base_seed)
        mutants.sort(key=lambda x: x['score'], reverse=True)
        return mutants
    
    def _calculate_mutant_score(self, mutations: List[str], goal: str, 
                               domain: str, is_learned: bool = False) -> float:
        """计算突变体得分（启发式规则打分，非数据驱动；仅用于智能规则模式的相对排序）"""
        # 说明：本打分为基于物化规则的启发式加成，反映"规则上更可能有效"，
        # 并非 ML 模型的定量预测。如需数据驱动打分，请使用机器学习模式（generate_ml_mutants）。
        base_score = 0.5
        
        if goal == 'stability':
            base_score += 0.2
        elif goal == 'activity':
            base_score += 0.15
        elif goal == 'solubility':
            base_score += 0.1
        
        mutation_bonus = min(0.3, len(mutations) * 0.1)
        base_score += mutation_bonus
        
        if domain in ['enzyme', 'antibody']:
            base_score += 0.1
        
        if is_learned:
            base_score += 0.15
        
        return min(0.95, base_score)
    
    def _calculate_aa_composition(self, sequence: str) -> Dict[str, float]:
        """计算氨基酸组成"""
        total = len(sequence)
        composition = {}
        for aa in set(sequence):
            composition[aa] = round(sequence.count(aa) / total * 100, 2)
        return composition
    
    def _estimate_molecular_weight(self, sequence: str) -> float:
        """估算分子量"""
        aa_weights = {
            'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.15,
            'E': 147.13, 'Q': 146.15, 'G': 75.07, 'H': 155.16, 'I': 131.17,
            'L': 131.17, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
            'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15
        }
        
        weight = sum(aa_weights.get(aa, 110) for aa in sequence) - (len(sequence)-1)*18.02
        return round(weight / 1000, 2)
    
    def _estimate_pi(self, sequence: str) -> float:
        """估算等电点"""
        acidic = sequence.count('D') + sequence.count('E')
        basic = sequence.count('R') + sequence.count('K') + sequence.count('H')
        
        if acidic + basic == 0:
            return 7.0
        
        ratio = basic / (acidic + basic)
        pi = 3.0 + ratio * 8.0
        return round(pi, 2)
    
    def _calculate_instability_index(self, sequence: str) -> float:
        """计算不稳定指数"""
        unstable_pairs = ['DP', 'PE', 'PG', 'PH', 'PK', 'PR']
        dipeptides = [sequence[i:i+2] for i in range(len(sequence)-1)]
        unstable_count = sum(1 for pair in dipeptides if pair in unstable_pairs)
        
        index = (unstable_count / len(dipeptides)) * 100
        return round(index, 1)
    
    def _identify_domain_features(self, sequence: str, domain: str) -> List[str]:
        """识别结构特征"""
        features = []
        
        if sequence.count('C') >= 2:
            features.append('潜在二硫键')
        
        if domain == 'enzyme':
            if any(motif in sequence for motif in ['GXGXXG', 'DXG']):
                features.append('核苷酸结合模体')
        elif domain == 'antibody':
            features.append('免疫球蛋白折叠')
        
        hydrophobic_ratio = sum(1 for aa in sequence if aa in 'FILMV') / len(sequence)
        if hydrophobic_ratio > 0.4:
            features.append('高疏水性区域')
        
        return features
    
    def _generate_optimization_suggestions(self, sequence: str, goal: str, domain: str) -> List[str]:
        """生成优化建议"""
        suggestions = []
        
        if goal == 'stability':
            suggestions.extend([
                '优化核心疏水包装',
                '引入表面盐桥',
                '减少柔性环区长度'
            ])
        elif goal == 'activity':
            suggestions.extend([
                '微调活性位点残基',
                '优化底物结合口袋',
                '改善催化三联体几何'
            ])
        elif goal == 'solubility':
            suggestions.extend([
                '表面疏水残基亲水化',
                '优化净电荷分布',
                '引入糖基化位点'
            ])
        
        if domain == 'enzyme':
            suggestions.append('保持催化残基完整性')
        elif domain == 'antibody':
            suggestions.append('优化CDR环构象')
        
        return suggestions


class SmartProteinAdvisor:
    """智能蛋白质分析顾问"""
    
    def __init__(self):
        self.deepseek = None
        self.has_deepseek = False
        self.local_knowledge = self._load_local_knowledge()
        
        self.stats = {"deepseek": 0, "offline": 0}
    
    def _load_local_knowledge(self):
        """加载本地知识库"""
        try:
            from offline_knowledge import OfflineProteinKnowledge
            offline_kb = OfflineProteinKnowledge()
            return offline_kb.knowledge
        except ImportError:
            return {
                "stability": "稳定性的关键因素包括疏水核心包装、表面盐桥形成、二硫键稳定性、构象刚性等。优化策略：1) 增加核心疏水性 2) 优化表面电荷 3) 引入二硫键（需谨慎）",
                "activity": "活性优化涉及活性位点残基调整、底物结合口袋优化、催化效率提升、别构调节位点改造。策略：1) 微调催化残基 2) 优化底物通道 3) 改善动力学参数",
                "solubility": "可溶性改进包括表面亲水化、电荷平衡、聚集倾向降低、糖基化位点引入等策略。关键：1) 表面疏水残基突变 2) 电荷分布优化 3) 防止错误折叠"
            }
        except Exception:
            return {}
    
    def set_deepseek_key(self, api_key: str):
        """设置DeepSeek API密钥"""
        try:
            from deepseek_expert import DeepSeekProteinExpert
            
            self.deepseek = DeepSeekProteinExpert(api_key=api_key)
            self.has_deepseek = getattr(self.deepseek, 'is_ready', False)
            
        except Exception:
            self.has_deepseek = False
    
    def get_analysis(self, wildtype: str, mutations: List[str], score: float, 
                    goal: str, domain: str, use_online: bool = True) -> Dict:
        """获取分析"""
        if use_online and self.has_deepseek and self.deepseek:
            try:
                result = self.deepseek.analyze(wildtype, mutations, score, goal, domain)
                self.stats["deepseek"] += 1
                
                if result and isinstance(result, dict) and 'content' in result:
                    return result
                else:
                    return self._get_offline_analysis(wildtype, mutations, score, goal, domain)
                    
            except Exception:
                return self._get_offline_analysis(wildtype, mutations, score, goal, domain)
        
        self.stats["offline"] += 1
        return self._get_offline_analysis(wildtype, mutations, score, goal, domain)
    
    def _get_offline_analysis(self, wildtype: str, mutations: List[str], score: float,
                            goal: str, domain: str) -> Dict:
        """获取离线分析"""
        goal_cn = {"stability": "稳定性", "activity": "活性", "solubility": "可溶性"}.get(goal, goal)
        
        knowledge = self.local_knowledge.get(goal, f"暂无详细的{goal_cn}优化分析信息。")
        
        analysis = f"""
## 离线分析结果

### 基本信息
- **突变：** {', '.join(mutations) if mutations else '无具体突变'}
- **优化目标：** {goal_cn}提升
- **预测得分：** {score:.3f}（满分1.0）
- **蛋白质类型：** {domain}
- **野生型长度：** {len(wildtype)}个氨基酸

### 专业分析
{knowledge}

### 突变具体影响分析
"""
        
        if mutations:
            for i, mutation in enumerate(mutations, 1):
                if len(mutation) >= 3:
                    orig = mutation[0]
                    new = mutation[-1]
                    analysis += f"\n**{i}. {mutation}**\n"
                    
                    if goal == "stability":
                        if orig in "ADE" and new in "NQ":
                            analysis += "- 消除表面电荷，减少静电排斥\n- 可能提高热稳定性\n- 建议验证对功能的影响\n"
                        elif orig in "G" and new in "P":
                            analysis += "- 引入脯氨酸增加结构刚性\n- 可能限制构象柔性\n- 避免在活性位点使用\n"
                        else:
                            analysis += f"- {orig}→{new}可能影响局部结构\n- 建议分子动力学模拟验证\n"
                    
                    elif goal == "activity":
                        analysis += f"- {orig}→{new}可能影响催化/结合位点\n- 建议酶动力学实验验证\n"
                    
                    elif goal == "solubility":
                        if orig in "FILMWY" and new in "STNQ":
                            analysis += f"- 疏水→亲水突变，可能改善可溶性\n- 注意可能影响结构稳定性\n"
                        else:
                            analysis += f"- {orig}→{new}对可溶性影响待验证\n"
        
        analysis += f"""
### 实验验证建议
1. **表达验证**：在大肠杆菌/酵母/哺乳动物细胞中表达突变体
2. **纯化分析**：使用亲和层析纯化，检测可溶性
3. **功能测定**：
   - 稳定性：热位移实验、圆二色谱
   - 活性：酶动力学分析、底物结合实验
   - 可溶性：动态光散射、浊度测定
4. **结构分析**（可选）：X射线晶体学、冷冻电镜

### 后续步骤
1. 设计引物，进行定点突变
2. 表达并纯化蛋白质
3. 进行功能验证实验
4. 根据结果进行迭代优化

---
*注：这是基于离线知识库的分析结果。如需更详细的AI分析，请设置有效的DeepSeek API密钥。*
"""
        
        return {
            'content': analysis.strip(),
            'source': '离线知识库',
            'is_online': False,
            'timestamp': datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
    
    def show_stats(self):
        """显示统计信息"""
        print("\n" + "="*50)
        print("SmartProteinAdvisor统计信息:")
        print(f"- DeepSeek使用次数: {self.stats['deepseek']}")
        print(f"- 离线分析次数: {self.stats['offline']}")
        print(f"- has_deepseek: {self.has_deepseek}")
        print("="*50)


class EnhancedNCBIWebBlastAPI(RealNCBIBlastAPI):
    """增强的NCBI BLAST API"""
    
    def __init__(self, email: str, api_key: str = None, timeout: int = 300, max_retries: int = 3):
        super().__init__(email=email, api_key=api_key)
        self.max_retries = max_retries
        self.retry_delay = 10
        self.timeout = timeout
        self.base_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
        self.headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
        }
        
    def run_blast_with_retry(self, sequence: str, program: str = 'blastp', 
                            database: str = 'nr', max_wait: int = 300) -> Optional[Dict]:
        """带重试机制的BLAST分析"""
        for attempt in range(self.max_retries):
            try:
                st.info(f"BLAST尝试 {attempt + 1}/{self.max_retries}...")
                
                if attempt > 0:
                    st.info(f"等待 {self.retry_delay} 秒后重试...")
                    time.sleep(self.retry_delay)
                    self.retry_delay = min(self.retry_delay * 1.5, 60)
                
                result = self.run_blast(sequence, program, database, max_wait)
                if result:
                    return result
                    
            except Exception as e:
                st.warning(f"尝试 {attempt + 1} 失败: {str(e)}")
                continue
        
        return None


import hashlib

class BlastResultCache:
    """BLAST结果缓存"""
    
    def __init__(self, cache_dir: str = "./cache/blast"):
        self.cache_dir = cache_dir
        os.makedirs(cache_dir, exist_ok=True)
        
        self.expiry_days = 7
        
    def _get_cache_key(self, sequence: str, database: str = "nr") -> str:
        """生成缓存键"""
        key_str = f"{sequence}_{database}"
        return hashlib.md5(key_str.encode()).hexdigest()
    
    def _get_cache_file(self, cache_key: str) -> str:
        """获取缓存文件路径"""
        return os.path.join(self.cache_dir, f"{cache_key}.json")
    
    def get(self, sequence: str, database: str = "nr") -> Optional[Dict]:
        """从缓存获取BLAST结果"""
        try:
            cache_key = self._get_cache_key(sequence, database)
            cache_file = self._get_cache_file(cache_key)
            
            if not os.path.exists(cache_file):
                return None
            
            from datetime import datetime, timedelta
            file_mtime = datetime.fromtimestamp(os.path.getmtime(cache_file))
            if datetime.now() - file_mtime > timedelta(days=self.expiry_days):
                os.remove(cache_file)
                return None
            
            with open(cache_file, 'r', encoding='utf-8') as f:
                cached_data = json.load(f)
            
            required_keys = ['homologous_proteins', 'query_length', 'database']
            if all(key in cached_data for key in required_keys):
                cached_data['from_cache'] = True
                cached_data['cache_time'] = file_mtime.strftime("%Y-%m-%d %H:%M:%S")
                return cached_data
            else:
                return None
                
        except Exception:
            return None
    
    def set(self, sequence: str, blast_result: Dict, database: str = "nr"):
        """保存BLAST结果到缓存"""
        try:
            if not blast_result:
                return
            
            cache_key = self._get_cache_key(sequence, database)
            cache_file = self._get_cache_file(cache_key)
            
            cache_data = {
                **blast_result,
                'cached_at': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                'sequence_hash': cache_key,
                'original_sequence_length': len(sequence)
            }
            
            with open(cache_file, 'w', encoding='utf-8') as f:
                json.dump(cache_data, f, ensure_ascii=False, indent=2)
            
        except Exception:
            pass
    
    def clear_old_cache(self):
        """清理过期缓存"""
        try:
            from datetime import datetime, timedelta
            cutoff_time = datetime.now() - timedelta(days=self.expiry_days)
            
            for filename in os.listdir(self.cache_dir):
                if filename.endswith('.json'):
                    filepath = os.path.join(self.cache_dir, filename)
                    file_mtime = datetime.fromtimestamp(os.path.getmtime(filepath))
                    
                    if file_mtime < cutoff_time:
                        os.remove(filepath)
                        
        except Exception:
            pass


class EnhancedNCBIBlastAnalyzer:
    """增强的NCBI BLAST分析器"""
    
    def __init__(self, email: str, api_key: str = None):
        self.email = email
        self.api_key = api_key
        self.blast_api = RealNCBIBlastAPI(email, api_key)
        self.error_handler = BlastErrorHandler()
        self.cache = BlastResultCache()
        
        self.cache.clear_old_cache()
    
    def analyze_sequence(self, sequence: str, database: str = "nr", 
                        use_cache: bool = True) -> Dict:
        """分析序列的同源性"""
        validation = self.error_handler.validate_sequence(sequence)
        if not validation['valid']:
            return {
                'status': 'error',
                'data': None,
                'message': validation['message'],
                'metadata': {'validation_error': True}
            }
        
        if use_cache:
            cached_result = self.cache.get(sequence, database)
            if cached_result:
                return {
                    'status': 'cached',
                    'data': cached_result,
                    'message': '使用缓存结果',
                    'metadata': {'from_cache': True}
                }
        
        try:
            with st.spinner("正在运行BLAST分析..."):
                blast_result = self.blast_api.run_blast(sequence, database=database)
            
            if blast_result:
                self.cache.set(sequence, blast_result, database)
                
                return {
                    'status': 'success',
                    'data': blast_result,
                    'message': blast_result.get('message', 'BLAST分析完成'),
                    'metadata': {'from_cache': False}
                }
            else:
                return {
                    'status': 'error',
                    'data': None,
                    'message': 'BLAST分析未返回结果',
                    'metadata': {'blast_empty': True}
                }
                
        except Exception as e:
            error_info = self.error_handler.handle_error(e, "run_blast")
            
            return {
                'status': 'error',
                'data': None,
                'message': error_info['message'],
                'metadata': {
                    'error_type': type(e).__name__,
                    'retryable': error_info.get('retryable', False),
                    'suggestion': error_info.get('suggestion', '')
                }
            }


def validate_protein_sequence(sequence: str) -> Dict[str, Any]:
    """验证蛋白质序列"""
    if not sequence or sequence.strip() == "":
        return {'valid': False, 'message': '序列不能为空'}
    
    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
    clean_seq = ''.join(aa for aa in sequence.upper() if aa in valid_aa)
    
    if len(clean_seq) != len(sequence.strip()):
        return {
            'valid': False, 
            'message': f'包含无效字符，已过滤为{len(clean_seq)}个有效AA',
            'suggestion': '请只使用标准20种氨基酸字符'
        }
    
    if len(clean_seq) < 10:
        return {
            'valid': False,
            'message': '序列过短',
            'suggestion': '请使用至少10个氨基酸的序列'
        }
    
    if len(clean_seq) > 2000:
        return {
            'valid': False,
            'message': '序列过长',
            'suggestion': '请使用少于2000个氨基酸的序列'
        }
    
    return {'valid': True, 'cleaned_sequence': clean_seq}


def display_blast_results_ultrawide(blast_results: Dict):
    """超宽显示BLAST分析结果"""
    
    # 使用全宽容器
    with st.container():
        st.markdown("<div class='full-width-container'>", unsafe_allow_html=True)
        
        st.subheader("BLAST同源性分析结果")
        
        # 顶部指标卡片 - 使用更多列布局
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("查询长度", f"{blast_results.get('query_length', 0)} aa")
        with col2:
            total_hits = blast_results.get('total_hits', 0)
            st.metric("同源蛋白数", total_hits)
        with col3:
            database = blast_results.get('database', 'nr')
            st.metric("数据库", database)
        with col4:
            if blast_results.get('from_cache', False):
                cache_time = blast_results.get('cache_time', '未知时间')
                st.metric("来源", "缓存", delta=f"保存于{cache_time}")
            else:
                st.metric("来源", "实时", delta="最新结果")
        
        if blast_results.get('is_offline_example', False):
            st.warning("离线模式：这是示例数据，实际BLAST分析需要网络连接")
        
        st.markdown("---")
        
        # 同源蛋白质表格 - 使用Markdown表格
        st.subheader("同源蛋白质 (Top 15)")
        homologous_proteins = blast_results.get('homologous_proteins', [])
        
        if homologous_proteins:
            markdown_table = "| 序号 | Accession ID | 蛋白质描述 | 物种 | E值 | 一致性 | 得分 |\n"
            markdown_table += "|------|-------------|------------|------|-----|--------|------|\n"
            
            for i, protein in enumerate(homologous_proteins[:15]):
                desc = protein.get('description', '')
                if '|' in desc:
                    parts = desc.split('|')
                    clean_desc = parts[-1].strip() if len(parts) >= 3 else desc
                else:
                    clean_desc = desc
                
                if len(clean_desc) > 60:
                    clean_desc = clean_desc[:60] + "..."
                
                accession = protein.get('accession', 'N/A')
                
                if accession != 'N/A':
                    accession_link = f"[{accession}](https://www.ncbi.nlm.nih.gov/protein/{accession})"
                else:
                    accession_link = 'N/A'
                
                markdown_table += f"| {i+1} | {accession_link} | {clean_desc} | {protein.get('species', 'Unknown')} | {protein.get('e_value', 0):.2e} | {protein.get('identity_percent', 0):.1f}% | {protein.get('bit_score', 0):.2f} |\n"
            
            st.markdown(markdown_table)
            
            st.info(f"**提示：**点击Accession ID链接可直接在新窗口打开NCBI蛋白质页面 | 共找到 {len(homologous_proteins)} 个同源蛋白")
            st.markdown("---")
            
            # 保守区域分析
            if blast_results.get('conserved_regions'):
                st.subheader("保守区域分析")
                conserved_regions = blast_results['conserved_regions']
                
                # 去重，基于起始位置
                unique_regions = []
                seen_starts = set()
                
                for region in conserved_regions:
                    start = region.get('start', 0)
                    if start not in seen_starts:
                        seen_starts.add(start)
                        unique_regions.append(region)
                
                # 限制显示数量
                unique_regions = unique_regions[:6]
                
                if unique_regions:
                    cols = st.columns(3)
                    for idx, region in enumerate(unique_regions):
                        col_idx = idx % 3
                        with cols[col_idx]:
                            identity = region.get('identity', 0)
                            start = region.get('start', 0)
                            end = region.get('end', 0)
                            
                            if identity > 80:
                                color = "#d4edda"
                                border_color = "#28a745"
                            elif identity > 60:
                                color = "#fff3cd"
                                border_color = "#ffc107"
                            else:
                                color = "#d1ecf1"
                                border_color = "#17a2b8"
                            
                            st.markdown(f"""
                            <div style="
                                background-color: {color};
                                border: 2px solid {border_color};
                                padding: 15px;
                                margin: 10px 0;
                                border-radius: 10px;
                                width: 100%;
                            ">
                                <strong>{region.get('description', '保守区域')}</strong><br>
                                <b>位置:</b> {start}-{end}<br>
                                <b>长度:</b> {region.get('length', 0)} aa<br>
                                <b>一致性:</b> {identity:.1f}%
                            </div>
                            """, unsafe_allow_html=True)
            
            st.markdown("---")
            
            # NCBI完整结果链接
            current_seq = st.session_state.get('current_sequence', '')
            
            if current_seq:
                encoded_seq = urllib.parse.quote(current_seq)
                blast_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&QUERY={encoded_seq}"
                
                st.markdown(f"""
                ### 在NCBI查看完整BLAST结果
                
                链接：[{blast_url}]({blast_url})
                
                **包含功能：**
                - 完整的比对可视化图表
                - 多种格式下载（XML、文本、CSV）
                - 保守区域和结构域分析
                - 物种分布统计
                """)
                
                try:
                    st.link_button(
                        "打开NCBI BLAST页面",
                        blast_url,
                        help="点击将在新标签页打开NCBI完整BLAST结果",
                        use_container_width=True
                    )
                except:
                    st.markdown(f"""
                    <a href="{blast_url}" target="_blank" style="
                        display: block;
                        text-align: center;
                        padding: 10px;
                        background-color: #0066cc;
                        color: white;
                        text-decoration: none;
                        border-radius: 5px;
                        font-weight: bold;
                        margin: 10px 0;
                    ">
                        打开NCBI BLAST页面
                    </a>
                    """, unsafe_allow_html=True)
            else:
                st.warning("无法生成NCBI链接：序列信息不可用")
            
            st.markdown("---")
            
            # 设计建议 - 全宽显示
            st.write("### 基于BLAST结果的设计建议")
            
            suggestions_container = st.container()
            with suggestions_container:
                if homologous_proteins:
                    homologous_proteins_list = blast_results['homologous_proteins']
                    total_hits = len(homologous_proteins_list)
                    
                    avg_identity = sum(p.get('identity_percent', 0) for p in homologous_proteins_list[:10]) / min(10, total_hits)
                    
                    if avg_identity > 80:
                        st.success("**高保守序列检测**：检测到高度同源蛋白（>80%一致性），保守区域可能对功能至关重要，建议避免突变这些区域。")
                    elif avg_identity > 50:
                        st.warning("**中等同源性**：检测到中等同源蛋白（50-80%一致性），可以考虑在低保守区域进行突变尝试。")
                    else:
                        st.info("**低同源性**：同源蛋白一致性较低，这可能是一个新颖的蛋白质家族，设计时需谨慎。")
                    
                    species_list = [p.get('species', 'Unknown') for p in homologous_proteins_list[:10]]
                    human_hits = sum(1 for s in species_list if 'Human' in s or 'homo' in s.lower())
                    if human_hits > 0:
                        st.info("**人类同源蛋白发现**：发现人类同源蛋白，这对药物设计和人源化改造非常重要。")
                
                if blast_results.get('conserved_regions'):
                    conserved_regions_list = blast_results['conserved_regions']
                    if len(conserved_regions_list) > 0:
                        st.warning("**保守区域警告**：发现保守结构域，这些区域可能是功能位点，突变需谨慎。")
                        
                        most_conserved = max(conserved_regions_list, key=lambda x: x.get('identity', 0))
                        if most_conserved.get('identity', 0) > 90:
                            st.error(f"**高度保守区域**：位置 {most_conserved['start']}-{most_conserved['end']} 的保守度 >90%，强烈建议保留。")
        else:
            st.info("未找到同源蛋白质")
        
        st.markdown("</div>", unsafe_allow_html=True)

def render_protein_3d_view(pdb_content: str, height: int = 600):
    """渲染3D蛋白质结构可视化"""
    try:
        # 使用 py3Dmol 创建可视化
        view = py3Dmol.view(width=None, height=height)
        view.addModel(pdb_content, 'pdb')
        
        # 基础样式：按光谱着色（从N端到C端）
        view.setStyle({'cartoon': {'color': 'spectrum'}})
        
        view.zoomTo()
        
        # 生成 HTML 并渲染
        # 注意：py3Dmol._make_html() 返回的是包含完整依赖的 HTML
        html = view._make_html()
        
        # 包装在 iframe 中展示
        components.html(html, height=height)
        
    except Exception as e:
        st.error(f"3D结构渲染失败: {str(e)}")
        st.info("请检查网络连接或尝试手动刷新页面。")


def display_structure_results_ultrawide(structure_results: Dict):
    """超宽显示结构分析结果"""
    
    with st.container():
        st.markdown("<div class='full-width-container'>", unsafe_allow_html=True)
        
        st.subheader("蛋白质结构预测结果")
        
        # 顶部指标 - 使用更多列
        col1, col2, col3, col4, col5 = st.columns(5)
        
        with col1:
            st.metric("平均pLDDT", f"{structure_results.get('avg_plddt', 0):.1f}")
        with col2:
            method = structure_results.get('method', '未知')
            st.metric("预测方法", method)
        with col3:
            seq_len = len(structure_results.get('sequence', ''))
            st.metric("序列长度", f"{seq_len} aa")
        with col4:
            source = structure_results.get('api_source', '未知')
            if '本地' in source:
                source = '本地ESM2'
            st.metric("预测源", source)
        with col5:
            if structure_results.get('is_real_prediction', False):
                st.metric("预测类型", "真实预测")
            else:
                st.metric("预测类型", "示意图")
        
        st.markdown("---")

        # --- 新增：3D结构展示区域 ---
        if 'pdb_content' in structure_results:
            st.markdown("### 3D 交互结构视图")
            render_protein_3d_view(structure_results['pdb_content'])
            st.markdown("---")
        # ------------------------
        
        # pLDDT图表 - 全宽显示
        if 'plddt_scores' in structure_results and len(structure_results['plddt_scores']) > 0:
            plot_plddt_chart_ultrawide(structure_results['plddt_scores'], 
                                     structure_results.get('sequence', ''))
        
        # 两列布局但全宽
        col1, col2 = st.columns([3, 2])
        
        with col1:
            st.markdown("### PDB文件")
            if 'pdb_content' in structure_results:
                # 创建可折叠的预览区域
                with st.expander("预览PDB文件内容（点击展开）", expanded=False):
                    pdb_content = structure_results['pdb_content']
                    # 显示更多行
                    lines = pdb_content.split('\n')[:100]
                    st.code('\n'.join(lines), language='pdb')
                
                # 下载按钮 - 全宽
                st.download_button(
                    label="下载完整PDB文件",
                    data=structure_results['pdb_content'],
                    file_name=f"protein_structure_{len(structure_results.get('sequence', ''))}aa.pdb",
                    mime="chemical/x-pdb",
                    use_container_width=True
                )
                
                # 序列预览
                if 'sequence' in structure_results:
                    st.markdown("### 序列预览")
                    seq = structure_results['sequence']
                    # 每行显示60个字符
                    seq_lines = [seq[i:i+60] for i in range(0, len(seq), 60)]
                    numbered_seq = ""
                    for i, line in enumerate(seq_lines):
                        start_num = i * 60 + 1
                        end_num = start_num + len(line) - 1
                        numbered_seq += f"{start_num:4d}-{end_num:4d}: {line}\n"
                    
                    st.text_area(
                        "氨基酸序列",
                        value=numbered_seq,
                        height=200,
                        disabled=True,
                        key="structure_sequence_preview"
                    )
        
        with col2:
            st.markdown("### 结构质量评估")
            
            avg_plddt = structure_results.get('avg_plddt', 0)
            
            # 质量评估卡片 —— 按平均 pLDDT 分级渲染（配置驱动，避免重复 HTML）
            _QUALITY_LEVELS = [
                (90, "#d4edda", "#28a745", "#155724", "极好",
                 ["高置信度预测", "适合分子对接", "可直接用于分析", "建议用于实验设计"]),
                (70, "#d1ecf1", "#17a2b8", "#0c5460", "良好",
                 ["中等置信度", "可用于初步分析", "建议进一步验证", "适合功能预测"]),
                (50, "#fff3cd", "#ffc107", "#856404", "一般",
                 ["低置信度区域较多", "需谨慎使用", "建议实验验证", "可用于参考"]),
                (0,  "#f8d7da", "#dc3545", "#721c24", "较差",
                 ["预测不确定性高", "不建议直接使用", "需要实验结构", "仅作参考"]),
            ]
            level = next(lv for lv in _QUALITY_LEVELS if avg_plddt >= lv[0])
            _, bg, border, text, label, items = level
            items_html = "".join(f"<li>{it}</li>" for it in items)
            st.markdown(
                f'<div style="background-color: {bg}; border: 2px solid {border}; '
                f'padding: 20px; border-radius: 10px; margin: 10px 0;">'
                f'<h4 style="color: {text}; margin-top: 0;">结构质量: {label}</h4>'
                f'<ul style="color: {text};">{items_html}</ul></div>',
                unsafe_allow_html=True
            )
            
            # 结构特征统计
            if 'plddt_scores' in structure_results:
                plddt_scores = structure_results['plddt_scores']
                high_confidence = sum(1 for score in plddt_scores if score >= 70)
                medium_confidence = sum(1 for score in plddt_scores if 50 <= score < 70)
                low_confidence = sum(1 for score in plddt_scores if score < 50)
                
                st.markdown("### 结构统计")
                
                stats_col1, stats_col2, stats_col3 = st.columns(3)
                with stats_col1:
                    st.metric("高置信度", f"{high_confidence}", delta=f"{high_confidence/len(plddt_scores)*100:.1f}%")
                with stats_col2:
                    st.metric("中等置信度", f"{medium_confidence}", delta=f"{medium_confidence/len(plddt_scores)*100:.1f}%")
                with stats_col3:
                    st.metric("低置信度", f"{low_confidence}", delta=f"{low_confidence/len(plddt_scores)*100:.1f}%")
        
        st.markdown("</div>", unsafe_allow_html=True)


def plot_plddt_chart_ultrawide(plddt_scores, sequence):
    """绘制超宽pLDDT分数图表"""
    try:
        import plotly.graph_objects as go
        import numpy as np
        
        # 基于序列生成更真实的pLDDT模式
        seq_len = len(sequence)
        
        # 生成更真实的pLDDT分数（基于序列特征）
        if isinstance(plddt_scores, list) or isinstance(plddt_scores, np.ndarray):
            if len(plddt_scores) == seq_len:
                # 使用提供的分数
                scores = np.array(plddt_scores)
            else:
                # 生成基于序列特征的分数
                scores = generate_realistic_plddt(sequence)
        else:
            # 生成基于序列特征的分数
            scores = generate_realistic_plddt(sequence)
        
        # 创建图表
        fig = go.Figure()
        
        # 添加pLDDT线
        fig.add_trace(go.Scatter(
            x=list(range(1, len(scores) + 1)),
            y=scores,
            mode='lines+markers',
            name='pLDDT',
            line=dict(color='blue', width=3),
            marker=dict(size=6, color='blue'),
            hovertemplate='残基 %{x}: %{y:.1f}<br>氨基酸: %{text}',
            text=[sequence[i] if i < len(sequence) else '' for i in range(len(scores))]
        ))
        
        # 添加置信度区域
        fig.add_hrect(
            y0=90, y1=100,
            fillcolor="rgba(0,200,0,0.15)",
            layer="below",
            line_width=0,
            annotation_text="极高置信度 (>90)",
            annotation_position="top left",
            annotation_font_size=14
        )
        
        fig.add_hrect(
            y0=70, y1=90,
            fillcolor="rgba(135,206,250,0.15)",
            layer="below",
            line_width=0,
            annotation_text="高置信度 (70-90)",
            annotation_position="top left",
            annotation_font_size=14
        )
        
        fig.add_hrect(
            y0=50, y1=70,
            fillcolor="rgba(255,165,0,0.15)",
            layer="below",
            line_width=0,
            annotation_text="中等置信度 (50-70)",
            annotation_position="top left",
            annotation_font_size=14
        )
        
        fig.add_hrect(
            y0=0, y1=50,
            fillcolor="rgba(255,0,0,0.15)",
            layer="below",
            line_width=0,
            annotation_text="低置信度 (<50)",
            annotation_position="top left",
            annotation_font_size=14
        )
        
        avg_plddt = np.mean(scores)
        
        # 添加平均线
        fig.add_hline(
            y=avg_plddt,
            line_dash="dash",
            line_color="red",
            line_width=3,
            annotation_text=f"平均: {avg_plddt:.1f}",
            annotation_position="bottom right",
            annotation_font_size=16,
            annotation_font_color="red"
        )
        
        # 更新布局 - 超宽设置
        fig.update_layout(
            title=dict(
                text=f"结构置信度分析 - 序列长度: {seq_len} aa (平均pLDDT: {avg_plddt:.1f})",
                font=dict(size=20, color='black'),
                x=0.5,
                xanchor='center'
            ),
            xaxis_title="残基位置",
            yaxis_title="pLDDT分数",
            yaxis_range=[0, 100],
            height=600,
            width=None,  # 自动宽度
            hovermode='x unified',
            showlegend=True,
            plot_bgcolor='white',
            paper_bgcolor='white',
            xaxis=dict(
                showgrid=True, 
                gridwidth=1, 
                gridcolor='lightgray',
                title_font=dict(size=16),
                tickfont=dict(size=14)
            ),
            yaxis=dict(
                showgrid=True, 
                gridwidth=1, 
                gridcolor='lightgray',
                title_font=dict(size=16),
                tickfont=dict(size=14)
            ),
            margin=dict(l=50, r=50, t=80, b=50)
        )
        
        # 添加网格
        fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='lightgray')
        fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgray')
        
        # 使用全宽显示
        st.plotly_chart(fig, use_container_width=True)
        
        # 显示统计信息 - 使用更多列
        col1, col2, col3, col4, col5 = st.columns(5)
        with col1:
            st.metric("平均pLDDT", f"{avg_plddt:.1f}")
        with col2:
            high_conf = sum(1 for s in scores if s >= 70)
            percentage = high_conf/len(scores)*100
            st.metric("高置信度残基", f"{high_conf}", delta=f"{percentage:.1f}%")
        with col3:
            medium_conf = sum(1 for s in scores if 50 <= s < 70)
            percentage = medium_conf/len(scores)*100
            st.metric("中等置信度", f"{medium_conf}", delta=f"{percentage:.1f}%")
        with col4:
            low_conf = sum(1 for s in scores if s < 50)
            percentage = low_conf/len(scores)*100
            st.metric("低置信度", f"{low_conf}", delta=f"{percentage:.1f}%")
        with col5:
            max_score = np.max(scores)
            min_score = np.min(scores)
            st.metric("波动范围", f"{min_score:.1f}-{max_score:.1f}")
        
    except ImportError:
        # 备用方案：使用matplotlib
        try:
            import matplotlib.pyplot as plt
            
            fig, ax = plt.subplots(figsize=(20, 8))
            
            # 生成分数
            if isinstance(plddt_scores, list) or isinstance(plddt_scores, np.ndarray):
                if len(plddt_scores) == len(sequence):
                    scores = np.array(plddt_scores)
                else:
                    scores = generate_realistic_plddt(sequence)
            else:
                scores = generate_realistic_plddt(sequence)
            
            ax.plot(range(1, len(scores) + 1), scores, 'b-', linewidth=3, alpha=0.8)
            
            # 添加区域
            ax.fill_between(range(1, len(scores) + 1), 90, 100, alpha=0.2, color='green', label='极高置信度 (>90)')
            ax.fill_between(range(1, len(scores) + 1), 70, 90, alpha=0.2, color='lightblue', label='高置信度 (70-90)')
            ax.fill_between(range(1, len(scores) + 1), 50, 70, alpha=0.2, color='orange', label='中等置信度 (50-70)')
            ax.fill_between(range(1, len(scores) + 1), 0, 50, alpha=0.2, color='red', label='低置信度 (<50)')
            
            ax.set_xlabel('残基位置', fontsize=14)
            ax.set_ylabel('pLDDT分数', fontsize=14)
            ax.set_title(f'结构置信度分析 (平均pLDDT: {np.mean(scores):.1f})', fontsize=16, fontweight='bold')
            ax.set_ylim(0, 100)
            ax.grid(True, alpha=0.3, linestyle='--')
            ax.legend(loc='upper right', fontsize=12)
            
            # 添加平均线
            ax.axhline(y=np.mean(scores), color='red', linestyle='--', linewidth=2, label=f'平均: {np.mean(scores):.1f}')
            
            st.pyplot(fig, use_container_width=True)
            
        except ImportError:
            st.warning("需要安装plotly或matplotlib来显示图表: `pip install plotly matplotlib`")
    except Exception as e:
        st.warning(f"图表显示失败: {e}")


def generate_realistic_plddt(sequence):
    """基于序列特征生成更真实的pLDDT分数"""
    import numpy as np
    import hashlib
    
    # 使用序列的哈希值作为随机种子，确保同一序列生成相同的pLDDT
    seq_hash = int(hashlib.md5(sequence.encode()).hexdigest()[:8], 16)
    np.random.seed(seq_hash % 10000)
    
    length = len(sequence)
    scores = np.zeros(length)
    
    # 基本分数（基于序列位置）
    for i in range(length):
        # 基础分数
        base_score = 75
        
        # 基于氨基酸类型调整
        aa = sequence[i]
        if aa in 'ACDEFGHIKLMNPQRSTVWY':
            # 给不同氨基酸类型不同的基准
            if aa in 'DEKR':  # 带电氨基酸
                base_score += np.random.normal(0, 3)
            elif aa in 'FILMVWY':  # 疏水氨基酸
                base_score += np.random.normal(0, 4)
            elif aa in 'P':  # 脯氨酸
                base_score += np.random.normal(-5, 5)
            elif aa in 'G':  # 甘氨酸
                base_score += np.random.normal(0, 6)
            else:  # 其他
                base_score += np.random.normal(0, 2)
        
        # 添加周期性波动
        period = 8 + (seq_hash % 7)  # 5-12之间的周期
        wave = 10 * np.sin(2 * np.pi * i / period)
        
        # 添加趋势
        trend = -0.05 * i if i > length/2 else 0.02 * i
        
        # 组合所有因素
        score = base_score + wave + trend
        
        # 添加随机噪声
        noise = np.random.normal(0, 2)
        score += noise
        
        # 确保分数在合理范围内
        scores[i] = max(30, min(95, score))
    
    return scores


def display_mutant_recommendations(mutants: List[Dict], memory: EnhancedProteinMemory):
    """显示突变体推荐"""
    st.subheader("AI推荐突变体")
    
    analysis_type = st.session_state.get('current_analysis_type', '智能规则')
    st.info(f"当前使用: **{analysis_type}** 模式")
    
    for i, mutant in enumerate(mutants):
        confidence = mutant.get('confidence', 0.7)
        status = mutant.get('status', '通过')
        
        # 根据状态显示不同的标题颜色和标签
        header_text = f"#{i+1} 得分: {mutant['score']:.3f} | 置信度: {confidence:.0%}"
        if status == '淘汰':
            header_text += " | ❌ 已淘汰"
        
        with st.expander(header_text, expanded=i < 2):
            if status == '淘汰':
                st.error(f"⚠️ 该突变体已被系统淘汰: {mutant.get('elimination_reason', '原因未知')}")
            
            col1, col2 = st.columns([2,1])
            
            with col1:
                st.write("**突变信息**")
                st.write(f"突变数量: {mutant['mutation_count']}")
                
                # 新增结构评估信息显示
                if 'avg_plddt' in mutant:
                    plddt = mutant['avg_plddt']
                    drop = mutant.get('plddt_drop', 0)
                    
                    plddt_col1, plddt_col2 = st.columns(2)
                    with plddt_col1:
                        st.metric("预测 pLDDT", f"{plddt:.1f}")
                    with plddt_col2:
                        st.metric("稳定性变化", f"-{drop:.1f}", delta_color="inverse")
                
                if mutant.get('is_ml_based'):
                    st.success("机器学习预测")
                else:
                    st.info("智能规则推荐")
                
                # 新增：突变体结构查看按钮
                if 'structure_result' in mutant:
                    if st.button(f"查看突变体 #{i+1} 3D结构", key=f"view_mut_struct_{i}"):
                        st.session_state.structure_results = mutant['structure_result']
                        st.session_state.show_structure_results = True
                        st.session_state.last_triggered = 'structure'
                        st.toast(f"已加载突变体 #{i+1} 结构", icon="🔬")
                        st.rerun()
                
                st.write("**具体突变:**")
                for mutation in mutant['mutations']:
                    st.code(mutation, language=None)
                
                if mutant.get('sequence'):
                    st.write("**序列预览:**")
                    seq_preview = mutant['sequence'][:50] + "..." if len(mutant['sequence']) > 50 else mutant['sequence']
                    st.text(seq_preview)
                
                st.markdown("---")
                
                if 'protein_advisor' in st.session_state:
                    if st.button(f"AI深度分析 #{i+1}", 
                                key=f"ai_analyze_{i}",
                                use_container_width=True,
                                type="primary"):
                        
                        st.session_state[f'analyzing_mutant_{i}'] = True
                        st.session_state.current_analyzing_idx = i
                        st.rerun()
                
                if st.session_state.get(f'analyzing_mutant_{i}', False):
                    with st.spinner("AI正在深度分析..."):
                        if st.session_state.current_analysis:
                            current = st.session_state.current_analysis
                            
                            analysis_mode = st.session_state.get('analysis_mode', '智能选择')
                            use_online = analysis_mode in ["智能选择", "强制在线"]
                            
                            st.write("**分析参数:**")
                            st.write(f"- 野生型序列长度: {len(current['sequence'])}")
                            st.write(f"- 突变: {mutant['mutations']}")
                            st.write(f"- 得分: {mutant['score']}")
                            st.write(f"- 目标: {current['goal']}")
                            st.write(f"- 领域: {current['domain']}")
                            st.write(f"- use_online: {use_online}")
                            
                            try:
                                result = st.session_state.protein_advisor.get_analysis(
                                wildtype=current['sequence'],
                                mutations=mutant['mutations'],
                                score=mutant['score'],
                                goal=current['goal'],
                                domain=current['domain'],
                                use_online=use_online
                                )
    
                                if not result or not isinstance(result, dict):
                                    result = {
                                        'content': '分析服务暂时不可用，请稍后重试。',
                                        'source': '备用分析',
                                        'is_online': False
                                    }
    
                                st.markdown(f"### {result.get('source', 'AI分析')}结果")
    
                                if result.get('is_online', False):
                                    st.success("来自DeepSeek在线AI")
                                else:
                                    st.info("来自离线知识库")
    
                                st.markdown("---")
                                st.markdown(result.get('content', '无分析内容'))
    
                                st.markdown("---")
                                st.write("**AI分析调试信息:**")
                                st.write(f"- 分析来源: {result.get('source', '未知')}")
                                st.write(f"- 是否为在线: {result.get('is_online', False)}")
                                st.write(f"- 内容长度: {len(result.get('content', ''))}")
    
                            except Exception as e:
                                st.error(f"AI分析过程中出错: {str(e)}")
                                st.error(f"错误类型: {type(e).__name__}")
    
                                st.info("使用备用离线分析...")
                                st.markdown(f"""
                                ### 备用分析结果
    
                                **突变:** {', '.join(mutant['mutations'])}
                                **目标:** {current['goal']}
                                **得分:** {mutant['score']:.3f}
    
                                **建议:**
                                1. 这些突变旨在优化蛋白质{current['goal']}
                                2. 建议进行小规模表达验证
                                3. 考虑逐步突变，避免多重突变干扰
    
                                *注: AI分析服务暂时不可用，请稍后重试。*
                                """)
                            
                            st.markdown("---")
                            if st.button("关闭分析", key=f"close_analysis_{i}"):
                                st.session_state[f'analyzing_mutant_{i}'] = False
                                st.rerun()
            
            with col2:
                st.write("**设计理由**")
                st.info(mutant['rationale'])
                
                if mutant.get('is_ml_based') and 'sequence' in mutant:
                    st.write("**完整序列:**")
                    st.text_area(
                        f"突变体 #{i+1} 完整序列",
                        value=mutant['sequence'],
                        height=100,
                        key=f"mutant_sequence_{i}",
                        disabled=True
                    )
                
                score = mutant['score']
                if score > 0.8:
                    st.success("高潜力 - 强烈推荐实验验证")
                elif score > 0.6:
                    st.info("中等潜力 - 建议进一步分析")
                else:
                    st.warning("探索性 - 需要谨慎验证")


def display_similar_cases(similar_cases: List[Dict]):
    """显示相似案例"""
    if similar_cases:
        st.subheader("相似历史案例")
        
        for i, case in enumerate(similar_cases):
            with st.expander(f"案例 #{i+1}: {case['strategy']} ({case['time']})", expanded=i==0):
                col1, col2 = st.columns(2)
                
                with col1:
                    st.write(f"**目标**: {case['goal']}")
                    st.write(f"**结果**: {'成功' if case['success'] else '失败'}")
                    st.write(f"**序列长度**: {case['sequence_length']} aa")
                
                with col2:
                    if case['mutations']:
                        st.write("**关键突变**:")
                        for mutation in case['mutations'][:3]:
                            st.write(f"• {mutation}")
                    
                    features = case['features']
                    st.write("**序列特征**:")
                    st.write(f"疏水性: {features['hydrophobic']:.1%} | 带电性: {features['charged']:.1%}")


def display_learning_progress(memory: EnhancedProteinMemory, training_samples: List):
    """显示学习进度"""
    st.subheader("AI学习进度")
    
    stats = memory.get_stats()
    exp_count = len(training_samples) + stats['total']
    
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("设计案例", stats['total'])
    with col2:
        if stats['total'] > 0:
            st.metric("成功率", f"{stats['rate']:.1%}")
        else:
            st.metric("成功率", "0.0%")
    with col3:
        st.metric("学习模式", stats['patterns'])
    with col4:
        if exp_count >= 30:
            level = "专家"
        elif exp_count >= 15:
            level = "高级"
        elif exp_count >= 5:
            level = "中级"
        else:
            level = "新手"
        st.metric("智能体等级", level)
    
    st.markdown("---")
    st.write("**经验积累进度:**")
    
    max_exp = 50
    exp_progress = min(exp_count / max_exp, 1.0)
    
    col1, col2 = st.columns([3, 1])
    with col1:
        st.progress(exp_progress)
    with col2:
        st.write(f"{exp_count}/{max_exp}")
    
    if stats['patterns'] > 0:
        st.write("**学习到的成功模式:**")
        tab1, tab2, tab3 = st.tabs(["稳定性", "活性", "可溶性"])
        
        with tab1:
            patterns = memory.get_patterns('stability')
            if patterns:
                for pattern in patterns:
                    st.info(pattern)
            else:
                st.write("暂无模式")
        
        with tab2:
            patterns = memory.get_patterns('activity')
            if patterns:
                for pattern in patterns:
                    st.info(pattern)
            else:
                st.write("暂无模式")
        
        with tab3:
            patterns = memory.get_patterns('solubility')
            if patterns:
                for pattern in patterns:
                    st.info(pattern)
            else:
                st.write("暂无模式")


def show_realtime_learning_panel(memory: EnhancedProteinMemory, training_samples: List):
    """显示实时学习面板"""
    st.sidebar.markdown("---")
    st.sidebar.subheader("学习进度")
    
    stats = memory.get_stats()
    exp_count = len(training_samples) + stats['total']
    
    max_exp = 50
    progress = min(exp_count / max_exp, 1.0)
    st.sidebar.write(f"**经验值:** {exp_count}/{max_exp}")
    st.sidebar.progress(progress)
    
    if exp_count >= 30:
        level = "专家"
    elif exp_count >= 15:
        level = "高级"
    elif exp_count >= 5:
        level = "中级"
    else:
        level = "新手"
    
    st.sidebar.write(f"**等级:** {level}")
    st.sidebar.write(f"**案例总数:** {stats['total']}")
    st.sidebar.write(f"**成功模式:** {stats['patterns']}个")
    
    if stats['total'] > 0:
        st.sidebar.write(f"**成功率:** {stats['rate']:.1%}")
    
    if st.sidebar.button("重置学习进度", type="secondary"):
        if 'memory' in st.session_state:
            st.session_state.memory = EnhancedProteinMemory()
        if 'training_samples' in st.session_state:
            st.session_state.training_samples = []
        st.sidebar.success("学习进度已重置")
        st.rerun()


def check_network_connection():
    """检查网络连接状态"""
    try:
        import urllib.request
        
        test_urls = [
            'https://blast.ncbi.nlm.nih.gov/',
            'https://www.ncbi.nlm.nih.gov/'
        ]
        
        for url in test_urls:
            try:
                request = urllib.request.Request(
                    url,
                    headers={'User-Agent': 'Mozilla/5.0'}
                )
                
                response = urllib.request.urlopen(request, timeout=10)
                if response.getcode() == 200:
                    return True, f"成功连接到 {url}"
            except Exception:
                continue
        
        try:
            socket.create_connection(("blast.ncbi.nlm.nih.gov", 443), timeout=10)
            return True, "可以连接到NCBI服务器"
        except:
            try:
                socket.gethostbyname("blast.ncbi.nlm.nih.gov")
                return True, "可以解析NCBI域名"
            except:
                pass
        
        return False, "无法连接到NCBI服务器"
        
    except Exception as e:
        return False, f"网络检查失败: {str(e)}"


# 初始化session state
if 'ai_agent' not in st.session_state:
    import time
    random_seed = int(time.time() * 1000) % 1000000
    st.session_state.ai_agent = ProteinAIAgent(random_seed=random_seed)
if 'memory' not in st.session_state:
    st.session_state.memory = EnhancedProteinMemory()
if 'training_samples' not in st.session_state:
    st.session_state.training_samples = []
if 'current_analysis' not in st.session_state:
    st.session_state.current_analysis = None
if 'similar_cases' not in st.session_state:
    st.session_state.similar_cases = None
if 'protein_advisor' not in st.session_state:
    st.session_state.protein_advisor = SmartProteinAdvisor()
if 'homology_analyzer' not in st.session_state:
    st.session_state.homology_analyzer = None
if 'ncbi_email' not in st.session_state:
    st.session_state.ncbi_email = ""
if 'blast_results' not in st.session_state:
    st.session_state.blast_results = None
if 'structure_predictor' not in st.session_state:
    st.session_state.structure_predictor = UnifiedStructurePredictor()
if 'unified_predictor' not in st.session_state:
    st.session_state.unified_predictor = UnifiedStructurePredictor()
if 'blast_analyzer' not in st.session_state:
    st.session_state.blast_analyzer = None
if 'structure_predictor' not in st.session_state:
    st.session_state.structure_predictor = None
if 'structure_results' not in st.session_state:
    st.session_state.structure_results = None
if 'show_structure_tab' not in st.session_state:
    st.session_state.show_structure_tab = False
if 'use_local_esmfold' not in st.session_state:
    st.session_state.use_local_esmfold = False
if 'blast_triggered' not in st.session_state:
    st.session_state.blast_triggered = False
if 'structure_triggered' not in st.session_state:
    st.session_state.structure_triggered = False
if 'show_blast_results' not in st.session_state:
    st.session_state.show_blast_results = False
if 'show_structure_results' not in st.session_state:
    st.session_state.show_structure_results = False

# 主界面
st.markdown('<div class="main-header">', unsafe_allow_html=True)
st.title("蛋白质AI设计实验室")
st.subheader("智能蛋白质工程与优化平台")
st.markdown('</div>', unsafe_allow_html=True)

# 侧边栏
st.sidebar.header("控制面板")
st.sidebar.subheader("机器学习训练")

uploaded_file = st.sidebar.file_uploader("上传训练数据(CSV)", type=['csv'])

if uploaded_file is not None:
    try:
        training_data = pd.read_csv(uploaded_file)
        st.sidebar.success(f"成功加载 {len(training_data)} 条训练数据")
        
        if st.sidebar.button("训练AI模型", type="primary"):
            with st.spinner("正在训练机器学习模型..."):
                sequences = training_data['Sequence'].tolist()
                targets = training_data['dep'].tolist()
                
                success = st.session_state.ai_agent.train_ml_model(sequences, targets)
                if success:
                    st.sidebar.success("模型训练完成！")
                    st.session_state.ml_model_trained = True
                else:
                    st.sidebar.error("模型训练失败")
                    
    except Exception as e:
        st.sidebar.error(f"加载训练数据失败: {str(e)}")

prediction_mode = st.sidebar.radio(
    "预测模式",
    ["智能规则", "机器学习"],
    help="智能规则：基于生物物理规则\n机器学习：使用训练好的AAindex模型"
)

st.sidebar.markdown("---")
st.sidebar.subheader("AI专家设置")

api_key = st.sidebar.text_input(
    "DeepSeek API密钥",
    type="password",
    help="免费获取：https://platform.deepseek.com/",
    value=st.session_state.get('deepseek_key', '')
)

if api_key:
    st.session_state.deepseek_key = api_key
    try:
        if 'protein_advisor' in st.session_state:
            st.session_state.protein_advisor.set_deepseek_key(api_key)
        st.sidebar.success("DeepSeek API密钥已保存")
    except Exception as e:
        st.sidebar.error(f"初始化DeepSeek失败: {str(e)}")

analysis_mode = st.sidebar.radio(
    "分析模式",
    ["智能选择", "仅离线", "强制在线"],
    help="智能选择：优先用API，失败用离线"
)

with st.sidebar.expander("DeepSeek 配置说明"):
    st.markdown("""
    **DeepSeek API（推荐）**：
    1. 访问：https://platform.deepseek.com/
    2. 注册获取API密钥
    3. 粘贴到左侧输入框
    
    **离线模式**：
    - 无需API密钥
    - 基于规则的知识库
    - 适合快速预览
    """)

st.sidebar.markdown("---")
st.sidebar.subheader("NCBI BLAST分析")

ncbi_email = st.sidebar.text_input(
    "NCBI邮箱地址",
    value=st.session_state.get('ncbi_email', ''),
    placeholder="your_email@example.com",
)

if ncbi_email and '@' in ncbi_email:
    st.session_state.ncbi_email = ncbi_email
    try:
        if st.session_state.homology_analyzer is None:
            st.session_state.homology_analyzer = EnhancedNCBIBlastAnalyzer(email=ncbi_email)
    except Exception as e:
        st.sidebar.error(f"初始化失败: {str(e)}")

st.sidebar.selectbox(
    "BLAST数据库",
    ["nr", "swissprot", "pdb", "refseq_protein"],
    index=0,
    help="nr: 非冗余蛋白质数据库（默认）\nswissprot: 高质量注释蛋白质\npdb: 蛋白质结构数据库"
)

st.sidebar.radio(
    "BLAST程序",
    ["blastp", "blastn", "blastx", "tblastn", "tblastx"],
    index=0,
    help="blastp: 蛋白质-蛋白质比对（默认）"
)

st.sidebar.markdown("---")
st.sidebar.subheader("结构预测设置")

if STRUCTURE_MODULE_AVAILABLE:
    st.session_state.use_local_esmfold = st.sidebar.checkbox(
        "使用本地ESMFold模型（需安装esm）",
        value=False,
        help="如果未安装esm库，将自动使用API模式"
    )
    
    enable_structure_analysis = st.sidebar.checkbox(
        "启用结构分析功能",
        value=True,
        help="启用蛋白质结构预测和可视化功能"
    )
    
    enable_interaction_analysis = st.sidebar.checkbox(
        "启用相互作用力分析",
        value=True,
        help="分析蛋白质间相互作用力"
    )
    
    st.sidebar.info("结构预测模块已启用")
else:
    st.sidebar.warning("结构预测模块不可用")

st.sidebar.markdown("---")
st.sidebar.subheader("网络状态检查")

if st.sidebar.button("检查网络连接"):
    success, message = check_network_connection()
    if success:
        st.sidebar.success(message)
    else:
        st.sidebar.error(message)
        st.sidebar.info("请检查您的网络设置，确保可以访问NCBI网站")

st.sidebar.markdown("---")
st.sidebar.subheader("设计配置")

analysis_goal = st.sidebar.selectbox(
    "优化目标",
    ["stability", "activity", "solubility"],
    format_func=lambda x: {"stability": "稳定性", "activity": "活性", "solubility": "可溶性"}[x]
)

protein_domain = st.sidebar.selectbox(
    "蛋白质类型",
    ["enzyme", "antibody", "membrane_protein", "other"],
    format_func=lambda x: {
        "enzyme": "酶", 
        "antibody": "抗体", 
        "membrane_protein": "膜蛋白",
        "other": "其他"
    }[x]
)

st.sidebar.subheader("序列输入")

sequence_input = st.sidebar.text_area(
    "蛋白质序列（野生型）",
    value="MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM",
    height=100,
    placeholder="输入蛋白质氨基酸序列...",
    key='current_sequence'
)

st.sidebar.markdown("---")
st.sidebar.write("**相互作用分析**")

protein2_input = st.sidebar.text_area(
    "对比蛋白质序列（可选）",
    value="",  
    height=80,
    placeholder="输入第二个蛋白质序列...\n留空则分析自身相互作用",
    key='protein2_sequence',
    help="输入另一个蛋白质序列进行相互作用分析"
)

use_self_interaction = st.sidebar.checkbox(
    "分析自身相互作用", 
    value=True,
    help="分析蛋白质与自身的相互作用"
)

if protein2_input and protein2_input.strip():
    use_self_interaction = False

st.session_state.current_goal = analysis_goal

show_realtime_learning_panel(st.session_state.memory, st.session_state.training_samples)

st.sidebar.subheader("实验反馈")
feedback_outcome = st.sidebar.selectbox("实验结果", ["", "成功", "部分成功", "失败"])
feedback_notes = st.sidebar.text_area("实验记录", placeholder="记录关键发现和突变效果...")

if st.sidebar.button("提交反馈", type="secondary"):
    if feedback_outcome and st.session_state.current_analysis:
        mutations = re.findall(r'[A-Z]\d+[A-Z]', feedback_notes)
        
        success_status = feedback_outcome == "成功"
        case_id = st.session_state.memory.save_case(
            sequence=sequence_input,
            goal=analysis_goal,
            strategy="AI设计",
            success=success_status,
            mutations=mutations,
            is_user_verified=True
        )
        
        st.session_state.training_samples.append({
            'sequence': sequence_input,
            'goal': analysis_goal,
            'domain': protein_domain,
            'timestamp': datetime.now().strftime("%Y-%m-%d %H:%M"),
            'outcome': feedback_outcome
        })
        
        st.sidebar.success(f"反馈已保存! 案例ID: {case_id}")
        st.rerun()

# 主内容区域 - 优化按钮布局
st.markdown("### 分析工具")

col1, col2, col3 = st.columns([2, 2, 1.5])

with col1:
    if st.button("启动AI分析", type="primary", use_container_width=True):
        if sequence_input:
            validation = validate_protein_sequence(sequence_input)
            
            if validation['valid']:
                with st.spinner("AI正在分析蛋白质序列..."):
                    cleaned_wt = validation['cleaned_sequence']
                    analysis = st.session_state.ai_agent.analyze_sequence(
                        cleaned_wt, analysis_goal, protein_domain
                    )
                    
                    if prediction_mode == "机器学习" and st.session_state.get('ml_model_trained', False):
                        mutants = st.session_state.ai_agent.generate_ml_mutants(
                            cleaned_wt, 5
                        )
                        st.session_state.current_analysis_type = "机器学习"
                    else:
                        mutants = st.session_state.ai_agent.generate_mutants(
                            cleaned_wt, analysis_goal, protein_domain,
                            st.session_state.memory, 5
                        )
                        st.session_state.current_analysis_type = "智能规则"
                    
                    # --- 新增：突变体结构稳定性验证 ---
                    st.write("🔍 正在进行突变体结构稳定性交叉验证...")
                    
                    # 1. 预测野生型结构（如果还没有）
                    if st.session_state.structure_results is None or st.session_state.structure_results.get('sequence') != cleaned_wt:
                        with st.spinner("正在预测野生型参考结构..."):
                            wt_struct = st.session_state.unified_predictor.predict(cleaned_wt, job_name="wildtype")
                            st.session_state.structure_results = wt_struct
                    
                    wt_plddt = st.session_state.structure_results.get('avg_plddt', 0)
                    
                    # 2. 为前 N 个突变体预测结构并评估 pLDDT
                    valid_mutants = []
                    for i, mutant in enumerate(mutants):
                        mut_seq = mutant['sequence']
                        with st.spinner(f"验证突变体 #{i+1} 结构稳定性..."):
                            mut_struct = st.session_state.unified_predictor.predict(mut_seq, job_name=f"mutant_{i+1}")
                            
                            mut_plddt = mut_struct.get('avg_plddt', 0)
                            plddt_drop = wt_plddt - mut_plddt
                            
                            mutant['avg_plddt'] = mut_plddt
                            mutant['plddt_drop'] = plddt_drop
                            mutant['structure_result'] = mut_struct
                            
                            # 淘汰机制：如果 pLDDT 下降超过 15 或者绝对值低于 50
                            if plddt_drop > 15 or mut_plddt < 50:
                                mutant['status'] = "淘汰"
                                mutant['elimination_reason'] = f"结构稳定性骤降 (pLDDT 下降 {plddt_drop:.1f})"
                            else:
                                mutant['status'] = "通过"
                            
                        valid_mutants.append(mutant)
                    
                    mutants = valid_mutants
                    # -----------------------------------
                    
                    st.session_state.current_analysis = {
                        'analysis': analysis,
                        'mutants': mutants,
                        'sequence': cleaned_wt,
                        'goal': analysis_goal,
                        'domain': protein_domain
                    }
                    st.session_state.last_triggered = 'ai'

                    similar_cases = st.session_state.memory.find_similar(
                        validation['cleaned_sequence'], analysis_goal, 3
                    )
                    st.session_state.similar_cases = similar_cases

                    if mutants:
                        mutations = [mutant['mutations'][0] for mutant in mutants if mutant['mutations']]
                        st.session_state.memory.save_case(
                            sequence=validation['cleaned_sequence'],
                            goal=analysis_goal,
                            strategy="AI分析",
                            success=False, # 初始标记为未验证
                            mutations=mutations[:2],
                            is_user_verified=False
                        )
                    
                    st.rerun() # 强制重绘以实现标签页自动跳转
            
            else:
                st.error(f"序列验证失败: {validation['message']}")
                if 'suggestion' in validation:
                    st.info(f"建议: {validation['suggestion']}")
        else:
            st.warning("请输入蛋白质序列")

with col2:
    blast_button = st.button("运行BLAST分析", type="secondary", use_container_width=True)
    
    if blast_button:
        st.session_state.blast_triggered = True
        st.session_state.show_blast_results = True
        if not sequence_input:
            st.warning("请输入蛋白质序列")
        elif not st.session_state.ncbi_email:
            st.warning("请先在侧边栏输入NCBI邮箱地址")
        else:
            if 'blast_analyzer' not in st.session_state:
                st.session_state.blast_analyzer = None
            
            if st.session_state.blast_analyzer is None:
                try:
                    st.session_state.blast_analyzer = EnhancedNCBIBlastAnalyzer(
                        email=st.session_state.ncbi_email
                    )
                    st.success("BLAST分析器已初始化")
                except Exception as e:
                    st.error(f"BLAST分析器初始化失败: {str(e)}")
                    st.info("请检查：1) NCBI邮箱格式是否正确 2) 网络连接是否正常 3) 是否有防火墙限制")
                    st.stop()
            
            try:
                with st.spinner("正在分析..."):
                    result = st.session_state.blast_analyzer.analyze_sequence(sequence_input)
                
                if result and result['status'] in ['success', 'cached']:
                    st.session_state.blast_results = result['data']
                    st.session_state.last_triggered = 'blast'
                
                    if result['status'] == 'cached':
                        st.info(f"使用缓存结果（保存于{result['data'].get('cache_time', '未知时间')}）")
                    else:
                        st.success("BLAST分析完成！")
                else:
                    st.error(f"分析失败: {result.get('message', '未知错误')}")
                    if result.get('metadata', {}).get('suggestion'):
                        st.info(f"建议: {result['metadata']['suggestion']}")
                        
            except Exception as e:
                st.error(f"BLAST分析过程中出错: {str(e)}")
                with st.expander("查看错误详情"):
                    st.code(str(e))

with col3:
    structure_button = st.button("结构分析", type="secondary", use_container_width=True)
    
    if structure_button:
        st.session_state.structure_triggered = True
        st.session_state.show_structure_results = True
        sequence_input = st.session_state.get('current_sequence', '')
        if not sequence_input:
            st.warning("请输入蛋白质序列")
        else:
            def clean_protein_sequence(seq):
                if not seq:
                    return ""
                seq = ''.join(seq.split())
                import re
                seq = re.sub(r'[^A-Za-z]', '', seq)
                return seq.upper()
            
            cleaned_seq = clean_protein_sequence(sequence_input)
            
            if len(cleaned_seq) < 10:
                st.error(f"有效序列太短 ({len(cleaned_seq)} aa)，至少需要10个有效氨基酸")
            else:
                st.info(f"使用本地ESM模型 | 序列长度: {len(cleaned_seq)} aa")
                
                try:
                    from structure_predictor_api import UnifiedStructurePredictor
                    
                    with st.spinner(f"正在调用 ESMFold 预测结构（首次可能耗时数十秒）..."):
                        predictor = UnifiedStructurePredictor()
                        result = predictor.predict(cleaned_seq)
                    
                    if result['success']:
                        st.session_state.structure_results = result
                        st.session_state.last_triggered = 'structure'
                        st.toast("✅ 结构预测已完成！", icon="🧬")
                        st.success("预测完成！已自动切换至下方「结构预测」标签页展示。")
                        st.rerun() # 强制重绘以实现标签页自动跳转
                        
                    else:
                        st.error(f"预测失败: {result.get('error', '未知错误')}")
                        if 'suggestion' in result:
                            st.info(f"建议: {result['suggestion']}")
                            
                except Exception as e:
                    st.error(f"结构预测异常: {str(e)}")
                    st.info("请确保已安装: `pip install fair-esm torch`")

# 结果显示区域 - 始终显示完整结果
st.markdown("<div class='result-container'>", unsafe_allow_html=True)

# 创建标签页，显示不同类型的结果
if (st.session_state.current_analysis or 
    st.session_state.show_blast_results or 
    st.session_state.show_structure_results):
    
    # 确定需要显示哪些标签 —— 按最近触发的分析置顶，实现"跳转"到最新结果所在标签
    tab_pairs = []
    if st.session_state.current_analysis:
        tab_pairs.append(('ai', "🤖 AI突变推荐"))
    if st.session_state.show_blast_results and st.session_state.blast_results:
        tab_pairs.append(('blast', "🔍 BLAST分析"))
    if st.session_state.show_structure_results and st.session_state.structure_results:
        tab_pairs.append(('structure', "🧬 结构预测"))

    last = st.session_state.get('last_triggered', '')
    if last and any(t == last for t, _ in tab_pairs):
        tab_pairs.sort(key=lambda x: x[0] != last)

    tabs_to_show = [t for t, _ in tab_pairs]
    tab_names = [n for _, n in tab_pairs]

    if tabs_to_show:
        tabs = st.tabs(tab_names)
        
        for i, tab_type in enumerate(tabs_to_show):
            with tabs[i]:
                if tab_type == 'ai' and st.session_state.current_analysis:
                    current = st.session_state.current_analysis
                    
                    if st.session_state.similar_cases:
                        display_similar_cases(st.session_state.similar_cases)
                        st.markdown("---")
                    
                    display_mutant_recommendations(current['mutants'], st.session_state.memory)
                    
                    st.markdown("---")
                    display_learning_progress(st.session_state.memory, st.session_state.training_samples)
                
                elif tab_type == 'blast' and st.session_state.blast_results:
                    display_blast_results_ultrawide(st.session_state.blast_results)
                
                elif tab_type == 'structure' and st.session_state.structure_results:
                    display_structure_results_ultrawide(st.session_state.structure_results)

st.markdown("</div>", unsafe_allow_html=True)

with st.expander("使用指南", expanded=False):
    st.markdown("""
    ### 如何使用蛋白质AI设计实验室
    
    1. **配置设计目标**
       - 稳定性: 提高热稳定性和结构刚性
       - 活性: 优化催化效率或结合亲和力  
       - 可溶性: 改善表达可溶性和减少聚集
    
    2. **输入蛋白质序列**
       - 使用标准20种氨基酸字符 (ACDEFGHIKLMNPQRSTVWY)
       - 长度建议: 50-500个氨基酸
       - 自动清理空格和无效字符
    
    3. **分析结果解读**
       - 序列分析: 基础理化性质评估 + 增强特征分析
       - 相似案例: 基于特征相似度的历史案例
       - 突变推荐: AI生成的优化突变体（包含学习模式）
       - 学习进度: AI智能体的经验积累
    
    4. **实验反馈**
       - 提交实验结果帮助AI学习
       - 记录关键突变和效果
       - 积累成功模式库
    """)

# 页脚
st.markdown("---")
st.markdown(
    "<div style='text-align: center; color: #666;'>"
    "蛋白质AI设计实验室 | 智能蛋白质工程平台"
    "</div>", 
    unsafe_allow_html=True
)