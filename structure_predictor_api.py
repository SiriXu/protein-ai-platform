# structure_predictor_api.py
import numpy as np
from typing import Dict, Any, Optional
import hashlib
import math

class LocalESMPredictor:
    """本地ESM预测器 - 优化版本"""
    
    def __init__(self):
        self.model = None
        self.alphabet = None
        self.is_loaded = False
    
    def load_model(self):
        """加载本地ESM模型"""
        try:
            import esm
            import torch
            
            print("正在加载本地ESM2模型...")
            
            # 加载小型ESM2模型
            self.model, self.alphabet = esm.pretrained.esm2_t6_8M_UR50D()
            self.model = self.model.eval()
            
            # 使用GPU如果可用
            if torch.cuda.is_available():
                self.model = self.model.cuda()
                print("使用GPU加速")
            else:
                print("使用CPU模式")
            
            self.is_loaded = True
            print("本地ESM2模型加载成功")
            return True
            
        except ImportError as e:
            print(f"未安装fair-esm: {str(e)}")
            return False
        except Exception as e:
            print(f"模型加载失败: {str(e)}")
            return False
    
    def predict(self, sequence: str, job_name: str = "protein") -> Dict[str, Any]:
        """使用本地ESM模型预测 - 改进版本"""
        if not self.is_loaded:
            if not self.load_model():
                return {
                    'success': False, 
                    'error': '无法加载ESM模型',
                    'suggestion': '请安装fair-esm: pip install fair-esm'
                }
        
        try:
            import torch
            
            print(f"正在预测结构 ({len(sequence)} aa)...")
            
            # 获取ESM嵌入
            batch_converter = self.alphabet.get_batch_converter()
            data = [(job_name, sequence)]
            _, _, batch_tokens = batch_converter(data)
            
            if torch.cuda.is_available():
                batch_tokens = batch_tokens.cuda()
            
            with torch.no_grad():
                print("计算蛋白质嵌入...")
                results = self.model(batch_tokens, repr_layers=[6])
                embeddings = results["representations"][6][0].cpu().numpy()
            
            print("生成PDB结构...")
            
            # 生成更真实的pLDDT分数
            plddt_scores = self._generate_plddt_scores(sequence, embeddings)
            avg_plddt = np.mean(plddt_scores)
            
            # 生成PDB结构
            pdb_content = self._generate_pdb_structure(sequence, embeddings, plddt_scores)
            
            return {
                'success': True,
                'pdb_content': pdb_content,
                'plddt_scores': plddt_scores,
                'avg_plddt': float(avg_plddt),
                'sequence': sequence,
                'method': 'local_esm2',
                'api_source': '本地ESM2',
                'is_real_prediction': False,
                'note': f'基于ESM2嵌入生成的结构示意图（非真实折叠预测），序列长度: {len(sequence)} aa'
            }
            
        except Exception as e:
            return {
                'success': False, 
                'error': f'预测失败: {str(e)}',
                'suggestion': '请检查模型是否正确加载'
            }
    
    def _generate_plddt_scores(self, sequence: str, embeddings: np.ndarray) -> np.ndarray:
        """生成基于序列和嵌入的pLDDT分数"""
        length = len(sequence)
        scores = np.zeros(length)
        
        # 使用序列哈希确保可重复性
        seq_hash = int(hashlib.md5(sequence.encode()).hexdigest()[:8], 16)
        rng = np.random.RandomState(seq_hash % 10000)
        
        for i in range(length):
            if i < len(embeddings):
                # 基于嵌入的复杂性计算分数
                emb = embeddings[i]
                emb_norm = np.linalg.norm(emb)
                complexity = emb_norm / 100  # 归一化
                
                # 基础分数
                base_score = 70.0
                
                # 基于氨基酸类型的调整
                aa = sequence[i]
                if aa in 'DEKR':  # 带电氨基酸通常更难预测
                    base_score += rng.normal(-3, 2)
                elif aa in 'FILMVWY':  # 疏水氨基酸
                    base_score += rng.normal(2, 3)
                elif aa in 'P':  # 脯氨酸
                    base_score += rng.normal(-5, 4)
                elif aa in 'G':  # 甘氨酸
                    base_score += rng.normal(-2, 5)
                else:
                    base_score += rng.normal(0, 2)
                
                # 嵌入复杂性影响
                base_score += 5 * complexity
                
                # 位置效应
                if i < 10 or i > length - 10:  # 末端通常更难预测
                    base_score -= 5
                
                # 周期性波动
                period = 7 + (seq_hash % 6)
                wave = 8 * math.sin(2 * math.pi * i / period)
                
                # 趋势
                if i > length * 0.7:  # C端
                    trend = -0.03 * (i - length * 0.7)
                else:
                    trend = 0
                
                score = base_score + wave + trend
                
                # 添加随机噪声
                noise = rng.normal(0, 1.5)
                score += noise
                
                scores[i] = max(30, min(95, score))
            else:
                # 对于超出嵌入范围的位置使用默认值
                scores[i] = 65.0
        
        return scores
    
    def _generate_pdb_structure(self, sequence: str, embeddings: np.ndarray, plddt_scores: np.ndarray) -> str:
        """生成PDB结构文件"""
        pdb_lines = [
            "HEADER    ESM2 LOCAL PREDICTION",
            f"TITLE     Protein Structure - Generated by ESM2",
            "REMARK    Model: esm2_t6_8M_UR50D",
            "REMARK    Method: Embedding-based structure generation",
            "REMARK    Note: This is a simulated structure for visualization",
            "REMARK    "
        ]
        
        chain_id = 'A'
        atom_num = 1
        residue_num = 1
        
        # 生成螺旋-转角-螺旋的二级结构
        for i, aa in enumerate(sequence):
            if i >= len(embeddings):
                break
            
            # 基于嵌入生成坐标
            emb = embeddings[i]
            
            # 简化的二级结构模拟
            if i % 20 < 10:  # 螺旋区域
                radius = 4.5
                rise = 1.5
            else:  # 转角区域
                radius = 6.0
                rise = 0.8
            
            angle = 2 * math.pi * i / 3.6
            x = radius * math.cos(angle)
            y = radius * math.sin(angle)
            z = i * rise
            
            # 添加一些基于嵌入的扰动
            if len(emb) >= 3:
                x += 0.3 * (emb[0] / 10 - 0.5)
                y += 0.3 * (emb[1] / 10 - 0.5)
                z += 0.2 * (emb[2] / 10 - 0.5)
            
            plddt = plddt_scores[i] if i < len(plddt_scores) else 70.0
            
            # 添加主链原子
            atoms = [
                ("N", x, y, z, plddt),
                ("CA", x+1.5, y, z, plddt),
                ("C", x+1.5, y+1.2, z, plddt),
                ("O", x+0.8, y+1.8, z, plddt)
            ]
            
            for atom_name, x_, y_, z_, bfactor in atoms:
                line = f"ATOM  {atom_num:5d} {atom_name:^4s} {aa:3s} {chain_id}{residue_num:4d}    {x_:8.3f}{y_:8.3f}{z_:8.3f}  1.00{bfactor:6.2f}           {atom_name[0]}"
                pdb_lines.append(line)
                atom_num += 1
            
            # 添加CB原子（除甘氨酸外）
            if aa != 'G':
                cb_x = x + 2.0
                cb_y = y - 0.5
                cb_z = z
                line = f"ATOM  {atom_num:5d}  CB  {aa:3s} {chain_id}{residue_num:4d}    {cb_x:8.3f}{cb_y:8.3f}{cb_z:8.3f}  1.00{plddt:6.2f}           C"
                pdb_lines.append(line)
                atom_num += 1
            
            residue_num += 1
        
        # 添加连接信息
        for i in range(1, len(sequence)):
            pdb_lines.append(f"CONECT {i*5-4:5d} {i*5:5d}")
        
        pdb_lines.append("TER")
        pdb_lines.append("END")
        
        return "\n".join(pdb_lines)

class EsmFoldAPIPredictor:
    """ESMFold 在线 API 预测器 —— 调用 ESM Atlas 获取真实折叠结构"""

    API_URL = "https://api.esmatlas.com/v1/prediction/pdb"

    def predict(self, sequence: str, job_name: str = "protein", timeout: int = 120) -> Dict[str, Any]:
        """调用 ESMFold 在线 API 预测真实结构"""
        try:
            import requests
        except ImportError:
            return {
                'success': False,
                'error': '未安装 requests',
                'suggestion': 'pip install requests'
            }

        try:
            resp = requests.post(
                self.API_URL,
                data=sequence.encode('utf-8'),
                headers={"Content-Type": "text/plain"},
                timeout=timeout
            )
            if resp.status_code != 200:
                return {
                    'success': False,
                    'error': f'ESMFold API 返回 {resp.status_code}: {resp.text[:200]}',
                    'suggestion': 'API 可能限流或暂不可用，将回退到本地示意图'
                }

            pdb_content = resp.text
            plddt_scores = self._parse_plddt_from_pdb(pdb_content, len(sequence))
            avg_plddt = float(np.mean(plddt_scores)) if len(plddt_scores) > 0 else 0.0

            return {
                'success': True,
                'pdb_content': pdb_content,
                'plddt_scores': plddt_scores,
                'avg_plddt': avg_plddt,
                'sequence': sequence,
                'method': 'esmfold_api',
                'api_source': 'ESMFold 在线 API（真实预测）',
                'is_real_prediction': True,
                'note': f'ESMFold 真实结构预测，序列长度: {len(sequence)} aa，平均 pLDDT: {avg_plddt:.1f}'
            }
        except Exception as e:
            return {
                'success': False,
                'error': f'ESMFold API 调用失败: {str(e)}',
                'suggestion': '请检查网络连接，将回退到本地示意图'
            }

    @staticmethod
    def _parse_plddt_from_pdb(pdb_content: str, seq_len: int) -> np.ndarray:
        """从 PDB 的 B-factor 列解析每残基 pLDDT（取 CA 原子）"""
        scores = []
        seen = set()
        for line in pdb_content.splitlines():
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    resseq = int(line[22:26].strip())
                    bfactor = float(line[60:66].strip())
                    if resseq not in seen:
                        scores.append(bfactor)
                        seen.add(resseq)
                except (ValueError, IndexError):
                    continue
        if not scores:
            return np.full(max(seq_len, 1), 70.0)
        arr = np.array(scores, dtype=float)
        if len(arr) < seq_len:
            arr = np.concatenate([arr, np.full(seq_len - len(arr), arr[-1])])
        return arr[:seq_len]


class UnifiedStructurePredictor:
    """统一结构预测器 —— 优先 ESMFold 在线 API，失败回退本地示意图"""

    def __init__(self):
        self.api_predictor = EsmFoldAPIPredictor()
        self.local_predictor = LocalESMPredictor()

    def predict(self, sequence: str, priority: str = "api_first",
                job_name: str = "protein", show_info: bool = True) -> Dict[str, Any]:
        """
        预测结构：优先调用 ESMFold 在线 API 获取真实结构；
        API 不可用时回退到本地 ESM2 嵌入生成的结构示意图（明确标注非真实预测）。
        """
        api_result = self.api_predictor.predict(sequence, job_name)
        if api_result.get('success'):
            return api_result

        local_result = self.local_predictor.predict(sequence, job_name)
        local_result['fallback_reason'] = api_result.get('error', 'ESMFold API 不可用')
        return local_result