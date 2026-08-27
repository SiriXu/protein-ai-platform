from typing import List, Dict, Any, Optional
import json

class OfflineProteinKnowledge:
    """离线蛋白质知识库"""
    
    def __init__(self):
        self.knowledge = self._load_knowledge()
    
    def _load_knowledge(self):
        """加载知识库"""
        return {
            "mutation_rules": {
                "稳定性优化": [
                    "A→L/V/I: 增加疏水性，核心区域适用",
                    "D/E→N/Q: 消除表面电荷排斥",
                    "G→P: 限制构象柔性（避免活性位点）",
                    "表面K/R→Q/N: 减少表面电荷簇"
                ],
                "活性优化": [
                    "W/Y→A/G: 减小底物通道空间位阻",
                    "H→D/E: 调整催化残基pKa",
                    "表面D/E→K/R: 增强底物结合静电相互作用"
                ],
                "可溶性优化": [
                    "核心F/W→S/T: 核心疏水残基亲水化",
                    "表面K/R→E/D: 平衡表面电荷",
                    "引入N-X-S/T: N-糖基化位点"
                ]
            },
            "experiment_tips": [
                "热稳定性：差示扫描量热法(DSC)或热位移分析",
                "活性测定：使用天然底物，设阳性对照",
                "可溶性：动态光散射(DLS)检测聚集"
            ]
        }
    
    def explain(self, mutations: List[str], goal: str) -> str:
        """解释突变"""
        explanations = []
        goal_cn = {"stability": "稳定性", "activity": "活性", "solubility": "可溶性"}.get(goal, goal)
        
        for mutation in mutations:
            if len(mutation) >= 3:
                orig, new = mutation[0], mutation[-1]
                change = f"{orig}→{new}"
                
                # 根据目标提供解释
                if goal == "stability":
                    if orig in "ADE" and new in "NQ":
                        explanations.append(f"**{mutation}**：消除表面电荷，减少静电排斥，提高热稳定性")
                    elif orig in "FILMV" and new in "LIV":
                        explanations.append(f"**{mutation}**：优化疏水核心包装，增强van der Waals相互作用")
                
                elif goal == "activity":
                    if orig == "W" and new in "AG":
                        explanations.append(f"**{mutation}**：减小色氨酸的大侧链体积，可能改善底物通道")
                
                elif goal == "solubility":
                    if orig in "FILMWY" and new in "STNQ":
                        explanations.append(f"**{mutation}**：疏水残基亲水化，改善可溶性")
        
        if explanations:
            return "\n\n".join(explanations)
        return f"这些突变旨在优化蛋白质{goal_cn}，建议实验验证具体效果。"