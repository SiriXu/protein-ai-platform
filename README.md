# 🧬 Protein AI Design Lab (蛋白质智能设计实验室)

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://protein-ai-design-lab.streamlit.app/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 🌟 项目简介
这是一个基于人工智能的蛋白质设计与优化平台，旨在辅助生物科研人员高效地进行蛋白质突变体筛选、同源性分析及三维结构验证。

本项目集成了：
- **AAindex + FFT 编码**：从理化性质和频域特征维度提取蛋白质序列信息。
- **GBR 机器学习模型**：实现突变体效应的定量打分。
- **ESMFold (ESM-2)**：端到端的蛋白质 3D 结构预测与稳定性验证。
- **DeepSeek AI 专家**：引入大语言模型为突变建议提供生物学解释与风险评估。

## 🚀 核心功能
1. **智能突变设计**：基于规则或机器学习模型推荐高潜力的蛋白质突变体。
2. **结构稳定性筛选**：自动识别并过滤可能导致结构崩塌（pLDDT 骤降）的危险突变。
3. **同源性分析 (BLAST)**：集成 NCBI BLAST API，从进化保守性维度辅助决策。
4. **3D 交互可视化**：实时预览突变体与野生型的结构对比。
5. **经验记忆系统**：支持湿实验数据反馈，构建“设计-验证-学习”的研发闭环。

## 🛠️ 安装指南

### 环境要求
- Python 3.8+
- NVIDIA GPU (可选，用于加速结构预测)

### 安装步骤
1. 克隆仓库：
   ```bash
   git clone https://github.com/your-username/protein-prediction.git
   cd protein-prediction
   ```
2. 安装依赖：
   ```bash
   pip install -r requirements.txt
   ```

## 📖 使用说明
1. 启动应用：
   ```bash
   streamlit run app.py
   ```
2. 在侧边栏输入 **DeepSeek API Key** 和 **NCBI Email**。
3. 输入蛋白质野生型序列，点击“启动 AI 分析”。
4. 在“结构预测”标签页查看三维模型，在“AI 突变推荐”中获取详细建议。

## 📂 项目结构
- `app.py`: Web 界面与主逻辑入口。
- `AAindex.py`: 核心特征工程模块。
- `structure_predictor_api.py`: 结构预测引擎。
- `deepseek_expert.py`: AI 专家交互模块。
- `data/`: 存放 AAindex 编码表等基础数据。

## 🤝 贡献指南
欢迎提交 Issue 或 Pull Request 来完善本项目！

## 📄 开源协议
本项目采用 [MIT License](LICENSE) 开源协议。
