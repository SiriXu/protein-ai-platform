# deepseek_expert.py
import sys
import os
from openai import OpenAI
from typing import List, Dict
import time

class DeepSeekProteinExpert:
    """DeepSeek蛋白质AI专家"""
    
    def __init__(self, api_key: str = None):
        print(f"\n{'='*50}")
        print("初始化DeepSeekProteinExpert...")
        
        # 修复密钥格式
        self.api_key = self._fix_key_format(api_key)
        print(f"API密钥: {self.api_key[:15]}...")
        
        self.is_ready = False
        self.client = None
        
        if self.api_key and len(self.api_key) > 30:
            try:
                # 创建客户端
                self.client = OpenAI(
                    api_key=self.api_key,
                    base_url="https://api.deepseek.com",  # 使用不带/v1的端点
                    timeout=30.0
                )
                
                # 快速测试连接
                self.is_ready = self._quick_test()
                
                if self.is_ready:
                    print("DeepSeek API连接成功！")
                else:
                    print("DeepSeek API连接测试失败")
                    
            except Exception as e:
                print(f"DeepSeek初始化异常: {e}")
                self.is_ready = False
        else:
            print("API密钥无效或太短")
            self.is_ready = False
        
        print(f"初始化完成: is_ready={self.is_ready}")
        print(f"{'='*50}\n")
    
    def _fix_key_format(self, key):
        """修复API密钥格式"""
        if not key:
            return key
        
        key = key.strip()
        
        # 确保以'sk-'开头
        if not key.startswith('sk-'):
            if len(key) > 30:  # 可能是忘了加前缀
                key = 'sk-' + key
            else:
                return key
        
        return key
    
    def _quick_test(self):
        """快速连接测试"""
        try:
            # 发送一个非常简单的请求
            response = self.client.chat.completions.create(
                model="deepseek-chat",
                messages=[{"role": "user", "content": "测试连接，请回复'OK'"}],
                max_tokens=5,
                temperature=0.1
            )
            result = response.choices[0].message.content.strip()
            print(f"测试响应: '{result}'")
            return True
        except Exception as e:
            print(f"连接测试失败: {e}")
            return False
    
    def analyze(self, wildtype: str, mutations: List[str], 
               score: float, goal: str, domain: str) -> Dict:
        """分析蛋白质突变 - 增强版，带重试机制"""
        print(f"\n{'='*50}")
        print("开始DeepSeek分析...")
        
        if not self.is_ready or not self.client:
            print("DeepSeek未就绪，返回离线分析")
            return {
                'content': "DeepSeek API未就绪，请检查API密钥和网络连接。",
                'source': 'DeepSeek离线',
                'is_online': False
            }
        
        # 重试机制
        max_retries = 3
        retry_delay = 5  # 秒
        
        for attempt in range(max_retries):
            try:
                print(f"第 {attempt + 1}/{max_retries} 次尝试...")
                
                # 如果这不是第一次尝试，等待一段时间
                if attempt > 0:
                    print(f"等待 {retry_delay} 秒后重试...")
                    time.sleep(retry_delay)
                    retry_delay = min(retry_delay * 1.5, 15)  # 指数退避，最大15秒
                
                # 构建prompt
                prompt = self._build_analysis_prompt(wildtype, mutations, score, goal, domain)
                
                print("发送请求到DeepSeek API...")
                start_time = time.time()
                
                # 使用更短的超时时间
                timeout_seconds = 30  # 默认30秒超时
                
                try:
                    # 发送请求
                    response = self.client.chat.completions.create(
                        model="deepseek-chat",
                        messages=[
                            {"role": "system", "content": "你是专业的蛋白质工程专家，请提供专业、实用、可操作的分析建议。用中文回答。"},
                            {"role": "user", "content": prompt}
                        ],
                        temperature=0.7,
                        max_tokens=1000,  # 减少token数量以加快响应
                        timeout=timeout_seconds
                    )
                    
                    elapsed_time = time.time() - start_time
                    content = response.choices[0].message.content
                    print(f"收到DeepSeek响应！长度: {len(content)}字符，耗时: {elapsed_time:.2f}秒")
                    
                    return {
                        'content': content,
                        'source': 'DeepSeek AI专家',
                        'is_online': True,
                        'response_time': f"{elapsed_time:.2f}秒"
                    }
                    
                except Exception as api_error:
                    error_msg = str(api_error)
                    print(f"API请求错误: {error_msg}")
                    
                    # 如果是超时错误，继续重试
                    if "timeout" in error_msg.lower() or "timed out" in error_msg.lower():
                        if attempt < max_retries - 1:
                            continue  # 继续重试
                        else:
                            raise TimeoutError(f"API请求超时，已重试{max_retries}次")
                    else:
                        # 其他错误直接抛出
                        raise api_error
                        
            except TimeoutError as te:
                print(f"DeepSeek分析超时: {te}")
                # 最后一次尝试也超时了
                if attempt == max_retries - 1:
                    return {
                        'content': "DeepSeek API请求超时，服务器响应较慢。建议稍后重试或使用离线分析。",
                        'source': 'DeepSeek超时',
                        'is_online': False,
                        'suggestion': '尝试：1. 稍后重试 2. 检查网络连接 3. 使用离线分析'
                    }
            except Exception as e:
                print(f"DeepSeek分析失败: {e}")
                # 如果不是超时错误，直接返回错误信息
                if attempt == max_retries - 1:
                    return {
                        'content': f"DeepSeek在线分析失败: {str(e)[:100]}",
                        'source': 'DeepSeek API错误',
                        'is_online': False,
                        'error_type': type(e).__name__
                    }
        
        # 所有重试都失败
        return {
            'content': "DeepSeek API请求失败，请检查网络连接和API密钥。",
            'source': 'DeepSeek连接失败',
            'is_online': False,
            'suggestion': '建议：1. 检查API密钥 2. 检查网络连接 3. 使用离线分析'
        }
    
    def _build_analysis_prompt(self, wildtype, mutations, score, goal, domain):
        """构建分析提示词"""
        goal_cn = {"stability": "稳定性", "activity": "活性", "solubility": "可溶性"}.get(goal, goal)
        
        return f"""作为蛋白质工程专家，请分析以下突变对{goal_cn}的影响：

突变：{', '.join(mutations)}

蛋白质信息：
- 类型：{domain}
- 优化目标：提高{goal_cn}
- 野生型长度：{len(wildtype)} AA
- 预测得分：{score:.3f}/1.0

请提供：
1. 每个突变的生物物理效应分析
2. 对结构和功能的潜在影响
3. 实验验证建议
4. 风险评估

要求：专业、实用、可操作，用中文回答。"""