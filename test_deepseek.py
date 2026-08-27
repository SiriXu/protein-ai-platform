# test_deepseek.py
import sys
import os

# 你的API密钥（临时测试用）
API_KEY = "你的API密钥"

def test_direct():
    """直接测试DeepSeek API"""
    try:
        from openai import OpenAI
        
        print("测试1: 不带/v1后缀")
        try:
            client = OpenAI(
                api_key=API_KEY,
                base_url="https://api.deepseek.com",
                timeout=30
            )
            
            response = client.chat.completions.create(
                model="deepseek-chat",
                messages=[{"role": "user", "content": "Hello"}],
                max_tokens=10
            )
            print(f"✅ 测试1成功: {response.choices[0].message.content}")
        except Exception as e:
            print(f"❌ 测试1失败: {e}")
        
        print("\n测试2: 带/v1后缀")
        try:
            client2 = OpenAI(
                api_key=API_KEY,
                base_url="https://api.deepseek.com/v1",
                timeout=30
            )
            
            response2 = client2.chat.completions.create(
                model="deepseek-chat",
                messages=[{"role": "user", "content": "Hello"}],
                max_tokens=10
            )
            print(f"✅ 测试2成功: {response2.choices[0].message.content}")
        except Exception as e:
            print(f"❌ 测试2失败: {e}")
            
    except ImportError:
        print("❌ 请安装openai: pip install openai")
    except Exception as e:
        print(f"❌ 测试异常: {e}")

def test_with_requests():
    """使用requests直接测试"""
    import requests
    import json
    
    print("\n使用requests测试API端点...")
    
    url = "https://api.deepseek.com/chat/completions"
    headers = {
        "Authorization": f"Bearer {API_KEY}",
        "Content-Type": "application/json"
    }
    
    data = {
        "model": "deepseek-chat",
        "messages": [{"role": "user", "content": "Hello"}],
        "max_tokens": 10
    }
    
    try:
        response = requests.post(url, headers=headers, json=data, timeout=30)
        print(f"状态码: {response.status_code}")
        if response.status_code == 200:
            result = response.json()
            print(f"✅ 请求成功: {result['choices'][0]['message']['content']}")
        else:
            print(f"❌ 请求失败: {response.text}")
    except Exception as e:
        print(f"❌ requests测试失败: {e}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        API_KEY = sys.argv[1]
    
    if not API_KEY or API_KEY == "你的API密钥":
        print("请提供API密钥: python test_deepseek.py your_api_key_here")
        sys.exit(1)
    
    print(f"测试API密钥: {API_KEY[:10]}...")
    test_direct()
    test_with_requests()