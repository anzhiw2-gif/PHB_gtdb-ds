#!/usr/bin/env python3
"""分析 PAZy api-docs 页面内容，找 API 调用示例"""
import requests
import re

r = requests.get("https://www.pazy.eu/api-docs", timeout=20)
txt = r.text
# 找 JS 脚本 / 接口地址 / fetch 示例
scripts = re.findall(r'src="([^"]+)"', txt)
print("scripts:", scripts[:10])
for pat in [r'https?://[^"\'<> ]+', r'/api[a-zA-Z0-9/_-]*', r'fetch\([^)]{0,80}']:
    found = sorted(set(re.findall(pat, txt)))
    print(f"\nPATTERN {pat}:")
    for f in found[:25]:
        print("  ", f[:150])
