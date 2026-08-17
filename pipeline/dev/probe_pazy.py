#!/usr/bin/env python3
"""探测 PAZy API 文档与端点"""
import requests
import re

r = requests.get("https://www.pazy.eu/api-docs", timeout=20)
print("status:", r.status_code, "len:", len(r.text))
# 常见 swagger 位置
for pat in [r'"(/api[^"?]*)"', r"'(/api[^']*)'", r'url:\s*"([^"]+)"', r'"paths"\s*:\s*\{([^}]{0,500})']:
    m = re.findall(pat, r.text)
    if m:
        print("PATTERN", pat, "->", sorted(set(m))[:20])
# 尝试 swagger json
for u in ["https://www.pazy.eu/swagger.json", "https://www.pazy.eu/openapi.json", "https://www.pazy.eu/api/openapi.json", "https://www.pazy.eu/api/swagger.json"]:
    try:
        rr = requests.get(u, timeout=15)
        print(u, "->", rr.status_code, len(rr.text))
        if rr.status_code == 200 and rr.text.strip().startswith("{"):
            try:
                j = rr.json()
                print("  paths:", list(j.get("paths", {}).keys())[:25])
            except Exception as e:
                print("  json err", e)
    except Exception as e:
        print(u, "ERR", e)
