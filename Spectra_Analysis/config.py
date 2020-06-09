import json

CONFIG = {}

with open('config.json') as f:
    CONFIG = json.load(f)