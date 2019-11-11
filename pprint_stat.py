import json
import sys
with open(sys.argv[1], 'r') as handle:
    parsed = json.load(handle)
print(json.dumps(parsed, indent=4, sort_keys=True))
