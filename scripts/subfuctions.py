import sys 
import json

# load FEL results
def load_json_file(filename):
    try:
        json_dict = json.load(filename)
        return json_dict
    except Exception as e:
        print(e)
        sys.exit(1)