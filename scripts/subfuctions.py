import sys 
import json
import numpy as np
import re

# load FEL results
def load_json_file(file_handle):
    raw = re.sub(r'\binf\b', 'Infinity', file_handle.read())
    try:
        json_dict = json.loads(raw)
        return json_dict
    except Exception as e:
        print(e)
        sys.exit(1)
    
    
# create list of str when nargs='+' is used for arguments
def convert_nargs_to_list(input):
    input_list = [x.replace("[","").replace("]","").replace(",","") for x in input]
    return input_list