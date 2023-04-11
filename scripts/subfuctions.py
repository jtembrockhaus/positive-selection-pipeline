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
        
# create list of str when nargs='+' is used for arguments
def convert_nargs_to_list(input):
    input_list = [x.replace("[","").replace("]","").replace(",","") for x in input]
    return input_list