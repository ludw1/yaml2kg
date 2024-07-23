This script converts yaml files to a knowledge graph html. It utilizes the storage of the documentation in a context-aware data format seen in [this json file](json/rest_lessvariables.json).
## Usage
1. Install required python packages using pip install -r requirements.txt
2. Call the script with the yaml file as the positional argument and optionally the -o flag to specify the output html name: ```python yaml2kg.py myyaml.yaml -o kg.html```

## Notes
If there are TupleTools in the yaml that are not in the updated dictionary, only a node with the TupleTool name and applied options will be added.
