This script converts yaml files to a knowledge graph html. It utilizes the storage of the documentation in a context-aware data format seen in [this json file](https://github.com/ludw1/yaml2kg/blob/plotly_script/json/rest_lessvariables.json).
## Usage
1. Install using pip via ```pip install git+https://github.com/ludw1/yaml2kg```
2. Run the command ```yaml2kg myyaml.yaml -o output.html``` replacing myyaml and output with the desired input yaml file and desired output html file path

## Developing
To develop yaml2kg, simply clone the repository and install the required packages via the requirements.txt.

## Notes
If there are TupleTools in the yaml that are not in the updated dictionary, only a node with the TupleTool name and applied options will be added.
