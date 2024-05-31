This script converts yaml files to a knowledge graph html. It utilizes the storage of the documentation in a context-aware data format seen in [this json file](json/rest_lessvariables.json).
## Usage
After installing the required python packages using ```pip install -r requirements.txt``` calling the script with the specified yaml file will produce a html file in the output directory which can also be specified using the -o flag.

## To-Dos
- Add support for Loki Hybrid tupletool.
- Throw exception when there a tool thats not in the documentation is present in the yaml
- Decide on visual form of the knowledge graph


## Acknowledgements
Special thanks to Mr. Yo Mizutani for the visual.py, enabling visualization through pyvis: https://users.cs.utah.edu/~yos/2021/02/02/plotly-python.html.