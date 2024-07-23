#!/usr/bin/env python
import yaml
import json
import argparse
import networkx as nx
import uuid
from pyvis.network import Network
import gzip
import requests

restr_label = ["basic","head","charged"] # labels that indicate variables which are restricted to certain particles

# node colors
exp_color = "purple"
var_color = "green"
option_color = "orange"
option_value_color = "skyblue"
particle_color = "red"
ttcolor = "blue"
functor_color = "mint"

# file paths
var_file = r".\json\rest_lessvariables.json"
config_file = "config.yaml"
files = {}

# options
link_hover  = False # create nodes from hover data

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(prog = "yaml2kg", description="Create a graph from a yaml file.")
    parser.add_argument("file", type=str, help="The path to the yaml file.")
    parser.add_argument("-o", "--output", type=str, help="The path to the output file. Defaults to graph.html.")
    return parser.parse_args()

def get_metadata() -> bool:
    """
    Function that fetches the metadata from the backend.
    Uses the config file to get the filenames and the backend url.
    Data stored in global file object.
    return: True if successful, False otherwise."""
    with open(config_file,"r") as f:
        config = yaml.safe_load(f)
    for file in config['filenames']:
        response = requests.get(config['backend_url'] + list(file.values())[0] + ".json.gz", headers = {"Accept-Encoding": "gzip"})
        if response.status_code == 200:
            data = str(gzip.decompress(response.content), "utf-8") # decode gzip file
            files[list(file.keys())[0]] = json.loads(data) # turn string back into dict
        else:
            print(f"Error: {response.status_code}")
            return False
    return True
def get_mapped_list(decay_temp: str) -> list:
    """Get the mapped list from the decay file. This takes the str from the yaml.
    decay_temp: The decay template from the yaml file.
    return: The mapped list from the decay file."""
    data = files['decay_file']
    for i in list(data.values()):
        desc = i["descriptors"]
        if decay_temp == desc["template"]:
            return desc["mapped_list"]
                
def build_decay_graph(decay: list) -> nx.DiGraph:
    """
    This constructs a nx.DiGraph from a nested decay list. 
    param decay: A nested list of particles, where a nested list indicates a subdecay. Each particle has a branch and a particle name.
    return: A nx.DiGraph object representing the decay tree.""" 

    def get_graph_params(term: list | dict, graph: nx.DiGraph, parent_id: str) -> nx.DiGraph:
        """
        Recursive function that constructs the decay tree.
        term: A list of particles, where nesting indicates a subdecay.
        graph: The current nx.DiGraph object.
        parent_id: The name of the current decay head.
        return: The updated nx.DiGraph object."""
        if not isinstance(term, list): # if term / the particle is not a list then turn into a node and connect to the decay head
            branch,particle = term.values() # branch is unique, so use it for nodes
            new_node = f"{branch}"
            label = f"{particle}" # use particle name for label
            new_edge = "" if parent_id == 0 else (f"{parent_id}",f"{branch}") # dont connect the decay head to anything
            return {'new_node': new_node, 'new_edge': new_edge, 'label':label}

        for element in term:
            if isinstance(element, list): # if there are subdecays
                sub_parent_id = list(graph.nodes())[-1] # get the current decay head
                graph = get_graph_params(element, graph, sub_parent_id) # call decay tree function with the subarray
            else:
                params = get_graph_params(element, graph, parent_id)
                graph.add_node(params['new_node'],label = params['label'], color = particle_color, type = "particle")
                if params["new_edge"] != "":
                    graph.add_edge(*params['new_edge'], label = "decays")
                    
        return graph
    G = nx.DiGraph()
    graph = get_graph_params(decay, G, 0)
    return graph


def loop_varjson(
        G: nx.DiGraph, 
        tool: dict, 
        particle: str, 
        variables: dict, 
        label_dict: dict,
        tupletoolname: uuid.UUID
    ) -> dict:
    """
    Recursive function that loops through the documentation file and links the variables to the tuple tool.
    G: The nx.DiGraph object representing the decay tree. Will be modified in place so no need to return it.
    tool: The current tuple tool. Dict with the name of the tool as key and option:value pairs as values.
    particle: The current particle name.
    variables: The variable dict of the current tuple tool.
    label_dict: A dictionary that maps uuids to labels.
    tupletoolname: The UUID of the current tuple tool.
    return: The updated label_dict dictionary."""
    
    options = (list(tool.values())[0]) # options of the tupletool
    for var in variables.items(): # check all variables
        if isinstance(var[1],dict): # if the variable is a dict, it is either locked behind an option or a restriction
            if (var[0].split(",")[0].replace("!","") not in restr_label): # if it is not a restriction check if the option is true
                for option in options:
                    if option in list(variables.keys()):
                        if options[option]:
                            loop_varjson(G,tool,particle,variables[option],label_dict,tupletoolname) # recursively go through 
            else:
                restriction = var[0] # the variable name
                pconfrom = 1
                for res in restriction.split(","): #check if each restriction is fulfilled
                    pprop = G._node[particle][res] if "!" not in res else not G._node[particle][res.replace("!","")]                        
                    pconfrom = pconfrom and pprop
                if pconfrom:
                    loop_varjson(G,tool,particle,variables[var[0]],label_dict,tupletoolname) # if yes, go through the variables
        else: # if the variable is a string
            varname = var[0].replace("head",particle+"_"+options["ExtraName"]) # change the variable name to match output
            for option in list(options.keys()): # check if there are variable names which depend on options
                if option in var[0]:
                    varname = varname.replace("head",particle+"_"+options["ExtraName"]).replace(option,str(options[option]))
            varexname = variables[var[0]] # get the explanation of the variable
            G.add_node(varname, color = var_color, label = varname, expl = varexname, type = "variable")
            G.add_edge(tupletoolname,varname)
            if (link_hover): # create a node for the hover data
                varexname = uuid.uuid4()
                label_dict[varexname] = "{}".format(variables[var[0]])
                G.add_node(varexname, color = exp_color, label = varexname, type = "explanation")
                G.add_edge(varname,varexname)
    return label_dict
def link_loki(G: nx.DiGraph, tool: dict, variables: dict, label_dict: dict) -> dict:
    """
    This function is used to link the LoKi functors to the graph.
    G: The usual nx.DiGraph that contains everything.
    tool: The current tuple tool dict.
    variables: The variable dict of the current tuple tool.
    label_dict: A dictionary that maps uuids to labels.
    return: The updated label_dict dictionary.
    """
    options = (list(tool.values())[0])
    for option in options:
        if "Variables" in option:
            optionname = list(label_dict.keys())[list(label_dict.values()).index(option)]
            func_dict = options[option]
            for func in func_dict: # every user defined variable that will be written out
                funcname = uuid.uuid4()
                label_dict[funcname] = "{}".format(func)
                # this creates a node for the user given name with hover data for the loki functor
                G.add_node(funcname, color = var_color, label = func, expl = func_dict[func], type = "loki_variable")
                G.add_edge(optionname,funcname)
                # now match loki functors to documentation
                split = "".join([x if x.isalpha() else " " for x in func_dict[func]]).split(" ")
                for lokifunc in split:
                    if lokifunc:
                        for var in variables.items():
                            if lokifunc == var[0].split(",")[0]:
                                varexplname = uuid.uuid4()
                                label_dict[varexplname] = "{}".format(lokifunc)
                                G.add_node(varexplname, color = functor_color, label = lokifunc, expl = var[1], type = "loki_expl")
                                G.add_edge(funcname,varexplname)
                                

    return label_dict
def link_var(G: nx.DiGraph, tool: dict, particle: str) -> dict: 
    """
    Check if the tupletool should be linked to the particle and calls the function that loops through var.json.
    G: The usual nx.DiGraph that contains everything.
    tool: The current tuple tool dict.
    particle: The current particle name.
    return: The updated label_dict dictionary.
    """
    if "LoKi" in list(tool.keys())[0]: # if the tool is a LoKi tool, handle it differently
        label_dict = create_tuple_tool(G, tool, particle)
        tupletoolname = list(label_dict.keys())[list(label_dict.values()).index(list(tool.keys())[0])]
        with open(var_file) as f:
            data = json.load(f)
            variables = data["LoKi__Hybrid__TupleTool"]["variables"]
            label_dict = label_dict | link_loki(G, tool, variables, label_dict)
    else:
        with open(var_file) as f:
            data = json.load(f)
            if list(tool.keys())[0].split("/")[0] not in data.keys(): # if the tool is not in the var.json, only add the tuple tool node and the options
                label_dict = create_tuple_tool(G, tool, particle) # create the tuple tool node
                return label_dict
            variables = data[list(tool.keys())[0].split("/")[0]]["variables"]
            if isinstance(list(variables.values())[0],dict): # this means that there are no variables which every particle can have
                        restriction = list(variables.keys())[0] # check if the current particle fulfills restrictions
                        pconfrom = 1
                        for res in restriction.split(","): # for multiple restrictions seperated by ","
                            pprop = G._node[particle][res] if "!" not in res else not G._node[particle][res.replace("!","")] # ! = not (basic e.g.)
                            pconfrom = pconfrom and pprop
                        if pconfrom: # if yes, move on as usual
                            label_dict = create_tuple_tool(G,tool,particle)
                            tupletoolname = list(label_dict.keys())[list(label_dict.values()).index(list(tool.keys())[0])]
                        else: # else dont link the tupletool and pass back
                            return {}
            else:
                label_dict = create_tuple_tool(G, tool, particle) # create the tuple tool node
                tupletoolname = list(label_dict.keys())[list(label_dict.values()).index(list(tool.keys())[0])] # get the tupletool name written in the var.json
            label_dict = label_dict | loop_varjson(G, tool, particle, variables, label_dict, tupletoolname) # link variables to tuple tool
    return label_dict
    


def create_tuple_tool(G: nx.DiGraph, tool: dict, particle: str) -> dict:
    """Creates tuple tool node and connects it with the options.
    return: A dictionary that maps uuids to labels."""
    tupletoolname = list(tool.keys())[0]
    options = (list(tool.values())[0])

    uuid_mapping = {}
    # create uuid for every tupletool so they appear distinct in the graph
    ttname = uuid.uuid4()
    G.add_node(ttname, color = ttcolor, type = "tupletool", lokitool = "LoKi" in tupletoolname)
    G.add_edge(particle, ttname)
    uuid_mapping[ttname] = "{}".format(tupletoolname)

    for key in options.keys():
        keyname = uuid.uuid4()
        uuid_mapping[keyname] = "{}".format(key)
        G.add_node(keyname, color = option_color, optval = str(options[key]), type = "option", loki = "LoKi" in tupletoolname)
        G.add_edge(ttname,keyname)
        if(link_hover):
            optionname = uuid.uuid4()
            uuid_mapping[optionname] = "{}".format(str(options[key]))
            G.add_node(optionname, color = option_value_color, type = "option_value") 
            G.add_edge(keyname, optionname)
    return uuid_mapping

def link_all(graph: nx.DiGraph, config: dict, apl_tools: list) -> dict:
    """
    The main function that links all tupletools, variables and explanations to the decay tree.
    It goes backwards through the yaml and links tupletools that have unique tuplenames and extranames.
    graph: The usual nx.DiGraph that contains everything.
    config: The yaml file as a dict.
    apl_tools: A list of all applied tools.
    return: The updated label_dict dictionary.
    """
    label_dict = {} # keeps track of uuids 
    for group in list(config["groups"].keys())[::-1]: # go through groups in reverse order
        for tool in config["groups"][group]["tools"][::-1]:
            toolname = list(tool.keys())[0] + "{}".format(tool[list(tool.keys())[0]].get('ExtraName',"")) # this ensures that the tool has a different name and extra name 
            for particle in group.split(","):
                if toolname in apl_tools[particle]: # if the tool was already applied, skip it
                    pass
                else: # otherwise link it and mark that is has been applied
                    label_dict = label_dict | link_var(graph,tool,particle)
                    apl_tools[particle].append(toolname)
    for particles in list(config["branches"].items()):
        particle_name = particles[0]
        for tool in particles[1]["tools"][::-1]:
            toolname = toolname = list(tool.keys())[0] + "{}".format(tool[list(tool.keys())[0]].get('ExtraName',""))
            if toolname in apl_tools[particle_name]:
                pass
            else:
                label_dict = label_dict | link_var(graph,tool,particle_name)
                apl_tools[particle_name].append(toolname)
    for tool in config["tools"]:
        toolname = list(tool.keys())[0] + "{}".format(tool[list(tool.keys())[0]].get('ExtraName',""))
        if "Event" in toolname and toolname not in apl_tools[particle]:
            graph.add_node("Event", color = "red", label = "Event", type = "particle")
            graph.add_edge("Event",list(graph.nodes)[0])
            label_dict = label_dict | link_var(graph,tool,"Event")
        else:
            for particle in list(apl_tools.keys()):
                if toolname in apl_tools[particle]:
                    pass
                else:
                    label_dict = label_dict | link_var(graph,tool,particle)
                    apl_tools[particle].append(toolname)
    return label_dict

def style_graph(G: nx.DiGraph, label_dict: dict):
    """
    Makes the graph pretty.
    G: The nx.DiGraph object representing the decay tree.
    label_dict: A dictionary that maps uuids to labels.
    """
    maxi = 0 # maximum depth of the graph
    for i, layer in enumerate(nx.topological_generations(G)):
        for n in layer:
            G.nodes[n]["level"] = i
            if i > maxi:
                maxi = i
    for node, label in label_dict.items(): # needed for pyvis as it cant use the label dict
        G.nodes[node]["label"] = str(label)

    if not link_hover: # set the data on hover
        for node in G.nodes():
            if G.nodes[node]["type"] in ["variable", "loki_variable", "loki_expl"]:
                G.nodes[node]["title"] = G.nodes[node]["expl"]
            elif G.nodes[node]["type"] == "option":
                G.nodes[node]["title"] = G.nodes[node]["optval"]
            elif G.nodes[node]["type"] == "particle" and G.nodes[node]["label"] != "Event":
                G.nodes[node]["title"] = f"{G.nodes[node]['label']} \n Basic: {G.nodes[node]['basic']} \n Head: {G.nodes[node]['head']} \n Charge: {G.nodes[node]['charge']} \n Mass: {round(float(G.nodes[node]['mass']), 2)} MeV"
            else:
                G.nodes[node]["title"] = G.nodes[node]["label"]
    
def main(): 
    """Opens NTuple Wizrad yaml file, calls every other function.
    """
    args = parse_args()
    output = args.output if args.output else r"graph.html"
    yaml_file = args.file
    if not get_metadata():
        print("Error: Could not fetch metadata.")
        return
    with open(yaml_file,"r") as f:
        config = yaml.safe_load(f)
    decayl = get_mapped_list(config["descriptorTemplate"])

    applied_tupletools = {}
    graph = build_decay_graph(decayl) # construct decay graph
    data = files['prop_file'] # edit particle properties
    for node in graph.nodes():
        if node == list(graph.nodes)[0]:
            graph._node[node]["head"] = 1
        else:
            graph._node[node]["head"] = 0
        graph._node[node]["charged"] = 1 if data[graph._node[node]["label"]]["charge"] != 0 else 0
        graph._node[node]["charge"] = data[graph._node[node]["label"]]["charge"]
        graph._node[node]["mass"] = data[graph._node[node]["label"]]["mass"]
        graph._node[node]["basic"] = 0 if graph.out_degree(node) != 0 else 1
        applied_tupletools[node] = []
    files.clear() # clear the files to save memory, as decay file is ~50 MB

    label_dict = link_all(graph,config,applied_tupletools) # link all tupletools, variables etc.
    style_graph(graph,label_dict)

    # draw graph
    G = nx.convert_node_labels_to_integers(graph)
    

    net = Network(
        notebook=False,
        # bgcolor="#1a1a1a",
        cdn_resources="remote",
        height="900px",
        width="100%",
        select_menu=True,
        # font_color="#cccccc",
        filter_menu=True
    )

    net.from_nx(G)
    node_l = [node for node in list(net.nodes) if node["type"] == "particle"]
    counts = {}
    
    for node in node_l:
        counts[node["level"]] = counts.get(node["level"], 0) + 1
        node["physics"] = False # disable physics for particles
        node["size"] = 15
    net.force_atlas_2based(central_gravity=0.01, gravity=-31,overlap = 1)
    for level in sorted(counts): # this is so that the particles are aligned in a decay tree
            i = 0
            for node in net.nodes:
                    if node["color"] == particle_color:
                        if node["level"] == level:
                                node["y"] = node["level"]*200-1000
                                node["x"] = counts[node["level"]]*20*i-counts[node["level"]]*20
                                i += 1
    net.save_graph(output)

if __name__ == "__main__":
    main()