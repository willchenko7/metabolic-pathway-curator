

#pathways = findAllPathsBetweenTwoMetabolites('pyruvate','ethanol',metabolic_graph)

def findAllPathsBetweenTwoMetabolites(precursor,target,metabolic_graph,reaction_dict):
    paths = find_all_paths(metabolic_graph,precursor,target,path=[])
    pathways = []
    for path in paths:
        path_pathways = findAllPathwaysFromPath(path,reaction_dict)
        pathways.append(path_pathways)
    return pathways

def find_all_paths(graph, start, end, path=[]):
        path = path + [start]
        if start == end:
            return [path]
        if start not in graph:
            return []
        paths = []
        for node in graph[start]:
            if node not in path:
                newpaths = find_all_paths(graph, node, end, path)
                for newpath in newpaths:
                    paths.append(newpath)
        return paths

def findReactionsFromNodes(met1,met2,reaction_dict):
    rxns = []
    rxn_indices = []
    [rxn_list, rxn_index_list] = createRxnIndex(reaction_dict)
    for rxn in reaction_dict['reactions']:
        if met1 in reaction_dict[rxn]['nodes'] and met2 in reaction_dict[rxn]['nodes']:
            rxns.append(rxn)
            rxn_indices.append(rxn_list.index(rxn))
    return rxns, rxn_indices

def findAllPathwaysFromPath(path,reaction_dict):
    pathways = [[]]
    for i in range(1,len(path)):
        print(i)
        met1 = path[i]
        met2 = path[i-1]
        [rxns, rxn_indices] = findReactionsFromNodes(met1,met2,reaction_dict)
        if len(rxns) == 1:
            for j in range(0,len(pathways)):
                pathways[j].extend(rxns)
        elif len(rxns) > 1:
            olenp = len(pathways)
            pathways = duplicateLoLNumber(pathways,len(rxns))
            start = 0
            stop = olenp
            for rxn in rxns:
                print(rxn)
                for k in range(start,stop):
                    print(k)
                    pathways[k] = pathways[k] + [rxn]
                start = start + olenp
                stop = stop + olenp
    return pathways

def copyLoL(lol):
    nlol = []
    for i in range(0,len(lol)):
        nlol.append(lol[i])
    return nlol

def duplicateLoLNumber(lol,number):
    nlol = copyLoL(lol)
    for j in range(0,number-1):
        nlol.extend(lol)
    return nlol
    
def createRxnIndex(reaction_dict):
    rxn_list = reaction_dict['reactions']
    rxn_index_list = list(range(0,len(rxn_list)))
    return rxn_list, rxn_index_list
