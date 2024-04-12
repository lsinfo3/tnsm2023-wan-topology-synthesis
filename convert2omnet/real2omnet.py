# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 18:28:03 2021

@author: katha
"""


import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import pickle

def graph2omnet(G, name, plot = False):
    
    if plot:
        plt.figure(figsize=(10,5))
        ax = plt.gca()
        ax.set_title(name)
        pos=nx.spring_layout(G)
        weights = list(nx.get_edge_attributes(G,'weight').values())
        weights = [np.log10(x) + 1 for x in weights]
        nx.draw(G,pos,width=weights,with_labels=True)
        plt.show()
        
    n = len(G.nodes)
    
    mst = nx.minimum_spanning_edges(G, algorithm='kruskal', weight='weight', keys=True, data=True, ignore_nan=False)
    
    with open('skeletons/openflow_ned.txt') as f:
        lines_of = f.readlines()
    lines_of.append("\n")   
    

    
    for i in range (n):
        lines_of.append("ofs_"+str(i)+": Open_Flow_Domain;\n")
 
    lines_of.append("\n")    
    lines_of.append("connections allowunconnected:\n")
    
    
    for (u,v,w) in G.edges(data=True):
        lines_of.append("ofs_"+str(u)+".gateDPlane++ <--> DistanceChannel {  distance = "+str(w.get('weight'))+"km; } <--> ofs_"+str(v)+".gateDPlane++;\n")

    lines_of.append("\n")

    for (u,v,w) in list(mst):
        lines_of.append("ofs_"+str(u)+".gateCPlane++ <--> DistanceChannel {  distance = "+str(w.get('weight'))+"km; } <--> ofs_"+str(v)+".gateCPlane++;\n")
       
    lines_of.append("\n")   
    
    numControllers = [1]
      
    numPlacements = 1
    placement_id = 1
    
    with open('skeletons/openflow_ini.txt') as skele_ini:
        lines_ini_of = skele_ini.readlines()
    lines_ini_of.append("\n")
    
    
    for numC in numControllers:  # these two loops dont really do anything, since both only loop one time and "controller1" is hardcoded below (artifact from trying out distributed architectures)
        for j in range(numPlacements):
            lines_of.append("\n")
            lines_ini_of.append("\n")

            
            lines_ini_of.append("[Config cfg"+str(placement_id)+"]\n")
            lines_ini_of.append("**.placementID = "+str(placement_id)+"\n")
            lines_ini_of.append("**.numControllers = "+str(numC)+"\n")
            
            closeness = nx.closeness_centrality(G, distance="weight")
            closest_node = max(closeness, key=closeness.get)
 

            lines_of.append("ofs_"+str(closest_node)+".gateCPlane++ <--> backboneline <--> open_flow_controller1"+".ethg++ if (numControllers == "+str(numC)+" && placementID == "+str(placement_id)+");\n")

            for node in G.nodes:   
                    lines_ini_of.append("**.ofs_"+str(node)+".open_flow_switch*.OF_Switch.connectAddress = \"open_flow_controller1"+"\"\n")

            
    lines_of.append("}")

    
    with open("random_networks/"+str(name)+'.ned', 'w') as ned:
        for line in lines_of:
            line = line.replace("PLACEHOLDER_NETWORK", name)
            ned.write("%s" % line)
            
    with open("random_ini/"+str(name)+'.ini', 'w') as ini:
        for line_ini in lines_ini_of:
            line_ini = line_ini.replace("PLACEHOLDER_NETWORK", "openflow.networks.synthetic."+name)
            ini.write("%s" % line_ini)
            
    

names = ["GTSSLOVAKIA", "BREN", "GEANT2001", "BTNORTHAMERICA"] # one of GEANT2001, BREN, BTNORTHAMERICA, GTSSLOVAKIA
names = ["GEANT2001"]

for name in names:
    path = "C:\\cnsm2022-wan-topology-synthesis\\adj_matrices\\" + name + "_weighted.txt"
    real = np.loadtxt(path)
    real_graph = nx.from_numpy_array(real)
    graph2omnet(real_graph, "REAL_"+name, plot = False)