# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 16:36:17 2022

@author: katha
"""

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pickle
import scipy 
import scipy.stats
import pandas as pd
import copy

def propagate_data(G, L):
    for (a,b,c) in G.edges(data =True):
        weight = c['weight']
        nx.set_node_attributes(L, values = {(a,b):weight}, name='weight')

def get_edge_attributes(G, name): # https://stackoverflow.com/questions/62564983/how-to-efficiently-get-edge-weights-from-a-networkx-graph-nx-graph
    edges = G.edges(data=True)
    return dict( (x[:-1], x[-1][name]) for x in edges if name in x[-1] )
    
def get_node_attributes(G, name):
    nodes = G.nodes(data=True)
    return dict( (x[:-1], x[-1][name]) for x in nodes if name in x[-1] )

def get_bin(value, edges):
    if value <= edges[0]:
        return 0
    
    for i in range(0, len(edges)-1):
        if (edges[i] <= value)and (value <= edges[i+1]):
            bin_idx = np.where(edges == edges[i])
            return bin_idx[0][0]
        
    if (value >= edges[len(edges)-1]):
            return len(edges)-2
        
def get_clostest_bin(value, theoretical_edges, real_edges):
    diffs = np.abs(real_edges - value)
    mindiff_idx = np.argmin(diffs)
    bin_val = real_edges[mindiff_idx]
    
    # if the closest bin edge is smaller than the value to bin, then the bin edge is the upper limit -> lower index by one to assign the correct bin
    if bin_val < value:
        bin_idx = np.where(theoretical_edges == bin_val)[0][0]-1
    # else the closest bin edge is bigger, than it's the lower bin edge -> index already correct
    else:
        bin_idx = np.where(theoretical_edges == bin_val)[0][0]
    return bin_idx


name = "BREN"
metrics = ["EBC","BC","CC","DC"]
adj_matrix_real = np.loadtxt("C:\\cnsm2022-wan-topology-synthesis\\adj_matrices\\"+name+"_weighted.txt")
path = "C:\\cnsm2022-wan-topology-synthesis\\"+name+"_naive\\"

for metric in metrics:
    np.random.seed(0)

    G = nx.Graph(np.squeeze(adj_matrix_real),data = True)
    if metric == "BC" or metric == "CC" or metric == "DC":
        L=nx.line_graph(G)
        propagate_data(G,L)
        G=L
    
    weights_all = adj_matrix_real[adj_matrix_real != 0] 
    bin_edges_weights = np.sort(np.unique(weights_all)) # we do not bin multiple edge weights together, as they are too different, we utilize all unique values as bin edges
    
    # raw values of edge weights
    if metric == "BC" or metric == "CC" or metric == "DC":
        weights_edges = list(get_node_attributes(G, 'weight').values()) # from node attr for line graph approach
    if metric == "EBC":
        weights_edges = list(get_edge_attributes(G, 'weight').values()) # normal if edge betweenness
        
    # get the raw values for the respective centrality metric of the edges (or the nodes for the line graphs)
    if metric == "BC":
        centrality_edges = np.array(list(nx.betweenness_centrality(G).values()))
    if metric == "CC":
        centrality_edges = np.array(list(nx.closeness_centrality(G).values()))
    if metric == "DC":
        centrality_edges = np.array(list(nx.degree_centrality(G).values()))
    if metric == "EBC":
        centrality_edges = np.array(list(nx.edge_betweenness_centrality(G).values()))    
        
    for numbins in np.flip([1,2,4, 8, 16, 32, 64, 128]):
        #bin_egdes_bcw = np.linspace( np.min(bcsw_edges), np.max(bcsw_edges), numbins + 1)
        
        # bin the whole possible value range (0 to 1 for normalized centralities)
        bin_egdes_centrality = np.linspace(0, 1, numbins + 1)
        
        # contains bin edges etc., we utilize "statistic='mean'" to get the "mean" value of a joint bin for the edge weights, which is essentially one of the real weights, so we can assign it later
        stats = scipy.stats.binned_statistic_2d(weights_edges, centrality_edges, weights_edges,statistic='mean', bins=[bin_edges_weights, bin_egdes_centrality],expand_binnumbers=True)
        
        # we utilize "statistic='count'" to calculate the joint probability to sample from later
        jointProbs = scipy.stats.binned_statistic_2d(weights_edges, centrality_edges, None,statistic='count', bins=[bin_edges_weights, bin_egdes_centrality],expand_binnumbers=True)[0]
        
        # as some bins may be empty, we figure out the "real" bins that contain values
        c = np.where(jointProbs.any(axis=0))[0]
        a_upper = bin_egdes_centrality[c]
        a_lower =  bin_egdes_centrality[c+1]
        a = np.unique(np.concatenate((a_lower,a_upper)))
    
        means_cc = []
        means_bc = []
        means_weights = []
    
        graph_metrics=[]
    
        # for k in range(100):
        #     for s in range(10):
        for k in range(0,1000):
                with open(path +"synth_sample_2K_"+str(k)+"_.pkl","rb") as f:
        #       with open(path +"synth_sample_"+str(k)+"_"+str(s)+"_"+"BW"+"_"+str(int(True))+"_.pkl","rb") as f:
                    L_synth = pickle.load(f)
                    
                if metric == "BC" or metric == "CC" or metric == "DC":
                    G_synth = nx.line_graph(L_synth)
                    propagate_data(L_synth,G_synth)
                    
                if metric == "EBC":
                    G_synth = copy.deepcopy(L_synth)
                    
                # calculate centralities for the synthetic network
                if metric == "BC":
                    centralities = nx.betweenness_centrality(G_synth)
                if metric == "CC":
                    centralities = nx.closeness_centrality(G_synth)
                if metric == "DC":
                    centralities = nx.degree_centrality(G_synth)
                if metric == "EBC":
                    centralities = nx.edge_betweenness_centrality(G_synth)
                 
                if metric == "BC" or metric == "CC" or metric == "DC":
                        
                    for u,w in G_synth.nodes(data = True):
                        centrality = centralities[(u)] # specific centrality for a specific node
                        bin_id = get_bin(centrality, stats[2]) # get the theoretical bin it belongs to
                        rel_freq = jointProbs[:,bin_id] # get the probability corresponding to the bin
                        probs = rel_freq/rel_freq.sum() # convert count to actual probability
                    
                        if (np.isnan(probs.sum())): # if theoretical bin is empty            
                            bin_id = get_clostest_bin(centrality, bin_egdes_centrality, a) # get the closet real bin
                            rel_freq = jointProbs[:,bin_id] # repeat above
                            probs = rel_freq/rel_freq.sum() # repeat above
        
                        means = stats[0][:,bin_id] # get the edge weight corresponding to the mean bin value of the weights
                        new_weight = np.random.choice(means, p=probs) # sample according to given joint probability
                        w['weight'] = new_weight # assign weight to line graph
                        
                        L_synth.edges[(u)]["weight"] = new_weight # backpropagate data to the real graph (not the line graph)
        
                if metric == "EBC":
                    for u,v,w in G_synth.edges(data = True):
                        centrality = centralities[(u,v)]
                        bin_id = get_bin(centrality, stats[2])
                        rel_freq = jointProbs[:,bin_id]
                        probs = rel_freq/rel_freq.sum()
                        
                        if (np.isnan(probs.sum())):                 
                            bin_id = get_clostest_bin(centrality, bin_egdes_centrality, a)
                            rel_freq = jointProbs[:,bin_id]
                            probs = rel_freq/rel_freq.sum()

                        
                        means = stats[0][:,bin_id]
                        new_weight = np.random.choice(means, p=probs)
                        w['weight'] = new_weight
                        L_synth = G_synth
                        
                with open("sampling_2k/synth_sample_"+str(k)+"_2K_"+metric+"_b"+str(numbins)+".pkl","wb") as f:
                            pickle.dump(L_synth,f)     
                # aaaaaaaand calculate the (weighted) centrality metrics for comparison later
                ccw = np.mean(list(nx.closeness_centrality(L_synth, distance='weight').values()))
                bcw = np.mean(list(nx.betweenness_centrality(L_synth, weight='weight').values()))
                
                means_weights = means_weights + [sum(list(get_edge_attributes(L_synth, 'weight').values()))]
                means_bc= means_bc+[bcw] 
                means_cc= means_cc+[ccw]
        
                graph_metrics =  graph_metrics + [[ k, sum(list(get_edge_attributes(L_synth, 'weight').values())),bcw,ccw]]
        
        graph_metrics_df = pd.DataFrame(graph_metrics, columns=["i", "Weights", "BCW", "CCW"])    
        graph_metrics_df.to_csv("graph_metrics_node_"+metric+"_sampling_"+name+"_2K_"+str(numbins)+".csv",index=False,sep=";")
