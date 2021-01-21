# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 20:05:29 2020

@author: gosl241
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 20:05:29 2020
@author: gosl241
"""

import argparse
import sys

import hyphalnet.hypha as hyp
from hyphalnet.hypha import hyphalNetwork
import hyphalnet.proteomics as prot
import hyphalnet.hyphEnrich as hyEnrich
import hyphalnet.hyphaeStats as hyStats
import pandas as pd
import synapseclient
from synapseclient import File, Table, table
import pickle
from scipy import stats

syn = synapseclient.Synapse()
syn.login()
data = syn.tableQuery("SELECT * FROM syn22172602").asDataFrame()


class kvdictAppendAction(argparse.Action):
    """
    argparse action to split an argument into KEY=VALUE form
    on the first = and append to a dictionary.
    """
    def __call__(self, parser, args, values, option_string=None):
        assert(len(values) == 1)
        values = values[0].split(',')
        for val in values:
            try:
                (k, v) = val.split("=", 2)
            except ValueError as ex:
                raise argparse.ArgumentError(self, \
                                             "could not parse argument \"{values[0]}\" as k=v format")
            d = getattr(args, self.dest) or {}
            d[k] = v
            setattr(args, self.dest, d)

#Parser information for command line
parser = argparse.ArgumentParser(description="""Get data from the proteomic \
                                 data commons and build community networks""")
parser.add_argument('--enrich', dest='doEnrich', action='store_true',\
                    default=False, help='Flag to do GO enrichment')
parser.add_argument('--saveGraphs', dest='toFile', action='store_true',\
                    default=False, help='Flag to save networks to file')
parser.add_argument('--getDistances', dest='getDist', action='store_true',\
                    default=False, help='Get and save distances')
parser.add_argument('--fromFile', dest='fromFile', nargs=1,\
                    action=kvdictAppendAction,metavar='KEY=VALUE',\
                    help='Key/value params for extra files')


#def tumor_genes(data_frame, group, subgroup, value):
#    significant_mut = data_frame[abs(data_frame[value]) > 0]
#    mutation_dictionary = (significant_mut.groupby(group).apply(lambda x: dict(zip(x[subgroup], x[value]))).to_dict())
#    return mutation_dictionary


def significant_prots(data_frame, group, subgroup, value,quantThresh=0.01):
    '''
    Gets top most expressed proteins in each patient
    '''
    pquants = pd.DataFrame({'thresh':data_frame.groupby(group)[value].quantile(1.0-quantThresh)})
    tdat = data_frame.merge(pquants, on=group)
    tdat = tdat.assign(topProt=tdat[value] > tdat['thresh'])
    significant = tdat[tdat['topProt']]

    gene_dictionary = (significant.groupby(group).apply(lambda x: dict(zip(x[subgroup], x[value]))).to_dict())
    return gene_dictionary

def loadFromFile(file_name_dict):
    hyphae = dict()
    for key, fname in file_name_dict.items():
        hyphae[key] = hyp.load_from_file(fname)
    return hyphae

def main():
    ''''
    Very basic scripts that pull the panCan hyphae and computes distances to store
    '''
    panCan = hyp.load_from_file(syn.get('syn22392951').path)
    g = panCan.interactome #interactome is already in the object!!!

    args = parser.parse_args()
    beta = 0.5
    proteomics_dictionary = significant_prots(data, 'AML sample', 'Gene', 'LogFoldChange')
    #gene_dictionary = tumor_genes(data, 'AML sample', 'Gene', 'Tumor VAF')

    hyphae = dict()
    hyphae['panCan'] = panCan
    #hyphae['mutations'] = hyphalNetwork(gene_dictionary, g.copy(), beta)
    hyphae['patients'] = hyphalNetwork(proteomics_dictionary, g.copy(), beta)

    #now compute graph distances to ascertain fidelity
    res = hyStats.compute_all_distances(hyphae)
    #res.to_csv('amlNetworkdistances.csv')
    tab = table.build_table("AML PanCan Network Distances", 'syn22128879',res)
    syn.store(tab)


if __name__ == '__main__':
    main()
