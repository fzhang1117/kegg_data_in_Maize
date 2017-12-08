## kegg_pathway_geneget.py
## Get all the genes in a pathway and trans them from Entrez ID to Gramene ID
## Zhang Fei
## 2017-11-14

import re
from Bio.KEGG import REST

def dic_cpd_build():
    dic = {}
    cpd_list = REST.kegg_list('compound')
    contains = cpd_list.read()
    cpd_list.close()
    contains = contains.strip('\n').split('\n')
    for line in contains:
        line = line.split('\t')
        key = line[0].split(':')[1]
        value = line[1]
        if dic.get(key) is None:
            dic[key] = value
    return dic

def compound_ask(request, my_dic):
    compound = REST.kegg_link('cpd' ,request)
    contain_cpd = compound.read()
    compound.close()
    contain_cpd = contain_cpd.strip('\n').split('\n')
    fl_out = './result/' + request + '.txt'
    fh_out = open(fl_out, 'a')
    formula, exact_mass, mol_weight = "-", "-", "-"

    title = ['ENTRY', 'NAME', 'FORMULA', 'EXACT_MASS', 'MOL_WEIGHT']
    fh_out.writelines('\t'.join(title))
    fh_out.writelines('\n')
    for line in contain_cpd:
        line = line.split('\t')
        if len(line) >= 2:
            cpd = line[1]
            entry = line[1][4: ]
            name = my_dic[entry]
            cpd_info = REST.kegg_get(cpd).read()
            cpd_info = cpd_info.strip('\n').split('\n')
            for line2 in cpd_info:
                if line2[0: 7] == 'FORMULA':
                    formula = ' '.join(line2.split()).split(' ')[1]
                if line2[0: 10] == 'EXACT_MASS':
                    exact_mass = ' '.join(line2.split()).split(' ')[1]
                if line2[0: 10] == 'MOL_WEIGHT':
                    mol_weight = ' '.join(line2.split()).split(' ')[1]
            output = [entry, name, formula, exact_mass, mol_weight]
            fh_out.writelines('\t'.join(output))
            fh_out.writelines('\n')
    fh_out.close()
    return 0

def pathway_lists():
    result = []
    path_list = REST.kegg_list("pathway/zma")
    contains = path_list.read()
    path_list.close()
    contains = contains.strip('\n')
    contains = contains.split('\n')
    for line in contains:
        line = line.split('\t')
        pathway = line[0].split(':')[1]
        result.append(pathway)
    return result

def main():
    zma_path = pathway_lists()
    my_dic = dic_cpd_build()
    for pathway in zma_path:
        request = 'map' + pathway[3: ]
        print request
        if request != 'map00196' and request != 'map00511' and request != 'map00514' and request != 'map00601' and request != 'map00603' and request != 'map00604' and request != 'map03008' and request != 'map03010':
            compound_ask(request, my_dic)
    return 0

main()
