## kegg_pathway_geneget.py
## Get all the genes in a pathway and trans them from Entrez ID to Gramene ID
## Zhang Fei
## 2017-11-14

import re, pickle
from Bio.KEGG import REST
from Bio import SeqIO
from Bio.KEGG import Enzyme
from Bio.KEGG.KGML import KGML_parser

def my_dic_build():
    fh_dic = open("./Zea_mays_gene_info.dic", 'rb')
    my_dic = pickle.load(fh_dic)
    fh_dic.close()
    return my_dic

def pathway_ask(request, my_dic):
    pathway = REST.kegg_get(request)
    contains = pathway.read()
    pathway.close()
    outfl = './result/genes_' + request + '.txt'
    fh_output = open(outfl, 'a')
    contains = contains.split('\n')
    for line in contains:
        if line[0: 4] == 'NAME':
            print line
        if line[0: 11] == 'DESCRIPTION':
            print line
        if line[0: 4] == 'GENE' and re.search('\[EC', line) is not None:
            fh_output.writelines('\t'.join(['Entrez_id', 'Maize_v3', 'Maize_v4', 'KO', 'EC', 'Gene_description']))
            fh_output.writelines('\n')
            KO = line.split('[KO:')[1].split(']')[0]
            EC = 'EC' + line.split('[EC:')[1].split(']')[0]
            Entrez_id = ' '.join(line.split()).split(' ')[1]
            if my_dic.get(Entrez_id) is not None:
                Gene_v3 = my_dic[Entrez_id][2]
                Gene_v4 = my_dic[Entrez_id][1]
                Gene_description = my_dic[Entrez_id][7]
            else:
                Gene_v3, Gene_v4, Gene_description = 'not found', 'not found', 'not found'
            output = [Entrez_id, Gene_v3, Gene_v4, KO, EC, Gene_description]
            print '\t'.join(output)
            fh_output.writelines('\t'.join(output))
            fh_output.writelines('\n')
        elif line[0: 4] == 'GENE' and re.search('\[EC', line) is None:
            fh_output.writelines('\t'.join(['Entrez_id', 'Maize_v3', 'Maize_v4', 'KO', 'EC', 'Gene_description']))
            fh_output.writelines('\n')
            KO = line.split('[KO:')[1].split(']')[0]
            EC = '-'
            Entrez_id = ' '.join(line.split()).split(' ')[1]
            if my_dic.get(Entrez_id) is not None:
                Gene_v3 = my_dic[Entrez_id][2]
                Gene_v4 = my_dic[Entrez_id][1]
                Gene_description = my_dic[Entrez_id][7]
            else:
                Gene_v3, Gene_v4, Gene_description = 'not found', 'not found', 'not found'
            output = [Entrez_id, Gene_v3, Gene_v4, KO, EC, Gene_description]
            print '\t'.join(output)
            fh_output.writelines('\t'.join(output))
            fh_output.writelines('\n')

        if re.search('\[KO:', line) is not None and line[0: 4]!= 'GENE' and re.search('\EC:', line) is not None:
            KO = line.split('[KO:')[1].split(']')[0]
            EC = 'EC' + line.split('[EC:')[1].split(']')[0]
            Entrez_id = ' '.join(line.split()).split(' ')[0]
            if my_dic.get(Entrez_id) is not None:
                Gene_v3 = my_dic[Entrez_id][2]
                Gene_v4 = my_dic[Entrez_id][1]
                Gene_description = my_dic[Entrez_id][7]
            else:
                Gene_v3, Gene_v4, Gene_description = 'not found', 'not found', 'not found'
            output = [Entrez_id, Gene_v3, Gene_v4, KO, EC, Gene_description]
            print '\t'.join(output)
            fh_output.writelines('\t'.join(output))
            fh_output.writelines('\n')
        elif re.search('\[KO:', line) is not None and line[0: 4]!= 'GENE' and re.search('\EC:', line) is None:
            KO = line.split('[KO:')[1].split(']')[0]
            EC = '-'
            Entrez_id = ' '.join(line.split()).split(' ')[0]
            if my_dic.get(Entrez_id) is not None:
                Gene_v3 = my_dic[Entrez_id][2]
                Gene_v4 = my_dic[Entrez_id][1]
                Gene_description = my_dic[Entrez_id][7]
            else:
                Gene_v3, Gene_v4, Gene_description = 'not found', 'not found', 'not found'
            output = [Entrez_id, Gene_v3, Gene_v4, KO, EC, Gene_description]
            print '\t'.join(output)
            fh_output.writelines('\t'.join(output))
            fh_output.writelines('\n')
    fh_output.close()
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
    my_dic = my_dic_build()
    for pathway in zma_path:
        pathway_ask(pathway, my_dic)
    return 0

main()
"""
# by open a kgml object, we can extract genes and compounds in a pathway

request = REST.kegg_get('zma00020/kgml')
contains = KGML_parser.read(request, 'r')

for gene in contains.maps:
    print gene

print(len(contains.entries))
print(len(contains.reactions))
print(len(contains.genes))
print contains
"""

"""
fh_out = open('zma00020.gb', 'w')
fh_out.write(request.read())
fh_out.close

contains = SeqIO.read("zma00020.gb", "genbank")

print contains
"""
#contains = contains.strip('\n')
#contains = contains.split('\n')
#print contains
"""
fh_test = open("test.txt", 'a')
fh_test.write(contains)
fh_test.close()
"""
