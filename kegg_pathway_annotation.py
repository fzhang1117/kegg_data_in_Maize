## kegg_pathway_annotation.py ##
## Annotate kegg pathway ##
## Zhang Fei ##
## 2017-12-27 ##

import sys
from Bio.KEGG import REST

def dic_pathway_build():
    dic = {}
    pathway_list = REST.kegg_list('pathway', org = 'zma')
    contains = pathway_list.read()
    print contains
    pathway_list.close()

    contains = contains.strip('\n').split('\n')
    for line in contains:
        line = line.split('\t')
        key = line[0].split(':')[1]
        value = line[1]
        if dic.get(key) is None:
            dic[key] = value
    return dic

def pathway_ask(dic, fh_query, fh_out):
    for line in fh_query:
        line = line.strip('\r\n')
        if dic.get(line) is None:
            output = [line, '-']
            fh_out.writelines('\t'.join(output))
            fh_out.writelines('\n')
        else:
            output = [line, dic[line]]
            fh_out.writelines('\t'.join(output))
            fh_out.writelines('\n')
    return 0

def main():
    fh_query = open(sys.argv[1], 'r')
    fh_out = open(sys.argv[2], 'a')
    dic_pathway = dic_pathway_build()
    pathway_ask(dic_pathway, fh_query, fh_out)
    return 0

main()
