import os
import glob
import argparse
from Bio import SeqIO


############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-a",
                    help="ARG file of genome name", type=str, default='Genome_with_ARG.txt',
                    metavar='your ARG file of genome name')
parser.add_argument("-i",
                    help="Int file of genome name", type=str, default='Genome_with_IntI1.txt',
                    metavar='your Int file of genome name')
parser.add_argument("-g",
                    help="genome file of genome information", type=str,
                    default='assembly_summary_genbank_taxon.txt.normalized2',
                    metavar='your genome file of genome information')
parser.add_argument('--r',
                    default="Network", action='store', type=str, metavar='Network',
                    help="Set the result directory for edge and node table (default: Network)")


################################################## Definition ########################################################
args = parser.parse_args()
ARG_dir, ARG_file = os.path.split(os.path.abspath(args.a))
Int_dir, Int_file = os.path.split(os.path.abspath(args.i))
anno_dir, anno_file = os.path.split(os.path.abspath(args.g))
try:
    os.mkdir(args.r)
except OSError:
    pass


################################################### Function #########################################################
def node_add(node,shape,lable,lablecolor,color,length):
    f1=open(os.path.join(args.r,str(anno_file).replace('.txt','.node.txt')),'ab')
    f1.write(str(node)+'\t'+str(shape)+'\t'+str(lable)+'\t'+str(lablecolor)+'\t'+str(color)+'\t'+str(length)+'\n')


def edge_add(node1,node2,edgetype,color,length):
    f1 = open(os.path.join(args.r, str(anno_file).replace('.txt','.edge.txt')), 'ab')
    f1.write(str(node1) + '\t' + str(node2) + '\t' + str(edgetype) + '\t' + str(color)+ '\t' + str(length) + '\n')


def normalize(f1,f2):
    for line in f1:
        newline='\t'.join(str(line).split('\t')[0:3])
        for i in range(3,8+1):
            if str(line).split('\t')[i]!='NA' or str(line).split('\t')[i]!='':
                newline+='\t'+str(line).split('\t')[i]
            else:
                for j in reversed(range(3, i)):  # replace 'NA' or empty annotation
                    if str(line).split('\t')[j] !='NA' and str(line).split('\t')[j]!='':
                        newline += '\t' +str(line).split('\t')[j] + '_' + str(i)
                        break
        if str(line).split('\t')[9] =='NA' or str(line).split('\t')[9] =='' or 'environmental samples' in str(line).split('\t')[9] :
            node2=str(line).split('\t')[1]
        else:
            node2 =str(line).split('\t')[9]
        newline += '\t' + node2+'\t'+'\t'.join(str(line).split('\t')[10:len(str(line).split('\t'))])
        f2.write(newline)
    f1.close()
    f2.close()


def normalize2(f3):
    Taxon = dict()
    for lines in f3:
        for i in range(3, 7 + 1):
            node2 = str(lines).split('\t')[i + 1]
            node1 = str(lines).split('\t')[i]
            if node2 not in Taxon:
                Taxon.setdefault(node2, [[node1 + '--' + node2, 1]])
            else:
                Set = 0
                for key in Taxon[node2]:
                    if node1 + '--' + node2 == key[0]:
                        key[1] += 1
                        Set = 1
                        break
                if Set == 0:
                    Taxon[node2].append([node1 + '--' + node2, 1])
    Taxon2 = dict()
    for node2 in Taxon:
        newkey = ''
        max = 0
        for key in Taxon[node2]:
            if key[1] > max:
                max = key[1]
                newkey = key
        Taxon2.setdefault(node2, newkey)
    f3.close()
    return Taxon2


def normalize3(f3,f4,Taxon2):
    for line in f3:
        newline = '\t'.join(str(line).split('\t')[0:3])
        node2=str(line).split('\t')[8]
        templine=[str(line).split('\t')[9],node2]
        for i in reversed(range(3, 7+1)):
            node1=Taxon2[node2][0].split('--')[0]
            templine.append(node1)
            node2=node1
        for newnode in reversed(templine):
            newline+='\t'+newnode
        newline += '\t' + '\t'.join(str(line).split('\t')[10:len(str(line).split('\t'))])
        f4.write(newline)
    f3.close()
    f4.close()


def normalize4(edgefile):
    Taxon = dict()
    for lines in edgefile:
        node1 = str(lines).split('\t')[0]
        node2 = str(lines).split('\t')[1]
        if node2 not in Taxon:
            Taxon.setdefault(node2, [[node1 + '--' + node2, 1]])
        else:
            Set = 0
            for key in Taxon[node2]:
                if node1 + '--' + node2 == key[0]:
                    key[1] += 1
                    Set = 1
                    break
            if Set == 0:
                Taxon[node2].append([node1 + '--' + node2, 1])
    for node2 in Taxon:
        if len(Taxon[node2])>1:
            print node2
            print Taxon[node2]
    edgefile.close()


def network(line,nodetable,edgetable,Color,ARG,Int):
    lastnewlable = 0
    for i in range(3,9+1): #taxon nodes and edges
        if i ==3:
            if 'Bacteria' + '-' + str(line).split('\t')[i] not in edgetable:
                edge_add('Bacteria', str(line).split('\t')[i], 'taxon', 'None', 1)
                edgetable.append('Bacteria' + '-' + str(line).split('\t')[i])
        if str(line).split('\t')[i] not in Nalist:
            node1=str(line).split('\t')[i]
            if len(str(line).split('\t')[i].split(' '))>=2:
                lableshort=str(line).split('\t')[i].split(' ')[-2]+' '+str(line).split('\t')[i].split(' ')[-1]
            else:
                lableshort=str(line).split('\t')[i]
            if node1 not in nodetable and i < 9:
                node_add(node1, 'taxon', lableshort,'', 'white', 100) #node color as pathogen
                nodetable.append(node1)
        if i<=8:#add edge
            if str(line).split('\t')[i+1] not in Nalist or (i==8 and str(line).split('\t')[i + 1] in Nalist):
                node2=str(line).split('\t')[i+1]
                edgecolor='taxon'
                if i==8:
                    if (node2 =='NA' or 'environmental samples' in node2):
                        node2=str(line).split('\t')[0]+'_'+str(line).split('\t')[1]+'.sp'
                    else:
                        node2=str(line).split('\t')[0]+'_'+node2+'.sp'
                    ARGInt = ''
                    lablecolor = Color.get(str(line).split('\t')[-1].split('\r')[0].split('\n')[0], 'None')
                    if node2 not in nodetable:
                        if str(line).split('\t')[0] in ARG and str(line).split('\t')[0] in Int:
                            ARGInt = 'ARGInt'
                        elif str(line).split('\t')[0] in ARG and str(line).split('\t')[0] not in Int:
                            ARGInt = 'ARG'
                        elif str(line).split('\t')[0] not in ARG and str(line).split('\t')[0] in Int:
                            ARGInt = 'Int'
                    if node2 not in ARGIntlist:
                        ARGIntlist.setdefault(node2, [ARGInt, lablecolor])
                    elif ARGInt != '' and ARGIntlist[node2][0] != '':
                        if ARGInt != ARGIntlist[node2][0]:
                            ARGIntlist[node2] = ['ARGInt', lablecolor]
                    elif ARGInt != '' and ARGIntlist[node2][0] == '':
                        ARGIntlist[node2] = [ARGInt, lablecolor]
                    edgecolor=ARGInt
                if lastnewlable == 0:
                    if str(line).split('\t')[i] in Nalist:
                        node1=''
                        for j in reversed(range(3,i)):#add 'NA' node
                            if str(line).split('\t')[j] not in Nalist:
                                node1=str(line).split('\t')[j]+'_'+str(i)
                                if len(str(line).split('\t')[j].split(' ')) >= 2:
                                    lableshort = str(line).split('\t')[j].split(' ')[-2] + ' ' + \
                                                 str(line).split('\t')[j].split(' ')[-1]
                                else:
                                    lableshort = str(line).split('\t')[j]
                                node_add(node1, 'taxon', lableshort, '',
                                         'white', 100) #node color as pathogen
                                nodetable.append(node1)
                                edge_add(str(line).split('\t')[j], str(node1), 'taxon', '', 1)
                                edgetable.append(str(line).split('\t')[j] + '-' + str(node1))
                                break
                    elif node2==str(line).split('\t')[i]: # replicate taxon names add lable
                        newlable = i
                        node_add(node2 + str(newlable), 'taxon', lableshort, '', lablecolor,
                                 100)
                        node2 = node2 + str(newlable)
                        lastnewlable=newlable
                    else:
                        node1=str(line).split('\t')[i]
                else:
                    node1 = str(line).split('\t')[i]+str(lastnewlable)
                    lastnewlable=0
                if str(node1)+'-'+str(node2) not in edgetable and node1!='':
                    edge_add(str(node1), str(node2),'taxon',edgecolor, 1)
                    edgetable.append(str(node1)+'-'+str(node2))


def netword_add(Repli):
    for Repli_node in Repli:
        if Repli[Repli_node][0] > 1:
            node_add(Repli_node, 'Frequency', '*'+str(Repli[Repli_node][0]), 'yellow', 'Frequency', 500+Repli[Repli_node][0]*10)
            edge_add(Repli[Repli_node][1], Repli_node, 'Frequency', 'Frequency', 0)


################################################### Programme #######################################################
nodetable=[]
edgetable=[]
# output column names
node_add('Node','Shape','Lable','Lablecolor','Color','Length')
node_add('Bacteria','taxon','Bacteria','white','Bacteria',100)
edge_add('Node1','Node2','Type','Color','Length')
# set up edge colors for pathogen (1) and non-pathogen (0)
Color=dict()
Color.setdefault('0','white')
Color.setdefault('1','red')

# input genomes with ARGs
ARG=[]
for line in open(os.path.join(ARG_dir, ARG_file),'rb'):
    ARG.append(str(line).split('\r')[0].split('\n')[0])
# input genomes with integrases
Int=[]
for line in open(os.path.join(Int_dir, Int_file),'rb'):
    Int.append(str(line).split('\r')[0].split('\n')[0])


# normalize the taxonomy metadata
Nalist=['NA','environmental samples']
#normalize(open(os.path.join(anno_dir, anno_file),'rb'),open(os.path.join(anno_dir, anno_file+'.normalized'),'wb'))
#Taxon2=normalize2(open(os.path.join(anno_dir, anno_file + '.normalized'), 'rb'))
#normalize3(open(os.path.join(anno_dir, anno_file + '.normalized'), 'rb'),
#           open(os.path.join(anno_dir, anno_file),'wb'),Taxon2)


# set up the skeleton of phylogenetic tree
ARGIntlist=dict()
for line in open(os.path.join(anno_dir, anno_file),'rb'):
    if 'assembly_accession' not in str(line).split('\t')[0]:
        network(str(line),nodetable,edgetable,Color,ARG,Int)


# add the information of ARGs and integrases to the phylogenetic tree
for node1 in ARGIntlist:
    if node1 not in nodetable:
        node_add(node1, ARGIntlist[node1][0], '', '', ARGIntlist[node1][1], 150)  # node color as pathogen
        nodetable.append(node1)
        #node_add(node1 + ARGIntlist[node1], ARGIntlist[node1],
        #                     '', ''
        #                     , ARGIntlist[node1], 150)
        #edge_add(node1, node1 + ARGIntlist[node1], ARGIntlist[node1], '', 1)
normalize4(open(os.path.join(args.r, str(anno_file).replace('.txt','.edge.txt')), 'rb'))