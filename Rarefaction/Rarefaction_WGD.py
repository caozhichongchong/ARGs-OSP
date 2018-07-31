import glob
import os
import random
import argparse


############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-r',
                    default=0, action='store', type=int, metavar='0 to 10',
                    help="number of reps")
parser.add_argument('-e',
                    default=0, action='store', type=int, metavar='0 or 1',
                    help="0 for false, 1 for true")


################################################## Definition ########################################################
args = parser.parse_args()
# total number of genomes
Totalsp=54718
# initialization
if args.e==0:
    try:
        os.mkdir('rarefaction')
    except OSError:
        pass
# input ARG list
    ARGlist=dict()
    for lines in open('ARGlist.txt','rb'):
        ARGlist.setdefault(lines.split('\t')[0],int(lines.split('\t')[-1].split('\r')[0].split('\n')[0]))
    ARGtable=dict()
    for lines in open('WGD_ARG.mothertable.rarefaction','rb'):
        try:
            spnum=int(float(lines.split('\t')[-1].split('\r')[0].split('\n')[0]))
            ARGnum=int(float(lines.split('\t')[-2]))
            if spnum not in ARGtable:
                ARGtable.setdefault(spnum,[ARGnum])
            else:
                ARGtable[spnum].append(ARGnum)
        except ValueError:
            pass
    TotalARG = 0
    ARGin = []
    i=0
    Raretable3 = range(0, 100)
    reps=args.r
    f1 = open('rarefaction/ARG_rarefaction.WGD.' + str(reps) + '.txt', 'wb')
    # main calculation
    while Raretable3!=[]:
        Randomrange = random.choice(Raretable3)
        if Randomrange in Raretable3:
            Raretable3.remove(Randomrange)
            Raretable = range(int((float(Randomrange))/100.0*float(Totalsp)+1.0),
                              int((float(Randomrange)+1.0)/100.0*float(Totalsp))+1)
            print [reps, Randomrange, int((float(Randomrange))/100.0*float(Totalsp)+1.0),
                              int((float(Randomrange)+1.0)/100.0*float(Totalsp))]
            while Raretable!=[] and TotalARG < 4049:
                Randomseq = random.choice(Raretable)
                if Randomseq in Raretable:
                    Raretable.remove(Randomseq)
                    i += 1
                    if Randomseq not in ARGtable:
                        pass
                    else:
                        ARGset = ARGtable[Randomseq]
                        for ARG in ARGset:
                            if ARG not in ARGin:
                                TotalARG += 1
                                ARGin.append(ARG)
                            else:
                                pass
                    f1.write(str(i) + '\t' + str(TotalARG) + '\n')
    f1.close()
else:
    # merge the results
    Raretable=0
    Raretable2=0
    f1 = open('rarefaction/ARG_rarefaction.WGD.average.txt', 'ab')
    for k in range(0,Totalsp,1):
        i+=1
        Num=0
        for reps in range(0,10,1):
            os.system('sed -n '+str(i)+','+str(i)+'p rarefaction/ARG_rarefaction.WGD.' + str(reps) + '.txt > rarefaction/temp.txt')
            for lines in open('rarefaction/temp.txt', 'rb'):
                try:
                    Num+=float(lines.split('\t')[-1].split('\n')[0])
                except ValueError:
                    print lines
        f1.write(str(i) + '\t' + str(float(Num)/10.0) + '\n')
    f1.close()
    
