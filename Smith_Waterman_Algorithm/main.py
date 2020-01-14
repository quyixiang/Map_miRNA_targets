import os
from itertools import islice
from smith_waterman import *
from quit import quit_func

header_l=[]
seq_l=[]
with open("sequence_miRNA.txt","r") as f:
    for line in f:
        if line[0]=='>':
            header=line.split(' ')
            header=str(header[0])
            header=header.strip('>')
            header_l.append(header)
        else:
            seq=str(line)
            seq=seq.strip('\n')
            seq_l.append(seq)

miRNA_dic=dict(zip(header_l,seq_l))

file_output= open("compare_result.csv","w")
file_output.write("miRNA" + "," + "mRNA" + "," + "result" + "\n")


def output(mrna_name,plus,minus,zero,ratio):
    fi=open("result_total.csv","a")
    fi.write(str(mrna_name)+","+str(plus)+","+str(minus)+","+str(zero)+","+str(ratio)+"\n")
    fi.close()

if __name__ == '__main__':
    prepare()
    
    for key in miRNA_dic:
        rna=[]
        rna_seq_l=[]
        path=os.getcwd()+'/result/'+str(key)
        filedir = path+'/sequence_result'
        if os.path.exists(filedir):
            glist=os.listdir(filedir)
            glist_new=[]
            for item in glist:
                filepath=filedir+'/'+item
                if os.stat(filepath).st_size!=0 and os.stat(filepath).st_size<1048576:
                    with open(filepath) as ff:
                        for line in islice(ff, 1, None):
                            line=line.strip()
                            rna.append(line)
                        rna_seq=''.join(rna)
                        if rna_seq.find("N")==-1:
                            rna_seq_l.append(rna_seq)
                            glist_new.append(item)
                else:
                    continue

            resultdir=path+'/result'
            if os.path.exists(resultdir):
                pass
            else:
                os.mkdir(resultdir)
            mi=miRNA_dic[key]
            miname=key.strip('hsa-')
            total_results={"plus":0,"minus":0,"0":0}
            #glist_new是每个miRNA中存在的mRNA名字
            for i in range(len(glist_new)):
                mR=rna_seq_l[i]
                mRname=glist_new[i]
                print(mRname)
                water_resolution(mi,mR,miname,mRname,resultdir)
                outp=analyze_results(miname,mRname,file_output)
                total_results[outp]+=1
                plus_ratio=total_results["plus"]/(total_results["plus"]+total_results["minus"]+total_results["0"])
            output(key,total_results["plus"],total_results["minus"],total_results["0"],plus_ratio)

        
        else:
            continue

file_output.close()