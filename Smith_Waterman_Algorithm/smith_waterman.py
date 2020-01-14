from __future__ import division, print_function
from compareResults import compare
import numpy as np

import time
import random


def prepare():
    global S
    # Read the score matrix from file
    S = np.genfromtxt('ScoreMatrix.csv',delimiter=',')
    global base_val
    base_val = {'AA':0,'AG':1,'AC':2,'AU':3,'GA':4,'GG':5,'GC':6,'GU':7,'CA':8,'CG':9,'CC':10,'CU':11,'UA':12,'UG':13,'UC':14,'UU':15}
    global pt
    pt ={'gap': -10}


def mch(alpha,beta):
    global S
    return S[base_val[alpha]][base_val[beta]]


def tracebk(i,j,s1,s2):
    # Traceback
    global align1,align2,H,T
    while T[i][j] != 0:
        if T[i][j] == 3:
            a1 = s1[i-1:i+1]
            a2 = s2[j-1:j+1]
            i -= 1
            j -= 1
        elif T[i][j] == 2:
            a1 = '--'
            a2 = s2[j-1:j+1]
            j -= 1
        elif T[i][j] == 1:
            a1 = s1[i-1:i+1]
            a2 = '--'
            i -= 1
        a1 = a1[::-1]
        a2 = a2[::-1]
        align1 += a1
        align2 += a2
    align1 = align1[::-1]
    align2 = align2[::-1]


def output_short(a,max_list_element):
    global align1, align2
    align1 = align1[::-1]
    align2 = align2[::-1]
    sym, al1, al2 = '', '', ''
    al1 += align1[0:2]
    al2 += align2[0:2]
    for i in range(3,len(align1),2):
        al1 += align1[i]
        al2 += align2[i]
    al1 = al1[::-1]
    al2 = al2[::-1]
    print('No.',a,"\ttarget:",al1,"\tmiRseq:",al2,'\tx:',max_list_element[1],'\tscore:',max_list_element[0])

def output_file(a,max_list_element,out_file):
    global align1, align2
    align1 = align1[::-1]
    align2 = align2[::-1]
    sym, al1, al2 = '', '', ''
    al1 += align1[0:2]
    al2 += align2[0:2]
    for i in range(3,len(align1),2):
        al1 += align1[i]
        al2 += align2[i]
    al1 = al1[::-1]
    al2 = al2[::-1]
    out_file.write('No.'+str(a)+"\ttarget:"+al1+"\tmiRseq:"+al2+'\tx:'+str(max_list_element[1])+'\tscore:'+str(max_list_element[0])+"\n")

def output_csv(mirna, mrna, result, f_out):
    print("output_csv")
    #f_out.write("miRNA" + "," + "mRNA" + "," + "result" + "\n")
    print("head")
    f_out.write(str(mirna) + "," + str(mrna) + "," + str(result) + "\n")

def output_full(a,max_score):
    # Print the result, to see the full result
    global align1, align2
    print("*" * 50)
    print("No.",a)
    print('\n')
    # print(align1)
    # print(align2,'\n')
    align1 = align1[::-1]
    align2 = align2[::-1]
    sym, al1, al2 = '', '', ''
    al1 += align1[0:2]
    al2 += align2[0:2]
    for i in range(3,len(align1),2):
        al1 += align1[i]
        al2 += align2[i]
    al1 = al1[::-1]
    al2 = al2[::-1]

    al2n = ''
    for k in range(len(al2)):
        complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}
        a2 = al2[k]
        if a2 == '-':
            al2n += '-'
        else:
            al2n += complement[a2]
        
    iden = 0
    for j in range(len(al1)):
        a1 = al1[j]
        a2 = al2n[j]
        if a1 == a2:                
            sym += '|'
            iden += 1
        elif (a1,a2) in [('G','A'),('U','C'),('U','T')]:
            sym += ':'
        else: 
            sym += ' '

    identity = iden / len(al1) * 100
    print('Identity = %f percent' % identity)    
    print('Score =', max_score)
    print('\n')
    print(al1)
    print(sym)
    print(al2)
    print('\n')

    
def ptmax(H,x,i,j):
    for (i2,j2) in [(i-1,j),(i+1,j)]:
        try:
            if x < H[i2][j2]:
                return False
        except IndexError:
            pass
    return True


def water(s1, s2, f_parameter = 10.0, out_file = None):
    m, n = len(s1), len(s2)
    global H,T
    H = np.zeros((m, n))
    T = np.zeros((m, n))    

    max_list = []
    # Score, Pointer Matrix
    for i in range(1, m):
        for j in range(1, n):
            if j in range(2,9):                
                sc_diag = H[i-1][j-1] + f_parameter*mch(s1[i-1:i+1], s2[j-1:j+1])
                sc_up = H[i][j-1] + f_parameter*pt['gap']
                sc_left = H[i-1][j] + f_parameter*pt['gap']
                H[i][j] = max(sc_left,sc_up,sc_diag)
            else:
                sc_diag = H[i-1][j-1] + mch(s1[i-1:i+1], s2[j-1:j+1])
                sc_up = H[i][j-1] + pt['gap']
                sc_left = H[i-1][j] + pt['gap']
                H[i][j] = max(sc_left, sc_up, sc_diag)
            trace_value = {0:0 , sc_left:1 , sc_up:2 , sc_diag:3}
            T[i][j] = trace_value[H[i][j]]
    for i in range(1, m):
        for j in range(1, n):
            if ptmax(H,H[i][j],i,j) == True and H[i][j] > min(m,n)*45.0 and j == n-1:
                max_list=max_list+[(H[i,j],i,j)]
    max_list = sorted(max_list, key = lambda v:v[0], reverse=True)
    
    global align1, align2 

    if out_file != None:
        for a in range(len(max_list)):
            align1, align2 = '', ''
            i,j=max_list[a][1],max_list[a][2]
            tracebk(i,j,s1,s2)
            output_file(a+1,max_list[a],out_file)
    else:     
        for a in range(len(max_list)):
            align1, align2 = '', ''
            i,j=max_list[a][1],max_list[a][2]
            tracebk(i,j,s1,s2)
            output_short(a+1,max_list[a])
            #output_full(a+1,max_list[a][0])


def water_shuffle(ref,qu,mRname,miname,dir,repeat=50):
    global iterate_shuffle_number
    global result_target, result_shuffle

    iterate_shuffle_number = 0
    result_target=dir+'/rt_'+miname+'_'+mRname+'.txt'
    fout = open(result_target,"w")
    water(ref,qu,out_file = fout)
    fout.close()

    result_shuffle=dir+'/rs_'+miname+'_'+mRname+'.txt'
    fout = open(result_shuffle,"w")
    shuff_seq = qu
    print("intial miRNA",shuff_seq)
    for count in range(repeat):
        shuff_seq = ''.join(random.sample(list(qu),len(qu)))
        iterate_shuffle_number += 1
        print("No.",iterate_shuffle_number)
        print("New_Seq.",shuff_seq)
        print("progress:{0}%".format(round((count + 1) * 100 / repeat)), end="\r")
        qu = shuff_seq
    water(ref,shuff_seq,out_file = fout)
    fout.close()


def water_resolution(mi,mR,miname,mRname,dir):
    start = time.clock()

    qu = (mi.upper()).replace('T','U')
    ref = (mR.upper()).replace('T','U')
    water_shuffle(ref,qu,mRname,miname,dir)
    

    print("Time used:",(time.clock() - start),'second\n')

def analyze_results(miname,mRname,file_output):
    print("analyze_results")
    result_compare = compare(result_target, result_shuffle)
    output_csv(miname, mRname, result_compare,file_output)
    return result_compare

def quit_func():
    quit_info = "\n<Press Enter to Quit>"
    print(quit_info.center(20))
    raw_input()
    quit()