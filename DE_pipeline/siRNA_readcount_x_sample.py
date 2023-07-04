from operator import add
import pandas as pd
import collections

def spot_doubles(dict4):
        groups=[]
        dict1=dict(dict4)
        for b in list(dict1.values()):
                el=[k for k,v in dict1.items() if collections.Counter(v) == collections.Counter(b) ]
                for r in el:
                        del dict1[r]
                if el:
                        groups.append(el)
        return groups

with open("siRNA_read_count.tab") as r:
        file=[a.strip() for a in r]
dict3={}

for a in file:
        dict3[a.split('\t')[0]]=a.split('\t')[1:7]



df = pd.DataFrame()

for ll in range(20,31):
        b1=ll-2
        print(str(ll))
        dict2={}
        with open("DaRe_"+str(ll)+"_overlap_"+str(b1)+"_table_fast.fasta") as r:
                file2=[a.strip() for a in r.readlines() if a[0] == ">"]
        with open("DaRe_"+str(ll)+"_overlap_"+str(b1)+"_table_fast.fasta") as r:
                file2_1=[a.strip() for a in r.readlines() if a[0] != ">"]

        for a in range(0,len(file2)):
                if file2[a].split("|")[0]+file2[a].split("|")[1] in dict2.keys():
                        dict2[file2[a].split("|")[0]+file2[a].split("|")[1]]+=[file2_1[a]]
                else:
                        dict2[file2[a].split("|")[0]+file2[a].split("|")[1]]=[file2_1[a]]
        gp=spot_doubles(dict2)
        print(dict2.keys())
        gp2=[a[0] for a in gp ]
        for c in dict2.keys():
                if len(dict2[c]) > 1:
                        b=[]            
                        for a in dict2[c]:
                                b+=list(map(int, dict3[a])),
                        tot=[sum(i) for i in zip(*b)]
                else:
                        tot=list(map(int, dict3[dict2[c][0]]))
                if c in gp2:
                    new_row={"location":c+"_"+str(ll)+"_siRNA", "sp1":tot[0], "sp2":tot[1], "sp3":tot[2], "sp4":tot[3], "sp5":tot[4], "sp6":tot[5] }
                    df=df.append(new_row, ignore_index=True)
        with open("siRNA_doubles.tab","a+") as r:
                r.writelines('\t'.join([a+"_"+str(ll)+"_siRNA" for a in ciao]) + '\n' for ciao in gp )


print(df)
df.to_csv("siRNA_readcount_x_sample.tab",sep="\t",index=False)
