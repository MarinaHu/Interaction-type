import numpy as np
import os
import argparse

#Possible chains in the protein
chain = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N']
#Molecule types include ligand, solvent and residues from the protein. A residue have a main chain (mc) and a side chian (sc).
restype = ['mc','sc','ligand','solvent']
#The atoms included in the main chain residue
mc = ['N', 'CA', 'C', 'O', 'H', 'HA', 'OXT', 'HA2', 'HA3', 'H?', 'W']
#Interaction types between the ligand and the residues: wc, wide contact with an interactomic gap distance >= 0.25 Angstrom; cc, close contact with an interactomic gap distance < 0.25 Angstrom; bo, big overlap with overlapping van der Waal radii >= 0.4 Angstrom; so, small overlap with overlapping van der Waal radii < 0.4 Angstrom; hb, hydrogen bond; all, include all the interaction types mentioned previously.
intypes = ['wc','cc','bo','so','hb','all']

#Organize the format for the seed
def gen_seed(seed_raw):
    seed_raw = seed_raw.replace(';',' ')
    seed_raw = seed_raw.replace(':',' ')
    seed_raw = seed_raw.split()
    seed = []
    for i in range(len(seed_raw)//2):
        seed_id = []
        seed_id = seed_raw[2*i+1].replace(',',' ')
        seed_id = seed_id.split()
        for j in range(len(seed_id)):
            seed.append(seed_raw[2*i] + '-' + seed_id[j])
    return seed

#Generating res_atom.dat files based on interaction types in the order of the interaction counting
def gen_int(filename, seed, solvent): 
    with open(filename, 'r') as fp: 
        data = fp.readlines()
        res_chain_dict = {}
        res_atom_dict = {}
        res_intfreq_dict = {}
        solv = []
        for line in data:
            intype = line[6:8].replace(' ','')
            chain1 = line[9:11].replace(' ','')
            res1 = int(line[11:16].strip())
            code1 = line[16:20].replace(' ','')
            atom1 = line[20:25].replace(' ','')
            chain2 = line[25:28].replace(' ','')
            chain2 = chain2.replace(':','')
            res2 = line[28:33].replace(' ','')
            res2 = res2.replace(chain2,'')
            res2 = int(res2)
            code2 = line[32:36].replace(' ','')
            atom2 = line[37:41].replace(' ','')
            a = chain1+'-'+ str(res1)
            b = chain2+'-'+ str(res2)
            if (a in seed) and (b not in seed):
                res_chain_dict[res1] = chain1
                res_chain_dict[res2] = chain2
                if code2 == solvent:
                    if res2 not in solv:
                        solv.append(res2)                                     
                if a not in res_intfreq_dict.keys():
                    res_atom_dict[a] = [[[],[],[],[],[],[]],[[],[],[],[],[],[]]]
                    res_intfreq_dict[a] = [[0,0,0,0,0,0],[0,0,0,0,0,0]]
                    for i in range(len(intypes)-1):
                        if intype == intypes[i]:
                            res_atom_dict[a][0][i].append(atom1)
                            res_atom_dict[a][0][5].append(atom1)
                            res_intfreq_dict[a][0][i] += 1
                            res_intfreq_dict[a][0][5] += 1                         
                else:
                    for i in range(len(intypes)-1):
                        if intype == intypes[i]:
                            res_intfreq_dict[a][0][i] += 1
                            res_intfreq_dict[a][0][5] += 1                                
                            if atom1 not in res_atom_dict[a][0][i]:
                                res_atom_dict[a][0][i].append(atom1)
                                if atom1 not in res_atom_dict[a][0][5]:
                                    res_atom_dict[a][0][5].append(atom1)
                if b not in res_intfreq_dict.keys():
                    res_atom_dict[b] = [[[],[],[],[],[],[]],[[],[],[],[],[],[]]]
                    res_intfreq_dict[b] = [[0,0,0,0,0,0],[0,0,0,0,0,0]]
                    if (res2 in solv) or (atom2 in mc):
                        for i in range(len(intypes)-1):
                            if intype == intypes[i]:
                                res_intfreq_dict[b][0][i] += 1
                                res_atom_dict[b][0][i].append(atom2)
                                res_atom_dict[b][0][5].append(atom2)
                                res_intfreq_dict[b][0][5] += 1
                    else:
                        for i in range(len(intypes)-1):
                            if intype == intypes[i]:
                                res_intfreq_dict[b][1][i] +=1
                                res_atom_dict[b][1][i].append(atom2)
                                res_atom_dict[b][1][5].append(atom2) 
                                res_intfreq_dict[b][1][5] += 1   
                else:
                    if (res2 in solv) or (atom2 in mc):
                        for i in range(len(intypes)-1):
                            if intype == intypes[i]:
                                res_intfreq_dict[b][0][i] += 1
                                res_intfreq_dict[b][0][5] += 1 
                                if atom2 not in res_atom_dict[b][0][i]:
                                    res_atom_dict[b][0][i].append(atom2)
                                    if atom2 not in res_atom_dict[b][0][5]:
                                        res_atom_dict[b][0][5].append(atom2)
                    else:
                        for i in range(len(intypes)-1):
                            if intype == intypes[i]:
                                res_intfreq_dict[b][1][i] += 1
                                res_intfreq_dict[b][1][5] += 1 
                                if atom2 not in res_atom_dict[b][1][i]:
                                    res_atom_dict[b][1][i].append(atom2)
                                    if atom2 not in res_atom_dict[b][1][5]:
                                        res_atom_dict[b][1][5].append(atom2)
            if (a not in seed) and (b in seed):
                res_chain_dict[res1] = chain1
                res_chain_dict[res2] = chain2
                if code1 == solvent:
                    if res1 not in solv:
                        solv.append(res1)                                     
                if b not in res_intfreq_dict.keys():
                    res_atom_dict[b] = [[[],[],[],[],[],[]],[[],[],[],[],[],[]]]
                    res_intfreq_dict[b] = [[0,0,0,0,0,0],[0,0,0,0,0,0]]
                    for i in range(len(intypes)-1):
                        if intype == intypes[i]:
                            res_atom_dict[b][0][i].append(atom2)
                            res_atom_dict[b][0][5].append(atom2)
                            res_intfreq_dict[b][0][i] += 1
                            res_intfreq_dict[b][0][5] += 1                         
                else:
                    for i in range(len(intypes)-1):
                        if intype == intypes[i]:
                            res_intfreq_dict[b][0][i] += 1
                            res_intfreq_dict[b][0][5] += 1                                
                            if atom2 not in res_atom_dict[b][0][i]:
                                res_atom_dict[b][0][i].append(atom2)
                                if atom2 not in res_atom_dict[b][0][5]:
                                    res_atom_dict[b][0][5].append(atom2)
                if a not in res_intfreq_dict.keys():
                    res_atom_dict[a] = [[[],[],[],[],[],[]],[[],[],[],[],[],[]]]
                    res_intfreq_dict[a] = [[0,0,0,0,0,0],[0,0,0,0,0,0]]
                    if (res1 in solv) or (atom1 in mc):
                        for i in range(len(intypes)-1):
                            if intype == intypes[i]:
                                res_intfreq_dict[a][0][i] += 1
                                res_atom_dict[a][0][i].append(atom1)
                                res_atom_dict[a][0][5].append(atom1)
                                res_intfreq_dict[a][0][5] += 1
                    else:
                        for i in range(len(intypes)-1):
                            if intype == intypes[i]:
                                res_intfreq_dict[a][1][i] +=1
                                res_atom_dict[a][1][i].append(atom1)
                                res_atom_dict[a][1][5].append(atom1) 
                                res_intfreq_dict[a][1][5] += 1   
                else:
                    if (res1 in solv) or (atom1 in mc):
                        for i in range(len(intypes)-1):
                            if intype == intypes[i]:
                                res_intfreq_dict[a][0][i] += 1
                                res_intfreq_dict[a][0][5] += 1 
                                if atom1 not in res_atom_dict[a][0][i]:
                                    res_atom_dict[a][0][i].append(atom1)
                                    if atom1 not in res_atom_dict[a][0][5]:
                                        res_atom_dict[a][0][5].append(atom1)
                    else:
                        for i in range(len(intypes)-1):
                            if intype == intypes[i]:
                                res_intfreq_dict[a][1][i] += 1
                                res_intfreq_dict[a][1][5] += 1 
                                if atom1 not in res_atom_dict[a][1][i]:
                                    res_atom_dict[a][1][i].append(atom1)
                                    if atom1 not in res_atom_dict[a][1][5]:
                                        res_atom_dict[a][1][5].append(atom1)            
    
    return res_chain_dict, res_atom_dict, res_intfreq_dict, solv

def w_res_atom_int(res_chain_dict, res_atom_dict, res_intfreq_dict, solv, seed):
    res_chain_freqall_dict = {}
    for k,v in res_chain_dict.items():
        x = v + '-' + str(k)
        res_chain_freqall_dict[k] = [v, (res_intfreq_dict[x][0][5] + res_intfreq_dict[x][1][5])]  
    sorted_res_chain_freqall_dict_by_freq = sorted(res_chain_freqall_dict.items(), key=lambda x:x[1][1], reverse=True)
    res_chain_freqall_dict = dict(sorted_res_chain_freqall_dict_by_freq) 
    keys = list(res_chain_freqall_dict.keys())
    for i in keys:
        keyname = res_chain_dict[i] +'-'+ str(i) 
        for j in range(len(intypes)): 
            with open('res_atom_dict_' + intypes[j] +'.dat', 'a') as f:
                if keyname in seed:
                    if res_intfreq_dict[keyname][0][j] != 0: 
                        f.write('%-4s'%res_chain_dict[i] + '%-8s'%i + '%-8s'%res_intfreq_dict[keyname][0][j])
                        for h in range(len(res_atom_dict[keyname][0][j])): 
                            f.write('%-6s'%str(res_atom_dict[keyname][0][j][h]))
                        f.write('\n') 
                    else:
                        pass
                elif i in solv:
                    if res_intfreq_dict[keyname][0][j] != 0: 
                        f.write('%-4s'%res_chain_dict[i] + '%-8s'%i+ '%-8s'%res_intfreq_dict[keyname][0][j])
                        for h in range(len(res_atom_dict[keyname][0][j])): 
                            f.write('%-6s'%str(res_atom_dict[keyname][0][j][h]))
                        f.write('\n') 
                else:
                    if (res_intfreq_dict[keyname][0][j] != 0) or (res_intfreq_dict[keyname][1][j] != 0):
                        sum = res_intfreq_dict[keyname][0][j] + res_intfreq_dict[keyname][1][j]  
                        if res_atom_dict[keyname][0][j] != 0: 
                            f.write('%-4s'%res_chain_dict[i]+ '%-8s'%i+ '%-8s'%str(sum)) 
                            for h in range(len(res_atom_dict[keyname][0][j])): 
                                f.write('%-6s'%res_atom_dict[keyname][0][j][h])
                        else:
                            f.write('%-4s'%res_chain_dict[i]+ '%-8s'%i+ '%-8s'%str(sum)) 
                        if res_atom_dict[keyname][1][j] != 0:
                            for h in range(len(res_atom_dict[keyname][1][j])): 
                                f.write('%-6s'%res_atom_dict[keyname][1][j][h])
                            f.write('\n')
                        else:
                            f.write('\n') 
    return

def w_res_atom_freq(res_chain_dict, res_intfreq_dict, solv, seed):
    res_chain_freqall_dict = {}
    for k,v in res_chain_dict.items():
        x = v + '-' + str(k)
        res_chain_freqall_dict[k] = [v, (res_intfreq_dict[x][0][5] + res_intfreq_dict[x][1][5])]  
    sorted_res_chain_freqall_dict_by_freq = sorted(res_chain_freqall_dict.items(), key=lambda x:x[1][1], reverse=True)
    res_chain_freqall_dict = dict(sorted_res_chain_freqall_dict_by_freq) 
    keys = list(res_chain_freqall_dict.keys())
    for i in keys:
        keyname = res_chain_dict[i] +'-'+ str(i) 
        with open ('res_intfreq_dict.dat', 'a') as f:
            if keyname in seed:
                f.write('%-8s'%keyname) 
                for j in range(len(res_intfreq_dict[keyname][0])):
                    f.write('%-8s'%str(res_intfreq_dict[keyname][0][j]))
                f.write('\n')
            elif i in solv:
                f.write('%-8s'%keyname) 
                for j in range(len(res_intfreq_dict[keyname][0])):
                    f.write('%-8s'%str(res_intfreq_dict[keyname][0][j]))
                f.write('\n')
            else:
                f.write('%-8s'%keyname) 
                for j in range(len(res_intfreq_dict[keyname][0])):
                    f.write('%-8s'%str(res_intfreq_dict[keyname][0][j]))
                f.write('\n')
                f.write('%-8s'%keyname) 
                for j in range(len(res_intfreq_dict[keyname][1])):
                    f.write('%-8s'%str(res_intfreq_dict[keyname][1][j]))
                f.write('\n')
    return





if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Trim large PDB file according to res_atoms.dat, write trimmed pdb in working directory')
    parser.add_argument('-probe', dest='probe', default = 'None', help='protonated pdbfile')
    parser.add_argument('-s', dest='seed', default = 'None', help='Chain:Resid1,Resid2;Chain:Resid1,Resid2')
    parser.add_argument('-slv', dest='solvent', default = 'None', help='water is represented by HOH')
    args = parser.parse_args()

    probe = args.probe
    seed = args.seed
    solvent = args.solvent
    
    seed = gen_seed(seed) 
    res_chain_dict, res_atom_dict, res_intfreq_dict, solv = gen_int(probe, seed, solvent)
    w_res_atom_int(res_chain_dict, res_atom_dict, res_intfreq_dict, solv, seed)
    w_res_atom_freq(res_chain_dict, res_intfreq_dict, solv, seed) 