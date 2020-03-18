"""
Updated on 03/18/20
@author: Tatiana Lenskaia
"""

import time as tm
import sys
import os

if len(sys.argv) == 2:
    path = sys.argv[1]+"\\"
    cur_dir = sys.argv[1] + " directory"
else:
    path = ""
    cur_dir = "the current directory"



if not os.path.exists(path+'Output_files'):
    os.makedirs(path+'Output_files')



def CheckInputFiles(path):
    if os.path.exists(path+"list.txt"):
        fList = open(path+"list.txt")
        
        flag = 0
        
        for line in fList:
            line = line.strip()
            if line != "":
                if (os.path.exists(path+"Input_files\\FASTA\\"+line+".fasta")==False) or (os.path.exists(path+"Input_files\\CRISPRCasdb\\"+line+".csv") == False):
                    flag = 1
                    break
        fList.close()
        
        if flag == 1:
            print("Not all input files are supplied")
            print("Please check that .fasta and .csv files are provided for each replicon accession number from list.txt")
    else:
        
        
        print("The list of replicons is missing in", cur_dir)
        print("Please provide a list of replicon accession numbers (list.txt)")
    return
CheckInputFiles(path)

#=============================================

def CountGC(seq):
    seq = seq.lower();
    n_seq = len(seq)
    
    n_a = seq.count("a");
    n_c = seq.count("c");
    n_g = seq.count("g");
    n_t = seq.count("t");
    
    n_bases = n_a + n_c + n_g + n_t
    
    n_other = n_seq - n_bases;
    
    if n_bases != 0:
        gc = round(1.0*(n_c+n_g)/n_bases*100,2)
    else:
        gc = "N/A"
    
    return [gc, n_other]



def GetListFromFile(fInName):
    t = []
    fIn = open(fInName, "r")
    for line in fIn:
        line = line.strip()
        if line != "":
            if line not in t:
                t.append(line)
    #print(len(t))
    fIn.close()
    return t


def GetDictFromFile(fInName, sep, header):
    fIn = open(fInName,"r")
    lines = fIn.readlines()
    if header == 1:
        header_line = lines[0]
        lines = lines[1:]
    
    d = {}
    t = []
    
    n = len(lines[0].split(sep))
    
    for line in lines:
        line = line.strip()
        if line != "":
            t_line = line.split(sep)
            if len(t_line) != n:
                print("Check format!", line, t_line)
            else:
                if t_line[0] not in d:
                    t.append(t_line[0])
                    d[t_line[0]] = t_line[1:]
                else:
                    print("First colum has non-unique values!")
    fIn.close()
    return [t,d]


#Updated: May 4, 2019    
def GetText(finName):
    '''Extracts text from a single fasta file'''
    fin = open(finName, 'r')
    text = ''
    for line in fin:
        line = line.strip()
        if (line != "") and (line[0] != '>'):
            text = text + line
    fin.close()
    return text



def Rev_cmp_upper(st):
    '''Computes the reverse complement of a string'''
    st = st.upper()
    cmp_st = st.translate(str.maketrans("ACGT","TGCA"))
    rev_cmp_st = cmp_st[::-1]
    return rev_cmp_st;

def CreateDictLocD_upper(text, m, gtp = "l"):
    if (m <= 0) or (m > len(text)):
        print("(n = "+str(m)+") is not a valid window size for this genome!!!");
        return {}
	
    d_g = dict()
    nn = len(text);
    gtype = gtp.lower();
    gtype = gtype[0];
	
    if gtype == "c":
        text = text + text[0:(m-1)];
        lastpos = nn;
    elif gtype == "l":
        lastpos = nn-m+1;
    else:
        print("Is this genome linear or circular?");
        return d_g;
		
    for ii in range (lastpos):
        bl = text[ii:(ii+m)]; 
        bl = bl.upper();
        if bl in d_g :
            dd = d_g[bl];
            dd[ii] = ""
            d_g[bl] = dd;
        else:
            dd = {}
            dd[ii] = ""
            d_g[bl] = dd;
    return d_g;	


def SearchOrganism_plm_cas(fRname, name,kind, count, ddd, d_rep_sps, d_fname):
    '''Analyzes self-targeting events in a given organism'''
    ct_chr = 0
    ct_plm = 0
    ct_other = 0
       
    d_sps = {}

    for seq in ddd:
        
        if "chromo" in ddd[seq].lower():
            ct_chr = ct_chr + 1
        elif "plasmid" in ddd[seq].lower():
            ct_plm = ct_plm + 1
        else:
            ct_other = ct_other+1
        
        
        if seq in d_rep_sps:
 
            for sp in d_rep_sps[seq]:
                nn = len(sp)
                if nn not in d_sps:
                    dd = {}
                    dd[sp] = [seq]
                    d_sps[nn] = dd
                else:
                    dd = d_sps[nn]
                    if sp in dd:
                        dd[sp].append(seq)
                    else:
                        dd[sp] = [seq]
                    d_sps[nn] = dd
                  
    # d_sps contains a dictionary of spacers for a genome
    #print(d_sps,"\n")
    
    d_res = {}
    d_gcas = {}
    d_gcrs = {}
    
    sep2 = ","
    
    for rep in d_gens[name] :
        fIn1 = open(path+"Input_files\\CRISPRCasdb\\"+rep+".csv","r")
        lines1 = fIn1.readlines()
        d_crs = {}
        d_cas = {}
        for ln in lines1:
            ln = ln.strip()
            if (ln != ""):
                t_ln = ln.split(sep2)
        
                r_id = t_ln[0]
                r_class = t_ln[1]
                r_name = t_ln[2]
                
        
                if r_class.lower() == "locuscrispr":
                    if r_id not in d_crs:
                        d_crs[r_id] = t_ln[2:]
                        
                if r_class.lower() == "clustercas":
                    if r_id not in d_cas:
                        d_cas[r_id] = t_ln[2:]
                            
                        
        fIn1.close()
    
        d_gcas[rep] = d_cas
        d_gcrs[rep] = d_crs

    
    
    for nn in d_sps:
        
        for rep in d_gens[name] :
            
            d_crs = d_gcrs[rep]
            d_cas = d_gcas[rep]
            
            fInName = path+"Input_files\\FASTA\\"+rep+".fasta"
            
            text = GetText(fInName)
            
            if len(text) != int(d_fname[rep][1]):
                print("Mismatch!!!")
            else:
                d_b = CreateDictLocD_upper(text,nn,d_fname[rep][2])

           
                for sp in d_sps[nn]:
                    sp_rc = Rev_cmp_upper(sp)
                    sp_len = len(sp)
                                  
                    # sp spacer
                    # nn length of a spaccer                  
                    
                    c_dir = 0
                    c_rev = 0
                    
                    c_in = 0
                    c_out = 0
                    
                    
                    
                    if sp in d_b:
                        c_dir = len(d_b[sp])
                        
                        t = sorted(d_b[sp].keys())
                        
                        
                        for pos in t:
                            fl = 0
                            
                            for it in d_crs:
                                ttt = d_crs[it]
                                if (pos >= int(ttt[1])) and (pos+sp_len <= int(ttt[1])+int(ttt[2])):
                                    fl = 1
                                    break
                
                            if fl == 1:
                                d_b[sp][pos] = "F"+"|"+ttt[0]+"_x"
                                c_in = c_in+1
                            else:
                                d_b[sp][pos] = "F"+"|"+"-"
                                c_out = c_out+1
                           
                    
                    if sp_rc in d_b:
                        c_rev = len(d_b[sp_rc])
                        
                        t = sorted(d_b[sp_rc].keys())


                    
                        for pos in t:
                            fl = 0
                            
                            for it in d_crs:
                                ttt = d_crs[it]
                                
                                if (pos >= int(ttt[1])) and (pos+sp_len <= int(ttt[1])+int(ttt[2])):
                                    fl = 1
                                    break
                
                            if fl == 1:
                                d_b[sp_rc][pos] = "R"+"|"+ttt[0]+"_x"
                                c_in = c_in+1
                            else:
                                d_b[sp_rc][pos] = "R"+"|"+"-"
                                c_out = c_out+1

                    
                    if sp in d_res:
                        if rep not in d_res[sp]:
                            d_res[sp][rep] = [c_dir,c_rev, c_in, c_out]
                        else:
                            print("Duplicate seq", rep)
                    else:
                        d_res[sp] = {}
                        if rep not in d_res[sp]:
                            d_res[sp][rep] = [c_dir,c_rev, c_in, c_out]
                        else:
                            print("Duplicate seq", rep)
                d_b.clear()
                
                
    
    num_sp = len(d_res)
    num_stsp = 0
    
    
    for sp in d_res:

        sp_total = 0
        sp_in = 0
        sp_out = 0
        sp_dir = 0
        sp_rev = 0
        
        sp_len = len(sp)
        
        for rep in d_res[sp]:
            total = d_res[sp][rep][0]+d_res[sp][rep][1]
            sp_dir = sp_dir + d_res[sp][rep][0]
            sp_rev = sp_rev + d_res[sp][rep][1]
            sp_total = sp_total+total
            sp_in = sp_in + d_res[sp][rep][2]
            sp_out = sp_out +d_res[sp][rep][3]

        
        if sp_out > 0:
            num_stsp = num_stsp + 1
        
        
        
        

        cas_ind = ""
        st_ind = ""
        plm_ind = ""
        
        if sp_len in d_sps:
            d_len_sps = d_sps[sp_len]
    

            
            
            t_cas = []
            fl_plm = 0
            for it in d_len_sps[sp]:
                
                d_cas = d_gcas[it]
                t_cas.append(d_cas)
                if "plasmid" in d_gens[name][it].lower():
                    fl_plm = 1
            
            if fl_plm == 0:
                plm_ind = "No"
            else:
                plm_ind = "Yes"
                
                
             
            flag = 0
            for jj in range(len(t_cas)):
                if len(t_cas[jj]) > 0:
                    flag = 1
                    
            if flag ==1:
                cas_ind = "Yes"
            else:
                cas_ind = "No"
        
        if sp_out > 0:
            st_ind = "Yes"
        else:
            st_ind = "No"
            
        
        fOutSP = open(fRname,"a")
        print(kind,name,sp,len(sp), sp_total, sp_in, sp_out, st_ind,plm_ind,cas_ind, file = fOutSP, sep = "\t")
        fOutSP.close()
        

    return [ct_chr, ct_plm,ct_other, num_sp, num_stsp]
    
    



#============================================
#-------------partA (summary file)

# Summary for replicons in genomes
fInName = path+"Input_files\\FASTA\\summary.txt"
res = GetDictFromFile(fInName,"\t",0)
d_fname = res[1]
print('Read-in',fInName.rsplit("\\",1)[1])


#---------------part B (wildcard spacers are discarded)

fInName = path+"Input_files\\CRISPRCasdb\\20190618_spacer_34.fasta"
fIn = open(fInName,"r")
lines = fIn.readlines()
fIn.close()

print("Read-in",fInName.rsplit("\\",1)[1])

d = {}

ct_wildcard = 0

fOut = open(path+"Output_files\\wildcard_spacers.txt","w")
fOut.write("seq_id(s)"+"\t"+"spacer"+"\t"+"spacer_length"+"\t"+"n_other"+"\n")

# d links spacers to the sequences where these spacers were found
# keys = spacers
# values = lists of sequence accession numbers
t_lines = [] 
for line in lines:
    line = line.strip()
    if line[0] == ">":
        if t_lines != []:
            headline = t_lines[0].strip()
            headline = headline[1:]
            t_headline = headline.split(" ",1)
            name = t_headline[0]
            
            t_name = name.strip().split("|")
            
            sp = t_lines[1]
            sp_other = CountGC(sp)[1]
            
            if len(t_lines) == 2:
                if (sp_other == 0):
                    if (sp not in d):
                        d[sp] = t_name
                    else:
                        print("double entry!", sp)
                else:
                    ct_wildcard = ct_wildcard + 1
                    fOut.write(t_lines[0]+"\t"+ t_lines[1]+"\t"+str(len(sp))+"\t"+str(sp_other)+"\n")
                    
        t_lines = [line]
        
            
    else:
        if line.strip() != "":
            t_lines.append(line)


if t_lines != []:
    headline = t_lines[0].strip()
    headline = headline[1:]
    t_headline = headline.split(" ",1)
    name = t_headline[0]
    
    t_name = name.strip().split("|")

    sp = t_lines[1]
    sp_other = CountGC(sp)[1]    
   
    if len(t_lines) == 2:
        if (sp_other == 0):
            if (sp not in d):
                d[sp] = t_name
            else:
                print("double entry!", sp)
        else:
            ct_wildcard = ct_wildcard + 1
            fOut.write(t_lines[0]+"\t"+ t_lines[1]+"\t"+str(len(sp))+"\t"+str(sp_other))

fOut.close()    
num = len(d)

print("\n",ct_wildcard, "wildcard spacers were discarded from futher analysis")
print("(See wildcard_spacers.txt for details)","\n")

#d_rep_sps stores spacers for each sequence accession number
# keys = sequence accession number
# values = spacers found in the sequence


d_rep_sps = {}        
for sp in d:
    t_st = d[sp]
    for it in t_st:    
        if it not in d_rep_sps:
            d_rep_sps[it] = [sp]
        else:
            tt = d_rep_sps[it]
            tt.append(sp)
            d_rep_sps[it] = tt
            
#---------------------end of part B

#-------------------Main part-------------------

# The analysis of self-targeting events in the organisms of interest listed in export.txt
fInName = path+"Input_files\\CRISPRCasdb\\export.txt"
print("Read-in",fInName.rsplit("\\",1)[1],"\n")

fIn = open(fInName,"r")
lines = fIn.readlines()
fIn.close()
header = lines[0]
header = header.strip()
t_header = header.split("\t")
#print(header)
lines = lines[1:]
n = len(lines)
#print(fInName)

print("Analysis is in progress...","\n")


n_lines = len(lines)


sep ="\t"

d = {}
d_ids = {}
d_gens = {}
d_reps = {}




# Organism level statistics
fResName = path+"Output_files\\organisms.txt"
fRes = open(fResName,"w")
print("Order","Super kindom","Strain name", "Release date", "SEQ count","CAS count","CRISPR count","Assembly description","Number of spacers","Number of self-targeting spacers",file = fRes, sep = "\t")
fRes.close()

# d_gens links organisms to its replicons:
# keys = organism names
# values = the corresponding replicons

# d_reps links replicons to the corresponding organism
# keys = replicon accession numbers
# values = the corresponding organism name


# Spacer level statistics
fRname = path+"Output_files\\spacers.txt"
fOutSP = open(fRname,"w")
print("Super kindom","Strain name","Spacer","Spacer_length","total_cp","in_cp(spacer duplicates)","out_cp(targets)", "Self-targeting","On_plasmid","Cas_genes", file = fOutSP, sep = "\t")
fOutSP.close()


print("Num | Organism |Time,sec.")

ct = 0
ct_dup = 0
c = 0

for line in lines:
    line = line.strip()
    if line != "":
        ct = ct+1
        
        
        t_line = line.split(sep)
    
        kind = t_line[0]
        name = t_line[1]
        info = t_line[6]
        
        
        
        
        t_info = info.split("|")
        
        if name != t_info[3]:
            c = c+1
            print(name,"\n", t_info[3])
            
        if name not in d:
            d[name] = t_info[4:]
            tt = t_info[4:]
            nn = len(tt)
    
            if nn % 2 == 0:
                dd = {}
                for i in range(0,nn,2):
                    if tt[i] not in d_ids:
                        d_ids[tt[i]] = 1
                    else:
                        print("duplicate id!", tt[i])
                    if tt[i] not in dd:
                        dd[tt[i]] = tt[i+1]
                    if tt[i] not in d_reps:
                        d_reps[tt[i]] = name
                    else:
                        print("duplicate id!",tt[i])
                        
                
                t_st = tm.time()
                
                d_gens[name] = dd
        
                res = SearchOrganism_plm_cas(fRname, name,kind,ct, d_gens[name], d_rep_sps, d_fname)
            
                s_line = sep.join(line.split(sep)[0:(-2)])  
                
                fRes = open(fResName,"a")
                fRes.write(str(ct)+sep+s_line+sep+str(res[3])+sep+str(res[4])+"\n")
                fRes.close()
                
                print(ct,"|",name,"|", round(tm.time()-t_st,2),"sec.")

            
            
        else:
            print("Duplicate organism!",line)
            ct_dup = ct_dup+1
        

print("\n")
   
print("Number of organisms analyzed :", len(d_gens))
print("Number of replicons processed :",len(d_reps),"\n")

#-----------end of main part 

