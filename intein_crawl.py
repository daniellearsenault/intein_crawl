import os
cwd=os.getcwd()
#input = 'intein_crawl_TEST_input.fst'
input = raw_input('Input FASTA filename:\n')
Iout=''
Eout=''
#determine correct starting index to present first entry
start = 0
I_log_lines=[]
E_log_lines=[]

#getting intein and extein output filenames
if (input+'INTEIN_CRAWL_INT_LOG.txt') and (input+'INTEIN_CRAWL_EXT_LOG.txt') in os.listdir(cwd):
    I_log = open((input+'INTEIN_CRAWL_INT_LOG.txt'),'r')
    I_log_lines=I_log.readlines()
    int_decision = raw_input("Would you still like to name your INTEINS ONLY output file "+str(I_log_lines[0][:-1])+"?\n Type YES to keep this filename or NO to provide a new filename.\n")
    if (int_decision == 'YES')or(int_decision == 'Yes')or(int_decision == 'yes'):
        Iout = str(I_log_lines[0])
    if (int_decision == 'NO')or(int_decision == 'No')or(int_decision == 'no'):
        Iout = raw_input("Please provide a new INTEINS ONLY output file name:\n")
    E_log = open((input+'INTEIN_CRAWL_EXT_LOG.txt'),'r')
    E_log_lines=E_log.readlines()
    ext_decision = raw_input("Would you still like to name your EXTEINS ONLY output file "+str(E_log_lines[0][:-1])+"?\n Type YES to keep this filename or NO to provide a new filename.\n")
    if (ext_decision == 'YES')or(ext_decision == 'Yes')or(ext_decision == 'yes'):
        Eout = str(E_log_lines[0])
    if (ext_decision == 'NO')or(ext_decision == 'No')or(ext_decision == 'no'):
        Eout = raw_input("Please provide a new EXTEINS ONLY output file name:\n")
    I_log.close()
    E_log.close()

if (input+'INTEIN_CRAWL_INT_LOG.txt') not in os.listdir(cwd) and (input+'INTEIN_CRAWL_EXT_LOG.txt') not in os.listdir(cwd):
    Iout = raw_input("What would you like to name the INTEINS ONLY output file? ex. inteins_only.fst\n")
    I_log = open((input+'INTEIN_CRAWL_INT_LOG.txt'),'a')
    I_log.write(str(Iout))
    I_log.write('\n')
    I_log.close()
    Eout = raw_input("What would you like to name the EXTEINS ONLY output file? ex. exteins_only.fst\n")
    E_log = open((input+'INTEIN_CRAWL_EXT_LOG.txt'),'a')
    E_log.write(str(Eout))
    E_log.write('\n')
    E_log.close()

def get_last_log_entry(log):
    annos = []
    for line in log:
        if '>' in line:
            annos.append(log.index(line))
    last_i = annos[-1]
    last = log[last_i]
    anno_start_i = last.find('__')
    anno_start = last[1:anno_start_i]
    pickup = int(anno_start)
    return pickup+1

I_log = open((input+'INTEIN_CRAWL_INT_LOG.txt'),'r')
I_log_lines=I_log.readlines()
E_log = open((input+'INTEIN_CRAWL_EXT_LOG.txt'),'r')
E_log_lines=E_log.readlines()
if len(I_log_lines)>1 or len(E_log_lines)>1:
    #input has been analyzed previously, start with proper index
    #ex. if the last entry split was entry at index 7 in entries, start with index 8 in entries
    start = get_last_log_entry(I_log_lines)
I_log.close()
E_log.close()


#read input file
fst1 = open(input,'r')
fst = fst1.readlines()
#establish intein deliniating characters
Nterm = raw_input("What few residues (~2-4) denote the beginning of the N-terminal of your intein? ex. CL\n")
Cterm = raw_input("What residues (~2-4) denote the end of the C-terminal of your intein? ex. TGN\n")

counter = 0
entries=[]
begins = []

#parse entries from input file into array called entries, with index after every >
for line in fst:
    if '>' in line:
        begins.append(fst.index(line))
begins.append(len(fst))
for i in range(len(begins)-1):
    oldanno = fst[begins[i]]
    rest = oldanno[1:]
    newanno = ('>'+str(counter)+'____'+rest)
    hold = []
    hold.append(newanno)
    counter+=1
    for n in range(begins[i]+1,begins[i+1]):
        hold.append(fst[n])
    entries.append(hold)

#now index of entry in entries is easily retrieved from the annotation line


#-------------------------------functions-------------------------------------

def two_coords(N_line,N_index,C_line,C_index,entry):
    I_log = open((input+'INTEIN_CRAWL_INT_LOG.txt'),'a')
    E_log = open((input+'INTEIN_CRAWL_EXT_LOG.txt'),'a')

    C_index+=(len(Cterm)-1)
    #cut entry and append to intein and extein logs
    I_log.write(str(entry[0]))
    E_log.write(str(entry[0]))
    int_range = [x for x in range(N_line,C_line+1)]
    if N_line!=0:
        for n in range(1,N_line):
            E_log.write(str(entry[n]))
    E_log.write(str(entry[N_line][:N_index]))
    I_log.write(str(entry[N_line][N_index:]))
    for n in range(int_range[1],int_range[-1]):
        I_log.write(str(entry[n]))
    I_log.write(str(entry[C_line][:C_index+1]))
    E_log.write(str(entry[C_line][C_index+1:]))
    if C_line!=len(entry)-1:
        for n in range(int_range[-1]+1,len(entry)):
            E_log.write(str(entry[n]))

    I_log.write('\n')
    I_log.close()

    E_log.write('\n')
    E_log.close()

def get_hits(line,term,line_index):
    hits = []
    temp=line
    add_on=0
    while term in temp:
        i = temp.find(term)
        hit = str('line: '+str(line_index)+', index: '+str(i+add_on)+', '+str(temp[i:i+6]))
        hits.append(hit)
        add_on+=(temp.find(term)+1)
        temp=line[add_on+1:]
    return hits
#---------------------------------------------------------------------


i = start
while i < len(entries):
    entry = entries[i]
    #present the user with the entry
    print('*** ENTRY: ***')
    for l in entry:
        send = str('line'+str(entry.index(l))+':'+str(l[:-1]))
        print(send)
    #find N and C intein term hits
    N_hits = []
    C_hits =[]
    for line in entry:
        if entry.index(line)!=0:
            if len(get_hits(line,Nterm,entry.index(line)))!=0:
                N_hits.append(get_hits(line,Nterm,entry.index(line)))
            if len(get_hits(line,Cterm,entry.index(line)))!=0:
                C_hits.append(get_hits(line,Cterm,entry.index(line)))
    #report hits to user
    if len(N_hits)!=0:
        print('The intein N terminal marker was found in the following lines at the provided indices:')
        for hit in N_hits:
            print(hit)
    if len(C_hits)!=0:
        print('The intein C terminal marker was found in the following lines at the provided indices:')
        for hit in C_hits:
            print(hit)
    if len(N_hits)==0:print('No intein N terminal markers found for this entry.')
    if len(C_hits)==0:print('No intein C terminal markers found for this entry.')
    #obtain intein bounds from user
    print('Use the following prompts to choose where you would like to place the correct intein boundary markers:\n')
    print('\n')
    N_line = int(raw_input("What is the LINE of the GOOD N-term marker? ex. 6 -- Type -1 if None were retrieved to append FULL entry to both output files. Type -2 if you want to REMOVE this entry from your set (will remain in original input file, will not be added to output files.)\n"))
    N_index = int(raw_input("What is the INDEX within this line of the GOOD N-term marker? ex. 22 -- Type -1 if None were retrieved to append FULL entry to both output files. Type -2 if you want to REMOVE this entry from your set (will remain in original input file, will not be added to output files.)\n"))
    print('\n')
    C_line = int(raw_input("What is the LINE of the GOOD C-term marker? ex. 6 -- Type -1 if None were retrieved to append FULL entry to both output files. Type -2 if you want to REMOVE this entry from your set (will remain in original input file, will not be added to output files.)\n"))
    C_index = int(raw_input("What is the INDEX within this line of the GOOD C-term marker? ex. 22 -- Type -1 if None were retrieved to append FULL entry to both output files. Type -2 if you want to REMOVE this entry from your set (will remain in original input file, will not be added to output files.)\n"))
    #if bounds are unclear append entire entry to both logs
    if N_line==-1 or N_index==-1 or C_line==-1 or C_index==-1:
        print("Incomplete or no bounds provided. Adding full entry to both intein and extein files for manual editing by user.")
        I_log = open((input+'INTEIN_CRAWL_INT_LOG.txt'),'a')
        E_log = open((input+'INTEIN_CRAWL_EXT_LOG.txt'),'a')
        for line in entry:
            I_log.write(str(line))
            E_log.write(str(line))
        I_log.write('\n')
        E_log.write('\n')
        I_log.close()
        E_log.close()
    #if bounds are clear, append appropriate edited entries to both logs
    if N_line>-1 and N_index>-1 and C_line>-1 and C_index>-1 :
        two_coords(N_line,N_index,C_line,C_index,entry)

    i+=1

Int_out = open(Iout, 'w')
Ext_out = open(Eout, 'w')

I_log = open((input+'INTEIN_CRAWL_INT_LOG.txt'),'r')
E_log = open((input+'INTEIN_CRAWL_EXT_LOG.txt'),'r')
I_log_lines = I_log.readlines()
for line in I_log_lines:
    if I_log_lines.index(line)!=0:
        Int_out.write(line)
E_log_lines = E_log.readlines()
for line in E_log_lines:
    if E_log_lines.index(line)!=0:
        Ext_out.write(line)

I_log.close()
Int_out.close()
E_log.close()
Ext_out.close()

print('*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~*~*')
print('Inteins Only file has been written to '+str(Iout))
print('Exteins Only file has been written to '+str(Eout))
print('')
print('Thank you for using the intein_crawl script.')
print('Enjoy your files :)')
print('*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~*~*')

