# intein_crawl
Walk through MSFASTA of intein-containing entries with consistent identifying Nterm and Cterm boundary characters to yield intein-only and extein-only output files with original input file order conserved.

**version 0.0 April 22, 2022**
This preliminary version is only designed to handle entries with **single** insertions. I am working on a second version that can handle situations such as genes invaded by multiple inteins, stay tuned :)


**INPUT:** MSFASTA file (prot. ideally) of entries that contain a nested element, such as an intein. The element must have a set of ~2-4 identifying characters for both the N-terminal and C-terminal that appear in the majority of the instances of the element. (ex. an intein that nearly always starts with 'CL' and ends with 'TGN')

**OUTPUT:** 2 MSFASTA files: an edited internal/interrupting element-only set of input entries with original order retained (**ex. intein only**), and an edited external/host element-only set of input entries with original order retained (**ex. extein only**)

**HOW TO USE:** Save script and input file in one directory, and set as working directory in terminal/terminal-equivalent. Use the command python intein-crawl.py to begin. The user will be greeted with different prompts and walked through the process.

**Flow:**
1) Use command python intein_crawl.py to begin.  
<img width="1052" alt="image" src="https://user-images.githubusercontent.com/56440050/164769548-819e8dff-41a7-4d3a-8c3f-ba642f8c3e8d.png"> 

2) The user will be prompted for the name of their input file (ex. **test_ints_and_exts.fst**), what they would like to name their output files (**ex. test_inteins_only.fst and test_extein_only.fst**), and the characters that appear at the N-terminal and C-terminal ends of their intein (**ex. if the intein being searched for typically follows a 'CL[...]TGN' sequence structure, 'CL' would designate the N-terminal and 'TGN' would designate the C-terminal**).  
<img width="1051" alt="image" src="https://user-images.githubusercontent.com/56440050/164770371-81e42a43-311e-4cf6-bb0f-e03eb77038f2.png"> 

3) The user will then be walked through each entry of their input file, will be presented with the (line,index) pair of each N-term marker and C-term marker found in the entry, and will ask the user to identify the correct boundaries to splice the entry at. The full entry will be shown to the user with each round as well, so they may see where the provided hits fall in the full entry to better inform bound selection. If the bounds are unclear, simply submit -1 when prompted for the bounds, and the full original entry will be appended to both the intein and extein output files for the user to manually inspect afterwards.  

3a) Clear bounds:
![image](https://user-images.githubusercontent.com/56440050/164780984-59ec636c-bd48-4eb0-bcb1-66abe9cc3c5e.png)
3b) Unclear bounds:
![image](https://user-images.githubusercontent.com/56440050/164781204-b7882886-43aa-4f08-a7a6-59b71975257e.png)

4) Since this script will likely be used on larger data sets, the user may terminate the run at any point, and any entries the user has already sorted through will be stored in appropriate log files. The next time the script is run a directory with the input and log files present, the run will pick back up at the entry the user left off at after the previous run. At the beginning of a secondary run, the user will be asked if they would like to keep the names of the output files they originally submitted, or if they would like to change them.   
Once all entries have been sorted through, the intein-only and extein-only files will be generated and saved in the working directory.  

4a) Beginning of second run on same input:
![image](https://user-images.githubusercontent.com/56440050/164781840-2cf9679a-d068-4091-988b-b4e3f495dbe1.png)
4b) Second run starts where previous run left off:
![image](https://user-images.githubusercontent.com/56440050/164782138-cc75fc50-9665-4a8b-bc86-c4b1a1b46c4a.png)


**Provided testing materials:** 

**intein_crawl.py** - script 

**test_ints_and_exts.fst** - input file 

Save these two files in a working directory, and use python intein_crawl.py to run, prompts will guide from there. 

**NOTE:** The N-terminal marker for the test data is **CL** , and the C-terminal marker for the test data is **TGN.** Enter these characters when prompted for N-term and C-term characters.  

<img width="680" alt="image" src="https://user-images.githubusercontent.com/56440050/164783586-c22ef902-fe7a-4ab3-8dda-f9802ce0305a.png">




p.s.
Script will be most efficient on **protein-level** sequence data.  

p.p.s
This script can be used for any multi-entry file with '>' annotation deliniations above each entry (ex. FASTA format) where each entry contains a nested element with relatively consistent characters at the N and C terminals, though the verbeage is intein-centric as that was the original purpose of the script. If using for non intein purposes, consider the intein-extein analog features in your data set: intein = the nested element with identifiers; extein = the host element surrounding the nested element.  

p.p.p.s. if thats a thing:
version 0.0 (April 22, 2022) is designed for use on single-intein entries (1 intein/interrupting element per entry). Currently working on a new version that can handle **multiply** invaded entries, which may take me some time to work out :)
