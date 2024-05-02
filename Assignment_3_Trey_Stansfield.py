#imports the nesscary packages
import numpy as np
import pandas as pd

#This function is for ptt files
#It takes a dataframe from a ptt file and an output file name as input
def operon_maker_ptt(ptt,operon_file):
    #To offer more control over my iterations I am using a while loop
    #I need to set the counter varible beforehand
    x = 0
    #Since I am comparing the current index to its neighbor 
    #I need to the loop early to avoid a key error
    while x < len(ptt)-1:
        #flag to say if we are currently in operon
        in_operon = False
        #These varibles collect the strand of the current index and its next nieghbor
        #They also collect the trailing distance of the current
        #and the starting of the next nieghbor
        strand = ptt.loc[x,'Strand']
        next_strand = ptt.loc[x+1,'Strand']
        loc = int(ptt.loc[x,'Location'].split('..')[1])
        next_loc = int(ptt.loc[x+1,'Location'].split('..')[0])
        
        #This checks if the index and its next neighor have the same strand
        #and are within 50 bp of each other
        if strand == next_strand and next_loc - loc <= 50:
            #if true a list is created of the index gene and its neighbor
            #I used the Synonym as the name as Halobacterium had many null gene names
            operon = [ptt.loc[x,'Synonym']]
            operon.append(ptt.loc[x+1,'Synonym'])
            #this sets the  in_operon flag to true
            in_operon = True
        
        #This iterates to the next index
        x += 1
        
        #If we have two genes in operon this loop checks to see if the subsequent
        #genes are also in the same operon
        while in_operon:
            
            #same logic as above
            strand = ptt.loc[x,'Strand']
            next_strand = ptt.loc[x+1,'Strand']
            loc = int(ptt.loc[x,'Location'].split('..')[1])
            next_loc = int(ptt.loc[x+1,'Location'].split('..')[0])
            
            if strand == next_strand and next_loc - loc <= 50:
                #if the genes continue to be in operon we keep adding to this 
                #operon list and iterate
                #We only add the next gene into operon since we added the index
                #and its neighbor were added previously
                operon.append(ptt.loc[x+1,'Synonym'])
                x += 1
            
            else:
                #once the genes are no longer in one operon we switch the
                #flag to false
                in_operon = False
                #the operon is then written to a new line on the output
                #file
                with open(operon_file,'a') as write:
                    write.write(','.join(operon) + '\n')
                x += 1

#this function works the same as the above but it is for gff files                
def operon_maker_gff(gff,operon_file):
    x = 0
    
    while x < len(gff)-1:         
        in_operon = False
        #the main change is here 
        #All of the data comes from different columns
        strand = gff.loc[x,6]
        next_strand = gff.loc[x+1,6]
        loc = int(gff.loc[x,4])
        next_loc = int(gff.loc[x+1,3])
        
        if strand == next_strand and next_loc - loc <= 50:
            #Gff file is based on contigs so contig names used over 
            #gene names
            operon = [gff.loc[x,0]]
            operon.append(gff.loc[x+1,0])
            in_operon = True
            
        x += 1
        
        while in_operon:
           
            strand = gff.loc[x,6]
            next_strand = gff.loc[x+1,6]
            loc = int(gff.loc[x,4])
            next_loc = int(gff.loc[x+1,3])
            
            if strand == next_strand and next_loc - loc <= 50:
                operon.append(gff.loc[x+1,0])
                x += 1
                
            else:
                in_operon = False
                print(operon)
                with open(operon_file,'a') as write:
                    write.write(','.join(operon) + '\n')
                x += 1
                

#running for e.coli ptt
e_coli = pd.read_csv('E_coli_K12_MG1655.ptt',sep='\t', skiprows=2)

#Creates the output file and writes a header to it
with open('operon_list_ecoli.txt','a') as write:
    write.write('E.Coli Operon List:\n')

operon_maker_ptt(e_coli,'operon_list_ecoli.txt')

#running for b.subtilis ptt
b_subtilis = pd.read_csv('B_subtilis_168.ptt',sep='\t', skiprows=2)

with open('operon_list_bsubtilis.txt','a') as write:
    write.write('B.subtilis Operon List:\n')
    
operon_maker_ptt(b_subtilis,'operon_list_bsubtilis.txt')

#running for halobacterium ptt
halobacterium = pd.read_csv('Halobacterium_NRC1.ptt',sep='\t', skiprows=2)

with open('operon_list_halobacterium.txt','a') as write:
    write.write('Halobacterium Operon List:\n')
    
operon_maker_ptt(halobacterium,'operon_list_halobacterium.txt')

#running for synechocystis ptt
synechocystis = pd.read_csv('Synechocystis_PCC6803_uid159873.ptt',sep='\t', 
                            skiprows=2)
with open('operon_list_synechocystis.txt','a') as write:
    write.write('Synechocystis Operon List:\n')
    
operon_maker_ptt(synechocystis,'operon_list_synechocystis.txt')

#running for hoatzin gff
hoatzin = pd.read_csv('2088090036.gff',sep='\t', skiprows=1, header=None)

with open('operon_list_hoatzin.txt','a') as write:
    write.write('Crop microbiome from Hoatzin Operon List:\n')
    
operon_maker_gff(hoatzin, 'operon_list_hoatzin.txt')

