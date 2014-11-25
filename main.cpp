#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#include "pretest.h"

#include <pthread.h>
#include <time.h>
#include "bitarray.h"
#include "readset.h"
#include "kmerset.h"
#include "graph.h"
#include "contigmerge.h"
#include "constructcontigset.h"
#include "scaffolding.h"
#include "fillgap.h"


using namespace std;
 
int main(int argc, char *argv[])
{  
 
    int i = 0;
    int j = 0;
    if((argc-3)%3 != 0){
        cout<<"Please Input Correct Parameters!"<<endl;
        exit(0);
    }
    
    long int setNumber = (argc-3)/3;
    long int kmerLength = atoi(argv[argc-2]);
    long int threadNumber = atoi(argv[argc-1]);
    
    if(kmerLength<32){
        globalKmerLength = kmerLength;
    }else{
        globalKmerLength = 32;
    }
    
    //cout<<kmerLength<<"--"<<globalKmerLength<<endl;
    
    minKmerLength = kmerLength - 4;
    maxKmerLength = kmerLength + 4;
    
    ReadSet * readSet = new ReadSet[setNumber];
    KmerSet * kmerSet = new KmerSet[setNumber];
    
    for(i = 0; i<setNumber;i++){
        readSet[i].address = new char[100];
        kmerSet[i].address = new char[100];
        strncpy(readSet[i].address, argv[i*3+1],100);
        char * address = new char[100];
        strcpy(address, "kmerFrequency");
        char * tempAddress = new char[100];
        sprintf(tempAddress, "%d", i+1);
        strcat(address, tempAddress);   
        strcpy(kmerSet[i].address, address);
        
        readSet[i].insertSize = atoi(argv[i*3+2]);
        readSet[i].var = atoi(argv[i*3+3]);
        
        delete [] tempAddress;
        delete [] address;       
    }
    
    
    
    InitReadSet(readSet, setNumber, threadNumber);
    
    KmerSetHashTableHead * kmerSetHashTableHead = new KmerSetHashTableHead;
    kmerSetHashTableHead->kmerLength = globalKmerLength;
    
    InitKmerSet(readSet,kmerSet, kmerSetHashTableHead, globalKmerLength, setNumber, threadNumber);
       
    DBGraphHead * deBruijnGraphHead = new DBGraphHead;
    
    InitDBGraph(readSet, kmerLength, kmerSetHashTableHead,deBruijnGraphHead,threadNumber);
    
    for(i = 0; i<setNumber; i++){
        GetAvgKmerNumberAndGapProblity(deBruijnGraphHead->deBruijnGraph, readSet, i, kmerSetHashTableHead->kmerSetHashTable, kmerSetHashTableHead->kmerSetHashTableCount, kmerLength);
        //cout<<(readSet+i)->insertSize<<"--"<<(readSet+i)->var<<"--"<<(readSet+i)->gapProblity<<"--"<<(readSet+i)->avgKmerNumber<<"--"<<(readSet+i)->gapAvg<<"--"<<(readSet+i)->avgReadNumber<<endl;
    }
    
    ContigSet * contigSetHead = new ContigSet;
    ConstructContigSet(deBruijnGraphHead, readSet, contigSetHead, setNumber, kmerSetHashTableHead->kmerSetHashTable, kmerSetHashTableHead->kmerSetHashTableCount, kmerLength, threadNumber, extendCutOff);
    
         
    char str[30];
    //strcpy(str,"contigSet.fa");
    char str1[30];
    //strcpy(str1,"contigSetLong.fa");
    
    ContigSet * temp = contigSetHead->next;
    
    
    temp = contigSetHead->next;
    //WriteContigSet(temp, str);
    //WriteContigSetLong(temp, str1);
    
    contigSetHead->next = ContigMerge(contigSetHead->next, deBruijnGraphHead->deBruijnGraph, kmerLength, readSet, setNumber);
    contigSetHead->next = ContigMerge(contigSetHead->next, deBruijnGraphHead->deBruijnGraph, kmerLength, readSet, setNumber);
    //exit(0);
    temp = contigSetHead->next;
    i = 0;
    while(temp!=NULL){
        i++;
        temp = temp->next; 
    }
    
    
    
    strcpy(str,"contigSet.fa");
    strcpy(str1,"contigSetLong.fa");
    
    temp = contigSetHead->next;
    WriteContigSet(temp, str);
    WriteContigSetLong(temp, str1);
    
    
    ScaffoldSetHead * scaffoldSetHead = ScaffoldingContigSet(temp, readSet, setNumber, 3*kmerLength, kmerLength, threadNumber);
    
    
    //char * address0 = new char[30];
    //strcpy(address0, "scaffold0.fa");
    //WriteScaffoldSet(scaffoldSetHead, address0);
    FillGap(scaffoldSetHead, readSet, setNumber, kmerSetHashTableHead, kmerLength, 1);
    
    
    contigSetHead = GetContigSetFromScaffoldSetHead(scaffoldSetHead);
    DeleteScaffoldSetHead(scaffoldSetHead);
    temp = contigSetHead->next;
    
    scaffoldSetHead = ScaffoldingContigSet(temp, readSet, setNumber, kmerLength, kmerLength, threadNumber);    
    FillGap(scaffoldSetHead, readSet, setNumber, kmerSetHashTableHead, kmerLength, threadNumber); 
        
    
    
    
}
