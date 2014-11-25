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
    
    
    cout<<"0--Getting Read Library Start."<<endl;
    InitReadSet(readSet, setNumber, threadNumber);
    cout<<"0--Getting Read Library End."<<endl;
    
    
    KmerSetHashTableHead * kmerSetHashTableHead = new KmerSetHashTableHead;
    kmerSetHashTableHead->kmerLength = globalKmerLength;
    cout<<"1--Getting K-mer set Start."<<endl;
    InitKmerSet(readSet,kmerSet, kmerSetHashTableHead, globalKmerLength, setNumber, threadNumber);
    cout<<"1--Getting K-mer set End."<<endl;
    
    DBGraphHead * deBruijnGraphHead = new DBGraphHead;
    
    cout<<"2--Getting De Bruijn Graph Start."<<endl;
    InitDBGraph(readSet, kmerLength, kmerSetHashTableHead,deBruijnGraphHead,threadNumber);
    cout<<"2--Getting De Bruijn Graph End."<<endl;
    
    for(i = 0; i<setNumber; i++){
        GetAvgKmerNumberAndGapProblity(deBruijnGraphHead->deBruijnGraph, readSet, i, kmerSetHashTableHead->kmerSetHashTable, kmerSetHashTableHead->kmerSetHashTableCount, kmerLength);
        //cout<<(readSet+i)->insertSize<<"--"<<(readSet+i)->var<<"--"<<(readSet+i)->gapProblity<<"--"<<(readSet+i)->avgKmerNumber<<"--"<<(readSet+i)->gapAvg<<"--"<<(readSet+i)->avgReadNumber<<endl;
    }
    
    ContigSet * contigSetHead = new ContigSet;
    
    cout<<"3--Construct Contig Set Start."<<endl;
    ConstructContigSet(deBruijnGraphHead, readSet, contigSetHead, setNumber, kmerSetHashTableHead->kmerSetHashTable, kmerSetHashTableHead->kmerSetHashTableCount, kmerLength, threadNumber, extendCutOff);
    cout<<"3--Construct Contig Set End."<<endl;
         
    char str[30];
    //strcpy(str,"contigSet.fa");
    char str1[30];
    //strcpy(str1,"contigSetLong.fa");
    
    ContigSet * temp = contigSetHead->next;
    
    
    temp = contigSetHead->next;
    //WriteContigSet(temp, str);
    //WriteContigSetLong(temp, str1);
    cout<<"4--Merge Contigs Start."<<endl;
    contigSetHead->next = ContigMerge(contigSetHead->next, deBruijnGraphHead->deBruijnGraph, kmerLength, readSet, setNumber);
    contigSetHead->next = ContigMerge(contigSetHead->next, deBruijnGraphHead->deBruijnGraph, kmerLength, readSet, setNumber);
    cout<<"4--Merge Contigs End."<<endl;
    //exit(0);
    temp = contigSetHead->next;
    i = 0;
    while(temp!=NULL){
        i++;
        temp = temp->next; 
    }
    
    
    //cout<<"Output Contig Set."<<endl;
    strcpy(str,"contigSet.fa");
    strcpy(str1,"contigSetLong.fa");
    
    temp = contigSetHead->next;
    WriteContigSet(temp, str);
    WriteContigSetLong(temp, str1);
    
    cout<<"5--Scaffold and Fill gap Start."<<endl;
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
    cout<<"5--Scaffold and Fill gap End."<<endl;    
    
    
    
}
