#ifndef KMERSET_HEAD
#define KMERSET_HEAD

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "readset.h"
#include "bitarray.h"
#include <cstring>
#include "common.h"

using namespace std;

#pragma pack(2)
typedef struct KmerHashTable{
    unsigned long long int kmer;
    long int frequency;
    KmerHashTable(){
        kmer = 0;
        frequency = 0;
    }
}KmerHashTable;
#pragma pack ()

#pragma pack(2)
typedef struct KmerSetHashTable{
    unsigned long long int kmer;
    long int * frequency;
    KmerSetHashTable(){
        kmer = 0;
        frequency = NULL;
    }
}KmerSetHashTable;
#pragma pack ()

#pragma pack(2)
typedef struct KmerSetHashTableHead{
    KmerSetHashTable * kmerSetHashTable;
    unsigned long long int kmerSetHashTableCount;
    double * avgKmerFrequency;
    long int setNumber;
    long int kmerLength;
    KmerSetHashTableHead(){
        kmerSetHashTable = NULL;
        kmerSetHashTableCount = 0;
        avgKmerFrequency = NULL;
        setNumber = 0;
        kmerLength = 0;
    }
}KmerSetHashTableHead;
#pragma pack ()

#pragma pack(2)
typedef struct KmerSet{
    KmerHashTable * kmerHashTable;
    unsigned long long int kmerNumber;
    unsigned int kmerLength;
    unsigned long long int kmerHashTableCount;
    char * address;
    int addressNumber;
    KmerSet(){
        kmerHashTable = NULL;
        kmerNumber = 0;
        kmerLength = 0;
        kmerHashTableCount = 0;
        address = NULL;
        addressNumber = 0;
    }
}KmerSet;
#pragma pack ()

unsigned long long int KmerSetHashTableSearch(char * value, KmerSetHashTable * kmerSetHashTable, unsigned int kmerLength, unsigned long int kmerSetHashTableCount){
    
    unsigned long int i = 0;
    unsigned long int j = 0;
    int p = 0;
    i = Hash(value, kmerLength, kmerSetHashTableCount);
    unsigned long long int tempValue = 0;
    SetBit(&tempValue, 0, kmerLength, value);
    tempValue++;
    
    while(kmerSetHashTable[i].kmer!=0){   
        if(kmerSetHashTable[i].kmer == tempValue){
            return (i + 1);
        }
        i = (i+1)%kmerSetHashTableCount; 
    }
       
    return 0;
      
}

int CopyToDisc1(KmerSetHashTable * kmerSetHashTable, unsigned int kmerLength, unsigned long long int kmerSetHashTableCount, char * address, int freNumber){

    char * tempAddress = new char[1000];
    sprintf(tempAddress, "%d", 123);
    char * newAddress = new char[strlen(tempAddress)+strlen(address)+1];
    AppendRight(newAddress, address, tempAddress, 1);
          
    ofstream ocout;
    ocout.open(newAddress,ios::app);

    unsigned long int i = 0;
    char * temp11 = new char[kmerLength + 1];
    unsigned long long int temp = 0;
    
    for(i = 0; i < kmerSetHashTableCount; i++){
        if(kmerSetHashTable[i].kmer!=0){                 
            temp = kmerSetHashTable[i].kmer - 1;               
            GetBit(&temp,0,kmerLength, temp11);
            ocout<<temp11<<",";
            for(int j = 0; j<freNumber; j++){
                ocout<<kmerSetHashTable[i].frequency[j]<<",";
            }
            ocout<<endl;

        }
    }
    ocout.close();
}

int CopyToDisc(KmerHashTable * kmerHashTable, unsigned int kmerLength, unsigned long int kmerHashTableCount, char * address, int index){

    FILE * fp = NULL;
    if(index!=-1){
        char * tempAddress = new char[1000];
        sprintf(tempAddress, "%d", index);
        char * newAddress = new char[strlen(tempAddress)+strlen(address)+1];
        AppendRight(newAddress, address, tempAddress, 1);
        
        fp = fopen(newAddress, "w+");
        
        delete [] tempAddress;
        delete [] newAddress;
    }
    if(index==-1){
        fp = fopen(address, "a+");
    }
    
    unsigned long int i = 0;
    char * temp = new char[kmerLength + 1];
    int kmerLengthB = sizeof(unsigned long long int);
    int frequencyL= sizeof(long int);
    for(i = 0; i < kmerHashTableCount; i++){
        if(kmerHashTable[i].kmer>0){                 
            fwrite(&kmerHashTable[i].kmer, kmerLengthB, 1, fp);
            fwrite(&kmerHashTable[i].frequency, frequencyL, 1, fp);
            kmerHashTable[i].kmer = 0;
            kmerHashTable[i].frequency = 0;
        }
    }

    fclose(fp);
}

int CopyToDisc2(KmerSetHashTable * kmerSetHashTable, int setNumber, unsigned int kmerLength, unsigned long int kmerSetHashTableCount, char * address){

    FILE * fp = NULL;

    fp = fopen(address, "a+");
      
    unsigned long long int i = 0;
    unsigned long int j = 0;
    char * temp = new char[kmerLength + 1];
    int kmerLengthB = sizeof(unsigned long long int);
    int frequencyL= sizeof(long int);
    for(i = 0; i < kmerSetHashTableCount; i++){
        if(kmerSetHashTable[i].kmer>0){
            fwrite(&i, kmerLengthB, 1, fp);                 
            fwrite(&kmerSetHashTable[i].kmer, kmerLengthB, 1, fp);
            for(j =0;j<setNumber; j++){
                fwrite(&kmerSetHashTable[i].frequency[j], frequencyL, 1, fp);
            }
        }
    }

    fclose(fp);
}

void GetKmerSetHashTableFromAddress(KmerSetHashTable * kmerSetHashTable, int setNumber, char * address){
    int kmerLengthB = sizeof(unsigned long long int);
    int frequencyL= sizeof(long int);
    unsigned long long int i = 0;
    unsigned long int j = 0;
    FILE * fp = fopen(address, "rb+");
    
    unsigned long long int kmer = 0;
    long int frequency = 0;
    while(fread(&i, kmerLengthB, 1, fp) == 1){
        fread(&kmer, kmerLengthB, 1, fp);
        kmerSetHashTable[i].kmer = kmer;
        kmerSetHashTable[i].frequency = new long int[setNumber];
        for(j = 0; j<setNumber; j++){
            fread(&frequency, frequencyL, 1, fp);
            kmerSetHashTable[i].frequency[j] = frequency;
        }
    }
}

void KmerHashTableInsert(char * kmer, long int frequency, unsigned int kmerLength, KmerHashTable * kmerHashTable, unsigned long int kmerHashTableCount, unsigned long long int * kmerNumber){
     
    int j = 0;

    unsigned long int i = Hash(kmer, kmerLength, kmerHashTableCount);
    long int previousIndex = i;
    
    unsigned long long int tempKmer = 0;
    SetBit(&tempKmer, 0, kmerLength, kmer);
    tempKmer++;
    unsigned long long int cur = 0;
    do{
        previousIndex = i;
        i = (i + 1) % (kmerHashTableCount);
        cur = __sync_val_compare_and_swap(&kmerHashTable[previousIndex].kmer, 0, tempKmer);
        
    }while(cur!=0&&cur!=tempKmer);
    
    if(cur == 0){
        __sync_fetch_and_add(kmerNumber, 1);
    }
    
    __sync_fetch_and_add(&kmerHashTable[previousIndex].frequency, frequency);
      
}

void KmerSetHashTableInsert(char * kmer, long int frequency, int setIndex, unsigned int kmerLength, KmerSetHashTable * kmerHashTable, unsigned long int kmerHashTableCount, unsigned long long int * kmerNumber){
     
    int j = 0;

    unsigned long int i = Hash(kmer, kmerLength, kmerHashTableCount);
    long int previousIndex = i;
    
    unsigned long long int tempKmer = 0;
    SetBit(&tempKmer, 0, kmerLength, kmer);
    tempKmer++;
    unsigned long long int cur = 0;
    do{
        previousIndex = i;
        i = (i + 1) % (kmerHashTableCount);
        cur = __sync_val_compare_and_swap(&kmerHashTable[previousIndex].kmer, 0, tempKmer);
        
    }while(cur!=0&&cur!=tempKmer);
    
    if(cur == 0){
        __sync_fetch_and_add(kmerNumber, 1);
    }
    
    __sync_fetch_and_add(&kmerHashTable[previousIndex].frequency[setIndex], frequency);
          
}

#pragma pack(2)
typedef struct GetKmerSetP{
    ReadSet * readSet;
    KmerSet * kmerSet;
    int threadIndex;
    int totalThreadNumber;
    int * threadSignal;
    int * addressIndex;
}GetKmerSetP;
#pragma pack ()

int CompareKmerHashTable(const void * a, const void * b){
    KmerHashTable * ha = (KmerHashTable *)a;
    KmerHashTable * hb = (KmerHashTable *)b;
        
    if(ha->kmer==0&&hb->kmer>0)return -1;
    if(ha->kmer>0&&hb->kmer==0)return 1;
    if(ha->kmer==0&&hb->kmer==0)return 0;
    
    unsigned long long int temp = 7;
    return (((ha->kmer) - 1)&temp) - (((hb->kmer) - 1)&temp);

}

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

void * GetKmerSetThread(void * arg){
    GetKmerSetP * getKmerSetP = (GetKmerSetP *)arg;
    unsigned long int i = 0;
    unsigned long int p = 0;
    int j = 0;
	
	int kmerLength = getKmerSetP->kmerSet->kmerLength;
	int readLength = getKmerSetP->readSet->readLength;
	char * kmer = new char[kmerLength + 1];
	p = readLength * getKmerSetP->threadIndex;
	int offset = readLength * getKmerSetP->totalThreadNumber;
	i = getKmerSetP->threadIndex;
	
	//cout<<getKmerSetP->kmerSet->kmerHashTableCount<<endl;
	getKmerSetP->threadSignal[getKmerSetP->threadIndex] = 0;
	while(i<getKmerSetP->readSet->readNumber){                         
            for(j=0;j<readLength-kmerLength+1;j++){
                GetBit(getKmerSetP->readSet->readSet, p+j, kmerLength, kmer);
                KmerHashTableInsert(kmer, 1, kmerLength, getKmerSetP->kmerSet->kmerHashTable, getKmerSetP->kmerSet->kmerHashTableCount, &getKmerSetP->kmerSet->kmerNumber);
            }
            getKmerSetP->threadSignal[getKmerSetP->threadIndex] = 1;
            
            pthread_mutex_lock(&mutex);
            if(getKmerSetP->kmerSet->kmerNumber>(0.75*getKmerSetP->kmerSet->kmerHashTableCount)){
                while(1){
                    int t = 0;
                    for(t = 0; t < getKmerSetP->totalThreadNumber; t++){
                        if(getKmerSetP->threadSignal[t] == 0){
                            break;
                        }
                    }
                    if(t == getKmerSetP->totalThreadNumber){
                        break;
                    }                    
                }

                qsort(getKmerSetP->kmerSet->kmerHashTable, getKmerSetP->kmerSet->kmerHashTableCount, sizeof(KmerHashTable), CompareKmerHashTable);

                CopyToDisc(getKmerSetP->kmerSet->kmerHashTable, kmerLength, getKmerSetP->kmerSet->kmerHashTableCount, getKmerSetP->kmerSet->address, getKmerSetP->kmerSet->addressNumber);

                allkmernum = allkmernum + getKmerSetP->kmerSet->kmerNumber;
                getKmerSetP->kmerSet->kmerNumber = 0;
                getKmerSetP->kmerSet->addressNumber ++;
                           
            }else{
                getKmerSetP->threadSignal[getKmerSetP->threadIndex] = 0;
            }
            pthread_mutex_unlock(&mutex);
            p = p + offset; 
            i = i + getKmerSetP->totalThreadNumber;                                  
		
	}
	getKmerSetP->threadSignal[getKmerSetP->threadIndex] = 1;

	
	delete [] kmer;

}

void GetKmerSet(ReadSet * readSet, KmerSet * kmerSet, int totalThreadNumber){
    
    pthread_t tid[totalThreadNumber];
    
    int i = 0;
    GetKmerSetP * getKmerSetP = new GetKmerSetP[totalThreadNumber];
    int * threadSignal = new int[totalThreadNumber]; 
    allkmernum = 0;
    
    for(i = 0; i< totalThreadNumber; i++){
        getKmerSetP[i].readSet = readSet;
        getKmerSetP[i].kmerSet = kmerSet;
        getKmerSetP[i].threadIndex = i;
        getKmerSetP[i].totalThreadNumber = totalThreadNumber;
        getKmerSetP[i].threadSignal = threadSignal;
        threadSignal[i] = 0;
             
        if(pthread_create(&tid[i], NULL, GetKmerSetThread, (void *)&getKmerSetP[i])){
            cout<<"create thread wrong!"<<endl;
            exit(1);
        }      
    }
    
    for(i = 0; i < totalThreadNumber; i++){
        pthread_join(tid[i], NULL);
    }

    qsort(kmerSet->kmerHashTable, kmerSet->kmerHashTableCount, sizeof(KmerHashTable), CompareKmerHashTable);
    CopyToDisc(kmerSet->kmerHashTable, kmerSet->kmerLength, kmerSet->kmerHashTableCount, kmerSet->address, kmerSet->addressNumber);
    allkmernum = allkmernum + kmerSet->kmerNumber;

    kmerSet->kmerNumber = 0;
    kmerSet->addressNumber ++;
}
#pragma pack(2)
typedef struct MergeKmerSetP{
    KmerSet * kmerSet;
    int threadIndex;
    int totalThreadNumber;
    int * threadSignal;
    int * addressIndex;
    int * copyToDiscToken;
}MergeKmerSetP;
#pragma pack ()

void * MergeKmerSetThread(void * arg){
    MergeKmerSetP * mergeKmerSetP = (MergeKmerSetP *)arg;
    unsigned long int i = 0;
    unsigned long int p = 0;
    int j = 0;
    int q = 0;
	char * kmer = NULL;
	int kmerLength = mergeKmerSetP->kmerSet->kmerLength;
	
	FILE * fp[mergeKmerSetP->kmerSet->addressNumber];
	
	
	long int tempFrequency = 0;
	int kmerLengthB = sizeof(unsigned long long int);
	unsigned long long int tempKmer = 0;
	int frequencyL = sizeof(long int);
	
	for(q = 0; q < mergeKmerSetP->kmerSet->addressNumber; q++){
            char * tempAddress = new char[1000];
            sprintf(tempAddress, "%d", q);
            char * newAddress = new char[strlen(tempAddress)+strlen(mergeKmerSetP->kmerSet->address)+1];
            AppendRight(newAddress, mergeKmerSetP->kmerSet->address, tempAddress, 1);        
            fp[q] = fopen(newAddress, "a+");
            if(fp[q] == NULL){
                //cout<<"ddddddddddd"<<endl;
                exit(0);
            }
            fseek(fp[q], (kmerLengthB+frequencyL)*mergeKmerSetP->threadIndex, 1);   
            delete [] tempAddress;
            delete [] newAddress;     
    }
    
	
	
	int offset = (kmerLengthB+frequencyL)*(mergeKmerSetP->totalThreadNumber-1);
	char * tempKmer1 = new char[kmerLength + 1];
	
	for(j = 0; j < 8; j++){

        mergeKmerSetP->threadSignal[mergeKmerSetP->threadIndex] = 0;
        for(q = 0; q < mergeKmerSetP->kmerSet->addressNumber; q++){
            while(fread(&tempKmer, kmerLengthB, 1, fp[q])== 1){              
                tempKmer--;
                unsigned long long int temp = (tempKmer & 7);
                if((int)temp != j){

                    break;
                }              
                GetBit(&tempKmer, 0, kmerLength, tempKmer1);                
                fread(&tempFrequency, frequencyL, 1, fp[q]);
                KmerHashTableInsert(tempKmer1, tempFrequency, kmerLength, mergeKmerSetP->kmerSet->kmerHashTable, mergeKmerSetP->kmerSet->kmerHashTableCount, &mergeKmerSetP->kmerSet->kmerNumber);
                if(fseek(fp[q], offset, 1)!=0){
                    break;
                }                          
            }
            fseek(fp[q], -1*kmerLengthB, 1);         
        }
        mergeKmerSetP->threadSignal[mergeKmerSetP->threadIndex] = 1;
        
        pthread_mutex_lock(&mutex);
        if(*(mergeKmerSetP->copyToDiscToken)==j){
            (*(mergeKmerSetP->copyToDiscToken)) ++;                                     
            while(1){
                int t = 0;
                for(t = 0; t < mergeKmerSetP->totalThreadNumber; t++){
                    if(mergeKmerSetP->threadSignal[t] == 0){
                        break;
                    }
                }
                if(t == mergeKmerSetP->totalThreadNumber){
                    break;
                }                    
            }

            CopyToDisc(mergeKmerSetP->kmerSet->kmerHashTable, kmerLength, mergeKmerSetP->kmerSet->kmerHashTableCount, mergeKmerSetP->kmerSet->address, -1);
            allkmernum = allkmernum + mergeKmerSetP->kmerSet->kmerNumber;
            mergeKmerSetP->kmerSet->kmerNumber = 0;

        }
        pthread_mutex_unlock(&mutex);
        
    }
    
    for(q = 0; q < mergeKmerSetP->kmerSet->addressNumber; q++){           
        fclose(fp[q]);
    }

}

void MergeKmerSet(KmerSet * kmerSet, int totalThreadNumber){
    
    pthread_t tid[totalThreadNumber];
    
    allkmernum = 0;
    
    int i = 0;
    MergeKmerSetP * mergeKmerSetP = new MergeKmerSetP[totalThreadNumber];
    int * threadSignal = new int[totalThreadNumber];
    int * copyToDiscToken = new int; 
    * copyToDiscToken = 0;
    for(i = 0; i< totalThreadNumber; i++){
        mergeKmerSetP[i].kmerSet = kmerSet;
        mergeKmerSetP[i].threadIndex = i;
        mergeKmerSetP[i].totalThreadNumber = totalThreadNumber;
        mergeKmerSetP[i].threadSignal = threadSignal; 
        mergeKmerSetP[i].copyToDiscToken = copyToDiscToken; 
             
        if(pthread_create(&tid[i], NULL, MergeKmerSetThread, (void *)&mergeKmerSetP[i])){
            cout<<"create thread wrong!"<<endl;
            exit(1);
        }      
    }
    
    for(i = 0; i < totalThreadNumber; i++){
        pthread_join(tid[i], NULL);
    }
    

    
    for(i = 0; i < kmerSet->addressNumber; i++){
        char * tempAddress = new char[1000];
        sprintf(tempAddress, "%d", i);
        char * newAddress = new char[strlen(tempAddress)+strlen(kmerSet->address)+1];
        AppendRight(newAddress, kmerSet->address, tempAddress, 1);
        
        remove(newAddress);
        delete [] tempAddress;
        delete [] newAddress;        
    }
    
    
}

#pragma pack(2)
typedef struct GetKmerNumberP{
    char * address;
    int kmerLength;
    unsigned long int number;
    int threadIndex;
    int totalThreadNumber;
}GetKmerNumberP;
#pragma pack ()

void * GetKmerNumberThread(void * arg){
    GetKmerNumberP * getKmerNumberP = (struct GetKmerNumberP *)arg;  	
	FILE * fp = fopen(getKmerNumberP->address,"a+");
	
	long int tempFrequency = 0;
    int kmerLengthB = sizeof(unsigned long long int);
	int frequencyL = sizeof(long int);
	fseek(fp, (kmerLengthB + frequencyL)*getKmerNumberP->threadIndex + kmerLengthB, 1);
	int offset = (kmerLengthB + frequencyL)*getKmerNumberP->totalThreadNumber - frequencyL;
	while(1==fread(&tempFrequency, frequencyL, 1, fp)){   
		if(tempFrequency>minKmerFrequency){		    
            getKmerNumberP->number ++;
		}
		if(fseek(fp, offset, 1)!=0){
            break;
        }
	}
	fclose(fp);
	unsigned long int * ret = new unsigned long int;
	* ret = getKmerNumberP->number;
	return (void *)ret;
}

long int GetKmerNumber(char * address, int kmerLength, int totalThreadNumber){
    int i = 0;
    pthread_t tid[totalThreadNumber];
    GetKmerNumberP getKmerNumberP[totalThreadNumber];
    long int num = 0;
    
    for(i = 0; i < totalThreadNumber; i++){
        getKmerNumberP[i].address = address;
        getKmerNumberP[i].kmerLength = kmerLength; 
        getKmerNumberP[i].number = 0;
        getKmerNumberP[i].threadIndex = i;
        getKmerNumberP[i].totalThreadNumber = totalThreadNumber;
        if(pthread_create(&tid[i], NULL, GetKmerNumberThread, (void *)&getKmerNumberP[i])){
            cout<<"create thread wrong!"<<endl;
            exit(1);
        }
    }
    void * ret[totalThreadNumber];
    for(i = 0; i < totalThreadNumber; i++){
        pthread_join(tid[i], &ret[i]);
    }
    for(i = 0; i < totalThreadNumber; i++){
        num = num + (*(unsigned long int *)ret[i]);
    }
    //cout<<num<<endl;
    return num;
}

#pragma pack(2)
typedef struct GetKmerSetAddressP{
    char * address;
    KmerSetHashTable * kmerSetHashTable;
    unsigned long long int kmerSetHashTableCount;
    unsigned long long int kmerNumber;
    int kmerLength;
    int freNumber;
    int freIndex;
    int threadIndex;
    int totalThreadNumber;
}GetKmerSetAddressP;
#pragma pack ()

void * GetKmerSetAddressThread(void * arg){
    GetKmerSetAddressP * getKmerSetAddressP = (GetKmerSetAddressP *)arg;
    unsigned long int i = 0;
    unsigned long int p = 0;
    int j = 0;
	char * kmer = NULL;
	int kmerLength = getKmerSetAddressP->kmerLength;
	int freNumber = getKmerSetAddressP->freNumber;
	int freIndex = getKmerSetAddressP->freIndex;
	
	FILE *fp = fopen(getKmerSetAddressP->address, "rb+");
	int kmerLengthB = sizeof(unsigned long long int);
	int frequencyL = sizeof(long int);
	unsigned long long int tempKmer = 0;
	long int tempFrequency = 0;
	char * tempKmer1 = new char[kmerLength + 1];
	
	fseek(fp, (kmerLengthB + frequencyL)*getKmerSetAddressP->threadIndex + kmerLengthB , 1);
	int offset = (kmerLengthB + frequencyL)*(getKmerSetAddressP->totalThreadNumber - 1) + kmerLengthB;
	int offset1 = (kmerLengthB + frequencyL)*(getKmerSetAddressP->totalThreadNumber);
		
	while(fread(&tempFrequency, frequencyL, 1, fp)==1){
        if(tempFrequency>minKmerFrequency){
            fseek(fp, -1*(kmerLengthB + frequencyL), 1);
            fread(&tempKmer, kmerLengthB, 1, fp);
            tempKmer--;
            GetBit(&tempKmer, 0, kmerLength, tempKmer1);
            KmerSetHashTableInsert(tempKmer1, tempFrequency, freIndex, kmerLength, getKmerSetAddressP->kmerSetHashTable, getKmerSetAddressP->kmerSetHashTableCount, &getKmerSetAddressP->kmerNumber);
            if(fseek(fp, offset1, 1)!=0){
                break;
            }
            continue;
        }
        if(fseek(fp, offset, 1)!=0){
            break;
        }
	}
}

unsigned long int GetKmerSetAddress(char * address, KmerSetHashTable * kmerSetHashTable, unsigned long int kmerSetHashTableCount, int kmerLength, int setNumber, int setIndex, int totalThreadNumber){
    
    pthread_t tid[totalThreadNumber];
    
    int i = 0;
    GetKmerSetAddressP * getKmerSetAddressP = new GetKmerSetAddressP[totalThreadNumber];
    unsigned long int num = 0;
    
    for(i = 0; i< totalThreadNumber; i++){
        getKmerSetAddressP[i].address = address;
        getKmerSetAddressP[i].kmerSetHashTable = kmerSetHashTable;
        getKmerSetAddressP[i].kmerSetHashTableCount = kmerSetHashTableCount;
        getKmerSetAddressP[i].kmerNumber = 0;
        getKmerSetAddressP[i].kmerLength = kmerLength;
        getKmerSetAddressP[i].freNumber = setNumber;
        getKmerSetAddressP[i].freIndex = setIndex;
        getKmerSetAddressP[i].threadIndex = i;
        getKmerSetAddressP[i].totalThreadNumber = totalThreadNumber;
             
        if(pthread_create(&tid[i], NULL, GetKmerSetAddressThread, (void *)&getKmerSetAddressP[i])){
            cout<<"create thread wrong!"<<endl;
            exit(1);
        }      
    }
    
    for(i = 0; i < totalThreadNumber; i++){
        pthread_join(tid[i], NULL);
    }
    
    for(i = 0; i<totalThreadNumber; i++){
        num = num + getKmerSetAddressP[i].kmerNumber;
    }
    
    remove(address);

    return num;
}

double * GetAvgKmerFrequency(KmerSetHashTableHead * kmerSetHashTableHead){
    long int i = 0;
    long int j = 0;
    long int num = 0;
    
    long int all = 0;
    
    double * avg = new double[kmerSetHashTableHead->setNumber];
    
    for(i = 0; i<kmerSetHashTableHead->setNumber; i++){
        avg[i] = 0;
        all = 0;
        num = 0;
        for(j = 0; j<kmerSetHashTableHead->kmerSetHashTableCount; j++){
            if(kmerSetHashTableHead->kmerSetHashTable[j].frequency!=NULL&&kmerSetHashTableHead->kmerSetHashTable[j].frequency[i]!=0){
                all = all + kmerSetHashTableHead->kmerSetHashTable[j].frequency[i];
                num++;
            }
        }
        avg[i] = (double)all/(double)num;
    }
    
    return avg;
}

double GetAvgKmerFrequency(char * contig, long int readSetIndex, long int start, long int end, KmerSetHashTableHead * kmerSetHashTableHead){
    
    if(end<0){
        return 0;
    }
    
    long int i = 0;
    long int j = 0;
    long int kmerLength = kmerSetHashTableHead->kmerLength;
    char * tempKmer = new char[kmerLength + 1];
    char * tempKmerRC = new char[kmerLength + 1];
    
    
    double avgKmerNum = 0;
    long int matchKmerNumber = 0;
    
    for(j=start;j<end;j++){
        
        SubContig(tempKmer, contig, j, j+kmerLength);
        ReverseComplement(tempKmer, tempKmerRC);

        long int hashP = KmerSetHashTableSearch(tempKmer,kmerSetHashTableHead->kmerSetHashTable,kmerLength,kmerSetHashTableHead->kmerSetHashTableCount);
        long int hashP1 = KmerSetHashTableSearch(tempKmerRC,kmerSetHashTableHead->kmerSetHashTable,kmerLength,kmerSetHashTableHead->kmerSetHashTableCount);

        if(hashP!=0){
            avgKmerNum = avgKmerNum + 
                kmerSetHashTableHead->kmerSetHashTable[hashP-1].frequency[readSetIndex];
            matchKmerNumber++;
        } 
        if(hashP1!=0){
            avgKmerNum = avgKmerNum + 
                kmerSetHashTableHead->kmerSetHashTable[hashP1-1].frequency[readSetIndex];
            matchKmerNumber++;
        }       
    }
    avgKmerNum = avgKmerNum/matchKmerNumber;
    return avgKmerNum;
    
}

void InitKmerSet(ReadSet * readSet, KmerSet * kmerSet, KmerSetHashTableHead * kmerSetHashTableHead, int kmerLength, int setNumber, unsigned int totalThreadNumber){
    int i = 0;
    unsigned long long int j = 0;
    unsigned long int allReadNumber = 0;
    allkmernum = 0;
    for(i = 0; i<setNumber; i++){
        allReadNumber = allReadNumber + (readSet[i].readLength - kmerLength + 1)*readSet[i].readNumber;
    }
    KmerHashTable * kmerHashTable = new KmerHashTable[allReadNumber];

    for(i = 0; i<setNumber; i++){

        kmerSet[i].kmerLength = kmerLength;
        kmerSet[i].kmerHashTableCount = allReadNumber;
        kmerSet[i].kmerHashTable = kmerHashTable;
        GetKmerSet(readSet+i, kmerSet+i, totalThreadNumber);
        MergeKmerSet(kmerSet+i, totalThreadNumber);
    }
    

    
    delete [] kmerHashTable;
    unsigned long long int kmerSetHashTableCount = 0;
    for(i = 0; i<setNumber; i++){
        kmerSet[i].kmerNumber = GetKmerNumber(kmerSet[i].address, kmerLength, totalThreadNumber);
        kmerSetHashTableCount = kmerSetHashTableCount + kmerSet[i].kmerNumber;
        kmerSet[i].kmerHashTable = NULL;
    }
    
    kmerSetHashTableCount = (long int)(1.2*kmerSetHashTableCount);
    
    //ofstream ocout;
    //char * temp11 = new char[30];
    //strcpy(temp11,"kmersethashcount.fa");
    //ocout.open(temp11);


    KmerSetHashTable * kmerSetHashTable = new KmerSetHashTable[kmerSetHashTableCount];
    
    for(j = 0; j<kmerSetHashTableCount; j++){
        kmerSetHashTable[j].frequency = new long int[setNumber];
        for(i = 0; i<setNumber; i++){
            kmerSetHashTable[j].frequency[i] = 0;
        }  
    }
    
    for(i = 0; i<setNumber; i++){
        GetKmerSetAddress(kmerSet[i].address, kmerSetHashTable, kmerSetHashTableCount, kmerLength, setNumber, i, totalThreadNumber); 
    }
    
    
    char * temp = new char[30];
    strcpy(temp,"allKmerFrequency.fa");
    CopyToDisc2(kmerSetHashTable, setNumber, kmerLength, kmerSetHashTableCount, temp);
    
    
    KmerSetHashTable * kmerSetHashTable1 = new KmerSetHashTable[kmerSetHashTableCount];
    GetKmerSetHashTableFromAddress(kmerSetHashTable1, setNumber, temp);
    
    remove(temp); 
    //strcpy(temp,"allkmer.fa");
    //CopyToDisc1(kmerSetHashTable1, kmerLength, kmerSetHashTableCount, temp, setNumber);
    
    
    kmerSetHashTableHead->kmerSetHashTable = kmerSetHashTable;
    kmerSetHashTableHead->kmerSetHashTableCount = kmerSetHashTableCount;
    kmerSetHashTableHead->setNumber = setNumber;
    kmerSetHashTableHead->avgKmerFrequency = GetAvgKmerFrequency(kmerSetHashTableHead);
    
}



#endif
