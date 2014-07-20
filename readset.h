#ifndef READSET_H_INCLUDED 
#define READSET_H_INCLUDED 
#include "fstream"
#include <cstring>
#include <iostream>
#include "stdio.h"
#include "stdlib.h"
#include <pthread.h>

#include "bitarray.h"
#include "common.h"

using namespace std;

#pragma pack(2)
typedef struct ReadHashTable{
    unsigned long int keyValue;
    unsigned long int index;
    ReadHashTable(){
        keyValue = 0;
        index = 0;
    }
}ReadHashTable;
#pragma pack ()

#pragma pack(2)
typedef struct ReadSet{
    char * readSet;
    unsigned int readLength;
    unsigned long int readNumber;
    unsigned int insertSize;
    double var;
    double gapAvg;
    double gapSD;
    double gapProblity;
    double avgKmerNumber;
    double avgReadNumber;
    char * address;
    ReadHashTable * readHashTable;
}ReadSet; 
#pragma pack ()

#pragma pack(2)
typedef struct ReadMate{
    char * readMate;
    struct ReadMate * next;
    ReadMate(){
        readMate = NULL;
        next = NULL;
    }
}ReadMate;
#pragma pack()

#pragma pack(2)
typedef struct LongRead{
    char * read;
    long int readLength;
    LongRead(){
        read = NULL;
        readLength = 0;
    }
}LongRead;
#pragma pack ()

#pragma pack(2)
typedef struct LongReadSet{
    LongRead * longReadSet;
    long int readNumber;
    char * address;
    LongReadSet(){
        longReadSet = NULL;
        readNumber = 0;
        address = NULL;
    }
}LongReadSet;
#pragma pack ()

long int SearchReadNumber(char * read, ReadSet * readSet){
    
    long int readNumber = 0;
    
    unsigned long int i = 0;
    long int j = 0;
    unsigned long int t = 0;
    int index = 0;
    //cout<<"tt00:"<<endl;
    i = Hash(read, readSet->readLength, readSet->readNumber);
    //cout<<"tt01:"<<endl;
    long int start = 0;
    long int end = readSet->readNumber-1;
    long int mid = 0;
    char temp[readSet->readLength + 1];
    temp[readSet->readLength] = '\0';
    //cout<<"tt:"<<i<<endl;
    while(start <= end){       
        mid = (start + end)/2;
        if(i == readSet->readHashTable[mid].keyValue){
            break;
        }else if(i > readSet->readHashTable[mid].keyValue){
            start = mid + 1;
        }else if(i < readSet->readHashTable[mid].keyValue){
            end = mid - 1;
        }
    }
    //cout<<"tt1:"<<start<<"--"<<end<<endl;
    if(start>end){
        return readNumber;
    }

    for(j = mid; j < readSet->readNumber && i == readSet->readHashTable[j].keyValue; j++){
        
        start = readSet->readLength * (readSet->readHashTable[j].index);
        //cout<<"tt2:"<<start<<endl;
        index = 0;
        for(index = 0; index<readSet->readLength; index++){
            if(read[index] != GetBit(readSet->readSet, start+index)){
                index = -1;
                break;
            }
        }
        
        if(index != -1){
            readNumber++;
        }
         
    }
    
    for(j = mid-1; j >= 0 && i == readSet->readHashTable[j].keyValue; j--){
        start = readSet->readLength * (readSet->readHashTable[j].index);
        index = 0;
        //cout<<"tt3:"<<start<<endl;
        for(index = 0; index<readSet->readLength; index++){
            if(read[index] != GetBit(readSet->readSet, start+index)){
                index = -1;
                break;
            }
        }
        
        if(index != -1){
            readNumber++;
        }
    }
    //cout<<"tt4:"<<endl;
    return readNumber;
    
}

long int SearchRead(char * read, ReadSet * readSet){
    unsigned long int i = 0;
    long int j = 0;
    unsigned long int t = 0;
    int index = 0;
    //cout<<"tt00:"<<endl;
    i = Hash(read, readSet->readLength, readSet->readNumber);
    //cout<<"tt01:"<<endl;
    long int start = 0;
    long int end = readSet->readNumber-1;
    long int mid = 0;
    char temp[readSet->readLength + 1];
    temp[readSet->readLength] = '\0';
    //cout<<"tt:"<<i<<endl;
    while(start <= end){       
        mid = (start + end)/2;
        if(i == readSet->readHashTable[mid].keyValue){
            break;
        }else if(i > readSet->readHashTable[mid].keyValue){
            start = mid + 1;
        }else if(i < readSet->readHashTable[mid].keyValue){
            end = mid - 1;
        }
    }
    //cout<<"tt1:"<<start<<"--"<<end<<endl;
    if(start>end){
        return -1;
    }

    for(j = mid; j < readSet->readNumber && i == readSet->readHashTable[j].keyValue; j++){
        
        start = readSet->readLength * (readSet->readHashTable[j].index);
        //cout<<"tt2:"<<start<<endl;
        index = 0;
        for(index = 0; index<readSet->readLength; index++){
            if(read[index] != GetBit(readSet->readSet, start+index)){
                index = -1;
                break;
            }
        }
        
        if(index != -1){
            return j;
        }
         
    }
    
    for(j = mid-1; j >= 0 && i == readSet->readHashTable[j].keyValue; j--){
        start = readSet->readLength * (readSet->readHashTable[j].index);
        index = 0;
        //cout<<"tt3:"<<start<<endl;
        for(index = 0; index<readSet->readLength; index++){
            if(read[index] != GetBit(readSet->readSet, start+index)){
                index = -1;
                break;
            }
        }
        
        if(index != -1){
            return j;
        }
    }
    //cout<<"tt4:"<<endl;
    return -1;
    
}

ReadMate * SearchRightReadMate(char * read, ReadSet * readSet){
    unsigned long int i = 0;
    unsigned long int j = 0;
    long int t = 0;
    
    unsigned long int start = 0;
    i = Hash(read, readSet->readLength, readSet->readNumber);
    long int mid = SearchRead(read, readSet);
    
    if(mid == -1){
        return NULL;
    }

    ReadMate * readMateHead = new ReadMate; 
    ReadMate * mate = readMateHead;
    
    char temp[readSet->readLength+1];
    
    for(j = mid; j < readSet->readNumber && i == readSet->readHashTable[j].keyValue; j++){
        
        start = readSet->readLength * (readSet->readHashTable[j].index);  
        for(t=0; t<readSet->readLength; t++){
            if(read[t]!=GetBit(readSet->readSet,start+t)){
                break;
            }
        }
        if(t!=readSet->readLength){
            continue;
        }
        
        ReadMate * tempMate = new ReadMate;
        tempMate->readMate = new char[readSet->readLength + 1];
        
        if(readSet->readHashTable[j].index%2 == 0){
            start = start + 2*readSet->readLength;
        }
        
        for(t=0;t<readSet->readLength;t++){
            if(GetBit(readSet->readSet, start -1 - t) == 'A'){
                tempMate->readMate[t] = 'T';
            }else if(GetBit(readSet->readSet, start -1 - t) == 'T'){
                tempMate->readMate[t] = 'A';
            }else if(GetBit(readSet->readSet, start -1 - t) == 'G'){
                tempMate->readMate[t] = 'C';
            }else if(GetBit(readSet->readSet, start -1 - t) == 'C'){
                tempMate->readMate[t] = 'G';
            } 
        }
        tempMate->readMate[t] = '\0';  
        
        tempMate->next = readMateHead->next;
        readMateHead->next = tempMate;
    }

    for(j = mid-1; j >= 0 && i == readSet->readHashTable[j].keyValue; j--){
        start = readSet->readLength * (readSet->readHashTable[j].index);  
        for(t=0; t<readSet->readLength; t++){
            if(read[t]!=GetBit(readSet->readSet,start+t)){
                break;
            }
        }
        if(t!=readSet->readLength){
            continue;
        }
        
        ReadMate * tempMate = new ReadMate;
        tempMate->readMate = new char[readSet->readLength + 1];
        
        if(readSet->readHashTable[j].index%2 == 0){
            start = start + 2*readSet->readLength;
        }
        
        for(t=0;t<readSet->readLength;t++){
            if(GetBit(readSet->readSet, start -1 - t) == 'A'){
                tempMate->readMate[t] = 'T';
            }else if(GetBit(readSet->readSet, start -1 - t) == 'T'){
                tempMate->readMate[t] = 'A';
            }else if(GetBit(readSet->readSet, start -1 - t) == 'G'){
                tempMate->readMate[t] = 'C';
            }else if(GetBit(readSet->readSet, start -1 - t) == 'C'){
                tempMate->readMate[t] = 'G';
            } 
        }
        tempMate->readMate[t] = '\0';  
        
        tempMate->next = readMateHead->next;
        readMateHead->next = tempMate;
    }

    return mate;
    
}

ReadMate * SearchLeftReadMate(char * read, ReadSet * readSet){
    unsigned long int i = 0;
    unsigned long int j = 0;
    unsigned long int t = 0;
    
    unsigned long int start = 0;
    char readRC[readSet->readLength+1];
    ReverseComplement(read, readRC);
    i = Hash(readRC, readSet->readLength, readSet->readNumber);
    
    long int mid = SearchRead(readRC, readSet);
    if(mid == -1){
        return NULL;
    }

    ReadMate * readMateHead = new ReadMate; 
    ReadMate * mate = readMateHead;
    
    char temp[readSet->readLength+1];
    
    for(j = mid; j < readSet->readNumber && i == readSet->readHashTable[j].keyValue; j++){
        start = readSet->readLength * (readSet->readHashTable[j].index);  
        for(t=0; t<readSet->readLength; t++){
            if(readRC[t]!=GetBit(readSet->readSet,start+t)){
                break;
            }
        }
        if(t!=readSet->readLength){
            continue;
        }
        
        ReadMate * tempMate = new ReadMate;
        tempMate->readMate = new char[readSet->readLength + 1];
        
        if(readSet->readHashTable[j].index%2 != 0){
            start = start - readSet->readLength;
        }else{
            start = start + readSet->readLength;
        }
        
        for(t=0;t<readSet->readLength;t++){    
            tempMate->readMate[t] = GetBit(readSet->readSet, start + t);   
        }
        tempMate->readMate[t] = '\0';  
        
        tempMate->next = readMateHead->next;
        readMateHead->next = tempMate;
    }

    for(j = mid-1; j >= 0 && i == readSet->readHashTable[j].keyValue; j--){
        start = readSet->readLength * (readSet->readHashTable[j].index);  
        for(t=0; t<readSet->readLength; t++){
            if(readRC[t]!=GetBit(readSet->readSet,start+t)){
                break;
            }
        }
        if(t!=readSet->readLength){
            continue;
        }
        
        ReadMate * tempMate = new ReadMate;
        tempMate->readMate = new char[readSet->readLength + 1];
        
        if(readSet->readHashTable[j].index%2 != 0){
            start = start - readSet->readLength;
        }else{
            start = start + readSet->readLength;
        }
        
        for(t=0;t<readSet->readLength;t++){    
            tempMate->readMate[t] = GetBit(readSet->readSet, start + t);   
        }
        tempMate->readMate[t] = '\0';  
        
        tempMate->next = readMateHead->next;
        readMateHead->next = tempMate;
    }

    return mate;
    
}

int ReadHashTableInsert(ReadSet * readSet, unsigned long int pos, char * key)
{
    unsigned long int i = Hash(key, readSet->readLength, readSet->readNumber);
    readSet->readHashTable[pos].keyValue = i;
    readSet->readHashTable[pos].index = pos;
    return 1;
}

#pragma pack(2)
typedef struct GetReadNumberP{
    char * address;
    long int number;
    int threadIndex;
    int totalThreadNumber;
}GetReadNumberP;
#pragma pack ()

void * GetReadNumberThread(void * arg){
    GetReadNumberP * getReadNumberP = (struct GetReadNumberP *)arg;
    unsigned long int i = 0;
	char temp[1000];   	
	FILE * fp = fopen(getReadNumberP->address,"a+");
	
	while(NULL!=fgets(temp,1000,fp)){
		if(i%2!=0&&((unsigned long int)(i/2)%getReadNumberP->totalThreadNumber == getReadNumberP->threadIndex)){		    
		    getReadNumberP->number ++;
		}
		i++;
	}
	fclose(fp);
	
	unsigned long int * ret = new unsigned long int;
	* ret = getReadNumberP->number;
	return (void *)ret;
}

long int GetReadNumber(char * address, int totalThreadNumber){
    int i = 0;
    pthread_t tid[totalThreadNumber];
    GetReadNumberP getReadNumberP[totalThreadNumber];
    long int num = 0;
    
    for(i = 0; i < totalThreadNumber; i++){
        getReadNumberP[i].address = address;
        getReadNumberP[i].number = 0;
        getReadNumberP[i].threadIndex = i;
        getReadNumberP[i].totalThreadNumber = totalThreadNumber;
        if(pthread_create(&tid[i], NULL, GetReadNumberThread, (void *)&getReadNumberP[i])){
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
typedef struct GetReadSetP{
    ReadSet * readSet;
    int threadIndex;
    int totalThreadNumber;
}GetReadSetP;
#pragma pack ()

#pragma pack(2)
typedef struct GetLongReadSetP{
    LongReadSet * longReadSet;
    int threadIndex;
    int totalThreadNumber;
}GetLongReadSetP;
#pragma pack ()

pthread_mutex_t mutex3 = PTHREAD_MUTEX_INITIALIZER;

void * GetLongReadSetThread(void * arg){
    GetLongReadSetP * getLongReadSetP = (struct GetLongReadSetP *)arg;
    long int i = 0;
	char * temp = new char[20000];
	FILE * fp = fopen(getLongReadSetP->longReadSet->address,"a+"); 

	while(NULL!=fgets(temp,20000,fp)){
        
		if(i%2 == 1 && (i/2)%getLongReadSetP->totalThreadNumber == getLongReadSetP->threadIndex){
            
            getLongReadSetP->longReadSet->longReadSet[i/2].readLength = strlen(temp);
            getLongReadSetP->longReadSet->longReadSet[i/2].read = new char[(getLongReadSetP->longReadSet->longReadSet[i/2].readLength/4) + 1];

            SetBit(getLongReadSetP->longReadSet->longReadSet[i/2].read, 0, getLongReadSetP->longReadSet->longReadSet[i/2].readLength, temp);
                
		}	
		i++;	
	}
	
	fclose(fp);
	return NULL;
	
}

unsigned long int GetLongReadSet(LongReadSet * longReadSet, unsigned int totalThreadNumber){

    pthread_t tid[totalThreadNumber];
    
    int i = 0;
    GetLongReadSetP * getLongReadSetP = new GetLongReadSetP[totalThreadNumber];
    longReadSet->longReadSet = new LongRead[longReadSet->readNumber];
    
    for(i = 0; i< totalThreadNumber; i++){
        getLongReadSetP[i].longReadSet = longReadSet;
        getLongReadSetP[i].threadIndex = i;
        getLongReadSetP[i].totalThreadNumber = totalThreadNumber;
             
        if(pthread_create(&tid[i], NULL, GetLongReadSetThread, (void *)&getLongReadSetP[i])){
            cout<<"create thread wrong!"<<endl;
            exit(1);
        }      
    }
    
    for(i = 0; i < totalThreadNumber; i++){
        pthread_join(tid[i], NULL);
    }
    
}

void * GetReadSetThread(void * arg){
    GetReadSetP * getReadSetP = (struct GetReadSetP *)arg;
    long int i = 0;
	int readLength = getReadSetP->readSet->readLength;
	char temp[readLength + 2];
	FILE * fp = fopen(getReadSetP->readSet->address,"a+"); 

    unsigned long int p = readLength * (getReadSetP->threadIndex);
	while(NULL!=fgets(temp,1000,fp)){
		if(i%2 == 1 && (i/2)%getReadSetP->totalThreadNumber == getReadSetP->threadIndex){
            
            pthread_mutex_lock(&mutex3);
            SetBit(getReadSetP->readSet->readSet, p, readLength, temp);
            pthread_mutex_unlock(&mutex3);
            
            ReadHashTableInsert(getReadSetP->readSet, i/2, temp);

            p = p + readLength * getReadSetP->totalThreadNumber;
                
		}	
		i++;	
	}
	
	fclose(fp);
	return NULL;
	
}

unsigned long int GetReadSet(ReadSet * readSet, unsigned int totalThreadNumber){

    pthread_t tid[totalThreadNumber];
    
    int i = 0;
    GetReadSetP * getReadSetP = new GetReadSetP[totalThreadNumber];
    readSet->readSet = new char[((readSet->readLength * readSet->readNumber)/4) + 1];
    readSet->readHashTable = new ReadHashTable[readSet->readNumber];
    
    for(i = 0; i< totalThreadNumber; i++){
        getReadSetP[i].readSet = readSet;
        getReadSetP[i].threadIndex = i;
        getReadSetP[i].totalThreadNumber = totalThreadNumber;
             
        if(pthread_create(&tid[i], NULL, GetReadSetThread, (void *)&getReadSetP[i])){
            cout<<"create thread wrong!"<<endl;
            exit(1);
        }      
    }
    
    for(i = 0; i < totalThreadNumber; i++){
        pthread_join(tid[i], NULL);
    }
    
}

int CompareReadHashTable(const void * a, const void * b){
    ReadHashTable * ha = (ReadHashTable *)a;
    ReadHashTable * hb = (ReadHashTable *)b;

    return (ha->keyValue - hb->keyValue);
}

void WriteReadHashTable(ReadSet * readSet){
    long int i = 0;
    long int j = 0;
    
    char temp[20];
    strcpy(temp, "readHashTable.fa");
    
    ofstream ocout;
    ocout.open(temp);
    
    char read[readSet->readLength+1];
    
    for(i=0; i<readSet->readNumber;i++){
        GetBit(readSet->readSet, readSet->readLength*readSet->readHashTable[i].index, readSet->readLength, read);
        ocout<<i<<"--"<<read<<","<<readSet->readHashTable[i].index<<","<<readSet->readHashTable[i].keyValue<<endl;
    }
    
    ocout.close();
}

void InitReadSet(ReadSet * readSet, int readSetNumber, unsigned int totalThreadNumber){
    int i = 0;
    unsigned long int allReadNumber = 0;
    for(i = 0; i<readSetNumber; i++){
        readSet[i].readNumber = GetReadNumber(readSet[i].address, totalThreadNumber);
        readSet[i].readLength = GetReadOrKmerLength(readSet[i].address);
        //cout<<readSet[i].readLength<<endl;
    }
       
    for(i = 0; i<readSetNumber; i++){
        GetReadSet(readSet + i, totalThreadNumber);
        qsort(readSet[i].readHashTable, readSet[i].readNumber, sizeof(ReadHashTable), CompareReadHashTable);
        
    }
}



#endif
