#ifndef FILLGAP_H_INCLUDED
#define FILLGAP_H_INCLUDED

#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "readset.h"
#include "kmerset.h"
#include "bitarray.h"
#include "common.h"
#include "graph.h"

using namespace std;

#pragma pack(2)
typedef struct CreateSubDBGraphP{
    char * contig;
    DBGraphHead * deBruijnGraphHead;
    ReadSet * readSet;
    int kmerLength;
    int direction;
    int threadIndex;
    int totalThreadNumber;
    CreateSubDBGraphP(){
        contig = NULL;
        deBruijnGraphHead = NULL;
        readSet = NULL;
        kmerLength = 0;
        threadIndex = 0;
        totalThreadNumber = 0;
        direction = -1;
    }
}CreateSubDBGraphP;
#pragma pack ()

void * CreateSubDBGraphThread(void * arg){
    
    CreateSubDBGraphP * createSubDBGraphP = (CreateSubDBGraphP *)arg;
    
    long int readLength = createSubDBGraphP->readSet->readLength;
    int kmerLength = createSubDBGraphP->kmerLength;
    char * contig = createSubDBGraphP->contig;
    ReadSet * readSet = createSubDBGraphP->readSet;
    int direction = createSubDBGraphP->direction;
    long int insertSize = createSubDBGraphP->readSet->insertSize;
    long int var = createSubDBGraphP->readSet->var;
    
    long int contigLength = strlen(contig);
    long int len = contigLength;
    
    char * tempRead = new char[readLength + 1];
    char * tempReadRC = new char[readLength + 1];
    char * tempKmer = new char[kmerLength + 1];

    long int i =0;
    long int j = 0;
    
    long int q = -1;
    
    long int start = 0;
    long int end = 0;
    //cout<<"bb1"<<endl;
    if(direction = 1){
        end = contigLength - createSubDBGraphP->threadIndex * (len/createSubDBGraphP->totalThreadNumber)  + readLength;
        start = end - (len/createSubDBGraphP->totalThreadNumber);    
    }else if(direction = -1){
        start = createSubDBGraphP->threadIndex * (len/createSubDBGraphP->totalThreadNumber);  
        end = start + (len/createSubDBGraphP->totalThreadNumber) + readLength;
    }
    //cout<<"bb2"<<endl;
    for(i = start; i < end && i < contigLength - readLength; i++){
        //cout<<"tt:"<<i<<endl;
        SubContig(tempKmer, contig, i, i + kmerLength);
        //cout<<q<<"--"<<tempKmer<<"--"<<kmerLength<<"--"<<createSubDBGraphP->deBruijnGraphHead->nodeNumber<<"-start-"<<start<<"-end-"<<end<<"-contigLength-"<<contigLength<<endl;
        q = DBGraphInsertNode(q, tempKmer, kmerLength, createSubDBGraphP->deBruijnGraphHead->deBruijnGraph, createSubDBGraphP->deBruijnGraphHead->nodeNumber); 
        //cout<<"ss"<<endl;
        SubContig(tempRead, contig, i, i+readLength); 
        //ReverseComplement(tempRead, tempReadRC);
        
        ReadMate * tempReadMate = NULL;
        ReadMate * tempReadMateRC = NULL;
        if(createSubDBGraphP->direction == 1){
            tempReadMate = SearchRightReadMate(tempRead, readSet);
            //tempReadMateRC = SearchLeftReadMate(tempReadRC, readSet);
        }else{
            tempReadMate = SearchLeftReadMate(tempRead, readSet);
            //tempReadMateRC = SearchRightReadMate(tempReadRC, readSet);
        }
        //cout<<"tt1"<<endl;       
        ReadMate * temp1 = NULL;
        
        if(tempReadMate!=NULL){
            temp1 = tempReadMate;
            tempReadMate = tempReadMate->next;
            delete temp1;
            temp1 = NULL;
        }
        /*
        if(tempReadMateRC!=NULL){
            temp1 = tempReadMateRC;
            tempReadMateRC = tempReadMateRC->next;
            delete temp1;
            temp1 = NULL;
        }
        */
        while(tempReadMate!=NULL){
            
            long int p = -1;
            //cout<<"readmate:"<<p<<"--"<<tempReadMate->readMate<<endl;
            for(int t = 0; t<readLength - kmerLength + 1; t++){
                SubContig(tempKmer, tempReadMate->readMate, t, t+kmerLength);
                //cout<<"p:"<<p<<endl;
                p = DBGraphInsertNode(p, tempKmer, kmerLength, createSubDBGraphP->deBruijnGraphHead->deBruijnGraph, createSubDBGraphP->deBruijnGraphHead->nodeNumber);
            }
            //cout<<"kk"<<endl;      
            temp1 = tempReadMate;
            tempReadMate = tempReadMate->next;
            delete [] temp1->readMate;
            delete temp1; 
            temp1 = NULL;      
        } 
       
        
    }
    delete [] tempRead;
    delete [] tempReadRC;
    delete [] tempKmer;
    
}


void CreateSubDBGraph(char * contig, ReadSet * readSet, int kmerLength, int setNumber, DBGraphHead * deBruijnGraphHead, int direction, int totalThreadNumber){
    
    pthread_t tid[totalThreadNumber];
    
    int i = 0;
    unsigned long int allKmerNumber = 0;
    CreateSubDBGraphP * createSubDBGraphP = new CreateSubDBGraphP[totalThreadNumber];
    
    for(i = 0; i< totalThreadNumber; i++){
        createSubDBGraphP[i].contig = contig;
        createSubDBGraphP[i].deBruijnGraphHead = deBruijnGraphHead;
        createSubDBGraphP[i].readSet = readSet;
        createSubDBGraphP[i].kmerLength = kmerLength;
        createSubDBGraphP[i].direction = direction;
        createSubDBGraphP[i].threadIndex = i;
        createSubDBGraphP[i].totalThreadNumber = totalThreadNumber;
             
        if(pthread_create(&tid[i], NULL, CreateSubDBGraphThread, (void *)&createSubDBGraphP[i])){
            cout<<"create thread wrong!"<<endl;
            exit(1);
        }      
    }
    
    for(i = 0; i < totalThreadNumber; i++){
        pthread_join(tid[i], NULL);
    }
    
}

double WeightRightGapRegion(char * contigLeft, char * contigRight, char * right, ReadSet * readSet, long int readSetIndex,long int kmerLength, long int start, long int length, long int gapDistance, KmerSetHashTable * kmerSetHashTable, unsigned long int kmerSetHashTableCount){
    long int i = 0;
    long int j = 0;
    long int readLength = readSet[readSetIndex].readLength;
    long int insertSize = readSet[readSetIndex].insertSize;
    long int var = readSet[readSetIndex].var;
    
    long int nullNumber = 0;
    long int positionMatchNumber = 0;
    long int positionNonExist = 0;
    
    long int AllMatchNumber = 0;
  
    
    long int min = 0;

    char * tempContigLeft = contigLeft;
    long int tempLen = strlen(contigLeft);

    char * tempContigRight = new char[readLength + length + 1];

    for(i =0; i<readLength+length; i++){
        if(i<readLength-1){
            tempContigRight[i] = contigLeft[tempLen-readLength + 1 + i];
        }else{
            tempContigRight[i] = contigRight[kmerLength-1 + j];
            j++;
        }
    }
    tempContigRight[i] = '\0';
    
    char * temp = new char[readLength + 1];
    char * tempRC = new char[readLength + 1];
    
    //cout<<"uuuuuuuuuuuu"<<endl;    
    for(j=start;j<start+length;j++){

        SubContig(temp, tempContigRight, j, j + readLength);

        ReadMate * tempReadMate = SearchRightReadMate(temp, readSet+readSetIndex);
        ReverseComplement(temp, tempRC);
        long int a = SearchRead(tempRC, readSet + readSetIndex);
        
        ReadMate * temp1 = NULL;
        
        if(tempReadMate==NULL && a ==-1){
            positionNonExist++;
            continue;
        }
        
        if(tempReadMate!=NULL){
            temp1 = tempReadMate;
            tempReadMate = tempReadMate->next;
            delete temp1;
            temp1 = NULL;
        }  
        
        if(tempReadMate==NULL){
            nullNumber++;
            continue;
        }
        
        
        bool token = false;
        
        while(tempReadMate!=NULL){
            long int index1 = KMPIndexOfContigOfMisMatch(right,tempReadMate->readMate);
            long int tempDD = index1 + readLength;
            if(index1!=-1 && tempDD<insertSize+lambda*var){
                if(token!=true){
                    positionMatchNumber++;
                    token = true;
                }
            }       
            temp1 = tempReadMate;
            tempReadMate = tempReadMate->next;
            delete [] temp1->readMate;
            temp1->next = NULL;
            delete temp1;
            temp1 = NULL;            
            
        }
    }


    
    delete [] temp;
    delete [] tempRC;
    delete [] tempContigRight;  
    
    
    if(positionMatchNumber<=2){
        return 0;
    }
    
    double a = (double)(positionMatchNumber)/(double)(length - nullNumber - positionNonExist);
    
    return a;
    
    
}

int FillingGapRegion(char * contigLeft, long int contigLeftLength, char * contigRight, long int leftLast, long int rightFirst, long int storeLength, long int gapDistance, long int fillGapDistance, DBGraph * deBruijnGraph, long int graphNodeNumber, ReadSet * readSet, long int readSetCount, long int length, long int kmerLength, KmerSetHashTableHead *kmerSetHashTableHead, ContigSet * contigSet, int * repeateNumber){
    
    if(globalFillingPathNum>3 || wrongFillingPathNum >20){
    //if(globalFillingPathNum>20){
        return 0;
    }
    
    if(leftLast == rightFirst && fillGapDistance -kmerLength + 1 >= (long int)(0.5*gapDistance)){
        globalFillingPathNum++;
        return 1;
    }
        
    long int i = 0;
    long int j = 0;
    
    
    int * nodeRepeateNumber = new int[graphNodeNumber];
    for(i = 0;i<graphNodeNumber;i++){
        nodeRepeateNumber[i] = repeateNumber[i];
    }

    //cout<<"hh"<<endl;   
    long int outCount = 0;
    long int inCount = 0;
    
    long int pos = leftLast;
    
    outCount = NodeCount(deBruijnGraph[pos].outNode);
    if(outCount==0){
        //cout<<"first outCount == 0"<<endl;
    }

    while(outCount!=0){
        if(outCount == 1){
            //cout<<"node0"<<endl;           
            if(deBruijnGraph[pos].outNode->index==rightFirst && fillGapDistance -kmerLength + 1 >= (long int)(0.5*gapDistance)){  
                globalFillingPathNum++;
                delete [] nodeRepeateNumber;
                return 1;
            }

            pos = deBruijnGraph[pos].outNode->index;
            nodeRepeateNumber[pos]++;
            if(nodeRepeateNumber[pos]>fillRepeateNumberMax){
                //cout<<"fillRepeateNumberMax!"<<endl;
                wrongFillingPathNum ++;
                delete [] nodeRepeateNumber;
                return 0;
            }     
            //cout<<"node0:"<<pos<<endl;                
            if(contigLeftLength + strlen(deBruijnGraph[pos].contig)-kmerLength + 1 > storeLength){
                
                storeLengthPathNum ++;
                wrongFillingPathNum ++;
                delete [] nodeRepeateNumber;
                return 0;
            }
            contigLeft[contigLeftLength-kmerLength + 1] = '\0';
            strcat(contigLeft, deBruijnGraph[pos].contig);
            contigLeftLength = contigLeftLength + strlen(deBruijnGraph[pos].contig)-kmerLength + 1; 
            fillGapDistance = fillGapDistance + strlen(deBruijnGraph[pos].contig)-kmerLength + 1;
            outCount = NodeCount(deBruijnGraph[pos].outNode);
            //cout<<"node00--"<<outCount<<endl;
            
        }
        if(outCount>1){
            //cout<<"node1"<<endl;
            int singleMatchNumber = 0;
            int matchNumber[outCount];
            double tempWeight[outCount];
            for(i=0;i<outCount;i++){
                tempWeight[i] = 0;
                matchNumber[i] = 0;
            }
            long int tempAdjNode[outCount];
            for(i=0;i<outCount;i++){
                tempAdjNode[i] = -1;
            }
            long int tempPrevious = -1;
            i = 0;
            DFSNode * dfsNode;
            dfsNode = DFSTranverse(deBruijnGraph,pos,kmerLength,length + kmerLength);
            
            DFSNode * tempDFSNode = dfsNode; 
            DFSNode * previousDFSNode = dfsNode; 
            while(dfsNode!=NULL){
                singleMatchNumber = 0;
                if(strlen(dfsNode->contig)<=length + kmerLength){                    
                    if(dfsNode->adjNode!=tempPrevious){
                        if(tempPrevious!=-1){
                            i++;
                        }
                        tempAdjNode[i] = dfsNode->adjNode;
                        tempPrevious = dfsNode->adjNode;                  
                    }                                                       
                    dfsNode = dfsNode->next;
                    continue;
                }else{
                    double weight = 0;
                    double weightNum = 0;
                    for(int p = 0;p<readSetCount;p++){ 
                        //cout<<dfsNode->contig<<endl;                         
                        double tempDD = WeightRightCandidateContig(contigLeft,dfsNode->contig,strlen(deBruijnGraph[pos].contig),readSet,p,kmerLength,0,length,kmerSetHashTableHead->kmerSetHashTable,kmerSetHashTableHead->kmerSetHashTableCount);                       
                        if(tempDD!=-1){
                            weight = weight + tempDD;
                            weightNum++;
                        }
                        
                        if(tempDD>0){
                            singleMatchNumber++;
                        }
                        
                        if(gapDistance - fillGapDistance > readSet[i].insertSize - readSet[i].var){
                            continue;
                        }
                        
                        double tempDD1 = WeightRightGapRegion(contigLeft,dfsNode->contig,contigRight,readSet,p,kmerLength,0,length,storeLength - contigLeftLength, kmerSetHashTableHead->kmerSetHashTable,kmerSetHashTableHead->kmerSetHashTableCount); 
                        if(tempDD1<0.2){
                            weight = 0;
                            break;
                        }
                        
                                   
                    }                   
                    if(weightNum!=0){
                        weight = weight/(double)weightNum;
                    }
                                        
                    if(dfsNode->adjNode!=tempPrevious){
                        if(tempPrevious!=-1){
                            i++;
                        }
                        tempWeight[i] = weight;
                        tempAdjNode[i] =dfsNode->adjNode;  
                        tempPrevious = dfsNode->adjNode;
                        matchNumber[i] = singleMatchNumber;                   
                    }else{
                        if(tempWeight[i]<weight){
                            tempWeight[i] = weight;
                        }
                        if(matchNumber[i]<singleMatchNumber){
                            matchNumber[i] = singleMatchNumber;
                        }
                    }                        
                }
                dfsNode = dfsNode->next;
            }
            //cout<<"node2"<<endl;
            j = 0;
            double maxWeight = 0;
            long int maxPosition = -1;
            
            for(i=0;i<outCount-1;i++){
                for(j=i;j<outCount;j++){
                    if(tempWeight[i]<tempWeight[j]){
                       maxWeight = tempWeight[j];
                       tempWeight[j] = tempWeight[i];
                       tempWeight[i] = maxWeight;
                       maxPosition = tempAdjNode[j];
                       tempAdjNode[j] = tempAdjNode[i];
                       tempAdjNode[i] = maxPosition;
                       singleMatchNumber = matchNumber[j];
                       matchNumber[j] = matchNumber[i];
                       matchNumber[i] = singleMatchNumber;
                    }
                }
            }
                       
            if(matchNumber[0]<matchNumber[1]){
                maxWeight = tempWeight[0];
                tempWeight[0] = tempWeight[1];
                tempWeight[1] = maxWeight;
                maxPosition = tempAdjNode[0];
                tempAdjNode[0] = tempAdjNode[1];
                tempAdjNode[1] = maxPosition;
                singleMatchNumber = matchNumber[0];
                matchNumber[0] = matchNumber[1];
                matchNumber[1] = singleMatchNumber;
            }
            
            //cout<<"1:"<<maxWeight[0]<<"---"<<maxWeight[1]<<endl;
            
            if((tempAdjNode[0]==rightFirst || tempAdjNode[1]==rightFirst) && fillGapDistance - kmerLength + 1>= (long int)(0.5*gapDistance)){  
                globalFillingPathNum++;
                delete [] nodeRepeateNumber;
                return 1;                
            }

            if(tempWeight[0]<=0.2){
                //cout<<"maxWeight[0]<=0.2---"<<maxWeight[0]<<endl;
                delete [] nodeRepeateNumber;
                wrongFillingPathNum ++;
                return 0;
            }
            //cout<<"node3"<<endl;
            //if(matchNumber[0] > matchNumber[1] || tempWeight[0] - tempWeight[1]>0.2){
            if(tempWeight[0] - tempWeight[1]>0.2){             
                //cout<<"node4"<<endl;
                pos = tempAdjNode[0];
                if(contigLeftLength + strlen(deBruijnGraph[pos].contig)-kmerLength + 1>storeLength){
                    
                    storeLengthPathNum ++;
                    wrongFillingPathNum ++;
                    delete [] nodeRepeateNumber;
                    return 0;
                }
                nodeRepeateNumber[pos]++;
                if(nodeRepeateNumber[pos]>fillRepeateNumberMax){
                    //cout<<"fillRepeateNumberMax!"<<endl;
                    wrongFillingPathNum ++;
                    delete [] nodeRepeateNumber;
                    return 0;
                } 
                //cout<<"node1:"<<pos<<endl; 
                contigLeft[contigLeftLength-kmerLength + 1] = '\0';
                strcat(contigLeft, deBruijnGraph[pos].contig);
                contigLeftLength = contigLeftLength + strlen(deBruijnGraph[pos].contig)-kmerLength + 1;     
                fillGapDistance = fillGapDistance + strlen(deBruijnGraph[pos].contig)-kmerLength + 1;           
                outCount = NodeCount(deBruijnGraph[pos].outNode); 
                //cout<<"node5"<<endl;
                              
            }else{
                /*
                cout<<"node6"<<endl;
                delete [] nodeRepeateNumber;
                return 0;
                */
                if(tempAdjNode[1]!=-1&&tempAdjNode[0]!=-1){
                    if(tempWeight[0]>0.2){
                        pos = tempAdjNode[0];
                        if(contigLeftLength + strlen(deBruijnGraph[pos].contig)-kmerLength + 1>storeLength){
                            
                            storeLengthPathNum ++;
                            wrongFillingPathNum ++;
                            delete [] nodeRepeateNumber;
                            return 0;
                        }
                        nodeRepeateNumber[pos]++;
                        if(nodeRepeateNumber[pos]>fillRepeateNumberMax){
                            //cout<<"fillRepeateNumberMax!"<<endl;
                            delete [] nodeRepeateNumber;
                            wrongFillingPathNum ++;
                            return 0;
                        } 
                        //cout<<"qq:"<<storeLength<<"--"<<contigLeftLength + strlen(deBruijnGraph[pos].contig)-kmerLength + 1<<endl;
                        char * tempContigLeft = new char[storeLength+1];
                        CopyContig(tempContigLeft, contigLeft);
                        tempContigLeft[contigLeftLength - kmerLength + 1] = '\0';
                        strcat(tempContigLeft, deBruijnGraph[pos].contig);
                        contigLeftLength = contigLeftLength + strlen(deBruijnGraph[pos].contig)-kmerLength + 1;
                        fillGapDistance = fillGapDistance + strlen(deBruijnGraph[pos].contig)-kmerLength + 1;
                        //cout<<"tt--"<<pos<<endl; 
                        int tt = FillingGapRegion(tempContigLeft, contigLeftLength, contigRight,pos, rightFirst, storeLength, gapDistance, fillGapDistance, deBruijnGraph, graphNodeNumber, readSet, readSetCount, length, kmerLength, kmerSetHashTableHead, contigSet, nodeRepeateNumber);
                        if(tt == 1){
                            ContigSet * tempContigSet = new ContigSet;
                            tempContigSet->contig = tempContigLeft;
                            tempContigSet->next = contigSet->next;
                            contigSet->next = tempContigSet;
                        }else{
                            delete [] tempContigLeft;
                        }
                    }
                    
                    if(tempWeight[1]>0.2){        
                        pos = tempAdjNode[1];
                    
                        if(contigLeftLength + strlen(deBruijnGraph[pos].contig)-kmerLength + 1>storeLength){
                            
                            storeLengthPathNum ++;
                            wrongFillingPathNum ++;
                            delete [] nodeRepeateNumber;
                            return 0;
                        }
                        nodeRepeateNumber[pos]++;
                        if(nodeRepeateNumber[pos]>fillRepeateNumberMax){
                            //cout<<"fillRepeateNumberMax!"<<endl;
                            delete [] nodeRepeateNumber;
                            wrongFillingPathNum ++;
                            return 0;
                        }
                        char * tempContigLeft = new char[storeLength+1];
                        CopyContig(tempContigLeft, contigLeft);
                        tempContigLeft[contigLeftLength - kmerLength + 1] = '\0';
                        strcat(tempContigLeft, deBruijnGraph[pos].contig);
                        contigLeftLength = contigLeftLength + strlen(deBruijnGraph[pos].contig)-kmerLength + 1;
                        fillGapDistance = fillGapDistance + strlen(deBruijnGraph[pos].contig)-kmerLength + 1;
                        //cout<<"tt1--"<<pos<<endl; 
                        int tt = FillingGapRegion(tempContigLeft,contigLeftLength,contigRight,pos, rightFirst, storeLength, gapDistance, fillGapDistance, deBruijnGraph, graphNodeNumber, readSet, readSetCount, length, kmerLength, kmerSetHashTableHead, contigSet, nodeRepeateNumber);
                        if(tt == 1){
                            ContigSet * tempContigSet = new ContigSet;
                            tempContigSet->contig = tempContigLeft;
                            tempContigSet->next = contigSet->next;
                            contigSet->next = tempContigSet;
                        }else{
                            delete [] tempContigLeft;
                        }
                        
                    }                   
                }
                //cout<<"path many"<<endl;
                delete [] nodeRepeateNumber;
                
                return 0;
            }
        }
        if(outCount == 0){
            //cout<<"outCount == 0"<<endl;
        }
    }
    delete [] nodeRepeateNumber;
    return 0;    
}

char * GetOptimalGapRegion(ContigSet * contigSet, long int gapDistance){
    long int i = 0;
    long int j = 0;
    
    ContigSet * temp = contigSet;
    
    long int min = labs(strlen(temp->contig) - gapDistance);
    
    char * contig = temp->contig;
    
    while(temp!=NULL){
        if(min > labs(strlen(temp->contig) - gapDistance)){
            min = labs(strlen(temp->contig) - gapDistance);
            j = i;
            contig = temp->contig;            
        }
        i++;
        temp = temp->next;
    }
    
    char * result = new char[strlen(contig) + 1];
    CopyContig(result, contig);
    return result;
}


#pragma pack(2)
typedef struct GapRegionContig{
    ContigSet * contigSet;
    char * right;
    long int length;
    int kmerLength;
    long int overlap;
    long int gapDistance;
    GapRegionContig(){
        contigSet = NULL;
        right = NULL;
        length = 0;
        kmerLength = 0;
        overlap = 0;
        gapDistance = 0;
    }
}GapRegionContig;
#pragma pack ()

#pragma pack(2)
typedef struct GapRegionContigSet{
    GapRegionContig * gapRegionContig;
    GapRegionContigSet(){
        gapRegionContig = NULL;
    }
}GapRegionContigSet;
#pragma pack()

#pragma pack(2)
typedef struct FillGapP{
    
    ScaffoldSetHead * scaffoldSetHead;
    ReadSet * readSet;
    int setNumber;
    KmerSetHashTableHead * kmerSetHashTableHead;
    int kmerLength;
    GapRegionContigSet * gapRegionContigSet;
    int threadIndex;
    int totalThreadNumber;   
    
}FillGapP;
#pragma pack ()

GapRegionContig * SubFillGap(Scaffold * leftScaffold, Scaffold * rightScaffold, ReadSet * readSet, int setNumber, KmerSetHashTableHead * kmerSetHashTableHead, int kmerLength, int totalThreadNumber){

    long int i = 0;
    long int j = 0;
    
    char * leftContig = leftScaffold->contig;
    char * rightContig = rightScaffold->contig;
    long int gapDistance = leftScaffold->gapDistance;
    long int var = 0;
    long int insertSize = 0;
    long int readSetIndex = -1;
    
    long int nodeNumber = 0;
    for(i = 0; i<setNumber; i++){
        nodeNumber = nodeNumber + readSet[i].insertSize + lambda*readSet[i].var;
        if(gapDistance<readSet[i].insertSize - lambda*readSet[i].var && readSetIndex == -1){
            readSetIndex = i;
        }
        if(insertSize<readSet[i].insertSize){
            insertSize = readSet[i].insertSize;
            var = readSet[i].var;
        }
    }
    nodeNumber = 150*nodeNumber;
    
    KmerSetHashTable * kmerSetHashTable = kmerSetHashTableHead->kmerSetHashTable;
    long int kmerSetHashTableCount = kmerSetHashTableHead->kmerSetHashTableCount;
    
    long int presentKmerLength = minKmerLength;
    DBGraphHead * deBruijnGraphHead = new DBGraphHead;
    deBruijnGraphHead->deBruijnGraph = new DBGraph[nodeNumber];
    deBruijnGraphHead->nodeNumber = nodeNumber;
    int * graphNodeRepeateNumber = new int[nodeNumber];
    for(i = 0; i<nodeNumber; i++){
        graphNodeRepeateNumber[i] = 0;
    }
    
    long int leftContigLength = strlen(leftContig);
    long int len = leftContigLength;
    if(len>insertSize+2*lambda*var){
        len = insertSize+2*lambda*var;
    }
    long int storeLength = len + (long int)(2*leftScaffold->gapDistance) + 2*cut + var;
    char * tempContigLeft1 = new char[storeLength + 1];
    SubContig(tempContigLeft1, leftContig, leftContigLength - len, leftContigLength - cut);
    leftContigLength = len - cut;
    //cout<<"leftLength:"<<leftContigLength<<endl;
    long int rightContigLength = strlen(rightContig);
    len = rightContigLength;
    if(len>insertSize+2*lambda*var){
        len = insertSize+2*lambda*var;
    }
    char * tempContigRight1 = new char[len + 1];
    SubContig(tempContigRight1, rightContig, cut, len); 
    rightContigLength = len - cut;
    //cout<<"rightLength:"<<rightContigLength<<endl;
    ContigSet * contigSet = new ContigSet;  
    long int length = globalKmerLength;
    
    long int leftLengthFillGap = leftContigLength;
    long int rightLengthFillGap = rightContigLength;
    
    char * tempRightCut = NULL;
    char * tempContigLeft = NULL;
    char * tempContigRight = NULL;
    
    long int overlap = 0;
    
    long int fillGapDistance = 0;
    
    for(presentKmerLength = maxKmerLength; presentKmerLength>=minKmerLength; presentKmerLength = presentKmerLength - 2){
        fillGapDistance = 0;
        leftLengthFillGap = leftContigLength;
        rightLengthFillGap = rightContigLength;
        tempContigLeft = new char[storeLength + 1];
        CopyContig(tempContigLeft, tempContigLeft1);    
        tempContigRight = new char[rightContigLength + 1];
        CopyContig(tempContigRight, tempContigRight1);       
        //cout<<"jj"<<endl;
        for(i=0;i<setNumber;i++){
            CreateSubDBGraph(tempContigLeft, readSet + i, presentKmerLength, setNumber, deBruijnGraphHead, 1, 1);
            //cout<<"ff0"<<endl;
            CreateSubDBGraph(tempContigRight, readSet + i, presentKmerLength, setNumber, deBruijnGraphHead, -1, 1);
            //cout<<"ff"<<endl;
        }       
        //cout<<"jj1"<<endl;
        for(i=0;i<5;i++){
            RemoveTip(deBruijnGraphHead->deBruijnGraph,presentKmerLength,deBruijnGraphHead->nodeNumber);
            SimplePathMerge(deBruijnGraphHead->deBruijnGraph,presentKmerLength,deBruijnGraphHead->nodeNumber);
        }
        RemoveCycle(deBruijnGraphHead->deBruijnGraph,presentKmerLength,deBruijnGraphHead->nodeNumber, kmerSetHashTableHead);
        RemoveCycle1(deBruijnGraphHead->deBruijnGraph,presentKmerLength,deBruijnGraphHead->nodeNumber);
        SimplePathMerge(deBruijnGraphHead->deBruijnGraph,presentKmerLength,deBruijnGraphHead->nodeNumber);
        AdjustDBG(deBruijnGraphHead->deBruijnGraph,presentKmerLength,deBruijnGraphHead->nodeNumber);
        //deBruijnGraphHead->deBruijnGraph = OptimizeDBGraphSpace(deBruijnGraphHead->deBruijnGraph, deBruijnGraphHead->nodeNumber);
        
        //char * tempAddress = new char[30];
        //strcpy(tempAddress,"subDBGraph.fa");
        //WriteDBGraph(deBruijnGraphHead->deBruijnGraph, deBruijnGraphHead->nodeNumber, tempAddress);
        
                
        //cout<<"--"<<tempContigLeft<<endl;
        //cout<<"aa:"<<presentKmerLength<<"--"<<leftContigLength - presentKmerLength<<"--"<<leftContigLength<<endl;
        char * tempLeft = new char[presentKmerLength + 1];
        SubContig(tempLeft, tempContigLeft, leftContigLength - presentKmerLength, leftContigLength);
        //cout<<"index:"<<tempLeft<<endl;
        long int * leftLast = SearchDBGraph(deBruijnGraphHead, tempLeft);
        //cout<<"aa1"<<endl;
        char * tempRight = new char[presentKmerLength + 1];
        SubContig(tempRight, tempContigRight, 0, presentKmerLength);
        //cout<<"index1:"<<tempRight<<endl;
        long int * rightFirst = SearchDBGraph(deBruijnGraphHead, tempRight);
        //cout<<"aa2"<<endl;
        if(leftLast == NULL || rightFirst == NULL){
            
            DeleteDeBruijnGraph(deBruijnGraphHead);
            continue;
        }        
        //cout<<"aa3"<<endl;
        
        //if(leftLast[0] == rightFirst[0] && leftScaffold->gapDistance < var){
        if(leftLast[0] == rightFirst[0]){
            
            if(leftLast[1]+ presentKmerLength > rightFirst[1]){
                overlap = leftLast[1] + presentKmerLength - rightFirst[1];
                contigSet->contig = tempContigLeft;
                tempContigLeft = NULL;
            }else{
                tempRightCut = new char[presentKmerLength];
                SubContig(tempRightCut, deBruijnGraphHead->deBruijnGraph[leftLast[0]].contig, rightFirst[1] - presentKmerLength + 1, rightFirst[1]);
                char * tempCut = new char[rightFirst[1] - leftLast[1] - presentKmerLength + 1];
                SubContig(tempCut, deBruijnGraphHead->deBruijnGraph[leftLast[0]].contig, leftLast[1] + presentKmerLength, rightFirst[1]);
                tempContigLeft = (char *)realloc(tempContigLeft, leftContigLength + rightFirst[1] - leftLast[1] - presentKmerLength + 1);
                strcat(tempContigLeft, tempCut);
                contigSet->contig = tempContigLeft;
                tempContigLeft = NULL;
                delete [] tempCut;
                tempCut = NULL;
            }
            DeleteDeBruijnGraph(deBruijnGraphHead);
            break;            
        }
                
        //cout<<"kk--"<<leftLast[2]<<"--"<<leftLast[1]<<"--"<<leftLast[0]<<"--"<<presentKmerLength<<"--"<<leftContigLength<<endl;
        if(leftLast[1] + presentKmerLength != leftLast[2]){
            //cout<<"kk--"<<leftLast[2] - leftLast[1] - presentKmerLength + 1<<endl;
            fillGapDistance = fillGapDistance + leftLast[2] - leftLast[1] - presentKmerLength;
            char * tempCut = new char[leftLast[2] - leftLast[1] - presentKmerLength + 1];
            //cout<<"kk1"<<endl;
            SubContig(tempCut, deBruijnGraphHead->deBruijnGraph[leftLast[0]].contig, leftLast[1] + presentKmerLength, leftLast[2]);
            //cout<<"kk2"<<endl;
            if(leftContigLength + leftLast[2] - leftLast[1] - presentKmerLength>storeLength){
                tempContigLeft = (char *)realloc(tempContigLeft, leftContigLength + leftLast[2] - leftLast[1] - presentKmerLength + var);
            }
            strcat(tempContigLeft, tempCut);
            //cout<<"kk3"<<endl;  
            leftLengthFillGap = leftContigLength + leftLast[2] - leftLast[1] - presentKmerLength;          
        }     
        
        //cout<<"aa4"<<endl;
        
        //cout<<"kk1--"<<rightFirst[2]<<"--"<<rightFirst[1]<<"--"<<rightFirst[0]<<"--"<<presentKmerLength<<"--"<<rightContigLength<<endl;
        if(rightFirst[1] != 0){
            fillGapDistance = fillGapDistance + rightFirst[1];
            char * tempCut = new char[rightFirst[1] + rightContigLength + 1];
            tempRightCut = new char[rightFirst[1] + 1];
            SubContig(tempRightCut, deBruijnGraphHead->deBruijnGraph[rightFirst[0]].contig, 0, rightFirst[1]);
            SubContig(tempCut, deBruijnGraphHead->deBruijnGraph[rightFirst[0]].contig, 0, rightFirst[1]);
            strcat(tempCut, tempContigRight); 
            delete [] tempContigRight;
            tempContigRight = tempCut;
            tempCut = NULL;
            rightLengthFillGap = rightContigLength + rightFirst[1];           
        }    
        
        //cout<<"aa5--"<<leftLast[0]<<"--"<<rightFirst[0]<<endl;
        globalFillingPathNum = 0;
        fillGapDistance = 0;
        storeLengthPathNum = 0;
        wrongFillingPathNum = 0;
        int a = FillingGapRegion(tempContigLeft, leftLengthFillGap, tempContigRight, leftLast[0], rightFirst[0], storeLength, gapDistance + 2*cut, fillGapDistance, deBruijnGraphHead->deBruijnGraph, deBruijnGraphHead->nodeNumber,readSet, setNumber, length, presentKmerLength, kmerSetHashTableHead, contigSet, graphNodeRepeateNumber);
        //cout<<"tt"<<endl;
        DeleteDeBruijnGraph(deBruijnGraphHead);
        //cout<<"aa6--"<<a<<endl;       
        
        if(a == 1){
            contigSet->contig = tempContigLeft;
            tempContigLeft = NULL;
        }
        if(contigSet->next!=NULL || contigSet->contig!=NULL){           
            
            break;
            
        }
        //cout<<"tt1"<<endl;
        delete [] tempContigLeft;
        tempContigLeft = NULL;
        delete [] tempContigRight;
        tempContigRight = NULL;
        delete [] tempRightCut;
        tempRightCut = NULL;
        //cout<<"tt2"<<endl;
                 
    }
    
    //exit(0);
    //cout<<"tt3"<<endl;
    GapRegionContig * result = NULL;
    if(contigSet->contig!=NULL){
        result = new GapRegionContig;
        result->contigSet = contigSet;
        result->right = tempRightCut;
        result->length = leftContigLength;
        result->kmerLength = presentKmerLength;
        result->overlap = overlap;
        result->gapDistance = gapDistance;
    }else if(contigSet->next!=NULL){
        result = new GapRegionContig;
        result->contigSet = contigSet->next;
        result->right = tempRightCut;
        result->length = leftContigLength;
        result->kmerLength = presentKmerLength;
        result->overlap = overlap;
        result->gapDistance = gapDistance;
    }
    //cout<<"tt4"<<endl;
    delete [] deBruijnGraphHead->deBruijnGraph;
    delete [] tempContigLeft;
    delete [] tempContigRight;
    delete [] tempContigLeft1;
    delete [] tempContigRight1;
    //cout<<"tt5"<<endl;
    return result;    
    
}


void * FillGapThread(void * arg){
    
    FillGapP * fillGapP = (FillGapP *)arg;
    
    ScaffoldSet * scaffoldSet = fillGapP->scaffoldSetHead->scaffoldSet;
    long int scaffoldSetNumber = fillGapP->scaffoldSetHead->scaffoldSetNumber;
    long int gapNumber = fillGapP->scaffoldSetHead->gapNumber;
    ReadSet * readSet = fillGapP->readSet;
    int setNumber = fillGapP->setNumber;
    KmerSetHashTableHead * kmerSetHashTableHead = fillGapP->kmerSetHashTableHead;
    int kmerLength = fillGapP->kmerLength;
    GapRegionContigSet * gapRegionContigSet = fillGapP->gapRegionContigSet;
    int threadIndex = fillGapP->threadIndex;
    int totalThreadNumber = fillGapP->totalThreadNumber;
    
    //cout<<"--thread:"<<threadIndex<<endl;
    
    long int i = 0;
    long int j = 0;
    long int count = 0;
    
    for(i = 0; i < scaffoldSetNumber; i++){
        Scaffold * leftScaffold = scaffoldSet[i].scaffold;
        Scaffold * rightScaffold = leftScaffold->next; 
        while(rightScaffold!=NULL){
            //cout<<"gapCount:"<<count<<"--"<<threadIndex<<endl;
            if(count%totalThreadNumber != threadIndex){
                leftScaffold = rightScaffold;
                rightScaffold = rightScaffold->next;
                count++;
                continue;
            }
            gapRegionContigSet[count].gapRegionContig = 
                SubFillGap(leftScaffold, rightScaffold, readSet, setNumber, kmerSetHashTableHead, kmerLength, totalThreadNumber);
            count++;
            leftScaffold = rightScaffold;
            rightScaffold = rightScaffold->next;
        }       
    }
    //cout<<"--threadttttttttt:"<<threadIndex<<endl;
}

char * MergeGapContig(char * left, char * mid, char * rightCut, char * right, long int overlap, long int overlap1){
     
     long int i = 0;
     long int j = 0;
     
     long int leftLength = strlen(left);
     long int midLength = strlen(mid);
     long int rightLength = strlen(right);
     
     char * result = NULL;
     
     if(rightCut == NULL){
         result = (char *)realloc(left, leftLength + midLength + rightLength - overlap - overlap1 - 2*cut + 1);
         result[leftLength - overlap - cut] = '\0';
         strcat(result, mid);
         long int start = leftLength + midLength - overlap - overlap1 - cut;
         for(i=cut;i<rightLength;i++){
             result[start] = right[i];
             start++;
         }
         result[start] = '\0';
     }else{
         long int rightCutLength = strlen(rightCut);
         result = (char *)realloc(left, leftLength + midLength + rightCutLength + rightLength - overlap - overlap1 - 2*cut + 1);
         result[leftLength - overlap - cut] = '\0';
         strcat(result, mid);
         result[leftLength + midLength - overlap - overlap1 - cut] = '\0';
         rightCut = (char *)realloc(rightCut, rightCutLength + rightLength - cut + 1);
         long int start = rightCutLength;
         for(i=cut;i<rightLength;i++){
             rightCut[start] = right[i];
             start++;
         }
         rightCut[start] = '\0';
         strcat(result,rightCut);
     }    
     
     return result;
         
}

void MergeGapRegion(ScaffoldSetHead * scaffoldSetHead, GapRegionContigSet * gapRegionContigSet){
    
    long int i = 0;
    long int j = 0;
    
    for(i = 0; i < scaffoldSetHead->scaffoldSetNumber; i++){
        //cout<<"dd--"<<i<<"--"<<scaffoldSetHead->scaffoldSetNumber<<endl;
        Scaffold * temp = scaffoldSetHead->scaffoldSet[i].scaffold;
        while(temp->next!=NULL){
            if(gapRegionContigSet[j].gapRegionContig!=NULL){
                //cout<<"ddee--"<<j<<endl;
                long int length = gapRegionContigSet[j].gapRegionContig->length;
                //long int kmerLength = gapRegionContigSet[j].gapRegionContig->kmerLength;
                //char * mid = GetOptimalGapRegion(gapRegionContigSet[j].gapRegionContig->contigSet, length + temp->gapDistance + kmerLength -1);
                long int presentKmerLength = gapRegionContigSet[j].gapRegionContig->kmerLength;
                ContigSet * lastContigSet = gapRegionContigSet[j].gapRegionContig->contigSet;
                while(lastContigSet->next!=NULL){
                    lastContigSet = lastContigSet->next;
                }
                char * mid = lastContigSet->contig;
                //cout<<"dd2--"<<length<<"--"<<presentKmerLength<<endl;
                if(gapRegionContigSet[j].gapRegionContig->overlap > 0){
                    long int overlap = gapRegionContigSet[j].gapRegionContig->overlap;
                    temp->contig = (char *)realloc(temp->contig, strlen(temp->contig) + strlen(temp->next->contig) - overlap - 2*cut + 1);
                    temp->contig[strlen(temp->contig) - cut - overlap] = '\0';
                    char * tempContig = new char[strlen(temp->next->contig) - cut + 1];
                    SubContig(tempContig, temp->next->contig, cut, strlen(temp->next->contig));
                    strcat(temp->contig, tempContig);
                    delete [] tempContig;
                    tempContig = NULL;                   
                    
                }else{
                    temp->contig = MergeGapContig(temp->contig, mid, gapRegionContigSet[j].gapRegionContig->right, temp->next->contig, length, presentKmerLength-1);
                }               
                //cout<<"dd3"<<endl;
                Scaffold * temp1 = temp->next;
                temp->next = temp->next->next;
                delete [] temp1->contig;
                temp1->contig = NULL;
                temp1->next = NULL;
                delete temp1;
                j++;
                //cout<<"dd4"<<endl;
                continue;                
            }
            //cout<<"dd5"<<endl;
            j++;
            temp = temp->next;            
        }
    }  
}

void MergeGapRegionN(ScaffoldSetHead * scaffoldSetHead, GapRegionContigSet * gapRegionContigSet){
    
    long int i = 0;
    long int j = 0;
    
    for(i = 0; i < scaffoldSetHead->scaffoldSetNumber; i++){
        //cout<<"dd--"<<i<<"--"<<scaffoldSetHead->scaffoldSetNumber<<endl;
        Scaffold * temp = scaffoldSetHead->scaffoldSet[i].scaffold;
        long int gapDistance = temp->gapDistance;
        while(temp->next!=NULL){
            if(gapRegionContigSet[j].gapRegionContig!=NULL){
                //cout<<"ddee--"<<j<<endl;
                long int length = gapRegionContigSet[j].gapRegionContig->length;
                //long int kmerLength = gapRegionContigSet[j].gapRegionContig->kmerLength;
                //char * mid = GetOptimalGapRegion(gapRegionContigSet[j].gapRegionContig->contigSet, length + temp->gapDistance + kmerLength -1);
                long int presentKmerLength = gapRegionContigSet[j].gapRegionContig->kmerLength;
                ContigSet * lastContigSet = gapRegionContigSet[j].gapRegionContig->contigSet;
                
                if(lastContigSet->next!=NULL){
                    char * tempN = new char[gapDistance + 1];
                    int ss = 0;
                    for(ss = 0; ss < gapDistance; ss++){
                        tempN[ss] = 'N';
                    }
                    tempN[ss] = '\0';
                    temp->contig = (char *)realloc(temp->contig, strlen(temp->contig) + strlen(temp->next->contig) + gapDistance +1);
                    strcat(temp->contig, tempN);
                    strcat(temp->contig, temp->next->contig);
                    delete [] tempN;
                    tempN = NULL;
                    
                }else{
                    while(lastContigSet->next!=NULL){
                        lastContigSet = lastContigSet->next;
                    }
                    char * mid = lastContigSet->contig;
                    //cout<<"dd2--"<<length<<"--"<<presentKmerLength<<endl;
                    if(gapRegionContigSet[j].gapRegionContig->overlap > 0){
                        long int overlap = gapRegionContigSet[j].gapRegionContig->overlap;
                        temp->contig = (char *)realloc(temp->contig, strlen(temp->contig) + strlen(temp->next->contig) - overlap - 2*cut + 1);
                        temp->contig[strlen(temp->contig) - cut - overlap] = '\0';
                        char * tempContig = new char[strlen(temp->next->contig) - cut + 1];
                        SubContig(tempContig, temp->next->contig, cut, strlen(temp->next->contig));
                        strcat(temp->contig, tempContig);
                        delete [] tempContig;
                        tempContig = NULL;                   
                    
                    }else{
                        temp->contig = MergeGapContig(temp->contig, mid, gapRegionContigSet[j].gapRegionContig->right, temp->next->contig, length, presentKmerLength-1);
                    }
                }               
                //cout<<"dd3"<<endl;
                //cout<<"dd4"<<endl;
                //continue;                
            }else{
                char * tempN = new char[gapDistance + 1];
                int ss = 0;
                for(ss = 0; ss < gapDistance; ss++){
                    tempN[ss] = 'N';
                }
                tempN[ss] = '\0';
                temp->contig = (char *)realloc(temp->contig, strlen(temp->contig) + strlen(temp->next->contig) + gapDistance +1);
                strcat(temp->contig, tempN);
                strcat(temp->contig, temp->next->contig);
                delete [] tempN;
                tempN = NULL;
            }
            Scaffold * temp1 = temp->next;
            gapDistance = temp1->gapDistance;
            temp->next = temp->next->next;
            delete [] temp1->contig;
            temp1->contig = NULL;
            temp1->next = NULL;
            delete temp1;
            j++;
            //cout<<"dd5"<<endl;
            //j++;
            //temp = temp->next;            
        }
    }  
}

ContigSet * GetContigSetFromScaffoldSetHead(ScaffoldSetHead * scaffoldSetHead){

    long int i = 0;
    long int j = 0;
    
    ContigSet * contigSet = new ContigSet;
    ContigSet * temp = NULL;
    
    for(i = 0; i<scaffoldSetHead->scaffoldSetNumber; i++){
        Scaffold * tempScaffold = scaffoldSetHead->scaffoldSet[i].scaffold;
        while(tempScaffold!=NULL){
            temp = new ContigSet;
            temp->contig = tempScaffold->contig;
            tempScaffold->contig = NULL;
            temp->next = contigSet->next;
            contigSet->next = temp;
            tempScaffold = tempScaffold->next;
        }
    }
    
    return contigSet;

}

void DeleteScaffoldSetHead(ScaffoldSetHead * scaffoldSetHead){

    long int i = 0;
    long int j = 0;
    
    for(i = 0; i<scaffoldSetHead->scaffoldSetNumber; i++){
        Scaffold * tempScaffold = scaffoldSetHead->scaffoldSet[i].scaffold;
        Scaffold * previous = tempScaffold;
        while(tempScaffold!=NULL){
            previous = tempScaffold;
            delete tempScaffold->contig;
            tempScaffold = tempScaffold->next;
            previous->next = NULL;
            delete previous;
            previous = NULL;
        }
        scaffoldSetHead->scaffoldSet[i].scaffold = NULL;
    }
    
    delete [] scaffoldSetHead->scaffoldSet;
    scaffoldSetHead->scaffoldSet = NULL;
    delete scaffoldSetHead;
    scaffoldSetHead = NULL;

}

void WriteGapRegionSet(GapRegionContigSet * gapRegionContigSet, int gapRegionNumber){
    
    long int i = 0;
    long int j = 0;
    
    ofstream ocout;
    char * address = new char[30];
    strcpy(address, "gapRegionSet.fa");
    ocout.open(address);
    
    for(i=0;i<gapRegionNumber;i++){
        if(gapRegionContigSet[i].gapRegionContig == NULL){
            continue;
        }
        ContigSet * temp = gapRegionContigSet[i].gapRegionContig->contigSet;
        while(temp!=NULL){
            ocout<<">"<<i<<endl;
            ocout<<temp->contig<<endl;
            temp = temp->next; 
        }
    }
    
    delete [] address;
    ocout.close();
    
}

int ComputeResult(ScaffoldSetHead * scaffoldSetHead){
    long int i = 0;
    long int j = 0;
    long int count = 0;
    long int tempLength = 0;
    long int max = 0;
    long int n50 = 0;
    long int n90 = 0;
    long int totalLength = 0;
    
    
    for(i = 0; i<scaffoldSetHead->scaffoldSetNumber; i++){
        Scaffold * tempScaffold = scaffoldSetHead->scaffoldSet[i].scaffold;
        while(tempScaffold!=NULL){
            count ++;
            tempScaffold = tempScaffold->next;
        }
    }
    
    if(count==0){
        return 0;
    }
    long int * contigSetLength = new long int[count];
    long int t = 0;
    
    for(i = 0; i<scaffoldSetHead->scaffoldSetNumber; i++){
        Scaffold * tempScaffold = scaffoldSetHead->scaffoldSet[i].scaffold;
        while(tempScaffold!=NULL){
            contigSetLength[j] = strlen(tempScaffold->contig);
            totalLength = totalLength + contigSetLength[j];
            tempScaffold = tempScaffold->next;
            j++;
        }
    }
    
    i = 0;
    j = 0;
    
    for(i = 0; i<count-1; i++){
        for(j = i+1;j<count;j++){
            if(contigSetLength[i]<contigSetLength[j]){
                tempLength = contigSetLength[j];
                contigSetLength[j] = contigSetLength[i];
                contigSetLength[i] = tempLength;
            }
        }
    }
    max = contigSetLength[0];
    tempLength = 0;
    for(i = 0;i<count;i++){
        tempLength = tempLength + contigSetLength[i];
        if(tempLength>totalLength/2){
            if(n50==0){
                n50 = contigSetLength[i];
            } 
        }
        if(tempLength>totalLength*0.9){
            if(n90==0){
                n90 = contigSetLength[i];
            }
        }
    }
    
    char * lastResultStatics = new char[100];
    strcpy(lastResultStatics, "lastResultStatics.fa");
    ofstream ocout;
    ocout.open(lastResultStatics,ios::app);
    cout<<"num:"<<count<<"--maxLength:"<<max<<"--N50:"<<n50<<"--N90:"<<n90<<"--averageLength:"<<tempLength/count<<"--allLength:"<<totalLength<<endl;
    ocout<<"num:"<<count<<"--maxLength:"<<max<<"--N50:"<<n50<<"--N90:"<<n90<<"--averageLength:"<<tempLength/count<<"--allLength:"<<totalLength<<endl;
     
}


void FillGap(ScaffoldSetHead * scaffoldSetHead, ReadSet * readSet, int setNumber, KmerSetHashTableHead * kmerSetHashTableHead, int kmerLength, int totalThreadNumber){
    
    pthread_t tid[totalThreadNumber];
    
    long int i = 0;
    long int j = 0;
    
    FillGapP * fillGapP = new FillGapP[totalThreadNumber];
    GapRegionContigSet * gapRegionContigSet = new GapRegionContigSet[scaffoldSetHead->gapNumber];
    
    for(i = 0; i<totalThreadNumber; i++){
        //cout<<"start--"<<i<<endl;
        fillGapP[i].scaffoldSetHead = scaffoldSetHead;
        fillGapP[i].readSet = readSet;
        fillGapP[i].setNumber = setNumber;
        fillGapP[i].kmerSetHashTableHead = kmerSetHashTableHead;
        fillGapP[i].kmerLength = kmerLength;
        fillGapP[i].gapRegionContigSet = gapRegionContigSet;
        fillGapP[i].threadIndex = i;
        fillGapP[i].totalThreadNumber = totalThreadNumber;
        
        if(pthread_create(&tid[i], NULL, FillGapThread, (void *)&fillGapP[i])){
            cout<<"create thread wrong!"<<endl;
            exit(1);
        } 

    }
    
    for(i = 0; i < totalThreadNumber; i++){
        pthread_join(tid[i], NULL);
    }
    
    
    //WriteGapRegionSet(gapRegionContigSet, scaffoldSetHead->gapNumber);
    
    MergeGapRegionN(scaffoldSetHead,gapRegionContigSet);
   
    char * address = new char[30];
    strcpy(address, "scaffold.fa");
    WriteScaffoldSet(scaffoldSetHead, address);
    strcpy(address, "scaffoldLong.fa");
    WriteScaffoldSetLong(scaffoldSetHead, address);
    
}




#endif
