#ifndef CONSTRUCTCONTIGSET_H_INCLUDED 
#define CONSTRUCTCONTIGSET_H_INCLUDED 

#include <cstring>
#include <iostream>
#include "stdio.h"
#include "stdlib.h"
#include <math.h>
#include <fstream>

#include "readset.h"
#include "kmerset.h"
#include "common.h"
#include "graph.h"
#include "dfs.h" 
#include "contigmerge.h"

using namespace std;

#pragma pack(2)
typedef struct Contig{
    char * contig;
    Contig(){
        contig = NULL;
    }
}Contig;
#pragma pack ()

#pragma pack(2)
typedef struct ContigSet{
    char * contig;
    struct ContigSet * next;
    ContigSet(){
        contig = NULL;
        next = NULL;
    }
}ContigSet;
#pragma pack ()

#pragma pack(2)
typedef struct ContainedNode{
    int index;
    long int length;
    ContainedNode(){
        index = -1;
        length = 0;;
    }
}ContainedNode;
#pragma pack ()

ContainedNode * InitContainedNode(DBGraphHead * deBruijnGraphHead){
    ContainedNode * result = new ContainedNode[deBruijnGraphHead->nodeNumber];
    long int i = 0;
    for(i = 0; i<deBruijnGraphHead->nodeNumber; i++){
        result[i].length = strlen(deBruijnGraphHead->deBruijnGraph[i].contig);
    }
    return result;
}

bool ContainedNodeMatch(DBGraphHead * deBruijnGraphHead, ContainedNode * containedNode, long int pos){

    long int i = 0;
    long int cur = containedNode[pos].index;
    
    if(cur == 1){
        return true;
    }
    
    for(i = 0; i<deBruijnGraphHead->nodeNumber; i++){
        if( pos!=i && containedNode[pos].length == containedNode[i].length){

            if(ReverseComplementMatch(deBruijnGraphHead->deBruijnGraph[pos].contig,deBruijnGraphHead->deBruijnGraph[i].contig)){
                cur = containedNode[i].index;
                if(cur == 1){
                    return true;
                }
            }
            
        }
    }

    return false;

}

#pragma pack(2)
typedef struct SequenceNode{
    GraphNode * sequence;
    struct SequenceNode * next;
    SequenceNode(){
        sequence = NULL;
        next = NULL;
    }
}SequenceNode;
#pragma pack ()

long int RemoveExtendContigCycle(GraphNode * first, GraphNode * last, DBGraph * deBruijnGraph,long int kmerLength){
    long int i = 0;
    long int j = 0;
    GraphNode * one = NULL;
    GraphNode * two = NULL;
    GraphNode * temp = first;
    long int len = 0;
    long int len1 = 0;
    long int length = 0;
    while(temp->next!=NULL){
        if(temp->index==last->index){
            if(two == NULL && one == NULL){
                one = temp;
                len = j;
            }else if(two == NULL && one != NULL){
                two = temp;
                len1 = j;
            }else{
                one = two;
                two = temp;
                len = len1;
                len1 = j;
            }
        }
        temp = temp->next;
        j++;
    }
    if(two==NULL||one==NULL||(len1-len)!=(j-len1)){
        return 0;
    }
    while(i<=len1-len){
        if(one->index!=two->index){
            return 0;
        }
        length = length + strlen(deBruijnGraph[one->index].contig)-kmerLength +1;
        one = one->next;
        two = two->next;
        i++;
    }
    return length;
}

long int RemoveExtendContigCycleLeft(GraphNode * first, GraphNode * last, DBGraph * deBruijnGraph,long int kmerLength){
    long int i = 0;
    long int j = 1;
    GraphNode * one = NULL;
    GraphNode * two = NULL;
    GraphNode * temp = first->next;
    long int len = 0;
    long int len1 = 0;
    long int length = 0;
    while(temp!=NULL){
        if(temp->index==first->index){
            if(two == NULL && one == NULL){
                one = temp;
                len = j;
            }else if(two == NULL && one != NULL){
                two = temp;
                len1 = j;
                break;
            }
        }
        temp = temp->next;
        j++;
    }
    if(two==NULL||one==NULL||len!=(len1-len)){
        //cout<<len<<"--"<<len1<<endl;
        return 0;
    }
    temp = first;
    while(i<=len1-len){
        if(temp->index!=one->index){
            //cout<<"aa"<<endl;
            return 0;
        }
        length = length + strlen(deBruijnGraph[temp->index].contig)-kmerLength +1;
        one = one->next;
        temp = temp->next;
        i++;
    }
    return length;
}

long int SequenceNodeAdd(long int pos, GraphNode * first, GraphNode * last, DBGraph * deBruijnGraph,long int kmerLength, long int index){
    
    GraphNode * tempGraphNode = new GraphNode;
    tempGraphNode->index = pos;        
    long int cut = 0;
    if(index == 0){
        last->next = tempGraphNode;
        last = tempGraphNode;
        cut = RemoveExtendContigCycle(first,last,deBruijnGraph,kmerLength);
    }else{
        tempGraphNode->next = first;
        first = tempGraphNode;
        cut = RemoveExtendContigCycleLeft(first,last,deBruijnGraph,kmerLength);
    }
    return cut;
}

ofstream ocout;

double WeightRightCandidateContig(char * contigLeft, char * contigRight, long int repeatLength, ReadSet * readSet, long int readSetIndex,long int kmerLength, long int start, long int length, KmerSetHashTable * kmerSetHashTable, unsigned long int kmerSetHashTableCount){
    long int i = 0;
    long int j = 0;
    long int readLength = readSet[readSetIndex].readLength;
    long int insertSize = readSet[readSetIndex].insertSize;
    long int var = readSet[readSetIndex].var;
    
    long int nullNumber = 0;
    long int AllMatchNumber = 0;
    long int positionMatchNumber = 0;
    long int positionNonMatchNumber = 0;
    long int positionNonExist = 0;
    double outDegree = 0;
    double avgKmerNum = 0;
    
    long int min = 0;

    char * tempContigLeft;
    long int tempLen = strlen(contigLeft);
    if(tempLen>insertSize+var){
        min = insertSize+lambda*var;
        if(tempLen<min){
            min = tempLen;
        }
        tempContigLeft = new char[min+lambda*var+length-insertSize+1 + readLength];
        SubContig(tempContigLeft,contigLeft,tempLen-min,tempLen - insertSize + lambda*var + length + readLength);
    }else{
        return -1;
    }

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
    
    long int realAvgKmerNum = 0;
    char * temp = new char[readLength + 1];
    char * tempReverseComplement = new char[readLength + 1];
    char * tempKmer = new char[globalKmerLength + 1];
    char * tempKmerReverseComplement = new char[globalKmerLength + 1];
    
    long int continuousNE = 0;
    long int allExistReadNumber = 0;
    
    for(j = 0; j<tempLen - readLength; j++){
        SubContig(temp, contigLeft, tempLen - readLength - j, tempLen - j);
        long int a = SearchRead(temp, readSet + readSetIndex);
        ReverseComplement(temp, tempReverseComplement);
        long int b = SearchRead(tempReverseComplement, readSet + readSetIndex);
        if(a == -1 && b == -1){
            continuousNE++;
        }else{
            break;
        }
    }
    
    long int maxContinuousNE = continuousNE;
    
    long int outDegreeVar[1000];
    for(j=0;j<1000;j++){
        outDegreeVar[j] = 0;
    }
        
    for(j=start;j<start+length;j++){

        SubContig(temp, tempContigRight, j, j + readLength);

        ReadMate * tempReadMate = SearchLeftReadMate(temp, readSet+readSetIndex);
        long int a = SearchRead(temp, readSet + readSetIndex);
        ReadMate * temp1;
        
        if(tempReadMate==NULL && a ==-1){
            positionNonExist++;
            continuousNE++;
            if(maxContinuousNE<continuousNE){
                maxContinuousNE = continuousNE;
            }
            continue;
        }else{
            continuousNE = 0;
        }
        
        if(tempReadMate!=NULL){
            temp1 = tempReadMate;
            tempReadMate = tempReadMate->next;
            delete temp1;
        }   
        
        if(tempReadMate==NULL){
            nullNumber++;
            realAvgKmerNum = realAvgKmerNum + j;
            continue;
        }
        
        bool token = false;
        ReadMate * tempReadMate1 = tempReadMate;

        while(tempReadMate!=NULL){
            //if(j<readLength-repeatLength-1){
                allExistReadNumber ++;
            //}
            long int index1 = KMPIndexOfContigOfMisMatch(tempContigLeft,tempReadMate->readMate);
            long int tempDD = j + min - index1;
            if(index1!=-1 && tempDD<insertSize+lambda*var && tempDD>insertSize-lambda*var){
                outDegreeVar[AllMatchNumber] = tempDD;
                AllMatchNumber++;
                realAvgKmerNum = realAvgKmerNum + j;
                if(token!=true){
                    positionMatchNumber++;
                    token = true;
                }
                //outDegree = outDegree + labs(insertSize - tempDD);
                outDegree = outDegree + tempDD;
            }
            tempReadMate = tempReadMate->next;
        } 
        if(token==false){
            positionNonMatchNumber++;
        }
        
        while(tempReadMate1!=NULL){
            temp1 = tempReadMate1;
            delete [] tempReadMate1->readMate;
            tempReadMate1 = tempReadMate1->next;
            delete temp1;
        }                   
    }
    
    long int matchKmerNumber = 0;
    
    for(j=start+readLength-kmerLength;j<start+readLength-kmerLength+length;j++){
        
        SubContig(tempKmer, tempContigRight, j, j+globalKmerLength);
        ReverseComplement(tempKmer, tempKmerReverseComplement);

        long int hashP = KmerSetHashTableSearch(tempKmer,kmerSetHashTable,globalKmerLength,kmerSetHashTableCount);
        long int hashP1 = KmerSetHashTableSearch(tempKmerReverseComplement,kmerSetHashTable,globalKmerLength,kmerSetHashTableCount);

        if(hashP!=0){
            avgKmerNum = avgKmerNum + 
                kmerSetHashTable[hashP-1].frequency[readSetIndex];
            matchKmerNumber++;
        } 
        if(hashP1!=0){
            avgKmerNum = avgKmerNum + 
                kmerSetHashTable[hashP1-1].frequency[readSetIndex];
            matchKmerNumber++;
        }       
    }
    avgKmerNum = avgKmerNum/matchKmerNumber;

    
    delete [] temp;
    delete [] tempReverseComplement;
    delete [] tempKmer;
    delete [] tempKmerReverseComplement;
    delete [] tempContigLeft;
    delete [] tempContigRight;  

    
    double readP = (readLength-repeatLength-1)*readSet[readSetIndex].avgReadNumber;
    double readPD = 0;
    if(readP>0){
        readPD = sqrt(readP);
    }
    
    if(positionMatchNumber<2){
        return 0;
    }
    
    double a = (double)(positionMatchNumber)/(double)(length - nullNumber - positionNonExist);
    double b = fabs((outDegree/(double)(AllMatchNumber)) - (double)insertSize)/(2*(double)var);
    double c = avgKmerNum/readSet[readSetIndex].avgKmerNumber;

    double result = 0;
    double var1 = (double)positionNonExist-readSet[readSetIndex].gapProblity*(double)length;
    double d = sqrt((double)length*readSet[readSetIndex].gapProblity*(1-readSet[readSetIndex].gapProblity))*(lambda + 1);
    double a1 = a;
    double a2 = (double)(positionMatchNumber)/(double)(length);
    
    double r = ComputeSD(outDegreeVar,insertSize,var);
    
    if(var1 - d > 0){
        result = a*pow(r,0.5) - (1-r)*c*(var1-d)/d;        
        
    }else{       
        result = a*pow(r,0.5);       
    }

    ocout.close();
    
    if(result<0){
        return 0;
    }
    
    return result;
    
}


double WeightLeftCandidateContig(char * contigLeft, char * contigRight, long int repeatLength, ReadSet * readSet, long int readSetIndex,long int kmerLength, long int start, long int length, KmerSetHashTable * kmerSetHashTable, unsigned long int kmerSetHashTableCount){
    long int i = 0;
    long int j = 0;
    long int readLength = readSet[readSetIndex].readLength;
    long int insertSize = readSet[readSetIndex].insertSize;
    long int var = readSet[readSetIndex].var;
    
    long int nullNumber = 0;
    long int AllMatchNumber = 0;
    long int positionMatchNumber = 0;
    long int positionNonExist = 0;
    long int positionNonMatchNumber = 0;
    double outDegree = 0;
    double avgKmerNum = 0;
    
    long int min = 0;

    char * tempContigRight;
    long int tempLen = strlen(contigRight);
    if(tempLen>insertSize+var){
        min = insertSize+lambda*var;
        if(tempLen<min){
            min = tempLen;
        }
        
        tempContigRight = new char[min-insertSize+lambda*var+length + 1 + readLength];
        if(insertSize - lambda*var - length - readLength<0){
            SubContig(tempContigRight, contigRight,0,min);
        }else{
            SubContig(tempContigRight, contigRight,insertSize - lambda*var - length - readLength,min);
        }
        
    }else{
        return -1;
    }
    
    long int tempLen1 = strlen(contigLeft);
    
    char * tempContigLeft = new char[readLength + length + 1];
    for(i = 0; i<readLength+length; i++){
        if(i<length){
            tempContigLeft[i] = contigLeft[tempLen1-kmerLength + 1 - length + i];
        }else{
            tempContigLeft[i] = contigRight[j];
            j++;
        }
    }
    tempContigLeft[i] = '\0';
    
    long int realAvgKmerNum = 0;
    char * temp = new char[readLength + 1];
    char * tempReverseComplement = new char[readLength + 1];
    char * tempKmer = new char[globalKmerLength + 1];
    char * tempKmerReverseComplement = new char[globalKmerLength + 1];
    
    long int continuousNE = 0;
    long int allExistReadNumber = 0;
    
    
    for(j = 0; j<tempLen - readLength; j++){
        SubContig(temp, contigRight, j, j + readLength);
        long int a = SearchRead(temp, readSet + readSetIndex);
        ReverseComplement(temp, tempReverseComplement);
        long int b = SearchRead(tempReverseComplement, readSet + readSetIndex);
        if(a == -1 && b == -1){
            continuousNE++;
        }else{
            break;
        }
    }
    
    long int maxContinuousNE = continuousNE;
    
    long int outDegreeVar[1000];
    for(j=0;j<1000;j++){
        outDegreeVar[j] = 0;
    }
        
    for(j=start;j<start+length;j++){
        
        SubContig(temp, tempContigLeft, j, j+readLength);
        ReverseComplement(temp, tempReverseComplement);
        ReadMate * tempReadMate = SearchRightReadMate(temp, readSet+readSetIndex);
        long int a = SearchRead(tempReverseComplement, readSet + readSetIndex);  
        ReadMate * temp1;
        if(tempReadMate==NULL && a==-1){
            positionNonExist++;
            continuousNE++;
            if(maxContinuousNE<continuousNE){
                maxContinuousNE = continuousNE;
            }
            continue;
        }else{
            continuousNE = 0;
        }
        
        if(tempReadMate!=NULL){
            temp1 = tempReadMate;
            tempReadMate = tempReadMate->next;
            delete temp1;
        }
             
        if(tempReadMate==NULL){
            realAvgKmerNum = realAvgKmerNum + j;
            nullNumber++;
            continue;
        }      
        ReadMate * tempReadMate1 = tempReadMate;
        
        bool token = false;
        while(tempReadMate!=NULL){
           
           //if(j>length + repeatLength - readLength){
               allExistReadNumber++;
           //}
           
           long int index1 = KMPIndexOfContigOfMisMatch(tempContigRight,tempReadMate->readMate);
           long int tempDD = index1 + insertSize - lambda*var - j;
           if(insertSize - lambda*var - length - readLength<0){
               tempDD = index1 + readLength + length -j;
           }
           if(index1!=-1 && tempDD<insertSize+lambda*var && tempDD>insertSize-lambda*var){
               outDegreeVar[AllMatchNumber] = tempDD;
               AllMatchNumber++;
               realAvgKmerNum = realAvgKmerNum + length - j;
               if(token!=true){
                   positionMatchNumber++;
                   token = true;
               }
               
               outDegree = outDegree + tempDD;
               
           }
           tempReadMate = tempReadMate->next;
        }
        if(token == false){
            positionNonMatchNumber++;
        }
       
        while(tempReadMate1!=NULL){
             temp1 = tempReadMate1;
             delete [] tempReadMate1->readMate;
             tempReadMate1 = tempReadMate1->next;
             delete temp1;
        }                     
    }
    
    long int matchKmerNumber = 0;
    
    for(j=start;j<start+length;j++){
        SubContig(tempKmer, tempContigLeft, j, j+globalKmerLength);
        ReverseComplement(tempKmer, tempKmerReverseComplement);
        long int hashP = KmerSetHashTableSearch(tempKmer,kmerSetHashTable,globalKmerLength,kmerSetHashTableCount);
        long int hashP1 = KmerSetHashTableSearch(tempKmerReverseComplement,kmerSetHashTable,globalKmerLength,kmerSetHashTableCount);
        if(hashP!=0){
            avgKmerNum = avgKmerNum + 
                kmerSetHashTable[hashP-1].frequency[readSetIndex];
            matchKmerNumber++;
        }
        if(hashP1!=0){
            avgKmerNum = avgKmerNum + 
                kmerSetHashTable[hashP1-1].frequency[readSetIndex];
            matchKmerNumber++;
        }
    }
    avgKmerNum = avgKmerNum/matchKmerNumber;
     
    delete [] temp;
    delete [] tempReverseComplement;
    delete [] tempKmer;
    delete [] tempKmerReverseComplement;
    delete [] tempContigLeft;
    delete [] tempContigRight;  
    
    
    
    
    
    double readP = (readLength-repeatLength-1)*readSet[readSetIndex].avgReadNumber;
    double readPD = 0;
    if(readP>0){
        readPD = sqrt(readP);
    }
    
    if(positionMatchNumber<2){
        return 0;
    }
    
    double a = (double)(positionMatchNumber)/(double)(length - nullNumber - positionNonExist); 
    double b = fabs((outDegree/(double)(AllMatchNumber)) - (double)insertSize)/((double)var*2);
    double c = avgKmerNum/readSet[readSetIndex].avgKmerNumber;

    
    
    double a1 = a;
    double a2 = (double)(positionMatchNumber)/(double)(length);
    
    double result = 0;
    double d = sqrt((double)length*readSet[readSetIndex].gapProblity*(1-readSet[readSetIndex].gapProblity))*(lambda + 1);
    double var1 = (double)positionNonExist-readSet[readSetIndex].gapProblity*(double)length;
    
    double r = ComputeSD(outDegreeVar,insertSize,var);
    
    if(var1 - d > 0){
        result = a*pow(r,0.5) - (1-r)*c*(var1-d)/d;        
    }else{       
        result = a*pow(r,0.5);
    }
    
    if(result<0){
        return 0;
    }

    return result;
    
}

pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;

pthread_mutex_t mutex2 = PTHREAD_MUTEX_INITIALIZER;

char * ExtendContigRight(DBGraph * deBruijnGraph, long int pos, ReadSet * readSet, int setNumber, KmerSetHashTable * kmerSetHashTable, unsigned long int kmerSetHashTableCount, int kmerLength, ContainedNode * containedNode){
    
    long int i = 0;
    long int j = 0;
    long int outCount = 0;
    long int dfsLength = 0;
    long int contigLength = 0;
    long int storeLength = 0;
    long int posLength = 0;
    
    
    
    
    contigLength = strlen(deBruijnGraph[pos].contig);
    storeLength = contigLength + storeContigExtendLength;
    char * contig = new char[storeLength + 1];
    strcpy(contig, deBruijnGraph[pos].contig);
    
    outCount = NodeCount(deBruijnGraph[pos].outNode);  
    dfsLength = dfsExtendLength + kmerLength;
    
    GraphNode * first = new GraphNode;
    first->index = pos;
    GraphNode * last = first;
    
    while(outCount!=0){
        if(outCount == 0){
            DeleteGraphNode(first);
            return contig;
            break;
        }
        if(outCount == 1){
            pos = deBruijnGraph[pos].outNode->index;
            posLength = strlen(deBruijnGraph[pos].contig);
            contig[contigLength - kmerLength + 1] = '\0';
            contigLength = contigLength + posLength - kmerLength + 1;
            while(contigLength > storeLength){   
                storeLength = contigLength + storeContigExtendLength;
                contig = (char *)realloc(contig, storeLength + 1);              
            }
            strcat(contig, deBruijnGraph[pos].contig);   
            containedNode[pos].index = 1;
            
            outCount = NodeCount(deBruijnGraph[pos].outNode);
            
            GraphNode * tempGraphNode = new GraphNode;
            tempGraphNode->index = pos;
            last->next = tempGraphNode;
            last = tempGraphNode;
            
            long int cycleCut = RemoveExtendContigCycle(first,last,deBruijnGraph,kmerLength);
            if(cycleCut!=0){
                char * tempContig = new char[strlen(contig)-cycleCut+1];
                SubContig(tempContig,contig,0,strlen(contig)-cycleCut);
                delete []contig;
                contig = NULL;
                DeleteGraphNode(first);
                return tempContig;
            } 
                    
        }
        if(outCount>1){
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
            dfsNode = DFSTranverse(deBruijnGraph,pos,kmerLength,dfsLength);        
            DFSNode * tempDFSNode = dfsNode; 
            DFSNode * previousDFSNode = dfsNode; 
            while(dfsNode!=NULL){
                singleMatchNumber = 0;
                if(strlen(dfsNode->contig)<=dfsLength){            
                    if(dfsNode->adjNode!=tempPrevious){
                        if(tempPrevious!=-1){
                            i++;
                        }
                        tempAdjNode[i] =dfsNode->adjNode;
                        tempPrevious = dfsNode->adjNode;                  
                    }                                                    
                    dfsNode = dfsNode->next;
                    continue;
                }else{
                    double weight = 0;
                    long int weightNum = 0;
                    for(int p = 0;p<setNumber;p++){
                        
                        double tempDD = WeightRightCandidateContig(contig,dfsNode->contig,strlen(deBruijnGraph[pos].contig),readSet,p,kmerLength,0,dfsLength-kmerLength,kmerSetHashTable,kmerSetHashTableCount);
                        
                    
                        if(tempDD>=0){
                            weight = weight + tempDD;
                            weightNum++;
                        }
                        if(tempDD>0){
                            singleMatchNumber++;
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
            
            
            if((matchNumber[0] > matchNumber[1]) || ((tempWeight[0] - tempWeight[1]>0.2) && !(tempWeight[0]>0.7 && tempWeight[1]>0.7))){
                
                pos = tempAdjNode[0];
                
                posLength = strlen(deBruijnGraph[pos].contig);
                contig[contigLength - kmerLength + 1] = '\0';
                contigLength = contigLength + posLength - kmerLength + 1;
                while(contigLength > storeLength){   
                    storeLength = contigLength + storeContigExtendLength;
                    contig = (char *)realloc(contig, storeLength + 1);              
                }
                strcat(contig, deBruijnGraph[pos].contig);
                
                containedNode[pos].index = 1;
                
                outCount = NodeCount(deBruijnGraph[pos].outNode); 
                
                
                GraphNode * tempGraphNode = new GraphNode;
                tempGraphNode->index = pos;
                last->next = tempGraphNode;
                last = tempGraphNode;   
                long int cycleCut = RemoveExtendContigCycle(first,last,deBruijnGraph,kmerLength);     
                
                if(cycleCut!=0){
                    char * tempContig = new char[strlen(contig)-cycleCut+1];
                    SubContig(tempContig,contig,0,strlen(contig)-cycleCut);
                    delete []contig;
                    contig = NULL;
                    DeleteDFSNode(tempDFSNode);
                    DeleteGraphNode(first);
                    return tempContig;
                }            
                
            }else{
            	DeleteDFSNode(tempDFSNode);
                DeleteGraphNode(first);
                return contig;
                break;
            }
            DeleteDFSNode(tempDFSNode);
            
        }
    }
    DeleteGraphNode(first);
    return contig;
    
    
}

char * ExtendContigLeft(char * contig, DBGraph * deBruijnGraph, long int pos, ReadSet * readSet, int setNumber, KmerSetHashTable * kmerSetHashTable, unsigned long int kmerSetHashTableCount, int kmerLength, ContainedNode * containedNode){
    
    long int i = 0;
    long int j = 0;
    long int inCount = 0;

    long int dfsLength = 0;
    
    inCount = NodeCount(deBruijnGraph[pos].inNode);
    long int readLength = 0;
    for(i = 0;i<setNumber;i++){
        if(readLength<readSet[i].readLength){
            readLength = readSet[i].readLength;
        }
    }
    
    dfsLength = dfsExtendLength + kmerLength;
    
    GraphNode * first = new GraphNode;
    first->index = pos;
    GraphNode * last = first;
    
    
    
    while(inCount!=0){
        if(inCount==0){
            DeleteGraphNode(first);
            return contig;
            break;
        }
        if(inCount == 1){
            
            char * tempContig = contig;
            contig = new char[strlen(deBruijnGraph[deBruijnGraph[pos].inNode->index].contig)+strlen(contig)-kmerLength+2]; 
            AppendRight(contig,deBruijnGraph[deBruijnGraph[pos].inNode->index].contig,tempContig,kmerLength);
            delete []tempContig;
            pos = deBruijnGraph[pos].inNode->index;         
            
            containedNode[pos].index = 1;
            inCount = NodeCount(deBruijnGraph[pos].inNode); 
            
            GraphNode * tempGraphNode = new GraphNode;
            tempGraphNode->index = pos;
            tempGraphNode->next = first;
            first = tempGraphNode;
            
            long int cycleCut = RemoveExtendContigCycleLeft(first,last,deBruijnGraph,kmerLength);
            if(cycleCut!=0){
                tempContig = new char[strlen(contig)-cycleCut+1];
                SubContig(tempContig,contig,cycleCut,strlen(contig));
                delete []contig;
                contig = NULL;
                DeleteGraphNode(first);
                return tempContig;
            }
              
            
        }
        if(inCount>1){
            int singleMatchNumber = 0;
            int matchNumber[inCount];
            double tempWeight[inCount];
            for(i=0;i<inCount;i++){
                tempWeight[i] = 0;
                matchNumber[i] = 0;
            }
            long int tempAdjNode[inCount];
            for(i=0;i<inCount;i++){
                tempAdjNode[i] = -1;
            }
            
            long int tempPrevious = -1;
            i = 0;
            
            DFSNode * dfsNode = DFSTranverseLeft(deBruijnGraph,pos,kmerLength,dfsLength);

            DFSNode * tempDFSNode = dfsNode; 
            DFSNode * previousDFSNode = dfsNode; 

            while(dfsNode!=NULL){
                singleMatchNumber = 0;
                if(strlen(dfsNode->contig)<=dfsLength){
                    
                    if(dfsNode->adjNode!=tempPrevious){
                        if(tempPrevious!=-1){
                            i++;
                        }
                        tempAdjNode[i] =dfsNode->adjNode;
                        tempPrevious = dfsNode->adjNode;                  
                    }
                    
                    dfsNode = dfsNode->next;                    
                    continue;
                }else{
                    double weight = 0;
                    long int weightNum = 0;
                    for(int p = 0;p<setNumber;p++){
                        
                        double tempDD = WeightLeftCandidateContig(dfsNode->contig,contig,strlen(deBruijnGraph[pos].contig),readSet,p,kmerLength,0,dfsLength-kmerLength,kmerSetHashTable,kmerSetHashTableCount);
                        
                        if(tempDD>=0){
                            weight = weight + tempDD;
                            weightNum++;
                        }
                        if(tempDD>0){
                            singleMatchNumber++;
                        }                       
                    }
                    if(weightNum!=0){
                        weight = weight/(double)weightNum;
                    }
                    
                    if(dfsNode->adjNode!=tempPrevious){
                        if(tempPrevious!=-1){
                            i++;
                        }
                        tempWeight[i]=weight; 
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
            
            double maxWeight = 0;
            long int maxPosition = -1;
            
            for(i=0;i<inCount-1;i++){
                for(j=i;j<inCount;j++){
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
            
            
            if((matchNumber[0] > matchNumber[1]) || ((tempWeight[0] - tempWeight[1]>0.2) && !(tempWeight[0]>0.7 && tempWeight[1]>0.7))){
           
                pos = tempAdjNode[0];
                
                
                char * tempContig = contig;
                contig = new char[strlen(deBruijnGraph[pos].contig)+strlen(contig)-kmerLength+2]; 
                AppendRight(contig,deBruijnGraph[pos].contig,tempContig,kmerLength);
                delete []tempContig;
                containedNode[pos].index = 1;
                
                inCount = NodeCount(deBruijnGraph[pos].inNode);
                
                GraphNode * tempGraphNode = new GraphNode;
                tempGraphNode->index = pos;
                tempGraphNode->next = first;
                first = tempGraphNode;
            
                long int cycleCut = RemoveExtendContigCycleLeft(first,last,deBruijnGraph,kmerLength);
                if(cycleCut!=0){
                    tempContig = new char[strlen(contig)-cycleCut+1];
                    SubContig(tempContig,contig,cycleCut,strlen(contig));
                    delete []contig;
                    contig = NULL;
                    DeleteDFSNode(tempDFSNode);
                    DeleteGraphNode(first);
                    return tempContig;
                }   
                
            }else{
                DeleteDFSNode(tempDFSNode);
                DeleteGraphNode(first);
                return contig;
                break;
            }
            
        }
    }
    DeleteGraphNode(first);
    return contig;
    
    
}

void WriteContigSet(ContigSet * contigSet, char * str){
    long int i = 0;
    long int j = 0;
    long int p = 0;
    char *tempContig; 
    ofstream ocout;
    ocout.open(str);
    long int rowLength = 2000;
    while(contigSet!=NULL){
        long int len = strlen(contigSet->contig);
        char * longContig = new char[len-49];
        SubContig(longContig,contigSet->contig,25,len-25);
        len = len - 50;
        j = (long int)(len/rowLength);
        for(i = 0; i<j+1;i++){
            if(i!=j){
                tempContig = new char[rowLength+1];
                SubContig(tempContig,longContig,i*rowLength,i*rowLength+rowLength);
            }else{
                tempContig = new char[len - i*rowLength+1];
                SubContig(tempContig,longContig,i*rowLength,len);
            }
            ocout<<">"<<p<<"--"<<i<<endl;
            ocout<<tempContig<<endl;
            delete []tempContig;
        }
        p++;
        contigSet = contigSet->next;
    }
}

void WriteContigSetLong(ContigSet * contigSet, char * str){
    long int p = 0;
    ofstream ocout;
    ocout.open(str);
    while(contigSet!=NULL){
        long int len = strlen(contigSet->contig);
        ocout<<">"<<p<<":"<<len<<endl;
        ocout<<contigSet->contig<<endl;
        p++;
        contigSet = contigSet->next;
    }
}

ContigSet * GetContigSetLong(char * str){
    long int i = 0;
    
    ifstream icin;
    icin.open(str);
    
    char * tempContig = new char[5000000];
    ContigSet * head = new ContigSet;
    ContigSet * contigSet = head;
    
    while(icin.getline(tempContig, 5000000)){
        if(i%2!=0){
            ContigSet * tempContigSet = new ContigSet;
            tempContigSet->contig = new char[strlen(tempContig) + 1];
            CopyContig(tempContigSet->contig, tempContig);
            contigSet->next = tempContigSet;
            contigSet = tempContigSet;
            
        }
        i++;
    }
    contigSet = head->next;
    head->next = NULL;
    delete head;
    return contigSet;
    
}


#pragma pack(2)
typedef struct ConstructContigSetP{
      DBGraphHead * deBruijnGraphHead;
      ReadSet * readSet;
      ContigSet * contigSetHead;
      int setNumber;
      KmerSetHashTable * kmerSetHashTable;
      unsigned long int kmerSetHashTableCount;
      int kmerLength;
      ContainedNode * containedNode;
      int threadIndex;
      int totalThreadNumber;
}ConstructContigSetP;
#pragma pack ()

void * ConstructContigSetThread(void * arg){

    ConstructContigSetP * constructContigSetpP = (ConstructContigSetP *)arg;
    
    unsigned long int i = constructContigSetpP->threadIndex;
    unsigned long int j = 0;
    
    DBGraph * deBruijnGraph = constructContigSetpP->deBruijnGraphHead->deBruijnGraph;
    long int graphNodeNumber = constructContigSetpP->deBruijnGraphHead->nodeNumber;
    ContainedNode * containedNode = constructContigSetpP->containedNode;
    
    
        
    while(i<graphNodeNumber){
        
        if(containedNode[i].length<extendCutOff||ContainedNodeMatch(constructContigSetpP->deBruijnGraphHead, containedNode, i)){
            i = i + constructContigSetpP->totalThreadNumber;
            continue;
        }        
        
        containedNode[i].index = 1;
        
        char * contigRight = ExtendContigRight(deBruijnGraph, i, constructContigSetpP->readSet, constructContigSetpP->setNumber, constructContigSetpP->kmerSetHashTable, constructContigSetpP->kmerSetHashTableCount, constructContigSetpP->kmerLength, constructContigSetpP->containedNode);
        
        char * tempContig = ExtendContigLeft(contigRight, deBruijnGraph, i, constructContigSetpP->readSet, constructContigSetpP->setNumber, constructContigSetpP->kmerSetHashTable, constructContigSetpP->kmerSetHashTableCount, constructContigSetpP->kmerLength, constructContigSetpP->containedNode);
        
        
        
        
        ContigSet * tempContigSet = new ContigSet;
        tempContigSet->contig = tempContig;
        pthread_mutex_lock(&mutex2);     
        tempContigSet->next = constructContigSetpP->contigSetHead->next;
        constructContigSetpP->contigSetHead->next = tempContigSet;
        pthread_mutex_unlock(&mutex2);
        
        
                
        i = i + constructContigSetpP->totalThreadNumber;
    }

}


ContigSet * ConstructContigSet(DBGraphHead * deBruijnGraphHead, ReadSet * readSet, ContigSet * contigSetHead, int setNumber, KmerSetHashTable * kmerSetHashTable, unsigned long int kmerSetHashTableCount, int kmerLength, int totalThreadNumber, int extnedCutOff){
    
    pthread_t tid[totalThreadNumber];
    int i = 0;
    
    ConstructContigSetP * constructContigSetP = new ConstructContigSetP[totalThreadNumber];
    ContainedNode * containedNode = InitContainedNode(deBruijnGraphHead);
    for(i = 0; i<totalThreadNumber; i++){
        constructContigSetP[i].deBruijnGraphHead = deBruijnGraphHead;
        constructContigSetP[i].readSet = readSet;
        constructContigSetP[i].contigSetHead = contigSetHead;
        constructContigSetP[i].setNumber = setNumber;
        constructContigSetP[i].kmerSetHashTable = kmerSetHashTable;
        constructContigSetP[i].kmerSetHashTableCount = kmerSetHashTableCount;
        constructContigSetP[i].kmerLength = kmerLength;
        constructContigSetP[i].containedNode = containedNode;
        constructContigSetP[i].threadIndex = i;
        constructContigSetP[i].totalThreadNumber = totalThreadNumber;
        
        if(pthread_create(&tid[i], NULL, ConstructContigSetThread, (void *)&constructContigSetP[i])){
            cout<<"create thread wrong!"<<endl;
            exit(1);
        }
        
    }
    for(i = 0; i < totalThreadNumber; i++){
        pthread_join(tid[i], NULL);
    }
    
      
}

int CompareContigSet(const void * a, const void * b){
    Contig * ha = (Contig *)a;
    Contig * hb = (Contig *)b;

    return (strlen(hb->contig) - strlen(ha->contig));
}

ContigSet * GetContigSet(char * str){
    unsigned long long int i = 0;
    long int j = 0;
    unsigned long long int count = 0;
    char *temp = new char[10000000];

	ifstream icin1;
	icin1.open(str);
	
	while(icin1.getline(temp,10000000)){
		count++;
	}
	count = count/2;
	Contig * contigSet = new Contig[count];
	icin1.close();
	ifstream icin;
	icin.open(str);
	while(icin.getline(temp,10000000)){
		if(i%2==0){
            i++;
            continue;
        }
        contigSet[i/2].contig = new char[strlen(temp)+1];
        CopyContig(contigSet[i/2].contig,temp);
        i++;
	}
	qsort(contigSet, count, sizeof(Contig), CompareContigSet);
	icin.close();
	delete []temp;
	
	ContigSet * result = new ContigSet;
	ContigSet * first = result;
	for(i = 0; i<count; i++){
        
        result->contig = contigSet[i].contig;
        contigSet[i].contig = NULL;
        
        if(i!=count-1){
            result->next = new ContigSet;
            result = result->next;
        }
        
    }
    delete [] contigSet;
	
	return first;
}

double nihe(int * gap, int length, double p){
    
    int i = 0;
    int t = 0;
    double s = 1;
    double temp = 0;
    long int all = 0;
    
    for(i = 0; i<length; i++){
        all = all + gap[i];
    }
    
    int j = 0;
    
    for(i = 0; i<=length; i++){
        s = 1;
        j = i;
        for(t=length;t>length-i;t--){
            s = s*(double)t/(double)j;
            j--;
        }
        //cout<<"------------"<<s<<"----"<<all<<endl;
        double temp1 = (double)all*s*pow(p,i)*pow((1-p),length-i);
        double temp2 = temp1  - (double)gap[i];
        temp = temp + (pow(temp2,2)/temp1);
        //cout<<"--------------------"<<temp1<<"--"<<temp2<<endl;
    }
    
    return temp;
    
    
}

int GetMateProblity(char * contig, int length, ReadSet * readSet, int readSetIndex){
    
    long int i = 0;
    long int j = 0;
    long int readLength = readSet[readSetIndex].readLength;
    long int insertSize = readSet[readSetIndex].insertSize;
    long int var = readSet[readSetIndex].var;
    int index1 = 0;
    
    long int s = 0;
    
    long int result[101];
    for(i = 0; i<101; i++){
        result[i] = 0;
    }
    i = 0;
    
    long int nullNumber = 0;
    long int positionMatchNumber = 0;
    long int positionNonExist = 0;
    long int allMatchNumber;
    
    long int outDegree = 0;
    
    long int contigLength = strlen(contig);
    char temp[readLength + 1];
    
    ReadMate * temp1;
    
    char * distanceAddress = new char[10];
    strcpy(distanceAddress, "b.fa");
    ofstream ocout;
    ocout.open(distanceAddress);
    
    
    for(j = insertSize + lambda*var; j<contigLength - 2*readLength - length; j++){
        //cout<<"--"<<j<<endl;
        index1 = 0;
        SubContig(temp, contig, j, j + readLength);

        ReadMate * tempReadMate = SearchLeftReadMate(temp, readSet+readSetIndex);
        long int a = SearchRead(temp, readSet + readSetIndex);
        
        if(tempReadMate==NULL && a ==-1){
            positionNonExist++;
            index1 = 1;
        }
        
        if(tempReadMate!=NULL){
            temp1 = tempReadMate;
            tempReadMate = tempReadMate->next;
            delete temp1;
        }   
        
        if(tempReadMate==NULL&&index1!=1){
            nullNumber++;
        }
        
        bool token = false;

        while(tempReadMate!=NULL){
            index1 = KMPIndexOfContigOfMisMatch(contig,tempReadMate->readMate, j-insertSize-lambda*var);
            long int tempDD = j - index1 + readLength;
            if(index1!=-1 && tempDD<insertSize+lambda*var && tempDD>insertSize-lambda*var){
                allMatchNumber ++ ;
                if(token!=true){
                    positionMatchNumber++;
                    token = true;
                }
                outDegree = outDegree + tempDD;
            }
            temp1 = tempReadMate;
            tempReadMate = tempReadMate->next;
            delete [] temp1->readMate;
            delete temp1;       
        } 
        
        i++;
        if(i == length){
            double a = 0;
            if(length - nullNumber - positionNonExist!=0){
                a = (double)(positionMatchNumber)/(double)(length - nullNumber - positionNonExist);
            }else{
                s ++;
            }
            double b = fabs(((double)outDegree/(double)allMatchNumber) - insertSize);
            ocout<<"a:"<<a<<"--pMa:"<<positionMatchNumber<<"--pNo:"<<positionNonExist<<"--nu:"<<nullNumber<<"--b:"<<b/(double)var<<endl;
            result[(int)(a*100)]++;
            //cout<<"re::"<<(int)(a*100)<<endl;
            nullNumber = 0;
            positionMatchNumber = 0;
            positionNonExist = 0;
            allMatchNumber = 0;
            outDegree = 0;
            i = 0; 
            //s ++;
            /*
            if(s>2000){
                break;
            }
            */
        }
    }
    
    char address[30];
    strcpy(address, "result.fa");
    ofstream ocout1;
    ocout1.open(address, ios::app);
    
    for(i = 0; i<101; i++){
        ocout1<<result[i]<<endl;;
    }
    ocout1<<"------------------"<<s<<endl;
    ocout1.close();
    
}

int GetGapProblityChromosome(ContigSet * contigSet, int length, ReadSet * readSet, int readSetIndex){
    long int i = 0;
    long int j = 0;
    long int max = 0;
    long int all = 0;
    long int all1 = 0;
    long int allLength = 0;
    long int allLength1 = 0;
    long int tempAll = 0;
    long int readLength = readSet[readSetIndex].readLength;
    char read[readLength+1];
    char readRC[readLength+1];
    long int gapNumber[length+1];
    long int gapNumber1[length+1];
    
    for(i=0;i<=length;i++){
        gapNumber[i] = 0;
        gapNumber1[i] = 0;
    }

    char * testSubContig = new char[readLength+length];
    
    long int s = 0;
    
    char address[30];
    strcpy(address, "result.fa");
    ofstream ocout;
    ocout.open(address, ios::app);
    
    int error[101];
    int errorPossion[101];
    for(i=0;i<101;i++){
        error[i]=0;
        errorPossion[i]=0;
    }
    while(contigSet!=NULL){
        
        //ocout<<"--"<<s<<endl;
        s++;
        
        long int contigLength = strlen(contigSet->contig);
        j = 0;
        /*
        if(contigLength<100*length){
            break;
        }
        */

        int * gapIndex = new int[contigLength-readLength];                   
        for(i = 0; i<contigLength - readLength; i++){
            SubContig(read,contigSet->contig,i,i+readLength);
            ReverseComplement(read, readRC);
            long int a = SearchRead(read, readSet+readSetIndex);
            //cout<<"ss:"<<a<<endl;
            long int b = SearchRead(readRC, readSet+readSetIndex);
            //cout<<"aa:"<<a<<":"<<b<<endl;
            if(a==-1&&b==-1){
                all++;
                j++;
                gapIndex[i] = 1;
            }else{
                all1++;
                gapIndex[i] = 0;
            }
            if(i>=length-1){
                gapNumber[j]++;
                j = j - gapIndex[i-length+1];
            }
            //cout<<"bb"<<endl;
        }
        j = 0;
        char * region = new char[readLength + length + 1];
        for(i = 0; i<contigLength - readLength; i++){
            SubContig(region,contigSet->contig,i,i+readLength + length);
            char tempChar = region[readLength-1];
            while(1){
                int k = rand();
                if(k%4==0&&tempChar!='A'){
                    region[readLength-1] = 'A';
                    break;
                }
                if(k%4==1&&tempChar!='T'){
                    region[readLength-1] = 'T';
                    break;
                }
                if(k%4==2&&tempChar!='G'){
                    region[readLength-1] = 'G';
                    break;
                }
                if(k%4==3&&tempChar!='C'){
                    region[readLength-1] = 'C';
                    break;
                }
            }
            j = 0;
            for(int t = 0; t<length; t++){
                
                SubContig(read, region, t, t+readLength);
                ReverseComplement(read, readRC);
                long int a = SearchRead(read, readSet+readSetIndex);
                long int b = SearchRead(readRC, readSet+readSetIndex);         
                if(a==-1&&b==-1){
                    j++;
                }
                
            }
            gapNumber1[j]++;           
            
            //cout<<"bb"<<endl;
        }
        
        
        allLength = contigLength - readLength;
        delete [] gapIndex;  
        //cout<<all<<"--"<<all1<<"--"<<allLength<<endl;
        double gap = (double)all/(double)allLength;
        double avg = length*gap;
        double sd = sqrt(length*gap*(1-gap));
        
        j = 0;
        long int j1 = 0;
        long int possionNumber = 0;
        for(i=0;i<=length;i++){
            //ocout<<gapNumber[i]<<"--";
            //possionNumber = possionNumber + i*gapNumber[i];
            if(i>avg+(lambda+1)*sd||i<avg-(lambda+1)*sd){
                j = j+gapNumber[i];
                j1 = j1 + gapNumber1[i];
            }
        }
        ocout<<endl;
        
        //double gapPossion = (double)possionNumber/(double)(allLength-length);
        //double sdPossion = sqrt(gapPossion);
        /*
        for(i=0;i<=length;i++){
            if(i>gapPossion+lambda*sdPossion||i<gapPossion-lambda*sdPossion){
                j1 = j1+gapNumber[i];
            }
        }
        */
        //ocout<<all<<"----"<<allLength<<endl;
        ocout<<s<<"、length:"<<allLength<<endl;
        ocout<<"--gapProblity:"<<gap<<"--avg:"<<avg<<"--sd:"<<sd<<"--error:"<<(double)j/(double)(allLength)<<"--"<<(double)j1/(double)(allLength)<<endl;;
        //ocout<<"--possionPro:"<<gapPossion<<"--sdPossion:"<<sdPossion<<"--error:"<<(double)j1/(double)(allLength-length)<<endl;   
        //error[(int)((double)j/(double)(allLength)*100)] ++;
        //errorPossion[(int)((double)j1/(double)(allLength -length)*100)]++;
        
        gap = 0;
        avg = 0;
        sd = 0;
        allLength = 0;
        all = 0;
        all1 = 0;
        
        for(i=0;i<=length;i++){
            gapNumber[i] = 0;
            gapNumber1[i] = 0;
        }
           
        contigSet = contigSet->next;        
    }
    /*
    for(i=0;i<101;i++){
        ocout<<error[i]<<"--";
    }
    ocout<<endl;
    
    for(i=0;i<101;i++){
        ocout<<errorPossion[i]<<"--";
    }
    ocout<<endl;
    */
    ocout.close();
    return 1;
}

int GetGapProblity(ContigSet * contigSet, int length, ReadSet * readSet, int readSetIndex){
    long int i = 0;
    long int j = 0;
    long int max = 0;
    long int all = 0;
    long int all1 = 0;
    long int allLength = 0;
    long int allLength1 = 0;
    long int tempAll = 0;
    long int readLength = readSet[readSetIndex].readLength;
    char read[readLength+1];
    char readRC[readLength+1];
    long int gapNumber[length+1];
    long int gapNumber1[length+1];
    
    for(i=0;i<=length;i++){
        gapNumber[i] = 0;
        gapNumber1[i] = 0;
    }

    char * testSubContig = new char[readLength+length];
    
    long int s = 0;
    
    long int percentageOfCorrect = 0;
    long int percentageOfError = 0;
    
    char address[30];
    strcpy(address, "result.fa");
    ofstream ocout;
    ocout.open(address, ios::app);
    
    int error[101];
    int errorPossion[101];
    for(i=0;i<101;i++){
        error[i]=0;
        errorPossion[i]=0;
    }
    while(contigSet!=NULL){
        s++;
        
        long int contigLength = strlen(contigSet->contig);
        j = 0;

        int * gapIndex = new int[contigLength-readLength];                   
        for(i = 0; i<contigLength - readLength; i++){
            SubContig(read,contigSet->contig,i,i+readLength);
            ReverseComplement(read, readRC);
            long int a = SearchRead(read, readSet+readSetIndex);
            long int b = SearchRead(readRC, readSet+readSetIndex);

            if(a==-1&&b==-1){
                all++;
                j++;
                gapIndex[i] = 1;
            }else{
                all1++;
                gapIndex[i] = 0;
            }
            if(i>=length-1){
                gapNumber[j]++;
                j = j - gapIndex[i-length+1];
            }
        }
        j = 0;
        char * region = new char[readLength + length + 1];
        for(i = 0; i<contigLength - readLength; i++){
            SubContig(region,contigSet->contig,i,i+readLength + length);
            char tempChar = region[readLength-1];
            while(1){
                int k = rand();
                if(k%4==0&&tempChar!='A'){
                    region[readLength-1] = 'A';
                    break;
                }
                if(k%4==1&&tempChar!='T'){
                    region[readLength-1] = 'T';
                    break;
                }
                if(k%4==2&&tempChar!='G'){
                    region[readLength-1] = 'G';
                    break;
                }
                if(k%4==3&&tempChar!='C'){
                    region[readLength-1] = 'C';
                    break;
                }
            }
            j = 0;
            for(int t = 0; t<length; t++){
                
                SubContig(read, region, t, t+readLength);
                ReverseComplement(read, readRC);
                long int a = SearchRead(read, readSet+readSetIndex);
                long int b = SearchRead(readRC, readSet+readSetIndex);         
                if(a==-1&&b==-1){
                    j++;
                }
                
            }
            gapNumber1[j]++; 
        }
        
        
        allLength = allLength + contigLength - readLength;
           
        contigSet = contigSet->next;        
    }  
    //cout<<all<<"--"<<all1<<"--"<<allLength<<endl;

    
    double gap = (double)all/(double)allLength;
    double avg = length*gap;
    double sd = sqrt(length*gap*(1-gap));
    
    for(i=0;i<=length;i++){

        if(i>avg+(lambda+1)*sd||i<avg-(lambda+1)*sd){
            percentageOfCorrect = percentageOfCorrect + gapNumber[i];
            percentageOfError = percentageOfError + gapNumber1[i];
        }
    } 
    
    ocout<<s<<"、length:"<<allLength<<endl;
    ocout<<"--gapProblity:"<<1-gap<<"--avg:"<<avg<<"--sd:"<<sd<<"--error:"<<1-((double)percentageOfCorrect/(double)(allLength))<<"--"<<(double)percentageOfError/(double)(allLength)<<endl;;
    
    //delete [] region;
    
    ocout.close();
    return 1;
}

double ComputeSD(GraphNode * distance, long int pos, long int length, long int insertSize,long int var){
    long int i = 0;
    long int j = 0;
    double avg = 0;
    double sd = 0;
    
    double allNumber = 0;
    
    long int outNumber[7];
    for(i = 0; i<7; i++){
        outNumber[i] = 0;
    }
    i = 0;
    //cout<<"ff3"<<endl;
    for(i = pos; i< pos+length; i++){
        GraphNode * temp = distance[i].next;
        while(temp!=NULL){
            long int index = (temp->index-insertSize)/var;
            if(index<-2){
                outNumber[0]++;
            }else if(index<-1){
                outNumber[1]++;
            }else if(index<0){
                outNumber[2]++;
            }else if(index<1){
                outNumber[3]++;
            }else if(index<2){
                outNumber[4]++;
            }else{
                outNumber[5]++;
            }
            temp = temp->next;
            allNumber++;
        }
    }    
    //cout<<"ff4"<<endl;
    
    double all = 0;
    for(j = 0; j<6; j++){
        all = all + pow(outNumber[j],2);
    }
    
    double avg2 = 0;
    avg2 = avg2 + outNumber[0]*allNumber*0.02275;
    avg2 = avg2 + outNumber[1]*allNumber*0.13595;
    avg2 = avg2 + outNumber[2]*allNumber*0.3413;
    avg2 = avg2 + outNumber[3]*allNumber*0.3413;
    avg2 = avg2 + outNumber[4]*allNumber*0.13595;
    avg2 = avg2 + outNumber[5]*allNumber*0.02275;
    
    avg2 = pow(avg2,2);
    
    double avg3 = 0;
    avg3 = avg3 + pow(allNumber*0.02275,2);
    avg3 = avg3 + pow(allNumber*0.13595,2);
    avg3 = avg3 + pow(allNumber*0.3413,2);
    avg3 = avg3 + pow(allNumber*0.3413,2);
    avg3 = avg3 + pow(allNumber*0.13595,2);
    avg3 = avg3 + pow(allNumber*0.02275,2);
       
    //cout<<avg2/(all*avg3)<<endl;
    return avg2/(all*avg3); 
    
}

void GetReferenceCD(ContigSet * contigSet, int length, ReadSet * readSet, int readSetIndex){
    
    long int i = 0;
    long int j = 0;
    long int s = 0;
    
    long int readLength = readSet[readSetIndex].readLength;
    long int insertSize = readSet[readSetIndex].insertSize;
    long int var = readSet[readSetIndex].var;
    
    char read[readLength + 1];
    char readRC[readLength + 1];
    
    long int RM[20];
    long int CD[20];
    for(i = 0; i<20; i++){
        RM[i] = 0;
        CD[i] = 0;
    }
        
    while(contigSet!=NULL){
        
        long int contigLength = strlen(contigSet->contig);
        
        long int * matchNumber = new long int[contigLength - readLength];
        long int * mateNumber = new long int[contigLength - readLength];
        GraphNode * distance = new GraphNode[contigLength - readLength];
        
        i = 0;
        j = 0;
        //cout<<"aa--"<<s<<endl;
        s++;
        
        for(i = 0; i<contigLength - readLength - insertSize - lambda*var; i++){
            //cout<<"i:"<<i<<"--"<<contigLength<<endl;
            SubContig(read, contigSet->contig, i, i+readLength);
            ReverseComplement(read, readRC);
            ReadMate * tempReadMate = SearchRightReadMate(read, readSet+readSetIndex);
            long int a = SearchRead(readRC, readSet + readSetIndex); 
            ReadMate * temp1;
            
            if(tempReadMate==NULL && a==-1){
                matchNumber[i] = -2;
                mateNumber[i] = 0;
                continue;    
            }
        
            if(tempReadMate!=NULL){
                temp1 = tempReadMate;
                tempReadMate = tempReadMate->next;
                delete temp1;
            }
             
            if(tempReadMate==NULL){
                matchNumber[i] = -1;
                mateNumber[i] = 0;
                continue;
            }      
            ReadMate * tempReadMate1 = tempReadMate;
            matchNumber[i] = 0;
            mateNumber[i] = 0;
            bool token = false;
            GraphNode * last = NULL;
            while(tempReadMate!=NULL){
                mateNumber[i]++;                
                long int index1 = KMPIndexOfContigOfMisMatch(contigSet->contig, i + insertSize - lambda*var - readLength, i + insertSize + lambda*var, tempReadMate->readMate);
                //cout<<"index"<<endl;
                if(index1!=-1){
                    
                    if(last==NULL){
                        GraphNode * tempGraphNode = new GraphNode;
                        tempGraphNode->index = index1 - i + readLength;
                        distance[i].next = tempGraphNode;
                        last = tempGraphNode;
                    }else{
                        GraphNode * tempGraphNode = new GraphNode;
                        tempGraphNode->index = index1 - i + readLength;
                        last->next = tempGraphNode;
                        last = tempGraphNode;
                    }
                    
                    matchNumber[i]++;  
                }
                tempReadMate = tempReadMate->next;
            }
       
            while(tempReadMate1!=NULL){
                 temp1 = tempReadMate1;
                 delete [] tempReadMate1->readMate;
                 tempReadMate1 = tempReadMate1->next;
                 delete temp1;
            }
   
        
        }
        i = 0;
        //cout<<"aa1"<<endl;
        
        double mate = 0;
        double match = 0;
        
        for(i = 0; i<contigLength - readLength - length - insertSize - lambda*var; i++){
        
            if(matchNumber[i] > 0){
                match++;
                mate ++;
            }else if(matchNumber[i] == 0){
                mate ++;
            }
            
            if(i>=length-1){
                if(mate!=0){
                    RM[(long int)(10*(match/(mate)))]++;
                    //cout<<match<<"--"<<mate<<"--"<<match/mate<<endl;
                }else{
                    RM[0]++;
                }
                if(matchNumber[i-length+1] > 0){
                    match--;
                    mate --;
                }else if(matchNumber[i-length+1] == 0){
                    mate--;
                }                
            }
            //cout<<"ff"<<endl;
            double cd = ComputeSD(distance, i, length, insertSize, var);
            //cout<<cd<<endl;
            CD[(long int)(10*cd)]++;     
            //cout<<"ff1"<<endl;    
        
        }
        
        delete []matchNumber;
        delete []mateNumber;
        
        contigSet = contigSet->next;
        
    } 
    //cout<<"aa2"<<endl;
    long int allRM = 0;
    long int allCD = 0;
    for(i = 0 ; i<11; i++){
        allRM = allRM + RM[i];
        allCD = allCD + CD[i];
    }   
    //cout<<"aa3"<<endl;
    char address[30];
    strcpy(address, "resultCD.fa");
    ofstream ocout;
    ocout.open(address,ios::app);
    double * rmP = new double[11];
    double * cdP = new double[11]; 
    for(i = 0 ; i<11; i++){
        rmP[i] = (double)RM[i]/(double)allRM;
        cdP[i] = (double)CD[i]/(double)allCD;
    }
    //cout<<"aa4"<<endl;
    ocout<<"RM:";
    for(i = 0 ; i<11; i++){
        double temp = 0;
        for(j = i; j<11; j++){
            temp = temp + rmP[j];
        }
        ocout<<temp<<"--";
    } 
    ocout<<endl;
    ocout<<"CD:";
    for(i = 0 ; i<11; i++){
        double temp = 0;
        for(j = i; j<11; j++){
            temp = temp + cdP[j];
        }
        ocout<<temp<<"--";
    } 
    ocout<<endl;
                         
        
    
}

double nihePossion(long int * gap, int length, double p){
    
    int i = 0;
    int t = 0;
    double s = 1;
    double temp = 0;
    long int all = 0;
    
    for(i = 0; i<length; i++){
        all = all + gap[i];
    }
    
    int j = 0;
    double e = 2.71828;
    
    for(i = 0; i<=length; i++){
        s = 1;
        j = i;
        for(t=i;t>0;t--){
            s = s*t;
        }
        double temp1 = all*((1/pow(e,p))*pow(p,i)/s);  
        double temp2 = pow(temp1 - gap[i],2);    
        temp = temp + (temp2/temp1); 
    }
    
    return temp;
    
    
}


int GetGapProblityPossion(char * contig, int length, ReadSet * readSet, int readSetIndex, double & gap, double & sd){
    long int i = 0;
    long int j = 0;
    long int readLength = readSet[readSetIndex].readLength;
    
    char read[readLength+1];
    char readRC[readLength+1];
    long int gapNumber[length+1];
    
    for(i=0;i<=length;i++){
        gapNumber[i] = 0;
    }   
            
    long int contigLength = strlen(contig);
    
    if(contigLength<readLength + length){
        return 0;
    }
    
    int * gapIndex = new int[contigLength-readLength];                    
        
    for(i = 0; i<contigLength - readLength; i++){
            
        SubContig(read,contig,i,i+readLength);
        ReverseComplement(read, readRC);
        long int a = SearchRead(read, readSet+readSetIndex);
        long int b = SearchRead(readRC, readSet+readSetIndex);
        if(a==-1&&b==-1){
            j++;
            gapIndex[i] = 1;
        }else{
            gapIndex[i] = 0;
        }
        if(i>=length-1){
            gapNumber[j]++;
            j = j - gapIndex[i-length+1];
        }
    }
        
        
    long int gapNumberPossion = 0;
    long int gapZone = 0;
    for(i=0;i<=length;i++){
        gapNumberPossion = gapNumberPossion + gapNumber[i]*i;
        gapZone = gapZone + gapNumber[i];
    }
        
    gap = (double)gapNumberPossion/(double)gapZone;
    sd = sqrt(gap);
    
    delete [] gapIndex;
    
    return 1;
}

int GetInsertSize(char * contig, ReadSet * readSet, int readSetIndex, int & insertSizeLength, int & sd){
    long int i = 0;
    long int j = 0;
    long int readLength = readSet[readSetIndex].readLength;
    
    char read[readLength+1];
    char readRC[readLength+1];
    
    long int len = strlen(contig);
    long int t = 0;
    long int insertLength = 10000;
    
    long int * insertSize = new long int[insertLength];
    for(i=0;i<insertLength;i++){
        insertSize[i]=-1;
    }
    
    for(i=0;i<len-readLength;i++){
        
        SubContig(read,contig,i,i+readLength);
        ReadMate * tempReadMate = SearchRightReadMate(read, readSet+readSetIndex);
        ReverseComplement(read, readRC);
        
        long int a = SearchRead(readRC, readSet+readSetIndex);
        
        if(a == -1 && tempReadMate == NULL){
            continue;
        }
        
        if(tempReadMate!=NULL){
            ReadMate * temp11 = tempReadMate;
            tempReadMate = tempReadMate->next;
            delete temp11;
        }
                      
        while(tempReadMate!=NULL){

            long int index1 = KMPIndexOfContigOfMisMatch(contig,tempReadMate->readMate);
            if(index1>i){
                if(t<insertLength){
                    insertSize[t] = index1 - i + readLength;
                    t++;
                }        
            }
            ReadMate * temp11 = tempReadMate;
            tempReadMate = tempReadMate->next;
            delete []temp11->readMate;
            delete temp11;
        }
        
        if(t>=insertLength){
            break;
        }
               
    }
    
    j = 0;
    long int tempInsertSize = 0;
    long int tempSD = 0;
    for(i=0;i<t;i++){  
        tempInsertSize = tempInsertSize + insertSize[i];
    }
    tempInsertSize = tempInsertSize/t;
    for(i = 0; i < t; i++){
        tempSD = (long int)pow((double)(insertSize[i]-tempInsertSize),2) + tempSD;
    }
    tempSD = (long int)sqrt(tempSD/t);
    
    
    long int tempInsertSize1 = 0;
    long int tempSD1 = 0;
    j = 0;
    for(i=0;i<t;i++){
        if(insertSize[i]<tempInsertSize+(lambda)*tempSD&&insertSize[i]>tempInsertSize-(lambda)*tempSD){
            j++;
            tempInsertSize1 = tempInsertSize1 + insertSize[i];
        }else{
            insertSize[i] = -1;
        }
    }
    tempInsertSize1 = tempInsertSize1/j;
    for(i=0;i<t;i++){
        if(insertSize[i]!=-1){
            tempSD1 = (long int)pow((double)(insertSize[i]-tempInsertSize1),2) + tempSD1; 
        }
    }
    tempSD1 = (long int)sqrt(tempSD1/j);
    
    insertSizeLength = tempInsertSize1;
    sd = tempSD1;
    
    delete [] insertSize;
    
    return 1;
    
    
}

ContigSet * PreTest(){
    
    long int i = 0;
    long int j = 0;
    
    char * ref = new char[30];
    strcpy(ref, "genome.fasta");
    
    ifstream icin;
    icin.open(ref);
    
    char * temp = new char[1000];
    icin.getline(temp, 1000);
    icin.getline(temp, 1000);
    long int len = strlen(temp);
    i++;
    ContigSet * contigSet = new ContigSet;
    long int * lineNumber = new long int[300];
    long int genomeNumber = 1;
    while(icin.getline(temp, 1000)){
        if(temp[0] == '>'){
            lineNumber[genomeNumber-1] = i;
            genomeNumber++;
            i = 0;
            continue;
        }
        i++;
    }
    lineNumber[genomeNumber-1] = i;
    icin.close();
    
    ifstream icin1;
    icin1.open(ref);
    long int t = 0;
    i = 0;
    long int len1 = 0;
    //cout<<"tttttttttt-ssssssss"<<endl;
    while(icin1.getline(temp, 1000)){
        if(temp[0] == '>'){          
            ContigSet * temp = new ContigSet;
            temp->contig = new char[lineNumber[i]*len+1];
            temp->next = contigSet->next;
            contigSet->next = temp;
            t = 0;
            i++;
            //cout<<"tttttttttt----------"<<i<<endl;
            continue;
        }
        len1 = strlen(temp);
        for(j = 0; j<len1; j++){
            contigSet->next->contig[t+j] = temp[j];
        }
        t = t + len1;
        contigSet->next->contig[t] = '\0';
    }
    
    icin1.close();
    return contigSet->next;

}

void GetExProbablity(char * contig, ReadSet * readSet, int readSetIndex){
    long int i = 0;
    long int j = 0;
    
    long int contigLength = strlen(contig);
    long int readLength = readSet[readSetIndex].readLength;
    
    if(contigLength < 1000000){
        return;
    }
    
    char read[readLength+1];
    char readRC[readLength+1];
    
    long int number[100];
    for(i = 0; i<100; i++){
        number[i] = 0;
    }
    
    long int continueousGap = 0;
    long int allExist = 0;
    long int allNonExist = 0;
    
    for(i = 0; i<contigLength - readLength; i++){
        SubContig(read, contig, i, i+readLength);
        ReverseComplement(read, readRC);
        
        long int a = SearchReadNumber(read, readSet+readSetIndex);
        long int b = SearchReadNumber(readRC, readSet+readSetIndex);
        
        if(a==0&&b==0){
            allNonExist++;
            continueousGap++;
        }else{
            allExist = allExist + a + b;
            number[continueousGap]++;
            continueousGap = 0;
        }
        
    }
    
    double b = (double)allExist/(double)(contigLength - readLength);
    double e = 2.71828;
    
    char * tempStr = new char[100];
    strcpy(tempStr, "contineousGap.fa");
    ofstream ocout;
    ocout.open(tempStr, ios::app);
    for(i = 0; i<100; i++){
        ocout<<i<<"--"<<number[i]<<"--"<<pow(e, -b*(i+1))<<"--"<<(double)number[i]/(double)(contigLength - readLength)<<endl;
    }
    
    ocout<<"allExist--"<<allExist<<"--allNonExist:"<<allNonExist<<"--"<<b<<endl;
    
    
}


#endif 
