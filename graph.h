#ifndef GRAPH_H_INCLUDED
#define GRAPH_H_INCLUDED

#include "fstream"
#include <cstring>
#include <iostream>
#include "stdio.h"
#include "stdlib.h"
#include <pthread.h>
#include <math.h>

#include "readset.h"
#include "kmerset.h"
#include "bitarray.h"
#include "common.h"

using namespace std;

#pragma pack(2)
typedef struct GraphNode{
    long int index;
    struct GraphNode * next;
    GraphNode(){
        index = -1;
        next = NULL;
    }
}GraphNode;
#pragma pack ()

#pragma pack(2)
typedef struct DBGraph{
    char * contig;
    GraphNode * outNode;
    GraphNode * inNode;
    DBGraph(){
        contig = NULL;
        outNode = NULL;
        inNode = NULL;
    }
}DBGraph;
#pragma pack ()

#pragma pack(2)
typedef struct DBGraphHead{
    DBGraph * deBruijnGraph;
    unsigned long int nodeNumber;
    DBGraphHead(){
        deBruijnGraph = NULL;
        nodeNumber = 0;
    }
}DBGraphHead;
#pragma pack ()

long int * SearchDBGraph(DBGraphHead * deBruijnGraphHead, char * temp){
    long int i = 0;
    long int j = 0;
    long int len = strlen(temp);
    long int *next = new long int[len+1];
    KMPGetNext(temp, next);
    
    while(i<deBruijnGraphHead->nodeNumber){
        if(deBruijnGraphHead->deBruijnGraph[i].contig==NULL){
            i++;
            continue;
        }
        j = KMPIndexOfContig(deBruijnGraphHead->deBruijnGraph[i].contig,temp, next,1);
        if(j!=-1){
            long int * temp = new long int[3];
            temp[0] = i;
            temp[1] = j;
            temp[2] = strlen(deBruijnGraphHead->deBruijnGraph[i].contig);
            return temp;
        }
        i++;
    }
    
    delete [] next;
    
    return NULL;
}

void DeleteDeBruijnGraphInAndOutNode(GraphNode * temp){
    GraphNode * previous = temp;
    while(temp!=NULL){
        temp = temp->next;
        previous->next = NULL;
        delete previous;
        previous = temp;
    }
}

void DeleteDeBruijnGraph(DBGraphHead * deBruijnGraphHead){
    long int i = 0;
    long int j = 0;
    while(i<deBruijnGraphHead->nodeNumber){
        /*
        if(deBruijnGraphHead->deBruijnGraph[i].contig==NULL){
            i++;
            continue;
        }
        */
        delete [] deBruijnGraphHead->deBruijnGraph[i].contig;
        deBruijnGraphHead->deBruijnGraph[i].contig = NULL;
        DeleteDeBruijnGraphInAndOutNode(deBruijnGraphHead->deBruijnGraph[i].inNode);
        deBruijnGraphHead->deBruijnGraph[i].inNode = NULL;
        DeleteDeBruijnGraphInAndOutNode(deBruijnGraphHead->deBruijnGraph[i].outNode); 
        deBruijnGraphHead->deBruijnGraph[i].outNode = NULL;       
        i++;
    }
}

#pragma pack(2)
typedef struct CreateDBGraphP{
    DBGraphHead * deBruijnGraphHead;
    ReadSet * readSet;
    int kmerLength;
    KmerSetHashTable * kmerSetHashTable;
    unsigned long int kmerSetHashTableCount;
    int setNumber;
    char * kmerHashTableBitArray;
    int threadIndex;
    int totalThreadNumber;
    CreateDBGraphP(){
        deBruijnGraphHead = NULL;
        kmerHashTableBitArray = NULL;
        setNumber = 0;
        threadIndex = 0;
        totalThreadNumber = 0;
    }
}CreateDBGraphP;
#pragma pack ()

void InsertLeftMateNode(DBGraph * deBruijnGraph, long int right, long int left){
    //cout<<"ee"<<endl;
    GraphNode * temp = deBruijnGraph[right].inNode;
    while(temp!=NULL){
        if(temp->index == left){
            return;
        }
        temp = temp->next;
    }
    //cout<<"ee1"<<endl;

    temp = new GraphNode;
    
    temp->index = left;
    //cout<<right<<endl;
    temp->next = deBruijnGraph[right].inNode;
    //cout<<"ee2:"<<left<<endl;
    deBruijnGraph[right].inNode = temp;
}

void InsertRightMateNode(DBGraph * deBruijnGraph, long int left, long int right){
    GraphNode * temp = deBruijnGraph[left].outNode;
    while(temp!=NULL){
        if(temp->index == right){
            return;
        }
        temp = temp->next;
    }
    temp = new GraphNode;
    temp->index = right;
    temp->next = deBruijnGraph[left].outNode;
    deBruijnGraph[left].outNode = temp;
}

pthread_mutex_t mutex4 = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex5 = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex6 = PTHREAD_MUTEX_INITIALIZER;

long int DBGraphInsertNode(long int previousPos, char * kmer, int kmerLength, DBGraph * deBruijnGraph, unsigned long int graphHashCount){
    unsigned long int i = 0;
    
    i = Hash(kmer, kmerLength, graphHashCount);
    int j = 0;
    //cout<<"hash:"<<i<<endl;
    previous:
    while(deBruijnGraph[i].contig!=NULL){
        //cout<<"rr:"<<deBruijnGraph[i].contig<<endl;
        for(j = 0; j<kmerLength; j++){
            if(kmer[j]!=deBruijnGraph[i].contig[j]){
                j = -1;
                break;
            }
        }
        //cout<<"rr1:"<<j<<endl;
        if(j!=-1){         
            if(previousPos == -1 || previousPos == i){
                return i;
            }else{
                pthread_mutex_lock(&mutex4);
                //cout<<"ff"<<endl;
                InsertLeftMateNode(deBruijnGraph, i, previousPos);
                //cout<<"ff1"<<endl;
                InsertRightMateNode(deBruijnGraph, previousPos, i);
                //cout<<"ff2"<<endl;
                pthread_mutex_unlock(&mutex4); 
                return i;
            }
        }
        //cout<<"rr2:"<<j<<endl;
        i = (i + 1)%graphHashCount;
    }
    //cout<<"hash1:"<<i<<endl;
    pthread_mutex_lock(&mutex5);
    if(deBruijnGraph[i].contig==NULL){
        deBruijnGraph[i].contig = new char[kmerLength + 1];
        CopyContig(deBruijnGraph[i].contig, kmer);
        if(previousPos == -1){
            pthread_mutex_unlock(&mutex5);
            return i;
        }else{
            pthread_mutex_lock(&mutex4);
            //cout<<"aa"<<i<<"--"<<previousPos<<endl;
            InsertLeftMateNode(deBruijnGraph, i, previousPos);
            //cout<<"aa1"<<i<<"--"<<previousPos<<endl;
            InsertRightMateNode(deBruijnGraph, previousPos, i);
            //cout<<"aa2"<<i<<"--"<<previousPos<<endl;
            pthread_mutex_unlock(&mutex4); 
            pthread_mutex_unlock(&mutex5);
            return i;
        }
    }else{
        pthread_mutex_unlock(&mutex5);
        goto previous;
    }
    return -1;
    
}

void * CreateDBGraphThread(void * arg){
    CreateDBGraphP * createDBGraphP = (CreateDBGraphP *)arg;
    
    long int readLength = createDBGraphP->readSet->readLength;
    long int kmerLength = createDBGraphP->kmerLength;
    
    long int q = 0;
    long int qR = 0;
    
    unsigned long int p = createDBGraphP->threadIndex * readLength;
    
    char * tempRead = NULL;
    char * tempKmer[readLength - kmerLength + 1];
    char * tempKmerReverseComplement[readLength - kmerLength + 1];
    
    char * tempGlobalKmer[readLength - globalKmerLength + 1];
    char * tempGlobalKmerRC[readLength - globalKmerLength + 1];
    
    int * index = new int[readLength - globalKmerLength + 1];
    int * indexRC = new int[readLength - globalKmerLength + 1];
    
    
    
    unsigned long int i =0;
    long int j = 0;
    int t = 0;
    
    for(i = 0; i<readLength - kmerLength + 1; i++){
        tempKmer[i] = new char[kmerLength + 1];
        tempKmerReverseComplement[i] = new char[kmerLength + 1];
    }
    
    for(i = 0; i<readLength - globalKmerLength + 1; i++){
        tempGlobalKmer[i] = new char[globalKmerLength + 1];
        tempGlobalKmerRC[i] = new char[globalKmerLength + 1];
        index[i] = 0;
        indexRC[i] = 0;
    }

    
    i = createDBGraphP->threadIndex;
    while(i < createDBGraphP->readSet->readNumber){
        //cout<<"xiangya:"<<i<<endl;
        /*
        for(j = 0; j<readLength - globalKmerLength + 1; j++){
            GetBit(createDBGraphP->readSet->readSet, p + j, globalKmerLength, tempGlobalKmer[j]);
            
            t = 0;
            if(KmerSetHashTableSearch(tempGlobalKmer[j], createDBGraphP->kmerSetHashTable, globalKmerLength, createDBGraphP->kmerSetHashTableCount)!=0){
                t = 1;
                continue;
            }
            ReverseComplement(tempGlobalKmer[j], tempGlobalKmerRC[j]);
            if(KmerSetHashTableSearch(tempGlobalKmerRC[j], createDBGraphP->kmerSetHashTable, globalKmerLength, createDBGraphP->kmerSetHashTableCount)!=0){
                t = 1;
                continue;
            }
            
            break;
                        
        }
        if(t != 1){
            i = i + createDBGraphP->totalThreadNumber;
            p = p + createDBGraphP->totalThreadNumber * readLength;
            continue;
        }
        */
        
        for(j = 0; j<readLength - kmerLength + 1; j++){
            GetBit(createDBGraphP->readSet->readSet, p + j, kmerLength, tempKmer[j]);
            ReverseComplement(tempKmer[j], tempKmerReverseComplement[readLength - kmerLength - j]);                    
        } 
        
        for(j = 0; j<readLength - globalKmerLength + 1; j++){
            GetBit(createDBGraphP->readSet->readSet, p + j, globalKmerLength, tempGlobalKmer[j]);
            ReverseComplement(tempGlobalKmer[j], tempGlobalKmerRC[readLength - globalKmerLength - j]); 
            index[j] = KmerSetHashTableSearch(tempGlobalKmer[j], createDBGraphP->kmerSetHashTable, globalKmerLength, createDBGraphP->kmerSetHashTableCount);
            indexRC[readLength - globalKmerLength - j] = KmerSetHashTableSearch(tempGlobalKmerRC[readLength - globalKmerLength - j], createDBGraphP->kmerSetHashTable, globalKmerLength, createDBGraphP->kmerSetHashTableCount);                  
        }      
        
        q = -1;
        qR = -1;
        for(t = 0; t<readLength - kmerLength + 1; t++){
            
            if(t<readLength - globalKmerLength && (index[t] + indexRC[readLength - globalKmerLength - t] != 0)){
                q = DBGraphInsertNode(q, tempKmer[t], kmerLength, createDBGraphP->deBruijnGraphHead->deBruijnGraph, 2*createDBGraphP->kmerSetHashTableCount);
            }else if(t==readLength-globalKmerLength && (index[t] + indexRC[readLength - globalKmerLength - t] == 0)){
                break;
            }else if(t==readLength-globalKmerLength && (index[t] + indexRC[readLength - globalKmerLength - t] != 0)){
                q = DBGraphInsertNode(q, tempKmer[t], kmerLength, createDBGraphP->deBruijnGraphHead->deBruijnGraph, 2*createDBGraphP->kmerSetHashTableCount);
            }else if(t>readLength - globalKmerLength){
                q = DBGraphInsertNode(q, tempKmer[t], kmerLength, createDBGraphP->deBruijnGraphHead->deBruijnGraph, 2*createDBGraphP->kmerSetHashTableCount);
            }else{
                q = -1;
            }        
        }
        
        for(t = 0; t<readLength - kmerLength + 1; t++){
            
            if(t<readLength - globalKmerLength && (indexRC[t] + index[readLength - globalKmerLength - t] != 0)){
                qR = DBGraphInsertNode(qR, tempKmerReverseComplement[t], kmerLength, createDBGraphP->deBruijnGraphHead->deBruijnGraph, 2*createDBGraphP->kmerSetHashTableCount);
            }else if(t==readLength-globalKmerLength && (indexRC[t] + index[readLength - globalKmerLength - t] == 0)){
                break;
            }else if(t==readLength-globalKmerLength && (indexRC[t] + index[readLength - globalKmerLength - t] != 0)){
                qR = DBGraphInsertNode(qR, tempKmerReverseComplement[t], kmerLength, createDBGraphP->deBruijnGraphHead->deBruijnGraph, 2*createDBGraphP->kmerSetHashTableCount);
            }else if(t>readLength - globalKmerLength){
                qR = DBGraphInsertNode(qR, tempKmerReverseComplement[t], kmerLength, createDBGraphP->deBruijnGraphHead->deBruijnGraph, 2*createDBGraphP->kmerSetHashTableCount);
            }else{
                qR = -1;
            }
        
        }
        i = i + createDBGraphP->totalThreadNumber;
        p = p + createDBGraphP->totalThreadNumber * readLength;
    }
    
    for(i = 0; i<readLength - kmerLength + 1; i++){
        delete [] tempKmer[i];
        delete [] tempKmerReverseComplement[i];
        delete [] tempGlobalKmer[i];
        delete [] tempGlobalKmerRC[i];
    }
    delete [] index;
    delete [] indexRC;

    
}

void CreateDBGraph(ReadSet * readSet, int kmerLength, KmerSetHashTable * kmerSetHashTable, unsigned long int kmerSetHashTableCount, int setNumber, DBGraphHead * deBruijnGraphHead, int totalThreadNumber){
    
    pthread_t tid[totalThreadNumber];
    
    int i = 0;
    unsigned long int allKmerNumber = 0;
    CreateDBGraphP * createDBGraphP = new CreateDBGraphP[totalThreadNumber];

    for(i = 0; i< totalThreadNumber; i++){
        createDBGraphP[i].deBruijnGraphHead = deBruijnGraphHead;
        createDBGraphP[i].readSet = readSet;
        createDBGraphP[i].kmerLength = kmerLength;
        createDBGraphP[i].kmerSetHashTable = kmerSetHashTable;
        createDBGraphP[i].kmerSetHashTableCount = kmerSetHashTableCount;
        createDBGraphP[i].setNumber = setNumber;
        createDBGraphP[i].threadIndex = i;
        createDBGraphP[i].totalThreadNumber = totalThreadNumber;
             
        if(pthread_create(&tid[i], NULL, CreateDBGraphThread, (void *)&createDBGraphP[i])){
            cout<<"create thread wrong!"<<endl;
            exit(1);
        }      
    }
    
    for(i = 0; i < totalThreadNumber; i++){
        pthread_join(tid[i], NULL);
    }
    
    
}

int DBGSearchNode(GraphNode * tempGraphNode, long int num, long int newNum){
    while(tempGraphNode!=NULL){
        if(tempGraphNode->index==num){
            tempGraphNode->index = newNum;
            return 1;
        }
        tempGraphNode = tempGraphNode->next;
    }
    return 0;
}

int DBGSearchNode(GraphNode * tempGraphNode, long int num){
    while(tempGraphNode!=NULL){
        if(tempGraphNode->index==num){
            return 1;
        }
        tempGraphNode = tempGraphNode->next;
    }
    return 0;
}

int NodeCount(GraphNode * temp){
    int i = 0;
    while(temp!=NULL){
        i++;
        temp = temp->next;
    }
    return i;
}

int DBGRemoveInNode(DBGraph * graph, long int num){
    GraphNode * tempGraphNode = graph->inNode;
    GraphNode * temp = NULL;
    while(tempGraphNode!=NULL){
        if(tempGraphNode->index==num){
            if(temp==NULL){
                graph->inNode = tempGraphNode->next;
                tempGraphNode->next = NULL;
                delete tempGraphNode;
                return 1;
            }
            temp->next = tempGraphNode->next;
            tempGraphNode->next = NULL;
            delete tempGraphNode; 
            return 1;
        }
        temp = tempGraphNode;
        tempGraphNode = tempGraphNode->next;
    }
    return 0;
}

int DBGRemoveOutNode(DBGraph * graph, long int num){
    GraphNode * tempGraphNode = graph->outNode;
    GraphNode * temp = NULL;
    while(tempGraphNode!=NULL){
        if(tempGraphNode->index==num){
            if(temp==NULL){
                graph->outNode = tempGraphNode->next;
                tempGraphNode->next = NULL;
                delete tempGraphNode;
                return 1;
            }
            temp->next = tempGraphNode->next;
            tempGraphNode->next = NULL;
            delete tempGraphNode; 
            return 1;
        }
        temp = tempGraphNode;
        tempGraphNode = tempGraphNode->next;
    }
    return 0;
}

void AdjustDBG(DBGraph *deBruijnGraph, int kmerLength,long int max){
    long int i = 0;
    long int j = 0;
    long int m = 0;
    long int n = 0;
    GraphNode * tempGraphNode = NULL;
    
    while(i<max){
        if(deBruijnGraph[i].contig==NULL){
            i++;
            j++;
            continue;
        }
        if(j>0){
            deBruijnGraph[i-j].contig = deBruijnGraph[i].contig;
            deBruijnGraph[i-j].inNode = deBruijnGraph[i].inNode;
            deBruijnGraph[i-j].outNode = deBruijnGraph[i].outNode;
            deBruijnGraph[i].contig = NULL;
            deBruijnGraph[i].inNode = NULL;
            deBruijnGraph[i].outNode = NULL;
            tempGraphNode = deBruijnGraph[i-j].inNode;
            while(tempGraphNode!=NULL){
                if(DBGSearchNode(deBruijnGraph[tempGraphNode->index].outNode,i,i-j)!=1){
                    //cout<<i<<"--"<<j<<"--"<<deBruijnGraph[i].contig<<"--"<<deBruijnGraph[tempGraphNode->index].contig<<endl;
                    cout<<"1SimplePathMerge Error!"<<endl;
                    exit(0);
                }
                tempGraphNode = tempGraphNode->next;
            }
            tempGraphNode = deBruijnGraph[i-j].outNode;
            while(tempGraphNode!=NULL){
                if(DBGSearchNode(deBruijnGraph[tempGraphNode->index].inNode,i,i-j)!=1){
                    //cout<<i<<"--"<<j<<"--"<<deBruijnGraph[i].contig<<"--"<<deBruijnGraph[tempGraphNode->index].contig<<endl;
                    cout<<"2SimplePathMerge Error!"<<endl;
                    exit(0);
                }
                tempGraphNode = tempGraphNode->next;
            }
        }
        i++;
    }
    //cout<<"AdjustDBG End!"<<endl;
}


void SimplePathMerge(DBGraph *deBruijnGraph, int kmerLength, long int max){
    long int i = 0;
    long int j = 0;
    long int m = 0;
    long int n = 0;
    long int t = 0;
    long int previousIndex = -1;
    GraphNode *tempGraphNode = NULL;
    while(i<max){
        
        if(deBruijnGraph[i].contig==NULL){
            i++;
            continue;
        }
        while(NodeCount(deBruijnGraph[i].inNode)==1){
            previousIndex = deBruijnGraph[i].inNode->index;
            if(NodeCount(deBruijnGraph[previousIndex].outNode)!=1||DBGSearchNode(deBruijnGraph[previousIndex].inNode,i)==1){
                break;
            }
            t = 0;
            n = strlen(deBruijnGraph[previousIndex].contig);
            m = strlen(deBruijnGraph[i].contig);
            //cout<<"m:"<<m<<"--n:"<<n<<"current:"<<i<<"--pre--"<<previousIndex<<"--"<<deBruijnGraph[previousIndex].contig<<endl;
            char *tempContig = new char[m+n-kmerLength + 2];
            AppendRight(tempContig,deBruijnGraph[previousIndex].contig, deBruijnGraph[i].contig, kmerLength);
            
            delete deBruijnGraph[i].contig;
            deBruijnGraph[i].contig = tempContig;
            delete deBruijnGraph[i].inNode;
            deBruijnGraph[i].inNode = deBruijnGraph[previousIndex].inNode;
            deBruijnGraph[previousIndex].inNode = NULL;
            delete deBruijnGraph[previousIndex].outNode;
            deBruijnGraph[previousIndex].outNode = NULL;
            delete deBruijnGraph[previousIndex].contig;
            deBruijnGraph[previousIndex].contig = NULL;

            tempGraphNode = deBruijnGraph[i].inNode;
            while(tempGraphNode!=NULL){
                if(DBGSearchNode(deBruijnGraph[tempGraphNode->index].outNode,previousIndex,i)!=1){
                    cout<<"0SimplePathMerge Error!"<<endl;
                    exit(0);
                }
                tempGraphNode = tempGraphNode->next;
            }
        }
        i++;
    }
         
    //cout<<"SimplePathMerge End!"<<endl;
    
}

void RemoveTip(DBGraph *deBruijnGraph, int kmerLength, long int graphNodeNumber){
    long int i = 0;
    long int j = 0;
    long int inNodeCount = 0;
    long int outNodeCount = 0;
    GraphNode * tempGraphNode = new GraphNode;
    GraphNode * first = tempGraphNode;
    while(i<graphNodeNumber){
        if(deBruijnGraph[i].contig==NULL || strlen(deBruijnGraph[i].contig)>=2*kmerLength){
            i++;
            continue;
        }
        inNodeCount = NodeCount(deBruijnGraph[i].inNode);
        outNodeCount = NodeCount(deBruijnGraph[i].outNode);
        if((inNodeCount==1&&outNodeCount==0)||(inNodeCount==0&&outNodeCount==1)||(inNodeCount==0&&outNodeCount==0)){
            GraphNode * newGraphNode = new GraphNode;
            newGraphNode->index = i;
            tempGraphNode->next = newGraphNode;
            tempGraphNode = newGraphNode;
        }
        i++;
    }
    first = first->next;
    while(first!=NULL){
        inNodeCount = NodeCount(deBruijnGraph[first->index].inNode);
        outNodeCount = NodeCount(deBruijnGraph[first->index].outNode);
        if(inNodeCount==1&&outNodeCount==0){
            j = deBruijnGraph[first->index].inNode->index;
            DBGRemoveOutNode(deBruijnGraph+j,first->index);
            delete deBruijnGraph[first->index].inNode;
            deBruijnGraph[first->index].inNode = NULL;
            delete deBruijnGraph[first->index].contig;
            deBruijnGraph[first->index].contig = NULL; 
        }
        if(inNodeCount==0&&outNodeCount==1){
            j = deBruijnGraph[first->index].outNode->index;
            DBGRemoveInNode(deBruijnGraph+j,first->index);
            delete deBruijnGraph[first->index].outNode;
            deBruijnGraph[first->index].outNode = NULL;
            delete deBruijnGraph[first->index].contig;
            deBruijnGraph[first->index].contig = NULL; 
        }
        if(inNodeCount==0&&outNodeCount==0){
            delete deBruijnGraph[first->index].contig;
            deBruijnGraph[first->index].contig = NULL; 
        }
        first = first->next;
    }
    //cout<<"RemoveTip End!"<<endl;
}

void RemoveCycle(DBGraph *deBruijnGraph, int kmerLength, long int graphNodeNumber, KmerSetHashTableHead * kmerSetHashTableHead){
    long int i = 0;
    long int j = 0;
    long int m = 0;
    long int n = 0;
    long int p = 0;
    long int q = 0;
    long int t = 0;
    int repeat = 2;
    long int inNodeCount = 0;
    long int outNodeCount = 0;
    while(i<graphNodeNumber){
        if(deBruijnGraph[i].contig==NULL){
            i++;
            continue;
        }
        inNodeCount = NodeCount(deBruijnGraph[i].inNode);
        outNodeCount = NodeCount(deBruijnGraph[i].outNode);
        if(inNodeCount==1&&outNodeCount==1){
                                                      
            j = deBruijnGraph[i].inNode->index;
            if(deBruijnGraph[i].inNode->index == deBruijnGraph[i].outNode->index){
                m = strlen(deBruijnGraph[i].contig);
                n = strlen(deBruijnGraph[j].contig);
                repeat = 2;
                double avg = 0;
                for(t=0; t<kmerSetHashTableHead->setNumber;t++){
                    double avgKmerFrequency = GetAvgKmerFrequency(deBruijnGraph[i].contig, t, 0, m-kmerSetHashTableHead->kmerLength+1, kmerSetHashTableHead);
                    avgKmerFrequency = avgKmerFrequency/kmerSetHashTableHead->avgKmerFrequency[t];
                    avg = avg + avgKmerFrequency;
                }
                avg = avg/kmerSetHashTableHead->setNumber;
                if(avg<3){
                    repeat = 1;
                }
                char *tempContig = NULL;
                if(repeat == 1){
                    tempContig = new char[2*n+m-2*kmerLength + 3];
                    q = 0;
                    t = 0;
                    for(p = 0;p<2*n+m-2*kmerLength + 2;p++){
                        if(p<n){
                            tempContig[p] = deBruijnGraph[j].contig[p];
                        }else if(p<m+n-kmerLength+1){
                            tempContig[p] = deBruijnGraph[i].contig[kmerLength-1+q];
                            q++;
                        }else{
                            tempContig[p] = deBruijnGraph[j].contig[kmerLength-1+t];
                            t++;
                        }
                    }
                    tempContig[p] = '\0';
                }else{
                    
                    tempContig = new char[3*n+2*m-4*kmerLength + 5];
                    AppendRight(tempContig,deBruijnGraph[j].contig,deBruijnGraph[i].contig,kmerLength);
                    AppendRight(tempContig,tempContig,deBruijnGraph[j].contig,kmerLength);
                    AppendRight(tempContig,tempContig,deBruijnGraph[i].contig,kmerLength);
                    AppendRight(tempContig,tempContig,deBruijnGraph[j].contig,kmerLength);
                }
                
                delete deBruijnGraph[i].contig;
                delete deBruijnGraph[i].outNode;
                delete deBruijnGraph[i].inNode;
                deBruijnGraph[i].contig = NULL;
                deBruijnGraph[i].outNode = NULL;
                deBruijnGraph[i].inNode = NULL;
                delete deBruijnGraph[j].contig;
                deBruijnGraph[j].contig = tempContig;
                DBGRemoveInNode(deBruijnGraph+j,i);
                DBGRemoveOutNode(deBruijnGraph+j,i);
                
            }
        }
        i++;
    }
    //cout<<"RemoveCycle End!"<<endl;
}

void RemoveCycle1(DBGraph *deBruijnGraph, int kmerLength, long int graphNodeNumber){
    long int i = 0;
    long int j = 0;
    long int m = 0;
    long int n = 0;
    long int p = 0;
    long int q = 0;
    long int t = 0;
    long int inNodeCount = 0;
    long int outNodeCount = 0;

    while(i<graphNodeNumber){
        if(deBruijnGraph[i].contig==NULL){
            i++;
            continue;
        }
        inNodeCount = NodeCount(deBruijnGraph[i].inNode);
        outNodeCount = NodeCount(deBruijnGraph[i].outNode);
        if(inNodeCount==2&&outNodeCount==1){
            j = deBruijnGraph[i].outNode->index;
            if(!(NodeCount(deBruijnGraph[j].inNode)==1&&NodeCount(deBruijnGraph[j].outNode)==2)){
                i++;
                continue;
            }
            if(DBGSearchNode(deBruijnGraph[i].inNode,j)==0||DBGSearchNode(deBruijnGraph[j].outNode,i)==0){
                i++;
                continue;
            }
           
            m = strlen(deBruijnGraph[i].contig);
            n = strlen(deBruijnGraph[j].contig);
            char *tempContig = new char[n+m-kmerLength + 2];
            AppendRight(tempContig,deBruijnGraph[i].contig,deBruijnGraph[j].contig,kmerLength);
            char * tempContig1 = new char[2*strlen(tempContig)-kmerLength + 2];
            AppendRight(tempContig1,tempContig,tempContig,kmerLength);
            delete []deBruijnGraph[i].contig;
            delete tempContig;
            deBruijnGraph[i].contig = tempContig1;
            
            delete deBruijnGraph[i].outNode;
            deBruijnGraph[i].outNode = NULL;
            DBGRemoveInNode(deBruijnGraph + i,j);
            delete deBruijnGraph[j].inNode;
            deBruijnGraph[j].inNode = NULL;
            DBGRemoveOutNode(deBruijnGraph + j,i);
            deBruijnGraph[i].outNode = deBruijnGraph[j].outNode;
            DBGSearchNode(deBruijnGraph[deBruijnGraph[i].outNode->index].inNode,j,i);
            deBruijnGraph[j].outNode = NULL;
            delete []deBruijnGraph[j].contig;
            deBruijnGraph[j].contig = NULL;
                        
        }
        i++;
    }
    //cout<<"RemoveCycle1 End!"<<endl;
}

DBGraph * OptimizeDBGraphSpace(DBGraph * deBruijnGraph, unsigned long int & graphNodeCount){
    long int i = 0;
    long int j = 0;
    graphNodeCount = 0;
    while(deBruijnGraph[graphNodeCount].contig!=NULL){
        graphNodeCount++;
    }

    DBGraph * newDeBruijnGraph = new DBGraph[graphNodeCount+1];
    for(i=0;i<graphNodeCount;i++){
        newDeBruijnGraph[i].contig = deBruijnGraph[i].contig;
        newDeBruijnGraph[i].inNode = deBruijnGraph[i].inNode;
        newDeBruijnGraph[i].outNode = deBruijnGraph[i].outNode;
        deBruijnGraph[i].contig = NULL;
        deBruijnGraph[i].inNode = NULL;
        deBruijnGraph[i].outNode = NULL;
    }
    delete [] deBruijnGraph;
    return newDeBruijnGraph;
}

void WriteDBGraph(DBGraph * deBruijnGraph, long int count, char * str){
    long int i = 0;
    long int j = 0;
    long int m = 0;
    long int n = 0;
    GraphNode * tempGraphNode;
    ofstream ocout;
    ocout.open(str);
    while(i<count){
        if(deBruijnGraph[i].contig == NULL){
            i++;
            continue;
        }
        tempGraphNode = deBruijnGraph[i].inNode;
        ocout<<">"<<i<<"-contig:inNodeCount:"<<NodeCount(tempGraphNode)<<"---";
        while(tempGraphNode!=NULL){
            ocout<<tempGraphNode->index<<",";
            tempGraphNode = tempGraphNode->next;
        }
        tempGraphNode = deBruijnGraph[i].outNode;
        ocout<<"outNodeCount:"<<NodeCount(tempGraphNode)<<"---";
        while(tempGraphNode!=NULL){
            ocout<<tempGraphNode->index<<",";
            tempGraphNode = tempGraphNode->next;
        }
        ocout<<"len:"<<strlen(deBruijnGraph[i].contig)<<endl;
        ocout<<deBruijnGraph[i].contig<<endl;
        i++;   
    }

}

void GetDBGraphFromAddress(DBGraphHead * deBruijnGraphHead, char * str){
    long int i = 0;
    long int j = 0;
    long int m = 0;
    long int n = 0;
    long int p = 0;
    long int q = 0;
    long int graphNodeCount = 0;
    char *temp = new char[10000000];   	
	ifstream icin;
	icin.open(str);
	while(icin.getline(temp,10000000)){
		i++;
	}
	icin.close();
	ifstream icin1;
	icin1.open(str);
	graphNodeCount = (long int)(i/2);
	deBruijnGraphHead->nodeNumber = graphNodeCount;
    DBGraph * deBruijnGraph = new DBGraph[graphNodeCount+1];
    deBruijnGraphHead->deBruijnGraph = deBruijnGraph;
    i = 0;
    while(icin1.getline(temp,10000000)){
        int tempLength = strlen(temp);
		if(i%2!=0){
            deBruijnGraph[j].contig = new char[tempLength+1];
            strcpy(deBruijnGraph[j].contig,temp);
            deBruijnGraph[j].contig[tempLength] = '\0';
            j++;
            
        }else{
            n = 0;
            for(m=0;m<tempLength;m++){
                if(temp[m]==':'&&n==0){
                    n = 1;
                    continue;
                }
                if(temp[m]==':'&&n==1){
                   if(temp[m+1]!='0'){
                       m = m + 5;
                       while(temp[m]!='o'){
                           q = 0;
                           char * tempNum = new char[10];
                           while(temp[m]!=','){
                               tempNum[q] = temp[m];
                               m++;
                               q++;
                           }
                           tempNum[q] = '\0';
                           long int num = atoi(tempNum);
                           InsertLeftMateNode(deBruijnGraph,j,num);
                           m++;
                           delete tempNum;
                       }
                       
                   }
                   n++;
                   continue;
                }
                if(temp[m]==':'&&n==2){
                   if(temp[m+1]!='0'){
                       m = m + 5;
                       while(temp[m]!='l'){
                           q = 0;
                           char * tempNum = new char[10];
                           while(temp[m]!=','){
                               tempNum[q] = temp[m];
                               m++;
                               q++;
                           }
                           tempNum[q] = '\0';
                           long int num = atoi(tempNum);
                           InsertRightMateNode(deBruijnGraph,j,num);
                           m++;
                           delete tempNum;
                       }
                       
                   }
                   n++;
                }
            }
        }
        i++;
	}
	delete [] temp;	

}

void GetAvgKmerNumberAndGapProblity(DBGraph* deBruijnGraph, ReadSet * readSet, int readSetIndex, KmerSetHashTable * kmerSetHashTable,long int kmerSetHashTableCount, int kmerLength){
    long int i = 0;
    long int j = 0;
    long int max = 0;
    long int all = 0;
    long int tempAll = 0;
    long int readLength = readSet[readSetIndex].readLength;
    char read[readLength+1];
    while(deBruijnGraph[i].contig!=NULL){
        if(strlen(deBruijnGraph[i].contig)>j){
            max = i;
            j = strlen(deBruijnGraph[i].contig);
        }
        i++;
    }
    long int len =strlen(deBruijnGraph[max].contig);
    
    long int * insertSize = new long int[10*len];
    for(i=0;i<10*len;i++){
        insertSize[i]=-1;
    }
    long int t = 0;
    
    char * temp = new char[kmerLength + 1];
    char * tempReverseComplement = new char[kmerLength + 1];
    long int kmerNumber = 0;
    for(i=0;i<len-kmerLength;i++){
        SubContig(temp,deBruijnGraph[max].contig,i,i+kmerLength);
        unsigned long long int ss = KmerSetHashTableSearch(temp,kmerSetHashTable,kmerLength,kmerSetHashTableCount);
        if(ss != 0){
           all = all + kmerSetHashTable[ss-1].frequency[readSetIndex];
           kmerNumber++;
        }
        ReverseComplement(temp, tempReverseComplement);
        ss = KmerSetHashTableSearch(tempReverseComplement,kmerSetHashTable,kmerLength,kmerSetHashTableCount);
        if(ss != 0){
            all = all + kmerSetHashTable[ss-1].frequency[readSetIndex];
            kmerNumber++;
        }   
        
    }

    double avg = (double)all/(double)(kmerNumber);

    char * temp1 = new char[readLength + 1];
    delete [] tempReverseComplement;
    tempReverseComplement = new char[readLength + 1];
    all = 0;
    
    int * gapIndex = new int[len-readLength];  
    long int length = 21;
    j = 0;
    long int gapNumber[length+1];   
    for(i=0;i<=length;i++){
        gapNumber[i] = 0;
    }
    
    long int allReadNumber = 0;
    
    for(i=0;i<len-readLength;i++){
        
        SubContig(temp1,deBruijnGraph[max].contig,i,i+readLength);
        ReadMate * tempReadMate = SearchRightReadMate(temp1, readSet+readSetIndex);
        ReverseComplement(temp1, tempReverseComplement);
        
        long int a = SearchReadNumber(temp1, readSet+readSetIndex);
        long int b = SearchReadNumber(tempReverseComplement, readSet + readSetIndex);
        
        if(a == 0&& b == 0){
            all ++;
            j++;
            gapIndex[i] = 1;
        }else{
            allReadNumber= allReadNumber + a + b;
            gapIndex[i] = 0;
        }
        
        if(i>=length-1){
            gapNumber[j]++;
            j = j - gapIndex[i-length+1];
        }
        
        if(tempReadMate!=NULL){
            ReadMate * temp11 = tempReadMate;
            tempReadMate = tempReadMate->next;
            delete temp11;
        }
                      
        while(tempReadMate!=NULL){

            long int index1 = KMPIndexOfContigOfMisMatch(deBruijnGraph[max].contig,tempReadMate->readMate);
            if(index1>i){
                if(t<10*len){
                    insertSize[t] = index1 - i + readLength;
                    t++;
                }        
            }
            ReadMate * temp11 = tempReadMate;
            tempReadMate = tempReadMate->next;
            delete []temp11->readMate;
            delete temp11;
        }
               
    }
    
    
    long int possionNumber = 0;
    for(i=0;i<=length;i++){
        possionNumber = possionNumber + i*gapNumber[i];     
    }
    
    double gapPossion = (double)possionNumber/(double)(len - readLength -length);
    double sdPossion = sqrt(gapPossion);
    

    double gap = (double)all/(double)(len-readLength);
    j = 0;
    long int tempInsertSize = 0;
    long int tempSD = 0;
    for(i=0;i<t;i++){  
        tempInsertSize = tempInsertSize + insertSize[i];
    }

    /*
    tempInsertSize = tempInsertSize/t;
    for(i=0;i<t;i++){
        tempSD = (long int)pow((double)(insertSize[i]-tempInsertSize),2) + tempSD;
    }
    tempSD = (long int)sqrt(tempSD/t);
    cout<<"1:"<<tempInsertSize<<"--2:"<<tempSD<<endl;
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
    */
    //(readSet+readSetIndex)->insertSize = tempInsertSize1;
    //(readSet+readSetIndex)->var = tempSD1;
    (readSet+readSetIndex)->gapAvg = gapPossion;
    (readSet+readSetIndex)->gapSD = sdPossion;
    (readSet+readSetIndex)->gapProblity = gap;
    (readSet+readSetIndex)->avgKmerNumber = avg;
    (readSet+readSetIndex)->avgReadNumber = (double)allReadNumber/(double)(len-readLength);
    
    
    delete []temp;
    delete []temp1;
    delete []insertSize;
}



void InitDBGraph(ReadSet * readSet, int kmerLength, KmerSetHashTableHead* kmerSetHashTableHead,DBGraphHead * deBruijnGraphHead, int totalThreadNumber){
    unsigned long int i = 0;
    unsigned long int j = 0;
    int num = 0;
    KmerSetHashTable* kmerSetHashTable = kmerSetHashTableHead->kmerSetHashTable;
    unsigned long int kmerSetHashTableCount = kmerSetHashTableHead->kmerSetHashTableCount;
    long int setNumber = kmerSetHashTableHead->setNumber;
    
    

    deBruijnGraphHead->nodeNumber = 0;
    long int graphNode = 8*kmerSetHashTableCount;
    DBGraph * deBruijnGraph = new DBGraph[graphNode];
    deBruijnGraphHead->deBruijnGraph = deBruijnGraph;
    for(num = 0; num < setNumber; num++){  
        CreateDBGraph(readSet + num, kmerLength, kmerSetHashTable, kmerSetHashTableCount, setNumber, deBruijnGraphHead,totalThreadNumber);
    } 
    char * temp = new char[30];
    strcpy(temp, "graph.fa");
    //WriteDBGraph(deBruijnGraph, graphNode, temp);
    SimplePathMerge(deBruijnGraph,kmerLength,graphNode); 
    char * temp1 = new char[30];
    strcpy(temp1, "graph1.fa");
    //WriteDBGraph(deBruijnGraph, graphNode, temp1);
    AdjustDBG(deBruijnGraph,kmerLength,graphNode);
    char * temp2 = new char[30];
    strcpy(temp2, "graph2.fa");
    //WriteDBGraph(deBruijnGraph, graphNode, temp2);
    for(i=0;i<5;i++){
        RemoveTip(deBruijnGraph,kmerLength, graphNode);
        SimplePathMerge(deBruijnGraph,kmerLength,graphNode);
        AdjustDBG(deBruijnGraph,kmerLength,graphNode);
    }
    char * temp3 = new char[30];
    strcpy(temp3, "graph3.fa");
    //WriteDBGraph(deBruijnGraph, graphNode, temp3);  
    RemoveCycle(deBruijnGraph,kmerLength,graphNode, kmerSetHashTableHead);
    SimplePathMerge(deBruijnGraph,kmerLength,graphNode);
    AdjustDBG(deBruijnGraph,kmerLength,graphNode);
    RemoveCycle1(deBruijnGraph,kmerLength,graphNode);
    SimplePathMerge(deBruijnGraph,kmerLength,graphNode);
    AdjustDBG(deBruijnGraph,kmerLength,graphNode);
    char * temp4 = new char[30];
    strcpy(temp4, "graph4.fa");
    //WriteDBGraph(deBruijnGraph, graphNode, temp4);
    deBruijnGraphHead->deBruijnGraph = OptimizeDBGraphSpace(deBruijnGraph, deBruijnGraphHead->nodeNumber);
    char * temp5 = new char[30];
    strcpy(temp5, "graph5.fa");
    WriteDBGraph(deBruijnGraphHead->deBruijnGraph, deBruijnGraphHead->nodeNumber, temp5);
    
     
}




#endif
