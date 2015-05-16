#ifndef DFS_H_INCLUDED
#define DFS_H_INCLUDED

#include <cstring>
#include <iostream>
#include "stdio.h"
#include "stdlib.h"

#include "readset.h"
#include "kmerset.h"
#include "common.h"
#include "graph.h"

using namespace std;

#pragma pack(2)
typedef struct DFSNode{
    char * contig;
    struct DFSNode * next;
    long int adjNode;
    long int deep;
    DFSNode(){
        contig = NULL;
        next = NULL;
        adjNode = -1;
        deep = 0;
    }
}DFSNode;
#pragma pack ()

void DeleteDFSNode(DFSNode * head){   
    DFSNode * temp = NULL;
    while(head !=NULL){
        delete [] head->contig;
        head->contig = NULL;
        temp = head;
        head = head->next;
        delete temp;
        temp = NULL;
    }
}

DFSNode * DFS(DBGraph * deBruijnGraph, DFSNode * dfsNode, long int start, long int pos, long int kmerLength, long int length,DFSNode * firstDFSNode,long int index){     
     
     if(start!=0){
         dfsNode->deep++;
         
         if(dfsNode->deep>maxDeep){
             dfsNode->deep--;
             return dfsNode;
         }
         
         char *tempContig = new char[strlen(dfsNode->contig)+strlen(deBruijnGraph[pos].contig)-kmerLength+2];
         AppendRight(tempContig, dfsNode->contig,deBruijnGraph[pos].contig,kmerLength);
         delete []dfsNode->contig;
         dfsNode->contig = tempContig;
         
         tempContig = NULL;
     }
     
     long int len = strlen(dfsNode->contig);
     if(len >length){
         DFSNode * newDFSNode = new DFSNode;
         newDFSNode->adjNode = dfsNode->adjNode; 
         newDFSNode->deep = dfsNode->deep - 1;  
         if(start!=0){
             newDFSNode->contig = new char[strlen(dfsNode->contig)-strlen(deBruijnGraph[pos].contig)+kmerLength];
             SubContig(newDFSNode->contig,dfsNode->contig,0,strlen(dfsNode->contig)-strlen(deBruijnGraph[pos].contig)+kmerLength-1);            
         } 
         dfsNode->next = newDFSNode;
         dfsNode = newDFSNode;
                  
     }else{
         GraphNode * tempGraphNode = deBruijnGraph[pos].outNode;
         while(tempGraphNode!=NULL){                                 
             dfsNode = DFS(deBruijnGraph,dfsNode,1,tempGraphNode->index,kmerLength,length,firstDFSNode, index);   
             tempGraphNode = tempGraphNode->next;
         }
         if(start!=0){
             char * tempChar = new char[strlen(dfsNode->contig)-strlen(deBruijnGraph[pos].contig)+kmerLength];
             SubContig(tempChar,dfsNode->contig,0,strlen(dfsNode->contig)-strlen(deBruijnGraph[pos].contig)+kmerLength-1);
             delete []dfsNode->contig;
             dfsNode->contig = tempChar;
             dfsNode->deep --;
             tempChar = NULL;
         }       
     }
     return dfsNode;    
}

DFSNode * DFSTranverse(DBGraph * deBruijnGraph, long int pos, long int kmerLength, long int length){
    long int i = 0;
    long int j = 0;
    if(NodeCount(deBruijnGraph[pos].outNode)!=0){       
        GraphNode * tempGraphNode = deBruijnGraph[pos].outNode;
        DFSNode *dfsNode = new DFSNode;
        dfsNode->adjNode = tempGraphNode->index;
        dfsNode->contig = new char[strlen(deBruijnGraph[tempGraphNode->index].contig)+1];
        dfsNode->deep = 1;
        CopyContig(dfsNode->contig,deBruijnGraph[tempGraphNode->index].contig);
        DFSNode * tempDFSNode = dfsNode;  
        while(tempGraphNode!=NULL){                             
            dfsNode = DFS(deBruijnGraph,dfsNode,0,tempGraphNode->index,kmerLength,length,tempDFSNode,tempGraphNode->index);             
            tempGraphNode = tempGraphNode->next;
            if(tempGraphNode!=NULL){
                delete []dfsNode->contig;
                dfsNode->contig = new char[strlen(deBruijnGraph[tempGraphNode->index].contig)+1];
                CopyContig(dfsNode->contig,deBruijnGraph[tempGraphNode->index].contig);
                dfsNode->adjNode = tempGraphNode->index;
                dfsNode->deep = 1;
            }
            i++;
        }
        DFSNode * tempDFSNode1 = tempDFSNode;
        DFSNode * previousDFSNode = tempDFSNode;
        long int num = 1;
        while(tempDFSNode1->next!=NULL){
            num++;
            previousDFSNode = tempDFSNode1;
            tempDFSNode1 = tempDFSNode1->next;
        }
        if(num==1){
            return NULL;
        }
        previousDFSNode->next=NULL;
        delete []dfsNode->contig;
        dfsNode->contig = NULL;
        delete dfsNode;
        dfsNode = NULL;
        return tempDFSNode;
    }else{
        return NULL;
    }    
}


DFSNode * DFSLeft(DBGraph * deBruijnGraph, DFSNode * dfsNode, long int start, long int pos, long int kmerLength, long int length,DFSNode * firstDFSNode,long int index){
              
     if(start!=0){
         dfsNode->deep++;
         
         if(dfsNode->deep>maxDeep){
             dfsNode->deep--;
             return dfsNode;
         }
         
         char *tempContig = new char[strlen(dfsNode->contig)+strlen(deBruijnGraph[pos].contig)-kmerLength+2];
         AppendRight(tempContig,deBruijnGraph[pos].contig,dfsNode->contig,kmerLength);
         delete []dfsNode->contig;
         dfsNode->contig = tempContig;
         tempContig = NULL;
     }
     
     long int len = strlen(dfsNode->contig);
     if(len >length){
         DFSNode * newDFSNode = new DFSNode;
         newDFSNode->adjNode = dfsNode->adjNode; 
         newDFSNode->deep = dfsNode->deep-1;  
         if(start!=0){
             newDFSNode->contig = new char[strlen(dfsNode->contig)-strlen(deBruijnGraph[pos].contig)+kmerLength];
             SubContig(newDFSNode->contig,dfsNode->contig,strlen(deBruijnGraph[pos].contig)-kmerLength+1,strlen(dfsNode->contig));            
         } 
         dfsNode->next = newDFSNode;
         dfsNode = newDFSNode;    
           
     }else{
         GraphNode * tempGraphNode = deBruijnGraph[pos].inNode;
         while(tempGraphNode!=NULL){                                 
             dfsNode = DFSLeft(deBruijnGraph,dfsNode,1,tempGraphNode->index,kmerLength,length,firstDFSNode,index);            
             tempGraphNode = tempGraphNode->next;
             
         }
         if(start!=0){
             char * tempChar = new char[strlen(dfsNode->contig)-strlen(deBruijnGraph[pos].contig)+kmerLength];
             SubContig(tempChar,dfsNode->contig,strlen(deBruijnGraph[pos].contig)-kmerLength+1,strlen(dfsNode->contig));
             delete []dfsNode->contig;
             dfsNode->contig = tempChar;
             dfsNode->deep--;
             tempChar = NULL;
         }       
     }
     return dfsNode;    
}

DFSNode * DFSTranverseLeft(DBGraph * deBruijnGraph, long int pos, long int kmerLength, long int length){
    long int i = 0;
    long int j = 0;
    if(NodeCount(deBruijnGraph[pos].inNode)!=0){       
        GraphNode * tempGraphNode = deBruijnGraph[pos].inNode;
        DFSNode *dfsNode = new DFSNode;
        dfsNode->adjNode = tempGraphNode->index;
        dfsNode->contig = new char[strlen(deBruijnGraph[tempGraphNode->index].contig)+1];
        dfsNode->deep = 1;
        CopyContig(dfsNode->contig,deBruijnGraph[tempGraphNode->index].contig);
        DFSNode * tempDFSNode = dfsNode;  
        while(tempGraphNode!=NULL){                             
            dfsNode = DFSLeft(deBruijnGraph,dfsNode,0,tempGraphNode->index,kmerLength,length,tempDFSNode,tempGraphNode->index);             
            tempGraphNode = tempGraphNode->next;
            if(tempGraphNode!=NULL){
                delete []dfsNode->contig;
                dfsNode->contig = new char[strlen(deBruijnGraph[tempGraphNode->index].contig)+1];
                CopyContig(dfsNode->contig,deBruijnGraph[tempGraphNode->index].contig);
                dfsNode->adjNode = tempGraphNode->index;
                dfsNode->deep = 1;
            }
            i++;
        }
        DFSNode * tempDFSNode1 = tempDFSNode;
        DFSNode * previousDFSNode = tempDFSNode;
        long int num = 1;
        while(tempDFSNode1->next!=NULL){
            num++;
            previousDFSNode = tempDFSNode1;
            tempDFSNode1 = tempDFSNode1->next;
        }
        if(num==1){
            return NULL;
        }
        previousDFSNode->next=NULL;
        delete []dfsNode->contig;
        dfsNode->contig = NULL;
        delete dfsNode;
        dfsNode = NULL;
        return tempDFSNode;
    }else{
        return NULL;
    }    
}

#endif
