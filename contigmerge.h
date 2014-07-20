#ifndef CONTIGMERGE_H_INCLUDED 
#define CONTIGMERGE_H_INCLUDED 

#include <cstring>
#include <iostream>
#include "stdio.h"
#include "stdlib.h"
#include <math.h>

#include "readset.h"
#include "constructcontigset.h"


double WeightContigMerge(char * contigLeft, char * contigRight, ReadSet * readSet, long int readSetIndex,long int kmerLength, long int overlap, long int length){
    
    long int i = 0;
    long int j = 0;
    long int readLength = readSet[readSetIndex].readLength;
    long int insertSize = readSet[readSetIndex].insertSize;
    long int var = readSet[readSetIndex].var;
    
    long int nullNumber = 0;
    long int positionMatchNumber = 0;
    long int positionNonExist = 0;
    
    long int min = 0;
    long int index1 = 0;

    char * tempContigLeft = contigLeft;
    long int tempLengthLeft = strlen(contigLeft);
    if(tempLengthLeft>=insertSize + lambda*var){
        min = insertSize + lambda*var;
        tempContigLeft = new char[min+1];
        SubContig(tempContigLeft,contigLeft,tempLengthLeft-min,tempLengthLeft);
        tempLengthLeft = min;
    }

    char * tempContigRight = contigRight;
    long int tempLengthRight = strlen(contigRight);
    if(tempLengthRight>=insertSize + lambda*var){
        min = insertSize + lambda*var;
        tempContigRight = new char[min+1];
        SubContig(tempContigRight,contigRight,0,min);
        tempLengthRight = min;
    }

    char * temp = new char[readLength + 1];
    char * tempRC = new char[readLength + 1];
    
    ReadMate * temp1;
    
    long int * gapIndex = new long int[tempLengthRight];
    long int * nullIndex = new long int[tempLengthRight];
    long int * matchIndex = new long int[tempLengthRight];
    
    for(j=overlap-kmerLength;j<tempLengthRight - readLength - length && j<insertSize - lambda*var - readLength && tempLengthLeft > insertSize;j++){
        index1 = 0;
        matchIndex[j] = 0;
        
        SubContig(temp, tempContigRight, j, j + readLength);
        
        ReadMate * tempReadMate = SearchLeftReadMate(temp, readSet+readSetIndex);
        long int a = SearchRead(temp, readSet + readSetIndex);
        
        if(tempReadMate==NULL && a ==-1){
            positionNonExist++;
            gapIndex[j] = 1;
            index1 = 1;
        }else{
            gapIndex[j] = 0;
        }
        
        if(tempReadMate!=NULL){
            temp1 = tempReadMate;
            tempReadMate = tempReadMate->next;
            delete temp1;
        }   
        
        if(tempReadMate==NULL&&index1!=1){
            nullIndex[j] = 1;
            nullNumber++;
        }else{
            nullIndex[j] = 0;
        }
        
        bool token = false;

        while(tempReadMate!=NULL){
            index1 = KMPIndexOfContigOfMisMatch(tempContigLeft,tempReadMate->readMate);
            long int tempDD = j + tempLengthLeft - index1 + readLength - overlap;
            if(index1!=-1 && tempDD<insertSize+lambda*var && tempDD>insertSize-lambda*var){
                if(token!=true){
                    positionMatchNumber++;
                    matchIndex[j] = 1;
                    token = true;
                    break;
                }
            }
            temp1 = tempReadMate;
            tempReadMate = tempReadMate->next;
            delete [] temp1->readMate;
            delete temp1;       
        }
        if(token!=true){
            matchIndex[j] = 0;
        } 
        
        i++;
        if(i >= length){
            if(positionMatchNumber<=1){
                
                if(tempLengthLeft == insertSize + lambda*var){
                    delete [] tempContigLeft;
                }
                if(tempLengthRight == insertSize + lambda*var){
                    delete [] tempContigRight;
                }
                delete [] temp;
                delete [] tempRC;
                delete [] gapIndex;
                delete [] nullIndex;
                delete [] matchIndex;
                
                return 0;
            }
            double a = 0;
            if(length - nullNumber - positionNonExist>1){
                a = (double)(positionMatchNumber)/(double)(length - nullNumber - positionNonExist);
            }else{
                a = 0.2;
            }
            
            if(a<0.2){
                
                if(tempLengthLeft == insertSize + lambda*var){
                    delete [] tempContigLeft;
                }
                if(tempLengthRight == insertSize + lambda*var){
                    delete [] tempContigRight;
                }
                delete [] temp;
                delete [] tempRC;
                delete [] gapIndex;
                delete [] nullIndex;
                delete [] matchIndex;
                
                return 0;
            }
            
            nullNumber = 0;
            positionMatchNumber = 0;
            positionNonExist = 0;
            i = 0;
        }
                       
    }

    i = 0;
    nullNumber = 0;
    positionMatchNumber = 0;
    positionNonExist = 0;
    
    for(j=overlap - kmerLength;j<tempLengthLeft - readLength - length && j<insertSize - lambda*var - readLength && tempLengthRight>insertSize;j++){
        index1 = 0;
        matchIndex[j] = 0;
        SubContig(temp, tempContigLeft, tempLengthLeft - j - readLength, tempLengthLeft - j);

        ReadMate * tempReadMate = SearchRightReadMate(temp, readSet+readSetIndex);
        ReverseComplement(temp, tempRC);
        long int a = SearchRead(tempRC, readSet + readSetIndex);
        
        if(tempReadMate==NULL && a ==-1){
            positionNonExist++;
            gapIndex[j] = 1;
            index1 = 1;
        }else{
            gapIndex[j] = 0;
        }
        
        if(tempReadMate!=NULL){
            temp1 = tempReadMate;
            tempReadMate = tempReadMate->next;
            delete temp1;
        }   
        
        if(tempReadMate==NULL&&index1!=1){
            nullNumber++;
            nullIndex[j] = 1;
        }else{
            nullIndex[j] = 0;
        }
        
        bool token = false;

        while(tempReadMate!=NULL){
            index1 = KMPIndexOfContigOfMisMatch(tempContigRight,tempReadMate->readMate);
            long int tempDD = index1 + 2*readLength + j - overlap;
            if(index1!=-1 && tempDD<insertSize+lambda*var && tempDD>insertSize-lambda*var){
                if(token!=true){
                    positionMatchNumber++;
                    matchIndex[j] = 1;
                    token = true;
                    break;
                }
            }
            temp1 = tempReadMate;
            tempReadMate = tempReadMate->next;
            delete [] temp1->readMate;
            delete temp1;       
        } 
        
        if(token!=true){
            matchIndex[j] = 0;
        }
        
        i++;
        if(i >= length){
            
            if(positionMatchNumber<=1){
            
                
                if(tempLengthLeft == insertSize + lambda*var){
                    delete [] tempContigLeft;
                }
                if(tempLengthRight == insertSize + lambda*var){
                    delete [] tempContigRight;
                }
                delete [] temp;
                delete [] tempRC;
                delete [] gapIndex;
                delete [] nullIndex;
                delete [] matchIndex;
                
                
                return 0;
            }
            
            double a = 0;
            if(length - nullNumber - positionNonExist>1){
                a = (double)(positionMatchNumber)/(double)(length - nullNumber - positionNonExist);
            }else{
                a = 0.2;
            }
            
            if(a<0.2){
                
                
                if(tempLengthLeft == insertSize + lambda*var){
                    delete [] tempContigLeft;
                }
                if(tempLengthRight == insertSize + lambda*var){
                    delete [] tempContigRight;
                }
                delete [] temp;
                delete [] tempRC;
                delete [] gapIndex;
                delete [] nullIndex;
                delete [] matchIndex;
                
                return 0;
            }
            
            
            nullNumber = 0;
            positionMatchNumber = 0;
            positionNonExist = 0;
            i = 0;
        }
                       
    }


    delete [] temp;
    delete [] tempRC;
    delete [] gapIndex;
    delete [] matchIndex;
    delete [] nullIndex;
    if(tempLengthRight==insertSize+ lambda*var){
        delete [] tempContigRight;
    }
    if(tempLengthLeft==insertSize+ lambda*var){
        delete [] tempContigLeft;
    }
      
    return 1;
    
}


int MergeSubContig(ContigSet * contigSet, long int kmerLength){
    long int i = 0;
    long int j = 0;
    int count = 0;
    ContigSet * temp;
    
    ContigSet * first = contigSet;

    while(contigSet!=NULL){
        if(contigSet->contig==NULL){
            contigSet = contigSet->next;
            i++;
            continue;
        }
        temp = first;
        j = 0;
        while(temp!=NULL){
            if(temp->contig==NULL || i==j){
                j++;
                temp = temp->next;
                continue;
            }
            int len = strlen(temp->contig);
            long int * next = new long int[len+1-kmerLength];
            char * pattern = new char[len+1-kmerLength];
            SubContig(pattern, temp->contig, kmerLength, len);
            long int p = KMPIndexOfContig(contigSet->contig,pattern,next);
            if(p!=-1){
                
                delete []temp->contig;
                temp->contig = NULL;
                delete []next;
                next = NULL;
                delete [] pattern;
                pattern = NULL;
                j++;
                temp = temp->next;
                continue;
            }
            
            SubContig(pattern, temp->contig, 0, len-kmerLength);
            p = KMPIndexOfContig(contigSet->contig,pattern,next);
            if(p!=-1){
                
                delete []temp->contig;
                temp->contig = NULL;
                delete []next;
                next = NULL;
                delete [] pattern;
                pattern = NULL;
                j++;
                temp = temp->next;
                continue;
            }
                       
      
            char * tempReverseContig = new char[len+1];
            ReverseComplement(temp->contig,tempReverseContig);
            SubContig(pattern, tempReverseContig, kmerLength, len);
            p = KMPIndexOfContig(contigSet->contig,pattern,next);
            if(p!=-1){
                
                delete []temp->contig;
                temp->contig = NULL;
                delete []next;
                next = NULL;
                delete []tempReverseContig;
                tempReverseContig = NULL;
                delete [] pattern;
                pattern = NULL;
                j++;
                temp = temp->next;
                continue;
            }
            
            SubContig(pattern, tempReverseContig, 0, len-kmerLength);
            p = KMPIndexOfContig(contigSet->contig,pattern,next);
            if(p!=-1){
                
                delete []temp->contig;
                temp->contig = NULL;
                
            }
            
            
            delete []tempReverseContig;
            delete []next;
            next = NULL; 
            delete [] pattern;
            pattern = NULL;
            temp = temp->next; 
            j++;          
        }
        i++;
        contigSet = contigSet->next;
    }
    return 1;
}


char * LengthOfOverlapBetweenContig(char * left, char * right, ReadSet * readSet, int readSetCount, int kmerLength){
    long int i = 0;
    long int j = 0;
    long int t = 0;
    long int index = 0;
    long int leftLength = (long int)strlen(left);
    long int rightLength = (long int)strlen(right);
    
    
    for(i=0;i<leftLength;i++){
        index = 0;
        for(j=0;j<rightLength&&i+j<leftLength;j++){
            if(left[i+j]!=right[j]){
                index = 1;
                break;
            }
        }
        if(index==0){
            break;
        }
    }
    
    if((i+j)==leftLength&&j>kmerLength){ 
        for(int p = 0;p<readSetCount;p++){
            
            if(readSet[p].insertSize - lambda*readSet[p].var<j||leftLength<readSet[p].insertSize- 2*readSet[p].readLength||rightLength<readSet[p].insertSize- 2*readSet[p].readLength){
                continue;
            }
            double temp = WeightContigMerge(left, right, readSet, p, kmerLength, j, kmerLength);
            
            if(temp!=1){
                return NULL;
            }    
        } 
        char * contig = new char[leftLength + rightLength - j + 1];
        AppendRight(contig,left,right,j+1);
        return contig;
    }
    return NULL;
}

ContigSet * ContigMerge(ContigSet * contigSet, DBGraph * deBruijnGraph, int kmerLength, ReadSet * readSet, int readSetCount){
    long int i = 0;
    long int j = 0;
    int count = 0;

    MergeSubContig(contigSet, kmerLength);
    
    ContigSet * temp;
    ContigSet * first = contigSet;
    ContigSet * previous = NULL;
    
    while(contigSet!=NULL){
        
        if(contigSet->contig==NULL){
            contigSet = contigSet->next;
            i ++;
            continue;
        }
        temp = first;
        j = 0;
        while(temp!=NULL){
            if(temp->contig==NULL || i==j){
                temp = temp->next;
                j++;
                continue;
            }
            char * tempContig = NULL;
            
            tempContig = LengthOfOverlapBetweenContig(temp->contig,contigSet->contig,readSet,readSetCount,kmerLength);
            
            if(tempContig!=NULL){
                
                delete []contigSet->contig;
                contigSet->contig = tempContig;
                delete []temp->contig;
                temp->contig = NULL;
                temp = temp->next;
                j++;
                continue;
            }
            char * tempReverseLeft = new char[strlen(temp->contig)+1];
            ReverseComplement(temp->contig,tempReverseLeft);
            
            tempContig = LengthOfOverlapBetweenContig(tempReverseLeft,contigSet->contig,readSet,readSetCount,kmerLength);
            
            if(tempContig!=NULL){
                
                delete []contigSet->contig;
                contigSet->contig = tempContig;
                delete []temp->contig;
                temp->contig = NULL;
                delete []tempReverseLeft;
                tempReverseLeft = NULL;
                j++;
                temp = temp->next;
                continue;
                
            }
            delete []tempReverseLeft;
            tempReverseLeft = NULL;
            char * tempReverseRight = new char[strlen(contigSet->contig)+1];
            ReverseComplement(contigSet->contig,tempReverseRight);
            
            tempContig = LengthOfOverlapBetweenContig(temp->contig,tempReverseRight,readSet,readSetCount,kmerLength);
            
            if(tempContig!=NULL){
                
                delete []contigSet->contig;
                contigSet->contig = tempContig;
                delete []temp->contig;
                temp->contig = NULL;
                delete []tempReverseRight;
                tempReverseRight = NULL;
                temp = temp->next;
                j++;
                continue;
                
            }
            temp = temp->next;
            j++;
            delete []tempReverseRight;
            tempReverseRight = NULL;
        }
        i++;
        contigSet = contigSet->next;
    }
    
    
    
    temp = first;
    while(first!=NULL){
        if(first->contig==NULL){
            if(previous != NULL){
                previous->next = first->next;
                delete first;
            }else{
                ContigSet * temp11 = first;
                first = first->next;
                temp11->next = NULL;
                temp = first;
                delete temp11;
                continue;
            }           
        }else{
            previous = first;
        }
        first = previous->next;
    }
    return temp;
}










#endif
