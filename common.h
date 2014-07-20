#ifndef COMMON_H_INCLUDED 
#define COMMON_H_INCLUDED 

#include <cstring>
#include <iostream>
#include "stdio.h"
#include "stdlib.h"
#include "bitarray.h"
#include "math.h"
#include <fstream>

using namespace std;

static unsigned int globalKmerLength = 21;
static unsigned int lambda = 3;
static long int maxDeep = 13;
static long int allkmernum = 0;
static long int extendCutOff = 300;
static long int dfsExtendLength = 21;
static long int storeContigExtendLength = 10000;
static long int globalFillingPathNum = 0;
static long int minKmerLength = 17;
static long int maxKmerLength = 25; 
static long int fillRepeateNumberMax = 4;
static int cut = 21;
static int minKmerFrequency = 1;
static int storeLengthPathNum = 0;
static int wrongFillingPathNum = 0;

unsigned long int Hash(char * str, unsigned int len, unsigned long int max)  
{  
   unsigned int hash = 0;  
   unsigned int i = 0;  
  
   for(i = 0; i < len; str++, i++) {  
      hash = (*str) + (hash << 6) + (hash << 16) - hash;  
   }  
  
   return hash % max;  
} 

void KMPGetNext(char * pattern, long int * next){
    next[0] = -1;
    long int k = -1;
    long int j = 0;
    while(pattern[j]!= '\0'){
        if(k!=-1  && pattern[k]!= pattern[j]){
            k=next[k];
        }        
        ++j;
        ++k;
        if(pattern[k]==pattern[j]){
            next[j]=next[k];
        }else{
            next[j]=k;
        }
    }
}

int KMPIndexOfContig(char * contig, char * pattern, long int * next,int a){
    long int i = 0;
    long int j = 0;
    long int index = 0;
    long int len1 = strlen(contig);
    long int len2 = strlen(pattern);
    if(len1<len2){
        return -1;
    }
    while(i<len1&&j<len2){
        if(contig[i]==pattern[j]){
            i++;
            j++;
        }else{
            index += j-next[j];
            if(next[j]!=-1){
               j=next[j];
            }else{
                j=0;
                ++i;
            }
        }
    }
    if(pattern[j]=='\0'){
       return index;
    }else{
       return -1;  
    }
}

int KMPIndexOfContig(char * contig, char * pattern, long int * next){
    long int i = 0;
    long int j = 0;
    long int index = 0;
    long int len1 = strlen(contig);
    long int len2 = strlen(pattern);
    if(len1<len2){
        return -1;
    }
    KMPGetNext(pattern, next);
    while(i<len1&&j<len2){
        if(contig[i]==pattern[j]){
            i++;
            j++;
        }else{
            index += j-next[j];
            if(next[j]!=-1){
               j=next[j];
            }else{
                j=0;
                ++i;
            }
        }
    }
    if(pattern[j]=='\0'){
       return index;
    }else{
       return -1;  
    }
} 

int KMPIndexOfContigOfMisMatch(char * contig, char * pattern){
    long int i = 0;
    long int j = 0;
    long int index = 0;
    long int token = 0;
    long int p = 0;
    long int len1 = strlen(contig);
    long int len2 = strlen(pattern);
    if(len1<len2){
        cout<<"IndexOfContig Wrong!"<<endl;
        cout<<pattern<<"--"<<contig<<endl;
        //return -1;
        exit(0);
    }
    for(i=0;i<=len1-len2;i++){
        p = i;
        token = 0;
        index = 0;
        for(j=0;j<len2;j++){
            if(contig[p]!=pattern[j]){
                if(token == 0){
                    token = 1;
                }else{
                    index = -1;
                    break;
                }
            }
            p++;
        }
        if(index==0){
            return i;
        }
    }
    return -1;
}

int KMPIndexOfContigOfMisMatch(char * contig, long int start, long int end, char * pattern){
    long int i = 0;
    long int j = 0;
    long int index = 0;
    long int token = 0;
    long int p = 0;
    //long int len1 = strlen(contig);
    long int len2 = strlen(pattern);
    if(end - start < len2){
        cout<<"IndexOfContig Wrong!"<<endl;
        cout<<pattern<<"--"<<contig<<endl;
        return -1;
        //exit(0);
    }
    for(i=start;i<=end - len2;i++){
        p = i;
        token = 0;
        index = 0;
        for(j=0;j<len2;j++){
            if(contig[p]!=pattern[j]){
                if(token == 0){
                    token = 1;
                }else{
                    index = -1;
                    break;
                }
            }
            p++;
        }
        if(index==0){
            return i;
        }
    }
    return -1;
}

int KMPIndexOfContigOfMisMatch(char * contig, char * pattern, long int start){
    long int i = 0;
    long int j = 0;
    long int index = 0;
    long int token = 0;
    long int p = 0;
    long int len1 = strlen(contig);
    long int len2 = strlen(pattern);
    if(len1<len2){
        cout<<"IndexOfContig Wrong!"<<endl;
        cout<<pattern<<"--"<<contig<<endl;
        return -1;
        //exit(0);
    }
    for(i=start;i<=len1-len2;i++){
        p = i;
        token = 0;
        index = 0;
        for(j=0;j<len2;j++){
            if(contig[p]!=pattern[j]){
                if(token == 0){
                    token = 1;
                }else{
                    index = -1;
                    break;
                }
            }
            p++;
        }
        if(index==0){
            return i;
        }
    }
    return -1;
}



bool ReverseComplement(char * temp1, char * temp2){
    long int len = strlen(temp1);
    long int i = 0;
    long int j = 0;
    for(i=0;i<len;i++){
        if(temp1[i]=='A'){
            temp2[len-1-i]='T';
        }else if(temp1[i]=='T'){
            temp2[len-1-i]='A';
        }else if(temp1[i]=='G'){
            temp2[len-1-i]='C';
        }else if(temp1[i]=='C'){
            temp2[len-1-i]='G';
        }else{
             return false;
        }
    }
    temp2[len]='\0';
    return true;
}

bool ReverseComplementMatch(char * temp1, char * temp2){
    long int len = strlen(temp1);
    if(len!=strlen(temp2))
        return false;
    long int i = 0;
    long int j = 0;
    for(i=0;i<len;i++){
        if((temp1[i]=='A'&&temp2[len-1-i]=='T')||(temp1[i]=='T'&&temp2[len-1-i]=='A')
            ||(temp1[i]=='G'&&temp2[len-1-i]=='C')||(temp1[i]=='C'&&temp2[len-1-i]=='G')){
            continue;
        }else{
            return false;
        }
    }
    return true;
}


int GetReadOrKmerLength(char * address){
    int i = 0;
    int j = 0;
	char temp[1000];   	
	ifstream icin;
	icin.open(address);
	icin.getline(temp, 1000);
	icin.getline(temp, 1000);
	while(1){
        if(temp[i]=='\0') break;
        if(temp[i]=='\n') break;
        if(temp[i]=='\r') break;
        i++;
    }
    
	icin.close();
	return i;
}

void AppendRight( char * temp, char * temp1, char * temp2,long int kmerLength){
     long int i = 0;
     long int j = 0;
     long int len1 = strlen(temp1);
     long int len2 = strlen(temp2);
     for(i =0;i<len1+len2-kmerLength+1;i++){
         if(i<len1){
             temp[i] = temp1[i];
         }else{
             temp[i] = temp2[kmerLength-1+j];
             j++;
         }
     }
     temp[i] = '\0';
}

void CopyContig(char *temp, char * temp1){
     long int i = 0;
     long int j = 0;
     long int len1 = strlen(temp1);
     for(i =0;i<len1;i++){
         temp[i] = temp1[i];
     }
     temp[i] = '\0';
}

void SubContig(char *temp,char *temp1,long int start,long int end){
    long int i = 0;
    for(i = 0; i<end-start;i++){
        temp[i]=temp1[start+i];
    }
    temp[i] = '\0';
}

void FillContig(char * temp, long int start, char * temp1, long int start1, long int length){
    long int i = 0;
    long int j = 0;
    for(i = 0; i<length; i++){
        temp[start + i] = temp1[start1 + i];
    }
}

int CompareStr(char * temp, char * temp1, int len){
    if(temp==NULL&&temp1==NULL) return 0;
    if(temp==NULL&&temp1!=NULL) return -1;
    if(temp1==NULL&&temp!=NULL) return 1;
    
    int i = 0;
    for(i = 0; i<len; i++){
        if(GetBit(temp, i)>GetBit(temp1, i)){
            return 1;
        }
        if(GetBit(temp, i)<GetBit(temp1, i)){
            return -1;
        }
    }
    return 0;
}


double ComputeSD(long * data,long int insertSize,long int var){
    long int i = 0;
    long int j = 0;
    double avg = 0;
    double sd = 0;
    
    
    long int outNumber[7];
    for(i = 0; i<7; i++){
        outNumber[i] = 0;
    }
    i = 0;

    while(data[i]!=0){
        
        long int index = (data[i]-insertSize)/var;
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
        
        i++;
    }

    for(j = 0; j<7; j++){
        //ocout<<"--"<<outNumber[j];
    }

    double allNumber = i;
    
    avg = avg + pow(outNumber[0] - allNumber*0.02275,2)/(allNumber*0.02275);
    avg = avg + pow(outNumber[1] - allNumber*0.13595,2)/(allNumber*0.13595);
    avg = avg + pow(outNumber[2] - allNumber*0.3413,2)/(allNumber*0.3413);
    avg = avg + pow(outNumber[3] - allNumber*0.3413,2)/(allNumber*0.3413);
    avg = avg + pow(outNumber[4] - allNumber*0.13595,2)/(allNumber*0.13595);
    avg = avg + pow(outNumber[5] - allNumber*0.02275,2)/(allNumber*0.02275);
    
    double avg1 = 0;
    avg1 = avg1 + pow(outNumber[0] - allNumber*0.02275,2);
    avg1 = avg1 + pow(outNumber[1] - allNumber*0.13595,2);
    avg1 = avg1 + pow(outNumber[2] - allNumber*0.3413,2);
    avg1 = avg1 + pow(outNumber[3] - allNumber*0.3413,2);
    avg1 = avg1 + pow(outNumber[4] - allNumber*0.13595,2);
    avg1 = avg1 + pow(outNumber[5] - allNumber*0.02275,2);
    
    double all = 0;
    double allNumber1 = 0;
    for(j = 0; j<6; j++){
        all = all + pow(outNumber[j],2);
        allNumber1 = allNumber1 + outNumber[j];
    }
    
    /*
    double avgJ = a11J/6;
    double d = 0;
    double d1 = 0;
    double d2 = 0;
    for(j = 0; j<6; j++){
        d = d + pow(outNumber[j] - avgJ, 3);
        d1 = d1 + pow(outNumber[j] - avgJ, 2);
        d2 = d2 + pow(outNumber[j] - avgJ, 4);
    }
    d = d/6;
    d1 = d1/6;
    d2 = d2/6;
    double s = d/pow(d1,1.5);
    double k = d2/pow(d1,2);
    double jb = pow(s,2) + 0.25*pow(k-3,2);
    */
    
    double fre[6];
    for(j = 0; j<6; j++){
        fre[j] = (double)outNumber[j]/allNumber1;
    }
    double maxF = 0;
    if(maxF < fabs(fre[0] - 0.02275)){
        maxF = fabs(fre[0] - 0.02275);
    }
    if(maxF < fabs(fre[1] - 0.13595)){
        maxF = fabs(fre[1] - 0.13595);
    }
    if(maxF < fabs(fre[2] - 0.3413)){
        maxF = fabs(fre[2] - 0.3413);
    }
    if(maxF < fabs(fre[3] - 0.3413)){
        maxF = fabs(fre[3] - 0.3413);
    }
    if(maxF < fabs(fre[4] - 0.13595)){
        maxF = fabs(fre[4] - 0.13595);
    }
    if(maxF < fabs(fre[5] - 0.02275)){
        maxF = fabs(fre[5] - 0.02275);
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
       
    
    double ks = 1;
    if(maxF>1.63/sqrt(allNumber1)){
        ks = 0.01;
    }else if(maxF>1.37/sqrt(allNumber1)){
        ks = 0.1;
    }else if(maxF>1.36/sqrt(allNumber1)){
        ks = 0.15;
    }else if(maxF>1.22/sqrt(allNumber1)){
        ks = 0.3;
    }else if(maxF>1.07/sqrt(allNumber1)){
        ks = 0.5;
    }else if(maxF>0.87/sqrt(allNumber1)){
        ks = 0.7;
    }
    
    //ocout<<"avg:"<<avg<<"--maxF:"<<maxF<<"--ks:"<<ks<<"--"<<avg2/(all*avg3)<<"--";
    return avg2/(all*avg3); 
    //return 1-maxF;
    //return ks;
    
    /*
    avg = avg/i;
    i = 0;
    while(data[i]!=0){
        sd = sd + pow(data[i]-insertSize,2);
        i++;
    }
    sd = sd/i;
    return sqrt(sd);
    */
}


/*
double Possion(int length, double p, int n){
    
    int i = 0;
    int j = 0;
    double e = 2.71828;
    long int t = 1;
    for(i = 1; i<=n; i++){
        t = t*i;
    }
    
    return (pow(e,-length*p)*pow(length*p,n))/(double)t;
    
    
}
*/

#endif
