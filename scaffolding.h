#ifndef SCAFFOLDING_H_INCLUDED 
#define SCAFFOLDING_H_INCLUDED 

#include <cstring>
#include <iostream>
#include "stdio.h"
#include "stdlib.h"
#include <math.h>

#include "readset.h"
#include "graph.h"

double * Orientation(char * contigLeft, char * contigRight, ReadSet * readSet, long int readSetIndex, long int start, long int length, int index){
    long int i = 0;
    long int j = 0;
    long int readLength = readSet[readSetIndex].readLength;
    long int insertSize = readSet[readSetIndex].insertSize;
    long int var = readSet[readSetIndex].var;   
    
    long int nullNumber = 0;
    long int positionMatchNumber = 0;
    long int allMatchNumber = 0;
    long int positionNonExist = 0;
    
    double score = 0;
    double avgScore = 0;
    long int gapDistance = 0;
    long int avgGapDistance = 0;
    long int scoreNumber = 0;
    
    double * result = new double[4];
    for(i = 0; i<4;i++){
        result[i] = 0;
    }

    char * tempContigLeft = contigLeft;
    char * tempContigRight = contigRight;
    long int tempLeftLength = strlen(contigLeft);
    long int tempRightLength = strlen(contigRight);
    
    if(tempLeftLength >= insertSize + lambda*var){     
        tempContigLeft = new char[insertSize + lambda*var + 1];
        if(index!=1){
            SubContig(tempContigLeft,contigLeft,tempLeftLength-insertSize-lambda*var,tempLeftLength);
        }else{
            char * contigRC = new char[insertSize + lambda*var + 1];
            SubContig(contigRC,contigLeft,0, insertSize+lambda*var);
            ReverseComplement(contigRC, tempContigLeft);
            delete [] contigRC;
        }
        tempLeftLength = insertSize+lambda*var;
    }else if(index==1){
        tempContigLeft = new char[tempLeftLength + 1];
        ReverseComplement(contigLeft, tempContigLeft);
    }
    
    if(tempRightLength >= insertSize + lambda*var){
        tempContigRight = new char[insertSize + lambda*var + 1];
        if(index!=2){
            SubContig(tempContigRight,contigRight,0,insertSize + lambda*var);
        }else{
            char * contigRC = new char[insertSize + lambda*var + 1];
            SubContig(contigRC,contigRight,tempRightLength-insertSize-lambda*var,tempRightLength);
            ReverseComplement(contigRC, tempContigRight);
            delete [] contigRC;
        }  
        tempRightLength = insertSize+lambda*var;
    }else if(index == 2){
        tempContigRight = new char[tempRightLength + 1];
        ReverseComplement(contigRight, tempContigRight);
    }
    
    char * temp = new char[readLength + 1];
    char * tempRC = new char[readLength + 1];
    ReadMate * temp1;   
    for(i=start; i<tempRightLength-length-readLength; i++){
        
        if(avgGapDistance!=0&&labs(insertSize - avgGapDistance) + i + 2*readLength + length > insertSize -lambda*var ){
            break;
        }
        bool token = false; 
        SubContig(temp, tempContigRight, i, i + readLength);
        ReadMate * tempReadMate = SearchLeftReadMate(temp, readSet+readSetIndex);
        long int a = SearchRead(temp, readSet + readSetIndex);
            
        if(tempReadMate==NULL && a ==-1){
            positionNonExist++;
            token = true;
        }        
        if(tempReadMate!=NULL){
            temp1 = tempReadMate;
            tempReadMate = tempReadMate->next;
            delete temp1;
        }         
        if(tempReadMate==NULL && token == false){
            nullNumber++;
        }      
        token = false; 
        while(tempReadMate!=NULL){
            long int index1 = KMPIndexOfContigOfMisMatch(tempContigLeft,tempReadMate->readMate);
            long int tempDD = i + tempLeftLength - index1 + readLength;
            //cout<<"index1--"<<index1<<"--";
            if(index1!=-1 && tempDD>0 && tempDD<insertSize+lambda*var){
                allMatchNumber++;
                if(token!=true){
                    positionMatchNumber++;
                    token = true;
                }
                gapDistance = gapDistance + tempDD;
            }
            temp1 = tempReadMate;
            tempReadMate = tempReadMate->next;
            delete [] temp1->readMate;
            delete temp1;
        }   
        j++;
        if(j==length){   
            if(positionMatchNumber<=2){
                delete [] result;
                
                if(tempLeftLength == insertSize + lambda*var || index==1){
                    delete [] tempContigLeft;
                }
                if(tempRightLength == insertSize + lambda*var||index==2){
                    delete [] tempContigRight;
                }
                
                
                //cout<<"null--"<<positionMatchNumber<<"--"<<scoreNumber<<"--";
                return NULL;
            }
            gapDistance = gapDistance/allMatchNumber;
        
            if(avgGapDistance!=0){
                avgGapDistance = (avgGapDistance + gapDistance)/2;
            }else{
                avgGapDistance = gapDistance;
            }            
            score = (double)(positionMatchNumber)/(double)(length - nullNumber - positionNonExist);   
            if(score < 0.2 || labs(avgGapDistance - insertSize) > insertSize - lambda*var){
                delete [] result;
                
                if(tempLeftLength == insertSize + lambda*var || index==1){
                    delete [] tempContigLeft;
                }
                if(tempRightLength == insertSize + lambda*var||index==2){
                    delete [] tempContigRight;
                }
                
                //cout<<"nulltt--"<<score<<"--"<<scoreNumber<<endl;
                return NULL;
            }
            avgScore = score + avgScore;            
            nullNumber = 0;
            positionMatchNumber = 0;
            positionNonExist = 0;
            allMatchNumber = 0;
            gapDistance = 0;
            j = 0;
            scoreNumber++;
            i = i + var;
        }
    }
    
    if(scoreNumber!=0 && avgScore/scoreNumber>0.7){
        //cout<<"scoreNumber:"<<scoreNumber<<"--";
        result[0] = avgScore/scoreNumber;
        result[1] = labs(avgGapDistance - insertSize);
    }else{
        return NULL;
    }
    
    nullNumber = 0;
    positionMatchNumber = 0;
    allMatchNumber = 0;
    positionNonExist = 0;
    
    j = 0;
    score = 0;
    avgScore = 0;
    avgGapDistance = 0;
    gapDistance = 0;
    scoreNumber = 0;
    
    for(i=start; i<tempLeftLength-length-readLength; i++){
    
        if(avgGapDistance!=0 && labs(insertSize - avgGapDistance) + i + 2*readLength + length > insertSize - lambda*var){
            break;
        }
        bool token = false; 
        SubContig(temp, tempContigLeft, tempLeftLength - i - readLength, tempLeftLength - i);
        ReadMate * tempReadMate = SearchRightReadMate(temp, readSet+readSetIndex);
        //cout<<"aa4--"<<tempLeftLength<<"--"<<length<<"--"<<readLength<<"--"<<i<<"--"<<temp<<endl;
        ReverseComplement(temp,tempRC);
        long int a = SearchRead(tempRC, readSet + readSetIndex);
            
        if(tempReadMate==NULL && a ==-1){
            positionNonExist++;
            token = true;
        }        
        if(tempReadMate!=NULL){
            temp1 = tempReadMate;
            tempReadMate = tempReadMate->next;
            delete temp1;
        }         
        if(tempReadMate==NULL && token == false){
            nullNumber++;
        }      
        token = false; 
        while(tempReadMate!=NULL){
            long int index1 = KMPIndexOfContigOfMisMatch(tempContigRight,tempReadMate->readMate);
            long int tempDD = i + index1 + 2*readLength;
            if(index1!=-1 && tempDD>0 && tempDD<insertSize+lambda*var){
                allMatchNumber++;
                if(token!=true){
                    positionMatchNumber++;
                    token = true;
                }
                gapDistance = gapDistance + tempDD;
            }
            temp1 = tempReadMate;
            tempReadMate = tempReadMate->next;
            delete [] temp1->readMate;
            delete temp1;
        } 
        
        j++;
        if(j==length){         
            if(positionMatchNumber<=2){
                delete [] result;
                
                if(tempLeftLength == insertSize + lambda*var || index==1){
                    delete [] tempContigLeft;
                }
                if(tempRightLength == insertSize + lambda*var||index==2){
                    delete [] tempContigRight;
                }
                
                //cout<<"null--"<<positionMatchNumber<<"--"<<scoreNumber<<endl;
                return NULL;
            }
            
            gapDistance = gapDistance/allMatchNumber;
        
            if(avgGapDistance!=0){
                avgGapDistance = (avgGapDistance + gapDistance)/2;
            }else{
                avgGapDistance = gapDistance;
            }
            
            score = (double)(positionMatchNumber)/(double)(length - nullNumber - positionNonExist);
            if(score < 0.2 || labs(avgGapDistance - insertSize) > insertSize - lambda*var){
                delete [] result;
                
                if(tempLeftLength == insertSize + lambda*var || index==1){
                    delete [] tempContigLeft;
                }
                if(tempRightLength == insertSize + lambda*var||index==2){
                    delete [] tempContigRight;
                }
                
                //cout<<"nulltt--"<<score<<"--"<<scoreNumber<<endl;
                return NULL;
            }
            avgScore = avgScore + score;            
            nullNumber = 0;
            positionMatchNumber = 0;
            positionNonExist = 0;
            allMatchNumber = 0;
            gapDistance = 0;
            j=0;
            scoreNumber++;
            i = i + var;
        }
    }
    
    if(scoreNumber!=0 && avgScore/scoreNumber>0.7){
        //cout<<"scoreNumber:"<<scoreNumber<<"--";
        result[2] = avgScore/scoreNumber;
        result[3] = labs(avgGapDistance - insertSize);
    }else{
        return NULL;
    }

    delete [] temp;
    delete [] tempRC;
    if(tempLeftLength == insertSize + lambda*var || index==1){
        delete [] tempContigLeft;
    }
    if(tempRightLength == insertSize + lambda*var||index==2){
        delete [] tempContigRight;
    }

    /*
    cout<<"Result:";
    for(i = 0; i< 4; i++){
        cout<<result[i]<<"--";
    }
    cout<<endl;
    */

    return result;
    

}

double * OrientationBetweenContig(char * contigLeft, char * contigRight, ReadSet * readSet, int setNumber, long int start, long int length, int index){
    
    long int i = 0;
    long int j = 0;
    long int p = 0;
    double result[setNumber][4];
    for(i=0;i<setNumber;i++){
        for(j=0;j<4;j++){
            result[i][j] = 0;
        }
    }
    for(i = 0; i < setNumber; i++){
        double * tempResult = Orientation(contigLeft,contigRight,readSet,i,start,length, index);
        if(tempResult == NULL){
            continue;
        }
        for(j=0;j<4;j++){
            result[i][j] = tempResult[j];
        }
        delete [] tempResult; 
    }
    for(i=0;i<setNumber;i++){
        if(result[i][0]==0){
            continue;
        }
        for(j = i+1; j<setNumber; j++){
            if(result[j][0]==0){
                return NULL;
            }
        }

        double * temp = new double[6];
        temp[0] = result[i][0];
        temp[1] = i;
        temp[2] = (result[i][1]+result[i][3])/2;
        temp[3] = result[i][2];
        temp[4] = i;
        temp[5] = temp[2];

        return temp;
    }
    return NULL;
}

#pragma pack(2)
typedef struct ScaffoldingContigSetP{
    
    ContigSet * contigSet;
    long int count;
    ReadSet * readSet;
    int setNumber;
    long int start;
    long int length;
    int threadIndex;
    int totalThreadNumber;
    double (* linkRight)[5];
    double (* linkLeft)[5];
    
    
}ScaffoldingContigSetP;
#pragma pack ()

void * ScaffoldingContigSetThread(void * arg){
       
    long int i = 0;
    long int j = 0;
    
    ScaffoldingContigSetP * scaffoldingContigSetP = (ScaffoldingContigSetP *)arg;
    
    ContigSet * contigSet = scaffoldingContigSetP->contigSet;
    ReadSet * readSet = scaffoldingContigSetP->readSet;
    int setNumber = scaffoldingContigSetP->setNumber;
    long int start = scaffoldingContigSetP->start;
    long int length = scaffoldingContigSetP->length;
    int threadIndex = scaffoldingContigSetP->threadIndex;
    int totalThreadNumber = scaffoldingContigSetP->totalThreadNumber;
    double (*linkRight)[5] = scaffoldingContigSetP->linkRight;
    double (*linkLeft)[5] = scaffoldingContigSetP->linkLeft;
    
    
    ContigSet * leftContig = contigSet;
    ContigSet * rightContig = contigSet;
    
    while(leftContig!=NULL){
        //cout<<"first:"<<i<<endl;
        if(i%totalThreadNumber != threadIndex){
            leftContig=leftContig->next;
            i++;
            continue;
        }
        j = 0;
        rightContig = contigSet;
        while(rightContig!=NULL){
            if(i==j){
                j++;
                rightContig = rightContig->next;
                continue;
            }
            //cout<<"second:"<<j<<"--"<<"len:"<<strlen(rightContig->contig);            
            double * temp = OrientationBetweenContig(leftContig->contig,rightContig->contig,readSet,setNumber,start,length,0);
            if(temp!=NULL){
                //cout<<"aa:"<<temp[0]<<"--"<<temp[1]<<"--"<<temp[2]<<"--"<<temp[3]<<"--"<<temp[4]<<"--"<<temp[5]<<endl;
            }
            double * temp1 = NULL;
            double * temp2 = NULL;
            if(temp==NULL){
                temp1 = OrientationBetweenContig(leftContig->contig,rightContig->contig,readSet,setNumber,start,length,2);
            }
            if(temp1!=NULL){
                //cout<<"bb:"<<temp1[0]<<"--"<<temp1[1]<<"--"<<temp1[2]<<"--"<<temp1[3]<<"--"<<temp1[4]<<"--"<<temp1[5]<<endl;
            }
            if(temp!=NULL&&*(linkRight[i]+1)<temp[0]){
                *(linkRight[i]+0) = j;
                *(linkRight[i]+1)=temp[0];
                *(linkRight[i]+2)=temp[1];
                *(linkRight[i]+3)=1;
                *(linkRight[i]+4)=temp[2];
                
            }
            if(temp!=NULL&&*(linkLeft[j]+1)<temp[3]){
                *(linkLeft[j]+0) = i;
                *(linkLeft[j]+1)=temp[3];
                *(linkLeft[j]+2)=temp[4];
                *(linkLeft[j]+3)=1;
                *(linkLeft[j]+4)=temp[5];
                
            }
            if(temp1!=NULL&&*(linkRight[i]+1)<temp1[0]){
                *(linkRight[i]+0) = j;
                *(linkRight[i]+1)=temp1[0];
                *(linkRight[i]+2)=temp1[1];
                *(linkRight[i]+3)=-1;
                *(linkRight[i]+4)=temp1[2];
            }
            if(temp1!=NULL&&*(linkRight[j]+1)<temp1[3]){
                *(linkRight[j]+0) = i;
                *(linkRight[j]+1)=temp1[0];
                *(linkRight[j]+2)=temp1[1];
                *(linkRight[j]+3)=-1;
                *(linkRight[j]+4)=temp1[2];
                
            }
            
            if(temp==NULL&&temp1==NULL){
                temp2 = OrientationBetweenContig(leftContig->contig,rightContig->contig,readSet,setNumber,start,length,1);
            }
            if(temp2!=NULL){
                //cout<<"cc:"<<temp2[0]<<"--"<<temp2[1]<<"--"<<temp2[2]<<"--"<<temp2[3]<<"--"<<temp2[4]<<"--"<<temp2[5]<<endl;
            }
            if(temp2!=NULL&&*(linkLeft[i]+1)<temp2[0]){
                *(linkLeft[i]+0) = j;
                *(linkLeft[i]+1)=temp2[3];
                *(linkLeft[i]+2)=temp2[4];
                *(linkLeft[i]+3)=-1;
                *(linkLeft[i]+4)=temp2[5];
            }
            if(temp2!=NULL&&*(linkLeft[j]+1)<temp2[3]){
                *(linkLeft[j]+0) = i;
                *(linkLeft[j]+1)=temp2[3];
                *(linkLeft[j]+2)=temp2[4];
                *(linkLeft[j]+3)=-1;
                *(linkLeft[j]+4)=temp2[5];
                
            }
            
            delete []temp;
            delete []temp1;
            delete []temp2;          
            
            rightContig = rightContig->next;
            j++;
        }

        leftContig = leftContig->next;
        i++; 
    }
    
    
}

#pragma pack(2)
typedef struct Scaffold{
    char * contig;
    long int gapDistance;
    struct Scaffold * next;
    Scaffold(){
        contig = NULL;
        gapDistance = 0;
        next = NULL;
    }
}Scaffold;
#pragma pack ()

#pragma pack(2)
typedef struct ScaffoldSet{
    Scaffold * scaffold;
    ScaffoldSet(){
        scaffold = NULL;
    }
}ScaffoldSet;
#pragma pack ()

#pragma pack(2)
typedef struct ScaffoldSetHead{
    ScaffoldSet * scaffoldSet;
    long int scaffoldSetNumber;
    long int gapNumber;
    ScaffoldSetHead(){
        scaffoldSet = NULL;
        scaffoldSetNumber = 0;
        gapNumber = 0;
    }
}ScaffoldSetHead;
#pragma pack ()

void WriteScaffoldSet(ScaffoldSetHead * scaffoldSetHead, char * temp){
    
    long int i = 0;
    long int j = 0;
    long int t = 0;
    long int p = 0;
    
    ofstream ocout;
    ocout.open(temp);
    
    long int rowLength = 500;
    char * tempContig;
    
    for(t = 0; t < scaffoldSetHead->scaffoldSetNumber; t++){
        Scaffold * tempScaffold = scaffoldSetHead->scaffoldSet[t].scaffold;
        while(tempScaffold!=NULL){
            long int len = strlen(tempScaffold->contig);
            j = (long int)(len/rowLength);
            for(i = 0; i<j+1;i++){
                if(i!=j){
                    tempContig = new char[rowLength+1];
                    SubContig(tempContig,tempScaffold->contig,i*rowLength,i*rowLength+rowLength);
                }else{
                    tempContig = new char[len - i*rowLength+1];
                    SubContig(tempContig,tempScaffold->contig,i*rowLength,len);
                }
                ocout<<">"<<p<<"--"<<i<<"--"<<tempScaffold->gapDistance<<endl;
                ocout<<tempContig<<endl;
                delete []tempContig;
            }
            p++; 
            tempScaffold = tempScaffold->next;         
        }
    }
    
}





void WriteScaffoldSetLong(ScaffoldSetHead * scaffoldSetHead, char * temp){
    
    long int i = 0;
    long int j = 0;
    long int t = 0;
    long int p = 0;
    
    ofstream ocout;
    ocout.open(temp);
    
    for(t = 0; t < scaffoldSetHead->scaffoldSetNumber; t++){
        Scaffold * tempScaffold = scaffoldSetHead->scaffoldSet[t].scaffold;
        while(tempScaffold!=NULL){
            ocout<<">"<<p<<"--"<<strlen(tempScaffold->contig)<<endl;
            ocout<<tempScaffold->contig<<endl;
            tempScaffold = tempScaffold->next; 
            p++;        
        }
    }    
}

void WriteScaffoldSetLongN(ScaffoldSetHead * scaffoldSetHead, char * temp){
    
    long int i = 0;
    long int j = 0;
    long int t = 0;
    long int p = 0;
    
    ofstream ocout;
    ocout.open(temp);
    
    long int rowLength = 500;
    char * tempContig;
        
    for(t = 0; t < scaffoldSetHead->scaffoldSetNumber; t++){
        Scaffold * tempScaffold = scaffoldSetHead->scaffoldSet[t].scaffold;
        ocout<<">"<<p<<endl;
        while(tempScaffold!=NULL){
            ocout<<tempScaffold->contig;
            if(tempScaffold->gapDistance>10){
                char * tempN = new char[tempScaffold->gapDistance + 1];
                int ss = 0;
                for(ss = 0; ss < tempScaffold->gapDistance; ss++){
                    tempN[ss] = 'N';
                }
                tempN[ss] = '\0';
                ocout<<tempN; 
            }else{
                char * tempN = new char[10 + 1];
                int ss = 0;
                for(ss = 0; ss < 10; ss++){
                    tempN[ss] = 'N';
                }
                tempN[ss] = '\0';
                ocout<<tempN; 
            }            
            tempScaffold = tempScaffold->next;         
        }
        ocout<<endl;
        p++;
    }   
}

int ComputeResult(ContigSet * contigSet, long int length){
    long int i = 0;
    long int j = 0;
    long int count = 0;
    long int temp = 0;
    long int max = 0;
    long int n50 = 0;
    long int n90 = 0;
    long int totalLength = 0;
    ContigSet * tempContigSet = contigSet;
    while(tempContigSet!=NULL){
        tempContigSet = tempContigSet->next;
        count++;
    }
    if(count==0){
        return 0;
    }
    long int * contigSetLength = new long int[count];
    long int t = 0;
    
    tempContigSet = contigSet;
    while(i<count){
        contigSetLength[i] = strlen(tempContigSet->contig);
        if(contigSetLength[i]<1000){
            t++;
        }else{
            totalLength = totalLength + contigSetLength[i];
        }
        i++;
        tempContigSet = tempContigSet->next;
    }
    //cout<<t<<endl;
    for(i = 0;i<count-1;i++){
        for(j = i+1;j<count;j++){
            if(contigSetLength[i]<contigSetLength[j]){
                temp = contigSetLength[j];
                contigSetLength[j] = contigSetLength[i];
                contigSetLength[i] = temp;
            }
        }
    }
    max = contigSetLength[0];
    temp = 0;
    for(i = 0;i<count - t;i++){
        temp = temp + contigSetLength[i];
        if(temp>totalLength/2){
            if(n50==0){
                n50 = contigSetLength[i];
            } 
        }
        if(temp>totalLength*0.9){
            if(n90==0){
                n90 = contigSetLength[i];
            }
        }
    }
    cout<<"num:"<<count-t<<"--maxLength:"<<max<<"--N50:"<<n50<<"--N90:"<<n90<<"--averageLength:"<<temp/(count-t)<<"--allLength:"<<totalLength<<endl;
     
}

ScaffoldSetHead * ScaffoldingContigSet(ContigSet * contigSet, ReadSet * readSet, int setNumber, long int start, long int length, int totalThreadNumber){
    
    pthread_t tid[totalThreadNumber];
    
    long int i = 0;
    long int j = 0;
    long int count = 0;
    
    ScaffoldingContigSetP * scaffoldingContigSetP = new ScaffoldingContigSetP[totalThreadNumber];
    

    ContigSet * tempContigSet = contigSet;
    
    while(tempContigSet!=NULL){
        count++;
        tempContigSet = tempContigSet->next;
    }
    
    double linkRight[count][5];
    double linkLeft[count][5];
    
    for(i=0;i<count;i++){
        linkRight[i][0] = -1;  //相邻结点 编号 
        linkRight[i][1] = -1;  //得分值 
        linkRight[i][2] = 0;   //相邻结点编号 
        linkRight[i][3] = 0;   //正方方向 
        linkRight[i][4] = 0;   //距离大小 
        linkLeft[i][0] = -1;
        linkLeft[i][1] = -1;
        linkLeft[i][2] = 0;
        linkLeft[i][3] = 0;
        linkLeft[i][4] = 0;
    }
        
    for(i = 0; i<totalThreadNumber; i++){
        scaffoldingContigSetP[i].contigSet = contigSet;
        scaffoldingContigSetP[i].count = count;
        scaffoldingContigSetP[i].readSet = readSet;
        scaffoldingContigSetP[i].setNumber = setNumber;
        scaffoldingContigSetP[i].start = start;
        scaffoldingContigSetP[i].length = length;
        scaffoldingContigSetP[i].threadIndex = i;
        scaffoldingContigSetP[i].totalThreadNumber = totalThreadNumber;
        scaffoldingContigSetP[i].linkRight = linkRight;
        scaffoldingContigSetP[i].linkLeft = linkLeft;
        
        
        if(pthread_create(&tid[i], NULL, ScaffoldingContigSetThread, (void *)&scaffoldingContigSetP[i])){
            cout<<"create thread wrong!"<<endl;
            exit(1);
        }
    }
    
    for(i = 0; i < totalThreadNumber; i++){
        pthread_join(tid[i], NULL);
    }
    
    
    
    
    for(i=0;i<count;i++){
        if(linkRight[i][0]==-1){
            continue;
        }
        
    }
   
    for(i=0;i<count;i++){
        if(linkLeft[i][0]==-1){
            continue;
        }
        
    }
    
    Scaffold * scaffold[count];
    char * tempContig[count];
    tempContigSet = contigSet;
    long int link[count][4];
    long int gapDistance[count][2];
    
    for(i = 0; i < count; i++){
        scaffold[i] = NULL;
        tempContig[i] = tempContigSet->contig;
        tempContigSet = tempContigSet->next;
        link[i][0] = -1;
        link[i][1] = 0;
        link[i][2] = -1;
        link[i][3] = 0;
        gapDistance[i][0] = 0;
        gapDistance[i][1] = 0;
    }
    
    
    
    for(i = 0; i < count; i++){
        
        if(linkRight[i][3] == 1 && linkLeft[(int)linkRight[i][0]][0] == i && linkLeft[(int)linkRight[i][0]][3] == 1){
            link[i][2] = (int)linkRight[i][0];
            link[i][3] = 1;
            gapDistance[i][1] = (int)linkRight[i][4];
        }
        if(linkRight[i][3] == -1 && linkRight[(int)linkRight[i][0]][0] == i && linkRight[(int)linkRight[i][0]][3] == -1){
            link[i][2] = (int)linkRight[i][0];
            link[i][3] = -1;
            gapDistance[i][1] = (int)linkRight[i][4];
            
        }
        if(linkLeft[i][3] == 1 && linkRight[(int)linkLeft[i][0]][0] == i && linkRight[(int)linkLeft[i][0]][3] == 1){
            link[i][0] = (int)linkLeft[i][0];
            link[i][1] = 1;
            gapDistance[i][0] = (int)linkLeft[i][4];
        }
        if(linkLeft[i][3] == -1 && linkLeft[(int)linkLeft[i][0]][0] == i && linkLeft[(int)linkLeft[i][0]][3] == -1){
            link[i][0] = (int)linkLeft[i][0];
            link[i][1] = -1;
            gapDistance[i][0] = (int)linkLeft[i][4];
        }      
    }
        
    
    for(i = 0; i < count; i++){
        
        int right = link[i][2];
        int left = link[i][0];  
        
        if(tempContig[i] == NULL){
            continue;
        }
        
        scaffold[i] = new Scaffold;
        scaffold[i]->contig = tempContig[i];
        scaffold[i]->gapDistance = gapDistance[i][1];
        tempContig[i] = NULL;
        
        
        if((right == -1 && left == -1)){         
            continue;
        }
        
        
        
        Scaffold * present = scaffold[i];
        int token = link[i][3];
        
        while(right != -1&& right!=i){
            if(token == 1){
                
                Scaffold * temp = new Scaffold;
                temp->contig = tempContig[right];
                temp->gapDistance = gapDistance[right][1];
                tempContig[right] = NULL;
                present->next = temp;
                present = temp;
                token = token * link[right][3];
                right = link[right][2];
            }
            if(token == -1 && right!=i){
                
                Scaffold * temp = new Scaffold;
                temp->contig = new char[strlen(tempContig[right])+1];
                ReverseComplement(tempContig[right], temp->contig);
                temp->gapDistance = gapDistance[right][0];
                delete [] tempContig[right];
                tempContig[right] = NULL;
                present->next = temp;
                present = temp;
                token = token * link[right][1];
                right = link[right][0];             
            }
        }
        if(right == i){
            break;
        }
        
        token = link[i][1];
        while(left != -1 && left!=i){
            if(token == 1){
                
                Scaffold * temp = new Scaffold;
                temp->contig = tempContig[left];
                temp->gapDistance = gapDistance[left][1];
                tempContig[left] = NULL;
                temp->next = scaffold[i];
                scaffold[i] = temp;
                token = token * link[left][1];
                left = link[left][0];               
            }
            if(token == -1 && left!=i){
                
                Scaffold * temp = new Scaffold;
                temp->contig = new char[strlen(tempContig[left])+1];
                ReverseComplement(tempContig[left], temp->contig);
                temp->gapDistance = gapDistance[left][0];
                delete [] tempContig[left];
                tempContig[left] = NULL;
                temp->next = scaffold[i];
                scaffold[i] = temp;
                token = token * link[left][3];
                left = link[left][2];
            }
        }
        //cout<<endl;
        
    }
    
    j = 0;
    for(i = 0; i<count; i++){
        if(scaffold[i] != NULL){
            j++;
        }
    }
    ScaffoldSet * result = new ScaffoldSet[j];
    j = 0;
    
    long int gapNumber = 0;
    for(i = 0; i<count; i++){
        if(scaffold[i] != NULL){
            result[j].scaffold = scaffold[i];
            j++;
            Scaffold * temp = scaffold[i];
            while(temp->next!=NULL){
                gapNumber++;
                temp = temp->next;
            }            
            scaffold[i] = NULL;
        }  
    }
    
    
    ScaffoldSetHead * scaffoldSetHead = new ScaffoldSetHead;
    scaffoldSetHead->scaffoldSet = result;
    scaffoldSetHead->scaffoldSetNumber = j;
    
    scaffoldSetHead->gapNumber = gapNumber;
    
    /*
    char * address0 = new char[30];
    strcpy(address0, "contig.fa");
    WriteScaffoldSet(scaffoldSetHead, address0);
    strcpy(address0, "contigLong.fa");
    WriteScaffoldSetLong(scaffoldSetHead, address0);
    */
    
    return scaffoldSetHead;
    
}


#endif
