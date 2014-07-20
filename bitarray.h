#ifndef BITS_H_INCLUDED 
#define BITS_H_INCLUDED 

#define BASE 2 
#define SHIFT 2 
#define MASK 3 
  
void SetBit(unsigned long long int * bit_array, unsigned int bit_number, char value){

    if(value == 'A')
    {
        unsigned long long int temp = 3;
        (*bit_array) &= ~(temp<<(2*bit_number)); 
    }
    if(value == 'T')
    {
        unsigned long long int temp = 3;
        (*bit_array) |= (temp<<(2*bit_number));
    }
    if(value == 'G')
    {
        unsigned long long int temp = 1;
        (*bit_array) |= (temp<<(2*bit_number));
        temp = 1;
        (*bit_array) &= ~(temp<<(2*bit_number+1));
    }
    if(value == 'C')
    {
        unsigned long long int temp = 1;
        (*bit_array) |= (temp<<(2*bit_number+1));
        temp = 1;
        (*bit_array) &= ~(temp<<(2*bit_number)); 
    }
    if(value == 'N')
    {
        unsigned long long int temp = 3;
        (*bit_array) &= ~(temp<<(2*bit_number)); 
    }

}

void SetBit(unsigned long long int * bit_array, unsigned int bit_number, int len, char * value){
    int i = 0;
    for(i = 0; i<len; i++){
        SetBit(bit_array, i, value[i]);
    }
}

char GetBit(unsigned long long int * bit_array, unsigned int bit_number){
    unsigned long long int temp = 3;
    int a = ((*bit_array)>>(2*bit_number)) & temp;
    if(a == 0){
        return 'A';
    }
    if(a == 3){
        return 'T';
    }
    if(a == 1){
        return 'G';
    }
    if(a == 2){
        return 'C';
    }
}

void GetBit(unsigned long long int * bit_array, unsigned int bit_number, int len, char * value){
    int i = 0;
    for(i = 0; i<len; i++){
        value[i] = GetBit(bit_array, i);
    }
    value[i] = '\0';
}

void SetBit(char bit_array[], unsigned int bit_number, char value) 
{ 
    unsigned int shift = 6 - BASE*(bit_number & MASK);
    
    if(value == 'A')
    {
        bit_array[bit_number >> SHIFT] &= ~ (3 << shift); 
    }
    if(value == 'T')
    {
        bit_array[bit_number >> SHIFT] |= 3 << shift; 
    }
    if(value == 'G')
    {
        bit_array[bit_number >> SHIFT] &= ~ (2 << shift);
        bit_array[bit_number >> SHIFT] |= 1 << shift;
    }
    if(value == 'C')
    {
        bit_array[bit_number >> SHIFT] &= ~ (1 << shift);
        bit_array[bit_number >> SHIFT] |= 2 << shift; 
    } 
    if(value == 'N')
    {
        bit_array[bit_number >> SHIFT] &= ~ (3 << shift); 
    }
} 

char GetBit(char bit_array[], unsigned int bit_number)
{
    unsigned int shift = 6 - BASE*(bit_number & MASK);
    char base = (bit_array[bit_number >> SHIFT] & (3 << shift)) >> shift; 
    if(base == 0)
    {
        return 'A';
    }
    if(base == 3)
    {
        return 'T';
    }
    if(base == 1)
    {
        return 'G';
    }
    if(base == 2)
    {
        return 'C';
    }
}

void SetBit(char bit_array[], unsigned int bit_start_number, unsigned int bit_length, char * value) 
{ 
    unsigned int i = 0;
    for(i = 0; i < bit_length; i++)
    {
        SetBit(bit_array, bit_start_number + i, value[i]);
    }
} 

void SetBit(char bit_array[], unsigned int bit_start_number, unsigned int bit_length, char * value, int index) 
{ 
    unsigned int i = 0;
    for(i = 0; i < bit_length; i++)
    {
        SetBit(bit_array, bit_start_number + i, value[i]);
    }
    int j = ((bit_length/4) + 1)*4;
    for(i = i; i<j; i++){
        SetBit(bit_array, bit_start_number + i, 'A');
    }
} 

char * GetBit(char bit_array[], unsigned int bit_start_number, unsigned int bit_length)
{
    unsigned int i = 0;
    char * value = new char[bit_length + 1];
    for(i = 0; i < bit_length; i++)
    {
        value[i] = GetBit(bit_array, bit_start_number + i);
    }
    value[i] = '\0';
    return value;
}  

void GetBit(char bit_array[], unsigned int bit_start_number, unsigned int bit_length, char * value)
{
    unsigned int i = 0;
    for(i = 0; i < bit_length; i++)
    {
        value[i] = GetBit(bit_array, bit_start_number + i);
    }
    value[i] = '\0';
} 


#endif 
