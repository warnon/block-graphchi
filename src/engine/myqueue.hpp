#include<stdlib.h>

public myqueue{
private:

int * ptr;
int start;
int end;
int size;
int len;

public :
myqueue(int nshards):start(0),end(0){
	size=nshards;
	ptr=(int*)malloc(sizeof(int)*size)
	start=end=0;
	len=0;
}

bool isEmpty(){
	return 0==len ;
}

bool isFull(){
	return len == num;
}
int getNext(){
	int val=-1;
	if(len != 0){
		val = ptr[start];
		start = (start+1)%size;
		len--;
	}
	return val;
}

bool  addNext(int index){
	if( !isFull() ){
		ptr[end] = index;
		end = (end+1)%size;
		len++;
		return true;	
	}
	return false;
}	
};
