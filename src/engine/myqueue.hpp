#include<cstdlib>
#include<cstdio>
#include<stdint.h>
#include<assert.h>
//#include "util/pthread_tools.hpp"
using namespace graphchi;
//namespace graphchi{

class myqueue{
	private:
		mutex lock;

		float * ptr_curr;
		float * ptr_next;
		int start;
		int end;
		int size;
		int len;

	public :
		myqueue(int nshards):start(0),end(0){
			size=nshards;
			ptr_curr = (float*)malloc(sizeof(float)*size);
				assert(ptr_curr != NULL);
			ptr_next = (float*)malloc(sizeof(float)*size);	
				assert(ptr_next != NULL);
			for(int i=0; i<size; i++){
				ptr_curr[i] = ptr_next[i] = 0;
			}
			start=end=0;
			len=0;
		}
		
		myqueue(){
			ptr_next = ptr_curr = NULL;	
			start = end = 0;	
		}

		void resize(int sz){
			if(sz != size){
				if(ptr_curr != NULL) free(ptr_curr);	
				if(ptr_next != NULL) free(ptr_next);
				size = sz;
				ptr_curr = (float*)malloc(sizeof(float)*size);
					assert(ptr_curr != NULL);
				ptr_next = (float*)malloc(sizeof(float)*size); 
				assert(ptr_next != NULL);
				for(int i=0; i<size; i++){
					ptr_curr[i] = 0;
					ptr_next[i] = 0;	
				}
				start=end=0;
			}
		}

		void clear(){
			for(int idx=0; idx<size; idx++){
				/*
				__sync_fetch_and_add(&ptr_curr[idx], 0);	
				__sync_fetch_and_add(&ptr_next[idx], 0);	
				*/
				ptr_curr[idx] = 0;
				ptr_curr[idx] = 0;
			}	
		}

		~myqueue(){
			if(ptr_curr != NULL)
				free(ptr_curr);
			ptr_curr = NULL;
			if(ptr_next != NULL)
				free(ptr_next);
			ptr_next = NULL;
		}

		bool is_empty(){
			return 0==len ;
		}

		bool is_full(){
			return len == size;
		}
		//add priority to subgraph idx;
		void add_value(int idx, float value){
			if(idx < size && value > 0){
				//__sync_fetch_and_add(ptr_next+idx, value);	
				//__sync_fetch_and_add(&ptr_next[idx], value);	
				lock.lock();
				ptr_next[idx] += value;	
				lock.unlock();
			}	
		}	
		//get next subgraph with biggest priority
		int get_next(){
			/*
			int val=-1;
			if(len != 0){
				val = ptr[start];
				start = (start+1)%size;
				len--;
			}
			return val;
			*/
			int idx = 0;	
			for(int i=1; i<size; i++){
				if(ptr_curr[idx] < ptr_curr[i])		
					idx = i;
			}
			//__sync_fetch_and_add(ptr_curr+idx, 0);		
			//__sync_fetch_and_add(&ptr_curr[idx], 0);		
			ptr_curr[idx] = 0;
			return idx;
		}

		void new_iteration(){
			float* tmp = ptr_curr;	
			ptr_curr = ptr_next;
			ptr_next = tmp;
		}
		/*
		bool add_next(int index){
			if( !isFull() ){
				ptr[end] = index;
				end = (end+1)%size;
				len++;
				return true;	
			}
			return false;
		}	
		*/
};
//};
