
/**
 * @file
 * @author  Aapo Kyrola <akyrola@cs.cmu.edu>
 * @version 1.0
 *
 * @section LICENSE
 *
 * Copyright [2012] [Aapo Kyrola, Guy Blelloch, Carlos Guestrin / Carnegie Mellon University]
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 
 *
 * @section DESCRIPTION
 *
 * Application for computing the connected components of a graph.
 * The algorithm is simple: on first iteration each vertex sends its
 * id to neighboring vertices. On subsequent iterations, each vertex chooses
 * the smallest id of its neighbors and broadcasts its (new) label to
 * its neighbors. The algorithm terminates when no vertex changes label.
 *
 * @section REMARKS
 *
 * This application is interesting demonstration of the asyncronous capabilities
 * of GraphChi, improving the convergence considerably. Consider
 * a chain graph 0->1->2->...->n. First, vertex 0 will write its value to its edges,
 * which will be observed by vertex 1 immediatelly, changing its label to 0. Nexgt,
 * vertex 2 changes its value to 0, and so on. This all happens in one iteration.
 * A subtle issue is that as any pair of vertices a<->b share an edge, they will
 * overwrite each others value. However, because they will be never run in parallel
 * (due to deterministic parallellism of graphchi), this does not compromise correctness.
 *
 * @author Aapo Kyrola
 */


#include <cmath>
#include <string>
#include <stdio.h>
#include <map>

#include "graphchi_basic_includes.hpp"
#include "util/labelanalysis.hpp"
#include "util/bfsanalysis.hpp"

using namespace graphchi;

int         iterationcount = 0;
bool        scheduler = false;
//vid_t Root = 0;
const int top = 20;
//vid_t Max_Level =(vid_t)-1; 
//vid_t Max_Level = 1000000; 

//vid_t* bfs= NULL;
//std::vector<vid_t> bfs; 
//float   diff_iter=0.0f;
/**
 * Type definitions. Remember to create suitable graph shards using the
 * Sharder-program. 
 */
/*
struct CCVertexinfo{
	vid_t ccid;
	CCVertexinfo():ccid(0){};
	CCVertexinfo(vid_t cid):ccid(cid){};
	vid_t get_value(){
		return ccid;
	}
	friend std::ostream& operator <<(std::ostream& out, CCVertexinfo& v){
		return out<<v.ccid;	
	}
};

bool operator < (const CCVertexinfo& v1, const CCVertexinfo& v2){
	return v1.ccid < v2.ccid;
}
bool operator == (const CCVertexinfo& v1, const CCVertexinfo& v2){
	return v1.ccid == v2.ccid;
}
bool operator != (const CCVertexinfo& v1, const CCVertexinfo& v2){
	return v1.ccid != v2.ccid;
}
bool operator > (const CCVertexinfo& v1, const CCVertexinfo& v2){
	return v1.ccid > v2.ccid;
}
*/
/*
struct Vertexinfo{
	vid_t ccid;
	vid_t level; 
	Vertexinfo():ccid(0xffffffff),level(Max_Level){};
	Vertexinfo(vid_t cid, vid_t lv):ccid(cid),level(lv){};
	void set_ccid(vid_t cid){
		ccid = cid;
	}
	void set_level(vid_t lv){
		level = lv;
	}
	vid_t get_value(){
		return level;
	}
	
	friend std::ostream& operator << (std::ostream& out, Vertexinfo& vertex){
		out<<vertex.level<<"\t"<<vertex.ccid;
		return out;
	}
};
*/

struct bidirectional_label {
    vid_t smaller_one;
    vid_t larger_one;
	float weight;
    vid_t & neighbor_label(vid_t myid, vid_t nbid) {
        //assert(larger_one != 0xffffffffu);
        //assert(smaller_one != 0xffffffffu);
        
        if (myid < nbid) {
            return larger_one;
        } else {
            return smaller_one;
        }
    }
    
    vid_t & my_label(vid_t myid, vid_t nbid) {
        //assert(larger_one != 0xffffffffu);
        //assert(smaller_one != 0xffffffffu);
        
        if (myid < nbid) {
            return smaller_one;
        } else {
            return larger_one;
        }
    }
    
	/*
    // Annoying hack
    bool deleted() {
        return smaller_one == 0xffffffffu;
    }
	*/
};

/*
bool operator < (const BFSVertexinfo& v1, const BFSVertexinfo& v2){
	return v1.level < v2.level;
}
bool operator == (const BFSVertexinfo& v1, const BFSVertexinfo& v2){
	return v1.level == v2.level;
}
bool operator != (const BFSVertexinfo& v1, const BFSVertexinfo& v2){
	return v1.level != v2.level;
}
bool operator > (const BFSVertexinfo& v1, const BFSVertexinfo& v2){
	return v1.level > v2.level;
}
*/
//typedef Vertexinfo VertexDataType;       // vid_t is the vertex id type
//typedef vid_t EdgeDataType;
//typedef Vertexinfo VertexDataType;       // vid_t is the vertex id type
typedef vid_t VertexDataType;       // vid_t is the vertex id type
typedef bidirectional_label EdgeDataType;
bool flag_weight = false;
//static void VARIABLE_IS_NOT_USED parse<EdgeDataType>(EdgeDataType& edata, const char* s){
static void  parse(EdgeDataType& edata, const char* s){
	flag_weight = true;
	edata.weight = atof(s);	
}
/*
struct level_pair{
	vid_t label;
	vid_t count;
	// constructor
	level_pair():label(0), count(0){};	
	level_pair(vid_t lb, vid_t size):label(lb), count(size){};
	level_pair(vid_t lb):label(lb), count(0){};
};

std::vector<level_pair>	cc_size; 
*/
//static uint32_t repeat_count = 0;
//int top = 0;
/*
vid_t getNewId(int level){
	assert(level < );
	lock.lock();
		
}
*/
mutex lock;
FILE* vfout = NULL;
FILE* efout = NULL;
int nshards = 0;
std::vector<int> prefix_sum;

struct ReMapProgram: public GraphChiProgram<VertexDataType, EdgeDataType> {
	bool converged;
	bool interval_converged;
	void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {
		if (gcontext.iteration == 0){
				//Vertexinfo vdata = vertex.get_data();  
				//id_th = vertex.id();	
				vid_t row = vertex.id() / nshards; 
				//vid_t new_id = row * nshards + prefix_sum[vertex.id() % nshards];	
				vid_t new_id = row + prefix_sum[vertex.id() % nshards];	
				vertex.set_data(new_id);
				//source vertex in each CC
				//vdata.level = 0;
				//vertex.set_data(vdata);	
				//vid_t new_id = getNewId(vdata.ccid, vdata.level);
				for(int i=0; i<vertex.num_edges(); i++){
					bidirectional_label edata = vertex.edge(i)->get_data();
					edata.my_label(vertex.id(), vertex.edge(i)->vertex_id()) = new_id;
					vertex.edge(i)->set_data(edata);
				}
				lock.lock();
					fprintf(vfout, "%u\t%u\n", new_id, vertex.id());
				lock.unlock();
		}else{
			for(int i=0; i<vertex.num_outedges(); i++){
				bidirectional_label edata = vertex.outedge(i)->get_data();	
				vid_t my_id = edata.my_label(vertex.id(), vertex.outedge(i)->vertex_id());		
				vid_t nb_id = edata.neighbor_label(vertex.id(), vertex.outedge(i)->vertex_id());
				if(my_id == nb_id){
					std::cout<<"my_id="<<vertex.id()<<"\tmy_label="<<my_id <<"\tnb_label="<<nb_id
							<<"\tnb_vid="<<vertex.outedge(i)->vertex_id()<<std::endl;
					assert(my_id != nb_id);
				}
				if(!flag_weight){
					lock.lock();
					fprintf(efout, "%u\t%u\n", my_id, nb_id);
					lock.unlock();
				}else{
					lock.lock();
					fprintf(efout, "%u\t%u\t%.3f\n", my_id, nb_id, edata.weight);
					lock.unlock();
				}
			}			
		}
	}    
	/**
	 * Called before an iteration starts.
	 */
	void before_iteration(int iteration, graphchi_context &info) {

	}

	/**
	 * Called after an iteration has finished.
	 */
	void after_iteration(int iteration, graphchi_context &ginfo) {
		/*	
			if (converged) {
			std::cout << "Converged!" << std::endl;
			ginfo.set_last_iteration(iteration);
			}
			*/
		if(iteration == 1)	
			ginfo.set_last_iteration(iteration);
	}
	bool repeat_updates(graphchi_context &gcontext){
		return false;
	} 
	/**
	 * Called before an execution interval is started.
	 */
	void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &ginfo) {        
		interval_converged = true;		
	}

	/**
	 * Called after an execution interval has finished.
	 */
	void after_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &ginfo) {        
	}
};

void prefix_sum_calculate(int shards, int nvertices){

	prefix_sum.resize(shards, 0);
	int interval_size = nvertices / shards;
	if(nvertices % shards == 0 ){
		for(int i=0; i<shards; i++)
			prefix_sum[i] = interval_size;
	}else{
		//separator index is  
		int separator = nvertices % shards;		
		for(int i=1; i<=shards; i++){
			if(i <= separator){
				prefix_sum[i-1] = interval_size+1;
			}else{
				prefix_sum[i-1] = interval_size;	
			}
		}		
	}
	int tmp_sum = 0;		
	for(int i=0; i<shards; i++){
		int tmp = tmp_sum;
		tmp_sum += prefix_sum[i];	
		prefix_sum[i] = tmp; 
	}
	assert(nvertices == tmp_sum);
	for(int i=0; i<shards; i++){
		std::cout<<i<<" , prefix_sum="<<prefix_sum[i]<<std::endl;
	}
}

int main(int argc, const char ** argv) {
    /* GraphChi initialization will read the command line 
     arguments and the configuration file. */
    graphchi_init(argc, argv);

    /* Metrics object for keeping track of performance counters
     and other information. Currently required. */
    metrics m("blockupdate BFS-Partition");
    
    /* Basic arguments for application */
    std::string filename = get_option_string("file");  // Base filename
    int niters           = get_option_int("niters", 2); // Number of iterations (max)
    scheduler            = get_option_int("scheduler", false);
    /* Process input file - if not already preprocessed */
    nshards             = (int) convert_if_notexists<EdgeDataType>(filename, get_option_string("nshards", "auto"));
    /* Run */
	vfout = fopen((filename+".hv").c_str(), "w+");	
	efout = fopen((filename+".he").c_str(), "w+");
	assert(vfout != NULL && efout != NULL);

	fprintf(vfout, "# new_vid  old_vid\n");
	fprintf(efout, "# new_src  new_dst\n");
	fflush(vfout);	
	fflush(efout);
	std::cout<<"remap is started"<<std::endl;
	graphchi_engine<VertexDataType, EdgeDataType> engine2(filename, nshards, scheduler, m); 
	// get num verticesfrom engine
	int num_vertices = engine2.num_vertices();
	assert(num_vertices > 0);

	prefix_sum_calculate(nshards, num_vertices);

	ReMapProgram remapprogram;				
	engine2.run(remapprogram, niters);
	fclose(vfout);
	fclose(efout);
	std::cout<<"hash based partitioning has been remap"<<std::endl;
	
	m.start_time("label-analysis");
    /* Report execution metrics */
    metrics_report(m);
    return 0;
}

