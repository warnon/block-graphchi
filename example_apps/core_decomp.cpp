
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
 * Template for GraphChi applications. To create a new application, duplicate
 * this template.
 */


#include <cmath>
#include <string>

#include "util/labelanalysis.hpp"
#include "util/toplist.hpp"
#include "graphchi_basic_includes.hpp"

#define MAX_LEVEL 100000
//#define MAXK	50

using namespace graphchi;

/**
  * Type definitions. Remember to create suitable graph shards using the
  * Sharder-program. 
  */
  
  
  
struct bidirectional_label {
    int smaller_degree;
	//int smaller_level;
    int larger_degree;
	bidirectional_label(){
		smaller_degree = larger_degree = -1;

	}
 	bidirectional_label( int x){
		smaller_degree = larger_degree = -1;
	} 

	int & neighbor_degree(vid_t myid, vid_t nbid) {
        //assert(larger_cid != 0xffffffffu);
        //assert(smaller_cid != 0xffffffffu);
        
        if (myid < nbid) {
            return larger_degree;
        } else {
            return smaller_degree;
        }
    }

    int & my_degree(vid_t myid, vid_t nbid) {
        //assert(larger_cid != 0xffffffffu);
        //assert(smaller_cid != 0xffffffffu)
        if (myid < nbid) {
            return smaller_degree;
        } else {
            return larger_degree;
        }
    }

};

struct Vertinfo {
    //int comp_id;
	//int level;
    bool confirm;
	int degree;
	int core_num;
    Vertinfo() : degree(-1),core_num(MAX_LEVEL) ,confirm(false){}
    Vertinfo(int deg) : degree(deg),core_num(MAX_LEVEL) ,confirm(false){}
    Vertinfo(vid_t deg, bool confirmed) : degree(deg),core_num(MAX_LEVEL) ,confirm(confirmed){}
   	
	int get_core_num(){
		return core_num;
	}
	int get_degree(){
		return degree;
	}
	// modify vertex degree when their neighbors were visited
	void set_degree(int deg){
		degree = deg;
	} 
	bool is_confirm(){
		return confirm;
	}
    friend std::ostream& operator<< (std::ostream &out, Vertinfo &scc) {
        out << scc.core_num;
        return out;
    }
    
};
/* Overloaded operators to help with labelanalysis.hpp */

bool operator<(const Vertinfo &a, const Vertinfo &b);
bool operator<(const Vertinfo &a, const Vertinfo &b) {
    return a.core_num < b.core_num;
}

bool operator==(const Vertinfo &a, const Vertinfo &b);
bool operator==(const Vertinfo &a, const Vertinfo &b) {
    return a.core_num == b.core_num;
}
bool operator!=(const Vertinfo &a, const Vertinfo &b);
bool operator!=(const Vertinfo &a, const Vertinfo &b) {
    return a.core_num != b.core_num;
}

typedef Vertinfo VertexDataType;
typedef bidirectional_label EdgeDataType;

int global_core = 0;
// global variables

bool scheduler = false; 

int         iterationcount = 0;
//bool        scheduler = false;
/**
  * GraphChi programs need to subclass GraphChiProgram<vertex-type, edge-type> 
  * class. The main logic is usually in the update function.
  */


struct Core_decomp : public GraphChiProgram<VertexDataType, EdgeDataType> {
	bool any_left;// =false; 
	bool increase_core;
	bool local_converged;
	/**
	 *  Vertex update function.
	 */
	void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {

		if (gcontext.iteration == 0)
		{
			/* On first iteration, initialize vertex (and its edges). This is usually required, because
			   on each run, GraphChi will modify the data files. To start from scratch, it is easiest
			   do initialize the program in code. Alternatively, you can keep a copy of initial data files. */
				VertexDataType vertexdata = vertex.get_data();
				vertexdata.confirm = false;
				vertexdata.degree = vertex.num_edges();
				vertexdata.core_num = MAX_LEVEL;
				//assert(vertexdata.comp_id == 0);
				//vertexdata.level = source_level;
				vertex.set_data(vertexdata);
				for(int i=0; i<vertex.num_edges(); i++)
				{
					graphchi_edge<EdgeDataType> * e = vertex.edge(i);
					bidirectional_label edgelabel = e->get_data();
					edgelabel.my_degree(vertex.id(), e->vertexid) =vertexdata.degree; //vertex.id();
					//edgelabel.my_level(vertex.id(), e->vertexid) = source_level;
					//edgelabel.enable_propagate();
					e->set_data(edgelabel);
					//if(scheduler)
					//std::cout<<"add vertex id========"<<vertex.edge(i)->vertexid<<"===to schedule"<<std::endl;
				}// first vertex(root) addition should always succeed
				if(vertexdata.degree == 0){
					vertexdata.core_num = global_core;	
					vertexdata.confirm = true;
					vertex.set_data(vertexdata);
				}/*else if(vertexdata.degree == 1){
						if(scheduler)
						gcontext.scheduler->add_task(vertex.id());
				}*/
		} else {
			/* Do computation */ 	
			VertexDataType vertexdata = vertex.get_data();	
			if(vertexdata.confirm)
				return;
			int deleted_edges = 0; //= vertexdata.degree;
			int real_degree = vertex.num_edges();
			for(int i=0; i < real_degree; i++)
			{
				graphchi_edge<EdgeDataType> * e = vertex.edge(i);
				bidirectional_label edgelabel = e->get_data();
				if(0 == edgelabel.neighbor_degree(vertex.id(), e->vertexid))
					deleted_edges++;
		
			}
			if(real_degree-deleted_edges <= global_core)
			{
				vertexdata.degree = real_degree - deleted_edges;//get_cid();//my_level+1;//my_cid;
				assert(vertexdata.degree >= 0);
				//assert(vertexdata.comp_id < 15);
				vertexdata.confirm = true;
				vertexdata.core_num = global_core;
				vertex.set_data(vertexdata);
				//modifyEdges(vertex);
				//changed = true;
				for(int i=0; i < real_degree; i++)
				{
					graphchi_edge<EdgeDataType>* e =  vertex.edge(i);
					bidirectional_label edgelabel = e->get_data();
					edgelabel.my_degree(vertex.id(), e->vertexid) = vertexdata.degree;//my_cid;
					//edgelabel.my_level(vertex.id(), e->vertexid) = vertexdata.level;
					//int nlevel = edgelabel.neighbor_level(vertex.id(), e->vertexid);
					e->set_data(edgelabel);
					/*
					if(vertexdata.level < nlevel)
					{
						if(scheduler)
							gcontext.scheduler->add_task(e->vertexid);
					}
					*/
				}
				increase_core = false;
				local_converged = false;
			}
			any_left = true;
		}
	}

	/**
	 * Called before an iteration starts.
	 */
	void before_iteration(int iteration, graphchi_context &gcontext) {
		//changed = gcontext.iteration == 0;
		increase_core = true;
		any_left = false;	
	}

	/**
	 * Called after an iteration has finished.
	 */
	void after_iteration(int iteration, graphchi_context &gcontext) {	
		if(iteration == 0)
			return;	
		if(!any_left){
			gcontext.set_last_iteration(iteration);				
		}else{
			if(increase_core)
				global_core++;
		}
		std::cout<<"iteration#"<<iteration<<"\t core_num="<<global_core<<"========================="<<global_core<<std::endl;
	}

    virtual bool repeat_updates(graphchi_context &gcontext) {
            return !local_converged;
    }

	/**
	 * Called before an execution interval is started.
	 */
	void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {        
			local_converged = gcontext.iteration > 0;
			
	}

	/**
	 * Called after an execution interval has finished.
	 */
	void after_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {        
		 //if()
	}

};

int main(int argc, const char ** argv) {
    /* GraphChi initialization will read the command line 
       arguments and the configuration file. */
    graphchi_init(argc, argv);
    
    /* Metrics object for keeping track of performance counters
       and other information. Currently required. */
    metrics m("Core_decomp-Report");
    
    /* Basic arguments for application */
    std::string filename = get_option_string("file");  // Base filename

    int niters           = get_option_int("niters", 100000); // Number of iterations

    scheduler       = get_option_int("scheduler", false); // Whether to use selective scheduling
    
	
    /* Detect the number of shards or preprocess an input to create them */
    int nshards          = convert_if_notexists<int,EdgeDataType>(filename, 
                                                            get_option_string("nshards", "auto"));


	graphchi_engine<VertexDataType, EdgeDataType> engine(filename, nshards, scheduler, m);
	
	Core_decomp decompprogram;
	engine.run(decompprogram, niters);

	m.start_time("result analysis");	
	//count_labels<VertexDataType>(filename,20,1);
	analyze_labels<VertexDataType>(filename, 50);
	m.start_time("result analysis");	
    metrics_report(m);
	
    return 0;
}
