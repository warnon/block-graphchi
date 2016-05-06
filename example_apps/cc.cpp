
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
#include<stdio.h>

#include "graphchi_basic_includes.hpp"
#include "util/labelanalysis.hpp"

using namespace graphchi;

int         iterationcount = 0;
bool        scheduler = false;
vid_t Root = 0;
//vid_t Max_Level =(vid_t)-1; 
vid_t Max_Level = 100000; 

//vid_t* bfs= NULL;
//std::vector<vid_t> bfs; 
//float   diff_iter=0.0f;
/**
 * Type definitions. Remember to create suitable graph shards using the
 * Sharder-program. 
 */
/*
struct Vertexinfo{
	vid_t level;
	vid_t degree;
	Vertexinfo():level(0), degree(0){};
	Vertexinfo(vid_t lv):level(lv), degree(0){};
	Vertexinfo(vid_t lv, vid_t deg):level(lv), degree(deg){};
	
	vid_t get_value(){
		return level;
	}
	friend std::ostream& operator << (std::ostream& out, Vertexinfo& vertex){
		out<<vertex.level;
		return out;
	}
};

bool operator < (const Vertexinfo& v1, const Vertexinfo& v2){
	return v1.level < v2.level;
}
bool operator == (const Vertexinfo& v1, const Vertexinfo& v2){
	return v1.level == v2.level;
}
bool operator != (const Vertexinfo& v1, const Vertexinfo& v2){
	return v1.level != v2.level;
}
bool operator > (const Vertexinfo& v1, const Vertexinfo& v2){
	return v1.level > v2.level;
}
*/
//typedef Vertexinfo VertexDataType;       // vid_t is the vertex id type
typedef vid_t VertexDataType;       // vid_t is the vertex id type
typedef vid_t EdgeDataType;

struct level_pair{
	vid_t label;
	vid_t count;
	// constructor
	level_pair():label(0), count(0){};	
	level_pair(vid_t lb, vid_t size):label(lb), count(size){};
	level_pair(vid_t lb):label(lb), count(0){};
};

static uint32_t repeat_count = 0;
std::vector<level_pair>	level_size; 


struct CCProgram: public GraphChiProgram<VertexDataType, EdgeDataType> {
	bool converged;
	bool interval_converged;
	void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {
		//if (scheduler) gcontext.scheduler->remove_tasks(vertex.id(), vertex.id());
		if (gcontext.iteration == 0) {
			vid_t min_nid = vertex.id();	
			//vertex.set_data(Vertexinfo(min_nid, (vid_t)vertex.num_edges()));	
			vertex.set_data(min_nid);	
			for(int i=0; i<vertex.num_edges(); i++){
				min_nid = std::min(min_nid, vertex.edge(i)->vertex_id());	
			}	
			//vertex.set_data(Vertexinfo(min_nid, (vid_t)vertex.num_edges()));	
			//vertex.set_data(min_nid);	
			for(int i=0; i<vertex.num_outedges(); i++){
				//vertex.outedge(i)->set_data(vertex.id());
				vertex.outedge(i)->set_data(min_nid);
				if(scheduler)
					gcontext.scheduler->add_task(vertex.outedge(i)->vertex_id());
			}
		}else{
			//Vertexinfo vdata = vertex.get_data();
			vid_t curmin = vertex.get_data();
			vid_t old_value =curmin ; //vertex.get_data();
			for(int i=0; i < vertex.num_edges(); i++) {
				vid_t nblabel = vertex.edge(i)->get_data();
				curmin = std::min(nblabel, curmin); 
			}
			if(old_value > curmin ){
				//vdata.level = curmin;	
				vertex.set_data(curmin);
				//vertex.modified=true;
				converged = false;
				interval_converged = false;
				for(int i=0; i<vertex.num_edges(); i++){
					graphchi_edge<EdgeDataType> *  e = vertex.edge(i);
					//if(e->get_data() > curmin){
						//vertex.edge(i)->set_data(curmin);		
						e->set_data(curmin);		
						if(scheduler)
							gcontext.scheduler->add_task(vertex.edge(i)->vertex_id());
					//}
				}
			}else{
				//vertex.modified=false;
			}	
		}
	}    

	void before_iteration(int iteration, graphchi_context &info) {
		iterationcount++;
		/////////////////////////////////////////////////
		converged = iteration > 0;
		////////////////////////////////////////////////

	}

	/**
	 * Called after an iteration has finished.
	 */
	void after_iteration(int iteration, graphchi_context &ginfo) {
		
		if (converged) {
			std::cout << "Converged!" << std::endl;
			ginfo.set_last_iteration(iteration);
		}
	
	}
	bool repeat_updates(graphchi_context &gcontext){
		if(gcontext.iteration == 0)
			return false;
		else{
			if(0 == repeat_count % 100 || interval_converged){
			std::cout<<"repeat_updates isconverged="<<interval_converged<<"\t iteration="
				<<gcontext.iteration<<"\trepeat_count="<<repeat_count<<std::endl;
			}
			repeat_count = interval_converged ? 0 : repeat_count+1;
			return !interval_converged;
			//return false;
		}
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



int main(int argc, const char ** argv) {
    /* GraphChi initialization will read the command line 
     arguments and the configuration file. */
    graphchi_init(argc, argv);

    /* Metrics object for keeping track of performance counters
     and other information. Currently required. */
    metrics m("blockupdate BFS-Partition");
    
    /* Basic arguments for application */
    std::string filename = get_option_string("file");  // Base filename
    int niters           = get_option_int("niters", 100000); // Number of iterations (max)
    scheduler            = get_option_int("scheduler", false);
	Root	= (vid_t)get_option_int("root", 0);
 //   diff_iter			=get_option_int("diff",0);
    /* Process input file - if not already preprocessed */
    int nshards             = (int) convert_if_notexists<EdgeDataType>(filename, get_option_string("nshards", "auto"));
    
        /* Run */
	graphchi_engine<VertexDataType, EdgeDataType> engine(filename, nshards, scheduler, m); 

	/* Run analysis of the connected components  (output is written to a file) */
	// first run cc to detect all communities in the graph
	CCProgram ccprogram;
    engine.run(ccprogram, niters);

	//BFSProgram bfsprogram;

    m.start_time("label-analysis");
    
    analyze_labels<vid_t>(filename, 20);
    //count_bfs_levels<vid_t, level_pair >(filename, Max_Level, level_size);
    //count_bfs_levels<VertexDataType, level_pair>(filename, Max_Level, level_size);
    
    m.stop_time("label-analysis");
	std::cout<<"partition bfs result size="<<level_size.size()<<std::endl;	
    /* Report execution metrics */
    metrics_report(m);
    return 0;
}

