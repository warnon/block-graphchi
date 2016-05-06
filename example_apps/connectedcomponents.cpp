
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
vid_t Max_Cid =(vid_t)-1; 

//vid_t* wcc= NULL;
//std::vector<vid_t> wcc; 
//float   diff_iter=0.0f;
/**
 * Type definitions. Remember to create suitable graph shards using the
 * Sharder-program. 
 */
typedef vid_t VertexDataType;       // vid_t is the vertex id type
typedef vid_t EdgeDataType;

static uint32_t repeat_count = 0;
/**
 * GraphChi programs need to subclass GraphChiProgram<vertex-type, edge-type> 
 * class. The main logic is usually in the update function.
 */
struct ConnectedComponentsProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {
   	bool interval_converged; 
    bool converged;
    /**
     *  Vertex update function.
     *  On first iteration ,each vertex chooses a label = the vertex id.
     *  On subsequent iterations, each vertex chooses the minimum of the neighbor's
     *  label (and itself). 
     */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {
       //std::cout<<"vid="<<vertex.id()<<std::endl; 
       if (scheduler) gcontext.scheduler->remove_tasks_now(vertex.id(), vertex.id());
		/*
		if(vertex.id() == 50){
			std::cout<<"degree="<<vertex.num_edges()<<std::endl;	
		}
		*/
		//assert();	
        if (gcontext.iteration == 0) {
            vertex.set_data(vertex.id());
            if (scheduler)  gcontext.scheduler->add_task(vertex.id(), false);
        }
        
        /* On subsequent iterations, find the minimum label of my neighbors */
        vid_t curmin = vertex.get_data();
        for(int i=0; i < vertex.num_edges(); i++) {
            vid_t nblabel = vertex.edge(i)->get_data();
			if (gcontext.iteration == 0){
				 nblabel = vertex.edge(i)->vertex_id();  // Note!
			}
            curmin = std::min(nblabel, curmin); 
        }

		vertex.set_data(curmin);
        vid_t label = vertex.get_data();

        if (gcontext.iteration > 0) {
            for(int i=0; i < vertex.num_edges(); i++) {
                if (label < vertex.edge(i)->get_data()){
                    vertex.edge(i)->set_data(label);
                    /* Schedule neighbor for update */
                   if (scheduler){
						vid_t edge_vid = vertex.edge(i)->vertexid;
						if(edge_vid >= gcontext.interval_st){
							//schedule edge_vid to update in this iteration
						 	gcontext.scheduler->add_task(vertex.edge(i)->vertex_id(), true);
						}else{
							//schedule edge_vid to update in next iteration
						 	gcontext.scheduler->add_task(vertex.edge(i)->vertex_id(), false);
						}
					}
                    converged = false;
					//vertex.modified=true;
					interval_converged = false;
                }
            }
        } else if (gcontext.iteration == 0) {
            for(int i=0; i < vertex.num_outedges(); i++) {
                vertex.outedge(i)->set_data(label);
            }
        }
    }    
	bool repeat_updates(graphchi_context &gcontext){
		if(gcontext.iteration == 0)
			return false;
		else{
			if(0 == repeat_count % 100 || interval_converged){
				std::cout<<"exec_interval="<<gcontext.exec_interval<<"\trepeat_updates isconverged="<<interval_converged<<"\t iteration="
					<<gcontext.iteration<<"\t repeat_count="<<repeat_count<<std::endl;
			}
			repeat_count = interval_converged ? 0 : repeat_count+1;
			return !interval_converged;
			//return false;
		}
	} 

    /**
     * Called before an iteration starts.
     */
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
    metrics m("connected-components");
    
    /* Basic arguments for application */
    std::string filename = get_option_string("file");  // Base filename
    int niters           = get_option_int("niters", 10000); // Number of iterations (max)
    scheduler            = get_option_int("scheduler", false);
	int terminate		 = get_option_int("term", 0);
 //   diff_iter			=get_option_int("diff",0);
    /* Process input file - if not already preprocessed */
    int nshards             = (int) convert_if_notexists<EdgeDataType>(filename, get_option_string("nshards", "auto"));
	if(0 != terminate)
   		return 0; 
        /* Run */
	graphchi_engine<VertexDataType, EdgeDataType> engine(filename, nshards, scheduler, m); 

	//bool inmemmode = engine.num_vertices()*sizeof(VertexDataType) < (size_t)engine.get_membudget_mb()*1024L*1024L;
	bool inmemmode = false;//engine.num_vertices()*sizeof(VertexDataType) < (size_t)engine.get_membudget_mb()*1024L*1024L;
	if(inmemmode){
		logstream(LOG_INFO)<<"\nRunning WCC by holding vertices in-memory mode!!!!"<<std::endl;
		//engine.set_modifies_outedges(false);
		//engine.set_disable_outedges(true);
		//engine.set_modifies_inedges(true);
		//engine.set_only_adjacency(true);	
		//wcc(engine.num_vertices(), Max_Cid);
		//wcc = (vid_t*)malloc(sizeof(vid_t)*(engine.num_vertices()+1));
		//assert(wcc != NULL);
		//memset(wcc,Max_Cid,engine.num_vertices());
		//ConnectedComponentsProgramInmem program(engine.num_vertices());
		//ConnectedComponentsProgramInmem program(engine.num_vertices());
        //engine.run(program, niters);
		//ConnectedComponentsProgramCount prog;
		//engine.run(prog, 1);
		//free(wcc);
	}else{
		//assert(false);
        ConnectedComponentsProgram program;
        engine.run(program, niters);
    }
    /* Run analysis of the connected components  (output is written to a file) */
    m.start_time("label-analysis");
    
    analyze_labels<vid_t>(filename);
    
    m.stop_time("label-analysis");
    

    /* Report execution metrics */
    metrics_report(m);
    return 0;
}

