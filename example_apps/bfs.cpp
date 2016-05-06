
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
typedef vid_t VertexDataType;       // vid_t is the vertex id type
typedef vid_t EdgeDataType;

static uint32_t repeat_count = 0;

/**
 * GraphChi programs need to subclass GraphChiProgram<vertex-type, edge-type> 
 * class. The main logic is usually in the update function.
 */
struct BfsProgram: public GraphChiProgram<VertexDataType, EdgeDataType> {

	bool converged;
	bool interval_converged;
	//std::vector<VertexDataType> bfs; 
	//BfsProgramInmem(int nvertices):bfs(nvertices, Max_Level){};	
	/**
	 *  Vertex update function.
	 *  On first iteration ,each vertex chooses a label = the vertex id.
	 *  On subsequent iterations, each vertex chooses the minimum of the neighbor's
	 *  label (and itself). 
	 */
	void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {

		if (scheduler) gcontext.scheduler->remove_tasks_now(vertex.id(), vertex.id());


		if (gcontext.iteration == 0) {
			//vertex.set_data(vertex.id());
			//std::cout<<"vid="<<vertex.id()<<"\t cid="<<bfs[vertex.id()]<<std::endl;
			//assert(bfs[vertex.id()] == Max_Level);
			if(vertex.id() == Root){
				//bfs[vertex.id()] = Root;
				vertex.set_data(Root);	
				for(int i=0; i<vertex.num_outedges(); i++){
					vertex.outedge(i)->set_data(Root);		
					if(scheduler)
						gcontext.scheduler->add_task(vertex.outedge(i)->vertex_id());	
				}
			}else{
				//bfs[vertex.id()] = Max_Level;
				vertex.set_data(Max_Level);	
				for(int i=0; i<vertex.num_outedges(); i++){
					vertex.outedge(i)->set_data(Max_Level);
					if(scheduler)
						gcontext.scheduler->add_task(vertex.outedge(i)->vertex_id());
					//vertex.outedge(i)->set_data(Max_Level);
				}
			}	
			//if (scheduler)  gcontext.scheduler->add_task(vertex.id());
		}else{
			/* On subsequent iterations, find the minimum label of my neighbors */
			vid_t curmin = vertex.get_data();
			//vid_t curmin = bfs[vertex.id()];
			vid_t old_value =curmin ; //vertex.get_data();
			for(int i=0; i < vertex.num_inedges(); i++) {
				vid_t nblabel = vertex.inedge(i)->get_data();
				curmin = std::min(nblabel+1, curmin); 
			}

			if(old_value > curmin){
				vertex.set_data(curmin);
				//vertex.modified=true;
				converged = false;
				interval_converged = false;
				for(int i=0; i<vertex.num_outedges(); i++){
					if( vertex.outedge(i)->get_data() > curmin){
						vertex.outedge(i)->set_data(curmin);		
						if(scheduler){
							vid_t edge_vid = vertex.outedge(i)->vertex_id();
							if(edge_vid > gcontext.interval_st){
								gcontext.scheduler->add_task(edge_vid, true);
							}else{
								gcontext.scheduler->add_task(edge_vid, false);
							}
						}
					}
				}
			}	
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
			//vid_t i=0;
			/*
			   FILE* fp;
			   fp = fopen((ginfo.filename+"ccc").c_str(),"w+");
			   for(i=0; i<bfs.size(); i++){
			   fprintf(fp, "%u\t%u\n", i, bfs[i]);
			   }
			   fclose(fp);
			   */
		}
	}
	bool repeat_updates(graphchi_context &gcontext){
		if(gcontext.iteration == 0)
			return false;
		else{
			if(0 == repeat_count % 100 || interval_converged){
			std::cout<<"exec_interval="<<gcontext.exec_interval<<"\trepeat_updates isconverged="<<interval_converged<<"\t iteration="
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

/*
struct ConnectedComponentsProgramCount : public GraphChiProgram<VertexDataType, EdgeDataType> {
	bool converged;
	
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {
				
			//assert(vertex.get_data() == bfs[vertex.id()]);	
			vertex.set_data(bfs[vertex.id()]);	
	} 
    void before_iteration(int iteration, graphchi_context &info) {
      	/////////////////////////////////////////////////
	  	converged = true;
		////////////////////////////////////////////////

    }

  void after_iteration(int iteration, graphchi_context &ginfo) {
        if (converged) {
            std::cout << "check pass!" << std::endl;
            ginfo.set_last_iteration(iteration);

        }
    }

};
*/

int main(int argc, const char ** argv) {
    /* GraphChi initialization will read the command line 
     arguments and the configuration file. */
    graphchi_init(argc, argv);

    /* Metrics object for keeping track of performance counters
     and other information. Currently required. */
    metrics m("blockupdate BFS");
    
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

	bool inmemmode =false;//= engine.num_vertices()*sizeof(VertexDataType) < (size_t)engine.get_membudget_mb()*1024L*1024L;
	if(inmemmode){
		logstream(LOG_INFO)<<"\nRunning BFS by holding vertices in-memory mode!!!!"<<std::endl;
		//engine.set_modifies_outedges(false);
		//engine.set_disable_outedges(true);
		//engine.set_modifies_inedges(true);
		//engine.set_only_adjacency(true);	
		//bfs(engine.num_vertices(), Max_Level);

		//bfs = (vid_t*)malloc(sizeof(vid_t)*(engine.num_vertices()+1));
		//assert(bfs != NULL);
		//memset(bfs,Max_Level,engine.num_vertices());
		//BfsProgramInmem program(engine.num_vertices());
		BfsProgram program;
        engine.run(program, niters);
		//ConnectedComponentsProgramCount prog;
		//engine.run(prog, 1);

		//free(bfs);
	}else{
		//assert(false);
       //ConnectedComponentsProgram program;
		BfsProgram program;
       engine.run(program, niters);
    }
    /* Run analysis of the connected components  (output is written to a file) */
    m.start_time("label-analysis");
    
    analyze_labels<vid_t>(filename, 50);
    //count_bfs_levels<vid_t>(filename, Max_Level);
    
    m.stop_time("label-analysis");

    /* Report execution metrics */
    metrics_report(m);
    return 0;
}

