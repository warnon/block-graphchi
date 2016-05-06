
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

/**
 * GraphChi programs need to subclass GraphChiProgram<vertex-type, edge-type> 
 * class. The main logic is usually in the update function.
 */
struct ConnectedComponentsProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {
    
    bool converged;
    /**
     *  Vertex update function.
     *  On first iteration ,each vertex chooses a label = the vertex id.
     *  On subsequent iterations, each vertex chooses the minimum of the neighbor's
     *  label (and itself). 
     */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {
        
       if (scheduler) gcontext.scheduler->remove_tasks(vertex.id(), vertex.id());
        

        if (gcontext.iteration == 0) {
            vertex.set_data(vertex.id());
            if (scheduler)  gcontext.scheduler->add_task(vertex.id());
        }
        
        /* On subsequent iterations, find the minimum label of my neighbors */
        vid_t curmin = vertex.get_data();

	/*	if(vertex.id() == 69697){
			logstream(LOG_INFO)<<"vertex 69697: ------------"<<vertex.get_data()<<std::endl;

		}*/
        vid_t old_value =curmin ; //vertex.get_data();
        for(int i=0; i < vertex.num_edges(); i++) {
            vid_t nblabel = vertex.edge(i)->get_data();
			if (gcontext.iteration == 0){
				 nblabel = vertex.edge(i)->vertex_id();  // Note!
            //     vertex.edge(i)->set_data(10000000);
			}
            curmin = std::min(nblabel, curmin); 
        }
		
      	/*if( curmin == 0 && vertex.id() != 0) {
			logstream(LOG_INFO)<<"vertex===========vid:"<<vertex.id()<<std::endl;
       	//	exit(1);
			}*/
		
	 /* Set my label */
	   // if(std::abs(old_value-curmin) > diff_iter){
      
	    if(old_value != curmin ){
		 	vertex.set_data(curmin);
		//	vertex.modified=true;
		}else{
			vertex.set_data(curmin);
			vertex.modified=false;
		}	
      
        /** 
         * Broadcast new label to neighbors by writing the value
         * to the incident edges.
         * Note: on first iteration, write only to out-edges to avoid
         * overwriting data (this is kind of a subtle point)
         */
        vid_t label = vertex.get_data();

        if (gcontext.iteration > 0) {
            for(int i=0; i < vertex.num_edges(); i++) {
                if (label < vertex.edge(i)->get_data()) {
                    vertex.edge(i)->set_data(label);
                    /* Schedule neighbor for update */
                   if (scheduler) gcontext.scheduler->add_task(vertex.edge(i)->vertex_id(), true);
                    converged = false;
			//	logstream(LOG_INFO)<<"gcontext.iteration==================="<<gcontext.iteration<<std::endl;
					vertex.modified=true;
                }
            }
        } else if (gcontext.iteration == 0) {
            for(int i=0; i < vertex.num_outedges(); i++) {
                vertex.outedge(i)->set_data(label);
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
        }
    }
    
    /**
     * Called before an execution interval is started.
     */
    void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &ginfo) {        
    }
    
    /**
     * Called after an execution interval has finished.
     */
    void after_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &ginfo) {        
    }
    
};

struct ConnectedComponentsProgramInmem : public GraphChiProgram<VertexDataType, EdgeDataType> {
    
    bool converged;
	bool interval_converged;
   	std::vector< VertexDataType > wcc; 
	ConnectedComponentsProgramInmem(int nvertices):wcc(nvertices, Max_Cid){};	
    /**
     *  Vertex update function.
     *  On first iteration ,each vertex chooses a label = the vertex id.
     *  On subsequent iterations, each vertex chooses the minimum of the neighbor's
     *  label (and itself). 
     */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {
        
       if (scheduler) gcontext.scheduler->remove_tasks(vertex.id(), vertex.id());
        

        if (gcontext.iteration == 0) {
            //vertex.set_data(vertex.id());
			//std::cout<<"vid="<<vertex.id()<<"\t cid="<<wcc[vertex.id()]<<std::endl;
			//assert(wcc[vertex.id()] == Max_Cid);
			wcc[vertex.id()] = vertex.id();
			vertex.set_data(wcc[vertex.id()]);	
            if (scheduler)  gcontext.scheduler->add_task(vertex.id());
        }else{
        
        /* On subsequent iterations, find the minimum label of my neighbors */
        //vid_t curmin = vertex.get_data();
        vid_t curmin = wcc[vertex.id()];


        vid_t old_value =curmin ; //vertex.get_data();
        for(int i=0; i < vertex.num_edges(); i++) {
            //vid_t nblabel = vertex.edge(i)->get_data();
            vid_t nblabel = wcc[vertex.edge(i)->vertex_id()];
		
            curmin = std::min(nblabel, curmin); 
        }
   
	    if(old_value > curmin ){
		// 	vertex.set_data(curmin);
			wcc[vertex.id()] = curmin;	
			//vertex.modified=true;
			//vertex.set_data(wcc[vertex.id()]);	
			converged = false;
			interval_converged = false;
		}/*else{
		//	vertex.set_data(curmin);
			vertex.modified=false;
		}*/	
      
        /** 
         * Broadcast new label to neighbors by writing the value
         * to the incident edges.
         * Note: on first iteration, write only to out-edges to avoid
         * overwriting data (this is kind of a subtle point)
         */
        //vid_t label = vertex.get_data();
        //vid_t label = wcc[vertex.id()];
		
		vertex.set_data(wcc[vertex.id()]);	
        if (gcontext.iteration > 0) {
			/*
            for(int i=0; i < vertex.num_edges(); i++) {
                if (label < wcc[vertex.edge(i)->vertex_id()]) {
					wcc[vertex.edge(i)->vertex_id()] = label;
                    //vertex.edge(i)->set_data(label);
                   if (scheduler) gcontext.scheduler->add_task(vertex.edge(i)->vertex_id(), true);
                    converged = false;
					interval_converged = false;
					//vertex.modified=true;
                }
            }
			*/
        } else if (gcontext.iteration == 0) {
			/*	
            for(int i=0; i < vertex.num_outedges(); i++) {
                //vertex.outedge(i)->set_data(label);
                //wcc[vertex.outedge(i)->vertex_id] = label;
                wcc[vertex.outedge(i)->vertex_id()] = label;
            }
			*/
        }
		/*
		if(gcontext.iteration == gcontext.num_iterations-1){
			vertex.set_data(wcc[vertex.id()]);	
		}
		*/
			//vertex.set_data(wcc[vertex.id()]);	
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
			for(i=0; i<wcc.size(); i++){
				fprintf(fp, "%u\t%u\n", i, wcc[i]);
			}
			fclose(fp);
			*/
        }
    }
   	bool repeat_updates(graphchi_context &gcontext){
		if(gcontext.iteration == 0)
			return false;
		else{
			std::cout<<"repeat_updates isconverged="<<interval_converged<<"\t iteration="<<gcontext.iteration<<std::endl;
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
				
			//assert(vertex.get_data() == wcc[vertex.id()]);	
			vertex.set_data(wcc[vertex.id()]);	
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
    metrics m("connected-components");
    
    /* Basic arguments for application */
    std::string filename = get_option_string("file");  // Base filename
    int niters           = get_option_int("niters", 200); // Number of iterations (max)
    scheduler            = get_option_int("scheduler", false);
 //   diff_iter			=get_option_int("diff",0);
    /* Process input file - if not already preprocessed */
    int nshards             = (int) convert_if_notexists<EdgeDataType>(filename, get_option_string("nshards", "auto"));
    
        /* Run */
	graphchi_engine<VertexDataType, EdgeDataType> engine(filename, nshards, scheduler, m); 

	bool inmemmode = engine.num_vertices()*sizeof(VertexDataType) < (size_t)engine.get_membudget_mb()*1024L*1024L;
	if(inmemmode){
		logstream(LOG_INFO)<<"\nRunning WCC by holding vertices in-memory mode!!!!"<<std::endl;
		engine.set_modifies_outedges(false);
		//engine.set_disable_outedges(true);
		//engine.set_modifies_inedges(true);
		engine.set_only_adjacency(true);	
		//wcc(engine.num_vertices(), Max_Cid);
		//wcc = (vid_t*)malloc(sizeof(vid_t)*(engine.num_vertices()+1));
		//assert(wcc != NULL);
		//memset(wcc,Max_Cid,engine.num_vertices());
		ConnectedComponentsProgramInmem program(engine.num_vertices());
		//ConnectedComponentsProgramInmem program;
        engine.run(program, niters);
		//ConnectedComponentsProgramCount prog;
		//engine.run(prog, 1);
		//free(wcc);
	}else{
		assert(false);
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

