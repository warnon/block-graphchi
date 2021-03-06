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
 * Simple pagerank implementation. Uses the basic vertex-based API for
 * demonstration purposes. A faster implementation uses the functional API,
 * "pagerank_functional".
 */

#include <string>
#include <fstream>
#include <cmath>

#define GRAPHCHI_DISABLE_COMPRESSION


#include "graphchi_basic_includes.hpp"
#include "util/toplist.hpp"

using namespace graphchi;
 
#define THRESHOLD 1e-2    
#define RANDOMRESETPROB 0.15


typedef float VertexDataType;
typedef float EdgeDataType;

float epsilon = 0.001;
int* array = NULL;
int repeat_count = 0;
bool scheduler = false;

struct PagerankProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {
   	bool converged;
	bool interval_converged; 
    /**
      * Called before an iteration starts. Not implemented.
      */
    void before_iteration(int iteration, graphchi_context &info) {
			converged = iteration > 0;
    }
    
    /**
      * Called after an iteration has finished. Not implemented.
      */
    void after_iteration(int iteration, graphchi_context &ginfo) {
		/*
		if(converged){
			logstream(LOG_INFO)<<"converged!"<<std::endl;
			ginfo.set_last_iteration(iteration);
		}
		*/
    }
    
    /**
      * Called before an execution interval is started. Not implemented.
      */
    void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &ginfo) {        
			interval_converged = true;
    }
    
  	bool repeat_updates(graphchi_context &gcontext){
		return false;
		/*
		if(gcontext.iteration == 0)
			return false;
		else{
			if(repeat_count % 100 == 0 || interval_converged)
			logstream(LOG_INFO)<<"exec_interval="<<gcontext.exec_interval<<"\trepeated_updates isconverged="<<interval_converged<<
					"\t iteration="<<gcontext.iteration<<"\t repeat_count="<<repeat_count<<std::endl;
			repeat_count = interval_converged ? 0 : repeat_count+1; 
			return !interval_converged;
		}
		*/
	}  
    /**
      * Incremental Pagerank update function.
      */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &v, graphchi_context &ginfo) {

		if(scheduler) ginfo.scheduler->remove_tasks_now(v.id(), v.id());			

        if (ginfo.iteration == 0) {
            /* On first iteration, initialize vertex and out-edges. 
               The initialization is important,
               because on every run, GraphChi will modify the data in the edges on disk. 
             */
            v.set_data(RANDOMRESETPROB); 
            for(int i=0; i < v.num_outedges(); i++) {
                graphchi_edge<float> * edge = v.outedge(i);
                //edge->set_data(1.0 / v.num_outedges());
                edge->set_data(RANDOMRESETPROB*(1-RANDOMRESETPROB)/ v.num_outedges());
            }
			if(scheduler) ginfo.scheduler->add_task(v.id(), false);
        } else {
        	float sum=0;
			float old_value = v.get_data();
            /* Compute the sum of neighbors' weighted pageranks by
               reading from the in-edges. */
            for(int i=0; i < v.num_inedges(); i++) {
                float val = v.inedge(i)->get_data();
                sum += val;                    
				if(val != .0f)
                	v.inedge(i)->set_data(.0f);
            }
            /* Compute my pagerank */
            //float pagerank = RANDOMRESETPROB + (1 - RANDOMRESETPROB) * sum;
            //float pagerank = sum;
            
            /* Write my pagerank divided by the number of out-edges to
               each of my out-edges. */
            if (v.num_outedges() > 0 && sum > 0) {
                float pagerankcont = (1 - RANDOMRESETPROB)*sum / v.num_outedges();
                for(int i=0; i < v.num_outedges(); i++) {
                    graphchi_edge<float> * edge = v.outedge(i);
					float edata = edge->get_data();
                    edge->set_data(pagerankcont + edata);
					if(scheduler){
						if(edge->vertexid >= ginfo.interval_st){
							ginfo.scheduler->add_task(edge->vertexid, true);
						}else{
							ginfo.scheduler->add_task(edge->vertexid, false);
						}
					}
                }
            }
            /* Keep track of the progression of the computation.
               GraphChi engine writes a file filename.deltalog. */
            //ginfo.log_change(std::abs(pagerank - v.get_data()));
            /* Set my new pagerank as the vertex value */
            v.set_data(old_value+sum); 
			if(sum > epsilon){
				converged = false;
				interval_converged = false;
			}
        }
    }
};

/**
  * Faster version of pagerank which holds vertices in memory. Used only if the number
  * of vertices is small enough.
  */
struct PagerankProgramInmem : public GraphChiProgram<VertexDataType, EdgeDataType> {
    
    std::vector<EdgeDataType> pr;
    PagerankProgramInmem(int nvertices) :   pr(nvertices, RANDOMRESETPROB) {}
    
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &v, graphchi_context &ginfo) {
        if (ginfo.iteration > 0) {
            float sum=0;
            for(int i=0; i < v.num_inedges(); i++) {
              sum += pr[v.inedge(i)->vertexid];
            }
            if (v.outc > 0) {
                pr[v.id()] = (RANDOMRESETPROB + (1 - RANDOMRESETPROB) * sum) / v.outc;
            } else {
                pr[v.id()] = (RANDOMRESETPROB + (1 - RANDOMRESETPROB) * sum);
            }
        } else if (ginfo.iteration == 0) {
            if (v.outc > 0) pr[v.id()] = 1.0f / v.outc;
        }
        if (ginfo.iteration == ginfo.num_iterations - 1) {
            /* On last iteration, multiply pr by degree and store the result */
            v.set_data(v.outc > 0 ? pr[v.id()] * v.outc : pr[v.id()]);
        }
    }
    
};

int main(int argc, const char ** argv) {
    graphchi_init(argc, argv);
    metrics m("pagerank");
    global_logger().set_log_level(LOG_DEBUG);

    /* Parameters */
    std::string filename    = get_option_string("file"); // Base filename
    int niters              = get_option_int("niters", 100000);
    scheduler         	 	= get_option_int("scheduler", false);                    // Non-dynamic version of pagerank.
    int ntop                = get_option_int("top", 20);
    epsilon	= get_option_float("epsilon", 0.001);
    /* Process input file - if not already preprocessed */
    int nshards             = convert_if_notexists<EdgeDataType>(filename, get_option_string("nshards", "auto"));
	assert(0 != nshards);
	//array = (int*)malloc(sizeof(int)*nshards);
    /* Run */
    graphchi_engine<float, float> engine(filename, nshards, scheduler, m); 
    //engine.set_modifies_inedges(false); // Improves I/O performance.
    
    bool inmemmode = false;//engine.num_vertices() * sizeof(EdgeDataType) < (size_t)engine.get_membudget_mb() * 1024L * 1024L;
    if (inmemmode) {
        logstream(LOG_INFO) << "Running Pagerank by holding vertices in-memory mode!" << std::endl;
        engine.set_modifies_outedges(false);
        engine.set_disable_outedges(true);
        engine.set_only_adjacency(true);
        PagerankProgramInmem program(engine.num_vertices());
        engine.run(program, niters);
    } else {
        PagerankProgram program;
        engine.run(program, niters);
    }
    
    /* Output top ranked vertices */
    std::vector< vertex_value<float> > top = get_top_vertices<float>(filename, ntop);
    std::cout << "Print top " << ntop << " vertices:" << std::endl;
    for(int i=0; i < (int)top.size(); i++) {
        std::cout << (i+1) << ". " << top[i].vertex << "\t" << top[i].value << std::endl;
    }
    
    metrics_report(m);    
    return 0;
}

