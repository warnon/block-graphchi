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

#include <time.h>
#include <stdlib.h>
#include <sstream>

#define GRAPHCHI_DISABLE_COMPRESSION


#include "graphchi_basic_includes.hpp"
#include "util/toplist.hpp"

using namespace graphchi;
 
#define THRESHOLD 1e-2    
#define RANDOMRESETPROB 0.15


typedef float VertexDataType;
typedef float EdgeDataType;

struct PR_Count{
int count;
int tasks;
PR_Count(){
	count = tasks = 0;
}
};

float* vdata_array = NULL;
float epsilon = 0.001;
//int* array = NULL;
PR_Count* array = NULL;
int repeat_count = 0;
int max_repeat = 999999999;
bool scheduler = false;
bool num_tasks_print = false;
vid_t sum_tasks= 0;
vid_t curr_sum = 0;
graphchi_engine<float, float>* engine_ptr; 
//vertex_data_store<VertexDataType>* vdata_handle;

struct PagerankProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {
   	bool converged;
	bool interval_converged; 
    /**
      * Called before an iteration starts. Not implemented.
      */
    void before_iteration(int iteration, graphchi_context &ginfo) {
			/*
			if(ginfo.iteration == 1){
				ginfo.scheduler->add_task_to_all();	
			}
			*/
			converged = iteration > 0;
    }
    
    /**
      * Called after an iteration has finished. Not implemented.
      */
    void after_iteration(int iteration, graphchi_context &ginfo) {
		if(converged){
			logstream(LOG_INFO)<<"converged!"<<std::endl;
			ginfo.set_last_iteration(iteration);
		}
    }
    
    /**
      * Called before an execution interval is started. Not implemented.
      */
	void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &ginfo) {        
		interval_converged = true;
		if(num_tasks_print && scheduler){
			//vid_t sum = ginfo.scheduler->total_tasks(window_st, window_en);
			curr_sum = ginfo.scheduler->total_tasks(window_st, window_en);
			logstream(LOG_INFO)<<"num of vertices being scheduled="<<curr_sum<<"/"<<window_en-window_st+1<<std::endl;		
			sum_tasks += curr_sum;
			array[ginfo.exec_interval].count++;
			if(scheduler) array[ginfo.exec_interval].tasks += curr_sum;	
			else	array[ginfo.exec_interval].tasks +=	window_en - window_st + 1;  
		}
	}
    
  	bool repeat_updates(graphchi_context &gcontext){
		if(gcontext.iteration == 0)
			return false;
		else{
			//if(repeat_count % 100 == 0 || interval_converged)
			bool count = (max_repeat < 999999999 && repeat_count == max_repeat);
			if(repeat_count % 100 == 0 || interval_converged || count)
			logstream(LOG_INFO)<<"interval="<<gcontext.exec_interval<<" iter="<<gcontext.iteration
				<<"\trepeat_count= "<<repeat_count<<std::endl;
			repeat_count = (interval_converged || count) ? 0 : repeat_count+1; 
			//repeat_count = (interval_converged) ? 0 : repeat_count+1; 

			if(num_tasks_print){
				sum_tasks = interval_converged ? 0 : sum_tasks; 
			}
			
				
			//return !interval_converged ;
			return !(interval_converged || count);
		}
	}  
    /**
      * Pagerank update function.
      */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &v, graphchi_context &ginfo) {
		//if(scheduler) ginfo.scheduler->remove_tasks_now(v.id(), v.id());
        //float sum=0;
        if (ginfo.iteration == 0) {
            /* On first iteration, initialize vertex and out-edges. 
               The initialization is important,
               because on every run, GraphChi will modify the data in the edges on disk. 
             */
            for(int i=0; i < v.num_outedges(); i++) {
                graphchi_edge<float> * edge = v.outedge(i);
                //edge->set_data( RANDOMRESETPROB/v.num_outedges());
                edge->set_data(1.0 / v.num_outedges());
            }
            v.set_data(RANDOMRESETPROB); 
			//schedule this vertex for next iteration
			if(scheduler) ginfo.scheduler->add_task(v.id(), false);
        } else {

			if(scheduler) ginfo.scheduler->remove_tasks_now(v.id(), v.id());
        	float sum=0;
			float old_value = v.get_data();
            /* Compute the sum of neighbors' weighted pageranks by
               reading from the in-edges. */
            for(int i=0; i < v.num_inedges(); i++) {
                sum += v.inedge(i)->get_data();
                //sum += val;                    
            }
            /* Compute my pagerank */
            float pagerank = RANDOMRESETPROB + (1 - RANDOMRESETPROB) * sum;
            
            /* Write my pagerank divided by the number of out-edges to
               each of my out-edges. */
			
            if (v.num_outedges() > 0) {
                float pagerankcont = pagerank / v.num_outedges();
                for(int i=0; i < v.num_outedges(); i++) {
                    graphchi_edge<float> * edge = v.outedge(i);
					edge->set_data(pagerankcont);
					if(scheduler){
						if(edge->vertexid >= ginfo.interval_st){
							if(!ginfo.scheduler->is_scheduled(edge->vertexid, true))
								ginfo.scheduler->add_task(edge->vertexid, true);
						}else{
							if(!ginfo.scheduler->is_scheduled(edge->vertexid, false))
								ginfo.scheduler->add_task(edge->vertexid, false);
						}	
					}
				}
            }
            /* Keep track of the progression of the computation.
               GraphChi engine writes a file filename.deltalog. */
            ginfo.log_change(std::abs(pagerank - v.get_data()));
            
            /* Set my new pagerank as the vertex value */
            v.set_data(pagerank); 
			if(std::abs(pagerank-old_value) > epsilon){
				if(converged)
					converged = false;
				if(interval_converged)
					interval_converged = false;
				/*if (v.num_outedges() > 0) {
					float pagerankcont = pagerank / v.num_outedges();
					for(int i=0; i < v.num_outedges(); i++) {
						graphchi_edge<float> * edge = v.outedge(i);
						edge->set_data(pagerankcont);
					}
				}
				v.set_data(pagerank);*/ 
			}
        }
    }
    
};


struct initPagerankProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {
   	bool converged;
	bool interval_converged; 
        
    /**
      * Called after an iteration has finished. Not implemented.
      */
    void after_iteration(int iteration, graphchi_context &ginfo) {
		ginfo.set_last_iteration(iteration);
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
		/*
		interval_converged = true;
		if(num_tasks_print && scheduler){
			//vid_t sum = ginfo.scheduler->total_tasks(window_st, window_en);
			curr_sum = ginfo.scheduler->total_tasks(window_st, window_en);
			logstream(LOG_INFO)<<"num of vertices being scheduled="<<curr_sum<<"/"<<window_en-window_st+1<<std::endl;		
			sum_tasks += curr_sum;
			array[ginfo.exec_interval].count++;
			if(scheduler) array[ginfo.exec_interval].tasks += curr_sum;	
			else	array[ginfo.exec_interval].tasks +=	window_en - window_st + 1;  
		}
		*/
	}
    
  	bool repeat_updates(graphchi_context &gcontext){
		return false;
	}  
    /**
      * Pagerank update function.
      */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &v, graphchi_context &ginfo) {
		//if(scheduler) ginfo.scheduler->remove_tasks_now(v.id(), v.id());
        //float sum=0;
        if (ginfo.iteration == 0) {
			/*
            for(int i=0; i < v.num_outedges(); i++) {
                graphchi_edge<float> * edge = v.outedge(i);
                //edge->set_data( RANDOMRESETPROB/v.num_outedges());
                edge->set_data(1.0 / v.num_outedges());
            }
			*/
			vdata_array[(int)v.id()] = v.get_data();
            //v.set_data(RANDOMRESETPROB); 
			//schedule this vertex for next iteration
        } 
    }
    
};

vid_t active_count = 0;
float l1_error = 0;
std::vector<int> active_vertices;
std::vector<float> error_sum;
std::vector<double> pagerank_sum;
double pr_sum = 0;

struct comparePagerankProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {
   	bool converged;
	bool interval_converged; 
    /**
      * Called before an iteration starts. Not implemented.
      */
    void before_iteration(int iteration, graphchi_context &ginfo) {
			/*
			if(ginfo.iteration == 1){
				ginfo.scheduler->add_task_to_all();	
			}
			*/
			l1_error = 0;
			active_count = 0;	
			pr_sum = 0;
			converged = iteration > 0;
    }
    
    /**
      * Called after an iteration has finished. Not implemented.
      */
    void after_iteration(int iteration, graphchi_context &ginfo) {
		active_vertices.push_back(active_count);
		error_sum.push_back(l1_error);
		pagerank_sum.push_back(pr_sum);
		if(converged){
			logstream(LOG_INFO)<<"converged!"<<std::endl;
			ginfo.set_last_iteration(iteration);
		}
    }
    
    /**
      * Called before an execution interval is started. Not implemented.
      */
	void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &ginfo) {        
		interval_converged = true;
		if(num_tasks_print && scheduler){
			//vid_t sum = ginfo.scheduler->total_tasks(window_st, window_en);
			curr_sum = ginfo.scheduler->total_tasks(window_st, window_en);
			logstream(LOG_INFO)<<"num of vertices being scheduled="<<curr_sum<<"/"<<window_en-window_st+1<<std::endl;		
			sum_tasks += curr_sum;
			array[ginfo.exec_interval].count++;
			if(scheduler) array[ginfo.exec_interval].tasks += curr_sum;	
			else	array[ginfo.exec_interval].tasks +=	window_en - window_st + 1;  
		}
	}
    
  	bool repeat_updates(graphchi_context &gcontext){
		int interval_active = 0;
		float diff_sum = 0;

		if(gcontext.iteration == 0 || interval_converged){
			vertex_data_store<VertexDataType>* vdata_handle = engine_ptr->get_vertex_data_handle();
			int interval_len = gcontext.interval_en - gcontext.interval_st +1;
			for(int i=0; i<interval_len; i++){
				vid_t vid =(vid_t) i + gcontext.interval_st;
				assert(vdata_handle != NULL);
				float vdata = *(vdata_handle->vertex_data_ptr(vid));
				float diff = std::abs(vdata - vdata_array[vid]);
				diff_sum += diff;
				pr_sum += vdata;	
				if(diff > epsilon){
					interval_active++;	
				}
			}
			l1_error += diff_sum;
			active_count += interval_active;
		}

		if(gcontext.iteration == 0)
			return false;
		else{
						
			//if(repeat_count % 100 == 0 || interval_converged)
			bool count = (max_repeat < 999999999 && repeat_count == max_repeat);
			if(repeat_count % 100 == 0 || interval_converged || count)
			logstream(LOG_INFO)<<"interval="<<gcontext.exec_interval<<" iter="<<gcontext.iteration
				<<"\trepeat_count= "<<repeat_count <<" #active= "<<interval_active<<std::endl;
			repeat_count = (interval_converged || count) ? 0 : repeat_count+1; 
			//repeat_count = (interval_converged) ? 0 : repeat_count+1; 

			if(num_tasks_print){
				sum_tasks = interval_converged ? 0 : sum_tasks; 
			}
			//return !interval_converged ;
			return !(interval_converged || count);
		}
	}  
    /**
      * Pagerank update function.
      */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &v, graphchi_context &ginfo) {
		//if(scheduler) ginfo.scheduler->remove_tasks_now(v.id(), v.id());
        //float sum=0;
        if (ginfo.iteration == 0) {
            /* On first iteration, initialize vertex and out-edges. 
               The initialization is important,
               because on every run, GraphChi will modify the data in the edges on disk. 
             */
            for(int i=0; i < v.num_outedges(); i++) {
                graphchi_edge<float> * edge = v.outedge(i);
                //edge->set_data( RANDOMRESETPROB/v.num_outedges());
                edge->set_data(1.0 / v.num_outedges());
            }
            v.set_data(RANDOMRESETPROB); 
			//schedule this vertex for next iteration
			if(scheduler) ginfo.scheduler->add_task(v.id(), false);
        } else {

			if(scheduler) ginfo.scheduler->remove_tasks_now(v.id(), v.id());
        	float sum=0;
			float old_value = v.get_data();
            /* Compute the sum of neighbors' weighted pageranks by
               reading from the in-edges. */
            for(int i=0; i < v.num_inedges(); i++) {
                sum += v.inedge(i)->get_data();
                //sum += val;                    
            }
            /* Compute my pagerank */
            float pagerank = RANDOMRESETPROB + (1 - RANDOMRESETPROB) * sum;
            
            /* Write my pagerank divided by the number of out-edges to
               each of my out-edges. */
			
            if (v.num_outedges() > 0) {
                float pagerankcont = pagerank / v.num_outedges();
                for(int i=0; i < v.num_outedges(); i++) {
                    graphchi_edge<float> * edge = v.outedge(i);
					edge->set_data(pagerankcont);
					if(scheduler){
						if(edge->vertexid >= ginfo.interval_st){
							if(!ginfo.scheduler->is_scheduled(edge->vertexid, true))
								ginfo.scheduler->add_task(edge->vertexid, true);
						}else{
							if(!ginfo.scheduler->is_scheduled(edge->vertexid, false))
								ginfo.scheduler->add_task(edge->vertexid, false);
						}	
					}
				}
            }
            /* Keep track of the progression of the computation.
               GraphChi engine writes a file filename.deltalog. */
            ginfo.log_change(std::abs(pagerank - v.get_data()));
            
            /* Set my new pagerank as the vertex value */
            v.set_data(pagerank); 
			if(std::abs(pagerank-old_value) > epsilon){
				if(converged)
					converged = false;
				if(interval_converged)
					interval_converged = false;
				/*if (v.num_outedges() > 0) {
					float pagerankcont = pagerank / v.num_outedges();
					for(int i=0; i < v.num_outedges(); i++) {
						graphchi_edge<float> * edge = v.outedge(i);
						edge->set_data(pagerankcont);
					}
				}
				v.set_data(pagerank);*/ 
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
    metrics m("block-pagerank");
    global_logger().set_log_level(LOG_DEBUG);

    /* Parameters */
    std::string filename    = get_option_string("file"); // Base filename
    int niters              = get_option_int("niters", 100000);
    scheduler          		= get_option_int("scheduler", false);                    // Non-dynamic version of pagerank.
    int ntop                = get_option_int("top", 50);
    epsilon					= get_option_float("epsilon", 0.001);
	num_tasks_print			= get_option_int("print", false);
	max_repeat				= get_option_int("repeat", 999999999); // max number of repeat allowed for each subgraph
    /* Process input file - if not already preprocessed */
    int nshards             = convert_if_notexists<EdgeDataType>(filename, get_option_string("nshards", "auto"));
	assert(0 != nshards);
	//array = (int*)malloc(sizeof(int)*nshards);
	if(num_tasks_print){
		array = (PR_Count*)malloc(sizeof(PR_Count)*nshards);
		assert(array != NULL);
		for(int i=0; i<nshards; i++){
			array[i].count = 0;
			array[i].tasks = 0;	
		}
	}
    /* Run */
    graphchi_engine<float, float> engine(filename, nshards, scheduler, m); 
    engine.set_modifies_inedges(false); // Improves I/O performance.
    
	vdata_array = (float*) malloc(sizeof(float)*(engine.num_vertices()));
	assert(vdata_array != NULL);
	memset(vdata_array, 0, sizeof(float)*(engine.num_vertices()));
	/*	
	/////////
	graphchi_engine<float, float> engine4(filename, nshards, scheduler, m); 
    engine4.set_modifies_inedges(false); // Improves I/O performance.
	engine_ptr = &engine4;
	//vdata_handle = engine4.get_vertex_data_handle();
    comparePagerankProgram compprogram1;
    engine4.run(compprogram1, niters);
	//////
	*/
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
		//vdata_handle = engine.get_vertex_data_handle();
        engine.run(program, niters);
    }
   // init array 
    graphchi_engine<float, float> engine2(filename, nshards, scheduler, m); 
    engine2.set_modifies_inedges(false); // Improves I/O performance.
    initPagerankProgram initprogram;
    engine2.run(initprogram, niters);

	graphchi_engine<float, float> engine3(filename, nshards, scheduler, m); 
    engine3.set_modifies_inedges(false); // Improves I/O performance.
	engine_ptr = &engine3;

	//vdata_handle = engine3.get_vertex_data_handle();
    comparePagerankProgram compprogram;
    engine3.run(compprogram, niters);

    /* Output top ranked vertices */
    std::vector< vertex_value<float> > top = get_top_vertices<float>(filename, ntop);
    std::cout << "Print top " << ntop << " vertices:" << std::endl;
    for(int i=0; i < (int)top.size(); i++) {
        std::cout << (i+1) << ". " << top[i].vertex << "\t" << top[i].value << std::endl;
    }
    
    metrics_report(m);    
	free(vdata_array);
	//vdata_handle = NULL;
	if(	num_tasks_print ){
		FILE* pr_out = NULL;
		if(scheduler)
			pr_out = fopen((filename+".schd_pr_ana").c_str(), "w+");
		else
			pr_out = fopen((filename+".pr_ana").c_str(), "w+");
		assert(pr_out != NULL);
		int avg_count = 0;
		int avg_tasks = 0;
		int nvertices = engine.num_vertices();
		for(int i=0; i<nshards; i++){
			fprintf(pr_out, "%d\t%d\t%d\n", i, array[i].count, array[i].tasks);		
			avg_count += array[i].count;
			avg_tasks += array[i].tasks;
		}		
		fprintf(pr_out, "avg\t%.4f\t%.4f\n", (float)avg_count/nshards, (float)avg_tasks/nvertices);	
		fclose(pr_out);
		free(array);
		printf("avg\t%.4f\t%.4f\n", (float)avg_count/nshards, (float)avg_tasks/nvertices);	
	}
	/*
	srand(time(NULL));
	int rd = rand() % 10000;
	ostringstream ss;
	ss << rd;
	*/
	//std::string str = //std::to_string(rd);
	//FILE* pr_active = fopen((filename+ss.str()+".pr_active").c_str(), "w+");
	FILE* pr_active = fopen((filename+".blk_pr_active").c_str(), "w+");
	assert(pr_active != NULL);
	fprintf(pr_active, "iters,#active,L1_error,PR_Sum\n");	
	for(int i=0; i<active_vertices.size(); i++){
		fprintf(pr_active, "%d,%d,%f,%lf\n", i, active_vertices[i], error_sum[i], pagerank_sum[i]);	
	}			
	fclose(pr_active);
	pr_active = NULL;
    return 0;
}

