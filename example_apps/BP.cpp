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

#define GRAPHCHI_DISABLE_COMPRESSION


#include "graphchi_basic_includes.hpp"
#include "util/toplist.hpp"
#include "util/UnaryFactor.hpp"

using namespace graphchi;
 
#define THRESHOLD 1e-3
#define RANDOMRESETPROB 0.15

//const int ColorLevel = 12;
int DeltaColor = 5;
double Bound = 1e-3;
double Damping = 0.1;
double Sigma = 4.0;
double Lambda = 5.0;

struct VertexInfo{
	UnaryFactor belief;
	UnaryFactor potential;
	//gray value
	int obs;	
	/*
	VertexInfo(){
		
	}*/	
	VertexInfo(const VertexInfo& vdata){
		belief = vdata.belief;
		potential = vdata.potential;
	}
	VertexInfo(int grayvalue){
		obs = grayvalue;
	}
	void VertexInit(){
		belief.resize(ColorLevel*2+1);	
		potential.resize(ColorLevel*2+1);	
		belief.uniform();
		potential.uniform();
		belief.normalize();
		potential.normalize();
		double SigmaSq = Sigma*Sigma;
		for(int i=0; i<=ColorLevel; i++){
			potential.data[i] =	-(i-ColorLevel)*(i-ColorLevel)/(2.0*SigmaSq);	 	
		}
		potential.normalize();
	}	
};

struct EdgeInfo{
	UnaryFactor msg;		
	int largeObs;	
	int smallObs;
	/*
	EdgeInfo(){
		msg.resize(ColorLevel*2+1);	
	}	
	*/
	EdgeInfo(){}
	EdgeInfo(const EdgeInfo& edata){
		largeObs = edata.largeObs;
		smallObs = edata.smallObs;
		msg = edata.msg;	
	}
	void EdgeInit(){
		msg.uniform();
		msg.normalize();	
	}

	int& my_label(vid_t myid, vid_t nbid){
		assert(myid != 0xffffffff && nbid != 0xffffffff);
		if(myid < nbid)	 return smallObs;	
		else return largeObs;
	}		
	
	int& nb_label(vid_t myid, vid_t nbid){
		assert(myid != 0xffffffff && nbid != 0xffffffff);
		if(myid < nbid) return largeObs;
		else return smallObs;
	}
	/*
	friend std::ostream& operator << (std::ostream& ost, EdgeInfo edata){
		ost<<edata.msg<<edata.largeObs<<edata.smallObs;
		return ost;	
	}
	*/
};

typedef VertexInfo VertexDataType;
typedef EdgeInfo EdgeDataType;


float epsilon = 0.001;
//int* array = NULL;
//PR_Count* array = NULL;
int repeat_count = 0;
//int max_repeat = 999999999;
bool scheduler = false;
//bool num_tasks_print = false;
vid_t sum_tasks= 0;
vid_t curr_sum = 0;

struct BPProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {
   	bool converged;
	bool interval_converged; 
    /**
      * Called before an iteration starts. Not implemented.
      */
    void before_iteration(int iteration, graphchi_context &ginfo) {
			if(iteration == 0)
				srand(time(NULL));	
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
	}
    
  	bool repeat_updates(graphchi_context &gcontext){
		if(gcontext.iteration == 0)
			return false;
		else{
			//if(repeat_count % 100 == 0 || interval_converged)
			if(repeat_count % 20 == 0 || interval_converged)
			logstream(LOG_INFO)<<"interval="<<gcontext.exec_interval<<" iter="<<gcontext.iteration
				<<"\trepeat_count= "<<repeat_count<<std::endl;
			repeat_count = (interval_converged) ? 0 : repeat_count+1; 
			return !(interval_converged);
		}
	}  
    /**
      * belief propagation update function.
      */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &ginfo) {
		//if(scheduler) ginfo.scheduler->remove_tasks_now(v.id(), v.id());
        //float sum=0;
        if (ginfo.iteration == 0) {
			//randomly select a grayvalue ranging from 0 to 256
			int gray_value = rand() % 256;
			VertexDataType vdata = vertex.get_data();
			vdata.VertexInit();
			vdata.obs = gray_value;	
			vertex.set_data(vdata);
			for(int i=0; i < vertex.num_edges(); i++) {
				graphchi_edge<EdgeDataType> * edge;
				edge = vertex.edge(i);	
				EdgeDataType edata = edge->get_data();
				edata.EdgeInit();	
				edata.my_label(vertex.id(), edge->vertex_id()) = gray_value;
				edge->set_data(edata);
				/*
				   if(i<vertex.num_inedges()){
				   edge = vertex.inedge(i);	
				   EdgeDataType edata = edge->get_data();
				   edata.EdgeInit();	
				   edata.dstObs = gray_value;
				   edge->set_data(edata);
				   }else{
				   int idx = i-vertex.num_inedges();
				   edge = vertex.outedge(i);
				   EdgeDataType edata = edge->get_data();		
				   edata.EdgeInit();	
				   edata.srcObs = gray_value;
				   edge->set_data(edata);
				   }
				   */
            }
            //v.set_data(RANDOMRESETPROB); 
			//schedule this vertex for next iteration
				
			if(scheduler) ginfo.scheduler->add_task(vertex.id(), false);
        } else {
			if(scheduler) ginfo.scheduler->remove_tasks_now(vertex.id(), vertex.id());
        	//float sum=0;
			//float old_value = v.get_data();
            /* Compute the sum of neighbors' weighted pageranks by
               reading from the in-edges. */
			VertexDataType vdata = vertex.get_data();	
			vdata.belief.copy(vdata.potential);
            for(int i=0; i < vertex.num_edges(); i++) {
				vdata.belief.times((vertex.edge(i)->get_data()).msg);	
                //sum += v.inedge(i)->get_data();
                //sum += val;                    
            }
			vdata.belief.normalize();
            /* Compute my pagerank */
            //float pagerank = RANDOMRESETPROB + (1 - RANDOMRESETPROB) * sum;
            
            /* Write my pagerank divided by the number of out-edges to
               each of my out-edges. */
			double max_residual = .0;	
            if (vertex.num_edges() > 0) {
				UnaryFactor cavity, tmpFactor; 
				tmpFactor.resize(ColorLevel*2+1);
                //float pagerankcont = pagerank / v.num_outedges();
                for(int i=0; i < vertex.num_edges(); i++) {
					cavity = vdata.belief;
                    graphchi_edge<EdgeDataType> * edge = vertex.edge(i);
					EdgeDataType edata = edge->get_data();
					cavity.divide(edata.msg);
					cavity.normalize();	
					int myobs = edata.my_label(vertex.id(), edge->vertex_id());		
					int nbobs = edata.nb_label(vertex.id(), edge->vertex_id());		
					tmpFactor.convolve(myobs, nbobs, DeltaColor, ColorLevel, Lambda, cavity);
					tmpFactor.normalize();			
					tmpFactor.damp(edata.msg, Damping);
					//for termination check
					double residual = tmpFactor.residual(edata.msg);
					max_residual = max_residual < residual ? residual : max_residual; 
					edata.msg = tmpFactor;	
					edge->set_data(edata);	
					//edge->set_data(pagerankcont);
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
            //ginfo.log_change(std::abs(pagerank - v.get_data()));
            ginfo.log_change(max_residual);
            
            /* Set my new pagerank as the vertex value */
            //v.set_data(pagerank); 
			if(max_residual > Bound){
				if(converged)
					converged = false;
				if(interval_converged)
					interval_converged = false;
			}
        }
    }
    
};

/**
  * Faster version of pagerank which holds vertices in memory. Used only if the number
  * of vertices is small enough.
  */
/*
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
            v.set_data(v.outc > 0 ? pr[v.id()] * v.outc : pr[v.id()]);
        }
    }
    
};
*/
int main(int argc, const char ** argv) {
    graphchi_init(argc, argv);
    metrics m("block-BP");
    global_logger().set_log_level(LOG_DEBUG);

    /* Parameters */
    std::string filename    = get_option_string("file"); // Base filename
    int niters              = get_option_int("niters", 100000);
    scheduler          		= get_option_int("scheduler", false);                    // Non-dynamic version of pagerank.
    //int ntop                = get_option_int("top", 50);
    epsilon					= get_option_float("epsilon", 0.001);
	//num_tasks_print			= get_option_int("print", false);
	//max_repeat				= get_option_int("repeat", 999999999); // max number of repeat allowed for each subgraph
    /* Process input file - if not already preprocessed */
    int nshards             = convert_if_notexists<EdgeDataType>(filename, get_option_string("nshards", "auto"));
	assert(0 != nshards);
	//array = (int*)malloc(sizeof(int)*nshards);

    /* Run */
    graphchi_engine<VertexDataType, EdgeDataType> engine(filename, nshards, scheduler, m); 
    engine.set_modifies_inedges(false); // Improves I/O performance.
    
    bool inmemmode = false;//engine.num_vertices() * sizeof(EdgeDataType) < (size_t)engine.get_membudget_mb() * 1024L * 1024L;
    if (inmemmode) {
		/*
        logstream(LOG_INFO) << "Running Pagerank by holding vertices in-memory mode!" << std::endl;
        engine.set_modifies_outedges(false);
        engine.set_disable_outedges(true);
        engine.set_only_adjacency(true);
        PagerankProgramInmem program(engine.num_vertices());
        engine.run(program, niters);
		*/
    } else {
        BPProgram program;
        engine.run(program, niters);
    }
    
	/*
    std::vector< vertex_value<float> > top = get_top_vertices<float>(filename, ntop);
    std::cout << "Print top " << ntop << " vertices:" << std::endl;
    for(int i=0; i < (int)top.size(); i++) {
        std::cout << (i+1) << ". " << top[i].vertex << "\t" << top[i].value << std::endl;
    }
   	*/ 
	metrics_report(m);    
	/*
	if(	num_tasks_print ){
		array = (PR_Count*)malloc(sizeof(PR_Count)*nshards);
		assert(array != NULL);
		for(int i=0; i<nshards; i++){
			array[i].count = 0;
			array[i].tasks = 0;	
		}
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
	*/
    return 0;
}

