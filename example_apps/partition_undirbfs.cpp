
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
#include "util/undirbfsanalysis.hpp"

using namespace graphchi;

int         iterationcount = 0;
bool        scheduler = false;
vid_t Root = 0;
const int top = 20;
//vid_t Max_Level =(vid_t)-1; 
vid_t Max_Level = 1000000; 

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
struct Vertexinfo{
	vid_t ccid;
	vid_t level; 
	/*
	vid_t degree;
	Vertexinfo():level(0), degree(0){};
	Vertexinfo(vid_t lv):level(lv), degree(0){};
	Vertexinfo(vid_t lv, vid_t deg):level(lv), degree(deg){};
	*/	
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


struct bidirectional_label {
    vid_t smaller_one;
    vid_t larger_one;

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
typedef Vertexinfo VertexDataType;       // vid_t is the vertex id type
typedef bidirectional_label EdgeDataType;

struct level_pair{
	vid_t label;
	vid_t count;
	// constructor
	level_pair():label(0), count(0){};	
	level_pair(vid_t lb, vid_t size):label(lb), count(size){};
	level_pair(vid_t lb):label(lb), count(0){};
};

static uint32_t repeat_count = 0;
std::vector<level_pair>	cc_size; 
//int top = 0;
/*
vid_t getNewId(int level){
	assert(level < );
	lock.lock();
		
}
*/
struct CCProgram: public GraphChiProgram<VertexDataType, EdgeDataType> {
//struct CCProgram: public GraphChiProgram<VertexDataType, EdgeDataType> {
	bool converged;
	bool interval_converged;
	void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {
		//if (scheduler) gcontext.scheduler->remove_tasks(vertex.id(), vertex.id());
		Vertexinfo vdata = vertex.get_data();
		if (gcontext.iteration == 0) {
			vid_t min_nid = vertex.id();	
			//vertex.set_data(VertexDataType(min_nid, (vid_t)vertex.num_edges()));	
			//vdata.ccid = min_nid;	
			//vertex.set_data(vdata);	
			for(int i=0; i<vertex.num_edges(); i++){
				min_nid = std::min(min_nid, vertex.edge(i)->vertex_id());	
			}	
			//vertex.set_data(VertexDataType(min_nid, (vid_t)vertex.num_edges()));	
			vdata.ccid = min_nid;	
			vertex.set_data(vdata);	
			for(int i=0; i<vertex.num_edges(); i++){
				//vertex.outedge(i)->set_data(vertex.id());
				bidirectional_label edata = vertex.edge(i)->get_data();;
				edata.my_label(vertex.id(), vertex.edge(i)->vertex_id()) = min_nid;
				vertex.edge(i)->set_data(edata);
				//vertex.outedge(i)->set_data(min_nid);
				if(scheduler)
					gcontext.scheduler->add_task(vertex.edge(i)->vertex_id());
			}
		}else{
			//Vertexinfo vdata = vertex.get_data();
			//vid_t curmin = vdata.level;
			vid_t curmin = vdata.ccid;
			vid_t old_value =curmin ; //vertex.get_data();
			for(int i=0; i < vertex.num_edges(); i++) {
				//vid_t nblabel = vertex.edge(i)->get_data();
				bidirectional_label edata = vertex.edge(i)->get_data();
				vid_t nblabel = edata.neighbor_label(vertex.id(), vertex.edge(i)->vertex_id());
				curmin = std::min(nblabel, curmin); 
			}
			if(old_value > curmin ){
				vdata.ccid = curmin;	
				vertex.set_data(vdata);
				//vertex.modified=true;
				converged = false;
				interval_converged = false;
				for(int i=0; i<vertex.num_edges(); i++){
					bidirectional_label edata = vertex.edge(i)->get_data();
					edata.my_label(vertex.id(), vertex.edge(i)->vertex_id()) = curmin;
					vertex.edge(i)->set_data(edata);		
					if(scheduler)
						gcontext.scheduler->add_task(vertex.edge(i)->vertex_id());
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

/**
 * GraphChi programs need to subclass GraphChiProgram<vertex-type, edge-type> 
 * class. The main logic is usually in the update function.
 */
struct BFSProgram: public GraphChiProgram<VertexDataType, EdgeDataType> {
	bool converged;
	bool interval_converged;
	void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {
		//if (scheduler) gcontext.scheduler->remove_tasks(vertex.id(), vertex.id());
		Vertexinfo vdata = vertex.get_data();  
		if (gcontext.iteration == 0){
			if(vertex.id() == vdata.ccid){
				if(vertex.num_edges() > 0){
				//source vertex in each CC
					vdata.level = 0;
					vertex.set_data(vdata);	
					for(int i=0; i<vertex.num_edges(); i++){
						bidirectional_label edata = vertex.edge(i)->get_data();
						edata.my_label(vertex.id(), vertex.edge(i)->vertex_id()) = 0;
						vertex.edge(i)->set_data(edata);
						//vertex.edge(i)->set_data(Root);		
						if(scheduler)
							gcontext.scheduler->add_task(vertex.edge(i)->vertex_id());	
					}
				}else{
					// isolated vertices
					//vdata.level = Max_Level;	
					vdata.level = 0xffffffff;	
					vertex.set_data(vdata);
				}
			}else{
				vdata.level = Max_Level;
				vertex.set_data(vdata);	
				//vertex.set_data(Vertexinfo(Max_Level));	
				for(int i=0; i<vertex.num_edges(); i++){
					bidirectional_label edata = vertex.edge(i)->get_data();
					edata.my_label(vertex.id(), vertex.edge(i)->vertex_id()) = Max_Level;
					vertex.edge(i)->set_data(edata);
				}
			}	
		}else{
			//Vertexinfo vdata = vertex.get_data();
			vid_t curmin = vdata.level;
			vid_t old_value =curmin ; //vertex.get_data();
	
			for(int i=0; i < vertex.num_edges(); i++) {
				bidirectional_label edata = vertex.edge(i)->get_data();
				vid_t nblabel = edata.neighbor_label(vertex.id(), vertex.edge(i)->vertex_id()); 
				curmin = std::min(nblabel+1, curmin); 
			}
			if(old_value > curmin ){
				vdata.level = curmin;	
				vertex.set_data(vdata);
				converged = false;
				interval_converged = false;
		
				for(int i=0; i<vertex.num_edges(); i++){
					bidirectional_label edata = vertex.edge(i)->get_data();
					edata.my_label(vertex.id(), vertex.edge(i)->vertex_id()) = curmin;
					vertex.edge(i)->set_data(edata);
					//vertex.edge(i)->set_data(curmin);		
					if(scheduler)
						gcontext.scheduler->add_task(vertex.edge(i)->vertex_id());
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
		}
	}
	bool repeat_updates(graphchi_context &gcontext){
		if(gcontext.iteration == 0)
			return false;
		else{
			if(0 == repeat_count % 100){
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

std::map<int, int> cc_to_index;
std::vector<int> prefix_sum;
std::vector<std::vector<int> > bfs_level(top+1);

mutex lock;
vid_t getNewId(vid_t cid, vid_t level){
	std::map<int, int>::iterator iters;
	iters = cc_to_index.find(cid);
	assert(iters != cc_to_index.end());
	int row_idx = iters->second; 
	// vertices beyond top are addressed casually
	level = row_idx == top ? 0 : level;
	assert(level < bfs_level[row_idx].size());

	vid_t new_id = 0;
	lock.lock();
	new_id = bfs_level[row_idx][level];	
	bfs_level[row_idx][level]++;
	lock.unlock();
	return new_id;
} 
FILE* vfout = NULL;
FILE* efout = NULL;

struct ReMapProgram: public GraphChiProgram<VertexDataType, EdgeDataType> {
	bool converged;
	bool interval_converged;
	void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {
		//if (scheduler) gcontext.scheduler->remove_tasks(vertex.id(), vertex.id());
		if (gcontext.iteration == 0){
			if(vertex.num_edges() > 0){
				Vertexinfo vdata = vertex.get_data();  
				//source vertex in each CC
				//vdata.level = 0;
				//vertex.set_data(vdata);	
				vid_t new_id = getNewId(vdata.ccid, vdata.level);
				for(int i=0; i<vertex.num_edges(); i++){
					bidirectional_label edata = vertex.edge(i)->get_data();
					edata.my_label(vertex.id(), vertex.edge(i)->vertex_id()) = new_id;
					vertex.edge(i)->set_data(edata);
				}
				lock.lock();
					fprintf(vfout, "%u\t%u\n", new_id, vertex.id());
				lock.unlock();
			}
		}else{
			for(int i=0; i<vertex.num_outedges(); i++){
				bidirectional_label edata = vertex.outedge(i)->get_data();	
				vid_t my_id = edata.my_label(vertex.id(), vertex.outedge(i)->vertex_id());		
				vid_t nb_id = edata.neighbor_label(vertex.id(), vertex.outedge(i)->vertex_id());
				lock.lock();
				fprintf(efout, "%u\t%u\n", my_id, nb_id);
				lock.unlock();
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
	//top 			 	 = get_option_int("ntop", 20);
	Root	             = (vid_t)get_option_int("root", 0);
 //   diff_iter			=get_option_int("diff",0);
    /* Process input file - if not already preprocessed */
    int nshards             = (int) convert_if_notexists<EdgeDataType>(filename, get_option_string("nshards", "auto"));
    
        /* Run */
	graphchi_engine<Vertexinfo, EdgeDataType> engine(filename, nshards, scheduler, m); 
	//graphchi_engine<Vertexinfo, vid_t> engine(filename, nshards, scheduler, m); 

	//graphchi_engine<BFSVertexDataType, EdgeDataType> engine1(filename, nshards, scheduler, m); 

	//bool inmemmode =false;//= engine.num_vertices()*sizeof(VertexDataType) < (size_t)engine.get_membudget_mb()*1024L*1024L;
    /* Run analysis of the connected components  (output is written to a file) */
	// first run cc to detect all communities in the graph
	CCProgram ccprogram;
    engine.run(ccprogram, niters);
    
    //analyze_labels<vid_t>(filename, 20);
    //count_bfs_levels<vid_t, level_pair >(filename, Max_Level, cc_size);
    count_cc_size<Vertexinfo, level_pair>(filename, cc_size);
	assert(cc_size.size() > 0);   	 
	int tmp_sum = 0;
	// remap the id of each vertex based on undirBFS
	for(int i=0; i<(int)cc_size.size(); i++){
		if(i < top){
			if(0 == i){
				tmp_sum = 0;
			}		
			//prefix_sum[i] = tmp_sum;	
			prefix_sum.push_back(tmp_sum);
			tmp_sum += cc_size[i].count;
			cc_to_index[cc_size[i].label] = i;	
		}else{
			// store only top20 CC(account for 98% vertices)
			if(i == top){
				prefix_sum.push_back(tmp_sum);	
			}
			cc_to_index[cc_size[i].label] = top;
		}
	}
	
    m.stop_time("label-analysis");
	std::cout<<"partition bfs result size="<<cc_size.size()<<std::endl;	
	std::cout<<"before bfs phase"<<std::endl;	
	graphchi_engine<Vertexinfo, EdgeDataType> engine1(filename, nshards, scheduler, m); 
	BFSProgram bfsprogram;
	engine1.run(bfsprogram, niters);
	/*	
	graphchi_engine<Vertexinfo, EdgeDataType> engine2(filename, nshards, scheduler, m); 
	CheckProgram checkprogram;				
	engine2.run(checkprogram, niters);
	std::cout<<"partition BFS has been checked"<<std::endl;
	*/	
	
	//use 2D array to store each CC(row) and corresponding each bfs level size(colunm)
	//std::vector<std::vector<int> > bfs_level(top+1);
    //int isolated = count_bfs_size<Vertexinfo>(filename, (vid_t)top, cc_to_index, bfs_level);
    count_undirbfs_size<Vertexinfo>(filename, (vid_t)top, cc_to_index, bfs_level);
	std::cout<<"before remap phase"<<std::endl;	
	//int sum = 0;
	for(int i=0; i<(int)bfs_level.size(); i++){
		int start = prefix_sum[i];
		int tmp = 0;
		for(int j=0; j<(int)bfs_level[i].size(); j++){
			//total += bfs_level[i][j];
			tmp = bfs_level[i][j];
			bfs_level[i][j] = start;
			start += tmp;
		}
		if((i+1) < (int)bfs_level.size()){
			assert(start == prefix_sum[i+1]);	
		}
		//sum += total;
		//std::cout<<"row="<<i<<"\tcc size="<<total<<" prefix_sum="<<prefix_sum[i]<<std::endl;
	}
	//std::cout<<"all prefix sum has checked"<<std::endl;
	//std::cout<<"\nisolated vertices= "<<isolated<<"\tnonisolated= "<<sum<<"\t sum="<<isolated+sum<<std::endl;
	//uv, ue denotes vertices and edges undirected BFS
	vfout = fopen((filename+".uv").c_str(), "w+");	
	efout = fopen((filename+".ue").c_str(), "w+");
	assert(vfout != NULL && efout != NULL);

	fprintf(vfout, "# new_vid  old_vid\n");
	fprintf(efout, "# new_src  new_dst\n");
	fflush(vfout);	
	fflush(efout);
	std::cout<<"remap is started"<<std::endl;
	graphchi_engine<Vertexinfo, EdgeDataType> engine2(filename, nshards, scheduler, m); 
	ReMapProgram remapprogram;				
	engine2.run(remapprogram, 2);
	fclose(vfout);
	fclose(efout);
	std::cout<<"undirected BFS based partitioning has been remap"<<std::endl;
	
	m.start_time("label-analysis");
    /* Report execution metrics */
    metrics_report(m);
    return 0;
}

