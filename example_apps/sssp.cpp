#include <stdlib.h>
#include <cmath>
#include <string>

#include "graphchi_basic_includes.hpp"
#include "util/labelanalysis.hpp"
#include "util/toplist.hpp"
//#include "util/vertexlist.hpp"

using namespace graphchi;

bool        scheduler = true;
vid_t		single_source = 0;
bool		reset_edge_value = false;
bool		num_tasks_print = false;

/**
 * Type definitions. Remember to create suitable graph shards using the
 * Sharder-program. 
 */
struct edgeWithSrcValue{
	float value;
	float srcValue;

	edgeWithSrcValue(){
	}

	edgeWithSrcValue(float v){
		value = v;
	}
	
	edgeWithSrcValue(float v,float sv){
		value = v;
		srcValue = sv;
	}

/*
	void updateSrcValue(float sv){
		srcValue = sv;
	}*/
};

typedef float VertexDataType;       // vid_t is the vertex id type
typedef edgeWithSrcValue EdgeDataType;

static void parse(EdgeDataType& edata, const char* s){
	edata.value = atof(s);
}


//typedef float InputEdgeDataType;

const float infinity=99999999.0;
static uint32_t repeat_count = 0;

bool converged ;
bool interval_converged;

/**
 * GraphChi programs need to subclass GraphChiProgram<vertex-type, edge-type> 
 * class. The main logic is usually in the update function.
 */
struct SSSPProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {
    

    /**
     *  Vertex update function.
     *  On first iteration ,each vertex chooses a label = the vertex id.
     *  On subsequent iterations, each vertex chooses the minimum of the neighbor's
     *  label (and itself). 
     */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {
        
        if (scheduler) gcontext.scheduler->remove_tasks_now(vertex.id(), vertex.id());
        

        if (gcontext.iteration == 0) {
			if( vertex.id() == single_source ){
				vertex.set_data( 0.0 );
			
				converged = false;
				for(int j=0;j<vertex.num_outedges();j++){
            		if (scheduler)
						gcontext.scheduler->add_task(vertex.outedge(j)->vertexid, false);
					
					EdgeDataType edata = vertex.outedge(j)->get_data();
					edata.srcValue = 0.0;
					vertex.outedge(j)->set_data(edata);
				}
			}else{
				vertex.set_data(infinity);
				for(int j=0;j<vertex.num_outedges();j++){
					if (scheduler)
						gcontext.scheduler->add_task(vertex.outedge(j)->vertexid, false);
					EdgeDataType edata = vertex.outedge(j)->get_data();
					edata.srcValue = infinity;
					vertex.outedge(j)->set_data(edata);
					//vertex.outedge(j)->set_data(edgeWithSrcValue(vertex.outedge(j)->get_data().value,vertex.get_data()));
				}
			}

			//set value of out-edges
			float edge_value = 0.0;
			if( reset_edge_value ){
				srand( time(0) );
				for( int i=0; i<vertex.num_outedges(); i++){
					graphchi_edge<edgeWithSrcValue> *edge = vertex.outedge(i);
					edge_value = (float)rand()/(float)RAND_MAX;
					edge->set_data( edgeWithSrcValue(edge_value,0) );
					if(vertex.id()==0) printf( "value of out edge %d is %f\n", i, edge_value );
				}
			}

        }else{
			//std::cout<<vertex.id()<<" into >>>>>>>>>>>>>"<<std::endl;
			float min = vertex.get_data();
			for( int i=0; i<vertex.num_inedges(); i++ ){
				//float tmp = ( neighbor_value(vertex.inedge(i)) + vertex.inedge(i)->get_data().value ); 
				EdgeDataType edata = vertex.inedge(i)->get_data();
				float tmp = edata.srcValue + edata.value; 
				if( tmp < min ){
					min = tmp;
				}
			}
			if(min < vertex.get_data()){
				//set_data( vertex, min);
				vertex.set_data(min);
				for(int j=0; j<vertex.num_outedges(); j++){
					vid_t out_edge_vid = vertex.outedge(j)->vertexid;
					if (scheduler){
						if(out_edge_vid >= gcontext.interval_st){	
							//vertex is updated in this iteration
							gcontext.scheduler->add_task(vertex.outedge(j)->vertexid, true);
						}else{
							//vertex in interval that is before this interval is scheduled next iteration
							gcontext.scheduler->add_task(vertex.outedge(j)->vertexid, false);
						}
					}
					EdgeDataType edata = vertex.outedge(j)->get_data();
					edata.srcValue = min;
					vertex.outedge(j)->set_data(edata);
				}  
				converged = false;
				interval_converged = false;
			}
		}
	}

	bool repeat_updates(graphchi_context& gcontext){
		if(gcontext.iteration == 0)
			return false;
		else{
			if(0 == repeat_count % 100 || interval_converged){
			logstream(LOG_INFO)<<"repeat_upates="<<interval_converged<<"\t iteration="
				<<gcontext.iteration<<"\t repeat_count="<<repeat_count<<std::endl;
			}
			repeat_count = interval_converged ? 0 : repeat_count+1;	
			if(interval_converged){
				assert(gcontext.scheduler->total_tasks(gcontext.interval_st, gcontext.interval_en) == 0);
			}
			return !interval_converged;
		}
	}
    /**
     * Called before an iteration starts.
     */
    void before_iteration(int iteration, graphchi_context &info) {
		
        converged = iteration > 0;
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
		if(num_tasks_print){
			vid_t sum = ginfo.scheduler->total_tasks(window_st, window_en);
			logstream(LOG_INFO)<<"num of vertices scheduled="<<sum<<"/"<<window_en-window_st+1<<std::endl;		
			//sum_tasks += sum;
		}	
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
    metrics m("blockupdate-SSSP");
    
    /* Basic arguments for application */
    std::string filename = get_option_string("file");  // Base filename
    int niters           = get_option_int("niters", 100000); // Number of iterations (max)
    scheduler            = get_option_int("scheduler", false);
	single_source		 = get_option_int("root", 0);
	reset_edge_value	 = get_option_int("reset_edge_value", false);
	num_tasks_print		 = get_option_int("print", false);
	//int ntop 			 = get_option_int("top",30);
//    
//	std::cout << "------------------------\tExecution Mode\t----------------" << std::endl;
//	if ( allow_race )
//		std::cout << "--The Execution allows RACE! (non-deterministic)" << std::endl;
//	else
//		std::cout << "--The Execution is Deterministic!" << std::endl;
//	if ( scheduler )
//		std::cout << "--The Execution uses scheduler. " << std::endl;
//	else
//		std::cout << "--The Execution uses no scheduler. " << std::endl;
//
	std::cout << "--Single source = " << single_source << std::endl;
	if( reset_edge_value )
		std::cout << "--Reset edge value = " << "TRUE" << std::endl;
	else
		std::cout << "--Reset edge value = " << "FALSE" << std::endl;

	std::cout << "------------------------------------------------------------" << std::endl;

    /* Process input file - if not already preprocessed */
    //int nshards             = (int) convert_if_notexists<InputEdgeDataType,EdgeDataType>(filename, get_option_string("nshards", "auto"));
    int nshards             = (int) convert_if_notexists<EdgeDataType>(filename, get_option_string("nshards", "auto"));
    
    /* Run */
    SSSPProgram program;
    graphchi_engine<VertexDataType, EdgeDataType> engine(filename, nshards, scheduler, m); 
	
    engine.run(program, niters);
	/*
	if( allow_race )
		engine.set_enable_deterministic_parallelism(false);
	else
		engine.set_enable_deterministic_parallelism(true);
    engine.run(program, niters);
   	*/ 
    /* Run analysis of the connected components  (output is written to a file) */
//    m.start_time("label-analysis");
//    
//    analyze_labels<vid_t>(filename);
//    
//    m.stop_time("label-analysis");
    
    /* Report execution metrics */
		
	//std::vector<vertex_value<float> > top = get_header_vertices<float>(filename,ntop);
	/*
	std::vector<vertex_value<float> > top = get_top_vertices<float>(filename,ntop);
	std::cout<<"SSSP Print head "<<ntop<<" vertices: "<< std::endl;
	for(int i=0;i<(int)top.size();i++)
		std::cout<<top[i].vertex<<"\t "<<top[i].value<<std::endl;
	*/	
	analyze_labels<VertexDataType>(filename, 30);
	//get_top_vertices<VertexDataType>(filename, 20);
    metrics_report(m);
    return 0;
}

