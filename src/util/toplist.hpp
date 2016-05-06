
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
 * Tools for listing the TOP K values from a verte data file.
 */

#ifndef DEF_GRAPHCHI_TOPLIST
#define DEF_GRAPHCHI_TOPLIST

#include <vector>
#include <algorithm>
#include <errno.h>
#include <assert.h>

#include "io/stripedio.hpp"
#include "logger/logger.hpp"
#include "util/merge.hpp"
#include "util/ioutil.hpp"
#include "util/qsort.hpp"
#include "api/chifilenames.hpp"
#include "engine/auxdata/vertex_data.hpp"

namespace graphchi {
  
    template <typename VertexDataType>
    struct vertex_value {
        vid_t vertex;
        VertexDataType value;
        vertex_value() {}
        vertex_value(vid_t v, VertexDataType x) : vertex(v), value(x) {}
    };
    
    template <typename VertexDataType>
    bool vertex_value_greater(const vertex_value<VertexDataType> &a, const vertex_value<VertexDataType> &b) {
        return a.value > b.value;
    }
     
    /**
      * Reads the vertex data file and returns top N values.
      * Vertex value type must be given as a template parameter.
      * This method has been implemented in a manner to consume very little
      * memory, i.e the whole file is not loaded into memory (unless ntop = nvertices).
      * @param basefilename name of the graph
      * @param ntop number of top values to return (if ntop is smaller than the total number of vertices, returns all in sorted order)
      * @param from first vertex to include (default, 0)
      * @param to last vertex to include (default, all)
      * @return a vector of top ntop values  
     */
    template <typename VertexDataType>
    std::vector<vertex_value<VertexDataType> > get_top_vertices(std::string basefilename, int ntop, vid_t from=0, vid_t to=0) {
        typedef vertex_value<VertexDataType> vv_t;
        
        /* Initialize striped IO manager */
        metrics m("toplist");
        stripedio * iomgr = new stripedio(m);
        
        /* Initialize the vertex-data reader */
        vid_t readwindow = 1024 * 1024;
        size_t numvertices = get_num_vertices(basefilename);
        vertex_data_store<VertexDataType> * vertexdata =
        new vertex_data_store<VertexDataType>(basefilename, numvertices, iomgr);
           
        if ((size_t)ntop > numvertices) {
            ntop = (int)numvertices;
        }
                
        /* Initialize buffer */
        vv_t * buffer_idxs = (vv_t*) calloc(readwindow, sizeof(vv_t));
        vv_t * topbuf = (vv_t*) calloc(ntop, sizeof(vv_t));
        vv_t * mergearr = (vv_t*) calloc(ntop * 2, sizeof(vv_t));
        
        /* Iterate the vertex values and maintain the top-list */
        size_t idx = 0;
        vid_t st = 0;
        vid_t en = numvertices - 1;

        int count = 0;
        while(st <= numvertices - 1) {
            en = st + readwindow - 1;
            if (en >= numvertices - 1) en = numvertices - 1;
            
            /* Load the vertex values */
            vertexdata->load(st, en);
            
            int nt = en - st + 1;
            int k = 0;
            VertexDataType minima = VertexDataType();
            if (count > 0) {
                minima = topbuf[ntop - 1].value; // Minimum value that should be even considered
            }
            for(int j=0; j < nt; j++) {
                VertexDataType& val = *vertexdata->vertex_data_ptr(j + st);
                if (count == 0 || (val > minima)) {
                    buffer_idxs[k] = vv_t((vid_t)idx + from, val);
                    k++;
                }
                idx++;
            }
            nt = k; /* How many were actually included */
            
            /* Sort buffer-idxs */
            quickSort(buffer_idxs, nt, vertex_value_greater<VertexDataType>);
            
            /* Merge the top with the current top */
            if (count == 0) {
                /* Nothing to merge, just copy */
                memcpy(topbuf, buffer_idxs, ntop * sizeof(vv_t));
            } else {
                // void merge(ET* S1, int l1, ET* S2, int l2, ET* R, F f) {
                merge<vv_t>(topbuf, ntop, buffer_idxs, std::min(ntop, nt), mergearr,  vertex_value_greater<VertexDataType>);
                memcpy(topbuf, mergearr, ntop * sizeof(vv_t));
            }
            
            count++;
            st += readwindow;
        }
                   
        /* Return */
        std::vector< vv_t > ret;
        for(int i=0; i < ntop; i++) {
            ret.push_back(topbuf[i]);
        }
        free(buffer_idxs);
        free(mergearr);
        free(topbuf);

        delete vertexdata;
        delete iomgr;

        return ret;
    }


//////////////////////////////////////////////////////////////////////////////////  
    template <typename VertexDataType>
    struct deg_vertex {
        int vertex;//vid
        int value;//vertex degree
        //deg_vertex() {}
        deg_vertex(vid_t v, vid_t x) : vertex(v), value(x) {}
    };
    
    template <typename VertexDataType>
    bool vertex_deg_greater(const deg_vertex<VertexDataType> &a, const deg_vertex<VertexDataType> &b) {
        return a.value > b.value;
    }
	
	template <typename VertexDataType>
    bool vertex_deg_less(const deg_vertex<VertexDataType> &a, const deg_vertex<VertexDataType> &b) {
        return a.value < b.value;
    }

	// get top K unvisited degree vertices
    template <typename VertexDataType>
    int get_top_deg_vertices(std::string basefilename, int* Rt, int ntop, vid_t from=0, vid_t to=0) {
        typedef deg_vertex<VertexDataType> vv_t;
       	assert(Rt != NULL); 
        /* Initialize striped IO manager */
        metrics m("Root init");
        stripedio * iomgr = new stripedio(m);
        /* Initialize the vertex-data reader */
        vid_t readwindow = 1024 * 1024;
        size_t numvertices = get_num_vertices(basefilename);
        vertex_data_store<VertexDataType> * vertexdata =
        new vertex_data_store<VertexDataType>(basefilename, numvertices, iomgr);
           
        if ((size_t)ntop > numvertices) {
            ntop = (int)numvertices;
        }
                
        /* Initialize buffer */
        vv_t * buffer_idxs = (vv_t*) calloc(readwindow, sizeof(vv_t));
        vv_t * topbuf = (vv_t*) calloc(ntop, sizeof(vv_t));
        vv_t * mergearr = (vv_t*) calloc(ntop * 2, sizeof(vv_t));
        
        /* Iterate the vertex values and maintain the top-list */
        size_t idx = 0;
        vid_t st = 0;
        vid_t en = numvertices - 1;
	    //int max_deg = 0;
		//int max_id =0; 

        int count = 0;
        while(st <= numvertices - 1) {
            en = st + readwindow - 1;
            if (en >= numvertices - 1) en = numvertices - 1;
            
            /* Load the vertex values */
            vertexdata->load(st, en);
            
            int nt = en - st + 1;
            int k = 0;
            //VertexDataType minima = VertexDataType();
			int mini_deg = 0;
            if (count > 0) {
                mini_deg = topbuf[ntop - 1].value; // Minimum value that should be even considered
            }
            for(int j=0; j < nt; j++) {
                VertexDataType& val = *vertexdata->vertex_data_ptr(j + st);
				if(!val.is_confirm())
				{// only  unvisited vertices can be roots
				/*
				if(max_deg < val.get_degree()){
					max_deg = val.get_degree();
					max_id = idx;
				}
				*/
                if (count == 0 || (mini_deg < val.get_degree())) {
                    buffer_idxs[k] = vv_t((int)idx + from, val.get_degree());
                    k++;
                }
				}
                idx++;
            }
            nt = k; /* How many were actually included */
            
            /* Sort buffer-idxs */
            quickSort(buffer_idxs, nt, vertex_deg_greater<VertexDataType>);
            
            /* Merge the top with the current top */
            if (count == 0) {
                /* Nothing to merge, just copy */
                memcpy(topbuf, buffer_idxs, ntop * sizeof(vv_t));
            } else {
                // void merge(ET* S1, int l1, ET* S2, int l2, ET* R, F f) {
                merge<vv_t>(topbuf, ntop, buffer_idxs, std::min(ntop, nt), mergearr, vertex_deg_greater<VertexDataType>);
                memcpy(topbuf, mergearr, ntop * sizeof(vv_t));
            }
            
            count++;
            st += readwindow;
        }
                   
        /* Return */
        /*
		std::vector< vv_t > ret;
        for(int i=0; i < ntop; i++) {
            ret.push_back(topbuf[i]);
        }
		*/
		static int icount =1;
		std::string degname = basefilename+".rdeg";
		FILE* degfile = NULL;
		
		degfile = fopen(degname.c_str(), "a+");
		std::cout<<"degree of Roots::"<<std::endl;
		fprintf(degfile, "iteration-%d\t", icount++);
		for(int i=0; i<ntop; i++)
		{
			Rt[i] = topbuf[i].vertex;
			std::cout<<"vid-"<<Rt[i]<<"deg:"<<topbuf[i].value<<'\t';
			fprintf(degfile,"vid-%d deg-%d\t", Rt[i], topbuf[i].value);
		}
		/*
		std::cout<<"assert max_deg="<<max_deg<<"\t topbuf[0].value="<<topbuf[0].value<<std::endl;
		assert(topbuf[0].value == max_deg);
		*/
		fprintf(degfile, "\n");
		fclose(degfile);
		std::cout<<std::endl;
        free(buffer_idxs);
        free(mergearr);
        free(topbuf);

        delete vertexdata;
        delete iomgr;
		
        return ntop;
    }
//};

// get top degree vertex
    template <typename VertexDataType>
    int get_top_deg_vertex(std::string basefilename, bool max = true, vid_t from=0, vid_t to=0) {
        //typedef deg_vertex<VertexDataType> vv_t;
       	//assert(Rt != NULL); 
        /* Initialize striped IO manager */
        metrics m("top degree vertex");
        stripedio * iomgr = new stripedio(m);
        /* Initialize the vertex-data reader */
        vid_t readwindow = 1024 * 1024;
        size_t numvertices = get_num_vertices(basefilename);
        vertex_data_store<VertexDataType> * vertexdata =
        new vertex_data_store<VertexDataType>(basefilename, numvertices, iomgr);
        /* 
        if ((size_t)ntop > numvertices) {
            ntop = (int)numvertices;
        }
        */     
        /* Initialize buffer */
		/*
        vv_t * buffer_idxs = (vv_t*) calloc(readwindow, sizeof(vv_t));
        vv_t * topbuf = (vv_t*) calloc(ntop, sizeof(vv_t));
        vv_t * mergearr = (vv_t*) calloc(ntop * 2, sizeof(vv_t));
       	*/ 
        /* Iterate the vertex values and maintain the top-list */
		
        //size_t idx = 0;
        vid_t st = 0;
        vid_t en = numvertices - 1;
	    //int max_deg = 0;
		//int max_id =0; 

        int count = 0;
		int extreme_id = 0;
		// get max or min degree vertex
		int extreme_deg = max ? 0 : 2100000000;
        while(st <= numvertices - 1) {
            en = st + readwindow - 1;
            if (en >= numvertices - 1) en = numvertices - 1;
            
            /* Load the vertex values */
            vertexdata->load(st, en);
            
            int nt = en - st + 1;
            int k = 0;
            //VertexDataType minima = VertexDataType();
			//int mini_deg = 0;
			/*
            if (count > 0) {
                mini_deg = topbuf[ntop - 1].value; // Minimum value that should be even considered
            }
			*/
			for(int j=0; j < nt; j++) {
				VertexDataType& val = *vertexdata->vertex_data_ptr(j + st);
				if(max){
					if(val.get_degree() > extreme_deg){
						extreme_id = j+st; 
						extreme_deg = val.get_degree();
					}
				}else{
					if(val.get_degree() > 0){
						if(val.get_degree() < extreme_deg){
							extreme_id = j+st;
							extreme_deg = val.get_degree();
						}
					}
				}
			}
            nt = k; /* How many were actually included */
            
            /* Sort buffer-idxs */
            //quickSort(buffer_idxs, nt, vertex_deg_greater<VertexDataType>);
           	/* 
            if (count == 0) {
                memcpy(topbuf, buffer_idxs, ntop * sizeof(vv_t));
            } else {
                merge<vv_t>(topbuf, ntop, buffer_idxs, std::min(ntop, nt), mergearr, vertex_deg_greater<VertexDataType>);
                memcpy(topbuf, mergearr, ntop * sizeof(vv_t));
            }
           	*/ 
            count++;
            st += readwindow;
        }
		std::cout<<"extreme_id="<<extreme_id<<" extreme_degree="<<extreme_deg<<std::endl;
        return extreme_id;
    }
};

#endif

