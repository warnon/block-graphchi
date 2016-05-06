
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
 * Analyses output of label propagation algorithms such as connected components
 * and community detection. Memory efficient implementation.
 *
 * @author Aapo Kyrola
 */


#include <vector>
#include <algorithm>
#include <errno.h>
#include <assert.h>
#include <fstream>
#include <iostream>
//#include<fiostream.h>

#include "io/stripedio.hpp"
#include "logger/logger.hpp"
#include "util/merge.hpp"
#include "util/ioutil.hpp"
#include "util/qsort.hpp"
#include "api/chifilenames.hpp"
#include "engine/auxdata/vertex_data.hpp"

#ifndef DEF_GRAPHCHI_BFSLABELANALYSIS
#define DEF_GRAPHCHI_BFSLABELANALYSIS

using namespace graphchi;
/*
template <typename LabelType>
struct labelcount_bfs {
    //int  label;
    LabelType  label;
    int count;  // Count excludes the vertex which has its own id as the label. (Important optimization)
    labelcount_bfs(LabelType la, int c) : label(la), count(c) {}
    labelcount_bfs(LabelType la) : label(la), count(0) {}
    //labelcount_tt(int la) : label(la), count(0) {}
    //labelcount_tt() {}
};

template <typename LabelType>
bool label_ct_greater(const labelcount_bfs<LabelType> &a, const labelcount_bfs<LabelType> &b) {
    return a.count > b.count;
}

template <typename LabelType>
bool label_greater(const labelcount_bfs<LabelType> &a, const labelcount_bfs<LabelType> &b) {
    return a.label.level > b.label.level;
}

template <typename LabelType>
bool label_less(const labelcount_bfs<LabelType> &a, const labelcount_bfs<LabelType> &b) {
    return a.label.level < b.label.level;
}



template <typename LabelType>
bool operator !=(const labelcount_bfs<LabelType> &a, const labelcount_bfs<LabelType> &b) {
    return a.label.level != b.label.level;
}
*/

//count the size of each bfs level but exclude isolated vertices (whose level is isolated_level)
//template <typename LabelType, typename Pair>
template <typename LabelType>
int count_undirbfs_size(std::string basefilename, int ntop, std::map<int, int>& ccid_to_idx, std::vector< std::vector<int> >& bfs_vector){    
    //typedef labelcount_bfs<LabelType> labelcount_t;
    /**
     * NOTE: this implementation is quite a mouthful. Cleaner implementation
     * could be done by using a map implementation. But STL map takes too much
     * memory, and I want to avoid Boost dependency - which would have boost::unordered_map.
     */
    metrics m("label_count");
    stripedio * iomgr = new stripedio(m);
    /* Initialize the vertex-data reader */
    vid_t readwindow = 1024 * 1024;
    vid_t numvertices = (vid_t) get_num_vertices(basefilename);
    vertex_data_store<LabelType> * vertexdata =
    new vertex_data_store<LabelType>(basefilename, numvertices, iomgr);
    
    //std::vector<labelcount_t> curlabels;
    //bool first = true;
	int isolated_vertices = 0;
	//int	visited_vertices = 0; 
    //labelcount_t * buffer = NULL;
    //LabelType * buffer = NULL;
    //buffer = (LabelType*) calloc(readwindow, sizeof(LabelType));
   	//assert(buffer != NULL); 
    /* Iterate the vertex values and maintain the top-list */
    vid_t st = 0;
    vid_t en = numvertices - 1;
	vid_t curvid = 0;
    while(st <= numvertices - 1){
        en = st + readwindow - 1;
        if (en >= numvertices - 1) en = numvertices - 1;
        
        /* Load the vertex values */
        vertexdata->load(st, en);
        
        int nt = en - st + 1;
        //int ki=0;
        /* Mark vertices with its own label with 0xffffffff so they will be ignored */
		for(int i=0; i < nt; i++) { 
			LabelType lab = *vertexdata->vertex_data_ptr(i + st);
			if(lab.level == 0xffffffff){
				isolated_vertices++;
			}else{	
				std::map<int, int>::iterator iter;
				iter = ccid_to_idx.find(lab.ccid); 
				assert(iter != ccid_to_idx.end());	
				// get the index in the 2-D array 
				int i_idx = iter->second;
				//vertices beyond top20 is marked casually
				if(i_idx == ntop){
					if(bfs_vector[i_idx].size() == 0){
						bfs_vector[i_idx].resize(1,0);
					}
					bfs_vector[i_idx][0]++;
				}else{
					if(lab.level >= bfs_vector[i_idx].size()){
						bfs_vector[i_idx].resize(lab.level+1, 0);		
					}
					bfs_vector[i_idx][lab.level]++;	
				}
			}
			curvid++;
		}
		/*
        nt = ki;
        quickSort(buffer, nt, label_less<LabelType>);// important modify!!        
        std::vector<labelcount_t> newlabels;
        newlabels.reserve(nt);
        LabelType lastlabel = LabelType(0xffffffff, 0);
		for(int i=0; i < nt; i++) {
			if(buffer[i].ccid != LabelType(0xffffffff, 0).ccid){
				if (buffer[i].ccid !=  lastlabel.ccid) {
					newlabels.push_back(labelcount_t(buffer[i], 1));
				} else {
					newlabels[newlabels.size() - 1].count++;
				}
				lastlabel = buffer[i];
			}
		}
        if (first) {
            for(int i=0; i < (int)newlabels.size(); i++) {
                curlabels.push_back(newlabels[i]);
            }
        } else {
            int cl = 0;
            int nl = 0;
            std::vector< labelcount_t > merged;
            merged.reserve(curlabels.size() + newlabels.size());
            while(cl < (int)curlabels.size() && nl < (int)newlabels.size()) {
                if (newlabels[nl].label.ccid == curlabels[cl].label.ccid) {
                    merged.push_back(labelcount_t(newlabels[nl].label, newlabels[nl].count + curlabels[cl].count));
                    nl++; cl++;
                } else {
                    if (newlabels[nl].label.ccid < curlabels[cl].label.ccid) {
                        merged.push_back(newlabels[nl]);
                        nl++;
                    } else {
                        merged.push_back(curlabels[cl]);
                        cl++;
                    }
                }
            }
            while(cl < (int)curlabels.size()) merged.push_back(curlabels[cl++]);
            while(nl < (int)newlabels.size()) merged.push_back(newlabels[nl++]);
            
            curlabels = merged;
        }
        
        first = false;
		*/
        st += readwindow;
    }
	/*	
    std::sort(curlabels.begin(), curlabels.end(), label_ct_greater<LabelType>);
		
    std::string outname = basefilename + ".cc_size";
    std::ofstream resf;
	int sum_vertices = 0;
    resf.open(outname.c_str());
    if (resf == NULL) {
        logstream(LOG_ERROR) << "Could not write label outputfile : " << outname << std::endl;
        return curlabels.size();
    }
	//cc_result.reserve(curlabels.size());
    for(int i=0; i < (int) curlabels.size(); i++) {
        resf << curlabels[i].label << "," << curlabels[i].count << std::endl;
		//cc_result.push_back(Pair(curlabels[i].label.ccid, curlabels[i].count));
		sum_vertices += curlabels[i].count;
		
    }
    resf.close();
	
	std::cout<<"valid vertices="<<sum_vertices<<"\tvisited vertices="<<visited_vertices<<std::endl;
   	sum_vertices += isolated_vertices; 
    std::cout << "Total number of different labels (components/communities): " << curlabels.size()
			 <<"\tisolated_vertices="<<isolated_vertices<<"\tsum_vertices="<<sum_vertices<< std::endl;
    std::cout << "List of labels was written to file: " << outname << std::endl;
	
	for(int i=0; i < (int)std::min((size_t)20, curlabels.size()); i++) {
        std::cout << (i+1) << ". label: " << curlabels[i].label << ", size: " << curlabels[i].count + 1 << std::endl;
    }
	*/
	//int size = curlabels.size();
	//std::cout<<"free buffer"<<std::endl;
    //free(buffer);
   	//buffer = NULL; 
	//std::cout<<"after free buffer"<<std::endl;
    delete vertexdata;
    delete iomgr;
	
	return isolated_vertices;	
	//return curlabels.size();	
}


#endif


