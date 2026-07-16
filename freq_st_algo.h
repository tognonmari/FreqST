#pragma once

#include <cmath>
#include <iostream>
#include <limits>
#include <list>
#include <vector>
#include <random>
#include <algorithm>
#include <random>
#include <boost/container/flat_set.hpp>
#include "roaring.hh"

#include "free_space_graph_free_axis.h"
#include "kdtree_range_search.h"
#include "metric_space.h"
#include "io.h"
#include "trajectory.h"
#include "canonical_pathlets.h"
#include "curve_simplification.h"
#include "grid_range_search.h"
#include "utility.h"
namespace frechet{

struct freq_subtrajectory_algo_output_config{

    bool maximal = false;
    bool keep_matching_ids = false;
    int min_length = 1;

};
template<metric_space m_space>
class freq_subtrajectory_sampler{

    public:
    using space = m_space;
    using point_t = space::point_t;
    using distance_function_t = space::distance_function_t;
    using distance_t = distance_function_t::distance_t;
    using trajectory_t = trajectory_collection<space>;
    using index_t = trajectory_t::index_t;
    using range_search_t = grid_range_search<space>;
    using subtrajectory_t = trajectory_t::subtrajectory_t;
    using id_t = trajectory_t::id_t;
    
    private: 
        struct transaction_and_points_pair{
            //Fields
            index_t traj_id;
            std::vector<std::pair<index_t, std::vector<index_t>>> closeby_points; //this structure represents the closeby beginnings, i need the extra vectors to count the number of pathlets in the nested case. 
            int capacity_ub;
            int num_closeby_pathlets() const {
                
                return std::accumulate(closeby_points.begin(),closeby_points.end(),size_t{0},[&](size_t acc, const auto& kv){return acc + kv.second.size();});
            }
            //I want to rank the more popular transactions first 
            bool operator<(const transaction_and_points_pair& other) const {
                int my_size = this->num_closeby_pathlets();
                int other_size = other.num_closeby_pathlets();
                if(my_size<other_size){
                    return false;
                }
                if(my_size>other_size){
                    return true;
                }
                return traj_id<other.traj_id;
            }

            bool operator==(const transaction_and_points_pair& other ) const{
                return traj_id == other->traj_id;
            }
            
        };
    public:
        freq_subtrajectory_sampler(const trajectory_t& trajectory, const trajectory_t& pathlets,
            float eps, 
            float del, 
            distance_t radius, int minimum_length, int random_seed, double grid_side_factor, bool thorough, bool inspect) : the_trajectory(trajectory), the_pathlets(pathlets), epsilon(eps), delta(del), distance_threshold(radius), 
                                                                                                                min_length(minimum_length), seed(random_seed), grid_side_factor(grid_side_factor), thorough(thorough) , inspect(inspect){
                the_pathlet_trees.reserve(the_pathlets.get_id_at(the_pathlets.get_actual_size()-1)+1);
                assert(the_pathlet_trees.capacity()>0);
                std::mt19937 seeded_generator(seed);
                this->mt = seeded_generator;
            }
        
        void fill_beginnings_of_pathlet_vector(){

            //Given the dataset (this->the_trajectory;) and the min_legth, 
            //we detect the beginning points of the long pathelets and insert them in the beggining vectors.
            //These beginnings will be used to instantiate the grid in the vc dimension
            std::cout <<std::format("I am about to fill the beginnings vector")<< std::endl;
            id_t last_trajectory = the_pathlets.get_id_at(the_pathlets.get_actual_size()-1);
            id_t current_visiting_id = the_pathlets.get_id_at(0);
            //Step one iterate through whole dataset and build a pathlet tree
            bool last_iter = false;
            while(true){
                
                if(current_visiting_id ==last_trajectory){
                    last_iter = true;
                }
                trajectory_t pathlet_mother = the_pathlets.slice_trajectory_by_id(current_visiting_id);
                LightWeightBinaryPathletTree pathlet_tree(pathlet_mother, current_visiting_id, floor(log2(pathlet_mother.total_size())) + 1,1);
                
                //std::cout<< "Pathelet tree for tid "<< current_visiting_id<< std::endl;
                //std::cout << pathlet_tree.toString()<< std::endl;
                //std::cout << std::format("pathlet tree vector has capacity {}.\n", the_pathlet_trees.capacity());
                if(current_visiting_id >= the_pathlet_trees.size()){

                    the_pathlet_trees.resize(current_visiting_id+1);

                }
                this->the_pathlet_trees[current_visiting_id] = (pathlet_tree);
                //For debugging purposes: print the tree
                
                std::vector<PathletNode<space>> min_length_pathlets = pathlet_tree.getMinLengthPathlets(this-> min_length);
                //std::cout<< "CALLING HERE for pm "<< current_visiting_id << std::endl;
                index_t trajectory_offset = the_pathlets.get_first_point_in_trajectory(current_visiting_id);
                for (auto& pathlet_node: min_length_pathlets){

                    index_t left_extreme = pathlet_node.pathlet.first;
                    this->pathlet_beginnings[left_extreme+ trajectory_offset].first = the_pathlets[left_extreme + trajectory_offset];
                    this->pathlet_beginnings[left_extreme+ trajectory_offset].second.push_back(trajectory_offset + pathlet_node.pathlet.second);
                    //this->pathlet_beginnings.insert({left_extreme + trajectory_offset, the_pathlets[left_extreme + trajectory_offset]});
                }

                //Move on to the next trajectory
                if(last_iter){
                    break;
                }

                current_visiting_id = the_pathlets.get_id_at(trajectory_offset + pathlet_mother.get_actual_size());

            }
            //For each pathlet tree, go down from the root and if there's a node @ that level of length >= min_length add it.
            //Stop when along the level all pathelts are shorter. 
             
            //For debugging purposes: print the pathlet beginnings 
            std::cout << "The pathlet beginnings are : "<< std::endl;
            //for (const auto& pt : this->pathlet_beginnings){
            //    std::cout<< pt.first <<"\n";
            //}



        }
        
        void fill_ends_of_pathlet_vector(){
            //Given the dataset (this->the_trajectory;) and the min_legth, 
            //we detect the beginning points of the long pathelets and insert them in the beggining vectors.
            //These beginnings will be used to instantiate the grid in the vc dimension
            id_t last_trajectory = the_pathlets.get_id_at(the_pathlets.get_actual_size()-1);
            id_t current_visiting_id = the_pathlets.get_id_at(0);
            //Step one iterate through whole dataset and build a pathlet tree
            bool last_iter = false;
            while(true){
                
                if(current_visiting_id ==last_trajectory){
                    last_iter = true;
                }
                trajectory_t pathlet_mother = the_pathlets.slice_trajectory_by_id(current_visiting_id);
                BinaryPathletTree pathlet_tree(pathlet_mother, current_visiting_id, floor(log2(pathlet_mother.total_size())) + 1,min_length);
                //For debugging purposes: print the tree
                //std::cout<< "Pathelet tree for tid "<< current_visiting_id<< std::endl;
                //std::cout << pathlet_tree.toString()<< std::endl;
                
                std::vector<PathletNode<space>> min_length_pathlets = pathlet_tree.getMinLengthPathlets(this-> min_length);
                //std::cout<< "CALLING HERE for pm "<< current_visiting_id << std::endl;
                index_t trajectory_offset = the_pathlets.get_first_point_in_trajectory(current_visiting_id);
                for (auto& pathlet_node: min_length_pathlets){

                    index_t right_extreme = pathlet_node.pathlet.second;
                    index_t left_extreme = pathlet_node.pathlet.first; //I use the first point of the pathelt as a key to the pathlet
                    this->pathlet_ends[right_extreme + trajectory_offset].first = the_pathlets[right_extreme + trajectory_offset];
                    this->pathlet_ends[right_extreme + trajectory_offset].second.push_back(left_extreme + trajectory_offset);
                    //this->pathlet_ends.insert({left_extreme + trajectory_offset, the_pathlets[right_extreme + trajectory_offset]});
                }

                //Move on to the next trajectory
                if(last_iter){
                    break;
                }

                current_visiting_id = the_pathlets.get_id_at(trajectory_offset + pathlet_mother.get_actual_size());

            }
            //For each pathlet tree, go down from the root and if there's a node @ that level of length >= min_length add it.
            //Stop when along the level all pathelts are shorter. 
             
            //For debugging purposes: print the pathlet beginnings 
            //std::cout << "The pathlet beginnings are : "<< std::endl;
            //for (const auto& pt : this->pathlet_ends){
            //    std::cout<< pt.first <<"\n";
            //}


        }
        
        void generate_chernoff_sample(){

            //Step 1: compute sample size according to Chernoff rule. 

            (this->sampled_trajs_ids).clear();

            int sample_size = (int) (3*4 / (epsilon * epsilon)) * log(2 * this->total_pathlet_number_respecting_ids() / delta);
            std::cout << "Chernoff sample size with espilon "<<epsilon << ",delta "<< delta <<" is: "<< sample_size <<std::endl;
            //Step 2: assert sampling is worthwhile
            if(sample_size > the_trajectory.num_trajectories()){

                std::cerr << "Chernoff Bound was too loose for your dataset."<< std::endl;

                std::exit(1);

            }
            
            //Step 3: Sample indexes without replacement
            
            this->sample_trajectories(sample_size); 

        }

        void generate_vc_sample(){

            //Step 1: compute sample size according to rule. 

            (this->sampled_trajs_ids).clear();

            int sample_size = (int) (2 / (epsilon * epsilon)) * (this->vc_dim()+ log(2 / delta));
            std::cout << "VCdim sample size with espilon "<<epsilon << ",  delta "<< delta <<", radius "<< distance_threshold<< " is: "<< sample_size <<std::endl;
            //Step 2: assert sampling is worthwhile
            
            if(sample_size > the_trajectory.num_trajectories()){

                std::cerr << "VC Bound was too loose for your dataset."<< std::endl;

                std::exit(1);

            }
            //Step 3: Sample indexes with replacement
            this->sample_trajectories(sample_size);
            
        }

        void generate_vc_no_erase_sample(){

            //Step 1: compute sample size according to rule. 

            (this->sampled_trajs_ids).clear();

            int sample_size = (int) (2 / (epsilon * epsilon)) * (this->vc_dim_no_erase()+ log(2 / delta));
            std::cout << "VCdim sample size with espilon "<<epsilon << ",  delta "<< delta <<", radius "<< distance_threshold<< " is: "<< sample_size <<std::endl;
            //Step 2: assert sampling is worthwhile
            
            if(sample_size > the_trajectory.num_trajectories()){

                std::cerr << "VC Bound was too loose for your dataset."<< std::endl;

                std::exit(1);

            }
            //Step 3: Sample indexes with replacement
            this->sample_trajectories(sample_size);
            
        }
        //VC dimension estimate is faster but coarse
        void generate_rough_vc_sample(){

            //Step 1: compute sample size according to rule. 
            (this->sampled_trajs_ids).clear();

            int sample_size = (int) (2 / (epsilon * epsilon)) * (this->rough_vc_dim() + log(2 / delta));
            std::cout << " Rough VCdim sample size with espilon "<<epsilon << ",  delta "<< delta <<", radius "<< distance_threshold<< " is: "<< sample_size <<std::endl;
            //Step 2: assert sampling is worthwhile
            if(sample_size > the_trajectory.num_trajectories()){

                std::cerr << "VC Bound was too loose for your dataset."<< std::endl;

                std::exit(1);

            }
            //Step 3: Sample indexes with replacement
            this->sample_trajectories(sample_size);

        }
        
        void generate_turbo_rough_vc_no_erase_sample(){

            //Step 1: compute sample size according to rule. 
            (this->sampled_trajs_ids).clear();

            int sample_size = (int) (2 / (epsilon * epsilon)) * (this->turbo_rough_vc_dim_no_erase() + log(2 / delta));
            std::cout << "Turbo Rough VCdim sample size with espilon "<<epsilon << ",  delta "<< delta <<", radius "<< distance_threshold<< " is: "<< sample_size <<std::endl;
            //Step 2: assert sampling is worthwhile
            if(sample_size > the_trajectory.num_trajectories()){

                std::cerr << "VC Bound was too loose for your dataset."<< std::endl;

                std::exit(1);

            }
            //Step 3: Sample indexes with replacement
            this->sample_trajectories(sample_size);

        }
        
        void generate_pathlets_vc_dim_sample(){
            //Step 1: compute sample size according to rule. 
            (this->sampled_trajs_ids).clear();

            int sample_size = (int) (2 / (epsilon * epsilon)) * (this->pathlets_vc_dim() + log(2 / delta));
            std::cout << " Pathlets VCdim sample size with espilon "<<epsilon << ",  delta "<< delta <<", radius "<< distance_threshold<< " is: "<< sample_size <<std::endl;
            //Step 2: assert sampling is worthwhile
            if(sample_size > the_trajectory.num_trajectories()){

                std::cerr << "VC Bound was too loose for your dataset."<< std::endl;

                std::exit(1);

            }
            //Step 3: Sample indexes with replacement
            this->sample_trajectories(sample_size);
        }
        void dump_sample_to_file(std::string filename){

            //std::cout << "Started dumping the sample to a file"<< std::endl;
            assert(!sampled_trajs_ids.empty());
            std::ofstream fout(filename);
            for (int j = 0; j< sampled_trajs_ids.size(); j++){
                //std::cout << "I am printing the sample "<< j <<std::endl;
                this->print_subtrajectory_to_file(fout, sampled_trajs_ids.at(j));

            }

            return;
        }

        void generate_fixed_size_sample(int size){

            (this->sampled_trajs_ids).clear();

            if(size > the_trajectory.num_trajectories()){

                std::cerr << "VC Bound was too loose for your dataset."<< std::endl;

                std::exit(1);

            }
            this->sample_trajectories(size);
        }

        void generate_rough_vc_no_erase_sample(){

            //Step 1: compute sample size according to rule. 
            (this->sampled_trajs_ids).clear();

            int sample_size = (int) (2 / (epsilon * epsilon)) * (this->rough_vc_dim_no_erase() + log(2 / delta));
            std::cout << " Rough VCdim sample size with espilon "<<epsilon << ",  delta "<< delta <<", radius "<< distance_threshold<< " is: "<< sample_size <<std::endl;
            //Step 2: assert sampling is worthwhile
            if(sample_size > the_trajectory.num_trajectories()){

                std::cerr << "VC Bound was too loose for your dataset."<< std::endl;

                std::exit(1);

            }
            //Step 3: Sample indexes with replacement
            this->sample_trajectories(sample_size);

        }
        //For debuggining
        std::map<index_t, std::pair<point_t, std::vector<index_t>>> get_pathlet_beginnings(){

            return this->pathlet_beginnings;
        }
        
        void generate_pathlet_aware_vc_dim_sample(){
            //Step 1: compute sample size according to rule. 
            (this->sampled_trajs_ids).clear();

            int sample_size = (int) (2 / (epsilon * epsilon)) * (this->pathlet_aware_vc_dim() + log(2 / delta));
            std::cout << " Pathlet-aware VCdim sample size with espilon "<<epsilon << ",  delta "<< delta <<", radius "<< distance_threshold<< " is: "<< sample_size <<std::endl;
            //Step 2: assert sampling is worthwhile
            if(sample_size > the_trajectory.num_trajectories()){

                std::cerr << "VC Bound was too loose for your dataset."<< std::endl;

                std::exit(1);

            }
            //Step 3: Sample indexes with replacement
            this->sample_trajectories(sample_size);
        }

    private:

        int pathlet_aware_vc_dim(){
            

            //Compute VC Dimension 
            range_search_t search{the_pathlets[0], grid_side_factor*distance_threshold};
        
            //Informative insertion 
            for (const auto& pair : this->pathlet_beginnings){
                index_t beginning_id = pair.first;
                point_t location = pair.second.first;
                for (auto end_id :  pair.second.second){
                    
                    search.insert(beginning_id, location, end_id);
                }
            }
            
            range_search_t search_ends{the_pathlets[0], grid_side_factor*distance_threshold};

            //if(thorough && min_length >=3){
                //Informative insertion 
            for (const auto& pair : this->pathlet_ends){
                point_t location = pair.second.first;
                index_t end_id = pair.first;
                search_ends.insert(end_id, location);
                
            }
                //assert(search_ends.num_elements() == search.num_elements());
            //}
            /*
            range_search_t search_seconds{the_pathlets[0], grid_side_factor*distance_threshold};
            if(thorough && min_length >=4){
                for (const auto& pair : this->pathlet_beginnings){

                search_seconds.insert(pair.first,the_pathlets[pair.first +1]);

                }
                assert(search_seconds.num_elements() == search.num_elements());
            }

            range_search_t search_thirds{the_pathlets[0], grid_side_factor*distance_threshold};
            if(thorough && min_length >=5){
                for (const auto& pair : this->pathlet_beginnings){

                search_thirds.insert(pair.first,the_pathlets[pair.first +2]);

                }
                assert(search_thirds.num_elements() == search.num_elements());
            }
            
            */
            assert(search.num_elements() == this->pathlet_beginnings.size());
            //Encoding a pathlet <beginning, end> as beginning*10^x +end, decoding as beignning = encoding / 10^x, end=encoding %10^x
            index_t maximum_extreme = the_pathlets.total_size();
            index_t temp = maximum_extreme;
            //std::cout<< std::format("Maximum extreme is : {}\n", maximum_extreme);
            int x = 1;
            while (true){
                
                maximum_extreme /= 10;
                //std::cout<< std::format("Maximum extreme is : {}\n", maximum_extreme);
                if(maximum_extreme == 0){
                    break;
                }
                x++;

            } 


            //std::cout << std::format("The maximum value for an extreme is {}, hence x is {}.\n", temp, x);
            std::vector<int> c;
            std::vector<std::pair<id_t, RangeRoaringBitmap64<space>>> inverted_index; //RRBm64 becuse it keeps pathlet ids.
            //FILL UP with empty sets to prevent weird behaviour with 0 or 1 beginning trajectories. 
            //Everything crashes if trajectories are not sorted or not consecutively named.
            inverted_index.reserve(the_trajectory.num_trajectories_not_consecutive()+1);
            for (unsigned int i = 0; i <the_trajectory.num_trajectories_not_consecutive()+1; i++){

                inverted_index.push_back(std::pair<id_t, RangeRoaringBitmap64<space>>{i,{}});

            }
            std::map<long unsigned int, RangeRoaringBitmap<space>> columns_map;
            index_t last_seen_trajectory = the_trajectory.get_id_at(0);
            std::map<index_t, std::vector<index_t>> traj_set_beginnings;
            //RangeRoaringBitmap<space> traj_set_seconds(std::set<id_t>({}));
            //RangeRoaringBitmap<space> traj_set_thirds(std::set<id_t>({}));
            RangeRoaringBitmap<space> traj_set_ends(std::set<id_t>({}));
            float squared_distance_threshold = distance_threshold *distance_threshold;
            for(index_t i =0; i<=the_trajectory.get_actual_size(); i++){
                //if(i%10000 == 0){
                    
                //std::cout<< "Processing point "<< i<< " to find the c bound" << std::endl;

                //}
                if(the_trajectory.get_id_at(i) == last_seen_trajectory){

                    point_t point = the_trajectory[i];
                    

                    for (const auto idx: search.search_and_return_associated_ids(the_trajectory[i], this->distance_threshold*distance_threshold, the_trajectory.get_id_at(i))) {

                        traj_set_beginnings.insert(idx);

                    }

                    //if(thorough && min_length >=3){

                    for (const auto idx: search_ends.search(the_trajectory[i], this->distance_threshold*distance_threshold, the_trajectory.get_id_at(i))) {
                        //std::cout << std::format("I am inserting end index {}.\n", idx);
                        traj_set_ends.add(idx);

                    }
                    //}
                }
                else{
                    
                    // Append the result up to now to c
                    //if(thorough &&min_length >=3)
                    if(last_seen_trajectory>=0){
                        //std::cout << std::format("For last seen trajectory {} trajset beginnings has {} entries \n", last_seen_trajectory, traj_set_beginnings.size());
                        //std::cout << std::format("For last seen trajectory {} trajset ends has {} entries \n", last_seen_trajectory, traj_set_ends.count());
                        std::set<std::pair<index_t,index_t>> surviving_pathlets;
                        surviving_pathlets = hierarchical_counting_and_reporting(traj_set_beginnings, traj_set_ends);
                        //std::cout << std::format("For Last Seen Trajectory {} the number of surviving pathlets is : {}\n",last_seen_trajectory, surviving_pathlets.size());
                        //std::cout.flush();
                        c.push_back(floor(log2(surviving_pathlets.size()) + 1));
                        //encode surviving pathlets and write them in the inverted index.
                        for (const auto& pathlet : surviving_pathlets){
                            auto pid = encode_pathlet(pathlet, x);
                            columns_map[pid].add(last_seen_trajectory);
                            inverted_index[last_seen_trajectory].second.add(pid);
                        }

                    }

                    //}
                    //else{
                        
                    //    int num_distinct_pathlets = std::accumulate(traj_set_beginnings.begin(),traj_set_beginnings.end(),size_t{0},[&](size_t acc, const auto& kv){return acc + kv.second.size();});
                        //std::cout << std::format("Trajectory {} is close to {} pathlet beginnings and in this case to {} distinct pathlets.\n", last_seen_trajectory, traj_set_beginnings.size(), num_distinct_pathlets);
                    //    c.push_back(floor(log2(num_distinct_pathlets) + 1));
                    //}
                    
                    //initialize the set again 
                    traj_set_beginnings.clear();
                    traj_set_ends.clear();
                    last_seen_trajectory = the_trajectory.get_id_at(i);
                    //add info for the current point
                    for (const auto idx: search.search_and_return_associated_ids(the_trajectory[i], this->distance_threshold*distance_threshold, the_trajectory.get_id_at(i))) {

                        traj_set_beginnings.insert(idx);

                    }

                    

                    for (const auto idx: search_ends.search(the_trajectory[i], this->distance_threshold*distance_threshold,the_trajectory.get_id_at(i))) {

                        traj_set_ends.add(idx);

                    }

                    
                }

            }
            //PRINT THE DATA STRUCTURES TO CHECK 
            // Print inverted_index
            
            
            //std::cout << "---------------------------------Printing inverted_index...\n";
            //for (id_t i = 0; i <inverted_index.size(); i++){

            //    std::cout << std::format("Trajectory with id {} has {} potentially matching pathlets, which are {}\n", inverted_index[i].first, inverted_index[i].second.count(), inverted_index[i].second.to_string());

            //}
            
            std::cout<< "I finished adding\n";
            
            std::sort(inverted_index.begin(), inverted_index.end(), [](const auto& a, const auto& b) {return b.second < a.second;});
            
            
            //std::cout << "---------------------------------Printing SORTED inverted_index...\n";
            //for (id_t i = 0; i <inverted_index.size(); i++){

            //    std::cout << std::format("Trajectory with id {} has {} potentially matching pathlets, which are {}\n", inverted_index[i].first, inverted_index[i].second.count(), inverted_index[i].second.to_string());

            //}
            
            
            

            // Print map data structure
            /*
            std::cout << "---------------------------------Printing columns map...\n";
            for (const auto& entry : columns_map){

                unsigned int pid = entry.first;
                RangeRoaringBitmap<space> column(entry.second.get_range()); //inefficient but this is just for debugging 
                std::cout << std::format("Pathelt with pid {} appears in {} trajectories, which are {}\n", pid, column.count(), column.to_string());
            }     
            
            */       
            

            std::sort(c.begin(),c.end(), std::greater<>());
            
            //std::cout << "---------------------------------Printing SORTED inverted_index...\n";
            //for (id_t i = 0; i <inverted_index.size(); i++){

            //    std::cout << std::format("Trajectory with id {} has {} potentially matching pathlets, which are {}\n", inverted_index[i].first, inverted_index[i].second.count(), inverted_index[i].second.to_string());

            //}            
            
            //Now I start with the algorithm 
            int vc_dim = 0;
            
            for(int i = 0 ; i < c.size(); i++){

                if(vc_dim < c.at(i)){

                    vc_dim++;

                }

            }

            // Okay now i start with the vc_dim
            int k = vc_dim;
            while(true){

                //I identify c_k as the deepest trajectory which contains at least 2^k -1 patterns 
                int c_k = 0;
                RangeRoaringBitmap<space> tids_for_a_shatterable_set{};
                for (int idx = 0; idx< inverted_index.size(); idx++){
                    if(inverted_index[idx].second.count()>= POWERS_OF_TWO[k-1]-1){ 
                        c_k = idx;
                        tids_for_a_shatterable_set.add(inverted_index[idx].first);
                    }
                    else{
                        break;
                    }
                }
                std::cout << std::format("I am trying to exclude VC-dim = {}:  c_k is {}.\n", k, c_k);
                std::cout << std::format("The candidate tids for this case are {}, which are {}\n", tids_for_a_shatterable_set.count(), tids_for_a_shatterable_set.to_string());
                std::cout << std::format("I now build the restricted pathlet map...\n");
                std::map<long unsigned int, RangeRoaringBitmap<space>> restricted_column_map;
            
                for (const auto& item : columns_map){
                    
                    RangeRoaringBitmap<space> intersection(columns_map[item.first].get_range() & tids_for_a_shatterable_set.get_range());

                    if((intersection).count() >0 ){
                        restricted_column_map[item.first] = intersection;
                    }

                }
                //std::cout << std::format( "PRINTING the restricted map for k = {}...\n",k);

                //for (const auto& entry : restricted_column_map){
                //    unsigned int pid = entry.first;
                //    RangeRoaringBitmap<space> column(entry.second.get_range()); //inefficient but this is just for debugging 
                //    std::cout << std::format("Pathlet with pid {} appears in {} trajectories, which are {}\n", pid, column.count(), column.to_string());
                //}

                //std::cout << std::format("Converting the map to vector. (Recall you cannot remove duplicate sets unless they are VERY MUCH REPEATED - i omit this for now to be extra conservative in the upper bound).\n");
                std::vector<RangeRoaringBitmap<space>> restricted_ranges_list; // I will be looking only at the cardinality 
                for (const auto entry : restricted_column_map){
                    restricted_ranges_list.push_back(entry.second); //inefficient but this is just for debugging 
                }
                std::sort(restricted_ranges_list.begin(), restricted_ranges_list.end(),[](const auto& a, const auto& b) {return b < a;});
                //std::cout << "-----------------PRINTING SORTED RANGES LIST \n";
                //for (const auto& item : restricted_ranges_list){
                //    std::cout << item.to_string()<< "\n";
                //}

                //Prune duplicate ranges from restricted ranges list.
                for (size_t i = 0; i < restricted_ranges_list.size(); ++i) {
                    int size_of_support = restricted_ranges_list[i].count();
                    size_t target_position = i+POWERS_OF_TWO[size_of_support]-2; 
                    if(size_of_support ==1){
                        target_position = i+1;
                    }
                    if(target_position>=restricted_ranges_list.size()){
                            continue;
                    }
                    //If i have more than 2^{size_of_support}-1 copies of the support the exceeding copies are useless. 
                    //std::cout <<std::format("I have a support of size {} at position {}, so I need to check from position {} onwards.\n", size_of_support, i,target_position);
                    //std::cout.flush();
                    while(restricted_ranges_list[target_position]==restricted_ranges_list[i]){
                        //std::cout << std::format("I am erasing support {}.\n", (*(restricted_ranges_list.begin()+target_position)).to_string());
                        restricted_ranges_list.erase(restricted_ranges_list.begin() + target_position );
                    }
                }
                
                //std::cout << "-----------------PRINTING PRUNED SORTED RANGES LIST \n";
                //for (const auto& item : restricted_ranges_list){
                //    std::cout << item.to_string()<< "\n";
                //}

                bool missing_pairs = true;
                bool failed_pair_check = false;
                //PRUNE VIA PAIRS: i need at least k choose 2 pairs that appear in at least 2^{k-2} supports in the restricted ranges list.
                if(k>=3 && c_k<50){
                    
                    std::map<std::pair<id_t,id_t>, int> candidate_pairs;
                    int actually_usable_pairs  = 0;
                    int presence_threshold = POWERS_OF_TWO[k-2];
                    //std::cout<< std::format("POWERS_OF_TWO[0]={}, POWERS_OF_TWO[{}]={}\n ", POWERS_OF_TWO[0], k-2, POWERS_OF_TWO[k-2]);
                    int num_usable_pairs_threshold = binom(k, 2);
                    
                    for (size_t idx = 0; idx < restricted_ranges_list.size() && missing_pairs; ++idx) {
                        const auto& s = restricted_ranges_list[idx];

                        std::vector<id_t> tids;
                        tids.reserve(s.get_range().cardinality());
                        for (auto tid : s.get_range())
                            tids.push_back(tid);

                        for (int i = 0; i < (int)tids.size() && missing_pairs; i++) {
                            for (int j = i + 1; j < (int)tids.size() && missing_pairs; j++) {

                                std::pair<id_t,id_t> pair{tids[i], tids[j]};
                                candidate_pairs[pair]++;
                                if(candidate_pairs[pair]==presence_threshold){
                                    actually_usable_pairs++;
                                }
                                if(actually_usable_pairs >=num_usable_pairs_threshold){
                                    missing_pairs = false;
                                }

                            }
                        }

                    }
                    if (missing_pairs){
                        std::cout << std::format("I failed the pair check. For actually usable pairs {} and required threshold {}, with presence threshold {}.\n", actually_usable_pairs, num_usable_pairs_threshold, presence_threshold);

                        failed_pair_check = true;
                    }

                }




                // Now the first items in restricted_ranges_list will be sorted in decreasing order. 
                // If the cardinalities can carry a set of size k...  
                int reserved_pathlets_for_current_visiting_size = 0;
                int current_visiting_size = k;
                int needed_pathlets_for_current_visting_size = binom(k, k-current_visiting_size);
                
                bool unshatterable = false;
                for (int i = 0; i< restricted_ranges_list.size() ; i++){
                    int c_val = restricted_ranges_list[i].count(); //std::floor(log2(restricted_ranges_list[i].count())-1);
                    if(c_val>= current_visiting_size){
                        //it is good to keep for the shattering
                        //std::cout << std::format("I have found a pathlet whose cardinality is {}, for current visiting size of {} i needed {}. \n", c_val, current_visiting_size,current_visiting_size );
                        needed_pathlets_for_current_visting_size--;
                    }
                    else{
                        unshatterable = true;
                        break;
                    }
                    if(needed_pathlets_for_current_visting_size==0){
                        
                        current_visiting_size--;
                        needed_pathlets_for_current_visting_size = binom(k, k-current_visiting_size);
                        if(current_visiting_size == 0){
                            //std::cout << std::format("According to cardinallites a set of size {} can be shattered.\n", k);
                            break;
                        }
                    }
                    if(unshatterable){
                        break;
                    }

                }
                
                if(unshatterable || c_k == 0 || current_visiting_size>0 || failed_pair_check){
                    
                    std::cout<< std::format("With the restricted pathlets I cannot shatter a set with cardinality as large as {}. Hence, I can try the lower value as a better upper bound.\n", k);
                    k--;
                }
                else{
                    std::cout <<std::format("According to cardinalities I can shatter sets of size {}, so I report it as an upper bound to the vc-dim.\n", k);
                    break;
                }


            }            

            //std::cout <<std::format( "I could not shatter anything of size larger than {}. The standard unaware upper bound was instead {}\n", k, vc_dim);
            

            std::cout<<"OLD VC DIM ESTIMATE WAS "<<vc_dim<< "\n";
            std::cout <<"VC DIM ESTIMATE IS "<< k <<"\n";
            return k;

        }
        
        int rough_vc_dim_no_erase(){

            range_search_t search{the_pathlets[0], grid_side_factor*distance_threshold};

            for (const auto& pair : this->pathlet_beginnings){

                point_t location = pair.second.first;
                index_t beginning_id = pair.first;
                for (auto end_id : pair.second.second){
                    search.insert(beginning_id, location, end_id);
                }
                
            }
            
            range_search_t search_ends{the_pathlets[0], grid_side_factor*distance_threshold};

            if(thorough && min_length >=3){
                for (const auto& pair : this->pathlet_ends){
                    point_t location = pair.second.first;
                    index_t beginning_id = pair.first;
                    for (auto end_id : pair.second.second){
                        search.insert(beginning_id, location, end_id);
                    }
                }
                //assert(search_ends.num_elements() == search.num_elements());
            }
            int d = 0;
            std::set<transaction_and_points_pair> c;
            index_t last_seen_trajectory = the_trajectory.get_id_at(0);
            int counter=0;
            std::vector<std::pair<index_t, std::vector<index_t>>> retrieved_points_with_duplicates;
            int counter_ends =0; 
            int total_distances = 0;
            float squared_distance_threshold = distance_threshold*distance_threshold;
            for(index_t i =0; i<=the_trajectory.get_actual_size(); i++){
                //if(i%10000 == 0){
                    
                //std::cout<< "Processing point "<< i<< " to find the c bound" << std::endl;

                //}
                if(the_trajectory.get_id_at(i) == last_seen_trajectory){
                    
                    
                    std::map<index_t, std::vector<index_t>> retrieved_points = search.search_and_return_associated_ids_no_erase(the_trajectory[i], this->distance_threshold*distance_threshold,the_trajectory.get_id_at(i));
                    counter +=std::accumulate(retrieved_points.begin(),retrieved_points.end(),size_t{0},[&](size_t acc, const auto& kv){return acc + kv.second.size();});
                    retrieved_points_with_duplicates.insert(retrieved_points_with_duplicates.end(), retrieved_points.begin(), retrieved_points.end());
                    if(thorough && min_length >=3){
                        std::map<index_t, std::vector<index_t>> retrieved_ends = search.search_and_return_associated_ids_no_erase(the_trajectory[i], this->distance_threshold*distance_threshold, the_trajectory.get_id_at(i));
                        counter_ends +=std::accumulate(retrieved_ends.begin(),retrieved_ends.end(),size_t{0},[&](size_t acc, const auto& kv){return acc + kv.second.size();});
                    }

                    

                }
                else{

                    int min_c = counter;
                    if(thorough && min_length >= 3){
                        
                        min_c = std::min(counter, counter_ends);
                    }
                    // Append the result up to now to c
                    if(floor(log2(min_c*(1.0-(1/min_c))) + 1) > d) {

                        c.insert(transaction_and_points_pair{last_seen_trajectory, retrieved_points_with_duplicates,floor(log2(min_c*(1.0-(1/min_c))) + 1)});
                        int l_prime = prev(c.end())->capacity_ub;
                        if(l_prime >d){
                            d++;
                        }
                        else{

                            c.erase(prev(c.end()));

                        }
                    }
                    //c.push_back(floor(log2(min_c*(1.0-(1/min_c))) + 1));
                    
                    //initialize the set again 
                    counter = 0;
                    counter_ends = 0;
                    retrieved_points_with_duplicates.clear();
                    last_seen_trajectory = the_trajectory.get_id_at(i);
                    std::map<index_t, std::vector<index_t>> retrieved_points = search.search_and_return_associated_ids_no_erase(the_trajectory[i], this->distance_threshold*distance_threshold, the_trajectory.get_id_at(i));
                    counter +=std::accumulate(retrieved_points.begin(),retrieved_points.end(),size_t{0},[&](size_t acc, const auto& kv){return acc + kv.second.size();});
                    retrieved_points_with_duplicates.insert(retrieved_points_with_duplicates.end(), retrieved_points.begin(), retrieved_points.end());
                    if(thorough && min_length >=3){
                        std::map<index_t, std::vector<index_t>> retrieved_ends = search.search_and_return_associated_ids_no_erase(the_trajectory[i], this->distance_threshold*distance_threshold, the_trajectory.get_id_at(i));
                        counter_ends +=std::accumulate(retrieved_ends.begin(),retrieved_ends.end(),size_t{0},[&](size_t acc, const auto& kv){return acc + kv.second.size();});
                    }
                    //add info for the current point


                }

            }  
            //std::cout<<total_distances << std::endl;
            //std::sort(c.begin(),c.end(), std::greater<>());
            //std::cout <<"############ Details: #################\n";
            /*
            for (int i=0; i<c.size();i++){

                std::cout<< " H index vector at position "<< i<< " "<< c.at(i)<<std::endl;

            }
            */
            
            int vc_dim = prev(c.end())->capacity_ub;
            /*
             for(int i = 0 ; i < c.size(); i++){

                if(vc_dim < c.at(i)){

                    vc_dim++;

                }

            }
            
            */
           
            std::cout <<"VC DIM ESTIMATE IS "<< vc_dim <<"\n";
            if(inspect){

                print_duplicate_statistics(c);

            }
            return vc_dim;

        }
        
        int turbo_rough_vc_dim_no_erase(){            
            range_search_t search{the_pathlets[0], grid_side_factor*distance_threshold};

            for (const auto& pair : this->pathlet_beginnings){

                point_t location = pair.second.first;
                index_t beginning_id = pair.first;
                for (auto end_id : pair.second.second){
                    search.insert(beginning_id, location, end_id);
                }
                
            }
            std::cout << std::format("I finished inserting the pathlet beginnings.\n");
            std::cout << std::format("Here is the grid: \n {}\n", search.get_grid().to_string());
            range_search_t search_ends{the_pathlets[0], grid_side_factor*distance_threshold};

            if(thorough && min_length >=3){
                for (const auto& pair : this->pathlet_ends){
                    point_t location = pair.second.first;
                    index_t beginning_id = pair.first;
                    for (auto end_id : pair.second.second){
                        search.insert(beginning_id, location, end_id);
                    }
                }
                //assert(search_ends.num_elements() == search.num_elements());
            }
            int d = 0;
            std::set<transaction_and_points_pair> c;
            index_t last_seen_trajectory = the_trajectory.get_id_at(0);
            int counter=0;
            int counter_ends =0; 
            std::vector<std::pair<long unsigned int, std::vector<long unsigned int> >> retrieved_points_with_duplicates;
            int total_distances = 0;
            float squared_distance_threshold = distance_threshold*distance_threshold;
            for(index_t i =0; i<=the_trajectory.get_actual_size(); i++){
                //if(i%10000 == 0){
                    
                //std::cout<< "Processing point "<< i<< " to find the c bound" << std::endl;

                //}
                if(the_trajectory.get_id_at(i) == last_seen_trajectory){
                    
                    std::cout << std::format("Processing point {}\n", i);
                    int retrieved_points = search.count_nearby_points_weighted_no_erase(the_trajectory[i], this->distance_threshold*distance_threshold,the_trajectory.get_id_at(i));
                    counter += retrieved_points;
                    if(thorough && min_length >=3){
                        int retrieved_ends = search_ends.count_nearby_points_weighted_no_erase(the_trajectory[i], this->distance_threshold*distance_threshold, the_trajectory.get_id_at(i));
                        counter_ends += retrieved_ends;
                    }

                }
                else{

                    int min_c = counter;
                    if(thorough && min_length >= 3){
                        
                        min_c = std::min(counter, counter_ends);
                    }
                    // Append the result up to now to c
                    if(floor(log2(min_c*(1.0-(1/min_c))) + 1) > d) {

                        c.insert(transaction_and_points_pair{last_seen_trajectory, retrieved_points_with_duplicates ,floor(log2(min_c*(1.0-(1/min_c))) + 1)});
                        int l_prime = prev(c.end())->capacity_ub;
                        if(l_prime >d){
                            d++;
                        }
                        else{

                            c.erase(prev(c.end()));

                        }
                    }
                    //c.push_back(floor(log2(min_c*(1.0-(1/min_c))) + 1));
                    
                    //initialize the set again 
                    counter = 0;
                    counter_ends = 0;
                    retrieved_points_with_duplicates.clear();
                    last_seen_trajectory = the_trajectory.get_id_at(i);
                    std::cout << std::format("PROCESSING TID {}\n", last_seen_trajectory);
                    std::cout << std::format("Processing point {}\n", i);
                    int retrieved_points = search.count_nearby_points_weighted_no_erase(the_trajectory[i], this->distance_threshold*distance_threshold, the_trajectory.get_id_at(i));
                    counter += retrieved_points;
                    
                    if(thorough && min_length >=3){
                        int retrieved_ends = search.count_nearby_points_weighted_no_erase(the_trajectory[i], this->distance_threshold*distance_threshold, the_trajectory.get_id_at(i));
                        counter_ends +=retrieved_ends;
                    }
                    //add info for the current point


                }

            }  
            //std::cout<<total_distances << std::endl;
            //std::sort(c.begin(),c.end(), std::greater<>());
            //std::cout <<"############ Details: #################\n";
            /*
            for (int i=0; i<c.size();i++){

                std::cout<< " H index vector at position "<< i<< " "<< c.at(i)<<std::endl;

            }
            */
            
            int vc_dim = prev(c.end())->capacity_ub;
            /*
             for(int i = 0 ; i < c.size(); i++){

                if(vc_dim < c.at(i)){

                    vc_dim++;

                }

            }
            
            */
           
            std::cout <<"VC DIM ESTIMATE IS "<< vc_dim <<"\n";
            
            return vc_dim;
        }
        
        int rough_vc_dim(){
            range_search_t search{the_pathlets[0], grid_side_factor*distance_threshold};

            for (const auto& pair : this->pathlet_beginnings){

                point_t location = pair.second.first;
                index_t beginning_id = pair.first;
                for (auto end_id : pair.second.second){
                    search.insert(beginning_id, location, end_id);
                }
            
            }
        
            
            range_search_t search_ends{the_pathlets[0], grid_side_factor*distance_threshold};

            if(thorough && min_length >=3){
                for (const auto& pair : this->pathlet_ends){
                    point_t location = pair.second.first;
                    index_t beginning_id = pair.first;
                    for (auto end_id : pair.second.second){
                        search.insert(beginning_id, location, end_id);
                    }
                }
                assert(search_ends.num_elements() == search.num_elements());
            }

            
            std::vector<int> c;
            index_t last_seen_trajectory = the_trajectory.get_id_at(0);
            int counter=0;
            int counter_ends = 0;
            int total_distances = 0;
            float squared_distance_threshold = distance_threshold*distance_threshold;
            for(index_t i =0; i<=the_trajectory.get_actual_size(); i++){
                //if(i%10000 == 0){
                    
                //std::cout<< "Processing point "<< i<< " to find the c bound" << std::endl;

                //}
                if(the_trajectory.get_id_at(i) == last_seen_trajectory){
                    
                    
                    int ss = search.search(the_trajectory[i], this->distance_threshold*distance_threshold, the_trajectory.get_id_at(i)).size();
                    counter +=ss;
                    if(thorough && min_length >=3){
                        int ss_ends = search_ends.search(the_trajectory[i], this->distance_threshold*distance_threshold, the_trajectory.get_id_at(i)).size();
                        counter_ends += ss_ends;
                    }

                }
                else{

                    // Append the result up to now to c+
                    int min_c = counter;

                    if(thorough && min_length >= 3){
                        
                        min_c = std::min(counter, counter_ends);
                    }
                    c.push_back(floor(log2(min_c*(1.0-(1/min_c))) + 1));
                    //initialize the set again 
                    counter = 0;
                    counter_ends= 0;
                    last_seen_trajectory = the_trajectory.get_id_at(i);
                    int ss = search.search(the_trajectory[i], this->distance_threshold*distance_threshold, the_trajectory.get_id_at(i)).size();
                    counter += ss;
                    if(thorough && min_length >=3){
                        int ss_ends = search_ends.search(the_trajectory[i], this->distance_threshold*distance_threshold, the_trajectory.get_id_at(i)).size();
                        counter_ends += ss_ends;
                    }
                    //add info for the current point


                }

            }  
            //std::cout<<total_distances << std::endl;
            std::sort(c.begin(),c.end(), std::greater<>());
            //std::cout <<"############ Details: #################\n";
            /*
            for (int i=0; i<c.size();i++){

                std::cout<< " H index vector at position "<< i<< " "<< c.at(i)<<std::endl;

            }
            */

            int vc_dim = 0;
            
            for(int i = 0 ; i < c.size(); i++){

                if(vc_dim < c.at(i)){

                    vc_dim++;

                }

            }
            std::cout <<"VC DIM ESTIMATE IS "<< vc_dim <<"\n";
            return vc_dim;

        }

        int vc_dim(){
            std::cout << "Beginning vc dim computation"<<std::endl;
            //Compute VC Dimension 
            range_search_t search{the_pathlets[0], grid_side_factor*distance_threshold};
        
            //pathlet beginnings caan be used with non informative search
            for (const auto& pair : this->pathlet_beginnings){

                point_t location = pair.second.first;
                index_t beginning_id = pair.first;
                for (auto end_id : pair.second.second){
                    search.insert(beginning_id, location);
                }
                
            }
            
            range_search_t search_ends{the_pathlets[0], grid_side_factor*distance_threshold};

            //if(thorough && min_length >=3){
            for (const auto& pair : this->pathlet_ends){
                point_t location = pair.second.first;
                index_t end_id = pair.first;
                search_ends.insert(end_id, location);
                
            }
                //assert(search_ends.num_elements() == search.num_elements());
            //}

            assert(search.num_elements() == this->pathlet_beginnings.size());
            
            std::vector<int> c;
            index_t last_seen_trajectory = the_trajectory.get_id_at(0);
            RangeRoaringBitmap64<space> traj_set_beginnings;
            RangeRoaringBitmap64<space> traj_set_ends;
            //std::vector<index_t> traj_set_beginnings;
            //std::vector<index_t> traj_set_ends;
            float squared_distance_threshold = distance_threshold *distance_threshold;
            for(index_t i =0; i<=the_trajectory.get_actual_size(); i++){
                //if(i%10000 == 0){
                    
                //std::cout<< "Processing point "<< i<< " to find the c bound" << std::endl;

                //}
                if(the_trajectory.get_id_at(i) == last_seen_trajectory){

                    point_t point = the_trajectory[i];
                    

                    for (const auto idx: search.search(the_trajectory[i], this->distance_threshold*distance_threshold, the_trajectory.get_id_at(i))) {

                        traj_set_beginnings.add(idx);

                    }

                    //if(thorough && min_length >=3){
                    if(min_length >1){
                        
                        for (const auto idx: search_ends.search(the_trajectory[i], this->distance_threshold*distance_threshold, the_trajectory.get_id_at(i))) {
                        //std::cout << std::format("I am inserting end index {}.\n", idx);
                            traj_set_ends.add(idx);

                        }

                    }
                    //}
                }
                else{
                    
                    // Append the result up to now to c
                    //if(thorough && min_length >=3){
                    int num_distinct_pathlets = 0;
                    if(min_length >1){

                        num_distinct_pathlets = hierarchical_counting(traj_set_beginnings, traj_set_ends);
                    }
                    else{
                        num_distinct_pathlets = hierarchical_counting(traj_set_beginnings, traj_set_beginnings);
                    }
                    c.push_back(floor(log2(num_distinct_pathlets) + 1));

                    //}
                    //else{
                        
                    //    int num_distinct_pathlets = std::accumulate(traj_set_beginnings.begin(),traj_set_beginnings.end(),size_t{0},[&](size_t acc, const auto& kv){return acc + kv.second.size();});
                        //std::cout << std::format("Trajectory {} is close to {} pathlet beginnings and in this case to {} distinct pathlets.\n", last_seen_trajectory, traj_set_beginnings.size(), num_distinct_pathlets);
                    //    c.push_back(floor(log2(num_distinct_pathlets) + 1));
                    //}
                    
                    //initialize the set again 
                    traj_set_beginnings.clear();
                    traj_set_ends.clear();
                    last_seen_trajectory = the_trajectory.get_id_at(i);
                    //add info for the current point
                    for (const auto idx: search.search(the_trajectory[i], this->distance_threshold*distance_threshold, the_trajectory.get_id_at(i))) {

                        traj_set_beginnings.add(idx);

                    }

                    
                    if(min_length> 1){

                        for (const auto idx: search_ends.search(the_trajectory[i], this->distance_threshold*distance_threshold, the_trajectory.get_id_at(i))) {

                            traj_set_ends.add(idx);

                        }

                    }
                    

                    
                }

            }
            std::sort(c.begin(),c.end(), std::greater<>());
            
            int vc_dim = 0;
            
            for(int i = 0 ; i < c.size(); i++){

                if(vc_dim < c.at(i)){

                    vc_dim++;

                }

            }
            std::cout <<"VC DIM ESTIMATE IS "<< vc_dim <<"\n";

            return vc_dim;
        }
        
        //TODO: adapt it to nested pathlets? Now it does not support it
        int pathlets_vc_dim(){
            
            //Compute VC Dimension 
            range_search_t search{the_pathlets[0], grid_side_factor*distance_threshold, true};
        
            //HERE i should insert points from the trajectory in the map data structure

            for (int j = 0; j< the_trajectory.total_size(); j++){

                search.insert(the_trajectory.get_id_at(j), the_trajectory[j]);
            }
            
            range_search_t search_ends{the_pathlets[0], grid_side_factor*distance_threshold, true};
            //HERE i should insert points from the trajectory in the map data structure 
            for (int j = 0; j< the_trajectory.total_size(); j++){

                search_ends.insert(the_trajectory.get_id_at(j), the_trajectory[j]);
            }
            /*
            range_search_t search_seconds{the_pathlets[0], grid_side_factor*distance_threshold};
            if(thorough && min_length >=4){
                for (const auto& pair : this->pathlet_beginnings){

                search_seconds.insert(pair.first,the_pathlets[pair.first +1]);

                }
                assert(search_seconds.num_elements() == search.num_elements());
            }

            range_search_t search_thirds{the_pathlets[0], grid_side_factor*distance_threshold};
            if(thorough && min_length >=5){
                for (const auto& pair : this->pathlet_beginnings){

                search_thirds.insert(pair.first,the_pathlets[pair.first +2]);

                }
                assert(search_thirds.num_elements() == search.num_elements());
            }
            
            */
            //assert(search.num_elements() == this->pathlet_beginnings.size());
            
            std::vector<int> c;
            index_t last_seen_trajectory = the_trajectory.get_id_at(0);
            std::set<index_t> traj_set;
            //std::set<index_t> traj_set_seconds;
            //std::set<index_t> traj_set_thirds;
            std::set<index_t> traj_set_ends;
            float squared_distance_threshold = distance_threshold *distance_threshold;
            //here i should loop over the pathlets and 
            for(const auto& pair : this->pathlet_beginnings){
                
                traj_set.clear();
                traj_set_ends.clear();
                auto pathlet_idx = pair.first;
                auto pathlet_point = pair.second.first;
                for (const auto idx: search.search(pathlet_point, this->distance_threshold*distance_threshold, pathlet_idx)) {

                    traj_set.insert(idx);

                }
                if(thorough && min_length >= 3){

                    
                    for (const auto idx: search_ends.search(pathlet_point, this->distance_threshold*distance_threshold, pathlet_idx)) {

                    traj_set_ends.insert(idx);

                }

                
                /*
                
                if( thorough && min_length >=3){
                    c.push_back(floor(log2(std::count_if(traj_set.begin(), traj_set.end(), [&](const auto& x){ return traj_set_ends.contains(x); })) + 1));

                }
                else{
                    c.push_back(floor(log2(traj_set.size()) + 1));
                }
                */
               
                c.push_back(std::count_if(traj_set.begin(), traj_set.end(), [&](const auto& x){ return traj_set_ends.contains(x); }) + 1);

                
            }
            std::sort(c.begin(),c.end(), std::greater<>());
            //std::cout <<"############ Details: #################\n";
            
            //for (int i=0; i<c.size();i++){

            //    std::cout<< " H index vector at position "<< i<< " "<< c.at(i)<<std::endl;

            //}

            int vc_dim = 1;

            /*
            
            for(int i = 0 ; i < c.size(); i++){

                if(vc_dim < c.at(i)){

                    vc_dim++;

                }

            }
            
            */
            while(true){
                int reserved_pathlets = 0;
                for (int j = 0; j<vc_dim; j++)
                {   
                    for (auto it=c.begin(); it!=c.begin() +reserved_pathlets+binom(vc_dim, j); it++){
                        int c_val = *it;
                        if (c_val< vc_dim -j){
                            std::cout << std::format("I cannot pass the test for vc-dim equal to {}\n",vc_dim);
                            std::cout <<"VC DIM ESTIMATE IS "<< vc_dim-1 <<"\n";
                            return vc_dim;
                        }
                    }
                    
                    reserved_pathlets += binom(vc_dim, j);

                }
                std::cout<< "Passed test for vc_dim "<< vc_dim<<std::endl;
                vc_dim++;
            }
        }
    }
        
        int vc_dim_no_erase(){

            //Compute VC Dimension 
            range_search_t search{the_pathlets[0], grid_side_factor*distance_threshold};
        
            //pathlet beginnings caan be used with non informative search
            for (const auto& pair : this->pathlet_beginnings){

                point_t location = pair.second.first;
                index_t beginning_id = pair.first;
                for (auto end_id : pair.second.second){
                    search.insert(beginning_id, location);
                }
                
            }
            
            range_search_t search_ends{the_pathlets[0], grid_side_factor*distance_threshold};

            
            for (const auto& pair : this->pathlet_ends){
                point_t location = pair.second.first;
                index_t end_id = pair.first;
                search_ends.insert(end_id, location);
                
            }
            //assert(search_ends.num_elements() == search.num_elements());
            

            assert(search.num_elements() == this->pathlet_beginnings.size());
            
            std::vector<int> c;
            index_t last_seen_trajectory = the_trajectory.get_id_at(0);
            RangeRoaringBitmap64<space> traj_set_beginnings;
            RangeRoaringBitmap64<space> traj_set_ends;
            float squared_distance_threshold = distance_threshold *distance_threshold;
            for(index_t i =0; i<=the_trajectory.get_actual_size(); i++){
                if(i%10000 == 0){
                    
                std::cout<< "Processing point "<< i<< " to find the c bound" << std::endl;

                }
                if(the_trajectory.get_id_at(i) == last_seen_trajectory){

                    point_t point = the_trajectory[i];
                    

                    for (const auto idx: search.search_no_erase(the_trajectory[i], this->distance_threshold*distance_threshold, the_trajectory.get_id_at(i))) {

                        traj_set_beginnings.add(idx);

                    }

                    

                    for (const auto idx: search_ends.search_no_erase(the_trajectory[i], this->distance_threshold*distance_threshold, the_trajectory.get_id_at(i))) {

                        traj_set_ends.add(idx);

                    }
                    
                }
                else{
                    
                    // Append the result up to now to c
                    //if(thorough && min_length >=3){
                    int num_distinct_pathlets = hierarchical_counting(traj_set_beginnings, traj_set_ends);
                    c.push_back(floor(log2(num_distinct_pathlets) + 1));

                    //}
                    //else{
                        
                    //    int num_distinct_pathlets = std::accumulate(traj_set_beginnings.begin(),traj_set_beginnings.end(),size_t{0},[&](size_t acc, const auto& kv){return acc + kv.second.size();});
                        //std::cout << std::format("Trajectory {} is close to {} pathlet beginnings and in this case to {} distinct pathlets.\n", last_seen_trajectory, traj_set_beginnings.size(), num_distinct_pathlets);
                    //    c.push_back(floor(log2(num_distinct_pathlets) + 1));
                    //}
                    
                    //initialize the set again 
                    traj_set_beginnings.clear();
                    traj_set_ends.clear();
                    last_seen_trajectory = the_trajectory.get_id_at(i);
                    //add info for the current point
                    for (const auto idx: search.search_no_erase(the_trajectory[i], this->distance_threshold*distance_threshold, the_trajectory.get_id_at(i))) {

                        traj_set_beginnings.add(idx);

                    }

                    

                    for (const auto idx: search_ends.search_no_erase(the_trajectory[i], this->distance_threshold*distance_threshold, the_trajectory.get_id_at(i))) {

                        traj_set_ends.add(idx);

                    }

                    
                }

            }
            std::sort(c.begin(),c.end(), std::greater<>());
            
            int vc_dim = 0;
            
            for(int i = 0 ; i < c.size(); i++){

                if(vc_dim < c.at(i)){

                    vc_dim++;

                }

            }
            std::cout <<"VC DIM ESTIMATE IS "<< vc_dim <<"\n";
            return vc_dim;
        }

    /*
    int aggressive_vc_dim(){

        range_search_t search{the_trajectory};
        std::vector<int> c;
        index_t last_seen_trajectory = the_trajectory.get_id_at(0);
        std::set<index_t> traj_set;
        boost::container::flat_set<id_t> id_set;
        float sq_dist = distance_threshold * distance_threshold;
        int counter = 0;
        for(index_t i =0; i<=the_trajectory.get_actual_size(); i++){
            if(i%10000 == 0){
                
            std::cout<< "Processing point "<< i<< " to find the c bound" << std::endl;

            }
            if(the_trajectory.get_id_at(i) == last_seen_trajectory){
                for (index_t jj =0; jj <the_trajectory.total_size(); jj++){

                    auto d_ij = distance_function_t{}(the_trajectory[i], the_trajectory[jj]); //FI
                    if (d_ij <=  sq_dist){
                        counter++;
                        traj_set.insert(jj);
                        
                    }

                }
                
                for (const auto idx: search.search(i, this->distance_threshold)) {

                    traj_set.insert(idx);

                }
                
                
            }
            else{
                
                // Append the result up to now to c
                c.push_back(floor(log2(traj_set.size()) + 1));
                
                
                
                //initialize the set again 
                
                traj_set.clear();
                last_seen_trajectory = the_trajectory.get_id_at(i);
                //add info for the current point
                for (index_t jj =0; jj <the_trajectory.total_size(); jj++){

                    auto d_ij = distance_function_t{}(the_trajectory[i], the_trajectory[jj]); //FI
                    if (d_ij <=  sq_dist){
                        counter++;
                        //traj_set.insert(jj);
                        
                    }

                }


            }

        }
        std::sort(c.begin(),c.end(), std::greater<>());
        int vc_dim = 0;
        
        for(int i = 0 ; i < c.size(); i++){

            if(vc_dim < c.at(i)){

                vc_dim++;

            }

        }

        return vc_dim;

    

    }
    */
    
        void print_subtrajectory_to_file(std::ofstream& fout, id_t& id){

            size_t n = the_trajectory.num_trajectories();
            index_t j = the_trajectory.get_first_point_in_trajectory(id%n);
            
            while( the_trajectory.get_id_at(j) == id%n){
                
                fout << the_trajectory[j].x()<< " "<< the_trajectory[j].y()<<" "<< id<<std::endl;
                j++;
            }

            return;
        }

        int total_pathlet_number_respecting_ids(){
            int total = 0;

            for (id_t i = 0; i < the_pathlets.num_trajectories(); i ++){

                total += 2 * ceil(( the_pathlets.get_trajectory_size(i) / min_length));

            }

            return total;

        }
    
        void sample_trajectories(int sample_size){

        if(! sampled_trajs_ids.empty()){

            sampled_trajs_ids.clear();

        }
        std::vector<id_t> indexes;
        int n = the_trajectory.num_trajectories();
        for (int i = 1; i< n; i++){

            indexes.push_back(i);

        }
        
        std::ranges::shuffle(indexes, this->mt);
        // Extract sampled ids
        for (int i = 0; i < sample_size; i++){

            sampled_trajs_ids.push_back(indexes.at(i));

        }

        std::sort(sampled_trajs_ids.begin(),sampled_trajs_ids.end());

        id_t last_read_trajectory = sampled_trajs_ids.at(0);
        int repetitions = 0;
        for( int j = 1; j< sample_size; j++){
            
            if (sampled_trajs_ids.at(j) == last_read_trajectory){

                repetitions++;
                sampled_trajs_ids.at(j) += (id_t) (repetitions * n);

            }
            else{
                last_read_trajectory = sampled_trajs_ids.at(j);
                repetitions = 0;
            }
            if(repetitions > 0){

                std::cout<< "I have repetitions despite the shuffle"<< std::endl;
            }
        }
        //std::cout<< "I have finished the method to get the ids"<< std::endl;
        return;
    }

        void print_duplicate_statistics(std::set<transaction_and_points_pair>& c){

            for (auto& pair : c){
                //sort the retrieved points for each participating transaction
                std::vector<std::pair<index_t, std::vector<index_t>>> replica = pair.closeby_points;
                std::sort(replica.begin(), replica.end(), std::greater<>()); // this is indifferent
                // Now for each pathlet count the #of occurrences...
                index_t last_seen = -1;
                int occurrences = 0;
                std::vector<int> counters;
                for (auto pair : replica){
                    auto pt = pair.first;
                    if (pt != last_seen){
                        if (occurrences !=0){
                            
                            counters.push_back(occurrences*pair.second.size());

                        }
                        occurrences = 1;
                        last_seen = pt;

                    }
                    else{
                        occurrences++;
                    }
                }
                std::sort(counters.begin(), counters.end(), std::greater<>());
                float mean = std::reduce(counters.begin(), counters.end()) / ((float)(counters.size()));
                float empirical_variance = std::transform_reduce(counters.begin(), counters.end(), 0.0, std::plus<>(), [&](int x) {return (x-mean)*(x-mean);})/ ((float)(counters.size()-1));
                std::cout << "********************************\n";
                std::cout << std::format("TRANSACTION : {}\n", pair.traj_id);
                std::cout << std::format("NUMBER OF REPORTED POINTS : {}\n", replica.size());
                std::cout << std::format("AVERAGE NUMBER OF REPETITIONS : {}\n", mean);

                std::cout << std::format("TOP 5 NUMBER OF REPETITIONS : {}, {}, {}, {}, {}\n", counters.at(0), counters.at(1), counters.at(2), counters.at(3), counters.at(4) );

            }

            return;
        }
        static long long binom(int n, int k) {
            if (k > n - k) k = n - k;
            long long res = 1;
            for (int i = 0; i < k; ++i) {
                res = res * (n - i) / (i + 1);
            }
            return res;
        }
        
        std::set<std::pair<index_t,index_t>> hierarchical_counting_and_reporting(const std::map<index_t, std::vector<index_t>>& traj_set_beginnings, const RangeRoaringBitmap<space>& traj_set_ends ){


            std::set<std::pair<index_t,index_t>> surviving_pathlets; 
            std::set<std::pair<index_t,index_t>> surviving_pathlets_by_traj; 
            //roaring::Roaring pathlets_tids;
            //for(const auto idx : traj_set_ends){

            //    pathlets_tids.add(the_pathlets.get_id_at(idx));

            //}
            RangeRoaringBitmap<space> visited_pathlet_mothers(std::set<id_t>{});

            for (const auto& pair: traj_set_beginnings){
                //std::cout<< std::format("For beginning {} I have {} potential endpoints.\n", pair.first, pair.second.size());
                //auto tid = the_pathlets.get_id_at(pair.first);
                auto beginning = pair.first;
                auto pathlet_mother = the_pathlets.get_id_at(beginning);
                if(!visited_pathlet_mothers.get_range().contains(pathlet_mother)){
                    //Get possible pathlets at the lower level 
                    trajectory_t sliced_trajectory = the_pathlets.slice_trajectory_by_id(pathlet_mother);
                    BinaryPathletTree pathlet_tree(sliced_trajectory, pathlet_mother, floor(log2(sliced_trajectory.total_size())) + 1,min_length);
                    //For debugging purposes: print the tree
                    //std::cout<< "Pathelet tree for tid "<< pathlet_mother<< std::endl;
                    //std::cout << pathlet_tree.toString()<< std::endl;
                    
                    std::vector<PathletNode<space>> min_length_pathlets_low_level = pathlet_tree.getMinLengthPathletsAtLowLevels(this-> min_length);
                    std::vector<PathletNode<space>> min_length_pathlets= pathlet_tree.getMinLengthPathlets(this-> min_length);
                    //std::cout<< "CALLING HERE for pm "<< current_visiting_id << std::endl;
                    index_t trajectory_offset = the_pathlets.get_first_point_in_trajectory(pathlet_mother);
                    //NOw for each low level node I look if it is present, otherwise I delete it from the min_length_pathlets. 
                    for (size_t i = 0; i< min_length_pathlets_low_level.size(); i++){
                        auto pathlet_beginning = min_length_pathlets_low_level[i].pathlet.first;
                        auto pathlet_end = min_length_pathlets_low_level[i].pathlet.second;
                        if(!(traj_set_beginnings.contains(pathlet_beginning + trajectory_offset) && traj_set_ends.get_range().contains(pathlet_end + trajectory_offset))){
                            std::erase_if(min_length_pathlets, [&](const PathletNode<space>& pn){return min_length_pathlets_low_level[i].is_contained_by(pn);});
                        }
                    }
                    for (const auto alive_pn : min_length_pathlets){
                        surviving_pathlets.insert({trajectory_offset + alive_pn.pathlet.first, trajectory_offset + alive_pn.pathlet.second});
                    }
                    visited_pathlet_mothers.add(pathlet_mother);
                }


            }
            return surviving_pathlets;
        }

        int hierarchical_counting(const RangeRoaringBitmap64<space>& traj_set_beginnings, const RangeRoaringBitmap64<space>& traj_set_ends ){

            std::set<std::pair<index_t,index_t>> surviving_pathlets; 
             
            //roaring::Roaring pathlets_tids;
            //for(const auto idx : traj_set_ends){

            //    pathlets_tids.add(the_pathlets.get_id_at(idx));

            //}
            RangeRoaringBitmap<space> visited_pathlet_mothers(std::set<id_t>{});

            for (const auto& beginning: traj_set_beginnings.get_range()){
                //std::cout<< std::format("For beginning {} I have {} potential endpoints.\n", pair.first, pair.second.size());
                //auto tid = the_pathlets.get_id_at(pair.first);
                
                auto pathlet_mother = the_pathlets.get_id_at(beginning);
                if(!visited_pathlet_mothers.get_range().contains(pathlet_mother)){
                    //Get possible pathlets at the lower level 
                    //trajectory_t sliced_trajectory = the_pathlets.slice_trajectory_by_id(pathlet_mother);
                    //BinaryPathletTree pathlet_tree = the_pathlet_trees[pathlet_mother];
                    //For debugging purposes: print the tree
                    //std::cout<< "Pathelet tree for tid "<< pathlet_mother<< std::endl;
                    //std::cout << pathlet_tree.toString()<< std::endl;
                    
                    std::vector<PathletNode<space>> min_length_pathlets_low_level =  the_pathlet_trees[pathlet_mother].getMinLengthPathletsAtLowLevels(this-> min_length);
                    std::vector<PathletNode<space>> min_length_pathlets=  the_pathlet_trees[pathlet_mother].getMinLengthPathlets(this-> min_length);
                    //std::cout<< "CALLING HERE for pm "<< current_visiting_id << std::endl;
                    index_t trajectory_offset = the_pathlets.get_first_point_in_trajectory(pathlet_mother);
                    //NOw for each low level node I look if it is present, otherwise I delete it from the min_length_pathlets. 
                    for (size_t i = 0; i< min_length_pathlets_low_level.size(); i++){
                        auto pathlet_beginning = min_length_pathlets_low_level[i].pathlet.first;
                        auto pathlet_end = min_length_pathlets_low_level[i].pathlet.second;
                        if(!(traj_set_beginnings.get_range().contains(pathlet_beginning + trajectory_offset) && traj_set_ends.get_range().contains(pathlet_end + trajectory_offset))){
                            std::erase_if(min_length_pathlets, [&](const PathletNode<space>& pn){return min_length_pathlets_low_level[i].is_contained_by(pn);});
                        }
                    }
                    for (const auto alive_pn : min_length_pathlets){
                        surviving_pathlets.insert({trajectory_offset + alive_pn.pathlet.first, trajectory_offset + alive_pn.pathlet.second});
                    }
                    visited_pathlet_mothers.add(pathlet_mother);
                }


            }
            return surviving_pathlets.size();
        }
        static long unsigned int encode_pathlet(std::pair<index_t, index_t> pathlet, int exp){

            long unsigned int encoding = pathlet.first * POWERS_OF_TEN[exp] + pathlet.second;
            return encoding;

        }
        static std::string pathlet_key_to_string(std::pair<index_t,index_t> pathlet){
            return std::format(" ({},{}) ", pathlet.first, pathlet.second);
        }
    std::mt19937 mt;
    std::map<index_t, std::pair<point_t, std::vector<index_t>>> pathlet_beginnings;
    std::map<index_t, std::pair<point_t, std::vector<index_t>>> pathlet_ends;
    std::vector<id_t> sampled_trajs_ids;
    trajectory_t the_trajectory;
    trajectory_t the_pathlets;
    std::vector<LightWeightBinaryPathletTree<space>> the_pathlet_trees; 
    distance_t distance_threshold;
    bool thorough;
    float epsilon;
    float delta;
    int min_length;
    bool performed_sampling;
    int seed;
    float grid_side_factor;
    bool inspect;
    int vc_dim_ub = std::numeric_limits<int>::max();
};

template<metric_space m_space>
class frequent_subtrajectory_algo{
    public:
        using space = m_space;
        using trajectory_t = trajectory_collection<space>;
        using subtrajectory_t = trajectory_t::subtrajectory_t;
        using index_t = trajectory_t::index_t;
        using free_space_graph_t = free_space_graph_free_axis<space>;
        
        using point_t = space::point_t;
        using distance_function_t = space::distance_function_t;
        using distance_t = distance_function_t::distance_t;
        using binary_pathlet_tree_t = BinaryPathletTree<space>;
        using range_search_t = kd_tree_range_search<space>;
        using id_t = trajectory_t::id_t;
        struct frequent_pathlet{

            std::pair<index_t,index_t> extremes;
            id_t pathlet_mother;
            float frequency;
            float efficacy=-1.0;
            roaring::Roaring supporting_trajectories;
            friend inline bool operator<(const frequent_pathlet& lhs, const frequent_pathlet& rhs){
                
                int lhs_length = lhs.extremes.second -lhs.extremes.first +1;
                int rhs_length = rhs.extremes.second -rhs.extremes.first +1;

                if(lhs_length > rhs_length){

                    return true;

                }
                if (lhs_length < rhs_length){

                    return false;

                }
                if (lhs.pathlet_mother < rhs.pathlet_mother){
                    return true;
                }
                if (lhs.pathlet_mother> rhs.pathlet_mother){

                    return false;
                }
                return lhs.extremes.first < rhs.extremes.first;


            }
        };
        
    public:
        std::vector<frequent_pathlet> freq_pathlets;
        
        frequent_subtrajectory_algo(trajectory_t& sampled_traj, range_search_t& search, std::string dataset_file, float frequency_threshold, distance_t distance_thresh, freq_subtrajectory_algo_output_config configs) : search(search) {
            this->sample = sampled_traj; //I keep the original sample, I will build the simplification later in the constructor
            this->dataset_location = dataset_file;
            this->output_config = configs;
            //std::cout << sampled_traj.num_trajectories_not_consecutive()<< std::endl;
            //std::cout<<"Frequency threshold is "<< frequency_threshold << std::endl;
            this->integer_frequency_threshold = ceil(frequency_threshold *((int)sampled_traj.num_trajectories_not_consecutive()));
            //std::cout << "THE INTEGER FREQ THRESHOLD IS "<< this->integer_frequency_threshold<<std::endl;
            this-> last_parsed_trajectory = -1;
            this-> distance_threshold = distance_thresh;
            for (int i=0; i<25; i++){
                this->POWERS_OF_TWO.push_back(int(pow(2,i)));
            }
        }
        //Computes all frequent pathlets and saves in in this->freq_pathlets
        /*
        
        void compute_frequent_pathlets_with_trajectory_slicing_and_limited_fsg_construtction(){
            std::ifstream input_stream(this->dataset_location); //input stream that reads trajectories upon which we build the pathlets

            int sample_size = this->sample.num_trajectories_not_consecutive();
            int chunk_size = 50;//int(this-> sample.num_trajectories_not_consecutive()/50);
            
            
            while(!input_stream.eof()){
                //One pathlet tree at a time
                trajectory_t pathlet_mother = this->read_next_transaction_from_file(input_stream);

                BinaryPathletTree pathlet_tree(pathlet_mother, pathlet_mother.get_id_at(0),floor(log2(pathlet_mother.total_size())) + 1,output_config.min_length);
                //std::cout <<"The transaction has id "<< pathlet_mother.get_id_at(0) <<std::endl;
                bool no_frequent_for_this_tree = false;
                int num_visited_trajectories = 0;
                id_t next_first_id_of_chunk= this->sample.get_id_at(0);
                subtrajectory_t chunk;
                chunk.second = 0;
                //std::cout << "Testing frequency for id "<< pathlet_mother.get_id_at(0)<< std::endl;

                //Test the pathlet against the trajectories in a "chunk" and then update the chunk
                while (chunk.second < sample.total_size()-1){
                    num_visited_trajectories +=chunk_size;
                    chunk = extract_chunked_slice(this->sample,next_first_id_of_chunk, chunk_size);
                    
                    
                    //slice.first = sample.get_first_point_in_trajectory(next_id);
                    //slice.second =  sample.get_first_point_in_trajectory(next_id) + sample.get_trajectory_size(next_id);
                    
                    free_space_graph_flexible_t fsg(0);

                    //POPULATE FSG
                    //std::cout << fsg.to_string(sample,chunk) <<std::endl; 
                    this->pruned_fsg_population(fsg, chunk, pathlet_mother, pathlet_tree); //Should receive pathlet tree
                    std::cout << fsg.to_string(sample,chunk) <<std::endl; 
                    this->query_and_update_counts_for_all_pathlets<free_space_graph_flexible_t>(this->sample,fsg, chunk, pathlet_tree);
                                      
                    
                    //ASSERT AT LEAST ONE FREQUENT EXISTS in the pathlet tree
                    if(!some_potentially_frequent_exists(pathlet_tree, num_visited_trajectories, sample_size)){
                        no_frequent_for_this_tree = true;
                        break;
                    }
                    
                    
                    if(chunk.second < sample.total_size()-1){

                        id_t second_id = sample.get_id_at(chunk.second);

                        next_first_id_of_chunk = sample.get_id_at(sample.get_first_point_in_trajectory(second_id) + sample.get_trajectory_size(second_id));

                    }
                }

                //COLLECT THE FREQUENT ONES 
                if (!no_frequent_for_this_tree){
                //std::cout<< "I have found some frequent"<<std::endl;
                this->collect_frequent_pathlets_from_single_tree(pathlet_tree);
                }
            }

        }
        */

        void compute_frequent_pathlets_with_trajectory_slicing(int suggested_chunk_size = -1, bool detailed_outcome = false){
        
            std::ifstream input_stream(this->dataset_location); //input stream that reads trajectories upon which we build the pathlets

            int sample_size = this->sample.num_trajectories_not_consecutive();
            int num_skipped_pathlets = 0;
            int chunk_size = 50;//int(this-> sample.num_trajectories_not_consecutive()/50);
            //If i have an input chunk size 
            if (suggested_chunk_size >0){
                chunk_size = suggested_chunk_size;
            }
            
            while(!input_stream.eof()){
                //One pathlet tree at a time
                trajectory_t pathlet_mother = this->read_next_transaction_from_file(input_stream);
            
                BinaryPathletTree pathlet_tree(pathlet_mother, pathlet_mother.get_id_at(0),floor(log2(pathlet_mother.total_size())) + 1,1);
                //std::cout <<"The transaction has id "<< pathlet_mother.get_id_at(0) <<std::endl;
                if (forbidden_pathlet_mothers.contains(pathlet_mother.get_id_at(0))){

                    continue; // I don't bother testing the pathlets if i know they are not frequent  

                }
                
                bool no_frequent_for_this_tree = false;
                int num_visited_trajectories = 0;
                id_t next_first_id_of_chunk= this->sample.get_id_at(0);
                subtrajectory_t chunk;
                chunk.second = 0;
                //std::cout << "Testing frequency for id "<< pathlet_mother.get_id_at(0)<< std::endl;

                //Test the pathlet against the trajectories in a "chunk" and then update the chunk
                while (chunk.second < sample.total_size()-1){
                    num_visited_trajectories +=chunk_size;
                    chunk = extract_chunked_slice(this->sample,next_first_id_of_chunk, chunk_size);
                    /*
                    
                    slice.first = sample.get_first_point_in_trajectory(next_id);
                    slice.second =  sample.get_first_point_in_trajectory(next_id) + sample.get_trajectory_size(next_id);
                    */
                    free_space_graph_t fsg(0);

                    //POPULATE FSG
                    //std::cout << fsg.to_string(sample,chunk) <<std::endl; 
                    this->populate_all_columns_with_labels_for_single_slice(fsg, this->sample, chunk, pathlet_mother, pathlet_tree); //Should receive pathlet tree
                    //std::cout << fsg.to_string(sample,chunk) <<std::endl; 
                    this->query_and_update_counts_for_all_pathlets<free_space_graph_t>(this->sample, fsg, chunk, pathlet_tree);
                                      
                    
                    //ASSERT AT LEAST ONE FREQUENT EXISTS in the pathlet tree
                    if(!some_potentially_frequent_exists(pathlet_tree, num_visited_trajectories, sample_size)){
                        no_frequent_for_this_tree = true;
                        break;
                    }
                    
                    if(suggested_chunk_size>0 && num_visited_trajectories<=chunk_size){
                        //Run the check if someone matched against the """net"""
                        if(!some_pathlet_appears_from_pathlet_mother(pathlet_tree)){
                            no_frequent_for_this_tree = true;
                            break;
                        }
                        if(detailed_outcome){

                            num_skipped_pathlets += num_spared_pathlets(pathlet_tree);
                        }

                    }
                    
                    
                    if(chunk.second < sample.total_size()-1){

                        id_t second_id = sample.get_id_at(chunk.second);

                        next_first_id_of_chunk = sample.get_id_at(sample.get_first_point_in_trajectory(second_id) + sample.get_trajectory_size(second_id));

                    }
                }

                //COLLECT THE FREQUENT ONES 
                if (!no_frequent_for_this_tree){
                //std::cout<< "I have found some frequent"<<std::endl;
                this->collect_frequent_pathlets_from_single_tree(pathlet_tree);
                }

            }

            if(suggested_chunk_size >0 && detailed_outcome){

                std::cout << std::format("SKIPPED PATHLETS : {}\n", num_skipped_pathlets);

            }
        }

        void compute_pathlet_filter_with_theta_net(std::string netfilename, bool detailed_outcome){

            trajectory_t net_sample = read_trajectory_from_file<space>(netfilename);

            std::ifstream input_stream(this->dataset_location); //input stream that reads trajectories upon which we build the pathlets

            int sample_size = net_sample.num_trajectories_not_consecutive();
            int chunk_size = 50;//int(this-> sample.num_trajectories_not_consecutive()/50);
            
            
            while(!input_stream.eof()){
                //One pathlet tree at a time
                trajectory_t pathlet_mother = this->read_next_transaction_from_file(input_stream);

                BinaryPathletTree pathlet_tree(pathlet_mother, pathlet_mother.get_id_at(0),floor(log2(pathlet_mother.total_size())) + 1,1);
                
                int virtual_count = 1;
                
                //std::cout <<"The transaction has id "<< pathlet_mother.get_id_at(0) <<std::endl;
                bool no_frequent_for_this_tree = false;
                int num_visited_trajectories = 0;
                id_t next_first_id_of_chunk= net_sample.get_id_at(0);
                subtrajectory_t chunk;
                chunk.second = 0;
                //std::cout << "Testing frequency for id "<< pathlet_mother.get_id_at(0)<< std::endl;

                //Test the pathlet against the trajectories in a "chunk" and then update the chunk
                while (chunk.second < net_sample.total_size()-1){
                    num_visited_trajectories +=chunk_size;
                    chunk = extract_chunked_slice(net_sample, next_first_id_of_chunk, chunk_size);
                    /*
                    
                    slice.first = sample.get_first_point_in_trajectory(next_id);
                    slice.second =  sample.get_first_point_in_trajectory(next_id) + sample.get_trajectory_size(next_id);
                    */
                    free_space_graph_t fsg(0);

                    //POPULATE FSG
                    //std::cout << fsg.to_string(sample,chunk) <<std::endl; 
                    this->populate_all_columns_with_labels_for_single_slice(fsg, net_sample, chunk, pathlet_mother, pathlet_tree); //Should receive pathlet tree
                    //std::cout << fsg.to_string(sample,chunk) <<std::endl; 
                    this->query_and_update_counts_for_all_pathlets<free_space_graph_t>(net_sample,fsg, chunk, pathlet_tree);
                    
                    if(chunk.second < net_sample.total_size()-1){

                        id_t second_id = net_sample.get_id_at(chunk.second);

                        next_first_id_of_chunk = net_sample.get_id_at(net_sample.get_first_point_in_trajectory(second_id) + net_sample.get_trajectory_size(second_id));

                    }
                }

                //COLLECT THE FREQUENT ONES 
                if (!some_pathlet_appears_from_pathlet_mother(pathlet_tree)){
                //std::cout<< std::format("I am adding {} to forbidden trajectories \n", pathlet_tree.getTrajectoryId());
                    this->forbidden_pathlet_mothers.add(pathlet_tree.getTrajectoryId());
                    if(detailed_outcome){
                        this->num_spared_pathlets += get_spared_pathlets(pathlet_tree);
                    }
                }
            }

            //for (auto item : forbidden_pathlet_mothers){
            //    std::cout << item <<std::endl;
            //}
            if (detailed_outcome){
                std::cout << std::format("PRUNED PATHLETS : {}\n", this->num_spared_pathlets);
            }
        }
        //FLUSH THE FREQUENT PATHLETS TO A FILE
        //Format for a line: start_idx end_idx pathlet_mother_id frequency
        void dump_collected_pathlets_to_file(std::string outputfilename){

            std::ofstream outfile(outputfilename);
            // Dump all freq pathlets in the data structure to a file  
            for(frequent_pathlet& fp : this->freq_pathlets){

                outfile << fp.extremes.first <<" "<< fp.extremes.second<< " ";
                outfile << fp.pathlet_mother << " "<< fp.frequency<< std::endl;

            }

            outfile.close();
        }
    
    
    private:
        //Returns true if at least one pathlet in pathlet_tree can still be frequent
        bool some_potentially_frequent_exists(binary_pathlet_tree_t& pathlet_tree, int num_visited_trajectories, int sample_size){

            //Interrupt if no pathlet has hope to be frequent in the remaining steps
            int missing_count = sample_size - num_visited_trajectories;
            //std::cout <<"SAMPLE SIZE IS "<< sample_size<< "WHILE NUM VISITED TRAJECTORIES IS "<< num_visited_trajectories<< " HENCE MISSING COUNT "<<missing_count<<std::endl;
            int d = pathlet_tree.getDepth(); 
            bool found_potential_frequent= false;
            int num_sampled_trajs = this->sample.num_trajectories_not_consecutive();

            for (int level = d; d>=0; d--){

                int level_beginning = int(POWERS_OF_TWO[d])-1;

                //Traverse the tree from left to right
                for (int offset = 0; offset <=level_beginning; offset++){

                    int position = level_beginning + offset;
                    auto& pn = pathlet_tree.getNodeAt(position);
                    if(pn.isNULL){
                        
                        continue;

                    }
                    int count = pn.frequency; //SF IS HERE

                                     
                    if(count >= this->integer_frequency_threshold || (count + missing_count >= this->integer_frequency_threshold))  {

                        found_potential_frequent =  true;
                    }
                    else{
                        pathlet_tree.setInfrequent(position);
                    }
                }

            }

            return found_potential_frequent; //No hope or no frequent 
        }
        //Extracts indices of the sample trajectory_t corresponding to chu nk_size distinct sampled trajectories
        subtrajectory_t extract_chunked_slice(trajectory_t& to_be_chunked, id_t last_visited_id, id_t chunk_size){

            subtrajectory_t chunk;
            index_t starting_point_for_new_chunk = to_be_chunked.get_first_point_in_trajectory(last_visited_id);
            chunk.first = starting_point_for_new_chunk;
            int chunked_trajs =1;

            id_t last_id = to_be_chunked.get_id_at(starting_point_for_new_chunk);   
            while (chunked_trajs< chunk_size && last_id!=to_be_chunked.get_id_at(to_be_chunked.total_size()-1)){

                last_id = to_be_chunked.get_id_at(to_be_chunked.get_first_point_in_trajectory(last_id) + to_be_chunked.get_trajectory_size(last_id));
                chunked_trajs +=1;
            }


            chunk.second = to_be_chunked.get_first_point_in_trajectory(last_id) + to_be_chunked.get_trajectory_size(last_id) -1;
            
            //std::cout << "Computed a chunk that contains trajs from "<< to_be_chunked.get_id_at(chunk.first)<< " and "<<to_be_chunked.get_id_at(chunk.second)<<std::endl;
            return chunk;
        }
        //Reads next pathlet_mother from file 
        trajectory_t read_next_transaction_from_file(std::ifstream& file){

            double x,y;
            id_t id;
            trajectory_t pathlet_mother;

            file >> x >> y >> id;

            pathlet_mother.push_back({x,y}, id);
            this-> last_parsed_trajectory = id;
            
            std::streampos sp = file.tellg();
            while(file >> x>> y >> id){

                if(id != last_parsed_trajectory){
                    //std::cout<< "i have just found the beginning of trajectory "<< last_parsed_trajectory<< std::endl;
                    //ripristina pointer alla riga precedente
                    file.seekg(sp);
                    break;

                }
                pathlet_mother.push_back({x,y}, id);
                //update pointer
                sp = file.tellg();
            }
            
            //std::cout << "trajectory at the end is "<< pathlet_mother.get_id_at(1)<< std::endl;
            //std::cout << "actual size "<< pathlet_mother.get_actual_size()<<std::endl;
            //std::cout << "num trajectories "<< pathlet_mother.num_trajectories()<< std::endl;
            
            return pathlet_mother;

        }
        //After a pathlet_tree has been queried against the sample, it extracts the frequent pathlets and saves them into this->freq_pathlets
        void collect_frequent_pathlets_from_single_tree(binary_pathlet_tree_t& pathlet_tree){
            if (output_config.maximal) {
                this->collect_maximal_frequent_pathlets_from_single_tree(pathlet_tree);
                return;
            }

            int d = pathlet_tree.getDepth(); 

            int num_sampled_trajs = this->sample.num_trajectories_not_consecutive();

            for (int level = d; d>=0; d--){

                int level_beginning = int(POWERS_OF_TWO[d])-1;

                //Traverse the tree from left to right
                for (int offset = 0; offset <=level_beginning; offset++){

                    int position = level_beginning + offset;
                    auto& pn = pathlet_tree.getNodeAt(position);
                    //std::cout << "I am visitingq querying pathlet "<< pn.getPathlet().first <<" "<< pn.getPathlet().second<< std::endl;
                    if(pn.isNULL ||!(pn.isFrequent())){
                        
                        continue;

                    }

                    int count = pn.frequency; 
                    if(count >= this->integer_frequency_threshold)  {
                        
                        pathlet_tree.setEstimatedFrequency(position, ((double)count / num_sampled_trajs));
                        frequent_pathlet just_found;
                        just_found.extremes = pn.getPathlet();
                        just_found.pathlet_mother = (pathlet_tree.getTrajectoryId());
                        just_found.frequency = ((double) count / num_sampled_trajs);
                        just_found.supporting_trajectories = pn.getSupportingTrajectories();
                        if ((just_found.extremes.second-just_found.extremes.first +1 )>=output_config.min_length){
                            
                            freq_pathlets.push_back(just_found);
                        }
                    }    

                }
            }
        }
        //Queries the pathlets in pathlet_tree against the free space graph fsg, updating their counts when a match is found
        template <typename fsg_type>
        void query_and_update_counts_for_all_pathlets(trajectory_t& unsliced_traj,fsg_type& fsg, subtrajectory_t& slice, binary_pathlet_tree_t& pathlet_tree){
            
            int d = pathlet_tree.getDepth(); 

            for (int level = d; d>=0; d--){

                int level_beginning = int(POWERS_OF_TWO[d])-1;

                //Traverse the tree from left to right
                for (int offset = 0; offset <=level_beginning; offset++){

                    int position = level_beginning + offset;
                    auto& pn = pathlet_tree.getNodeAt(position);
                    if(pn.isNULL ||!(pn.isFrequent())){
                        if (!pn.isFrequent()){
                        //std::cout<< "SEARCH PRUNING"<<std::endl;
                        }
                        continue;

                    }
                    
                    std::set<id_t> matching_ids = fsg.query_one_pathlet_over_the_sample_with_labels_by_slice(unsliced_traj, pn.getPathlet(),1, slice); //SF IS HERE
                    
                    int count = matching_ids.size();
                    if (this->output_config.keep_matching_ids){

                        for (id_t idd : matching_ids){
                            //if(pathlet_tree.getTrajectoryId() ==951 ||pathlet_tree.getTrajectoryId() ==1854 || pathlet_tree.getTrajectoryId() ==2237){
                            //    if (idd ==951 ||idd ==1854 || idd ==2237){
                                    //std::cout << "Pathlet "<< pn.getPathlet().first <<" "<< pn.getPathlet().second<<"of trajectory "<<pathlet_tree.getTrajectoryId()<<" has matches with tid "<< idd<< std::endl;
                            //    }
                            //}
                            pn.addId(idd);
                        }

                    }
                    
                    //std::cout << "Pathlet "<< pn.getPathlet().first <<" "<< pn.getPathlet().second<<"has matches with "<< slice.first<<" "<< slice.second<< std::endl;
                    
                    //pn.frequency +=count;
                    pathlet_tree.setEstimatedFrequency(position, pathlet_tree.getNodeAt(position).getFrequency() + count);
                    
                }


            }
            
            

        }
        //After a pathlet_tree has been queried against the sample, it extracts the MAXIMAL frequent pathlets and saves them into this->freq_pathlets
        void collect_maximal_frequent_pathlets_from_single_tree(binary_pathlet_tree_t& pathlet_tree){
            int num_sampled_trajs = this->sample.num_trajectories_not_consecutive();

            if(pathlet_tree.getNodeAt(0).frequency>=this->integer_frequency_threshold){
                PathletNode pn = pathlet_tree.getNodeAt(0);
                frequent_pathlet just_found;
                just_found.extremes = pn.getPathlet();
                just_found.pathlet_mother = (pathlet_tree.getTrajectoryId());
                just_found.frequency = pn.frequency/num_sampled_trajs;
                
                freq_pathlets.push_back(just_found);
                return;
            }
            //TODO: clean up beacuse this is horrible
            std::queue<int> nodes_to_visit;
            PathletNode pn =pathlet_tree.getNodeAt(binary_pathlet_tree_t::left_child_idx(0));
            if(!pn.isNULL && pn.getLength() > 1){
                nodes_to_visit.push(binary_pathlet_tree_t::left_child_idx(0));
            }

            PathletNode pn1=pathlet_tree.getNodeAt(binary_pathlet_tree_t::right_child_idx(0));

            if(!pn1.isNULL && pn1.getLength() > 1){

                nodes_to_visit.push(binary_pathlet_tree_t::right_child_idx(0));

            }

            while(!nodes_to_visit.empty()){

                int pathlet_idx = nodes_to_visit.front();
                nodes_to_visit.pop();

                if(pathlet_tree.getNodeAt(pathlet_idx).frequency>=this->integer_frequency_threshold){

                    PathletNode pn = pathlet_tree.getNodeAt(pathlet_idx);
                    frequent_pathlet just_found;
                    just_found.extremes = pn.getPathlet();
                    just_found.pathlet_mother = (pathlet_tree.getTrajectoryId());
                    just_found.frequency = pn.frequency/num_sampled_trajs;
                    freq_pathlets.push_back(just_found);                    

                }
                else{

                    PathletNode pn =pathlet_tree.getNodeAt(binary_pathlet_tree_t::left_child_idx(pathlet_idx));
                    if(!pn.isNULL && pn.getLength() > 1){
                        nodes_to_visit.push(binary_pathlet_tree_t::left_child_idx(pathlet_idx));
                    }

                    PathletNode pn1 =pathlet_tree.getNodeAt(binary_pathlet_tree_t::right_child_idx(pathlet_idx));

                    if(!pn1.isNULL && pn1.getLength() > 1){

                        nodes_to_visit.push(binary_pathlet_tree_t::right_child_idx(pathlet_idx));

                    }

                }

            }


        }
        //Populates a single column of the free space graph 
        void populate_column_with_labels_for_single_slice(free_space_graph_t &fsg, trajectory_t& unsliced_traj, const point_t& point, const distance_t& query_distance, index_t column_index, subtrajectory_t& slice){
            int zeroes = 0;
            float sq_dist = query_distance * query_distance;
            //std::cout << "populating at distance  "<< query_distance << std::endl; 
            int highest_index = 0;
            int iterations =0;
            //ARTIGIANALE:
            for (int i = slice.first; i<=slice.second; i++){
                iterations++;
                auto d_ij = distance_function_t{}(unsliced_traj[i], point); //FI
                if (d_ij <= sq_dist ){
                    //std::cout << "d_ij for i "<< i<<" and j "<< column_index << " is "<< d_ij << " while sq_dist is "<< sq_dist<<std::endl;
                    fsg.add_zero_with_id_respecting_labels(i, unsliced_traj);
                    zeroes++;

                } 
            }
            //std::cout << zeroes<< std::endl;
            //std::cout<< "i create columns without sf"<< std::endl;
            //assert(iterations == sample.get_actual_size());
        }
        /*
        
        roaring::Roaring populate_column_with_labels_for_single_slice(free_space_graph_flexible_t &fsg, const point_t& point, const distance_t& query_distance, index_t column_index, subtrajectory_t& slice, roaring::Roaring& potentially_matching_trajs){
            int zeroes = 0;
            float sq_dist = query_distance * query_distance;
            //std::cout << "populating at distance  "<< query_distance << std::endl; 
            int highest_index = 0;
            int iterations =0;
            
            //ARTIGIANALE:
            roaring::Roaring matching_trajectories{};
            if(fsg.is_intermediate(column_index)){
                for(auto tid: potentially_matching_trajs){
                    
                    index_t trajectory_offset = sample.get_first_point_in_trajectory(tid);
                    int i =0;
                    while (i<sample.get_trajectory_size(tid)-1){
                        auto d_ij = distance_function_t{}(sample[trajectory_offset + i], point); //FI
                        if (d_ij <= sq_dist ){
                            //std::cout << "d_ij for i "<< i<<" and j "<< column_index << " is "<< d_ij << " while sq_dist is "<< sq_dist<<std::endl;
                            bool potential_connection = fsg.add_zero_with_id_respecting_labels(i,column_index, this->sample);
                            if(potential_connection){
                                
                                matching_trajectories.add(sample.get_id_at(i));

                            }

                        } 
                        i++;
                    }

                }
            }
            else{
                //TODO: retrieve from a grid
                for (int i = slice.first; i<=slice.second; i++){
                iterations++;
                auto d_ij = distance_function_t{}(sample[i], point); //FI
                if (d_ij <= sq_dist ){
                    //std::cout << "d_ij for i "<< i<<" and j "<< column_index << " is "<< d_ij << " while sq_dist is "<< sq_dist<<std::endl;
                    fsg.add_zero_with_id_respecting_labels(i, column_index, this->sample);
                    matching_trajectories.add(sample.get_id_at(i));
                    zeroes++;

                } 
                }
            }
            //std::cout << zeroes<< std::endl;
            //std::cout<< "i create columns without sf"<< std::endl;
            //assert(iterations == sample.get_actual_size());

            //possibly, these are less than the subset i am selecting them from 
            return matching_trajectories;
        }
        
        */

        //Populates all columns of the free space graph (from left to right) for a single pathlet_tree, against a specific slice of the sample
        void populate_all_columns_with_labels_for_single_slice(free_space_graph_t& fsg, trajectory_t& unsliced_traj, subtrajectory_t& slice, trajectory_t& pathlet_mother, binary_pathlet_tree_t& pathlet_tree){

            int num_col = pathlet_mother.get_actual_size();
            
            for (int j = 0; j< num_col; j++){
                // If a point is not frequent and i already know it, skip the distance computations 
                //get the node associated to that point
                //PathletNode pn = pathlet_tree.getPointPathletNode(j);
                //if(pn.isFrequent()){
                populate_column_with_labels_for_single_slice(fsg, unsliced_traj, pathlet_mother[j], this->distance_threshold, j,slice);
                //}
                if(j < num_col -1){

                    fsg.new_column();

                }
            }

            return;


        }
        
        bool some_pathlet_appears_from_pathlet_mother(binary_pathlet_tree_t& pathlet_tree){

            int d = pathlet_tree.getDepth(); 
            //int num_sampled_trajs = this->sample.num_trajectories_not_consecutive();

            for (int level = d; d>=0; d--){

                int level_beginning = int(POWERS_OF_TWO[d])-1;

                //Traverse the tree from left to right
                for (int offset = 0; offset <=level_beginning; offset++){

                    int position = level_beginning + offset;
                    auto& pn = pathlet_tree.getNodeAt(position);
                    //std::cout << "I am visitingq querying pathlet "<< pn.getPathlet().first <<" "<< pn.getPathlet().second<< std::endl;
                    if(pn.isNULL || pn.getLength()<this->output_config.min_length){
                        
                        continue;

                    }

                    int count = pn.frequency; 
                    //std::cout << std::format("frequency is {}\n", count);
                    if ((count>=1)){
                        
                        return true;

                    }

                }
            }

            return false;
        }

        int get_spared_pathlets(binary_pathlet_tree_t& pathlet_tree){

            // Count the number of spared_pathlets
            int d = pathlet_tree.getDepth(); 

            //int num_sampled_trajs = this->sample.num_trajectories_not_consecutive();
            int spared = 0;

            for (int level = d; d>=0; d--){

                int level_beginning = int(POWERS_OF_TWO[d])-1;

                //Traverse the tree from left to right
                for (int offset = 0; offset <=level_beginning; offset++){

                    int position = level_beginning + offset;
                    auto& pn = pathlet_tree.getNodeAt(position);
                    //std::cout << "I am visitingq querying pathlet "<< pn.getPathlet().first <<" "<< pn.getPathlet().second<< std::endl;
                    if(pn.isNULL || pn.getLength()<this->output_config.min_length){
                        
                        continue;

                    }

                    spared++;
                }
            }

            return spared;
        }

        
        std::vector<int> POWERS_OF_TWO;
        trajectory_t sample;
        range_search_t& search;
        std::string dataset_location;
        id_t last_parsed_trajectory;
        distance_t distance_threshold;
        int integer_frequency_threshold;
        freq_subtrajectory_algo_output_config output_config;
        roaring::Roaring forbidden_pathlet_mothers;
        int  num_spared_pathlets =0;
};
};

