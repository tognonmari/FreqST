#pragma once

#include <cmath>
#include <iostream>
#include <limits>
#include <list>
#include <vector>
#include <random>
#include <algorithm>
#include <random>

#include "free_space_graph_incremental.h"
#include "free_space_graph_free_axis.h"
#include "kdtree_range_search.h"
#include "metric_space.h"
#include "io.h"
#include "trajectory.h"
#include "canonical_pathlets.h"
#include "curve_simplification.h"
#include "freq_st_algo.h"
namespace frechet{

template<metric_space m_space>
class frequent_subtrajectory_algo_simplified{
    public:
        using space = m_space;
        using trajectory_t = trajectory_collection<space>;
        using index_t = trajectory_t::index_t;
        using free_space_graph_t = free_space_graph_free_axis<space>;
        using point_t = space::point_t;
        using distance_function_t = space::distance_function_t;
        using distance_t = distance_function_t::distance_t;
        using binary_pathlet_tree_t = BinaryPathletTree<space>;
        using range_search_t = kd_tree_range_search<space>;
        using subtrajectory_t = trajectory_t::subtrajectory_t;
        using curve_simplification_cluster_summary_t = frechet::internal::curve_simplification<space>::cluster_summary_t;
        using frequent_pathlet = typename frequent_subtrajectory_algo<space>::frequent_pathlet;
        /*
        struct frequent_pathlet{

            std::pair<index_t,index_t> extremes;
            id_t pathlet_mother;
            float frequency;
            
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
    */
    public:
        std::vector<frequent_pathlet> freq_pathlets;

        
        frequent_subtrajectory_algo_simplified(trajectory_t& sampled_traj, range_search_t& search, std::string dataset_file, float frequency_threshold, distance_t& distance_thresh, distance_t& simplification_factor) : sample(sampled_traj), distance_threshold(distance_thresh),curve_simplification_factor(simplification_factor), simplification{sample, distance_threshold*distance_threshold,curve_simplification_factor}, search(search) {
            this->sample = sampled_traj; 
            this->dataset_location = dataset_file;
            std::cout << sample.num_trajectories_not_consecutive()<< std::endl;
            std::cout<<"Frequency threshold is "<< frequency_threshold << std::endl;
            this->integer_frequency_threshold = ceil(frequency_threshold *((int)sample.num_trajectories_not_consecutive()));
            std::cout << "THE INTEGER FREQ THRESHOLD IS "<< this->integer_frequency_threshold<<std::endl;
            this-> last_parsed_trajectory = -1;
            this-> distance_threshold = distance_thresh;
            this->curve_simplification_factor = curve_simplification_factor;
            //frechet::internal::curve_simplification<space> temporary_simplification(sampled_traj, distance_threshold * distance_threshold, curve_simplification_factor);
            std::cout << "Simplified traj has size : "<< this->simplification.trajectory().total_size()<< std::endl;

        }


        void compute_all_frequent_pathlets(){

            //open full dataset file 

            std::ifstream input_stream(this->dataset_location);

            //std::cout <<"Starting reading the transactions."<<std::endl;

            while(!input_stream.eof()){

                trajectory_t pathlet_mother = this->read_next_transaction_from_file(input_stream);
                frechet::internal::curve_simplification<space> pathlet_mother_simplification(pathlet_mother, this->distance_threshold * this->distance_threshold, this->curve_simplification_factor);
                //std::cout <<"Parsed a transaction."<<std::endl;
                //std::cout<<" The transaction has ID "<<pathlet_mother.get_id_at(pathlet_mother.get_actual_size()-1)<<std::endl;
                //std::cout <<" I have this many points : "<< pathlet_mother.get_actual_size()<<std::endl;
                trajectory_t simplified_mother_trajectory = pathlet_mother_simplification.trajectory();
                BinaryPathletTree pathlet_tree(simplified_mother_trajectory, pathlet_mother.get_id_at(0),floor(log2(pathlet_mother_simplification.trajectory().total_size())) + 1,1);
                
                //INITIALIZE FREE SPACE DIAGRAM 
                free_space_graph_t fsg(0);
                
                //POPULATE FREE SPACE DIAGRAM : THE SAMPLE IS ALONG THE Y-AXIS, THE TRANSACTION WHICH BAERS THE PATHLETS FOR THE CURRENT TREE IS ALONG THE X-AXIS

                this->populate_all_columns_with_labels(fsg, pathlet_mother_simplification.trajectory());

                //std::cout <<"Populated the columns."<<std::endl;
                
                //COLLECT THE FREQUENT PATHLETS GENERATED BY THE CURRENT PATHLET MOTHER

                this->collect_all_frequent_pathlets(fsg, pathlet_tree);

            }



        }


        void compute_maximal_frequent_pathlets(){

            //open full dataset file 
            std::ifstream input_stream(this->dataset_location);
            //std::cout <<"Starting reading the transactions."<<std::endl;
            while(!input_stream.eof()){

                trajectory_t pathlet_mother = this->read_next_transaction_from_file(input_stream);
                frechet::internal::curve_simplification<space> pathlet_mother_simplification(pathlet_mother, this->distance_threshold * this->distance_threshold, this->curve_simplification_factor);
                //std::cout <<"Parsed a transaction."<<std::endl;
                //std::cout<<" The transaction has ID "<<pathlet_mother.get_id_at(pathlet_mother.get_actual_size()-1)<<std::endl;
                //std::cout <<" I have this many points : "<< pathlet_mother.get_actual_size()<<std::endl;
                trajectory_t simplified_mother_trajectory = pathlet_mother_simplification.trajectory();
                BinaryPathletTree pathlet_tree(simplified_mother_trajectory, pathlet_mother.get_id_at(0),floor(log2(pathlet_mother_simplification.trajectory().total_size())) + 1,1);
                
                //INITIALIZE FREE SPACE DIAGRAM 
                free_space_graph_t fsg(0);
                
                //POPULATE FREE SPACE DIAGRAM : THE SAMPLE IS ALONG THE Y-AXIS, THE TRANSACTION WHICH BAERS THE PATHLETS FOR THE CURRENT TREE IS ALONG THE X-AXIS

                this->populate_all_columns_with_labels(fsg, pathlet_mother_simplification.trajectory());

                //std::cout <<"Populated the columns."<<std::endl;
                
                //COLLECT THE FREQUENT PATHLETS GENERATED BY THE CURRENT PATHLET MOTHER

                this->collect_maximal_frequent_pathlets(fsg, pathlet_tree);

            }


        }

        //FLUSH THE FREQUENT PATHLETS TO A FILE

        void dump_collected_pathlets_to_file(std::string outputfilename){
            //TODO: capire se voglio dumparli simplified or unsimplified.
            std::ofstream outfile(outputfilename);
            // Dump all freq pathlets in the data structure to a file  
            for(frequent_pathlet& fp : this->freq_pathlets){

                outfile << fp.extremes.first <<" "<< fp.extremes.second<< " ";
                outfile << fp.pathlet_mother << " "<< fp.frequency<< std::endl;

            }

            outfile.close();
        }
        
        void unsimplify_collected_pathlets(trajectory_t& original_traj){
            frechet::internal::curve_simplification<space> cs(original_traj, distance_threshold * distance_threshold, curve_simplification_factor);
            for (auto& p: this->freq_pathlets){
                
                //from frequent pathlets to indexes in trajectories 
                subtrajectory_t offsets = p.extremes;
                index_t initial_point = cs.trajectory().get_first_point_in_trajectory(p.pathlet_mother);
                //assert(cs.trajectory().get_first_point_in_trajectory(p.pathlet_mother)==simplification.trajectory().get_first_point_in_trajectory(p.pathlet_mother));
                std::cout << "Initial point for trajectory "<< p.pathlet_mother << "is" << simplification.trajectory().get_first_point_in_trajectory(2)<< std::endl;
                subtrajectory_t st{initial_point+offsets.first, initial_point+ offsets.second};
                //Now unsimplify st
                std::optional<curve_simplification_cluster_summary_t> c = curve_simplification_cluster_summary_t{0,0,0.0,st.first,st.second};
                std::cout <<"Unsimplifying Pathlet "<<p.extremes.first << " "<< p.extremes.second << " with mother "<< p.pathlet_mother<< std::endl;
                std::cout << "This corresponds to subtrajectory "<< st.first << " "<< st.second << std::endl;
                std::optional<curve_simplification_cluster_summary_t> temp = cs.unsimplify(c);
                p.extremes.first = temp.value().left_column - original_traj.get_first_point_in_trajectory(p.pathlet_mother);
                p.extremes.second = temp.value().right_column - original_traj.get_first_point_in_trajectory(p.pathlet_mother);

            }


        }

    private:

        // The method queries the pathlets present in the BST and saves the frequent ones in the freq_pathlets data structure
        void collect_maximal_frequent_pathlets(free_space_graph_t& fsg, binary_pathlet_tree_t& pathlet_tree){

            int d = pathlet_tree.getDepth(); //Last filled level
            //std::cout<< "tree has depth "<< d <<std::endl; // assertion for my toy dataset
            int num_sampled_trajs = this->sample.num_trajectories();

            //BOTTOM-UP VISIT OF THE TREE (starts from the deepest level with at least one valid pathlet, can access levels directly)

            for (int level = d; d>=0; d--){

                //Compute offset for direct access (no traversal each time)

                int level_beginning = int(pow(2,d))-1;

                //Traverse the level from left to right
                for (int offset = 0; offset <=level_beginning; offset++){

                    int position = level_beginning + offset;

                    PathletNode pn = pathlet_tree.getNodeAt(position);


                    //if the current node in the binary pathlet tree is either improper or infrequent, skip it.

                    if(pn.isNULL ||!(pn.isFrequent())){
                        if (!pn.isFrequent()){
                        //std::cout<< "SEARCH PRUNING"<<std::endl;
                        }
                        continue;

                    }
                    
                    //FIND THE NUMBER OF TRAJECTORIES IN THE SAMPLE (i.e. along the Y-AXIS), WHICH MATCH AGAINST THE PATHLET (count)
                    int count = fsg.query_one_pathlet_over_the_sample_with_labels(sample, pn.getPathlet()); 
                    
                    //IF THE PATHLET IS FREQUENT, SAVE ITS FRERQUENCY. OTHERWISE MARK IT (AND ALL OF ITS ANCESTORS) AS INFREQUENT
                    if(count < this->integer_frequency_threshold){
                        
                        pathlet_tree.setInfrequent(position);

                    }
                    else{

                        pathlet_tree.setEstimatedFrequency(position, ((float) count / num_sampled_trajs));

                    }

                }


            }

            if(pathlet_tree.getNodeAt(0).isFrequent()){
                PathletNode pn = pathlet_tree.getNodeAt(0);
                frequent_pathlet just_found;
                just_found.extremes = pn.getPathlet();
                just_found.pathlet_mother = (pathlet_tree.getTrajectoryId());
                just_found.frequency = pn.frequency;
                
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

                if(pathlet_tree.getNodeAt(pathlet_idx).isFrequent()){

                    PathletNode pn = pathlet_tree.getNodeAt(pathlet_idx);
                    frequent_pathlet just_found;
                    just_found.extremes = pn.getPathlet();
                    just_found.pathlet_mother = (pathlet_tree.getTrajectoryId());
                    just_found.frequency = pn.frequency;
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
            /*
            
            if (curve_simplification_factor>0){

                //Unsimplify the pathlets
                for (auto& entry : freq_pathlets){
                    entry.extremes.first = simplification.get_original_index(entry.extremes.first);
                    entry.extremes.second = simplification.get_original_index(entry.extremes.second);
                }
            }
            
            */
            
        }
        

        void collect_all_frequent_pathlets(free_space_graph_t& fsg, binary_pathlet_tree_t& pathlet_tree){

            int d = pathlet_tree.getDepth(); 

            int num_sampled_trajs = this->sample.num_trajectories();

            for (int level = d; d>=0; d--){

                int level_beginning = int(pow(2,d))-1;

                //Traverse the tree from left to right
                for (int offset = 0; offset <=level_beginning; offset++){

                    int position = level_beginning + offset;
                    PathletNode pn = pathlet_tree.getNodeAt(position);
                    //std::cout << "I am visitingq querying pathlet "<< pn.getPathlet().first <<" "<< pn.getPathlet().second<< std::endl;
                    if(pn.isNULL ||!(pn.isFrequent())){
                        if (!pn.isFrequent()){
                        //std::cout<< "SEARCH PRUNING"<<std::endl;
                        }
                        continue;

                    }

                    int count = fsg.query_one_pathlet_over_the_sample_with_labels(sample, pn.getPathlet()); //SF IS HERE
                    
                    if(count < this->integer_frequency_threshold){
                        
                        pathlet_tree.setInfrequent(position);

                    }
                    else{
                        
                        pathlet_tree.setEstimatedFrequency(position, ((float)count / num_sampled_trajs));
                        assert(pathlet_tree.getNodeAt(position).frequency - ((float)count / num_sampled_trajs) < 0.0001 );
                        frequent_pathlet just_found;
                        just_found.extremes = pn.getPathlet();
                        just_found.pathlet_mother = (pathlet_tree.getTrajectoryId());
                        just_found.frequency = ((float) count / num_sampled_trajs);
                        freq_pathlets.push_back(just_found);
                    }

                }


            }


        }
        void populate_all_columns_with_labels(free_space_graph_t& fsg, const trajectory_t& pathlet_mother){
            int num_col = pathlet_mother.get_actual_size();

            for (int j = 0; j< num_col; j++){

                populate_column_with_labels(fsg, pathlet_mother[j], this->distance_threshold, j);
                if(j < num_col -1){

                    fsg.new_column();

                }
            }

            return;

        }
        //MATERIALIZES THE NON-ZERO NODES IN THE FSG + THE EDGES TO TRAVERSE IT 
        void populate_all_columns(free_space_graph_t& fsg, const trajectory_t& pathlet_mother){

            int num_col = pathlet_mother.get_actual_size();

            for (int j = 0; j< num_col; j++){

                populate_column_with_labels(fsg, pathlet_mother[j], this->distance_threshold, j);
                if(j < num_col -1){

                    fsg.new_column();

                }
            }

            return;
        }
        void populate_column_with_labels(free_space_graph_t &fsg, const point_t& point, const distance_t& query_distance, index_t column_index ){

            
            int zeroes = 0;
            float sq_dist = query_distance * query_distance;
            int highest_index = 0;
            int iterations =0;
            //ARTIGIANALE:
            for (int i = 0; i< sample.get_actual_size(); i++){
                iterations++;
                auto d_ij = distance_function_t{}(sample[i], point); //FI
                if (d_ij <= sq_dist ){

                    fsg.add_zero_with_id_respecting_labels(i, this->sample);
                    zeroes++;

                } 
            }
            //std::cout << zeroes<< std::endl;
            //std::cout<< "i create columns without sf"<< std::endl;
            //assert(iterations == sample.get_actual_size());

        }
        void populate_column_with_labels_and_range_search(free_space_graph_t &fsg, const point_t& point, const distance_t& query_distance, index_t column_index ){

            
            int zeroes = 0;
            //float sq_dist = query_distance * query_distance;
            int highest_index = 0;
            int iterations =0;
            //ARTIGIANALE:
            
            for (const auto idx: search.search_by_point(point, query_distance*query_distance)) {

                fsg.add_zero_with_id_respecting_labels(idx, this->sample);
                zeroes++;
            }
            /*
            for (int i = 0; i< sample.get_actual_size(); i++){
                iterations++;
                auto d_ij = distance_function_t{}(sample[i], point); //FI
                if (d_ij <= sq_dist ){

                    fsg.add_zero_with_id_respecting_labels(i, this->sample);
                    zeroes++;

                } 
            }
            
            */
            
            //std::cout<< "i create columns without sf"<< std::endl;
            //assert(iterations == sample.get_actual_size());

        }
        // MATERIALIZES A SINGLE COLUMN OF THE FSG, i.e. the nodes for which 
        void populate_column(free_space_graph_t &fsg, const point_t& point, const distance_t& query_distance, index_t column_index ){

            
            int zeroes = 0;
            float sq_dist = query_distance * query_distance;
            int highest_index = 0;
            int iterations =0;
            //ARTIGIANALE:
            for (int i = 0; i< sample.get_actual_size(); i++){
                iterations++;
                auto d_ij = distance_function_t{}(sample[i], point);
                if (d_ij <= sq_dist ){

                    fsg.add_zero(i);
                    zeroes++;

                } 
            }
            //assert(iterations == sample.get_actual_size());

        }


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

        trajectory_t& sample;
        range_search_t& search;
        std::string dataset_location;
        id_t last_parsed_trajectory;
        distance_t& distance_threshold;
        int integer_frequency_threshold;
        distance_t& curve_simplification_factor = 0;
        const frechet::internal::curve_simplification<space> simplification;
};


}