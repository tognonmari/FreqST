#pragma once

#include <cmath>
#include <iostream>
#include <limits>
#include <list>
#include <memory>
#include <optional>
#include <vector>
#include <random>
#include <chrono>
#include "free_space_graph_free_axis.h"
#include "frechet_distance.h"
#include "free_space_graph.h"
#include "kdtree_range_search.h"
#include "metric_space.h"
#include "profiling.h"
#include "io.h"
#include "subtrajectory_cluster.h"
#include "subtrajectory_routine_bbgll.h"
#include "trajectory.h"
#include "canonical_pathlets.h"
namespace chrono = std::chrono;
namespace frechet{

template<metric_space m_space>
class freq_subtrajectory_sampler{

    public:
    using space = m_space;
    using range_search_t = kd_tree_range_search<space>;

    using point_t = space::point_t;
    using distance_function_t = space::distance_function_t;
    using distance_t = distance_function_t::distance_t;
    
    using trajectory_t = trajectory_collection<space>;
    using index_t = trajectory_t::index_t;
    using subtrajectory_t = trajectory_t::subtrajectory_t;

    using subtrajectory_cluster_t = subtrajectory_cluster<space>;


    using id_t = trajectory_t::id_t;
    
    
    public:
        freq_subtrajectory_sampler(const trajectory_t& trajectory,
            float eps, 
            float del, 
            distance_t radius, int minimum_length, int random_seed) : the_trajectory(trajectory), epsilon(eps), delta(del), distance_threshold(radius), min_length(minimum_length), seed(random_seed){

                std::mt19937 seeded_generator(seed);
                this->mt = seeded_generator;
            }

        void generate_chernoff_sample(){

            //Step 1: compute sample size according to Chernoff rule. 
            (this->sampled_trajs_ids).clear();
            int sample_size = (int) (3 / (epsilon * epsilon)) * log(2 * this->total_pathlet_number_respecting_ids() / delta);
            std::cout << "Chernoff sample size with espilon "<<epsilon << ",delta "<< delta <<" is: "<< sample_size <<std::endl;
            //Step 2: assert sampling is worthwhile
            if(sample_size > the_trajectory.num_trajectories()){

                std::cerr << "Chernoff Bound was too loose for your dataset."<< std::endl;

                std::exit(1);

            }
            
            //Step 3: Sample indexes with replacement
            
            this->sample_trajectories(sample_size); //FILLS IN CLASS VARIABLE SAMPLE

        }

        void generate_vc_sample(){

            //Step 1: compute sample size according to Chernoff rule. 
            (this->sampled_trajs_ids).clear();
            int sample_size = (int) (2 / (epsilon * epsilon)) * (this->vc_dim() + log(1 / delta));
            std::cout << "VCdim sample size with espilon "<<epsilon << ",  delta "<< delta <<", radius "<< distance_threshold<< " is: "<< sample_size <<std::endl;
            //Step 2: assert sampling is worthwhile
            if(sample_size > the_trajectory.num_trajectories()){

                std::cerr << "VC Bound was too loose for your dataset."<< std::endl;

                std::exit(1);

            }
            //Step 3: Sample indexes with replacement
            this->sample_trajectories(sample_size);

        }

        void generate_rough_vc_sample(){

            //Step 1: compute sample size according to Chernoff rule. 
            (this->sampled_trajs_ids).clear();
            int sample_size = (int) (2 / (epsilon * epsilon)) * (this->rough_vc_dim() + log(1 / delta));
            std::cout << " Rough VCdim sample size with espilon "<<epsilon << ",  delta "<< delta <<", radius "<< distance_threshold<< " is: "<< sample_size <<std::endl;
            //Step 2: assert sampling is worthwhile
            if(sample_size > the_trajectory.num_trajectories()){

                std::cerr << "VC Bound was too loose for your dataset."<< std::endl;

                std::exit(1);

            }
            //Step 3: Sample indexes with replacement
            this->sample_trajectories(sample_size);

        }

        void dump_sample_to_file(std::string filename){
            std::cout << "Started dumping the sample to a file"<< std::endl;
            assert(!sampled_trajs_ids.empty());
            std::ofstream fout(filename);
            for (int j = 0; j< sampled_trajs_ids.size(); j++){
                std::cout << "I am printing the sample "<< j <<std::endl;
                this->print_subtrajectory_to_file(fout, sampled_trajs_ids.at(j));

            }
            return;
        }



    private:
    int rough_vc_dim(){
        range_search_t search{the_trajectory};
        std::vector<int> c;
        index_t last_seen_trajectory = the_trajectory.get_id_at(0);
        int counter;
        
        for(index_t i =0; i<=the_trajectory.get_actual_size(); i++){
            if(i%10000 == 0){

                
                std::cout<< "Processing point "<< i<< " to find the c bound" << std::endl;

            }
            if(the_trajectory.get_id_at(i) == last_seen_trajectory){

                counter += search.search(i, this->distance_threshold).size();

            }
            else{

                // Append the result up to now to c
                c.push_back(floor(log2(counter) + 1));
                //initialize the set again 
                counter = 0;
                last_seen_trajectory = the_trajectory.get_id_at(i);
                counter += search.search(i, this->distance_threshold).size();
                //add info for the current point


            }

        }
        std::cout << "Started sorting "<< std::endl;
        std::sort(c.begin(),c.end(), std::greater<>());
        std::cout << "Finished sorting "<< std::endl;
        int vc_dim = 0;
        
        for(int i = 0 ; i < c.size(); i++){

            if(vc_dim < c.at(i)){

                vc_dim++;

            }

        }

        return vc_dim;

    }
    int vc_dim(){

        //Compute VC Dimension 
        range_search_t search{the_trajectory};
        std::vector<int> c;
        index_t last_seen_trajectory = the_trajectory.get_id_at(0);
        std::set<index_t> traj_set;
        
        for(index_t i =0; i<=the_trajectory.get_actual_size(); i++){
            if(i%10000 == 0){

                
                std::cout<< "Processing point "<< i<< " to find the c bound" << std::endl;

            }
            if(the_trajectory.get_id_at(i) == last_seen_trajectory){

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
                for (const auto idx: search.search(i, this->distance_threshold)) {

                    traj_set.insert(idx);

                }


            }

        }
        std::cout << "Started sorting "<< std::endl;
        std::sort(c.begin(),c.end(), std::greater<>());
        std::cout << "Finished sorting "<< std::endl;
        int vc_dim = 0;
        
        for(int i = 0 ; i < c.size(); i++){

            if(vc_dim < c.at(i)){

                vc_dim++;

            }

        }

        return vc_dim;
    }

    void print_subtrajectory_to_file(std::ofstream& fout, id_t& id){

        size_t n = the_trajectory.num_trajectories();
        index_t j = the_trajectory.get_first_point_in_trajectory(id%n);
        
        while( the_trajectory.get_id_at(j) == id%n){
            j++;
            fout << the_trajectory[j].x()<< " "<< the_trajectory[j].y()<<" "<< id<<std::endl;
        }

        return;
    }

    int total_pathlet_number_respecting_ids(){
        int total = 0;

        for (id_t i = 0; i < the_trajectory.num_trajectories(); i ++){

            total += 2 * ceil(log2( the_trajectory.get_trajectory_size(i) / min_length));

        }

        return total;

    }
    //return type should be void-> TODO:correct this!!!
    std::vector<id_t> sample_trajectories(int sample_size){

        if(! sampled_trajs_ids.empty()){

            sampled_trajs_ids.clear();

        }
        int n = the_trajectory.num_trajectories();
        // Extract sampled ids
        for (int i = 0; i < sample_size; i++){

            sampled_trajs_ids.push_back((id_t)(mt())% n);

        }
        std::cout << "Started sorting "<< std::endl;
        std::sort(sampled_trajs_ids.begin(),sampled_trajs_ids.end());
        std::cout << "Finished sorting "<< std::endl;
        id_t last_read_trajectory = sampled_trajs_ids.at(0);
        int repetitions = 0;
        //ri scorriamo il vettore. ad ogni ripetizione i dello stesso id assegnamo un nuovo id 
        for( int j = 1; j< sample_size; j++){
            
            if (sampled_trajs_ids.at(j) == last_read_trajectory){

                repetitions++;
                sampled_trajs_ids.at(j) += (id_t) (repetitions * n);

            }
            else{
                last_read_trajectory = sampled_trajs_ids.at(j);
                repetitions = 0;
            }
        }
        std::cout<< "I have finished the method to get the ids"<< std::endl;
        return sampled_trajs_ids;
    }

    std::mt19937 mt;
    std::vector<id_t> sampled_trajs_ids;
    trajectory_t the_trajectory;
    distance_t distance_threshold;
    float epsilon;
    float delta;
    int min_length;
    bool performed_sampling;
    int seed;
    
};

template<metric_space m_space>
class frequent_subtrajectory_algo{
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

    private:
        struct frequent_pathlet{

            std::pair<index_t,index_t> extremes;
            id_t pathlet_mother;
            float frequency;

        };

    public:

        frequent_subtrajectory_algo(trajectory_t sampled_traj, std::string dataset_file, float frequency_threshold, distance_t distance_thresh) : search(sampled_traj){
            this->sample = sampled_traj;
            this->dataset_location = dataset_file;
            this->integer_frequency_threshold = ceil(frequency_threshold * this->sample.num_trajectories());
            std::cout << "THE INTEGER FREQ THRESHOLD IS "<< this->integer_frequency_threshold<<std::endl;
            this-> last_parsed_trajectory = -1;
            this-> distance_threshold = distance_thresh;
            //this->populate_range_search_tree_with_sample_points();
        }

        void populate_range_search_tree_with_sample_points(){

            range_search_t rs(sample);
            this->search = rs;
            return;
        }

        void compute_all_frequent_pathlets(){

            //open full dataset file 
            std::ifstream input_stream(this->dataset_location);
            //std::cout <<"Starting reading the transactions."<<std::endl;
            while(!input_stream.eof()){
                auto start = chrono::high_resolution_clock::now();
                trajectory_t pathlet_mother = this->read_next_transaction_from_file(input_stream);
                //std::cout <<"Parsed a transaction."<<std::endl;
                std::cout<<" The transaction has ID "<<pathlet_mother.get_id_at(pathlet_mother.get_actual_size()-1)<<std::endl;
                std::cout <<" I have this many points : "<< pathlet_mother.get_actual_size()<<std::endl;
                BinaryPathletTree pathlet_tree(pathlet_mother, pathlet_mother.get_id_at(0),floor(log2(pathlet_mother.total_size())) + 1,1);
                
                free_space_graph_t fsg(0);
                
                //for all the columns of the bst populate the column
                this->populate_all_columns(fsg, pathlet_mother);
                //std::cout <<"Populated the columns."<<std::endl;
                //maybe i need to rewrite the kd tree to access with the coordinates directly
                this->collect_all_frequent_pathlets(fsg, pathlet_tree);
                auto stop = chrono::high_resolution_clock::now();
                auto duration = duration_cast<chrono::milliseconds>(stop - start);
                std::cout<< "TIME : "<< duration.count()<< std::endl;
            }



        }


        void compute_maximal_frequent_pathlets(){

            //open full dataset file 
            std::ifstream input_stream(this->dataset_location);
            std::cout <<"Starting reading the transactions."<<std::endl;
            while(!input_stream.eof()){

                trajectory_t pathlet_mother = this->read_next_transaction_from_file(input_stream);
                //std::cout <<"Parsed a transaction."<<std::endl;
                //std::cout<<" The transaction has ID "<<pathlet_mother.get_id_at(pathlet_mother.get_actual_size()-1)<<std::endl;
                std::cout <<" I have this many points : "<< pathlet_mother.get_actual_size()<<std::endl;
                BinaryPathletTree pathlet_tree(pathlet_mother, pathlet_mother.get_id_at(0),floor(log2(pathlet_mother.total_size())) + 1,1);
                
                free_space_graph_t fsg(0);
                
                //for all the columns of the bst populate the column
                this->populate_all_columns(fsg, pathlet_mother);
                std::cout <<"Populated the columns."<<std::endl;
                //maybe i need to rewrite the kd tree to access with the coordinates directly
                this->collect_maximal_frequent_pathlets(fsg, pathlet_tree);

            }



        }

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

        void collect_maximal_frequent_pathlets(free_space_graph_t& fsg, binary_pathlet_tree_t& pathlet_tree){

            int d = pathlet_tree.getDepth(); //Last filled level
            //std::cout<< "tree has depth "<< d <<std::endl; // assertion for my toy dataset
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

                    int count = fsg.query_one_pathlet_over_the_sample(sample, pn.getPathlet()); //SF IS HERE
                    
                    if(count < this->integer_frequency_threshold){
                        
                        pathlet_tree.setInfrequent(position);

                    }
                    else{

                        pathlet_tree.setEstimatedFrequency(position, ((float) count / num_sampled_trajs));

                    }

                }


            }

            //TOP DOWN VISIT TO COLLECT THE MAXIMAL FREQUENT PATHLETS ON  THE BST
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


        }
        

        void collect_all_frequent_pathlets(free_space_graph_t& fsg, binary_pathlet_tree_t& pathlet_tree){

            int d = pathlet_tree.getDepth(); //Last filled level
            //std::cout<< "tree has depth "<< d <<std::endl; // assertion for my toy dataset
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

                    int count = fsg.query_one_pathlet_over_the_sample(sample, pn.getPathlet()); //SF IS HERE
                    
                    if(count < this->integer_frequency_threshold){
                        
                        pathlet_tree.setInfrequent(position);

                    }
                    else{
                        //TODO: CLEANUP
                        pathlet_tree.setEstimatedFrequency(position, ((float)count / num_sampled_trajs));
                        assert(pathlet_tree.getNodeAt(position).frequency - ((float)count / num_sampled_trajs) < 0.0001 );
                        frequent_pathlet just_found;
                        just_found.extremes = pn.getPathlet();
                        just_found.pathlet_mother = (pathlet_tree.getTrajectoryId());
                        just_found.frequency = ((float) count / num_sampled_trajs);
                        //std::cout<<"FREQUENCY: "<<pathlet_tree.getNodeAt(position).frequency<< std::endl;
                        freq_pathlets.push_back(just_found);
                    }

                }


            }


        }

        void populate_all_columns(free_space_graph_t& fsg, const trajectory_t& pathlet_mother){

            int num_col = pathlet_mother.get_actual_size();

            for (int j = 0; j< num_col; j++){

                populate_column(fsg, pathlet_mother[j], this->distance_threshold, j);
                if(j < num_col -1){

                    fsg.new_column();

                }
            }

            return;
        }

        void populate_column(free_space_graph_t &fsg, const point_t& point, const distance_t& query_distance, index_t column_index ){

            TIME_BEGIN();
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
            assert(iterations == sample.get_actual_size());
            //assert(zeroes > 0); NOT NEEDED 
            TIME_END(populate_column);

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

        trajectory_t sample;
        range_search_t search;
        std::string dataset_location;
        id_t last_parsed_trajectory;
        distance_t distance_threshold;
        int integer_frequency_threshold;
        std::vector<frequent_pathlet> freq_pathlets;
};


}