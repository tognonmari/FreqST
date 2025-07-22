#pragma once

#include <cmath>
#include <iostream>
#include <limits>
#include <list>
#include <memory>

#include "frechet_distance.h"
#include "free_space_graph.h"
#include "free_space_graph_incremental.h"
#include "kdtree_range_search.h"
#include "kcluster_detail.h"
#include "metric_space.h"
#include "profiling.h"
#include "subtrajectory_cluster.h"
#include "subtrajectory_routine_bbgll.h"
#include "subtrajectory_routine_rightstep.h"
#include "trajectory.h"
#include "freq_st_algo.h"
namespace frechet {
    
template<metric_space m_space>
class subtrajectory_clustering_algo {

public:
    using space = m_space;
    using free_space_graph_t = sparse_free_space_graph<space>;
    using range_search_t = kd_tree_range_search<space>;
    using free_space_graph_incremental_t = sparse_free_space_graph_incremental<space>;
    using point_t = space::point_t;
    using distance_function_t = space::distance_function_t;
    using distance_t = distance_function_t::distance_t;

    using trajectory_t = trajectory_collection<space>;
    using index_t = trajectory_t::index_t;
    using subtrajectory_t = trajectory_t::subtrajectory_t;

    using subtrajectory_cluster_t = subtrajectory_cluster<space>;

    using efficacy_factor_t = internal::efficacy_factor_vec<distance_t>;

    using frequent_pathlet = frechet::frequent_subtrajectory_algo<space>::frequent_pathlet;
    using frequent_subtrajectory_algo_t = frechet::frequent_subtrajectory_algo<space>;
    using cluster_quality_t = sparse_free_space_graph_incremental<space>::cluster_quality;
private:
    using k_cluster_detail = internal::k_cluster_detail<space>;

    using fixed_d_cluster = internal::fixed_distance_clustering<space>;

public:
    // set min_distance or max_distance to -1 to have the algorithm compute it.
    subtrajectory_clustering_algo(const trajectory_t &trajectory,
            distance_t min_distance = -1,
            distance_t max_distance = -1,
            const efficacy_factor_t &efficacy_factors = {1, 1, 1, false},
            const rightstep_config &config = {}) :
                    trajectory(trajectory),
                    range_search(this->trajectory),  // Using `this->trajectory`,
                                                     // otherwise `range_search` has a reference
                                                     // to a different trajectory than we actually use.
                    distance_limit(max_distance),
                    efficacy_factors(efficacy_factors),
                    config(config)
    {
        initialize_distance_limits(min_distance, max_distance);
    }

    // Perform clustering for $EVAL = CENTER$.
    // If `efficacy_factors.ignore_point_clusters` is true, clusters with empty reference trajectory are ignored when estimating efficacy;
    // instead, the points in these clustered are treated as un-clustered
    void perform_center_clustering() {
        // For k-centers we can directly maximize the coverage, skipping the (slow) cost computation.
        config.cost_per_pathlet = 0;
        std::vector<std::unique_ptr<fixed_d_cluster>> clustering_algos;
        for (const auto &dist: sq_distances) {
            clustering_algos.emplace_back(std::make_unique<fixed_d_cluster>(trajectory, dist));
        }
        #pragma omp parallel for
        for (size_t i = 0; i < clustering_algos.size(); ++i) {
            clustering_algos[i]->perform_clustering_rightstep(config);
            clustering_algos[i]->drop_inefficient_clusters_center(efficacy_factors);
            std::cout << "Clustered at distance " << clustering_algos[i]->get_sq_distance() << std::endl;
        }
        distance_t best_efficacy = std::numeric_limits<distance_t>::infinity();
        fixed_d_cluster* algo_with_best_result = nullptr;
        std::cout << "Finding the best clustering...\n";
        for (auto &algo: clustering_algos) {
            auto efficacy = k_cluster_detail::compute_efficacy_center(trajectory, algo->get_clusters(), efficacy_factors);
            std::cout << " Best so far: " << best_efficacy << ", current: " << efficacy << ", at distance " << algo->get_sq_distance() << "\n";
            if (efficacy < best_efficacy) {
                best_efficacy = efficacy;
                algo_with_best_result = algo.get();
            }
        }
        std::cout << "Efficacy of best found clustering: " << best_efficacy << ", at distance " << algo_with_best_result->get_sq_distance() << "\n";
        assert(algo_with_best_result != nullptr);
        pathlets = algo_with_best_result->get_clusters();
    }






    void perform_aided_centers_clustering(std::string infilename, std::string samplefilename){
        int first_uncovered_pathlet_by_distance = 0;
        std::vector<frequent_pathlet> max_freq_pathlets_by_distance;
        this->config.cost_per_pathlet = 0;
        std::vector<std::unique_ptr<fixed_d_cluster>> clustering_algos;
        for (const auto &dist: sq_distances) {
            clustering_algos.emplace_back(std::make_unique<fixed_d_cluster>(trajectory, dist));
        }
        #pragma omp parallel for
        for (size_t i = 0; i < clustering_algos.size(); ++i) {
            const auto &dist = sq_distances.at(i);
            trajectory_t sample = read_trajectory_from_file<space>(samplefilename); //samplefile needs to be passed as input
            std::cout << "Read the sample, now mining frequent pathlets"<< std::endl;
            // Compute the frequent pathlets, then swap the values of the freq pathlets vector into the new one before deletion
            float radius = std::sqrt(dist);
            range_search_t rs(sample);
            float frequency_threshold = 0.1;
            frequent_subtrajectory_algo_t algo(sample, rs, infilename, frequency_threshold, radius); //infilename needs to be passed as parameter to the clustering algorithm
            algo.compute_maximal_frequent_pathlets();
            std::swap(algo.freq_pathlets, max_freq_pathlets_by_distance);
            //TODO:this should not be a vector, but a hashmap, must be modified in pathlet tree etc... Unluckily simply changing this ds: Keep both of them for now
            std::sort(max_freq_pathlets_by_distance.begin(), max_freq_pathlets_by_distance.end(), pathlet_sorter_by_length);
            //auto rng = std::default_random_engine {};
            //std::ranges::shuffle(max_freq_pathlets_by_distance.at(i), rng);
            first_uncovered_pathlet_by_distance = (max_freq_pathlets_by_distance.size() - 1); //initialize with the last pathlet
            if(i==sq_distances.size() -1){

                std::cout << "Found this longest frequent pathlet " << max_freq_pathlets_by_distance.at(first_uncovered_pathlet_by_distance).extremes.first << " " << max_freq_pathlets_by_distance.at(first_uncovered_pathlet_by_distance).extremes.second << " with mother " << max_freq_pathlets_by_distance.at(first_uncovered_pathlet_by_distance).pathlet_mother <<  std::endl;
            }
            subtrajectory_cluster_t output_cluster;
            while(first_uncovered_pathlet_by_distance >=0){
                std::pair<index_t, index_t> offsets = max_freq_pathlets_by_distance.at(first_uncovered_pathlet_by_distance).extremes;
                id_t pathlet_mother = max_freq_pathlets_by_distance.at(first_uncovered_pathlet_by_distance).pathlet_mother;
                subtrajectory_t indexes_in_trajectory = from_freq_pathlet_to_subtrajectory(offsets, pathlet_mother);
                //Now build a free_space_graph_free_axis, and extract the cluster as they do.
                cluster_quality_t cl_qual(0,0,sq_distances[i], indexes_in_trajectory.first, indexes_in_trajectory.second);
                trajectory_t temp = clustering_algos[i]->get_trajectory();
                subtrajectory_clustering_rightstep<space> clustering_algos_rightstep(temp, std::vector<index_t>(temp.total_size(), index_t{1}), range_search, config);
                output_cluster = clustering_algos_rightstep.cluster_from_candidate(sq_distances[i], cl_qual);
            
                clustering_algos[i]->establish_cluster(output_cluster);

                //Update max_freq_pathlet_by_dist and the cursor
                
                if(!max_freq_pathlets_by_distance.empty()){
                    bool found = false;
                    for(int j = first_uncovered_pathlet_by_distance; j>=0; j--){

                        //convert pathlet to trajectory interval:
                        std::pair<index_t, index_t> offsets = max_freq_pathlets_by_distance[j].extremes;
                        id_t pathlet_mother = max_freq_pathlets_by_distance[j].pathlet_mother;
                        subtrajectory_t indexes_in_trajectory = from_freq_pathlet_to_subtrajectory(offsets, pathlet_mother);

                        //Now, check whether any point in the frequent pathlet is covered.
                        bool covered = false;
                        trajectory_t traj = clustering_algos[i]->get_trajectory();
                        for(index_t h = indexes_in_trajectory.first; h <= indexes_in_trajectory.second; h++){
                            if(traj.is_point_deleted(h)){
                                covered = true;
                                break;
                            }
                        }

                        if(!covered){
                            //If not, we can use this pathlet as center.
                            first_uncovered_pathlet_by_distance = j;
                            found = true;
                            //std::cout << "First uncovered pathlet for distance " << std::sqrt(sq_distances[i]) << " is at index " << j << "\n";
                            break;
                        }                 

                    }
                    if(!found){

                        first_uncovered_pathlet_by_distance = -1; //no uncovered pathlet found

                    }

                }
            }
            
            clustering_algos[i]->perform_clustering_rightstep(config);
            clustering_algos[i]->drop_inefficient_clusters_center(efficacy_factors);
            std::cout << "Clustered at distance " << clustering_algos[i]->get_sq_distance() << std::endl;
            max_freq_pathlets_by_distance.clear();
            first_uncovered_pathlet_by_distance = 0;
        }

        distance_t best_efficacy = std::numeric_limits<distance_t>::infinity();
        fixed_d_cluster* algo_with_best_result = nullptr;
        std::cout << "Finding the best clustering...\n";
        for (auto &algo: clustering_algos) {
            auto efficacy = k_cluster_detail::compute_efficacy_center(trajectory, algo->get_clusters(), efficacy_factors);
            std::cout << " Best so far: " << best_efficacy << ", current: " << efficacy << ", at distance " << algo->get_sq_distance() << "\n";
            if (efficacy < best_efficacy) {
                best_efficacy = efficacy;
                algo_with_best_result = algo.get();
            }
        }
        std::cout << "Efficacy of best found clustering: " << best_efficacy << ", at distance " << algo_with_best_result->get_sq_distance() << "\n";
        assert(algo_with_best_result != nullptr);
        pathlets = algo_with_best_result->get_clusters();

    }

    
    // Perform clustering for $EVAL = CENTER$.
    // If `efficacy_factors.ignore_point_clusters` is true, clusters with empty reference trajectory are ignored when estimating efficacy;
    // instead, the points in these clustered are treated as un-clustered
    void perform_means_clustering() {
        // The distance isn't fixed yet, so we will multiply with it later.
        config.cost_per_pathlet = efficacy_factors.c_2 / efficacy_factors.c_1; 
        std::vector<std::unique_ptr<fixed_d_cluster>> clustering_algos;
        std::vector<std::pair<bool, subtrajectory_cluster_t> > candidate_clusters;
        std::vector<distance_t> gamma;
        for (const auto &dist: sq_distances) {
            clustering_algos.emplace_back(std::make_unique<fixed_d_cluster>(trajectory, dist));
            candidate_clusters.emplace_back();
            gamma.emplace_back();
        }
        while (clustering_algos.front()->count_remaining_points() > 0) {
            #pragma omp parallel for
            for (size_t i = 0; i < clustering_algos.size(); ++i) {
                auto &[success, output_cluster] = candidate_clusters[i];
                success = clustering_algos[i]->find_best_cluster_rightstep(output_cluster, config);
                assert(success);
                //std::cout <<" Started computing gamma, supposedly slow"<< std::endl;
                gamma[i] = clustering_algos[i]->compute_gamma(output_cluster, efficacy_factors);
            }
            // Pick best cluster as in Section 4.3 of Agarwal et. al, 2018.
            size_t best_i = std::max_element(gamma.begin(), gamma.end()) - gamma.begin();
            auto &best_cluster = candidate_clusters[best_i].second;
            k_cluster_detail::prune_inefficient_subtrajectories(trajectory, best_cluster, gamma[best_i], efficacy_factors);
            // On very dense data sets, this can happen when very few points remain, due to floating point inaccuracies.
            // The remaining points are better left unclustered, assuming c1 > 0.
            if (best_cluster.get_subtrajectories().empty()) {
                break;
            }
            std::cout << "Best cluster has "
                << "distance: " << std::sqrt(sq_distances[best_i])
                << ", gamma: " << gamma[best_i]
                << ", vertices: " << best_cluster.number_of_vertices() << "\n";

            // Each algo has it's own trajectory and range search object.
            // I think this is required for the parallel for loop.
            for (auto &algo : clustering_algos) {
                algo->establish_cluster(best_cluster);
            }
        }
        // all algorithms store the same clustering.
        clustering_algos.front()->drop_inefficient_clusters_means(efficacy_factors);
        pathlets = clustering_algos.front()->get_clusters();
    }
    void perform_aided_means_clustering_trying_all_pathlets(std::string infilename, std::string samplefilename){
        std::cout<< "Performing sample aided clustering..." << std::endl;
        // The distance isn't fixed yet, so we will multiply with it later.
        config.cost_per_pathlet = efficacy_factors.c_2 / efficacy_factors.c_1; 
        std::vector<std::unique_ptr<fixed_d_cluster>> clustering_algos;
        std::vector<std::pair<bool, subtrajectory_cluster_t> > candidate_clusters;
        std::vector<distance_t> gamma;
        
        std::vector<std::set<frequent_pathlet>> max_freq_pathlets_by_distance;
        std::vector<int> num_uncovered_pathlets;
        for (const auto &dist: sq_distances) {
            std::cout<< "Initializing for distance "<<dist << std::endl;
            clustering_algos.emplace_back(std::make_unique<fixed_d_cluster>(trajectory, dist));
            std::cout<< "Finished Initializing cl algos for distance "<<dist << std::endl;
            candidate_clusters.emplace_back();
            std::cout<< "Finished Initializing candidate clusters for distance "<<dist << std::endl;
            gamma.emplace_back();
            std::cout<< "Finished Initializing their ds for distance "<<dist << std::endl;
            max_freq_pathlets_by_distance.emplace_back();
            std::cout<< "Finished Initializing for distance "<<dist << std::endl;
        }
        std::cout<< "Initialized data structures for aided means clustering." << std::endl;
        //Do the sampling and initialization of the vector of frequent sts: we use the same chenoff sample?? For now, yes
        //Generate the sample
        trajectory_t sample = read_trajectory_from_file<space>(samplefilename); //samplefile needs to be passed as input
        std::cout << "Read the sample, now mining frequent pathlets"<< std::endl;
        float frequency_threshold = 0.1;
        for(int i = 0; i<sq_distances.size(); i++){
            const auto &dist = sq_distances.at(i);
            // Compute the frequent pathlets, then swap the values of the freq pathlets vector into the new one before deletion
            float radius = std::sqrt(dist);
            range_search_t rs(sample);
            
            frequent_subtrajectory_algo_t algo(sample, rs, infilename, frequency_threshold, radius); //infilename needs to be passed as parameter to the clustering algorithm
            algo.compute_maximal_frequent_pathlets();
            //std::swap(algo.freq_pathlets, max_freq_pathlets_by_distance.at(i));
            //TODO:this should not be a vector, but a hashmap, must be modified in pathlet tree etc... Unluckily simply changing this ds: Keep both of them for now
            std::set<frequent_pathlet> temp(algo.freq_pathlets.begin(), algo.freq_pathlets.end());

            std::swap(temp, max_freq_pathlets_by_distance.at(i));
            std::cout <<(max_freq_pathlets_by_distance[i].size())<< ","<< algo.freq_pathlets.size()<<std::endl;
            
        }
        std::cout << "Frequent pathlets computed for all distances.\n";
        //TODO: rewrite this because this is horrible
        while (clustering_algos.front()->count_remaining_points() > 0) {
            
            #pragma omp parallel for
            for (size_t i = 0; i < clustering_algos.size(); ++i) {
                
                auto &[success, output_cluster] = candidate_clusters[i];
                
                if(!max_freq_pathlets_by_distance[i].empty()){
                    //extract the frequent pathlet as center
                    success = true;
                    distance_t temp_gamma = std::numeric_limits<float>::max();
                    subtrajectory_cluster_t candidate_cluster;
                    for(auto& cluster_center : max_freq_pathlets_by_distance[i]){
                        
                        std::pair<index_t, index_t> offsets = cluster_center.extremes;
                        id_t pathlet_mother = cluster_center.pathlet_mother;
                        
                        
                        subtrajectory_t indexes_in_trajectory = from_freq_pathlet_to_subtrajectory(offsets, pathlet_mother);
                        //Now build a free_space_graph_free_axis, and extract the cluster as they do.
                        cluster_quality_t cl_qual(0,0,sq_distances[i], indexes_in_trajectory.first, indexes_in_trajectory.second);
                        trajectory_t temp = clustering_algos[i]->get_trajectory();
                        assert(!temp.is_point_deleted(indexes_in_trajectory.first));
                        subtrajectory_clustering_rightstep<space> clustering_algos_rightstep(temp, std::vector<index_t>(temp.total_size(), index_t{1}), range_search, config);
                        candidate_cluster = clustering_algos_rightstep.cluster_from_candidate(sq_distances[i], cl_qual);
                        distance_t cluster_gamma = clustering_algos[i]->compute_gamma(output_cluster, efficacy_factors);
                        if (cluster_gamma < temp_gamma){

                            temp_gamma = cluster_gamma;
                            std::swap(candidate_cluster, output_cluster);

                        }
                        
                    }
                }
                else{
                    
                // If no uncovered pathlet is available, we can use the clustering algorithm to find a cluster.
                
                success = clustering_algos[i]->find_best_cluster_rightstep(output_cluster, config);
                assert(success);

                }
                gamma[i] = clustering_algos[i]->compute_gamma(output_cluster, efficacy_factors);
            }
            // Pick best cluster as in Section 4.3 of Agarwal et. al, 2018.
            size_t best_i = std::max_element(gamma.begin(), gamma.end()) - gamma.begin();
            auto &best_cluster = candidate_clusters[best_i].second;
            k_cluster_detail::prune_inefficient_subtrajectories(trajectory, best_cluster, gamma[best_i], efficacy_factors);
            //TODO: Once the best cluster is selected, delete the unwanted pathlets from the frequent files. I should have a nice way of compputing these
            // On very dense data sets, this can happen when very few points remain, due to floating point inaccuracies.
            // The remaining points are better left unclustered, assuming c1 > 0.
            if (best_cluster.get_subtrajectories().empty()) {
                break;
            }
            std::cout << "Best cluster has "
                << "distance: " << std::sqrt(sq_distances[best_i])
                << ", gamma: " << gamma[best_i]
                << ", vertices: " << best_cluster.number_of_vertices() << "\n";
            
            // Each algo has it's own trajectory and range search object.
            // I think this is required for the parallel for loop.
            for (auto &algo : clustering_algos) {
                best_cluster = best_cluster;
                //std::cout << "Deleting from algo"<< std::endl;
                algo->establish_cluster(best_cluster);
            }
            std::cout << "UNCOVERED POINTS : "<< clustering_algos.front()->count_remaining_points() << std::endl;
            //Find first uncovered frequent pathlet for each distance
            for (size_t i = 0; i < max_freq_pathlets_by_distance.size(); i++){

                if(!max_freq_pathlets_by_distance[i].empty()){
                    bool found = false;
                    for(auto iter = max_freq_pathlets_by_distance[i].begin(); iter != max_freq_pathlets_by_distance[i].end();){
                        bool covered = false;
                        frequent_pathlet p= *(iter);
                        trajectory_t traj = clustering_algos[i]->get_trajectory();
                        std::pair<index_t, index_t> offsets = p.extremes;
                        id_t pathlet_mother = p.pathlet_mother;
                        subtrajectory_t indexes_in_trajectory = from_freq_pathlet_to_subtrajectory(offsets, pathlet_mother);
                        for(index_t h = indexes_in_trajectory.first; h <= indexes_in_trajectory.second; h++){
                            if(traj.is_point_deleted(h)){
                                covered = true;
                                break;
                            }
                        }
                        if(covered){
                            iter = max_freq_pathlets_by_distance[i].erase(iter);
                        }
                        else{
                            ++iter;
                        }
                    }

                    /*
                    for  (auto& p:max_freq_pathlets_by_distance[i]){
                        std::cout << "Iterating over the distance "<< i <<std::endl;
                        bool covered = false;
                        
                        trajectory_t traj = clustering_algos[i]->get_trajectory();
                        std::pair<index_t, index_t> offsets = p.extremes;
                        id_t pathlet_mother = p.pathlet_mother;
                        subtrajectory_t indexes_in_trajectory = from_freq_pathlet_to_subtrajectory(offsets, pathlet_mother);

                        for(index_t h = indexes_in_trajectory.first; h <= indexes_in_trajectory.second; h++){
                            if(traj.is_point_deleted(h)){
                                covered = true;
                                break;
                            }
                        }
                        if(covered){

                            max_freq_pathlets_by_distance[i].erase((p));

                        }
                        if(max_freq_pathlets_by_distance[i].empty()){

                            break;
                        }
                        
                    }
                    
                    */
                    
                }
            }            
        
        }
        // all algorithms store the same clustering.
        clustering_algos.front()->drop_inefficient_clusters_means(efficacy_factors);
        pathlets = clustering_algos.front()->get_clusters();




    }
    void perform_aided_means_clustering(std::string infilename, std::string samplefilename) {
        std::cout<< "Performing sample aided clustering..." << std::endl;
        // The distance isn't fixed yet, so we will multiply with it later.
        config.cost_per_pathlet = efficacy_factors.c_2 / efficacy_factors.c_1; 
        std::vector<std::unique_ptr<fixed_d_cluster>> clustering_algos;
        std::vector<std::pair<bool, subtrajectory_cluster_t> > candidate_clusters;
        std::vector<distance_t> gamma;
        
        std::vector<std::vector<frequent_pathlet>> max_freq_pathlets_by_distance;
        std::vector<int> first_uncovered_pathlet_by_distance;
        for (const auto &dist: sq_distances) {
            std::cout<< "Initializing for distance "<<dist << std::endl;
            clustering_algos.emplace_back(std::make_unique<fixed_d_cluster>(trajectory, dist));
            std::cout<< "Finished Initializing cl algos for distance "<<dist << std::endl;
            candidate_clusters.emplace_back();
            std::cout<< "Finished Initializing candidate clusters for distance "<<dist << std::endl;
            gamma.emplace_back();
            std::cout<< "Finished Initializing their ds for distance "<<dist << std::endl;
            max_freq_pathlets_by_distance.emplace_back();
            std::cout<< "Finished Initializing for distance "<<dist << std::endl;
        }
        std::cout<< "Initialized data structures for aided means clustering." << std::endl;
        //Do the sampling and initialization of the vector of frequent sts: we use the same chenoff sample?? For now, yes
        //Generate the sample
        trajectory_t sample = read_trajectory_from_file<space>(samplefilename); //samplefile needs to be passed as input
        std::cout << "Read the sample, now mining frequent pathlets"<< std::endl;
        for(int i = 0; i<sq_distances.size(); i++){
            const auto &dist = sq_distances.at(i);
            // Compute the frequent pathlets, then swap the values of the freq pathlets vector into the new one before deletion
            float radius = std::sqrt(dist);
            range_search_t rs(sample);
            float frequency_threshold = 0.05;
            frequent_subtrajectory_algo_t algo(sample, rs, infilename, frequency_threshold, radius); //infilename needs to be passed as parameter to the clustering algorithm
            algo.compute_maximal_frequent_pathlets();
            std::swap(algo.freq_pathlets, max_freq_pathlets_by_distance.at(i));
            //TODO:this should not be a vector, but a hashmap, must be modified in pathlet tree etc... Unluckily simply changing this ds: Keep both of them for now
            std::sort(max_freq_pathlets_by_distance[i].begin(), max_freq_pathlets_by_distance[i].end(), pathlet_sorter_by_length);
            std::cout <<(max_freq_pathlets_by_distance[i].size())<<std::endl;
            //auto rng = std::default_random_engine {};
            //std::ranges::shuffle(max_freq_pathlets_by_distance.at(i), rng);
            first_uncovered_pathlet_by_distance.push_back(max_freq_pathlets_by_distance[i].size() - 1); //initialize with the last pathlet
            if(i==sq_distances.size() -1){

                std::cout << "Found this longest frequent pathlet " << max_freq_pathlets_by_distance[i].back().extremes.first << " " << max_freq_pathlets_by_distance[i].back().extremes.second << " with mother " << max_freq_pathlets_by_distance[i].back().pathlet_mother <<  std::endl;
            }
        }
        std::cout << "Frequent pathlets computed for all distances.\n";
        //TODO: rewrite this because this is horrible
        while (clustering_algos.front()->count_remaining_points() > 0) {
            
            #pragma omp parallel for
            for (size_t i = 0; i < clustering_algos.size(); ++i) {
                
                auto &[success, output_cluster] = candidate_clusters[i];
                if(first_uncovered_pathlet_by_distance[i]>=0){
                    //extract the frequent pathlet as center
                    success = true;
                    std::pair<index_t, index_t> offsets = max_freq_pathlets_by_distance[i].at(first_uncovered_pathlet_by_distance[i]).extremes;
                    id_t pathlet_mother = max_freq_pathlets_by_distance[i].at(first_uncovered_pathlet_by_distance[i]).pathlet_mother;
                    subtrajectory_t indexes_in_trajectory = from_freq_pathlet_to_subtrajectory(offsets, pathlet_mother);
                    //Now build a free_space_graph_free_axis, and extract the cluster as they do.
                    cluster_quality_t cl_qual(0,0,sq_distances[i], indexes_in_trajectory.first, indexes_in_trajectory.second);
                    trajectory_t temp = clustering_algos[i]->get_trajectory();
                    subtrajectory_clustering_rightstep<space> clustering_algos_rightstep(temp, std::vector<index_t>(temp.total_size(), index_t{1}), range_search, config);
                    output_cluster = clustering_algos_rightstep.cluster_from_candidate(sq_distances[i], cl_qual);
                }
                else{
                    
                // If no uncovered pathlet is available, we can use the clustering algorithm to find a cluster.
                
                success = clustering_algos[i]->find_best_cluster_rightstep(output_cluster, config);
                assert(success);

                }
                gamma[i] = clustering_algos[i]->compute_gamma(output_cluster, efficacy_factors);
            }
            // Pick best cluster as in Section 4.3 of Agarwal et. al, 2018.
            size_t best_i = std::max_element(gamma.begin(), gamma.end()) - gamma.begin();
            auto &best_cluster = candidate_clusters[best_i].second;
            k_cluster_detail::prune_inefficient_subtrajectories(trajectory, best_cluster, gamma[best_i], efficacy_factors);
            //TODO: Once the best cluster is selected, delete the unwanted pathlets from the frequent files. I should have a nice way of compputing these
            // On very dense data sets, this can happen when very few points remain, due to floating point inaccuracies.
            // The remaining points are better left unclustered, assuming c1 > 0.
            if (best_cluster.get_subtrajectories().empty()) {
                break;
            }
            std::cout << "Best cluster has "
                << "distance: " << std::sqrt(sq_distances[best_i])
                << ", gamma: " << gamma[best_i]
                << ", vertices: " << best_cluster.number_of_vertices() << "\n";
            
            // Each algo has it's own trajectory and range search object.
            // I think this is required for the parallel for loop.
            for (auto &algo : clustering_algos) {
                //best_cluster = best_cluster;
                //std::cout << "Deleting from algo"<< std::endl;
                algo->establish_cluster(best_cluster);
            }
            std::cout << "UNCOVERED POINTS : "<< clustering_algos.front()->count_remaining_points() << std::endl;
            //Find first uncovered frequent pathlet for each distance
            for (size_t i = 0; i < max_freq_pathlets_by_distance.size(); i++){

                if(!max_freq_pathlets_by_distance[i].empty()){
                    bool found = false;
                    for(int j = first_uncovered_pathlet_by_distance[i]; j>=0; j--){

                        //convert pathlet to trajectory interval:
                        std::pair<index_t, index_t> offsets = max_freq_pathlets_by_distance[i][j].extremes;
                        id_t pathlet_mother = max_freq_pathlets_by_distance[i][j].pathlet_mother;
                        subtrajectory_t indexes_in_trajectory = from_freq_pathlet_to_subtrajectory(offsets, pathlet_mother);

                        //Now, check whether any point in the frequent pathlet is covered.
                        bool covered = false;
                        trajectory_t traj = clustering_algos[i]->get_trajectory();
                        for(index_t h = indexes_in_trajectory.first; h <= indexes_in_trajectory.second; h++){
                            if(traj.is_point_deleted(h)){
                                covered = true;
                                break;
                            }
                        }

                        if(!covered){
                            //If not, we can use this pathlet as center.
                            first_uncovered_pathlet_by_distance[i] = j;
                            found = true;
                            //std::cout << "First uncovered pathlet for distance " << std::sqrt(sq_distances[i]) << " is at index " << j << "\n";
                            break;
                        }                 

                    }
                    if(!found){

                        first_uncovered_pathlet_by_distance[i] = -1; //no uncovered pathlet found

                    }

                }
            }            
        
        }
        // all algorithms store the same clustering.
        clustering_algos.front()->drop_inefficient_clusters_means(efficacy_factors);
        pathlets = clustering_algos.front()->get_clusters();

    }
    void perform_aided_means_clustering_top_k(std::string infilename, std::string samplefilename, int k){

        std::cout<< "Performing sample aided clustering..." << std::endl;
        // The distance isn't fixed yet, so we will multiply with it later.
        config.cost_per_pathlet = efficacy_factors.c_2 / efficacy_factors.c_1; 
        std::vector<std::unique_ptr<fixed_d_cluster>> clustering_algos;
        std::vector<std::pair<bool, subtrajectory_cluster_t> > candidate_clusters;
        std::vector<distance_t> gamma;
        
        std::vector<std::set<frequent_pathlet>> max_freq_pathlets_by_distance;
        std::vector<int> num_uncovered_pathlets;
        for (const auto &dist: sq_distances) {
            std::cout<< "Initializing for distance "<<dist << std::endl;
            clustering_algos.emplace_back(std::make_unique<fixed_d_cluster>(trajectory, dist));
            std::cout<< "Finished Initializing cl algos for distance "<<dist << std::endl;
            candidate_clusters.emplace_back();
            std::cout<< "Finished Initializing candidate clusters for distance "<<dist << std::endl;
            gamma.emplace_back();
            std::cout<< "Finished Initializing their ds for distance "<<dist << std::endl;
            max_freq_pathlets_by_distance.emplace_back();
            std::cout<< "Finished Initializing for distance "<<dist << std::endl;
        }
        std::cout<< "Initialized data structures for aided means clustering." << std::endl;
        //Do the sampling and initialization of the vector of frequent sts: we use the same chenoff sample?? For now, yes
        //Generate the sample
        trajectory_t sample = read_trajectory_from_file<space>(samplefilename); //samplefile needs to be passed as input
        std::cout << "Read the sample, now mining frequent pathlets"<< std::endl;
        float frequency_threshold = 0.05;
        for(int i = 0; i<sq_distances.size(); i++){
            const auto &dist = sq_distances.at(i);
            // Compute the frequent pathlets, then swap the values of the freq pathlets vector into the new one before deletion
            float radius = std::sqrt(dist);
            range_search_t rs(sample);
            
            frequent_subtrajectory_algo_t algo(sample, rs, infilename, frequency_threshold, radius); //infilename needs to be passed as parameter to the clustering algorithm
            algo.compute_maximal_frequent_pathlets();
            //std::swap(algo.freq_pathlets, max_freq_pathlets_by_distance.at(i));
            //TODO:this should not be a vector, but a hashmap, must be modified in pathlet tree etc... Unluckily simply changing this ds: Keep both of them for now
            std::set<frequent_pathlet> temp(algo.freq_pathlets.begin(), algo.freq_pathlets.end());

            std::swap(temp, max_freq_pathlets_by_distance.at(i));
            std::cout <<(max_freq_pathlets_by_distance[i].size())<< ","<< algo.freq_pathlets.size()<<std::endl;
            
        }
        std::cout << "Frequent pathlets computed for all distances.\n";
        //TODO: rewrite this because this is horrible
        while (clustering_algos.front()->count_remaining_points() > 0) {
            
            #pragma omp parallel for
            for (size_t i = 0; i < clustering_algos.size(); ++i) {
                
                auto &[success, output_cluster] = candidate_clusters[i];
                
                if(!max_freq_pathlets_by_distance[i].empty()){
                    //extract the frequent pathlet as center
                    success = true;
                    distance_t temp_gamma = std::numeric_limits<float>::max();
                    subtrajectory_cluster_t candidate_cluster;
                    int counter = 0;
                    for(auto& cluster_center : max_freq_pathlets_by_distance[i]){
                        
                        counter++;
                        std::pair<index_t, index_t> offsets = cluster_center.extremes;
                        id_t pathlet_mother = cluster_center.pathlet_mother;
                        subtrajectory_t indexes_in_trajectory = from_freq_pathlet_to_subtrajectory(offsets, pathlet_mother);
                        //Now build a free_space_graph_free_axis, and extract the cluster as they do.
                        cluster_quality_t cl_qual(0,0,sq_distances[i], indexes_in_trajectory.first, indexes_in_trajectory.second);
                        trajectory_t temp = clustering_algos[i]->get_trajectory();
                        subtrajectory_clustering_rightstep<space> clustering_algos_rightstep(temp, std::vector<index_t>(temp.total_size(), index_t{1}), range_search, config);
                        candidate_cluster = clustering_algos_rightstep.cluster_from_candidate(sq_distances[i], cl_qual);
                        distance_t cluster_gamma = clustering_algos[i]->compute_gamma(output_cluster, efficacy_factors);
                        if (cluster_gamma < temp_gamma){

                            temp_gamma = cluster_gamma;
                            std::swap(candidate_cluster, output_cluster);

                        }
                        if(counter >=k){

                            break;
                        }
                    }
                }
                else{
                    
                // If no uncovered pathlet is available, we can use the clustering algorithm to find a cluster.
                
                success = clustering_algos[i]->find_best_cluster_rightstep(output_cluster, config);
                assert(success);

                }
                gamma[i] = clustering_algos[i]->compute_gamma(output_cluster, efficacy_factors);
            }
            // Pick best cluster as in Section 4.3 of Agarwal et. al, 2018.
            size_t best_i = std::max_element(gamma.begin(), gamma.end()) - gamma.begin();
            auto &best_cluster = candidate_clusters[best_i].second;
            k_cluster_detail::prune_inefficient_subtrajectories(trajectory, best_cluster, gamma[best_i], efficacy_factors);
            //TODO: Once the best cluster is selected, delete the unwanted pathlets from the frequent files. I should have a nice way of compputing these
            // On very dense data sets, this can happen when very few points remain, due to floating point inaccuracies.
            // The remaining points are better left unclustered, assuming c1 > 0.
            if (best_cluster.get_subtrajectories().empty()) {
                break;
            }
            std::cout << "Best cluster has "
                << "distance: " << std::sqrt(sq_distances[best_i])
                << ", gamma: " << gamma[best_i]
                << ", vertices: " << best_cluster.number_of_vertices() << "\n";
            
            // Each algo has it's own trajectory and range search object.
            // I think this is required for the parallel for loop.
            for (auto &algo : clustering_algos) {
                best_cluster = best_cluster;
                //std::cout << "Deleting from algo"<< std::endl;
                algo->establish_cluster(best_cluster);
            }
            std::cout << "UNCOVERED POINTS : "<< clustering_algos.front()->count_remaining_points() << std::endl;
            //Find first uncovered frequent pathlet for each distance
            for (size_t i = 0; i < max_freq_pathlets_by_distance.size(); i++){

                if(!max_freq_pathlets_by_distance[i].empty()){
                    bool found = false;
                    for(auto iter = max_freq_pathlets_by_distance[i].begin(); iter != max_freq_pathlets_by_distance[i].end();){
                        bool covered = false;
                        frequent_pathlet p= *(iter);
                        trajectory_t traj = clustering_algos[i]->get_trajectory();
                        std::pair<index_t, index_t> offsets = p.extremes;
                        id_t pathlet_mother = p.pathlet_mother;
                        subtrajectory_t indexes_in_trajectory = from_freq_pathlet_to_subtrajectory(offsets, pathlet_mother);
                        for(index_t h = indexes_in_trajectory.first; h <= indexes_in_trajectory.second; h++){
                            if(traj.is_point_deleted(h)){
                                covered = true;
                                break;
                            }
                        }
                        if(covered){
                            iter = max_freq_pathlets_by_distance[i].erase(iter);
                        }
                        else{
                            ++iter;
                        }
                    }

                    /*
                    for  (auto& p:max_freq_pathlets_by_distance[i]){
                        std::cout << "Iterating over the distance "<< i <<std::endl;
                        bool covered = false;
                        
                        trajectory_t traj = clustering_algos[i]->get_trajectory();
                        std::pair<index_t, index_t> offsets = p.extremes;
                        id_t pathlet_mother = p.pathlet_mother;
                        subtrajectory_t indexes_in_trajectory = from_freq_pathlet_to_subtrajectory(offsets, pathlet_mother);

                        for(index_t h = indexes_in_trajectory.first; h <= indexes_in_trajectory.second; h++){
                            if(traj.is_point_deleted(h)){
                                covered = true;
                                break;
                            }
                        }
                        if(covered){

                            max_freq_pathlets_by_distance[i].erase((p));

                        }
                        if(max_freq_pathlets_by_distance[i].empty()){

                            break;
                        }
                        
                    }
                    
                    */
                    
                }
            }            
        
        }
        // all algorithms store the same clustering.
        clustering_algos.front()->drop_inefficient_clusters_means(efficacy_factors);
        pathlets = clustering_algos.front()->get_clusters();





    }
    void perform_aided_means_clustering_multiple_tries_over_length(std::string infilename, std::string samplefilename) {
        std::cout<< "Performing sample aided clustering..." << std::endl;
        // The distance isn't fixed yet, so we will multiply with it later.
        config.cost_per_pathlet = efficacy_factors.c_2 / efficacy_factors.c_1; 
        std::vector<std::unique_ptr<fixed_d_cluster>> clustering_algos;
        std::vector<std::pair<bool, subtrajectory_cluster_t> > candidate_clusters;
        std::vector<distance_t> gamma;
        
        std::vector<std::set<frequent_pathlet>> max_freq_pathlets_by_distance;
        std::vector<int> num_uncovered_pathlets;
        for (const auto &dist: sq_distances) {
            std::cout<< "Initializing for distance "<<dist << std::endl;
            clustering_algos.emplace_back(std::make_unique<fixed_d_cluster>(trajectory, dist));
            std::cout<< "Finished Initializing cl algos for distance "<<dist << std::endl;
            candidate_clusters.emplace_back();
            std::cout<< "Finished Initializing candidate clusters for distance "<<dist << std::endl;
            gamma.emplace_back();
            std::cout<< "Finished Initializing their ds for distance "<<dist << std::endl;
            max_freq_pathlets_by_distance.emplace_back();
            std::cout<< "Finished Initializing for distance "<<dist << std::endl;
        }
        std::cout<< "Initialized data structures for aided means clustering." << std::endl;
        //Do the sampling and initialization of the vector of frequent sts: we use the same chenoff sample?? For now, yes
        //Generate the sample
        trajectory_t sample = read_trajectory_from_file<space>(samplefilename); //samplefile needs to be passed as input
        std::cout << "Read the sample, now mining frequent pathlets"<< std::endl;
        float frequency_threshold = 0.1;
        for(int i = 0; i<sq_distances.size(); i++){
            const auto &dist = sq_distances.at(i);
            // Compute the frequent pathlets, then swap the values of the freq pathlets vector into the new one before deletion
            float radius = std::sqrt(dist);
            range_search_t rs(sample);
            
            frequent_subtrajectory_algo_t algo(sample, rs, infilename, frequency_threshold, radius); //infilename needs to be passed as parameter to the clustering algorithm
            algo.compute_maximal_frequent_pathlets();
            //std::swap(algo.freq_pathlets, max_freq_pathlets_by_distance.at(i));
            //TODO:this should not be a vector, but a hashmap, must be modified in pathlet tree etc... Unluckily simply changing this ds: Keep both of them for now
            std::set<frequent_pathlet> temp(algo.freq_pathlets.begin(), algo.freq_pathlets.end());

            std::swap(temp, max_freq_pathlets_by_distance.at(i));
            std::cout <<(max_freq_pathlets_by_distance[i].size())<< ","<< algo.freq_pathlets.size()<<std::endl;
            
        }
        std::cout << "Frequent pathlets computed for all distances.\n";
        //TODO: rewrite this because this is horrible
        while (clustering_algos.front()->count_remaining_points() > 0) {
            
            #pragma omp parallel for
            for (size_t i = 0; i < clustering_algos.size(); ++i) {
                
                auto &[success, output_cluster] = candidate_clusters[i];
                
                if(!max_freq_pathlets_by_distance[i].empty()){
                    //extract the frequent pathlet as center
                    success = true;
                    distance_t temp_gamma = std::numeric_limits<float>::max();
                    subtrajectory_cluster_t candidate_cluster;
                    frequent_pathlet max_p = *(max_freq_pathlets_by_distance[i].begin());
                    int max_length = 1 + max_p.extremes.second - max_p.extremes.first;
                    for(auto& cluster_center : max_freq_pathlets_by_distance[i]){
                        
                        std::pair<index_t, index_t> offsets = cluster_center.extremes;
                        id_t pathlet_mother = cluster_center.pathlet_mother;
                        if( offsets.second - offsets.first +1 < max_length){
                            break;
                        }
                        
                        subtrajectory_t indexes_in_trajectory = from_freq_pathlet_to_subtrajectory(offsets, pathlet_mother);
                        //Now build a free_space_graph_free_axis, and extract the cluster as they do.
                        cluster_quality_t cl_qual(0,0,sq_distances[i], indexes_in_trajectory.first, indexes_in_trajectory.second);
                        trajectory_t temp = clustering_algos[i]->get_trajectory();
                        
                        subtrajectory_clustering_rightstep<space> clustering_algos_rightstep(temp, std::vector<index_t>(temp.total_size(), index_t{1}), range_search, config);
                        candidate_cluster = clustering_algos_rightstep.cluster_from_candidate(sq_distances[i], cl_qual);
                        distance_t cluster_gamma = clustering_algos[i]->compute_gamma(output_cluster, efficacy_factors);
                        if (cluster_gamma < temp_gamma){

                            temp_gamma = cluster_gamma;
                            std::swap(candidate_cluster, output_cluster);

                        }
                        
                    }
                }
                else{
                    
                // If no uncovered pathlet is available, we can use the clustering algorithm to find a cluster.
                
                success = clustering_algos[i]->find_best_cluster_rightstep(output_cluster, config);
                assert(success);

                }
                gamma[i] = clustering_algos[i]->compute_gamma(output_cluster, efficacy_factors);
            }
            // Pick best cluster as in Section 4.3 of Agarwal et. al, 2018.
            size_t best_i = std::max_element(gamma.begin(), gamma.end()) - gamma.begin();
            auto &best_cluster = candidate_clusters[best_i].second;
            k_cluster_detail::prune_inefficient_subtrajectories(trajectory, best_cluster, gamma[best_i], efficacy_factors);
            //TODO: Once the best cluster is selected, delete the unwanted pathlets from the frequent files. I should have a nice way of compputing these
            // On very dense data sets, this can happen when very few points remain, due to floating point inaccuracies.
            // The remaining points are better left unclustered, assuming c1 > 0.
            if (best_cluster.get_subtrajectories().empty()) {
                break;
            }
            std::cout << "Best cluster has "
                << "distance: " << std::sqrt(sq_distances[best_i])
                << ", gamma: " << gamma[best_i]
                << ", vertices: " << best_cluster.number_of_vertices() << "\n";
            
            // Each algo has it's own trajectory and range search object.
            // I think this is required for the parallel for loop.
            for (auto &algo : clustering_algos) {
                best_cluster = best_cluster;
                //std::cout << "Deleting from algo"<< std::endl;
                algo->establish_cluster(best_cluster);
            }
            std::cout << "UNCOVERED POINTS : "<< clustering_algos.front()->count_remaining_points() << std::endl;
            //Find first uncovered frequent pathlet for each distance
            for (size_t i = 0; i < max_freq_pathlets_by_distance.size(); i++){

                if(!max_freq_pathlets_by_distance[i].empty()){
                    bool found = false;
                    for(auto iter = max_freq_pathlets_by_distance[i].begin(); iter != max_freq_pathlets_by_distance[i].end();){
                        bool covered = false;
                        frequent_pathlet p= *(iter);
                        trajectory_t traj = clustering_algos[i]->get_trajectory();
                        std::pair<index_t, index_t> offsets = p.extremes;
                        id_t pathlet_mother = p.pathlet_mother;
                        subtrajectory_t indexes_in_trajectory = from_freq_pathlet_to_subtrajectory(offsets, pathlet_mother);
                        for(index_t h = indexes_in_trajectory.first; h <= indexes_in_trajectory.second; h++){
                            if(traj.is_point_deleted(h)){
                                covered = true;
                                break;
                            }
                        }
                        if(covered){
                            iter = max_freq_pathlets_by_distance[i].erase(iter);
                        }
                        else{
                            ++iter;
                        }
                    }

                    /*
                    for  (auto& p:max_freq_pathlets_by_distance[i]){
                        std::cout << "Iterating over the distance "<< i <<std::endl;
                        bool covered = false;
                        
                        trajectory_t traj = clustering_algos[i]->get_trajectory();
                        std::pair<index_t, index_t> offsets = p.extremes;
                        id_t pathlet_mother = p.pathlet_mother;
                        subtrajectory_t indexes_in_trajectory = from_freq_pathlet_to_subtrajectory(offsets, pathlet_mother);

                        for(index_t h = indexes_in_trajectory.first; h <= indexes_in_trajectory.second; h++){
                            if(traj.is_point_deleted(h)){
                                covered = true;
                                break;
                            }
                        }
                        if(covered){

                            max_freq_pathlets_by_distance[i].erase((p));

                        }
                        if(max_freq_pathlets_by_distance[i].empty()){

                            break;
                        }
                        
                    }
                    
                    */
                    
                }
            }            
        
        }
        // all algorithms store the same clustering.
        clustering_algos.front()->drop_inefficient_clusters_means(efficacy_factors);
        pathlets = clustering_algos.front()->get_clusters();



    }
    distance_t compute_means_efficacy() {
        return k_cluster_detail::compute_efficacy_means(trajectory, pathlets, efficacy_factors);
    }
    distance_t compute_center_efficacy() {
        return k_cluster_detail::compute_efficacy_center(trajectory, pathlets, efficacy_factors);
    }

    void print_pathlets() {
        k_cluster_detail::print_pathlets(pathlets);
    }

    // this->trajectory might be modified, e.g. because we deleted the covered points.
    // Thus, we need to pass the original trajectory as an argument.
    void print_clustering_spaced(trajectory_t &trajectory, std::ostream &stream) {
        k_cluster_detail::print_clustering_spaced(trajectory, pathlets, stream);
    }

    void print_clustering_csv(std::ostream &stream) {
        k_cluster_detail::print_clustering_csv(pathlets, stream);
    }

    const std::vector<subtrajectory_cluster_t>& get_clusters() const {
        return pathlets;
    }

private:
    trajectory_t trajectory;
    range_search_t range_search;
    std::vector<distance_t> sq_distances;
    std::vector<subtrajectory_cluster_t> pathlets;
    // Maximum distance at which to search for clusters
    distance_t distance_limit;
    // Weights for computing the efficacy
    const efficacy_factor_t efficacy_factors;
    rightstep_config config;
    
    void initialize_distance_limits(distance_t min_distance, distance_t max_distance){
        if(min_distance < 0 || max_distance < 0){
            const auto [trajectory_min_sq, trajectory_max_sq] = k_cluster_detail::compute_min_max_sq_distance(trajectory, range_search);
            if(min_distance < 0) {
                min_distance = std::sqrt(trajectory_min_sq);
                std::cout << "minimum distance: " << min_distance << "\n";
            }
            if(max_distance < 0) {
                max_distance = std::sqrt(trajectory_max_sq);
                std::cout << "maximum distance: " << max_distance << "\n";
            }
            // avoid weirdness if only one of min_distance, max_distance was computed and is now incompatible with the other (specified) one.
            assert(min_distance <= max_distance);
        }
        distance_limit = max_distance;
        k_cluster_detail::initialize_sq_distances(min_distance*min_distance, max_distance*max_distance, sq_distances);
    }
    //TODO : break ties by frequency????
    static inline bool pathlet_sorter_by_length(frequent_pathlet p1, frequent_pathlet p2){

        if(p1.extremes.second - p1.extremes.first < p2.extremes.second - p2.extremes.first ){

            return true;
        }
        
        return false;

    }
    //TODO. break ties by length????
    static inline bool pathlet_sorter_by_frequency(frequent_pathlet p1, frequent_pathlet p2){

        if(p1.frequency < p2.frequency){

            return true;
        }
        

        return false;

    }
    static inline bool pathlet_sorter_by_reverse_frequency(frequent_pathlet p1, frequent_pathlet p2){

        if(p1.frequency > p2.frequency){

            return true;
        }
        

        return false;

    }
    subtrajectory_t from_freq_pathlet_to_subtrajectory(std::pair<index_t,index_t> offsets, id_t pathlet_mother){

        index_t initial_point = trajectory.get_first_point_in_trajectory(pathlet_mother);
        //std::cout << "The pathelt mother starts at " << initial_point << std::endl;
        subtrajectory_t st = {initial_point+offsets.first, initial_point+ offsets.second};
        //std::cout << "The subtrajectory is from " << st.first << " to " << st.second << std::endl;
        return st;
    }

   
    subtrajectory_cluster_t extract_subtrajectory_cluster_from_pathlet(const subtrajectory_t columns, const distance_t &distance_max, std::unique_ptr<fixed_d_cluster>& clustering_algo ){

        
        free_space_graph_incremental_t free_space(columns.first, config.prefer_small_subtrajectories, config.cost_per_pathlet);
        for(index_t column =columns.first; column <= columns.second; ++column){
            if(column != columns.first) free_space.new_column();
            int zeroes = 0;
            if(range_search.search(column, distance_max).empty()) std::cout << "Empty range search for column " << column << "\n";
            for(const auto idx : range_search.search(column, distance_max)){
                // In rare cases, points fail to be deleted from the KD-tree.
                
                if (trajectory.is_point_deleted(idx)) continue;
                free_space.add_zero(my_spaced_index(idx));

                
            }
           
        }
        subtrajectory_cluster_t cluster;
        free_space.query_subtrajectories(cluster, my_spaced_index(columns.first), my_spaced_index(columns.second));
        std::cout << "    Actual cluster  without undoing the index" << cluster.number_of_vertices() << " points on " << cluster.size() << " trajectories\n";
        cluster.set_reference_trajectory({columns.first, columns.second});
        my_undo_index_spacing(cluster);
        std::cout << "Cluster center should have length "<< - columns.first + columns.second + 1 <<"\n";
        std::cout << "    Actual cluster has " << cluster.number_of_vertices() << " points on " << cluster.size() << " trajectories\n";
        return cluster;
    }

    
    
    // spacing out row indices ensures that subtrajectories stay within a trajectory.
    index_t my_spaced_index(index_t index){
        return index + trajectory.get_id_at(index);
    }

    // Transform: Subtrajectory (l + trajectory_id[l], r + trajectory_id[r]) to (l, r).
    // Assumes trajectories are non-overlapping and sorted in decreasing order
    void my_undo_index_spacing(subtrajectory_cluster_t &cluster){
        index_t orig_index = trajectory.total_size() - 1;
        // Only works if index is non-increasing in subsequent calls.
        auto unspace_index = [&](index_t& index){
            while (trajectory.is_point_deleted(orig_index) || my_spaced_index(orig_index) != index) {
                assert(orig_index >= 0);
                --orig_index;
            }
            index = orig_index;
        };
        for (index_t i = 0; i < cluster.size(); ++i){
            unspace_index(cluster[i].second);
            unspace_index(cluster[i].first);
        }
    }
};

} // End of namespace `frechet`

