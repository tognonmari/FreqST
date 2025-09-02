#pragma once

#include <cmath>
#include <iostream>

#include "frechet_distance.h"
#include "free_space_graph.h"
#include "kdtree_range_search.h"
#include "metric_space.h"
#include "profiling.h"
#include "subtrajectory_cluster.h"
#include "subtrajectory_routine_rightstep.h"
#include "trajectory.h"

namespace frechet {

namespace internal {

template<metric_space m_space>
class curve_simplification{

public:
    using space = m_space;
    using range_search_t = kd_tree_range_search<space>;

    using trajectory_t = trajectory_collection<space>;
    using index_t = trajectory_t::index_t;
    using subtrajectory_t = trajectory_t::subtrajectory_t;
    using distance_function_t = space::distance_function_t;
    using distance_t = distance_function_t::distance_t;
    using subtrajectory_cluster_t = subtrajectory_cluster<space>;
    using rightstep_algo_t = subtrajectory_clustering_rightstep<space>;
    using cluster_summary_t  = rightstep_algo_t::cluster_summary_t;

    using weights_t = std::vector<index_t>;

    curve_simplification(const trajectory_t &trajectory, const distance_t &sq_distance, const double simplification_factor) :
        original_trajectory(trajectory), original_distance(sq_distance), simplification_factor(simplification_factor),
        simplified_trajectory(), original_leftmost_index(), point_weights_(), range_search_(compute_simplification()) {

        }

    std::optional<cluster_summary_t> unsimplify(std::optional<cluster_summary_t> cluster) const {
        if (cluster) {
            std::cerr << "    Unsimplify [" << cluster->left_column << "," << cluster->right_column << "] to ";
            cluster->left_column = original_leftmost_index[cluster->left_column];
            cluster->right_column = original_leftmost_index[cluster->right_column] + point_weights_[cluster->right_column] - 1;
            std::cerr << "[" << cluster->left_column << "," << cluster->right_column << "]\n";
        }
        return cluster;
    }
    subtrajectory_cluster_t unsimplify_whole_cluster(subtrajectory_cluster_t& cluster){

        subtrajectory_cluster_t new_cluster;
        std::optional<cluster_summary_t> cl_qual = this->unsimplify(cluster_summary_t{0,0,0,cluster.get_reference_subtrajectory().first, cluster.get_reference_subtrajectory().second});
        new_cluster.set_reference_trajectory(subtrajectory_t{cl_qual->left_column, cl_qual->right_column});
        for (auto& cl_pathlet: cluster.get_subtrajectories()){

            cl_qual = this->unsimplify(cluster_summary_t{0,0,0,cl_pathlet.first, cl_pathlet.second});
            new_cluster.push_back(subtrajectory_t{cl_qual->left_column, cl_qual->right_column});

        }

        return new_cluster;
    }
    const index_t get_original_index(index_t simplified_index) const{

        return original_leftmost_index[simplified_index] + point_weights_[simplified_index] - 1;

    }
    const weights_t& point_weights() const {
        return point_weights_;
    }
    const trajectory_t& trajectory() const {
        return simplified_trajectory;
    }
    range_search_t& range_search() {
        return range_search_;
    }
    double sq_distance() const {
        // Reconstructing the actual cluster on the original curve only works
        // if we reduce the sq_distance on the simplified curve.
        const auto remaining_factor = 1.0 - simplification_factor;
        return original_distance * remaining_factor * remaining_factor;
    }

    trajectory_t& original(){

        return original_trajectory;
    }


    subtrajectory_t from_unsimplified_trajectory_to_simplified_subtrajectory(const subtrajectory_t offsets,  const id_t pathlet_mother){
        /*
        index_t pm_origin = simplified_trajectory.get_first_non_deleted_point_in_trajectory(pathlet_mother);
        subtrajectory_t indexes_in_simplified_trajectory{std::numeric_limits<index_t>::max(), std::numeric_limits<index_t>::max()};
        index_t aggregate_weights = 0;
        for (index_t t=0; t< simplified_trajectory.get_trajectory_size(pathlet_mother); t++){
            aggregate_weights += point_weights_[pm_origin +t];
            assert(simplified_trajectory.get_id_at(pm_origin + t) == pathlet_mother);
            if(indexes_in_simplified_trajectory.first == std::numeric_limits<index_t>::max() && aggregate_weights>= offsets.first ){

                indexes_in_simplified_trajectory.first = pm_origin+t;
            }

            if (indexes_in_simplified_trajectory.second == std::numeric_limits<index_t>::max() && aggregate_weights>= offsets.second){

                assert(indexes_in_simplified_trajectory.first != std::numeric_limits<index_t>::max());
                indexes_in_simplified_trajectory.second = pm_origin +t;
            } 
            std::cout << "Right extreme is "<< offsets.second << "and we are at "<< aggregate_weights <<std::endl;
            
            if(indexes_in_simplified_trajectory.second != std::numeric_limits<index_t>::max() && indexes_in_simplified_trajectory.first != std::numeric_limits<index_t>::max()){
                break;
            }
        }
        assert((indexes_in_simplified_trajectory.second != std::numeric_limits<index_t>::max() && indexes_in_simplified_trajectory.first != std::numeric_limits<index_t>::max()));
        
        
        return indexes_in_simplified_trajectory;
        
        */
        index_t lu = offsets.first;
        index_t ru = offsets.second;
        subtrajectory_t indexes_in_simplified_trajectory{std::numeric_limits<index_t>::max(), std::numeric_limits<index_t>::max()};
        if(simplified_trajectory.get_first_point_in_trajectory(pathlet_mother)+1 > simplified_trajectory.total_size()){

            return subtrajectory_t{simplified_trajectory.get_first_point_in_trajectory(pathlet_mother),simplified_trajectory.get_first_point_in_trajectory(pathlet_mother)};

        }
        for (index_t k =simplified_trajectory.get_first_point_in_trajectory(pathlet_mother); original_leftmost_index[k]<=ru; k++){
            
            if(indexes_in_simplified_trajectory.first == std::numeric_limits<index_t>::max() && original_leftmost_index[k]<=lu && original_leftmost_index[k+1]>lu){

                indexes_in_simplified_trajectory.first = k;

            }
            if(indexes_in_simplified_trajectory.second == std::numeric_limits<index_t>::max() && original_leftmost_index[k]<=ru && original_leftmost_index[k+1]>ru){

                indexes_in_simplified_trajectory.second = k;

            }
            if(indexes_in_simplified_trajectory.second != std::numeric_limits<index_t>::max() && indexes_in_simplified_trajectory.first != std::numeric_limits<index_t>::max()){
                break;
            }


        }
        assert((indexes_in_simplified_trajectory.second != std::numeric_limits<index_t>::max() && indexes_in_simplified_trajectory.first != std::numeric_limits<index_t>::max()));
        
        
        return indexes_in_simplified_trajectory;

    }
private:
    const trajectory_t& compute_simplification(){
        size_t l = 0;
        assert(simplified_trajectory.get_point_vector().empty());
        simplified_index.clear();
        for(size_t i=0; i <= original_trajectory.total_size(); ++i){
            if (i == original_trajectory.total_size() || !are_points_compatible(l, i)){
                // Compress the points in [l, i) to a single one.
                // Keep deleted points to remember that the trajectory was split into multiple pieces.
                id_t id = original_trajectory.get_id_at(l);
                simplified_trajectory.push_back(original_trajectory[l], id, true);
                //simplified_trajectory.push_back(original_trajectory[l], original_trajectory.get_id_at(l), true);
                point_weights_.push_back(i - l);
                original_leftmost_index.push_back(l);

                l = i;
            }
            if(i< original_trajectory.total_size()){

                simplified_index.push_back(simplified_trajectory.total_size());
            }
        }
        assert(simplified_index.size() == original_trajectory.total_size());
        return simplified_trajectory;
    }

    bool are_points_compatible(index_t i, index_t j){
        if (original_trajectory.get_id_at(i) != original_trajectory.get_id_at(j)) {
            return false;
        }
        auto const sq_distance = distance_function_t{}(original_trajectory[i], original_trajectory[j]);
        return sq_distance <= simplification_sq_distance();
    }
    double simplification_sq_distance() const {
        // If we simplify both curves by `distance * factor/2`, the maximum total error is `distance * factor`.
        const auto factor_halved = simplification_factor / 2;
        return original_distance * factor_halved * factor_halved;
    }

    const trajectory_t &original_trajectory;
    const distance_t original_distance;
    const double simplification_factor;
    weights_t simplified_index;
    trajectory_t simplified_trajectory;
    weights_t original_leftmost_index, point_weights_;

    range_search_t range_search_;
};

} // namespace internal

} // namespace frechet
