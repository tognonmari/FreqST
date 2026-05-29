#pragma once

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <algorithm>
#include <iterator>
#include <type_traits>
#include <vector>
#include <ranges>
#include <CGAL/Dimension.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Classification/Planimetric_grid.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_d.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Bbox_3.h>

#include "metric_space.h"
#include "trajectory.h"
#include "utility.h"
#include "generic_grid.h"
//This class should have the same interface as kd tree range search header file
namespace frechet{

//Class definition 
template<CGAL_metric_space_concept m_space> 
class grid_range_search{

    public:
    
        using space = m_space;
        using trajectory_t = trajectory_collection<space>;
        using index_t = trajectory_t::index_t;
        using grid_t = LowDimensionalGrid<space, index_t>;
    
        
    private: 
        using kernel = space::kernel;
        using point_t = space::point_t;
        using distance_function_t = space::distance_function_t;
        using distance_t = space::distance_function_t::distance_t;
        using property_map_t = trajectory_property_map_adapter<space>; // A map interface to access trajectory points by id
        using search_trait_base = space::search_traits;
        using search_traits = CGAL::Search_traits_adapter<typename property_map_t::key_type, property_map_t, search_trait_base>;
        using grid_range_t = std::vector<index_t>; 
        using point_id_t = index_t;
        using manual_map_t = std::map<point_id_t, point_t>; // I need this map to be able to erase pts out of range
        
    public:
        
        using result_t = std::vector<typename property_map_t::key_type>;
        using extra_informative_result_t = std::map<point_id_t, std::vector<point_id_t>>;
        using intermediate_result_t = std::array<std::map<point_id_t, std::pair<std::vector<point_t>, std::vector<point_id_t>>>, 2>;
        grid_range_search(const trajectory_t& trajectory, distance_t  grid_side) : point_map(trajectory), 
                                                                            grid(grid_side, trajectory[0]), inserting_trajectories_not_pathlets(false) {
            for (index_t i =0; i< trajectory.total_size(); i++){

                grid.insert(i, trajectory[i]);

            } 

        }
        grid_range_search(const trajectory_t& trajectory, distance_t  grid_side, bool inserting_trajectories_not_pathlets) : point_map(trajectory), 
                                                                            grid(grid_side, trajectory[0], inserting_trajectories_not_pathlets), inserting_trajectories_not_pathlets(inserting_trajectories_not_pathlets) {
            for (index_t i =0; i< trajectory.total_size(); i++){

                grid.insert(trajectory.get_id_at(i), trajectory[i]);

            } 

        }

        grid_range_search(const point_t center_point, distance_t grid_side, bool inserting_trajectories_not_pathlets) : point_map(manual_map_t{}),grid(grid_side, center_point), inserting_trajectories_not_pathlets(inserting_trajectories_not_pathlets) {

            manual_map_t new_map; // I need this map to be able to erase pts out of range 
            point_map = new_map;

        }

        grid_range_search(const point_t center_point, distance_t grid_side) : point_map(manual_map_t{}),grid(grid_side, center_point), inserting_trajectories_not_pathlets(false) {

            manual_map_t new_map; // I need this map to be able to erase pts out of range 
            point_map = new_map;

        }

        void insert(point_id_t id, point_t point){
            //If i insert points one by one and i don't want to create a uselessly chunked trajectory
            auto* map_ptr = std::get_if<manual_map_t>(&point_map);
            if (map_ptr == nullptr){
                return; //Cannot insert using the manula map 
            }
            map_ptr ->insert({id, point});
            grid.insert(id, point);

        }

        void insert(point_id_t id, point_t point, point_id_t extra_id_information){
            //If i insert points one by one and i don't want to create a uselessly chunked trajectory
            auto* map_ptr = std::get_if<manual_map_t>(&point_map);
            if (map_ptr == nullptr){
                return; //Cannot insert using the manula map 
            }
            map_ptr ->insert({id, point});
            grid.insert(id, point, extra_id_information);

        }

        extra_informative_result_t search_and_return_associated_ids(point_t point, distance_t squared_distance, id_t tid){
            auto search_distance_unsquared = std::sqrt(squared_distance);
            intermediate_result_t points_in_hypercube;
            //retrieve cell content
            points_in_hypercube = grid.search(point,search_distance_unsquared, tid);
            //delete points out of range
            if(inserting_trajectories_not_pathlets){
                erase_trajectories_out_of_range(points_in_hypercube[1], point, squared_distance);
            }
            else{
                
                erase_points_out_of_range(points_in_hypercube[1], point, squared_distance);

            }
            std::for_each(points_in_hypercube[1].begin(), points_in_hypercube[1].end(), [&](const auto& kv_pair) {points_in_hypercube[0].try_emplace(kv_pair.first);});
            //points_in_hypercube[0].insert(points_in_hypercube[0].end(), points_in_hypercube[1].begin(), points_in_hypercube[1].end());
            extra_informative_result_t result_map;
            for (const auto& item : points_in_hypercube[0]){

                result_map[item.first] = item.second.second;

            }
            //auto keys = points_in_hypercube[0] | std::views::keys;
            //informative_result_t result(keys.begin(), keys.end());
            return result_map;
            //return points_in_hypercube[0];


        }
        extra_informative_result_t search_and_return_associated_ids_no_erase(point_t point, distance_t squared_distance, id_t tid){
            auto search_distance_unsquared = std::sqrt(squared_distance);
            intermediate_result_t points_in_hypercube;
            //retrieve cell content
            points_in_hypercube = grid.search(point,search_distance_unsquared, tid);
            //delete points out of range
            std::for_each(points_in_hypercube[1].begin(), points_in_hypercube[1].end(), [&](const auto& kv_pair) {points_in_hypercube[0].try_emplace(kv_pair.first);});
            //points_in_hypercube[0].insert(points_in_hypercube[0].end(), points_in_hypercube[1].begin(), points_in_hypercube[1].end());
            extra_informative_result_t result_map;
            for (const auto& item : points_in_hypercube[0]){

                result_map[item.first] = item.second.second;

            }
            //auto keys = points_in_hypercube[0] | std::views::keys;
            //informative_result_t result(keys.begin(), keys.end());
            return result_map;
            //return points_in_hypercube[0];

        }
        result_t search(point_t point, distance_t squared_distance, id_t tid){
            
            auto search_distance_unsquared = std::sqrt(squared_distance);
            intermediate_result_t points_in_hypercube;
            //retrieve cell content
            points_in_hypercube = grid.search(point,search_distance_unsquared, tid);
            //delete points out of range
            if(inserting_trajectories_not_pathlets){
                erase_trajectories_out_of_range(points_in_hypercube[1], point, squared_distance);
            }
            else{
                
                erase_points_out_of_range(points_in_hypercube[1], point, squared_distance);

            }
            std::for_each(points_in_hypercube[1].begin(), points_in_hypercube[1].end(), [&](const auto& kv_pair) {points_in_hypercube[0].try_emplace(kv_pair.first);});
            //points_in_hypercube[0].insert(points_in_hypercube[0].end(), points_in_hypercube[1].begin(), points_in_hypercube[1].end());
            auto keys = points_in_hypercube[0] | std::views::keys;
            result_t result(keys.begin(), keys.end());
            return result;
            //return points_in_hypercube[0];
            
        }
        result_t search_no_erase(point_t point, distance_t squared_distance, id_t tid){
            
            auto search_distance_unsquared = std::sqrt(squared_distance);
            intermediate_result_t points_in_hypercube;
            //retrieve cell content
            points_in_hypercube = grid.search(point,search_distance_unsquared, tid);
            //DO NOT delete points out of range
            //erase_points_out_of_range(points_in_hypercube[1], index, squared_distance);
            //points_in_hypercube[0].insert(points_in_hypercube[0].end(), points_in_hypercube[1].begin(), points_in_hypercube[1].end());
            std::for_each(points_in_hypercube[1].begin(), points_in_hypercube[1].end(), [&](const auto& kv_pair) {points_in_hypercube[0].try_emplace(kv_pair.first);});
            //return points_in_hypercube[0].keys();
            auto keys = points_in_hypercube[0] | std::views::keys;
            result_t result(keys.begin(), keys.end());
            return result;
        }
        
        int count_nearby_points_weighted_no_erase(point_t point, distance_t squared_distance, id_t tid){

            auto search_distance_unsquared = std::sqrt(squared_distance);
            int intersecting_cells_points = grid.count(point, search_distance_unsquared, tid);

            return intersecting_cells_points;


        }

        grid_t& get_grid() {

            return this-> grid;
        }

        std::size_t num_elements(){
            auto* map_ptr = std::get_if<manual_map_t>(&point_map);
            if(map_ptr != nullptr){
                return map_ptr->size();
            }
            else{
                return std::get<property_map_t>(point_map).size();
                
            }
        }

    private:
        std::variant<manual_map_t, property_map_t> point_map; //reference to the original trajectory
        grid_t grid;
        const bool inserting_trajectories_not_pathlets;
        //Retrieve point either from map or adapter
        point_t get_point(point_id_t id){

            return std::visit([id](auto& m) { return m[id]; }, point_map);

        }

        //Search results clean up
        void erase_points_out_of_range(std::map<point_id_t, std::pair<std::vector<point_t>, std::vector<point_id_t>>>& map, point_t center, distance_t distance) {
            std::erase_if(map, [this, center, distance](const auto& kv) {
                return distance_function_t()(get_point(kv.first), center) > distance;
            });
        }

        void erase_trajectories_out_of_range(std::map<point_id_t, std::pair<std::vector<point_t>, std::vector<point_id_t>>>& map, point_t center, distance_t distance) {
            for (auto& pair : map){
                std::erase_if(pair.second.first, [this, center, distance, pair](const point_t p) {
                    return distance_function_t()(p, center) > distance;
                });
            }
            
            std::erase_if(map, [](const auto& kv) {
                return kv.second.first.empty();   
            });
        }

        //Bounding Box Methods
        static inline double find_x_min(const trajectory_t& trajectory) {

            double x_min = std::numeric_limits<double>::max();
            for (auto& pt : trajectory.get_point_vector()){
                if (pt.x() < x_min){
                    x_min = pt.x();
                }
            }
            return x_min;
        }
        static inline double find_x_max(const trajectory_t& trajectory) {

            double x_max = std::numeric_limits<double>::min();
            for (auto& pt : trajectory.get_point_vector()){
                if (pt.x() > x_max){
                    x_max = pt.x();
                }
            }
            return x_max;
        }
        static inline double find_y_min(const trajectory_t& trajectory) {

            double x_min = std::numeric_limits<double>::max();
            for (auto& pt : trajectory.get_point_vector()){
                if (pt.y() < x_min){
                    x_min = pt.y();
                }
            }
            return x_min;
        }
        static inline double find_y_max(const trajectory_t& trajectory) {

            double x_max = std::numeric_limits<double>::min();
            for (auto& pt : trajectory.get_point_vector()){
                if (pt.y() > x_max){
                    x_max = pt.y();
                }
            }
            return x_max;
        }
        static inline std::vector<index_t> idx_vector(const trajectory_t& trajectory){

            std::vector<index_t> idx_vec;
            for (index_t i = 0; i<trajectory.total_size(); i++){
                idx_vec.push_back(i);
            }

        }


};
};
