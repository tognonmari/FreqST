#pragma once

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <algorithm>
#include <iterator>
#include <type_traits>
#include <vector>

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
    
        
    private: 
        using kernel = space::kernel;
        using point_t = space::point_t;
        using distance_function_t = space::distance_function_t;
        using distance_t = space::distance_function_t::distance_t;
        using property_map_t = trajectory_property_map_adapter<space>; // A map interface to access trajectory points by id
        using search_trait_base = space::search_traits;
        using search_traits = CGAL::Search_traits_adapter<typename property_map_t::key_type, property_map_t, search_trait_base>;
        using grid_range_t = std::vector<index_t>; 
        using grid_t = LowDimensionalGrid<space, index_t>;
        using point_id_t = grid_t::point_id_t;
        using manual_map_t = std::map<point_id_t, point_t>; // I need this map to be able to erase pts out of range
    public:
        
        using result_t = std::vector<typename property_map_t::key_type>;
        
        grid_range_search(const trajectory_t& trajectory, distance_t  grid_side) : point_map(trajectory), 
                                                                            grid(grid_side, trajectory[0]) {
            for (index_t i =0; i< trajectory.total_size(); i++){

                grid.insert(i, trajectory[i]);

            } 

        }

        grid_range_search(const point_t center_point, distance_t grid_side) : point_map(manual_map_t{}),grid(grid_side, center_point){

            manual_map_t new_map; // I need this map to be able to erase pts out of range 
            point_map = new_map;

        }

        void insert(point_id_t id, point_t point){
            auto* map_ptr = std::get_if<manual_map_t>(&point_map);
            if (map_ptr == nullptr){
                return; //Cannot insert using the manula map 
            }
            map_ptr ->insert({id, point});
            grid.insert(id, point);

        }

        result_t search(point_t point, distance_t squared_distance){
            
            auto search_distance_unsquared = std::sqrt(squared_distance);
            std::array<result_t, 2> points_in_hypercube;
            //retrieve cell content
            points_in_hypercube = grid.search(point,search_distance_unsquared);
            //delete points out of range
            erase_points_out_of_range(points_in_hypercube[1], point, squared_distance);
            points_in_hypercube[0].insert(points_in_hypercube[0].end(), points_in_hypercube[1].begin(), points_in_hypercube[1].end());
            return points_in_hypercube[0];
            

        }
        result_t search_no_erase(point_t point, distance_t squared_distance){
            
            auto search_distance_unsquared = std::sqrt(squared_distance);
            std::array<result_t, 2> points_in_hypercube;
            //retrieve cell content
            points_in_hypercube = grid.search(point,search_distance_unsquared);
            //DO NOT delete points out of range
            //erase_points_out_of_range(points_in_hypercube[1], index, squared_distance);
            points_in_hypercube[0].insert(points_in_hypercube[0].end(), points_in_hypercube[1].begin(), points_in_hypercube[1].end());
            return points_in_hypercube[0];
            

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
        //Retrieve point either from map or adapter
        point_t get_point(point_id_t id){

            return std::visit([id](auto& m) { return m[id]; }, point_map);

        }

        //Search results clean up
        void erase_points_out_of_range(result_t &points, point_t center, distance_t distance) {
            std::erase_if(points, [this, center, distance](const index_t p) {
                return distance_function_t()(get_point(p), center) > distance;
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
