#pragma once
#include <queue>
#include <array>
#include <cassert>
#include <limits>
#include <list>
#include <unordered_map>
#include <ostream>
#include <utility>
#include <vector>

#include "ankerl/unordered_dense.h"

#include "profiling.h"
#include "recycling_object_pool.h"
#include "subtrajectory_cluster.h"
#include "trajectory.h"

namespace frechet {

template<metric_space m_space>
class free_space_graph_free_axis{

    public:
    using space = m_space;
    using point_t = m_space::point_t;
    using trajectory_t = trajectory_collection<space>;
    using index_t = trajectory_t::index_t;
    using subtrajectory_t = trajectory_t::subtrajectory_t;
    using subtrajectory_cluster_t = subtrajectory_cluster<space>;
    using id_t = trajectory_t::id_t;

public:
    using row_index_t = index_t;

private:
    struct vertex {
        using label_t = index_t;
        static constexpr label_t no_edge_label = std::numeric_limits<label_t>::max();

        vertex(const row_index_t &row_index) : row_index(row_index) {}

        row_index_t row_index;

        vertex* up = nullptr;

        vertex* left = nullptr;
        vertex* below_left = nullptr;
        vertex* below = nullptr;

        label_t label_left = no_edge_label;
        label_t label_below_left = no_edge_label;
        label_t label_below = no_edge_label;
        label_t min_label = no_edge_label;
    };

public:
    free_space_graph_free_axis(const index_t &init_right_column) : left_column(init_right_column),
                                                                right_column(init_right_column) {
        lowest_vertex_per_column.push_back(nullptr);
    }

    void new_column() {
        new_column(right_column + 1);
    }

    void new_column(index_t new_right_column) {
        if (!lowest_vertex_per_column.empty() && lowest_vertex_per_column.back() != nullptr) {
            candidate_for_left = lowest_vertex_per_column.back();
            candidate_for_below_left = candidate_for_left->below;
        } else {
            candidate_for_left = candidate_for_below_left = nullptr;
        }
        lowest_vertex_per_column.push_back(nullptr);
        highest_in_last_col = nullptr;
        highest_vertex_per_column.push_back(nullptr);
        right_column = new_right_column;
    }
    //Adds a new point to the fsg and connects it only to vertices with the same trajectory id in the_traj
    void add_zero_with_id_respecting_labels(row_index_t row_idx, const trajectory_t& the_traj){

        assert((highest_in_last_col == nullptr) || (highest_in_last_col->row_index < row_idx));
        auto *new_vertex = vertex_pool.construct(row_idx);
        if (lowest_vertex_per_column.back() == nullptr) {
            lowest_vertex_per_column.back() = new_vertex;
        }
        
        if(right_column == 0){
            
            new_vertex->min_label = 0;

        }
        if (right_column > 0) {
            
            advance_candidate_for_left(row_idx);
            //Add the link to candidat_for_left onlyi if they have the same id 
            if (candidate_for_left != nullptr && candidate_for_left->row_index == row_idx) {
                new_vertex->left = candidate_for_left;
                new_vertex->label_left = candidate_for_left->min_label;
            }
            if (candidate_for_below_left != nullptr && candidate_for_below_left->row_index == row_idx - 1 && the_traj.get_id_at(candidate_for_below_left->row_index)==the_traj.get_id_at(row_idx)) {
                new_vertex->below_left = candidate_for_below_left;
                new_vertex->label_below_left = candidate_for_below_left->min_label;
            }
        }
        if (highest_in_last_col != nullptr) {
            new_vertex->below = highest_in_last_col;
            new_vertex->label_below = (highest_in_last_col->row_index == row_idx - 1 && the_traj.get_id_at(highest_in_last_col->row_index)==the_traj.get_id_at(row_idx)) ? highest_in_last_col->min_label
                                                                                    : vertex::no_edge_label;
            highest_in_last_col->up = new_vertex;
        }
        highest_in_last_col = new_vertex;
        new_vertex->min_label = std::min({new_vertex->label_left,
                                          new_vertex->label_below_left,
                                          new_vertex->label_below,
                                          right_column});


    }

    // Prints the fsg. Must be read from bottom-right to top left. 
    //it is upside down w.r.t. the standard representations of free space graphs 
    std::string to_string(const trajectory_t &sample, subtrajectory_t& chunk){


        //Detect the width and the height of the fsg
        int width = lowest_vertex_per_column.size() -1;
        if (width <= 0){
            return "############## EMPTY FSG ##############\n";
        }
        index_t current_row = std::numeric_limits<unsigned>::max();
        index_t highest_height =0;
        std::vector<vertex*> highest_vertex_per_column;
        id_t current_visiting_trajectory;
        
        
        for (int j =0; j<lowest_vertex_per_column.size(); j++){

            if (lowest_vertex_per_column[j]==nullptr){
                continue;
            }
            
            auto it = lowest_vertex_per_column.at(j);
            
            row_index_t temp_index;
            while(it != nullptr){
                temp_index = it->row_index;
                if(temp_index > highest_height){
                    
                    highest_height = temp_index;
                    

                }
                
                it = it->up;
                
            }

        }
           
        
        //Find Lowest Row 
        std::vector<vertex*> next_to_be_visited;
        next_to_be_visited.reserve(lowest_vertex_per_column.size());
        for (int j=0; j <lowest_vertex_per_column.size(); j++){

            next_to_be_visited.push_back(lowest_vertex_per_column[j]);

            if (lowest_vertex_per_column[j]!=nullptr && lowest_vertex_per_column[j]->row_index<current_row){

                current_row = lowest_vertex_per_column[j]->row_index;
                current_visiting_trajectory = sample.get_id_at(lowest_vertex_per_column[j]->row_index);
                
            }
        }

        std::string s = "################## UPSIDE-DOWN FSG ( paths from bottom-right to upper-left )##################\n";
        //MISSING matching rows all below 
        index_t initial_row = chunk.first;
        current_visiting_trajectory = sample.get_id_at(chunk.first);
        for (unsigned i = chunk.first; i<=chunk.second; i++){
            if(current_visiting_trajectory !=sample.get_id_at(i)){
                current_visiting_trajectory=sample.get_id_at(i);
                s+= "-----------------------------\n";

            }
            s = s +std::string("TID : ") + std::to_string(sample.get_id_at(i))+std::string( " ROW ")+  std::to_string(i)+ std::string(" : ");
            for (int j=0; j< lowest_vertex_per_column.size(); j++){
                assert(next_to_be_visited[j]==nullptr || next_to_be_visited[j]->row_index>=i);
                if(next_to_be_visited[j]==nullptr || next_to_be_visited[j]->row_index>i){
                    s+=std::string("X");
                }
                else{
                    s+=std::string("0");
                    next_to_be_visited[j] =  next_to_be_visited[j]->up;
                }
            }
            
            s+= '\n';
        }
        return s;

    }

    
    // FOR FREQUENT SUBTRAJECTORIES 
    
    std::set<id_t> query_one_pathlet_over_the_sample_with_labels_by_slice(const trajectory_t&sample, const subtrajectory_t &pathlet, int integer_threshold, const subtrajectory_t& slice){
        
        int counter = 0;
        // It starts from the lowest 0 in the rightmost column of the pathlet, columns are indexed just like the points in the pathlet_mother. In the below fsg, if the pathlet is 0-2 it starts from the only zero along the lowest row.

        // 0 0 0 0
        // 0 1 0 1  
        // 1 1 0 1

        vertex* start_vertex = lowest_vertex_per_column.at(pathlet.second); 
        // Departure and Arrival column indexes
        this->left_column = pathlet.first;
        this->right_column  = pathlet.second;
        vertex* end_vertex = nullptr;
        std::set<id_t> matching_ids;
        //If start_vertex == nullptr it means that the right extreme of the pathlet is not close enough to any point in the sample -> we count 0. That column of the fsg is empty.
        if(start_vertex ==nullptr){

            return  matching_ids;
        }


        bool success = false;
        // I need the current trajectory in order to get to know how high i have to traverse to skip to the next one., or, equivalently, if I am querying valid matches.
        id_t current_visiting_trajectory = sample.get_id_at(start_vertex->row_index);
        
        // I need the last trajectory of the sample to know if I am finishing the visit of the fsg or if there is something weird going on.
        id_t last_trajectory_to_be_visited = (sample.get_id_at(slice.second));
        
        //std::cout<<"Started visiting trajectory "<< current_visiting_trajectory<<std::endl;
        //std::cout<<"Last trajectory "<< last_trajectory_to_be_visited<<std::endl;

        auto next_row = start_vertex->row_index;

        while (true) {
            //std::cout<<"Stuck here." <<std::endl;
            
            bool success = false;//find_match_with_pathlet_from_start_vertex_no_queues(sample, start_vertex,current_visiting_trajectory);
            if(start_vertex->min_label <=pathlet.first){
                success = true;
            }
            //std::cout << "Out of the matching function=> my segmentation fault is not there"<< std::endl;
            if(success){
                //std::cout<< "Success"<<std::endl;
                //std::cout<< "found match for the pathlet" << pathlet.first<< " "<< pathlet.second <<"at trajectory "<< current_visiting_trajectory << std::endl;
                counter++;
                matching_ids.insert(current_visiting_trajectory);
                //std::cout<< "Found a match in trajectory "<< current_visiting_trajectory<< "for pathlet "<< pathlet.first << " "<< pathlet.second<< std::endl;
                if(current_visiting_trajectory == last_trajectory_to_be_visited ){
                    
                    break;
                }
                start_vertex = find_next_eligible_vertex_after_success(start_vertex, sample, current_visiting_trajectory); 
                if(start_vertex == nullptr)//there's nothing above me 
                {
                    break;
                }
                current_visiting_trajectory =  sample.get_id_at(start_vertex->row_index);
                //std::cout << "after success i am moving ato traj "<< current_visiting_trajectory<< std::endl;
            }
            else{

                start_vertex = start_vertex->up;
                if(start_vertex == nullptr){

                    break;
                }
                current_visiting_trajectory =  sample.get_id_at(start_vertex->row_index);
            }
        }
        

        

        assert(counter==matching_ids.size());
        
        return matching_ids;


    }

    int query_one_pathlet_over_the_sample_with_labels(const trajectory_t &sample, const subtrajectory_t &pathlet, int integer_threshold){

        int counter = 0;
        // It starts from the lowest 0 in the rightmost column of the pathlet, columns are indexed just like the points in the pathlet_mother. In the below fsg, if the pathlet is 0-2 it starts from the only zero along the lowest row.

        // 0 0 0 0
        // 0 1 0 1  
        // 1 1 0 1

        vertex* start_vertex = lowest_vertex_per_column.at(pathlet.second); 
        // Departure and Arrival column indexes
        this->left_column = pathlet.first;
        this->right_column  = pathlet.second;
        vertex* end_vertex = nullptr;
        //If start_vertex == nullptr it means that the right extreme of the pathlet is not close enough to any point in the sample -> we count 0. That column of the fsg is empty.
        if(start_vertex ==nullptr){

            return counter;
        }


        bool success = false;
        // I need the current trajectory in order to get to know how high i have to traverse to skip to the next one., or, equivalently, if I am querying valid matches.
        id_t current_visiting_trajectory = sample.get_id_at(start_vertex->row_index);
        // I need the last trajectory of the sample to know if I am finishing the visit of the fsg or if there is something weird going on.
        id_t last_trajectory_to_be_visited = sample.get_id_at(sample.total_size()-1);
        
        //std::cout<<"Started visiting trajectory "<< current_visiting_trajectory<<std::endl;
        //std::cout<<"Last trajectory "<< last_trajectory_to_be_visited<<std::endl;

        auto next_row = start_vertex->row_index;

        while (true) {
            //std::cout<<"Stuck here." <<std::endl;
            
            bool success = false;//find_match_with_pathlet_from_start_vertex_no_queues(sample, start_vertex,current_visiting_trajectory);
            if(start_vertex->min_label <=pathlet.first){
                success = true;
            }
            //std::cout << "Out of the matching function=> my segmentation fault is not there"<< std::endl;
            if(success){
                //std::cout<< "Success"<<std::endl;
                //std::cout<< "found match for the pathlet" << pathlet.first<< " "<< pathlet.second <<"at trajectory "<< current_visiting_trajectory << std::endl;
                counter++;
                
                if(current_visiting_trajectory == last_trajectory_to_be_visited ){
                    break;
                }
                start_vertex = find_next_eligible_vertex_after_success(start_vertex, sample, current_visiting_trajectory); 
                if(start_vertex == nullptr)//there's nothing above me 
                {
                    break;
                }
                current_visiting_trajectory =  sample.get_id_at(start_vertex->row_index);
                //std::cout << "after success i am moving ato traj "<< current_visiting_trajectory<< std::endl;
            }
            else{

                start_vertex = start_vertex->up;
                if(start_vertex == nullptr){

                    break;
                }
                current_visiting_trajectory =  sample.get_id_at(start_vertex->row_index);
            }
        }

        return counter;

    }


    //=============================================================
    
    std::vector<vertex*> lowest_vertex_per_column;
    std::vector<vertex*> highest_vertex_per_column;
private:

    lost::recycling_object_pool<vertex> vertex_pool;
    // List of the lowest vertices in each column
    // `.front()` corresponds to the lowest vertex in the `left_column`
    // `.back()` corresponds to the lowest vertex in the `right_column`
       
    index_t left_column = 0;
    index_t right_column = 0;
    vertex* highest_in_last_col = nullptr;
    vertex* candidate_for_below_left = nullptr;
    vertex* candidate_for_left = nullptr;
    
    // Advance the vertices `candidate_for_left` such that it lies at least in the given row.
    // Sets `candidate_for_below_left` to be the next free vertex down from `candidate_for_left`.
    void advance_candidate_for_left(const row_index_t &row) {
        while (candidate_for_left != nullptr && candidate_for_left->row_index < row) {
            candidate_for_below_left = candidate_for_left;
            candidate_for_left = candidate_for_left->up;
        }
    }
    /*
    // Find the next vertex below `current_vertex` whose row index is at most `below_this`,
    // and which lies on a path that ends in a vertex in the column indicated by `left_column`.
    // Returns a pointer to the vertex if such a vertex was found, `nullptr` otherwise.
    vertex* find_eligible_row(vertex* current_vertex, index_t below_this) const {
        auto *v_ptr = current_vertex;
        while (v_ptr != nullptr &&
                (v_ptr->row_index > below_this || v_ptr->min_label > left_column)) {
            v_ptr = v_ptr->below;
        }
        return v_ptr;
    }
    */
    //Finds next vertex to start the search from after a match with current_visiting_trajectory is found. Returns null if there is not another eligible start vertex, ie. if there is not another start vertex above the current one.
    vertex* find_next_eligible_vertex_after_success(const vertex* start_vertex,const trajectory_t& sample, id_t& current_visiting_trajectory){

        
        vertex* next_vertex =start_vertex->up;
        //while I have not finished the fsg and i am still in the current trajectory
        while(next_vertex != nullptr && sample.get_id_at(next_vertex->row_index) == current_visiting_trajectory ){
            
            next_vertex =next_vertex->up;

        }
        
        return next_vertex;
    }
    /**
    // Extract the subtrajectory to `start_vertex`, write it into `output_trajectory` and return true,
    // unless it overlaps the subtrajectory defined by indices `left_column` and `right_column`.
    // If it overlaps this subtrajectory, return false.
    // The value of `output_trajectory` has no meaning if this function returns false.
    bool extract_trajectory(vertex *start_vertex,
                            subtrajectory_t &output_trajectory) const {
        TIME_BEGIN();

        if (left_column <= start_vertex->row_index && start_vertex->row_index <= right_column) {
            // We are already in "reference trajectory territory"
            TIME_END(extract_trajectory); return false;
        }
        auto current_column_idx = right_column;
        output_trajectory.first = start_vertex->row_index;
        output_trajectory.second = start_vertex->row_index;
        while (current_column_idx > left_column) {
            if (start_vertex->label_left <= left_column) {
                current_column_idx--;
                start_vertex = start_vertex->left;
            } else if (start_vertex->label_below_left <= left_column) {
                current_column_idx--;
                output_trajectory.first--;
                start_vertex = start_vertex->below_left;
            } else if (start_vertex->label_below <= left_column) {
                output_trajectory.first--;
                start_vertex = start_vertex->below;
            }
            if (left_column <= output_trajectory.first && output_trajectory.first <= right_column) {
                // We walked into "reference trajectory territory", so we report failure.
                TIME_END(extract_trajectory); return false;
            }
        }
        TIME_END(extract_trajectory); return true;
    }
*/
    // Extract the subtrajectory to `start_vertex`, write it into `output_trajectory` and return.
    // Writes the resulting vertex in the left row to `output_vertex`.
    // Returns `{true, next_row}` if successful, `{false, next_row}` if no valid subtrajectory was found.
    // Here, `next_row` is the row in which the algorithm should continue.
    // Non-valid subtrajectories are those overlapping the reference subtrajectory defined by `left_column` and `right_column`,
    // and those subtrajectories that cross input trajectory boundaries.
    // The value of `output_trajectory` and `output_vertex` has no meaning if this function returns false.
    std::tuple<bool, index_t, vertex*> extract_trajectory_respecting_ids(const trajectory_t &trajectory,
                                                               vertex *start_vertex,
                                                               subtrajectory_t &output_trajectory) const {
        TIME_BEGIN();

        if (left_column <= start_vertex->row_index && start_vertex->row_index <= right_column) {
            // We are already in "reference trajectory territory"
            
            TIME_END(extract_trajectory); return {false, left_column, nullptr};
        }
        auto current_column_idx = right_column;
        output_trajectory.first = start_vertex->row_index;
        output_trajectory.second = start_vertex->row_index;
        while (current_column_idx > left_column) {
            if (start_vertex->label_left <= left_column) {
                current_column_idx--;
                start_vertex = start_vertex->left;
            } else if (start_vertex->label_below_left <= left_column) {
                current_column_idx--;
                output_trajectory.first--;
                start_vertex = start_vertex->below_left;
            } else if (start_vertex->label_below <= left_column) {
                output_trajectory.first--;
                start_vertex = start_vertex->below;
            }
            if (left_column <= output_trajectory.first && output_trajectory.first <= right_column) {
                // We walked into "reference trajectory territory", so we report failure.
                TIME_END(extract_trajectory); return {false, left_column, nullptr};
            }
            if (trajectory.get_id_at(output_trajectory.first) != trajectory.get_id_at(output_trajectory.second)) {
                TIME_END(extract_trajectory); return {false, output_trajectory.second, nullptr};
            }
        }
        TIME_END(extract_trajectory); return {true, output_trajectory.second, start_vertex};
    }
    /*
    // Optimize the starting point of `subtrajectory` by walking down from end_vertex in the leftmost column,
    // without crossing the row of `next_start_vertex` or violating id constraints.
    void optimize_in_left_column(const trajectory_t &trajectory, subtrajectory_t &subtrajectory, vertex *end_vertex, vertex* const next_start_vertex) const {
        auto& row = subtrajectory.first;
        const auto trajectory_id = trajectory.get_id_at(row);
        assert(row == end_vertex->row_index);
        do {
            --row;
            end_vertex = end_vertex->below;
        } while ((!next_start_vertex || row > next_start_vertex->row_index) && end_vertex && end_vertex->row_index == row
            && trajectory.get_id_at(row) == trajectory_id && (row < left_column || row > right_column));
        ++row;
    }
    */


};
}
