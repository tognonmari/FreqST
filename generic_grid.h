#pragma once

#include "metric_space.h"
#include "utility.h"
#include <map>

//The equivalent of Kd_tree (employed by class kd_tree_range_search), but for grid_range_search
template<frechet::CGAL_metric_space_concept space, typename point_identifier_t> //types are the metric space (dimension is needed, along with a distance function i believe)
class LowDimensionalGrid{

    public: //exposed types
        using result_t = std::vector<point_identifier_t>;
    private:
        using kernel = space::kernel;
        using point_t = space::point_t;
        using distance_function_t = space::distance_function_t;
        using distance_t = space::distance_function_t::distance_t;
        using cell_content_t =std::vector<point_identifier_t>;
        static constexpr std::size_t dimension = space::dimension::value;
        using vector_t = std::conditional_t<(dimension==2), typename kernel::Vector_2, std::conditional_t<(dimension==3), typename kernel::Vector_3, void>>;


        template<std::size_t D>
        struct CellKey{ 
            
            std::array<long,D> discrete_cell_coordinates;         
            bool operator<(const CellKey& other) const {

                for (std::size_t i =0; i< dimension; i++){

                    if(this->discrete_cell_coordinates[i]<other.discrete_cell_coordinates[i])
                        return true;
                    if(this->discrete_cell_coordinates[i]> other.discrete_cell_coordinates[i]){
                        return false;
                    }
                }
                return false;
            }
             
            bool operator==(const CellKey& other) const {
                for (std::size_t i =0; i< dimension; i++){

                    if(this.discrete_cell_coordinates[i]!=other.discrete_cell_coordinates[i])
                        return false;

                }
                return true;
            }

            std::string to_string(){

                std::string s = "(";

                for (int i=0; i<dimension; i++)
                    s +=std::format(" {} ", discrete_cell_coordinates[i]);

                return s+")";

            }

            
        };
    public: 
        LowDimensionalGrid(distance_t grid_side, const point_t origin_point) : grid_side(grid_side), origin_point(origin_point) {


        }
        //WON'T CHECK FOR ID DuPLICATES
        void insert(point_identifier_t id, const point_t& point){

            CellKey<dimension> key = compute_cell_key(point);

            grid[key].push_back(id);

        }

        result_t search(point_t point, distance_t search_distance_unsquared){

            //Identify which cells I need to visit
            auto bounds = query_cell_delimiters(point, search_distance_unsquared);
            //std::cout << "CELL BOUNDS "<< bounds.first.to_string()<< " AND "<< bounds.second.to_string()<<std::endl;
            
            //Iterate over the visting cells 
            result_t intersecting_cells_points;
            //Define a lambda to pass this as argument
            auto query_cell_func = [this](CellKey<dimension>& cell, result_t& output,point_t& point,distance_t dist){
                this->query_cell(cell, point, dist, output);
            };
            this->get_points_in_ball_intersecting_cells<dimension>(bounds.first, bounds.second, bounds.first, query_cell_func, intersecting_cells_points, point, search_distance_unsquared, 0);
            return intersecting_cells_points;
        }

        std::size_t total_grid_cells(){

            return grid.size();
        }

        std::size_t total_stored_points(){

            int total = 0;
            for (auto it = grid.begin(); it != grid.end(); ++it){
                total += it->second.size();
            }
            return total;
        }


        std::string to_string(){
            std::string s;
            for (auto& pair : grid ){

                auto key = pair.first;
                auto value = pair.second;
                s += std::format("CELL {} LLV {}: {}\n", key.to_string(), print_point<point_t>(get_cell_lower_left_vertex(key)), get_points_id_string(value));

            }
            return s;

        }

    private:
        //class fields
        //The actual grid is a map with key = cell coordinates, value std::vector<point_t>
        std::map<CellKey<dimension>, cell_content_t> grid;
        const point_t origin_point;
        const double grid_side;

        //Verify whether the cell can have any intersection with the ball, then add points if some intersection exists
        void query_cell(CellKey<dimension>& cell, point_t& center, distance_t radius, result_t& output){
            //Assert the cell is not empty (i.e. assert it exists in the map)
            auto iter = grid.find(cell);
            if (iter==grid.end()){
                //std::cout << "Cell "<< cell.to_string()<<" is empty"<< std::endl;
                return;
            }
            //std::cout << "The cell "<<cell.to_string()<<" is not empty and i am indeed adding the points to the output"<< std::endl;
            
            //Assert some intersection exists: get cell vertices' list and check whether any of them is below the radius threshold
            for (point_t vertex : get_cell_vertices(cell)){
                //If some intersection is possible, append the cell content to the output.
                //Also, I add some fuzzyness for corner cases: better safe than sorry
                if (distance_function_t{}(vertex, center)<= radius*radius* 1.00001){
                    //std::cout << "The cell "<<cell.to_string()<<" is not empty and i am indeed adding the points to the output"<< std::endl;
                    output.insert(output.end(),grid[cell].begin(), grid[cell].end());
                    return;
                }

            }
            //std::cout << "Cell "<< cell.to_string() <<" doesn't intersect the range."<<std::endl;
            
        }
        //Dimension sensitive loop unfolding
        template<std::size_t dim, typename query_func>
        void  get_points_in_ball_intersecting_cells(const CellKey<dim>& min_cell, const CellKey<dim>& max_cell, CellKey<dim> current, query_func&& my_function, result_t& output, point_t center, distance_t radius, std::size_t d =0){

            if (d==dim){
                my_function(current, output,center, radius);
                return;
            }
            for (int i= min_cell.discrete_cell_coordinates[d]; i<=max_cell.discrete_cell_coordinates[d]; i++){

                current.discrete_cell_coordinates[d] = i;
                get_points_in_ball_intersecting_cells<dim>(min_cell,max_cell, current, my_function,output, center, radius,  d+1);
            }

        }

        
        std::vector<point_t>  get_cell_vertices(CellKey<dimension>& cell_id){

            point_t lower_left = get_cell_lower_left_vertex(cell_id);
            std::size_t num_vertices = 1<< dimension; //2^dimension
            std::vector<point_t> vertices;
            


            for (std::size_t mask =0; mask<num_vertices; ++mask) {

                std::array<double, dimension> coordinates;

                for (std::size_t d=0; d<dimension; d++) {

                    double offset = (mask & (1 <<d)) ? grid_side : 0.0;
                    
                    coordinates[d] = offset;

                }
                if constexpr( dimension ==2 ){
                    point_t p(coordinates[0]+ lower_left.x(), coordinates[1]+ lower_left.y());
                    vertices.push_back(p);

                }
                else if constexpr( dimension== 3){
                    point_t p(coordinates[0]+ lower_left.x(), coordinates[1]+ lower_left.y(), coordinates[2]+ lower_left.z());
                    vertices.push_back(p);


                }
                
            }
            /*
            
            std::cout << "Here are the Vertices of cell"<< cell_id.to_string()<< std::endl;
            for (auto vertex : vertices){
                std::cout <<print_point(vertex) <<std::endl;
            }
            */
            return vertices;

        }
        //Identifies where the point is stored
        CellKey<dimension> compute_cell_key(const point_t& point){

            CellKey<dimension> key;
            auto difference = point - origin_point;
            for (std::size_t i = 0; i<dimension; i++){

                key.discrete_cell_coordinates[i] = (long) floor(get_coord(difference, i)/grid_side);

            }

            return key;

        }

        point_t get_cell_lower_left_vertex(CellKey<dimension>& cell_id){

            point_t vertex = origin_point;
            std::array<distance_t, dimension> translation;

            for (size_t i =0; i<dimension; i++){

                translation[i] = cell_id.discrete_cell_coordinates[i] * grid_side;

            }
           
            return vertex + make_vector(translation);


        }

        static std::string get_points_id_string(cell_content_t& list){
            std::string s = " ";
            for (auto& idd : list){
                s += std::to_string(idd);
                s += " ";
            }

            return s;

        }

        std::pair<CellKey<dimension>, CellKey<dimension>> query_cell_delimiters(const point_t& point, distance_t radius) {

            std::pair<CellKey<dimension>, CellKey<dimension>> bounds;

            CellKey<dimension> ball_center_cell = compute_cell_key(point);
            //std::cout << "the point is located in cell "<< ball_center_cell.to_string()<< std::endl;
            long offset_cell_coord = (long) ceil(radius/grid_side); //side of the square "minimally" enclosing
            //std::cout << " The integer offset is "<< offset_cell_coord << std::endl;
            
            for (std::size_t i = 0;i< dimension; i++){

                bounds.first.discrete_cell_coordinates[i] = ball_center_cell.discrete_cell_coordinates[i] - offset_cell_coord;
                bounds.second.discrete_cell_coordinates[i] = ball_center_cell.discrete_cell_coordinates[i] + offset_cell_coord;
                
            }
            //std::cout << "The minimum cell is "<< bounds.first.to_string()<< ", The second is "<< bounds.second.to_string()<< std::endl;
            return bounds;

        }

        constexpr vector_t make_vector(std::array<distance_t, dimension> coordinates){

            if  constexpr (dimension==2){
                vector_t v(coordinates[0], coordinates[1]);
                return v;
            }
            else if constexpr (dimension==3){
                vector_t v(coordinates[0],coordinates[1], coordinates[2]);
                return v;
            }
            else{
                static_assert(dimension==2||dimension==3, "Unsupported dimension");
            }
        }


        template <typename PossiblePointType>
        static std::string print_point(const PossiblePointType& point){
            if constexpr (dimension == 2)
                return std::format(" {:.2f} {:.2f} ", point.x(), point.y());
            else if constexpr (dimension == 3)
                return std::format(" {:.2f} {:.2f} {:.2f} ", point.x(), point.y(), point.z());
            else
                return " Could not convert point to string. (High dimensional case)";
        }

        //Helper to support the various point_t: Point_2 and Point_3
        template <typename PossiblePointType>
        static double get_coord(const PossiblePointType& v, std::size_t i ){
            if constexpr (dimension == 2)
                return i==0 ? v.x() : v.y();
            else if constexpr (dimension == 3)
                return i==0 ? v.x() : (i==1 ? v.y() : v.z());
            else
                return v[i];
        }


};
