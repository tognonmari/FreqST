#include <CGAL/Dimension.h>
#include <CGAL/Simple_cartesian.h>
#include <iostream>
#include "CLI11.hpp"
#include <chrono>
#include "free_space_graph.h"
#include "io.h"
#include "freq_st_algo.h"
#include "metric_space.h"
#include "trajectory.h"
#include "utility"
#include "kdtree_range_search.h"
#include "free_space_graph_free_axis.h"
#include "vc_dim_extractor.h"

using namespace frechet;
namespace fs = std::filesystem;
namespace chrono = std::chrono;
 
using space = CGAL_metric_space<CGAL::Simple_cartesian<double>, CGAL::Dimension_tag<2>>;
using free_space_graph_t =  free_space_graph_free_axis<space>;
using frequent_subtrajectory_algo_t = frequent_subtrajectory_algo<space>;
using trajectory_t = trajectory_collection<space>;
using distance_function_t = space::distance_function_t;
using distance_t = distance_function_t::distance_t; 
using range_search_t = kd_tree_range_search<space>;

int main(int argc, char** argv){

    //Step 1: initialize config variables
    freq_subtrajectory_algo_output_config output_config;
    distance_t radius;  
    float frequency_threshold, epsilon, delta;
    int minimum_length = 1;
    std::string infilename, outfilename, pathlet_file_name;

    //Step 2: parse the input parameters
    
    CLI::App app{"Frequent Subtrajectory Extraction"};
    app.add_option("-r,--radius",
                   radius,
                   "The maximum radius for the Frechet distance. If not provided, will use the ")
                   ->required();  
    app.add_option("input",
                   infilename,
                   "The file with the dataset.")
                    ->required();

    CLI11_PARSE(app, argc, argv);   

    trajectory_t dataset = read_trajectory_from_file<space>(infilename);

    std::cout << "Num trajectories in the sample is "<< dataset.num_trajectories_not_consecutive()<<std::endl;
    
    range_search_t rs(dataset);
    
    vc_dim_extractor<space> extractor(dataset, rs, infilename, radius);

    extractor.compute_set_of_supports();
    std::cout << "Finished computing set of supports \n";
    extractor.convert_supports_to_bitsets();
    std::cout << "Finished conversion to bitsets \n";
    int vc_dim = extractor.compute_exact_vc_dimension();
    std::cout<< "Exact vc dim for dataset at distance "<< radius << " is: "<< vc_dim<< std::endl;

    return 0;
}

