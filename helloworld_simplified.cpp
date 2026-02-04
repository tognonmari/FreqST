#include <CGAL/Dimension.h>
#include <CGAL/Simple_cartesian.h>
#include <iostream>
#include "CLI11.hpp"
#include <chrono>
#include "free_space_graph.h"
#include "io.h"
#include "freq_st_algo_simplified.h"
#include "metric_space.h"
#include "trajectory.h"
#include "utility"
#include "kdtree_range_search.h"
#include "free_space_graph_free_axis.h"
#include "curve_simplification.h"
using namespace frechet;
namespace fs = std::filesystem;
namespace chrono = std::chrono;

using space = CGAL_metric_space<CGAL::Simple_cartesian<double>, CGAL::Dimension_tag<2>>;
using free_space_graph_t =  free_space_graph_free_axis<space>;
using frequent_subtrajectory_algo_t = frequent_subtrajectory_algo_simplified<space>;
using trajectory_t = trajectory_collection<space>;
using distance_function_t = space::distance_function_t;
using distance_t = distance_function_t::distance_t;
using range_search_t = kd_tree_range_search<space>;

int main(int argc, char** argv){

    //Step 1: initialize config variables

    distance_t radius;
    double frequency_threshold, epsilon, delta;
    int minimum_length = 1;
    std::string infilename, outfilename, pathlet_file_name;
    double curve_simplification_factor = 0.0;
    //Step 2: parse the input parameters
    
    CLI::App app{"Frequent Subtrajectory Extraction"};
    app.add_option("-r,--radius",
                   radius,
                   "The maximum radius for the Frechet distance. If not provided, will use the ")
                   ->required();  
    app.add_option("-f,--f",
                   frequency_threshold,
                   "The fraction of trajectories the subtrajectories must appear in in order to be considered frequent.")
                   ->required();
    app.add_option("input",
                   infilename,
                   "The file with the sample.")
                    ->required();
    app.add_option("pathlets",
                    pathlet_file_name,
                    "The file with the pathlets.")
                    ->required(); 
    app.add_option("outputfile",
                   outfilename,   
                   "The output file, where the frequent pathlets will be dumped. ")
                ->required();
    app.add_option("-s, --s",
                    curve_simplification_factor,
                    "The factor for curve simplification");                    
    CLI11_PARSE(app, argc, argv);   


    
    
    trajectory_t dataset = read_trajectory_from_file<space>(infilename);
    trajectory_t to_unsimplify = read_trajectory_from_file<space>(pathlet_file_name);
    std::cout << "Num trajectories in the sample is "<< dataset.num_trajectories()<<std::endl;
    //frechet::internal::curve_simplification<space> dataset_simplification(dataset,radius * radius, curve_simplification_factor);
                //std::cout <<"Parsed a transaction."<<std::endl;
    range_search_t rs(dataset);
    
    frequent_subtrajectory_algo_t algo(dataset, rs, pathlet_file_name, frequency_threshold, radius, curve_simplification_factor); 
    
    auto start = chrono::high_resolution_clock::now();
    algo.compute_all_frequent_pathlets_with_trajectory_slicing();
    auto stop =  chrono::high_resolution_clock::now();
    
    auto duration = duration_cast<chrono::milliseconds>(stop - start);
    std::cout<< "TIME : "<< duration.count()<< std::endl;
    algo.unsimplify_collected_pathlets(to_unsimplify);
    //algo.compute_all_frequent_pathlets();
   
    
    algo.dump_collected_pathlets_to_file(outfilename);

    return 0;
}
