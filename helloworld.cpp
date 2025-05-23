#include <CGAL/Dimension.h>
#include <CGAL/Simple_cartesian.h>

#include <omp.h>
#include <filesystem>
#include "CLI11.hpp"
#include <chrono>
#include "free_space_graph.h"
#include "io.h"
#include "canonical_pathlets.h"
#include "freq_st_algo.h"
#include "metric_space.h"
#include "trajectory.h"
#include "utility"
#include "validation.h"
#include "free_space_graph_free_axis.h"

using namespace frechet;
namespace fs = std::filesystem;
namespace chrono = std::chrono;
using space = CGAL_metric_space<CGAL::Simple_cartesian<double>, CGAL::Dimension_tag<2>>;

using bbgll_algo = subtrajectory_clustering_bbgll<space>;
using free_space_graph_t =  free_space_graph_free_axis<space>;
using frequent_subtrajectory_algo_t = frequent_subtrajectory_algo<space>;
using trajectory_t = trajectory_collection<space>;
using validation_t = validation<space>;
using index_t = trajectory_t::index_t;
using distance_function_t = space::distance_function_t;
using distance_t = distance_function_t::distance_t;
using range_search_t = kd_tree_range_search<space>;

enum class sampling_mode {
    chernoff = 0,
    vc_dim = 1
};

int main(int argc, char** argv){

    //Step 1: initialize config variables
    distance_t radius;
    float frequency_threshold, epsilon, delta;
    int minimum_length = 1;
    std::string infilename, outfilename, pathlet_file_name;
    sampling_mode mode;
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
  
    CLI11_PARSE(app, argc, argv);   


    
    
    trajectory_t dataset = read_trajectory_from_file<space>(infilename);
    std::cout << "NUm trajectories in the sample is "<< dataset.num_trajectories()<<std::endl;
    
    
    frequent_subtrajectory_algo_t algo(dataset, pathlet_file_name, frequency_threshold, radius); 
    
    auto start = chrono::high_resolution_clock::now();
    algo.compute_all_frequent_pathlets();
    auto stop =  chrono::high_resolution_clock::now();

    auto duration = duration_cast<chrono::seconds>(stop - start);

    std::cout << "TIME : "<< duration.count()<< std::endl;
    algo.dump_collected_pathlets_to_file(outfilename);
    return 0;
}
