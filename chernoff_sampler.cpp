#include <CGAL/Dimension.h>
#include <CGAL/Simple_cartesian.h>

#include <omp.h>

#include "CLI11.hpp"

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
    
    float frequency_threshold, epsilon, delta;
    int minimum_length = 1;
    int seed = 0;
    std::string infilename, outfiledir;
    sampling_mode mode;

    //Step 2: parse the input parameters
    
    CLI::App app{"Chernoff sampler for a trajectory dataset for frequent subtrajectory mining."};
    app.add_option("-e,--epsilon",
                 epsilon,
                 "The maximum allowed additive error on the frequency.")
        ->required();    
    app.add_option("-d,--delta",
                   delta,
                   "Confidence parameter.")
        ->required();
    app.add_option("input",
                   infilename,
                   "The trajectory file")
        ->required();
    app.add_option("-s, --seed",
                    seed,
                    "Seed for the random generator. Default is 0.");
    app.add_option("output",
                   outfiledir,   
                   "The directory where the sample will be dumped.")
        ->required();
  
    CLI11_PARSE(app, argc, argv);   


    
    
    trajectory_t dataset = read_trajectory_from_file<space>(infilename);
    freq_subtrajectory_sampler<space> sampler(dataset, epsilon, delta, 150, minimum_length, seed );
    sampler.generate_chernoff_sample();
    sampler.dump_sample_to_file(std::format("{}/chernoff_{}_{}_{}.txt", outfiledir, epsilon, delta, seed)); 
    sampler.generate_vc_sample();
    //frequent_subtrajectory_algo_t algo(dataset, infilename, 0.4, 50); 
    //algo.populate_range_search_tree_with_sample_points();

    //algo.compute_all_frequent_pathlets(); 
    //algo.dump_collected_pathlets_to_file("./prova.txt");
    return 0;
}
