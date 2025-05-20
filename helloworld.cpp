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
    distance_t radius;
    float frequency_threshold, epsilon, delta;
    int minimum_length = 1;
    std::string infilename, outfilename;
    sampling_mode mode;

    //Step 2: parse the input parameters
    
    CLI::App app{"Basic Frequent Subtrajectory Extraction via Sampling"};
    app.add_option("-r,--radius",
                   radius,
                   "The maximum radius for the Frechet distance. If not provided, will use the "); //If radius is not provided, it will be set as 
    app.add_option("-f,--f",
                   frequency_threshold,
                   "The fraction of trajectories the subtrajectories must appear in in order to be considered frequent.");
    app.add_option("-l,--length",
                 minimum_length,
                 "The minimum length of a frequent trajectory to be reported in output.");
    app.add_option("-e,--epsilon",
                 epsilon,
                 "The maximum allowed additive error on the frequency.")
        -> default_val(0.2);    
    app.add_option("-d,--delta",
                   delta,
                   "Confidence parameter.")
        ->default_val(0.1);
    app.add_option("-s,--s",
                mode,
                "Sampling bound to use.")
                ->default_val(0);
    app.add_option("input",
                   infilename,
                   "The trajectory file")
        ->required();
        //->check(CLI::ExistingFile);  

    app.add_option("output",
                   outfilename,   
                   "The output file. Leave blank to use stdout.");
  
    CLI11_PARSE(app, argc, argv);
    
    /*
     * The free space for this at distance 0.2 for toydataset looks like
     *
     *    0 1 2 3 4 5 6 7       trajectory id
     *
     * 7  1 1 1 0 1 1 1 0       2
     * 6  1 1 0 1 1 1 0 1       2
     * 5  1 0 1 1 1 0 1 1       1
     * 4  0 1 1 1 0 1 1 1       1
     * 3  1 1 1 0 1 1 1 0       0
     * 2  1 1 0 1 1 1 0 1       0
     * 1  1 0 1 1 1 0 1 1       0
     * 0  0 1 1 1 0 1 1 1       0
     *
     * When querying subtrajectories over reference [0,3],
     * we should only be getting the reference back when respecting ids,
     * even though [4,7] would be at distance < 0.2 when ignoring ids.
     */

    free_space_graph_t free_space{0};

    free_space.add_zero(0);
    free_space.add_zero(4);
    free_space.new_column();
    free_space.add_zero(1);
    free_space.add_zero(5);
    free_space.new_column();
    free_space.add_zero(2);
    free_space.add_zero(6);
    free_space.new_column(); 
    free_space.add_zero(3);

    free_space.add_zero(7);

    

    
    
    trajectory_t dataset = read_trajectory_from_file<space>(infilename);
    //freq_subtrajectory_sampler sampler(dataset, epsilon, delta, radius, minimum_length );
    //sampler.generate_chernoff_sample();
    //sampler.dump_sample_to_file("./prova.txt"); 
    
    frequent_subtrajectory_algo_t algo(dataset, infilename, 0.1, 150); 
    //algo.populate_range_search_tree_with_sample_points();

    algo.compute_all_frequent_pathlets();
    algo.dump_collected_pathlets_to_file("./prova.txt");
    return 0;
}
