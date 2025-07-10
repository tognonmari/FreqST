#include <CGAL/Dimension.h>
#include <CGAL/Simple_cartesian.h>

#include <omp.h>
#include <chrono>
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
#include <chrono>
using namespace frechet;
namespace chrono = std::chrono;
using space = CGAL_metric_space<CGAL::Simple_cartesian<double>, CGAL::Dimension_tag<2>>;


using free_space_graph_t =  free_space_graph_free_axis<space>;
using frequent_subtrajectory_algo_t = frequent_subtrajectory_algo<space>;
using trajectory_t = trajectory_collection<space>;
using validation_t = validation<space>;
using index_t = trajectory_t::index_t;
using distance_function_t = space::distance_function_t;
using distance_t = distance_function_t::distance_t;
using range_search_t = kd_tree_range_search<space>;
namespace chrono = std::chrono;
enum class sampling_mode {
    chernoff = 0,
    vc_dim = 1
};

int main(int argc, char** argv){

    //Step 1: initialize config variables
    
    float frequency_threshold, epsilon, delta;
    int minimum_length = 1;
    int seed = 0;
    distance_t radius;
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
    app.add_option("-r,--radius",
                    radius,
                    "The radius for the vc dimension.")
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


    std::cout << "your seed is : " << seed <<std::endl;
    
    trajectory_t dataset = read_trajectory_from_file<space>(infilename);
    freq_subtrajectory_sampler<space> sampler(dataset, epsilon, delta, radius, minimum_length, seed );
    //sampler.generate_fixed_size_sample(23000);
    //sampler.dump_sample_to_file(std::format("./berlin/fixed_size_sample_23000.txt"));
    //sampler.generate_chernoff_sample();
    //sampler.dump_sample_to_file(std::format("{}/chernoff_{}_{}_{}.txt", outfiledir, epsilon, delta, seed)); 
    //sampler.generate_vc_sample();
    auto start = chrono::high_resolution_clock::now();
    sampler.generate_rough_vc_sample();
    auto end = chrono::high_resolution_clock::now();
    auto duration = duration_cast<chrono::milliseconds>(end - start);
    std::cout << "TIME: "<<duration.count() << std::endl;
    //sampler.dump_sample_to_file(std::format("{}/vc_{}_{}_{}_{}.txt", outfiledir, epsilon, delta, seed, radius));
    //frequent_subtrajectory_algo_t algo(dataset, infilename, 0.4, 50); 
    //algo.populate_range_search_tree_with_sample_points();

    //algo.compute_all_frequent_pathlets(); 
    //algo.dump_collected_pathlets_to_file("./prova.txt");
    return 0;
}
