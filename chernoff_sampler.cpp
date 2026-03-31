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
    fixed_size = 0,
    chernoff = 1,
    vc = 2,
    rough_vc = 3,
    rough_vc_no_erase = 4
};

int main(int argc, char** argv){

    //Step 1: initialize config variables
    int fixed_size = 1;

    float epsilon = 0.1;
    float delta = 0.05;

    int minimum_length = 1;
    int seed = 0;
    int fixed_size= -1;
    distance_t radius;
    std::string infilename, outfiledir;
    sampling_mode mode;
    double grid_side_wrt_radius = 0.2;
    //Step 2: parse the input parameters
    
    CLI::App app{"Chernoff sampler for a trajectory dataset for frequent subtrajectory mining."};
    app.add_option("-e,--epsilon",
                 epsilon,
                 "The maximum allowed additive error on the frequency.");    
    app.add_option("-d,--delta",
                   delta,
                   "Confidence parameter.");
    app.add_option("-r,--radius",
                    radius,
                    "The radius for the vc dimension.");
    app.add_option("input",
                   infilename,
                   "The trajectory file")
        ->required();
    app.add_option("-c,--cardinality",
                   fixed_size,
                   "Fixed Size for fixed size sample generation.");
    app.add_option("-s, --seed",
                    seed,
                    "Seed for the random generator. Default is 0.");
    app.add_option("output",
                   outfiledir,   
                   "The directory where the sample will be dumped.")
        ->required();
    app.add_option("-m, --mode",
                   mode,
                   "Sampling mode: 0 for Fixed Size, 1 for Chernoff, 2 for VC, 3 for Rough VC Estimate.")
        ->required();
    app.add_option("-g, --grid_side",
                    grid_side_wrt_radius,
                    "The fraction of the radius to be used as grid side for the vc dimension.");                    
    CLI11_PARSE(app, argc, argv);   

    trajectory_t dataset = read_trajectory_from_file<space>(infilename);
    freq_subtrajectory_sampler<space> sampler(dataset, epsilon, delta, radius, minimum_length, seed, grid_side_wrt_radius);
    switch(mode){
        case sampling_mode::fixed_size:{
            sampler.generate_fixed_size_sample(fixed_size);
            sampler.dump_sample_to_file(std::format("{}/fixedsize_{}_{}.txt", outfiledir, fixed_size, seed));
            break;
        }
        case sampling_mode::chernoff:{
            sampler.generate_chernoff_sample();
            sampler.dump_sample_to_file(std::format("{}/chernoff_{}_{}_{}.txt", outfiledir, epsilon, delta, seed));
            break;
        }
        case sampling_mode::vc:{
            auto start = chrono::high_resolution_clock::now();
            sampler.generate_vc_sample();
            auto end = chrono::high_resolution_clock::now();
            auto duration = duration_cast<chrono::milliseconds>(end - start);
            std::cout << "TIME: "<<duration.count() << std::endl;
            
            //sampler.dump_sample_to_file(std::format("{}/vc_{}_{}_{}_{}.txt", outfiledir, epsilon, delta, seed, radius));
            break;
        }
        case sampling_mode::rough_vc:{
            auto start = chrono::high_resolution_clock::now();
            sampler.generate_rough_vc_sample();
            auto end = chrono::high_resolution_clock::now();
            auto duration = duration_cast<chrono::milliseconds>(end - start);
            std::cout << "TIME: "<<duration.count() << std::endl;
            //sampler.dump_sample_to_file(std::format("{}/roughvc_{}_{}_{}_{}.txt", outfiledir, epsilon, delta, seed, radius));
            break;
        }
        case sampling_mode::rough_vc_no_erase:{
            auto start = chrono::high_resolution_clock::now();
            sampler.generate_rough_vc_no_erase_sample();
            auto end = chrono::high_resolution_clock::now();
            auto duration = duration_cast<chrono::milliseconds>(end - start);
            std::cout << "TIME: "<<duration.count() << std::endl;
            //sampler.dump_sample_to_file(std::format("{}/roughvc_{}_{}_{}_{}.txt", outfiledir, epsilon, delta, seed, radius));
            break;
        }
    } 

    //std::cout << "your seed  is : " << seed <<std::endl;
    
    
    //freq_subtrajectory_sampler<space> sampler(dataset, epsilon, delta, radius, minimum_length, seed );
    //sampler.generate_fixed_size_sample(50);
    //sampler.dump_sample_to_file(std::format(std::format("{}/chernoff_{}_{}_{}.txt", outfiledir, epsilon, delta, seed)));
    //sampler.generate_chernoff_sample();
    //sampler.dump_sample_to_file(std::format("{}/chernoff_{}_{}_{}.txt", outfiledir, epsilon, delta, seed)); 
    //sampler.generate_rough_vc_sample();

    
    //auto start = chrono::high_resolution_clock::now();
    //sampler.generate_vc_sample();
    //auto end = chrono::high_resolution_clock::now();
    //auto duration = duration_cast<chrono::milliseconds>(end - start);
    //std::cout << "TIME: "<<duration.count() << std::endl;
    //sampler.dump_sample_to_file(std::format("{}/vc_{}_{}_{}_{}.txt", outfiledir, epsilon, delta, seed, radius));
    //start = chrono::high_resolution_clock::now();
    //sampler.generate_rough_vc_sample();
    //end = chrono::high_resolution_clock::now();
    //duration = duration_cast<chrono::milliseconds>(end - start);
    //std::cout << "TIME: "<<duration.count() << std::endl;
    
    //frequent_subtrajectory_algo_t algo(dataset, infilename, 0.4, 50); 
    //algo.populate_range_search_tree_with_sample_points();

    //algo.compute_all_frequent_pathlets(); 
    //algo.dump_collected_pathlets_to_file("./prova.txt");
    return 0;
}
