#include "freq_st_algo.h"
#include <bitset>
#include <boost/dynamic_bitset.hpp>
namespace frechet{
template <metric_space space>
class RangeBitset{

    public:
        using id_t = typename frequent_subtrajectory_algo<space>::id_t;
        RangeBitset(int num_trajectories, std::set<id_t>& range_ids): range(num_trajectories){
            
            for (id_t idd : range_ids){
                int bit_to_set = (int) idd;
                range.set(bit_to_set);
            }
        }
        RangeBitset(boost::dynamic_bitset<> bitset){

            range = bitset;

        }

        bool operator<(const RangeBitset<space>& b)const {
            if (this->range.count() == b.range.count()){
                //return lexicographical order
                for (std::size_t i = range.size(); i-- > 0; ) {   // from MSB to LSB
                    if (this->range[i] != b.range[i]){
                        return this->range[i] < b.range[i];
                    }
                }
                return false;
            }
            return this->range.count() < b.range.count();
        }

        bool operator==(const RangeBitset<space>& other) const{

            return this->range == other.range && this->range.size() == other.range.size();
        }

        int count() const{
            return this->range.count();
        }
        // When invoked returns a vector with all the subsets of this with this.count()-1 bits 
        std::set<RangeBitset> subsets_of_size_minus_one() const{
            std::set<RangeBitset> subsets;
            for (int pos = this->range.find_first(); pos != this->range.npos; pos= this->range.find_next(pos)){

                RangeBitset subset = this->range;
                subset.range.reset(pos);
                subsets.insert(subset);


            }
            return subsets;
        }
        //for debugging purposes
        std::string to_string() const {

            std::string s;
            boost::to_string(this->range, s);
            return s;

        }
        std::string to_string_readable() const{
            std::string s  = "";
            for (int pos = this->range.find_first(); pos != this->range.npos; pos= this->range.find_next(pos)){

                s+= std::format("{} ", pos);

            }
            return s;

        }

    private:
        boost::dynamic_bitset<> range;
    
};


template <metric_space space>
class vc_dim_extractor{
    private:
        struct steps_taken{
            bool computed_frequent_patterns = false;
            bool converted_supports_to_bitsets = false;
        };

    public:
        using frequent_subtrajectory_algo_t = frequent_subtrajectory_algo<space>;
        using trajectory_t = trajectory_collection<space>;
        using distance_function_t = space::distance_function_t;
        using distance_t = distance_function_t::distance_t;
        using range_search_t = kd_tree_range_search<space>;
        //takes same arguments as freqstalgo, as it needs to initialize it
        vc_dim_extractor(trajectory_t& sampled_traj, range_search_t& search, std::string dataset_file, distance_t distance_thresh):algo(sampled_traj,search, dataset_file,0.0,distance_thresh,output_config){
            this->dataset = sampled_traj;
        }

        void compute_set_of_supports(){

            if(!steps.computed_frequent_patterns){
                algo.compute_frequent_pathlets_with_trajectory_slicing();
                steps.computed_frequent_patterns = true;
            }

        }

        void convert_supports_to_bitsets(){
            this->compute_set_of_supports();
            if (!steps.converted_supports_to_bitsets){

                for (auto& fp : algo.freq_pathlets){
                    RangeBitset<space> rb(dataset.num_trajectories_not_consecutive()+1, fp.supporting_trajectories); // Adding +1 to support datasets with ids starting both at 0 and at 1
                    set_of_supports.insert(rb);

                }
                std::cout << "SUPPORTS "<< std::endl;
                int current_visiting_size = 1;
                
                for (auto&s : set_of_supports){
                    if (current_visiting_size != s.count()){
                        std::cout<< "SIZE "<< s.count()<<std::endl;
                        current_visiting_size = s.count();
                    }

                    std::cout << s.to_string_readable() <<std::endl;

                }
                
            
                
                steps.converted_supports_to_bitsets = true;
            }

        }
        
        int compute_exact_vc_dimension(){
            assert(steps.converted_supports_to_bitsets && steps.computed_frequent_patterns);
            int max_shattered=  2;
            std::set<RangeBitset<space>> shattered_subsets;
            //initialize the shattered subsets
            for (const auto& s: set_of_supports){
                if (s.count() ==1){
                    shattered_subsets.insert(s);
                }
                else{
                    break;
                }
            }

            int current_visiting_size =2;
            int largest_size_shattered = 1;
            bool found_shattered = false;
            for (const auto& s : set_of_supports){

                //Update the current visiting size
                if(s.count() > current_visiting_size && !found_shattered){
                    break;
                }
                if (s.count()> current_visiting_size){
                    current_visiting_size++;
                    found_shattered = false;
                }
                //Enumerate subsets of s of size current_visiting size -1 and check if they are shattered.
                std::set<RangeBitset<space>> subsets = s.subsets_of_size_minus_one();
                bool s_shattered = true;
                for (auto& smaller : subsets){

                    if( shattered_subsets.find(smaller) == shattered_subsets.end()){
                        s_shattered = false;
                        break;
                    }

                }
                if (s_shattered){
                    found_shattered = true;
                    shattered_subsets.insert(s);
                    if (s.count() >=2){
                        std::cout << "Shattered set "<<s.to_string_readable() << std::endl;
                    }
                    largest_size_shattered = s.count();
                }
            }

            return largest_size_shattered;
        }
            
        
    private:
    const freq_subtrajectory_algo_output_config output_config{
        .maximal = false,
        .keep_matching_ids = true,
        .min_length= 1
    };
    trajectory_t dataset;
    frequent_subtrajectory_algo_t algo;
    std::set<RangeBitset<space>> set_of_supports;
    steps_taken steps;
    
};

}


