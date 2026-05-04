#include "freq_st_algo.h"
#include <bitset>
#include <boost/dynamic_bitset.hpp>
#include "roaring.hh"
#include <omp.h>
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
class RangeSet{

    public:
        using id_t = typename frequent_subtrajectory_algo<space>::id_t;
        RangeSet(int num_trajectories, std::set<id_t> range_ids): range(range_ids){
            
        }
        
        RangeSet(std::set<id_t> range_ids){

            range = range_ids;

        }
        

        bool operator<(const RangeSet<space>& b)const {

            //IF THEY HAVE THE SAME SIZE:

            if (this->range.size() == b.range.size()){
                
                return std::lexicographical_compare(this->range.begin(), this->range.end(), b.range.begin(), b.range.end());
            }
            //IF THEY HAVE DIFFERENT SIZES
            return this->range.size() < b.range.size();
        }

        bool operator==(const RangeBitset<space>& other) const{

            return this->range.size() == other.range.size() && std::equal(this->range.begin(), this->range.end()), other.range.begin();
        }

        int count() const{
            return this->range.size();
        }
        // When invoked returns a vector with all the subsets of this with this.count()-1 bits 
        std::set<RangeSet> subsets_of_size_minus_one() const{
            std::set<RangeSet> subsets;
            for (auto& tid : this->range){
                RangeSet<space> subset(this->range);
                subset.remove(tid);
                subsets.insert(subset);
            }
            return subsets;
        }
        //for debugging purposes
        /*
        
        std::string to_string() const {

            std::string s;
            boost::to_string(this->range, s);
            return s;

        }
        */

        std::string to_string_readable() const{
            std::string s  = "";
            for (auto& tid : this->range){

                s+= std::format("{} ", tid);

            }
            return s;

        }

    private:
        void remove(id_t tid){
            this->range.erase(tid);
        }
        std::set<id_t> range;
    
};
template <metric_space space>
class RangeRoaringBitmap{

        public:
        using id_t = typename frequent_subtrajectory_algo<space>::id_t;

        
        RangeRoaringBitmap(int num_trajectories, std::set<id_t> range_ids): roaring_bitmap_id(std::numeric_limits<id_t>::max()) {
            for (id_t tid : range_ids){

                range.add(tid);
                last_added = tid;
            }
        }

        RangeRoaringBitmap(int num_trajectories, roaring::Roaring range_ids): roaring_bitmap_id(std::numeric_limits<id_t>::max()){
            range= range_ids;
        }

        RangeRoaringBitmap(roaring::Roaring range_ids) : roaring_bitmap_id(std::numeric_limits<id_t>::max()){

            range = range_ids;

        }
        
        RangeRoaringBitmap(const roaring::Roaring& range_ids, id_t roaring_bitmap_id ): roaring_bitmap_id(roaring_bitmap_id){

            range = range_ids;
            //this-> roaring_bitmap_id = roaring_bitmap_id;

        }

        RangeRoaringBitmap(std::set<id_t> range_ids, id_t roaring_bitmap_id){
            for (id_t tid : range_ids){

                range.add(tid);

            }
            this-> roaring_bitmap_id = roaring_bitmap_id;
        }

        RangeRoaringBitmap(std::set<id_t> range_ids) : roaring_bitmap_id(std::numeric_limits<id_t>::max()){
            for (id_t tid : range_ids){

                range.add(tid);
                last_added = tid;

            }
        }
        bool operator<(const RangeRoaringBitmap<space>& b) const {
            const auto ca = range.cardinality();
            const auto cb = b.range.cardinality();

            if (ca != cb)
                return ca < cb;

            const auto amin = range.minimum();   // roaring provides this
            const auto bmin = b.range.minimum();

            if (amin != bmin)
                return amin < bmin;

            return range.maximum() < b.range.maximum();
        }
        /*
        
        bool operator<(const RangeRoaringBitmap<space>& b)const {

            //IF THEY HAVE THE SAME SIZE:
            auto cthis = this->range.cardinality();
            auto cb = b.range.cardinality();
            if ( cthis == cb ){
                
                return std::lexicographical_compare(this->range.begin(), this->range.end(), b.range.begin(), b.range.end());
            }
            //IF THEY HAVE DIFFERENT SIZES
            return cthis < cb;
        }

        */
        bool operator==(const RangeRoaringBitmap<space>& other) const{

            return this->range== other.range; //built in method of roaring class
        }

        int count() const{
            return this->range.cardinality();
        }

        id_t get_range_id() const {

            return this->roaring_bitmap_id;
        }
        // When invoked returns a vector with all the subsets of this with this.count()-1 bits 
        std::set<RangeRoaringBitmap<space>> subsets_of_size_minus_one() const{
            std::set<RangeRoaringBitmap> subsets;
            for (auto tid : this->range){
                RangeRoaringBitmap<space> subset(this->range);
                subset.remove(tid);
                subsets.insert(subset);
            }
            return subsets;
        }
        //for debugging purposes
        /*
        
        std::string to_string() const {

            std::string s;
            boost::to_string(this->range, s);
            return s;

        }
        */

        RangeRoaringBitmap<space> operator|(const RangeRoaringBitmap& other ) const{

            return RangeRoaringBitmap(this->range | other.range);

        }

        RangeRoaringBitmap<space> operator&(const RangeRoaringBitmap& other ){

            return RangeRoaringBitmap(this->range & other.range);

        }

        RangeRoaringBitmap<space> operator&(const RangeRoaringBitmap& other ) const{

            return RangeRoaringBitmap(this->range & other.range);

        }

        std::string to_string_readable() const{
            
            std::string s  = "";
            for (auto& tid : this->range){

                s+= std::format("{} ", tid);

            }
            
            return s;

        }

        std::vector<id_t> to_vector() const{
            std::vector<id_t> v;

            for (auto tid: range){

                v.push_back(tid);
                last_added = tid;
            }

            return v;

        }

        roaring::Roaring get_range() const {
            return this->range;
        }

        std::string to_string() const {

            std::string s  = "";
            for (auto tid : this->range){

                s+= std::format("{} ", tid);

            }
            return s;


        }

        void add(id_t tid){
            this->range.add(tid);
            last_added = tid;
        }

        
        void remove(id_t tid){
            this->range.remove(tid);
        }
    private:

        

        roaring::Roaring range;
        id_t roaring_bitmap_id;
        id_t last_added;

};
template <metric_space space>
class vc_dim_extractor{
    private:
        struct steps_taken{
            bool computed_frequent_patterns = false;
            bool converted_supports_to_bitsets = false;
        };

        struct transaction_set_data{

            std::size_t max_cardinality_of_shattered_superset = 0;
            std::size_t generator_sets = 0;
            bool shattered = false;

        };

        struct TrieNode {
            // unique_ptr automatically deletes the child when the map is destroyed
            using children_map_t = std::map<id_t, std::unique_ptr<TrieNode>>;
            std::map<id_t, std::unique_ptr<TrieNode>> children; //Like std::map<int, Trienode*> but safe
            bool is_end = false;
            
        };

        class Trie {
            public:
                std::unique_ptr<TrieNode> root; //Like TrieNode* but safe

                Trie() : root(std::make_unique<TrieNode>()) {} //safe pointer initialization

                void insert(const RangeRoaringBitmap<space>& itemset) {
                    TrieNode* node = root.get();
                    for (auto item : itemset.get_range()) {
                        if (node->children.find(item) == node->children.end())
                            node->children[item] = std::make_unique<TrieNode>();
                        node = node->children[item].get();
                    }
                    node->is_end = true;
                }

                bool contains(const RangeRoaringBitmap<space>& itemset) const {
                    TrieNode* node = root.get();
                    for (auto item : itemset.get_range()) {
                        auto it = node->children.find(item);
                        if (it == node->children.end()) return false;
                        node = it->second.get();
                    }
                    return node->is_end;
                }

        
        };

        static bool all_subsets_exist(const Trie& trie, const RangeRoaringBitmap<space>& candidate){

            for (auto& subset : candidate.subsets_of_size_minus_one()){

                if(!trie.contains(subset)){

                    return false;

                }

            }
            return true;

        }

        static void generate_candidates_recursive(TrieNode* node, RangeRoaringBitmap<space>& prefix,int target_k,const Trie& trie,std::vector<RangeRoaringBitmap<space>>& candidates) {
            // If we're at depth k-1 we generate size k candidates
            if ((int)prefix.count() == target_k - 2) {
                //std::cout<< std::format("Generating candidates of size {} from prefix {}\n", target_k, prefix.to_string());
                std::vector<id_t> keys;

                for (auto& [item, _] : node->children)
                    keys.push_back(item);

                for (int i = 0; i < (int)keys.size(); ++i) {
                    for (int j = i + 1; j < (int)keys.size(); ++j) {
                        RangeRoaringBitmap candidate(prefix);
                        candidate.add(keys[i]);
                        candidate.add(keys[j]);

                        if (all_subsets_exist(trie, candidate)){
                            //std::cout<< std::format("All subsets of size {} exist for candidate {}.\n", target_k-1, candidate.to_string());
                            candidates.push_back(candidate);
                        }
                    }
                }
                return;
            }

            // Otherwise go deeper
            for (auto& [item, child] : node->children) {
                prefix.add(item);
                generate_candidates_recursive(child.get(), prefix, target_k, trie, candidates);
                prefix.remove(item);
            }
        };
        static void generate_candidates(const Trie& trie, int target_k, std::vector<RangeRoaringBitmap<space>>& candidates){

            std::vector<RangeRoaringBitmap<space>> global_candidates;
            std::vector<typename TrieNode::children_map_t::iterator> children_to_explore;
            children_to_explore.reserve(trie.root->children.size());
            for (auto it = trie.root->children.begin(); it != trie.root->children.end(); ++it){
                children_to_explore.push_back(it);
            }
            #pragma omp parallel num_threads(6)
            {
                std::vector<RangeRoaringBitmap<space>> local_candidates;
                #pragma omp for schedule(dynamic)
                for (int i = 0; i< children_to_explore.size(); i++){
                    auto& [key, value] = *(children_to_explore[i]);
                    RangeRoaringBitmap<space> prefix(std::set<id_t>{});
                    prefix.add(key);

                    generate_candidates_recursive(value.get(), prefix, target_k, trie, local_candidates);
                }

                #pragma omp critical
                {
                    global_candidates.insert(global_candidates.end(), local_candidates.begin(), local_candidates.end());
                }

                
            }

            candidates = std::move(global_candidates);
        }
        static bool is_shattered(const RangeRoaringBitmap<space>& candidate, const std::vector<RangeRoaringBitmap<space>>& ranges){

            std::vector<id_t> elems;
            elems.reserve(candidate.count());

            // Extract elements from Roaring bitmap
            roaring::Roaring trajectory_ids = candidate.get_range();
            for (auto it = trajectory_ids.begin(); it != trajectory_ids.end(); ++it) {
                elems.push_back(*it);
            }

            int n = elems.size();

            // Iterate over all non-empty subsets
            for (uint64_t mask = 1; mask < (1ULL << n); ++mask) {
                RangeRoaringBitmap<space> subset(std::set<id_t>{});
                for (int i = 0; i < n; ++i) {
                    if (mask & (1ULL << i)) {
                        subset.add(elems[i]);
                    }
                }
                //std::cout << std::format("I am searching for a support for set {} while trying to shatter set {}\n", subset.to_string(), candidate.to_string());
                bool found_a_range = false;
                //Check it is shattered
                for (const auto& range : ranges){
                    if((candidate & range )== subset){
                        //std::cout << std::format("Subset {} and range {} have intersection {}\n", subset.to_string(), range.to_string(), (subset & range ).to_string());
                        found_a_range = true;
                        break;
                    }
                }
                if( !found_a_range){
                    return false;
                }
            }
            //std::cout<<std::format("Set {} is shattered.\n", candidate.to_string());
            return true;
        };
        static void generate_pairs_by_enumeration(Trie& T_2, std::map<RangeRoaringBitmap<space>, transaction_set_data>& F_2,std::map<RangeRoaringBitmap<space>, transaction_set_data>& F_1,std::vector<roaring::Roaring> inverted_index){

            //Trivial generation of F_2 by enumeration
            //Bag all the non pruned trajectories
            std::vector<id_t> tids;
            for (const auto& [key, data] : F_1){
                tids.push_back(*(key.get_range().begin()));
            }
            #pragma omp parallel num_threads(6)
            {
                std::map<RangeRoaringBitmap<space>, transaction_set_data> local_F_2;
                #pragma omp for schedule(dynamic, 8) 
                for(int i = 0; i< (int) tids.size(); i++){
                    for (int j = i + 1; j < (int)tids.size(); j++) {
                        RangeRoaringBitmap<space> key(std::set<id_t>{});
                        key.add(tids[i]);
                        key.add(tids[j]);

                        auto [it, inserted] = local_F_2.try_emplace(key);
                        auto& value = it->second;
                        value.generator_sets++;
                        value.shattered = true;
                        
                        auto intersection = inverted_index[tids[i]] & inverted_index[tids[j]];

                        value.max_cardinality_of_shattered_superset = floor(log2(intersection.cardinality() + 1) + 1);

                        if (intersection.cardinality() == 0||intersection.cardinality() == inverted_index[tids[i]].cardinality() ||intersection.cardinality() == inverted_index[tids[j]].cardinality()) {
                            local_F_2.erase(key);
                        }
                    }

                }

                //merge into T_2 
                #pragma omp critical
                {
                    for (auto& [k, v] : local_F_2) {
                        T_2.insert(k);
                        F_2[k] = v;
                    }
                }
            }
        };
        static void generate_pairs_from_supports(Trie& T_2, std::map<RangeRoaringBitmap<space>, transaction_set_data>& F_2,const std::vector<RangeRoaringBitmap<space>>& set_of_supports, const std::vector<roaring::Roaring>& inverted_index){
            #pragma omp parallel num_threads(6)
            {
                int tid = omp_get_thread_num();
                int nthreads = omp_get_num_threads();

                #pragma omp critical
                std::cout << "Hello from thread " << tid<< " / " << nthreads << std::endl;
                std::map<RangeRoaringBitmap<space>, transaction_set_data> local_F_2;

                #pragma omp for schedule(dynamic, 8)
                for (size_t idx = 0; idx < set_of_supports.size(); ++idx) {
                    const auto& s = set_of_supports[idx];

                    std::vector<id_t> tids;
                    tids.reserve(s.get_range().cardinality());
                    for (auto tid : s.get_range())
                        tids.push_back(tid);

                    for (int i = 0; i < (int)tids.size(); i++) {
                        for (int j = i + 1; j < (int)tids.size(); j++) {

                            RangeRoaringBitmap<space> key(std::set<id_t>{});
                            key.add(tids[i]);
                            key.add(tids[j]);

                            auto [it, inserted] = local_F_2.try_emplace(key);
                            auto& value = it->second;
                            value.generator_sets++;
                            value.shattered = true;

                            auto intersection = inverted_index[tids[i]] & inverted_index[tids[j]];

                            value.max_cardinality_of_shattered_superset = floor(log2(intersection.cardinality() + 1) + 1);

                            if (intersection.cardinality() == inverted_index[tids[i]].cardinality() ||intersection.cardinality() == inverted_index[tids[j]].cardinality()) {
                                local_F_2.erase(key);
                            }
                        }
                    }

                }
                //merge into T_2
                #pragma omp critical
                {
                    for (auto& [k, v] : local_F_2) {
                        T_2.insert(k);
                        F_2[k] = v;
                    }
                }
            }
            

            /*
            
            //Building F_2 from the actually appearing pairs. 
            int iter = 1;
            for (const auto& s : set_of_supports){
                //std::cout << std::format("I  am looking at support set for a pathlet whose id is {}, which appears in trajectories {}\n",s.get_range_id(), s.to_string());
                iter++;
                if(iter%1000==0){
                    std::cout<< std::format("Examining support set {}/{}...\n", iter, set_of_supports.size()); 
                }
                std::vector<id_t> tids;
                for (auto tid : s.get_range()){
                    //std::cout << "-> I am pushing into the vector the trajectory id "<< tid <<std::endl;
                    tids.push_back(tid);
                }

                for (int i = 0; i< tids.size(); i++){

                    for (int j = i+1; j< tids.size(); j++){
                        
                        RangeRoaringBitmap<space> key(std::set<id_t>{tids[i], tids[j]});

                        //std::cout << "My key is "<< key.to_string()<< std::endl;
                        //std::cout << "Key0s cardinality is "<< key.count()<< "while it should be 2 \n";
                        auto [it, inserted] = F[2].try_emplace(key);
                        auto& value = it->second;
                        value.generator_sets++;
                        value.shattered = true;
                        //F[2][key].generator_sets++;
                        
                        roaring::Roaring intersection = inverted_index[tids[i]] & inverted_index[tids[j]];
                        //std::cout << std::format("Intersection of common pathlets for trajectories {} and {} is {}\n", tids[i], tids[j], RangeRoaringBitmap<space>(intersection).to_string() );
                        value.max_cardinality_of_shattered_superset = floor(log2((intersection.cardinality()+1))+1);
                        //std::cout<<std::format("Intersection of common pathlets for trajectories {} and {} has max_cardinality of a shattered superset equal to {}\n ", tids[i],tids[j], value.max_cardinality_of_shattered_superset);
                        if(intersection.cardinality() == inverted_index[tids[i]].cardinality() || intersection.cardinality() == inverted_index[tids[j]].cardinality()){
                            //std::cout << std::format("Trajectories {} and {} always share the same support, so I remove the pair from F_2\n", tids[i], tids[j]);
                            F[2].erase(key);

                        }
                        
                        

                    }

                }
            }
            
            */


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
            std::set<RangeRoaringBitmap<space>> temporary_set_of_supports;
            if (!steps.converted_supports_to_bitsets){

                for (auto& fp : algo.freq_pathlets){
                    RangeRoaringBitmap<space> temp(fp.supporting_trajectories); // Adding +1 to support datasets with ids starting both at 0 and at 1
                    //try to find the set temp in set_of_supports
                    //if not found insert the set with a proper unique id 
                    if(temporary_set_of_supports.find(temp) == temporary_set_of_supports.end()){

                        temporary_set_of_supports.insert(RangeRoaringBitmap<space>(fp.supporting_trajectories, temporary_set_of_supports.size())); // the first inserted set has id zero 

                    }

                }
                //std::cout << "SUPPORTS "<< std::endl;
                int current_visiting_size = 1;
                std::cout << "Set of supports size : "<< temporary_set_of_supports.size() << std::endl;
                int count = 0;
                for (auto&s : temporary_set_of_supports){

                    if (current_visiting_size != s.count()){
                    
                        current_visiting_size = s.count();
                        std::cout << std::format("I have {} supports of size {}\n", count, current_visiting_size-1);
                        count = 0;
                    }
                    count++;
                    //std::cout << s.to_string_readable() <<std::endl;

                }
                //set becomes a vector for parallelization later on 
                std::vector<RangeRoaringBitmap<space>> supports_vec(temporary_set_of_supports.begin(),temporary_set_of_supports.end());
                set_of_supports = supports_vec;
                std::cout << std::format("I have {} supports of size {}\n", count, current_visiting_size);
                
                steps.converted_supports_to_bitsets = true;
            }

        }
        /*
        
        int compute_exact_vc_dimension(){
            assert(steps.converted_supports_to_bitsets && steps.computed_frequent_patterns);
            int max_shattered=  2;
            std::set<RangeRoaringBitmap<space>> shattered_subsets;
            //initialize the shattered subsets
            for (const auto& s: set_of_supports){
                for (auto& id: s){
                    shattered_subsets.insert(RoaringRangeBitmap(roaring::Roaring{id}));
                }
            }
            assert(shattered_subsets.size() ==dataset.num_trajs() );
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
                //I have to check it here.


                std::set<RangeRoaringBitmap<space>> subsets = s.subsets_of_size_minus_one();
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
                        //std::cout << "Shattered set "<<s.to_string_readable() << std::endl;
                    }
                    largest_size_shattered = s.count();
                }
            }

            return largest_size_shattered;
        }
        
        */
        int compute_exact_vc_dimension(){
            //initializa the F vector with the families. 
            std::cout << std::format("Number of distinct supports is {}\n", set_of_supports.size());
            std::set<RangeRoaringBitmap<space>> shattered_subsets;
            std::vector<std::map<RangeRoaringBitmap<space>, transaction_set_data>> F;
            F.push_back(std::map<RangeRoaringBitmap<space>, transaction_set_data>{}); //F_0 = empty map 
            F.push_back(std::map<RangeRoaringBitmap<space>, transaction_set_data>{}); //F_1  = empty map 
            F.push_back(std::map<RangeRoaringBitmap<space>, transaction_set_data>{}); //F_2  = empty map 
            std::vector<Trie> T;
            T.emplace_back();
            T.emplace_back();
            //Build the inverted index data structure to access one trajectory's pathets.
            //Meanwhile initialize the shattered_subsets at F_1
            std::vector<roaring::Roaring> inverted_index;
            inverted_index.reserve(dataset.num_trajectories_not_consecutive()+1);
            //I am now adding the individual trajectories to F_1
            for (const auto& s: set_of_supports){

                for (const auto& id: s.get_range()){
                    F[1][RangeRoaringBitmap<space>(std::set<id_t>{id})].generator_sets++;
                    F[1][RangeRoaringBitmap<space>(std::set<id_t>{id})].shattered = true;
                    //std::cout <<  "TID is "<< id <<std::endl;
                    if(id>=inverted_index.size()){
                        inverted_index.resize(id+1);
                    }
                    inverted_index[id].add(s.get_range_id());  

                    F[1][RangeRoaringBitmap<space>(std::set<id_t>{id})].max_cardinality_of_shattered_superset = floor(std::log2(inverted_index[id].cardinality()+1) +1);              
                }

            }
            std::cout << "Built inverted index...\n";
            for (size_t tid = 0; tid< inverted_index.size(); tid++){
                bool found_negation = false;
                RangeRoaringBitmap<space> individual_set(std::set<id_t>{tid});
                for (auto& s : set_of_supports){
                    if ((individual_set.get_range() & s.get_range()).cardinality() == 0){
                        found_negation = true;
                        break;
                    }
                }
                if(!found_negation){
                    //Remove 
                    F[1].erase(individual_set);
                    std::cout << std::format("I am removing trajectory {} from F_1 because it always appears. \n", tid);
                }
                //std::cout << std::format("Pathlet ids which are contained in trajectory {} are : {}\n", tid, RangeRoaringBitmap<space>(inverted_index[tid]).to_string() );

            }
            //std::cout << "I manage to finish the insertions \n";

            assert(steps.converted_supports_to_bitsets && steps.computed_frequent_patterns);
            
            //std::cout << std::format(" Shattered sets have size: {}. \n The number of trajectories is {}\n ", shattered_subsets.size(), dataset.num_trajectories_not_consecutive() );
           
            int current_visiting_size =2;
            int largest_size_shattered = 1;


            bool found_shattered = false;
            int effective_pairs = 0;
            int possible_pairs = (dataset.num_trajectories_not_consecutive()-1) *dataset.num_trajectories_not_consecutive()/2;

            for(const auto& s: set_of_supports){
                effective_pairs +=s.count()*(s.count()-1) /2;
            }

            T.push_back(Trie());
            //Fill up F_2 with the appearing pairs.
            if(effective_pairs< possible_pairs){
                
                generate_pairs_from_supports(T[2], F[2], set_of_supports, inverted_index);
            }
            else{
                generate_pairs_by_enumeration(T[2], F[2], F[1], inverted_index);
            }
            
            
            std::cout << "Number of elements of F_2 is "<<F[2].size()<< std::endl;
            //std::cout << "Cleaning up pairs that always share the same support... \n";
            //std::erase_if(F[2], [](const auto& item) {auto const& [key, value] = item;  return value.max_cardinality_of_shattered_superset < 2;});
            //std::erase_if(F[2], [](const auto& item) {auto const& [key, value] = item;  return value.generator_sets < 2;});
            //std::cout << "Number of elements of F_2 is "<<F[2].size()<< std::endl;
            //I now  convert F_2 into a Trie to speed up generation of F_3
            
            /*
            
            std::cout<< "I am adding the pairs to the Trie...\n";
            for (const auto pair_and_data: F[2]){
                //std::cout << std::format("Inserting pair {} into the Trie\n", pair_and_data.first.to_string());
                T[2].insert(pair_and_data.first);

            }
            */

            int k = 3;
            while(true){
                std::cout << std::format("Generating candidates of size {}...\n", k);
                T.push_back(Trie());
                std::vector<RangeRoaringBitmap<space>> candidates;
                RangeRoaringBitmap<space> prefix(std::set<id_t>{});

                generate_candidates(T[k-1],  k, candidates);

                //I already check all subsets are present in the candidate generation phase. 
                int insertions = 0;

                //PARALLELIZE CANDIDATE EVALUATION
                #pragma omp parallel num_threads(6)
                {
                    std::set<RangeRoaringBitmap<space>> local_Fk;
                    int local_insertions = 0;
                    #pragma omp for schedule(dynamic, 8)
                    for (int i = 0; i<candidates.size(); i++){
                        auto c = candidates[i];
                        
                        RangeRoaringBitmap<space> intersection(inverted_index[*(c.get_range().begin())]);
                        for (const auto tid: c.get_range()){
                            intersection = intersection & inverted_index[tid];
                            if(intersection.count() == 0){
                                break;
                            }
                        }
                        if(intersection.count() == 0){
                            continue;
                        }
                        
                        bool a_subset_cant_be_separated_from_c = false;
                        
                        //se l'intesezione ha la stessa cardinalità dei sottinsiemi di taglia -1 it means it is not shatterable, so I won't bother inserting it.
                        for (const auto& subset: c.subsets_of_size_minus_one()){
                            RangeRoaringBitmap<space> intersection_subset(inverted_index[*(subset.get_range().begin())]);
                            for (const auto tid: subset.get_range()){
                                intersection_subset = intersection_subset & inverted_index[tid];
                            }
                            if(intersection_subset.count() == intersection.count()){
                                a_subset_cant_be_separated_from_c = true;
                                
                                break;
                            }
                            assert(intersection_subset.count()> intersection.count());

                        }
                        if(a_subset_cant_be_separated_from_c){
                            //std::cout << std::format("Set {} and its subset always share the same support, so I don't add the first set to T[{}]\n", c.to_string(), k);
                            continue;
                        }
                        if(!a_subset_cant_be_separated_from_c && is_shattered(c,set_of_supports)){
                            //std::cout << std::format("inserting the {}-ple {} into the Trie, as the trajectories' intersection is {}\n",k, c.to_string(), intersection.to_string());
                            local_Fk.insert(c);
                            local_insertions++;
                        }
                    }

                    //merge into Tk
                    #pragma omp critical
                    {
                        for (auto c : local_Fk){
                            T[k].insert(c);
                        }
                        insertions += local_insertions;
                    }


                }
                /*
                for (const auto& c : candidates){
                    //CHECK Shattering
                    RangeRoaringBitmap<space> intersection(inverted_index[*(c.get_range().begin())]);
                    for (const auto tid: c.get_range()){
                        intersection = intersection & inverted_index[tid];
                        if(intersection.count() == 0){
                            break;
                        }
                    }
                    if(intersection.count() == 0){
                        continue;
                    }
                    
                    bool a_subset_cant_be_separated_from_c = false;
                    
                    //se l'intesezione ha la stessa cardinalità dei sottinsiemi di taglia -1 it means it is not shatterable, so I won't bother inserting it.
                    for (const auto& subset: c.subsets_of_size_minus_one()){
                        RangeRoaringBitmap<space> intersection_subset(inverted_index[*(subset.get_range().begin())]);
                        for (const auto tid: subset.get_range()){
                            intersection_subset = intersection_subset & inverted_index[tid];
                        }
                        if(intersection_subset.count() == intersection.count()){
                            a_subset_cant_be_separated_from_c = true;
                            
                            break;
                        }
                        assert(intersection_subset.count()> intersection.count());

                    }
                    if(a_subset_cant_be_separated_from_c){
                        //std::cout << std::format("Set {} and its subset always share the same support, so I don't add the first set to T[{}]\n", c.to_string(), k);
                        continue;
                    }
                    if(!a_subset_cant_be_separated_from_c && is_shattered(c,set_of_supports)){
                        //std::cout << std::format("inserting the {}-ple {} into the Trie, as the trajectories' intersection is {}\n",k, c.to_string(), intersection.to_string());
                        T[k].insert(c);
                        insertions++;
                    }
                    
                */
                
                    
                    /*
                    
                    if(intersection.count() == 0 || intersection.count()==1){
                        continue;
                    }
                    bool a_subset_cant_be_separated_from_c = false;
                    
                    //se l'intesezione ha la stessa cardinalità dei sottinsiemi di taglia -1 it means it is not shatterable, so I won't bother inserting it.
                    for (const auto& subset: c.subsets_of_size_minus_one()){
                        RangeRoaringBitmap<space> intersection_subset(inverted_index[*(subset.get_range().begin())]);
                        for (const auto tid: subset.get_range()){
                            intersection_subset = intersection_subset & inverted_index[tid];
                        }
                        if(intersection_subset.count() == intersection.count()){
                            a_subset_cant_be_separated_from_c = true;
                            
                            break;
                        }
                        assert(intersection_subset.count()> intersection.count());

                    }
                    if(a_subset_cant_be_separated_from_c){
                        //std::cout << std::format("Set {} and its subset always share the same support, so I don't add the first set to T[{}]\n", c.to_string(), k);
                        continue;
                    }
                    
                    //TODO: CHECK SHATTERING BEFORE INSERTION or before giving up due to not usefulness at higher levels.
                
                    
                
                    
                */
                
                std::cout<< std::format("I have generated {} candidates of size {}, but due to lack of supports I am adding {}.\n", candidates.size(),k, insertions);
                
                if(insertions==0){
                    std::cout<< std::format("I could not fill the trie at size {}, so we break the apriori-like cycle.\n",k);
                    break;
                }
                k++;

            }
            return k-1;
        }
        
    private:
    const freq_subtrajectory_algo_output_config output_config{
        .maximal = false,
        .keep_matching_ids = true,
        .min_length= 1
    };
    trajectory_t dataset;
    frequent_subtrajectory_algo_t algo;
    std::vector<RangeRoaringBitmap<space>> set_of_supports;
    steps_taken steps;
    
};

}


