#pragma once

#include <boost/property_map/property_map.hpp>
#include <chrono>
#include <iostream>
#include <optional>
#include <type_traits>
#include "roaring.hh"
#include <ankerl/unordered_dense.h>

#include "metric_space.h"
#include "trajectory.h"

namespace frechet {

// A boost property map interface for `trajectory_collection<m_space>`.
template<metric_space m_space>
class trajectory_property_map_adapter {

public:
    using trajectory_t = trajectory_collection<m_space>;

    using value_type = m_space::point_t;
    using reference = const value_type&;
    using key_type = trajectory_t::index_t;
    using category = boost::lvalue_property_map_tag;

    trajectory_property_map_adapter(const trajectory_t &data) : data(data) {}

    reference operator[](const key_type &key) const {
        return data[key];
    }
    trajectory_t::index_t size() const{
        return data.get_actual_size();
    }
private:
    const trajectory_t &data;

};


template<typename T>
const T::value_type& get(const T &map, const typename T::key_type &key) {
    return map[key];
}


// Additional concepts

template<typename T>
concept has_CGAL_kernel = requires {
    typename T::kernel;
};

// Some more stuff

namespace internal {

template<typename pair_type>
struct pair_hash {
    using is_avalanching = void;

    [[nodiscard]] uint64_t operator()(const pair_type &s) const {
        using namespace ankerl::unordered_dense::detail::wyhash;
        const auto hash_a = hash(s.first);
        // TODO: might want to consider different ways to combine 
        return std::rotl(hash_a, 1) + hash_a + hash(s.second);
    }
};

template<typename T, typename Compare>
bool optimize(std::optional<T> &cur, T const&cand){
    if(!cur || Compare{}(cand, *cur)){
        cur = cand;
        return true;
    }
    return false;
}
template<typename T>
bool maximize(std::optional<T> &cur, T const&cand){
    return optimize<T, std::greater<>>(cur, cand);
}
template<typename T>
bool minimize(std::optional<T> &cur, T const&cand){
    return optimize<T, std::less<>>(cur, cand);
}

struct Timer{
    using clock_t = std::chrono::steady_clock;
    using time_point_t = decltype(clock_t::now());

    Timer() : start(clock_t::now()) {}

    double seconds() const {
        const auto now = clock_t::now();
        const auto nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(now - start).count();
        return nanoseconds * 1e-9;
    }

    time_point_t start;
};


template<typename T>
class prefix_sum{
public:
    using value_type = T;
    using storage_t = std::vector<T>;
    using index_t = std::size_t;

    prefix_sum(const storage_t &values) : data(values.size() + 1) {
        std::partial_sum(values.begin(), values.end(), data.begin() + 1);
    }

    // sum over [l, r)
    value_type sum(index_t l, index_t r) const {
        return data[r] - data[l];
    }

private:
    storage_t data;
};

} // namespace internal
template <metric_space space>
class RangeRoaringBitmap{

        public:
        using id_t = unsigned int;

        RangeRoaringBitmap(){
            roaring::Roaring temporary{};
            range = std::move(temporary);
        }
        
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
        /*
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
        
        */

        
        bool operator<(const RangeRoaringBitmap<space>& b) const {
            if (range.cardinality() != b.range.cardinality())
                return range.cardinality() < b.range.cardinality();

            auto ita = range.begin();
            auto itb = b.range.begin();

            for (; ita != range.end() && itb != b.range.end(); ++ita, ++itb) {
                if (*ita != *itb)
                    return *ita < *itb;
            }

            return false;
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

        void set_range_id(id_t range_id) {

            this->roaring_bitmap_id = range_id;

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

        void clear(){
            this->range = roaring::Roaring{};
        }

        //Risky bit
        roaring::Roaring range;
    private:

        

        
        id_t roaring_bitmap_id;
        id_t last_added;

};
template <metric_space space>
class RangeRoaringBitmap64{

        public:
        using id_t = long unsigned int;

        RangeRoaringBitmap64(){
            roaring::Roaring64Map temporary{};
            range = std::move(temporary);
        }
        
        RangeRoaringBitmap64(int num_trajectories, std::set<id_t> range_ids): roaring_bitmap_id(std::numeric_limits<id_t>::max()) {
            for (id_t tid : range_ids){

                range.add(tid);
                last_added = tid;
            }
        }

        RangeRoaringBitmap64(int num_trajectories, roaring::Roaring64Map range_ids): roaring_bitmap_id(std::numeric_limits<id_t>::max()){
            range= range_ids;
        }

        RangeRoaringBitmap64(roaring::Roaring64Map range_ids) : roaring_bitmap_id(std::numeric_limits<id_t>::max()){

            range = range_ids;

        }
        
        RangeRoaringBitmap64(const roaring::Roaring64Map& range_ids, id_t roaring_bitmap_id ): roaring_bitmap_id(roaring_bitmap_id){

            range = range_ids;
            //this-> roaring_bitmap_id = roaring_bitmap_id;

        }

        RangeRoaringBitmap64(std::set<id_t> range_ids, id_t roaring_bitmap_id){
            for (id_t tid : range_ids){

                range.add(tid);

            }
            this-> roaring_bitmap_id = roaring_bitmap_id;
        }

        RangeRoaringBitmap64(std::set<id_t> range_ids) : roaring_bitmap_id(std::numeric_limits<id_t>::max()){
            for (id_t tid : range_ids){

                range.add(tid);
                last_added = tid;

            }
        }
        /*
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
        
        */

        
        bool operator<(const RangeRoaringBitmap64<space>& b) const {
            if (range.cardinality() != b.range.cardinality())
                return range.cardinality() < b.range.cardinality();

            auto ita = range.begin();
            auto itb = b.range.begin();

            for (; ita != range.end() && itb != b.range.end(); ++ita, ++itb) {
                if (*ita != *itb)
                    return *ita < *itb;
            }

            return false;
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
        bool operator==(const RangeRoaringBitmap64<space>& other) const{

            return this->range== other.range; //built in method of roaring class
        }

        int count() const{
            return this->range.cardinality();
        }

        id_t get_range_id() const {

            return this->roaring_bitmap_id;
        }

        void set_range_id(id_t range_id) {

            this->roaring_bitmap_id = range_id;

        }
        // When invoked returns a vector with all the subsets of this with this.count()-1 bits 
        std::set<RangeRoaringBitmap64<space>> subsets_of_size_minus_one() const{
            std::set<RangeRoaringBitmap64> subsets;
            for (auto tid : this->range){
                RangeRoaringBitmap64<space> subset(this->range);
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

        RangeRoaringBitmap64<space> operator|(const RangeRoaringBitmap64& other ) const{

            return RangeRoaringBitmap64(this->range | other.range);

        }

        RangeRoaringBitmap64<space> operator&(const RangeRoaringBitmap64& other ){

            return RangeRoaringBitmap64(this->range & other.range);

        }

        RangeRoaringBitmap64<space> operator&(const RangeRoaringBitmap64& other ) const{

            return RangeRoaringBitmap64(this->range & other.range);

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

        roaring::Roaring64Map get_range() const {
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

        void clear(){
            this->range = roaring::Roaring64Map{};
        }

        //Risky bit
        roaring::Roaring64Map range;
    private:

        

        
        id_t roaring_bitmap_id;
        id_t last_added;

};
std::vector<int> POWERS_OF_TWO = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912};
std::vector<long unsigned int> POWERS_OF_TEN = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000};
} // namespace frechet

