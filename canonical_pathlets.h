#pragma once
#include <cmath>
#include <vector>
#include <queue>
#include "trajectory.h"
#include "metric_space.h"
namespace frechet{

template<metric_space m_space>
class PathletNode{
    public:
        using space = m_space;
        using index_t = std::size_t;
        using trajectory_t = trajectory_collection<space>;
        using subtrajectory_t = trajectory_t::subtrajectory_t;

    public:
        bool isNULL;
        bool frequent;
        float frequency;
        const subtrajectory_t pathlet;

        PathletNode(const int left, const int right) : pathlet({left, right}), frequent(true), isNULL(false){}
        inline int getLength(){ return this->pathlet.second - this-> pathlet.first + 1; }
        inline float getFrequency() { return this->frequency; }
        inline subtrajectory_t getPathlet(){ return this->pathlet; }
        inline bool isFrequent() { return this-> frequent; }
        inline void setNULL(){this->isNULL = true;}
        //inline bool isNULL() {return this->isNULL;}
};

//represents a collection of pathlets organized as a tree, based on a subtrajectory taken from a trajectory collection
template<metric_space m_space>
class BinaryPathletTree{
    
    public:
        using space = m_space;
        using index_t = std::size_t;
        using trajectory_t = trajectory_collection<space>;
        using subtrajectory_t = trajectory_t::subtrajectory_t;
        using cluster_t = std::vector<subtrajectory_t>;
        using id_t = trajectory_t::id_t;
        using PathletNode = frechet::PathletNode<space>;

    public:

        
        BinaryPathletTree(trajectory_t& supp, id_t traj_id, int max_depth, int min_length){
            this->support = supp;
            this-> n = support.get_trajectory_size(traj_id);
            this-> l = min_length;
            this-> d = max_depth;
            this->trajectory_id = traj_id;
            //assert(min_length <= n);

            assert(this->pathlet_collection.empty()); //assert flat_tree is empty

            this->build();
        }

        inline int getDepth() { return this-> d;}
        inline int getMinLength() { return this-> l; }
        inline int getTrajectoryLength() {return this-> n;}
        inline int getPathletCollectionLength(){return this->pathlet_collection.size();}
        

        void setInfrequent(int node_idx) {

            PathletNode& p = this-> pathlet_collection.at(node_idx);
            p.frequent = false;
            assert(p.frequent == this-> pathlet_collection.at(node_idx).isFrequent());
            int father = node_idx;
            while (father > 0) {
                
                father = get_father(father);
                PathletNode& ap = this-> pathlet_collection.at(father);
                if (!ap.frequent){ break; }
                ap.frequent = false;
            }
            //bottom-up visit of the tree, only where needed
        }

        void setEstimatedFrequency(int node_idx, float freq){

            pathlet_collection.at(node_idx).frequency = freq;


        }

        id_t getTrajectoryId(){ return this->trajectory_id; }

        /**
         * Retrieves ALL the pathlets currently marked as frequent in the tree. NO GUARANTEES ON THE CORRECTNESS IF THIS IS CALLED BEFORE PROCESSING.
        */
        /*cluster_t retrieveFrequentPathlets(){

            cluster_t freq_pathlets;
            std::queue<int> q;

            q.enqueue(0);

            while(!q.empty()){

                int current_visinting_index = q.dequeue();
                PathletNode p = pathlets.at(current_visinting_index);
                if(p == NULL){ continue; } // i am at a leaf
                if(p.frequent){

                    freq_pathlets.push_back(p.getPathlet());

                }
                else{

                    q.enqueue(left_child_idx(current_visinting_index));
                    q.enqueue(right_child_idx(current_visinting_index));

                }
            }
            
            return freq_pathlets;
            //top-down visit of the tree 

        }*/
        
        PathletNode getNodeAt(int node_idx){

            assert(node_idx < pathlet_collection.size());
            return pathlet_collection.at(node_idx);

        }



        private: //class variables + tree access methods

        std::vector<PathletNode> pathlet_collection; 
        int n;
        int d;
        int l;
        trajectory_t support;
        id_t trajectory_id;

        static inline int left_child_idx(int node_idx){

            return 2 * node_idx + 1;

        }

        static inline int right_child_idx(int node_idx){

            return 2* node_idx + 2;

        }

        static inline int get_father(int node_idx){

            assert(node_idx > 0);

            int father =  (node_idx % 2 ) ? (node_idx - 1) / 2 : (node_idx - 2) / 2; 

            return father;
        }

        void build(){
            PathletNode root(0, this-> n - 1);
            this-> pathlet_collection.push_back(root);
            assert(pathlet_collection.size() == 1);
            int level_beginning = 1;

            for (int level =1; level <= this->d; level++){

                level_beginning = (int( pow(2, level))) -1;
                
                for (int k = 0; k<= level_beginning; k++){

                    int position = level_beginning + k;
                    PathletNode parent  = this->pathlet_collection.at(get_father(position));
                    //std::cout << "Level "<< level<< "pos "<< position<< " father length "<<  pathlet_collection.at(get_father(position)).getPathlet().second- pathlet_collection.at(get_father(position)).getPathlet().first +1<< std::endl;
                    if( parent.isNULL || parent.getLength() <= this-> l){
                        
                        PathletNode nullKid(-1,-1);
                        nullKid.setNULL();
                        this-> pathlet_collection.push_back(nullKid);

                    }
                    else{

                        int left = parent.pathlet.first;
                        int right = parent.pathlet.second;
                        if( position % 2){
                            //left child
                            PathletNode kid(left, ((right - left + 1) /  2) - 1 + left);
                            this->pathlet_collection.push_back(kid);

                        }
                        else{
                            PathletNode kid( (right - left + 1) /  2 + left, right);
                            this->pathlet_collection.push_back(kid);

                        }

                    }

                }

                 //TODO:rephrase

            }

            //ASSERTIONS FOR TOY DS
            //assert(pathlet_collection.at(1).getPathlet().first == 0);
            //std::cout << "This is the second extreme of the 1st node of the second levedl (idx = 1)"<< pathlet_collection.at(1).getPathlet().second<< std::endl;
        }

        
};
}
