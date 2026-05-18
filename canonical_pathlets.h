#pragma once
#include <cmath>
#include <vector>
#include <queue>
#include "roaring.hh"

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
        using id_t = trajectory_t::id_t;
    public:
        bool isNULL;
        bool frequent; //Constructor sets it to true by default
        float frequency= 0.0;
        int father=-1;
        int left_child=-1;
        int right_child=-1;
        subtrajectory_t pathlet;
        roaring::Roaring supporting_trajectories;
    
        PathletNode(int left, int right) : pathlet({left, right}), frequent(true), isNULL(false){}
        PathletNode(int left, int right, int father_idx, int left_child_idx, int right_child_idx) : pathlet({left, right}), frequent(true), isNULL(false),father(father_idx), right_child(right_child_idx), left_child(left_child_idx){}

        inline int getLength(){ return this->pathlet.second - this-> pathlet.first + 1; }
        inline float getFrequency() { return this->frequency; }
        inline subtrajectory_t getPathlet(){ return this->pathlet; }
        inline bool isFrequent() { return this-> frequent; }
        inline void setNULL(){this->isNULL = true;}
        inline void addId(id_t id){this->supporting_trajectories.add(id);}
        //returns true if this is contained by other
        inline int getFather(){return this-> father;}
        inline int getRightChild(){return this->right_child;}
        inline int getLeftChild(){return this->left_child;}
        inline bool is_contained_by(const PathletNode& other){return other.pathlet.first<=this->pathlet.first && other.pathlet.second>=this->pathlet.second;}
        inline roaring::Roaring getSupportingTrajectories(){return this->supporting_trajectories;}
        inline bool isLeaf() {return left_child ==-1 && right_child == -1;}
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
        BinaryPathletTree(){
            
        }
        
        BinaryPathletTree(trajectory_t& supp, id_t traj_id, int max_depth, int min_length){
            this->support = supp;
            this-> n = support.get_trajectory_size(traj_id);
            this-> l = min_length;
            this-> d = max_depth;
            this->trajectory_id = traj_id;
            //assert(min_length <= n);

            assert(this->pathlet_collection.empty()); //assert flat_tree is empty

            this->build(); //Fills the pathlet_collection vector, which is then explored as a binary tree
        }

        //Getters
        inline int getDepth() { return this-> d;}
        inline int getMinLength() { return this-> l; }
        inline int getTrajectoryLength() {return this-> n;}
        inline int getPathletCollectionLength(){return this->pathlet_collection.size();}
        id_t getTrajectoryId(){ return this->trajectory_id; }

        //Setters
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
        
        PathletNode& getNodeAt(int node_idx){

            assert(node_idx < pathlet_collection.size());
            return pathlet_collection.at(node_idx); //Should be returning without the id set 

        }

        
        std::string toString(){
            
            std::string s = "Binary Pathlet Tree : ID "+ std::to_string(this->trajectory_id)+ "\n";
            for (int i=0; i< pathlet_collection.size(); i++){
                PathletNode p = pathlet_collection.at(i);
                if (p.isNULL){
                    s += std::format("Node {} Pathlet: NULL ", i);
                }
                else{
                    s += std::format("Node {} Pathlet: [ {}, {}]", i, p.getPathlet().first, p.getPathlet().second);
                }
                
                s += "\n";
            }
            s += "\n";
            return s;
        }
        //Given a point_idx in the support trajectory returns the node of the tree with the single-point pathlet containing it
        PathletNode getPointPathletNode(int point_idx){
        
            //std::cout << " ***** Searching for point pathlet "<< point_idx << " *****"<< std::endl;
            //top-down search
            bool found = false;
            PathletNode current_node = pathlet_collection.at(0); // root as the first node of the visit

            int current_idx = 0;
            if (current_node.getLength() ==1 && this->n ==1){
                return current_node;
            }
            //std::cout << "** Root pahtlet is "<< current_node.getPathlet().first << " "<< current_node.getPathlet().second<< " **"<<std::endl;
            //std::cout << " root pathelt should however be "<< pathlet_collection.at(0).getPathlet().first << " "<< pathlet_collection.at(0).getPathlet().second<< " **"<<std::endl;
            while(true){
                //exit conditions
                assert(!current_node.isNULL);
                if (current_node.getLength()==1 && (current_idx !=0 && this->n >1)){
                    //assert(current_node.getPathlet().first== point_idx);
                    //std::cout << "I am seeing a size one pathlet at point "<< current_node.getPathlet().first << " for point "<< point_idx<< std::endl;
                    return current_node;
                    
                }
                //traversal

                int left_child_location = left_child_idx(current_idx);
                PathletNode left_child = getNodeAt(left_child_location);
                //std::cout << "** Left Child pathlet is "<< left_child.getPathlet().first << " "<< left_child.getPathlet().second<< " **"<<std::endl;

                int right_child_location = right_child_idx(current_idx);
                PathletNode right_child = getNodeAt(right_child_location);
                //std::cout << "** Right Child pathlet is "<< right_child.getPathlet().first << " "<< right_child.getPathlet().second<< " **"<<std::endl;

                if (!left_child.isNULL && point_idx >= left_child.getPathlet().first && point_idx <= left_child.getPathlet().second){
                    current_idx = left_child_location;
                    current_node = getNodeAt(current_idx);
                }
                else{
                    //go left 
                    current_idx = right_child_location;
                    current_node = getNodeAt(current_idx);
                }

            }

            
        }
        std::vector<PathletNode> getMinLengthPathletsAtLowLevels(int min_length){
            //Top down visit of the tree which retrieves all pathlets above the length threshold.
            //Iterate through levels, from left to right, stop when the whole level if split would go under the threshold
            std::vector<PathletNode> min_length_pathlets;

            //Check the root
            if( pathlet_collection.at(0).getLength()>=min_length){
                min_length_pathlets.push_back(pathlet_collection.at(0));
            }
            else{
                return min_length_pathlets;
            }
            //Check the children
            int level_beginning = 1;
            for (int level =1; level <= this->d; level++){

                level_beginning = (int( pow(2, level))) -1;
                //Iterate through the level 
                int max_length_level = -1;
                for (int k = 0; k<= level_beginning; k++){

                    int position = level_beginning + k;
                    PathletNode current  = this->pathlet_collection.at(position);
                    if (current.isNULL){
                        continue;
                    }
                    if(current.getLength()>=min_length){
                        min_length_pathlets.push_back(current);

                    }
                    if(current.getLength()>=max_length_level){
                        max_length_level= current.getLength();
                    }
                    
                }

                if(max_length_level/2.0 < min_length -1){
                    break;
                }
            }
            //std::cout << std::format("For pathlet mother {} I have found {} pathlets at the lowest level of min length {}\n", this->getTrajectoryId(), min_length_pathlets.size(), min_length);
            std::reverse(min_length_pathlets.begin(), min_length_pathlets.end());
            for (size_t i = 0; i<min_length_pathlets.size()-1; i++){

                for (size_t j=i+1; j<min_length_pathlets.size(); j++){

                    if(min_length_pathlets[i].is_contained_by(min_length_pathlets[j])){
                        min_length_pathlets.erase(min_length_pathlets.begin() + j);
                    }

                }

            }
            return min_length_pathlets;
        }
        std::vector<PathletNode> getMinLengthPathlets(int min_length){
            //Top down visit of the tree which retrieves all pathlets above the length threshold.
            //Iterate through levels, from left to right, stop when the whole level if split would go under the threshold
            std::vector<PathletNode> min_length_pathlets;

            //Check the root
            if( pathlet_collection.at(0).getLength()>=min_length){
                min_length_pathlets.push_back(pathlet_collection.at(0));
            }
            else{
                return min_length_pathlets;
            }
            //Check the children
            int level_beginning = 1;
            for (int level =1; level <= this->d; level++){

                level_beginning = (int( pow(2, level))) -1;
                //Iterate through the level 
                int max_length_level = -1;
                for (int k = 0; k<= level_beginning; k++){

                    int position = level_beginning + k;
                    PathletNode current  = this->pathlet_collection.at(position);
                    if (current.isNULL){
                        continue;
                    }
                    if(current.getLength()>=min_length){
                        min_length_pathlets.push_back(current);

                    }
                    if(current.getLength()>=max_length_level){
                        max_length_level= current.getLength();
                    }
                    
                }

                if(max_length_level/2.0 < min_length -1){
                    break;
                }
            }
            //std::cout << std::format("For pathlet mother {} I have found {} pathlets of min length {}\n", this->getTrajectoryId(), min_length_pathlets.size(), min_length);
            return min_length_pathlets;
        }

        static inline int left_child_idx(int node_idx){

            return 2 * node_idx + 1;

        }

        static inline int right_child_idx(int node_idx){

            return 2* node_idx + 2;

        }

        private: 

        std::vector<PathletNode> pathlet_collection; 
        int n;
        int d;
        int l;
        trajectory_t support;
        id_t trajectory_id;

        static inline int get_father(int node_idx){

            assert(node_idx > 0);

            int father =  (node_idx % 2 ) ? (node_idx - 1) / 2 : (node_idx - 2) / 2; 

            return father;
        }
        /*
        void build(){

            pathlet_collection.clear();
            std::queue<int> nodes_to_process;
            
            //Create the root node:
            pathlet_collection.emplace_back(PathletNode(0, this->n-1) );
            queue.push(0);
            while(!q.empty()){

                int current = q.front();
                PathletNode& node = pathlet_collection[current];
                //If i don't need to expand the node, I continue
                if(node.getLength() <= this->l){
                    continue;
                }

                int mid = ((node.right - node.left + 1) /  2) - 1 + node.left;

                //Create left child - if it is admitted. 
                if(mid - node.left >=0){
                    PathletNode left_child(node.left, mid);
                    left_child.father = current;
                    int left_child_vector_idx = pathlet_collection.size();
                    node.left_child = left_child_vector_idx;
                    pathlet_collection.push_back(left_child);

                }

                //Create right child
                if(node.right - (mid+1) >=0){
                    PathletNode right_child(mid+1, node.right);
                    right_child.father = current;
                    int right_child_vector_idx = pathlet_collection.size();
                    node.rigth_child = right_child_vector_idx;
                    pathlet_collection.push_back(right_child);

                }

            }

        }
        
        */
        //Builds the binary tree of pathlets
        
        void build(){

            //Root is the first node of vector pathlet_collection
            PathletNode root(0, this-> n - 1);
            this-> pathlet_collection.push_back(root);
            assert(pathlet_collection.size() == 1);
            int level_beginning = 1;

            //Fill the vector level by level, from left to right
            for (int level =1; level <= this->d; level++){

                level_beginning = (int( pow(2, level))) -1;
                
                for (int k = 0; k<= level_beginning; k++){

                    int position = level_beginning + k;
                    PathletNode parent  = this->pathlet_collection.at(get_father(position));
                    //std::cout << "Level "<< level<< "pos "<< position<< " father length "<<  pathlet_collection.at(get_father(position)).getPathlet().second- pathlet_collection.at(get_father(position)).getPathlet().first +1<< std::endl;
                    //If the parent is a NULL pathlet (empty node of the vector) or is shorter than the mininmum length, then insert a NULL node
                    if( parent.isNULL || parent.getLength() <= this-> l){
                        
                        PathletNode nullKid(-1,-1);
                        nullKid.setNULL();
                        this-> pathlet_collection.push_back(nullKid);

                    }
                    else{ //The current parent node is valid, so create its kids

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
