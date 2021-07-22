/// \ingroup base
/// \class ttk::AssignmentExhaustive
/// \author Mathieu Pont <mathieu.pont@outlook.com>
///
/// Exhaustive Search for Unbalanced Assignement Problem
///   (TODO manage balanced problem)
///
/// The cost matrix in input has a size of (n + 1) x (m + 1)
/// - n is the number of jobs, m the number of workers
/// - the nth row contains the cost of not assigning workers
/// - the mth column is the same but with jobs
/// - the last cell (costMatrix[n][m]) is not used
///
/// An exhaustive search for the assignment problem has an exponential complexity
/// Use this algorithm only for small assigment problem


#ifndef _ASSIGNMENTEXHAUSTIVE_H
#define _ASSIGNMENTEXHAUSTIVE_H

namespace ttk {

  class AssignmentExhaustive : public Debug {

  public:
    AssignmentExhaustive(){};

    ~AssignmentExhaustive(){};

    template <typename dataType>
    int run(std::vector<matchingTuple> &matchings);
    
    template <typename dataType>
    dataType tryAssignment(std::vector<int> &asgn, 
                           std::vector<matchingTuple> &matchings);
    
    template<class dataType>
    std::string vectorToString(std::vector<dataType> &my_vector){
      std::stringstream result;
      std::copy(my_vector.begin(), my_vector.end(), std::ostream_iterator<dataType>(result, ""));
      return result.str();
    }
    
    //std::vector<std::vector<int>> constructAssignments(int min_dim, int max_dim);
    //std::vector<std::vector<int>> constructAssignments2(int min_dim, int max_dim);
    //std::vector<std::vector<int>> constructAssignments3(int min_dim, int max_dim);
    
    // This version creates a lot of duplicates (worst version)
    //std::vector<std::vector<int>> AssignmentExhaustive::constructAssignments(int min_dim, int max_dim){
    std::vector<std::vector<int>> constructAssignments(int min_dim, int max_dim){
      std::vector<std::vector<int>> allAsgn;
      std::queue<std::vector<int>> queue;
      std::set<std::string> unasgnAdded;
      
      std::vector<int> asgn_temp;
      queue.push(asgn_temp);
      
      while(!queue.empty()){
        std::vector<int> asgn = queue.front();
        queue.pop();
        int ind = asgn.size();
        for(int j = 0; j < max_dim; ++j){
          std::vector<int> new_asgn(asgn);
          // TODO avoid high complexity of find in vector
          std::vector<int>::iterator it = std::find(new_asgn.begin(), new_asgn.end(), j);
          if (it == new_asgn.end()){ // not found
            new_asgn.push_back(j);
          }else
            continue;
          
          if(ind == max_dim-1){
            // A new assignement without unassigned costs is made here
            //allAsgn.push_back(new_asgn);
            
            // Construct assignments with unassigned costs
            std::queue<std::tuple<std::vector<int>, int>> unassignedQueue;
            unassignedQueue.push(std::make_tuple(new_asgn, 0));
            while(!unassignedQueue.empty()){
              std::tuple<std::vector<int>, int> elem = unassignedQueue.front();
              unassignedQueue.pop();
              std::vector<int> elemVec = std::get<0>(elem);
              int elemInd = std::get<1>(elem);
              
              if(elemInd < min_dim){
                // Add unchanged assignment
                unassignedQueue.push(std::make_tuple(elemVec, elemInd+1));
                
                // Add new assignment (with unassigned cost)
                std::vector<int> new_unasgn(elemVec);
                new_unasgn.push_back(new_unasgn[elemInd]);
                new_unasgn[elemInd] = max_dim;
                unassignedQueue.push(std::make_tuple(new_unasgn, elemInd+1)); 
              }else{
                // A new assignment is made here            
                // TODO there is probably be a better way to avoid duplicates
                std::sort(elemVec.begin()+min_dim, elemVec.end());            
                std::string new_string = vectorToString(elemVec);
                auto it2 = std::find(unasgnAdded.begin(), unasgnAdded.end(), new_string);
                if (it2 == unasgnAdded.end()){
                  unasgnAdded.insert(vectorToString(elemVec));
                  allAsgn.push_back(elemVec);
                }//else
                  //std::cout << "collision " << vectorToString(elemVec) << std::endl;
              }
            }
            
          }else
            queue.push(new_asgn);
        }
      }
      
      return allAsgn;
    }
    
    //std::vector<std::vector<int>> AssignmentExhaustive::constructAssignments2(int min_dim, int max_dim){
    std::vector<std::vector<int>> constructAssignments2(int min_dim, int max_dim){
      std::vector<std::vector<int>> allAsgn;
      std::queue<std::tuple<std::vector<int>, std::vector<int>, bool>> queue;
      
      std::vector<int> asgn_temp, unasgn_temp;
      queue.push(std::make_tuple(asgn_temp, unasgn_temp, false));
      
      while(!queue.empty()){
        auto queueElem = queue.front();
        queue.pop();      
        std::vector<int> asgn = std::get<0>(queueElem);
        std::vector<int> unasgn = std::get<1>(queueElem);
        bool toUnasgn = std::get<2>(queueElem);
        int ind = asgn.size();
        for(int j = 0; j < max_dim+1; ++j){
          std::vector<int> new_asgn(asgn);
          std::vector<int> new_unasgn(unasgn);
          if(j < max_dim){
            // TODO avoid high complexity of find in vector
            std::vector<int>::iterator it = std::find(new_asgn.begin(), new_asgn.end(), j);
            std::vector<int>::iterator it2 = std::find(new_unasgn.begin(), new_unasgn.end(), j);
            if (it == new_asgn.end() and it2 == new_unasgn.end()){ // not found
              bool addIt = true;
              if(ind >= min_dim){
                for(unsigned int cpt = min_dim; cpt < new_asgn.size(); ++cpt)
                  if(new_asgn[cpt] > j)
                    addIt = false;
                for(unsigned int cpt = 0; cpt < new_unasgn.size(); ++cpt)
                  if(new_unasgn[cpt] > j)
                    addIt = false;
              }
              if(addIt){
                if(not toUnasgn)
                  new_asgn.push_back(j);
                else
                  new_unasgn.push_back(j);
              }else
                continue;
            }else
              continue;
            
            int new_ind = new_asgn.size();
            if(new_ind == max_dim){
              // A new assignement is made here
              for(auto new_unasgn_elem : new_unasgn)
                new_asgn.push_back(new_unasgn_elem);
              //std::sort(new_asgn.begin()+min_dim, new_asgn.end());
              allAsgn.push_back(new_asgn);
            }else
              queue.push(std::make_tuple(new_asgn, new_unasgn, false));
            
          }else{
            if(not toUnasgn and ind < min_dim){
              new_asgn.push_back(max_dim);
              queue.push(std::make_tuple(new_asgn, new_unasgn, true));
            }
          }
        }
      }
      
      return allAsgn;
    }
    
    //std::vector<std::vector<int>> AssignmentExhaustive::constructAssignments3(int min_dim, int max_dim){
    std::vector<std::vector<int>> constructAssignments3(int min_dim, int max_dim){
      std::vector<std::vector<int>> allAsgn;
      std::queue<std::tuple<std::vector<int>, std::vector<int>, bool,
                            std::vector<bool>, int>> queue;
      
      std::vector<int> asgn_temp, unasgn_temp;
      std::vector<bool> done_temp(max_dim, false);
      queue.push(std::make_tuple(asgn_temp, unasgn_temp, false, done_temp, -1));
      
      while(!queue.empty()){
        auto queueElem = queue.front();
        queue.pop();      
        std::vector<int> asgn = std::get<0>(queueElem);
        std::vector<int> unasgn = std::get<1>(queueElem);
        bool toUnasgn = std::get<2>(queueElem);
        std::vector<bool> done = std::get<3>(queueElem);
        int maxDone = std::get<4>(queueElem);
        
        int ind = asgn.size();
        for(int j = 0; j < max_dim+1; ++j){
          std::vector<int> new_asgn(asgn);
          std::vector<int> new_unasgn(unasgn);
          std::vector<bool> new_done(done);
          int new_maxDone = maxDone;
          if(j < max_dim){
            if(not done[j]){
              if(ind >= min_dim and maxDone > j)
                continue;
            }else
              continue;
            
            if(not toUnasgn){
              new_asgn.push_back(j);
              if(ind >= min_dim)
                new_maxDone = std::max(maxDone, j);
            }else{
              new_unasgn.push_back(j);
              new_maxDone = std::max(maxDone, j);
            }
            new_done[j] = true;          
            
            int new_ind = new_asgn.size();
            if(new_ind == max_dim){
              // A new assignement is made here
              for(auto new_unasgn_elem : new_unasgn)
                new_asgn.push_back(new_unasgn_elem);
              //std::sort(new_asgn.begin()+min_dim, new_asgn.end());
              allAsgn.push_back(new_asgn);
            }else
              queue.push(std::make_tuple(new_asgn, new_unasgn, false, new_done, new_maxDone));
            
          }else{
            if(not toUnasgn and ind < min_dim){
              new_asgn.push_back(max_dim);
              queue.push(std::make_tuple(new_asgn, new_unasgn, true, new_done, new_maxDone));
            }
          }
        }
      }
      
      return allAsgn;
    }

    inline void clear() {
      rowSize = 0;
      colSize = 0;
    }

    template <typename dataType>
    inline void clearMatrix() {
      std::vector<std::vector<dataType>> C
        = *((std::vector<std::vector<dataType>> *)Cptr);
      for(int r = 0, rS0 = rowSize; r < rS0; ++r)
        for(int c = 0, cS0 = colSize; c < cS0; ++c)
          C[r][c] = 0.0;
    }

    inline int setInput(int rowSize_, int colSize_, void *C_) {
      rowSize = rowSize_;
      colSize = colSize_;

      Cptr = C_;

      return 0;
    }
    
    template<class dataType>
    inline std::vector<std::vector<dataType>> getCostMatrix(){
      return *((std::vector<std::vector<dataType>> *)Cptr);
    }

    template<class dataType>
    void printTableVector(std::vector<std::vector<dataType>> &table){
      for(auto vecTemp : table){
        for(auto valTemp : vecTemp){
          std::cout << valTemp << " ";
        }
        std::cout << std::endl;
      }
      std::cout << " ------------------------------ " << std::endl;
    }
    
    void printAssignments(std::vector<std::vector<int>> &allAsgn){
      std::cout << " ------------------------------ " << std::endl;
      std::cout << "{";
      for(auto vecTemp : allAsgn){
        std::cout << "{";
        for(unsigned int i = 0; i < vecTemp.size(); ++i){
          auto valTemp = vecTemp[i];
          std::cout << valTemp;
          if(i != vecTemp.size() - 1) 
            std::cout << ",";
        }
        std::cout << "}, ";
      }
      std::cout << "}";
      std::cout << std::endl << " ------------------------------ " << std::endl;
    }
    
  private:
    void *Cptr;

    int rowSize = 0;
    int colSize = 0;
    
    // This attribute will store the computed assignments
    std::map<std::string, std::vector<std::vector<int>>> savedAsgn;
  };

  template <typename dataType>  
  dataType AssignmentExhaustive::tryAssignment(std::vector<int> &asgn, 
                                               std::vector<matchingTuple> &matchings){
    std::vector<std::vector<dataType>> costMatrix = getCostMatrix<dataType>();
    
    int nRows = costMatrix.size()-1;
    int nCols = costMatrix[0].size()-1;
    int max_dim = std::max(nRows, nCols);
    int min_dim = std::min(nRows, nCols);
    bool transpose = nRows > nCols;
    
    dataType cost = 0;
    for(int ind = 0; ind < asgn.size(); ++ind){
      int indMatrix = std::min(ind, min_dim);
      int i = (! transpose) ? indMatrix : asgn[ind];
      int j = (! transpose) ? asgn[ind] : indMatrix;
      cost += costMatrix[i][j];
      matchings.push_back(std::make_tuple(i, j, costMatrix[i][j]));
    }
    return cost;
  }
    
  template <typename dataType>
  int AssignmentExhaustive::run(std::vector<matchingTuple> &matchings){
    std::vector<std::vector<dataType>> costMatrix = getCostMatrix<dataType>();
    
    int nRows = costMatrix.size()-1;
    int nCols = costMatrix[0].size()-1;
    int max_dim = std::max(nRows, nCols);
    int min_dim = std::min(nRows, nCols);
    
    // --- Construct all possible assignments 
    std::vector<std::vector<int>> allAsgn;      
    // We hard write the basic assignments to avoid the call of constructAssignments
    // These assignments are, of course, automatically generated by constructAssignments
    if(min_dim == 1 and max_dim == 1)
      allAsgn = {{0}, {1,0}};
    else if(min_dim == 1 and max_dim == 2)
      allAsgn = {{0,1}, {2,0,1}, {1,0}};
    else if(min_dim == 1 and max_dim == 3)
      allAsgn = {{0,1,2}, {3,0,1,2}, {1,0,2}, {2,0,1}};
    else if(min_dim == 2 and max_dim == 2)
      allAsgn = {{0,1}, {0,2,1}, {2,1,0}, {2,2,0,1}, {1,0}, {1,2,0}, {2,0,1}};
    else{
      std::stringstream ss;
      ss << min_dim << "_" << max_dim;    
      std::string asgnName = ss.str();
      auto it = savedAsgn.find(asgnName);
      if (it == savedAsgn.end()){ // not found
        std::cout << asgnName << std::endl;
        //allAsgn = constructAssignments(min_dim, max_dim);
        //allAsgn = constructAssignments2(min_dim, max_dim);
        allAsgn = constructAssignments3(min_dim, max_dim);
        savedAsgn[asgnName] = allAsgn;
        std::cout << allAsgn.size() << " done" << std::endl;
      }else
        allAsgn = savedAsgn[asgnName];
    }
    
    // --- Try these assignments and get the better
    dataType bestCost = std::numeric_limits<dataType>::max();
    std::vector<matchingTuple> bestMatching;
    for(std::vector<int> asgn: allAsgn){
      std::vector<matchingTuple> tempMatching;
      dataType cost = tryAssignment<dataType>(asgn, tempMatching);
      if(bestCost > cost){
        bestCost = cost;
        bestMatching = tempMatching;
      }
    }
    matchings = bestMatching;
          
    return 0;
  }

} // namespace ttk

#endif
