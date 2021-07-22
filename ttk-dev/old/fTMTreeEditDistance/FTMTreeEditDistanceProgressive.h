/// \ingroup base
/// \class ttk::FTMTreeEditDistanceProgressive
/// \author Mathieu Pont <mathieu.pont@outlook.com>
///

#ifndef _FTMTREEEDITDISTANCEPROGRESSIVE_H
#define _FTMTREEEDITDISTANCEPROGRESSIVE_H

#include "FTMTreeEditDistanceBase.h"

namespace ttk {

  template <typename dataType>
  class FTMTreeEditDistanceProgressive : virtual public Debug, public FTMTreeEditDistanceBase {

  private:
    double addPercent = 0.5;
    
    std::vector<double> tree1Prices, tree2Prices;
    double t_assignment_time = 0;
    double epsilonDiviserMultiplier = 0;
    
  public:
    FTMTreeEditDistanceProgressive(){};

    ~FTMTreeEditDistanceProgressive(){};
    
    // ----------------------------------------
    // Tree Functions
    // ----------------------------------------
    void initTree(ftm::FTMTree_MT *tree, 
                  std::vector<std::tuple<SimplexId, SimplexId, dataType>> &treePairs,
                  std::vector<std::tuple<std::vector<ftm::idNode>, ftm::idNode>> &treeChildrenParent){
      for(int i = 0; i < treePairs.size(); ++i){
        auto pair = treePairs[i];
        SimplexId node1 = std::get<0>(pair);
        ftm::idNode parent1 = getParent(tree, node1);
        std::vector<ftm::idNode> children1 = getChildren(tree, node1);
        treeChildrenParent.push_back(std::make_tuple(children1, parent1));
        
        SimplexId node2 = std::get<1>(pair);
        ftm::idNode parent2 = getParent(tree, node2);
        std::vector<ftm::idNode> children2 = getChildren(tree, node2);
        treeChildrenParent.push_back(std::make_tuple(children2, parent2));
        
        //if(i < treePairs.size()/2){
          deleteNode(tree, node1);
          //if(node2 != tree->getNumberOfNodes()-1)
          deleteNode(tree, node2);
        //}
      }
    }
    
    template<class vecType>
    vecType getMinValueNotZero(std::vector<vecType> &vec){
      vecType minValue = std::numeric_limits<vecType>::max();
      bool foundOne = false;
      for(auto elem : vec)
        if(elem < minValue and elem != 0){
          minValue = elem;
          foundOne = true;
        }
      if(foundOne)
        return minValue;
      return 0;
    }
    
    template<class vecType>
    vecType getMaxValueNotZero(std::vector<vecType> &vec){
      vecType minValue = std::numeric_limits<vecType>::lowest();
      bool foundOne = false;
      for(auto elem : vec)
        if(elem > minValue and elem != 0){
          minValue = elem;
          foundOne = true;
        }
      if(foundOne)
        return minValue;
      return 0;
    }
    
    double initPrice(ftm::FTMTree_MT *tree, ftm::idNode node, bool isTree1){
      double newPrice = 0;
      auto childrens = getChildren(tree, node);
      for(auto children : childrens)
        newPrice += (isTree1) ? tree1Prices[children] : tree2Prices[children];
      newPrice /= (childrens.size() == 0) ? 1 : childrens.size();
      auto parent = getParent(tree, node);
      newPrice += (isTree1) ? tree1Prices[parent] : tree2Prices[parent];
      newPrice /= 2;
      newPrice += (isTree1) ? getMinValueNotZero(tree1Prices) : getMinValueNotZero(tree2Prices);
      //newPrice += (isTree1) ? getMaxValueNotZero(tree1Prices) : getMaxValueNotZero(tree2Prices);
      return newPrice;
    }
    
    void addNode(ftm::FTMTree_MT *tree, ftm::idNode node, bool isTree1,
                 std::tuple<std::vector<ftm::idNode>, ftm::idNode> &childrenParent){
      ftm::idNode parent = std::get<1>(childrenParent);
      if(node != parent)
        tree->makeSuperArc(node, parent);
      for(ftm::idNode children : std::get<0>(childrenParent)){
        parent = getParent(tree, children);
        if(children != parent){
          // Delete down arc from parent node
          for(ftm::idSuperArc i = 0; i < tree->getNode(parent)->getNumberOfDownSuperArcs(); ++i){
            ftm::idSuperArc arcId = tree->getNode(parent)->getDownSuperArcId(i);
            if(tree->getSuperArc(arcId)->getDownNodeId() == children)
              tree->getNode(parent)->removeDownSuperArc(arcId);
          } 
        }
        // Delete up arc from child nodes
        if(tree->getNode(children)->getNumberOfUpSuperArcs() != 0){
          ftm::idSuperArc arcId = tree->getNode(children)->getUpSuperArcId(0);
          tree->getNode(children)->removeUpSuperArc(arcId);          
        }
        // Add arc
        tree->makeSuperArc(children, node);
      }
      
      // Init price
      if(isTree1)
        tree1Prices[node] = initPrice(tree, node, isTree1);
      else
        tree2Prices[node] = initPrice(tree, node, isTree1);
    }
    
    void addPair(ftm::FTMTree_MT *tree, ftm::idNode node1, ftm::idNode node2, bool isTree1,
                 std::tuple<std::vector<ftm::idNode>, ftm::idNode> &childrenParent1,
                 std::tuple<std::vector<ftm::idNode>, ftm::idNode> &childrenParent2){
      addNode(tree, node1, isTree1, childrenParent1);
      //if(node2 != tree->getNumberOfNodes()-1)
      addNode(tree, node2, isTree1, childrenParent2);
    }
    
    void getNodesToUpdate(ftm::FTMTree_MT *tree, std::vector<ftm::idNode> nodes,
                          std::vector<ftm::idNode> &nodesToUpdate){
      std::vector<bool> nodesDone(tree->getNumberOfNodes(), false);
      std::stack<ftm::idNode> nodesStack;
      for(auto node : nodes)
        nodesStack.push(node);
      while(!nodesStack.empty()){
        ftm::idNode node = nodesStack.top();
        nodesStack.pop();
        nodesToUpdate.push_back(node);
        nodesDone[node] = true;
        ftm::idNode parent = getParent(tree, node);
        if(not nodesDone[parent]){
          nodesStack.push(parent);
        }
      }
      std::sort(nodesToUpdate.begin(), nodesToUpdate.end());
    }
    
    void updateTree(ftm::FTMTree_MT *tree, int iter, bool isTree1,
                    std::vector<std::tuple<SimplexId, SimplexId, dataType>> &treePairs,
                    std::vector<std::tuple<std::vector<ftm::idNode>, ftm::idNode>> &treeChildrenParent,
                    std::vector<ftm::idNode> &treeAddedNodes, std::vector<ftm::idNode> &nodesToUpdate){
      std::vector<ftm::idNode> nodesToSee;
      int size = std::max(treePairs.size()*addPercent, 1.);
      int noIter = (addPercent != 0) ? size : 1;
      if(addPercent != 0)
        iter *= size;
      for(int i = 0; i < noIter; ++i){
        if(iter < treePairs.size()){
          int treePairsIndex = treePairs.size()-1 - iter;
          ftm::idNode node1 = std::get<0>(treePairs[treePairsIndex]);
          ftm::idNode node2 = std::get<1>(treePairs[treePairsIndex]);
          int childrenParentIndex = treeChildrenParent.size()-1 - iter*2;
          
          addPair(tree, node1, node2, isTree1, treeChildrenParent[childrenParentIndex-1], 
                  treeChildrenParent[childrenParentIndex]);
          
          treeAddedNodes.push_back(node1);
          treeAddedNodes.push_back(node2);
          
          nodesToSee.push_back(node1);
        }
        iter += 1;
      }
      std::sort(treeAddedNodes.begin(), treeAddedNodes.end());
      getNodesToUpdate(tree, nodesToSee, nodesToUpdate);
    }
    
    // ----------------------------------------
    // Assignment Problem
    // ----------------------------------------        
    dataType forestAssignmentProblem(std::vector<std::vector<dataType>> &treeTable,
                                     std::vector<ftm::idNode> &children1,
                                     std::vector<ftm::idNode> &children2,
                                     std::vector<std::tuple<int, int>> &forestAssignment){
      // --- Create cost matrix
      int nRows = children1.size(), nCols = children2.size();
      std::vector<std::vector<dataType>> costMatrix(nRows+1, std::vector<dataType>(nCols+1));
      createCostMatrix(treeTable, children1, children2, costMatrix);
      //printTableVector(costMatrix);
      
      // --- Solve assignment problem
      std::vector<matchingTuple> matchings;
      
      // - Create prices
      std::vector<double> prices;
      for(int i = 0; i < children2.size(); ++i)
        prices.push_back(tree2Prices[children2[i]]);
      for(int i = 0; i < children1.size(); ++i)
        prices.push_back(tree1Prices[children1[i]]);
      
      nRows = costMatrix.size();
      nCols = costMatrix[0].size();
      
      AssignmentAuction<dataType> assignmentSolverAuction;
      assignmentSolverAuction.setInput(nRows, nCols, (void *)&costMatrix);
      assignmentSolverAuction.setBalanced(false);
      assignmentSolverAuction.setPrices(prices);
      assignmentSolverAuction.setNumberOfRounds(1);
      assignmentSolverAuction.setEpsilonDiviserMultiplier(epsilonDiviserMultiplier);
      assignmentSolverAuction.run(matchings);
      prices = assignmentSolverAuction.getPrices();
      
      /*Munkres assignmentSolverMunkres;
      assignmentSolverMunkres.setInput(nRows, nCols, (void *)&costMatrix);
      assignmentSolverMunkres.run<dataType>(matchings);*/
      
      // - Save prices
      for(int i = 0; i < prices.size(); ++i){
        if(i < children2.size())
          tree2Prices[children2[i]] = prices[i];
        else
          tree1Prices[children1[i-children2.size()]] = prices[i];
      }
      
      //printMatching(matchings);
      
      // --- Postprocess matching to create output assignment          
      dataType cost = postprocessAssignment<dataType>(matchings, children1, children2, forestAssignment);
      
      return cost;
    }
    
    void computeEquation13Progressive(int i, int j,
                     std::vector<std::vector<dataType>> &treeTable,
                     std::vector<std::vector<dataType>> &forestTable,
                     std::vector<std::vector<std::vector<std::tuple<int, int>>>> &forestBackTable,
                     std::vector<ftm::idNode> &children1, std::vector<ftm::idNode> &children2){
      if(children1.size() != 0 && children2.size() != 0){
        dataType forestTerm1, forestTerm2, forestTerm3;
        std::tuple<dataType, ftm::idNode> forestCoTerm1, forestCoTerm2;
        // Term 1
        forestCoTerm1 = computeTerm1_2<dataType>(children2, i, forestTable, true);                  
        forestTerm1 = forestTable[0][j] + std::get<0>(forestCoTerm1);
        
        // Term2
        forestCoTerm2 = computeTerm1_2<dataType>(children1, j, forestTable, false);
        forestTerm2 = forestTable[i][0] + std::get<0>(forestCoTerm2);
        
        // Term 3
        Timer t_assignment;
        std::vector<std::tuple<int, int>> forestAssignment;
        forestTerm3 = forestAssignmentProblem(treeTable, children1, children2, forestAssignment);
        t_assignment_time += t_assignment.getElapsedTime();
        
        // Compute table value
        forestTable[i][j] = std::min(std::min(forestTerm1, forestTerm2), forestTerm3);
        
        // Add backtracking information
        forestBackTable[i][j].clear();
        if(forestTable[i][j] == forestTerm3){
          forestBackTable[i][j] = forestAssignment;
        }else if(forestTable[i][j] == forestTerm2){
          forestBackTable[i][j].push_back(std::make_tuple(std::get<1>(forestCoTerm2), j));
        }else{
          forestBackTable[i][j].push_back(std::make_tuple(i, std::get<1>(forestCoTerm1)));
        }
      }else{
        // If one of the forest is empty we get back to equation 8 or 10
        forestTable[i][j] = (children1.size() == 0) ? forestTable[0][j] : forestTable[i][0];
      }
    }
    
    void globalEpsilonScaling(int iter, int maxIter){
      int middleIter = maxIter * 0.8;
      if(iter < middleIter ){
        epsilonDiviserMultiplier = 1/std::pow(5, middleIter - iter);
      }else{
        epsilonDiviserMultiplier = (epsilonDiviserMultiplier < 1) ? 1 : epsilonDiviserMultiplier*2;
      }
      std::cout << epsilonDiviserMultiplier << std::endl;
    }
    
    // ----------------------------------------
    // Main Functions
    // ----------------------------------------
    int computeEditDistanceProgressive(ftm::FTMTree_MT *tree1, ftm::FTMTree_MT *tree2,
                                std::vector<std::vector<dataType>> &treeTable,
                                std::vector<std::vector<dataType>> &forestTable,
                                std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
                                std::vector<std::vector<std::vector<std::tuple<int, int>>>> &forestBackTable
                                ){
      std::cout << "BEGIN PROGRESSIVE" << std::endl;
      
      Timer t_dyn;
      t_assignment_time = 0;
      
      tree1Prices.resize(tree1->getNumberOfNodes(), 0);
      tree2Prices.resize(tree2->getNumberOfNodes(), 0);
      
      std::vector<std::tuple<SimplexId, SimplexId, dataType>> tree1Pairs, tree2Pairs;
      tree1Pairs = computePersistencePairs<dataType>(tree1);
      tree2Pairs = computePersistencePairs<dataType>(tree2);
      //printPairs(tree1Pairs);
      //printPairs(tree2Pairs);      
      
      //tree1->printTree2();
      //tree2->printTree2();
      std::vector<std::tuple<std::vector<ftm::idNode>, ftm::idNode>> tree1ChildrenParent;
      std::vector<std::tuple<std::vector<ftm::idNode>, ftm::idNode>> tree2ChildrenParent;      
      initTree(tree1, tree1Pairs, tree1ChildrenParent);
      initTree(tree2, tree2Pairs, tree2ChildrenParent);
      //tree1->printTree2();
      //tree2->printTree2();
      
      std::vector<ftm::idNode> tree1AddedNodes, tree2AddedNodes;
      
      int maxIter = (addPercent != 0) ? 1/addPercent : std::max(tree1Pairs.size(), tree2Pairs.size());
      for(int iter = 0; iter < maxIter; ++iter){
        std::cout << iter << " / " << maxIter << std::endl;
        
        // Global epsilon-scaling
        globalEpsilonScaling(iter, maxIter);
        
        // Update Tree 1
        std::vector<ftm::idNode> tree1NodesToUpdate;
        updateTree(tree1, iter, true, tree1Pairs, tree1ChildrenParent, tree1AddedNodes, tree1NodesToUpdate);
          
        // Update Tree 2
        std::vector<ftm::idNode> tree2NodesToUpdate;
        updateTree(tree2, iter, false, tree2Pairs, tree2ChildrenParent, tree2AddedNodes, tree2NodesToUpdate);
          
        // Update Tables
        updateTables(tree1, tree2, tree1AddedNodes, tree2AddedNodes, tree1NodesToUpdate, tree2NodesToUpdate,
                     treeTable, forestTable, treeBackTable, forestBackTable);
        
        //
        /*dataType tree1MinValue = getMinValueNotZero(tree1Prices);
        for(int i = 0; i < tree1Prices.size(); ++i)
          if(tree1Prices[i] != 0)
            tree1Prices[i] -= tree1MinValue;
        dataType tree2MinValue = getMinValueNotZero(tree2Prices);
        for(int i = 0; i < tree2Prices.size(); ++i)
          if(tree2Prices[i] != 0)
            tree2Prices[i] -= tree2MinValue;*/
        
        //tree1->printTree2();
        //tree2->printTree2();
        /*printTableVector<dataType>(treeTable);
        printTableVector<dataType>(forestTable);*/
        
        //std::vector<std::tuple<ftm::idNode, ftm::idNode>> outputMatching;
        //computeMatching(treeBackTable, forestBackTable, outputMatching);
        //printOutputMatching<dataType>(outputMatching, tree1, tree2);
      }
      
      std::cout << "TIME DYNA.PROG. = " << t_dyn.getElapsedTime() << std::endl;
      std::cout << " - TIME ASSGNMT = " << t_assignment_time << std::endl;
      
      std::cout << "END PROGRESSIVE" << std::endl;
      return 0;
    }
    
    void updateTables(ftm::FTMTree_MT *tree1, ftm::FTMTree_MT *tree2,
                      std::vector<ftm::idNode> &tree1AddedNodes, std::vector<ftm::idNode> &tree2AddedNodes,
                      std::vector<ftm::idNode> &tree1NodesToUpdate, 
                      std::vector<ftm::idNode> &tree2NodesToUpdate,
                      std::vector<std::vector<dataType>> &treeTable,
                      std::vector<std::vector<dataType>> &forestTable,
                      std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
                      std::vector<std::vector<std::vector<std::tuple<int, int>>>> &forestBackTable){
      std::vector<ftm::idNode> tree2AddedNodesTemp;
      std::set_difference(tree2AddedNodes.begin(), tree2AddedNodes.end(), 
                          tree2NodesToUpdate.begin(), tree2NodesToUpdate.end(), 
                          std::inserter(tree2AddedNodesTemp, tree2AddedNodesTemp.begin()));

      for(auto node : tree1NodesToUpdate)
        processNode(tree1, tree2, node, true, tree1AddedNodes, tree2AddedNodesTemp, 
                    treeTable, forestTable, treeBackTable, forestBackTable);

      for(auto node : tree2NodesToUpdate)
        processNode(tree1, tree2, node, false, tree1AddedNodes, tree2AddedNodes,
                    treeTable, forestTable, treeBackTable, forestBackTable);
    }
    
    void processNode(ftm::FTMTree_MT *tree1, ftm::FTMTree_MT *tree2, ftm::idNode node, bool isTree1,
                     std::vector<ftm::idNode> &tree1AddedNodes, std::vector<ftm::idNode> &tree2AddedNodes,
                     std::vector<std::vector<dataType>> &treeTable,
                     std::vector<std::vector<dataType>> &forestTable,
                     std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
                     std::vector<std::vector<std::vector<std::tuple<int, int>>>> &forestBackTable){
      //std::cout << ((isTree1)?"tree1 ":"tree2 ") << node << std::endl;
      if(isTree1){
        int i = node+1;
        // --- Equation 8
        computeEquation8(tree1, node, i, treeTable, forestTable);
        
        // --- Equation 9
        computeEquation9(tree1, node, i, treeTable, forestTable);
        
        // --- Equation 12 and 13
        std::vector<ftm::idNode> children1 = getChildren(tree1, node);
        for(auto nodeJ : tree2AddedNodes){
          int j = nodeJ+1;
          std::vector<ftm::idNode> children2 = getChildren(tree2, nodeJ);
          // --- Equation 13
          computeEquation13Progressive(i, j, treeTable, forestTable, forestBackTable, children1, children2);
          
          // --- Equation 12
          computeEquation12(tree1, tree2, i, j, node, nodeJ, treeTable, 
                            forestTable, treeBackTable, children1, children2);
        }
      }else{
        int j = node+1;
        // --- Equation 10
        computeEquation10(tree2, node, j, treeTable, forestTable);
        
        // --- Equation 11
        computeEquation11(tree2, node, j, treeTable, forestTable);
        
        // --- Equation 12 and 13
        std::vector<ftm::idNode> children2 = getChildren(tree2, node);
        for(auto nodeI : tree1AddedNodes){
          int i = nodeI+1;
          std::vector<ftm::idNode> children1 = getChildren(tree1, nodeI);
          // --- Equation 13
          computeEquation13Progressive(i, j, treeTable, forestTable, forestBackTable, children1, children2);
          
          // --- Equation 12
          computeEquation12(tree1, tree2, i, j, nodeI, node, treeTable, 
                            forestTable, treeBackTable, children1, children2);
        }
      }
    }
    
    // ----------------------------------------
    // Utils
    // ----------------------------------------
    
  };

} // namespace ttk

#endif
