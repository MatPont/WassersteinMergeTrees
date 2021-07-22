/// \ingroup base
/// \class ttk::FTMTreeEditDistanceBase
/// \author Mathieu Pont <mathieu.pont@outlook.com>
///
/// This class is aim to be abstract (contains a pure virtual function)
/// The class to instance is FTMTreeEditDistance
/// TODO make "forestAssignmentProblem" or "computeEquation13" as pure virtual function

#ifndef _FTMTREEEDITDISTANCEBASE_H
#define _FTMTREEEDITDISTANCEBASE_H

#include <FTMTree.h>
#include <FTMNode.h>

#include "FTMTreeUtils.h"

using namespace ttk;

class FTMTreeEditDistanceBase : virtual public Debug {
  protected:
    int AssignmentSolver = 0;
    double EpsilonTree1 = 0;
    double EpsilonTree2 = 0;
    double Epsilon2Tree1 = 0;
    double Epsilon2Tree2 = 0;
    double PersistenceThreshold = 0;
    bool Parallelize = true;
    int NumberOfThreads = 0;
    bool BarycenterMergeTree = false;
    bool UseMinMaxPair = true;
    
    bool BranchDecomposition = false;
    int WassersteinPower = 2;
    bool NormalizedWasserstein = false;
    bool RescaledWasserstein = false;
    double NormalizedWassersteinReg = 0.;
    bool KeepSubtree = true;
    
    bool ProgressiveComputation = false;
    bool DistanceSquared = true;
    
    int verbose = 1;

  public:
    void setAssignmentSolver(int assignmentSolver){
      AssignmentSolver = assignmentSolver;
    }
    
    void setEpsilonTree1(double epsilon){
      EpsilonTree1 = epsilon;
    }
    
    void setEpsilonTree2(double epsilon){
      EpsilonTree2 = epsilon;
    }
    
    void setEpsilon2Tree1(double epsilon){
      Epsilon2Tree1 = epsilon;
    }
    
    void setEpsilon2Tree2(double epsilon){
      Epsilon2Tree2 = epsilon;
    }
    
    void setPersistenceThreshold(double pt){
      PersistenceThreshold = pt;
    }

    void setParallelize(bool para){
      Parallelize = para;
    }
    
    void setNumberOfThreads(int noThreads){
      NumberOfThreads = noThreads;
      if(NumberOfThreads == 0){
        NumberOfThreads = this->threadNumber_;
        std::cout << "no threads : " << NumberOfThreads << std::endl;
      }
    }
    
    void setBranchDecomposition(bool useBD){
      BranchDecomposition = useBD;
    }
    
    void setNormalizedWasserstein(bool normalizedWasserstein){
      NormalizedWasserstein = normalizedWasserstein;
    }
    
    void setRescaledWasserstein(bool rescaledWasserstein){
      RescaledWasserstein = rescaledWasserstein;
    }
    
    void setNormalizedWassersteinReg(double normalizedWassersteinReg){
      NormalizedWassersteinReg = normalizedWassersteinReg;
    }
    
    void setKeepSubtree(bool keepSubtree){
      KeepSubtree = keepSubtree;
    }
    
    void setProgressiveComputation(bool progressive){
      ProgressiveComputation = progressive;
      /*if(ProgressiveComputation)
        Preprocess = false;*/
    }
    
    void setBarycenterMergeTree(bool imt){
      BarycenterMergeTree = imt;
    }
    
    void setDistanceSquared(bool distanceSquared){
      DistanceSquared = distanceSquared;
    }
    
    void setUseMinMaxPair(bool useMinMaxPair){
      UseMinMaxPair = useMinMaxPair;
    }
    
    void setVerbose(int verb){
      verbose = verb;
    }

    // ----------------------------------------
    // Tree Preprocessing
    // ----------------------------------------
    template <class dataType>
    void mergeSaddle(ftm::FTMTree_MT *tree, double epsilon, 
                      std::vector<std::vector<ftm::idNode>> &treeNodeMerged,
                      bool mergeByPersistence = false){
      bool fullMerge = (epsilon == 100);
      
      if(mergeByPersistence)
        computePersistencePairs<dataType>(tree); // need to have the pairing (if merge by persistence)
      
      // Compute epsilon value
      dataType maxValue = tree->getValue<dataType>(0);
      dataType minValue = tree->getValue<dataType>(0);
      for(int i = 0; i < tree->getNumberOfNodes(); ++i){
        if(! isRoot(tree, i) and ! isLeaf(tree, i)){
          dataType iValue = tree->getValue<dataType>(i);
          if(mergeByPersistence){
            maxValue = (maxValue < iValue) ? iValue : maxValue;
            minValue = (minValue > iValue) ? iValue : minValue;
          }else{
            ftm::idNode parent = getParent(tree, i);
            dataType parentValue = tree->getValue<dataType>(parent);
            dataType tempMax = std::max(iValue, parentValue);
            dataType tempMin = std::min(iValue, parentValue);
            if((tempMax - tempMin) > (maxValue - minValue)){
              maxValue = tempMax;
              minValue = tempMin;
            }
          }
        }
      }
      epsilon = (maxValue - minValue) * epsilon / 100;   
      
      // Merge saddle
      std::stack<int> nodeStack;
      for(int i = tree->getNumberOfNodes()-1; i >= 0; --i)
        nodeStack.push(i);
      while(! nodeStack.empty()){
        ftm::idNode nodeId = nodeStack.top();
        nodeStack.pop();
        if(! isRoot(tree, nodeId) and ! isLeaf(tree, nodeId)){
          ftm::idNode parentNodeId = getParent(tree, nodeId);
          dataType nodeValue = tree->getValue<dataType>(nodeId);
          dataType parentNodeValue = tree->getValue<dataType>(parentNodeId);
          dataType diffValue = std::max(nodeValue, parentNodeValue) - 
                                std::min(nodeValue, parentNodeValue);
          if(diffValue <= epsilon){
            ftm::idNode nodeIdToDelete, nodeIdToKeep;
            if(mergeByPersistence){
              auto birthDeath1 = getBirthDeath<dataType>(tree, nodeId);
              auto birthDeath2 = getBirthDeath<dataType>(tree, parentNodeId);
              dataType pers1 = std::get<1>(birthDeath1) - std::get<0>(birthDeath1);
              dataType pers2 = std::get<1>(birthDeath2) - std::get<0>(birthDeath2);
              nodeIdToDelete = (pers1 > pers2) ? parentNodeId : nodeId;
              nodeIdToKeep = (pers1 > pers2) ? nodeId : parentNodeId;
              if(nodeIdToDelete == parentNodeId)
                nodeStack.push(nodeId);
              //parentNodeId = (pers1 > pers2) ? nodeId : parentNodeId;                    
            }else{
              nodeIdToDelete = nodeId;
              nodeIdToKeep = parentNodeId;
            }
            // Manage nodeMerged vector of vector
            for(auto node : treeNodeMerged[nodeIdToDelete])
              treeNodeMerged[nodeIdToKeep].push_back(node);
            treeNodeMerged[nodeIdToKeep].push_back(tree->getNode(nodeIdToDelete)->getOrigin());
            treeNodeMerged[nodeIdToDelete].clear();
            // Delete node
            deleteNode(tree, nodeIdToDelete); //tree->delNode(nodeIdToDelete);
          }
        }
      }
      
      // Reset pairing
      /*if(mergeByPersistence){
        for(int i = 0; i < tree->getNumberOfNodes(); ++i)
          tree->getNode(i)->setOrigin(ftm::nullNodes);
      }*/
      
      if(fullMerge){
        auto root = getRoot(tree);
        tree->getNode(root)->setOrigin(root);
      }
    }
    
    template <class dataType>
    void persistenceMerging(ftm::FTMTree_MT *tree, double epsilon){
      bool fullMerge = (epsilon == 100);
      epsilon /= 100;
      dataType maxPers = getMaximumPersistence<dataType>(tree);
      
      std::queue<ftm::idNode> queue;
      queue.push(getRoot(tree));
      while(!queue.empty()){
        ftm::idNode node = queue.front();
        queue.pop();
        if(!isRoot(tree, node))
          if(getNodePersistence<dataType>(tree, node) / maxPers < epsilon)
            setParent(tree, node, getParent(tree, getParent(tree, node)));
        for(auto child : getChildren(tree, node))
          queue.push(child);
      }
      
      if(fullMerge){
        auto root = getRoot(tree);
        if(tree->getNode(root)->getOrigin() != root){
          setParent(tree, tree->getNode(root)->getOrigin(), root);
          tree->getNode(root)->setOrigin(root);
        }
      }
    }
    
    template <class dataType>
    void persistenceThresholding(ftm::FTMTree_MT *tree, double persistenceThresholdT){
      dataType threshold = persistenceThresholdT / 100 * getMaximumPersistence<dataType>(tree);
      for(int i = 0; i < tree->getNumberOfNodes(); ++i){
        dataType nodePers = getNodePersistence<dataType>(tree, i);
        if((nodePers == 0 or nodePers <= threshold or not isNodeOriginDefined(tree, i)) and !isRoot(tree, i)){
//std::cout << "delete " << i << " " << tree->getNode(i)->getOrigin() << std::endl;
          deleteNode(tree, i);
          if(isNodeOriginDefined(tree, i) and tree->getNode(tree->getNode(i)->getOrigin())->getOrigin() == i)
            deleteNode(tree, tree->getNode(i)->getOrigin());
        }
      }
//std::cout << "==================== " << tree->getNumberOfNodes() << std::endl;
    }
    
    template <class dataType>
    void persistenceThresholding(ftm::FTMTree_MT *tree){
      persistenceThresholding<dataType>(tree, PersistenceThreshold);
    }
    
    template <class dataType>
    void preprocessTree(ftm::FTMTree_MT *tree, double epsilon, 
                        std::vector<std::vector<ftm::idNode>> &treeNodeMerged){
      // Manage inconsistent critical points
      // Critical points with same scalar value than parent
      for(int i = 0; i < tree->getNumberOfNodes(); ++i)
        if(!isNodeAlone(tree, i) and !isRoot(tree, i) and 
            tree->getValue<dataType>(getParent(tree, i)) == tree->getValue<dataType>(i))
          deleteNode(tree, i);
      // Valence 2 nodes
      for(int i = 0; i < tree->getNumberOfNodes(); ++i)
        if(tree->getNode(i)->getNumberOfUpSuperArcs() == 1 and 
            tree->getNode(i)->getNumberOfDownSuperArcs() == 1)
          deleteNode(tree, i);
      
      // Compute persistence pairs
      auto pairs = computePersistencePairs<dataType>(tree);
      //printPairs<dataType>(pairs);
        
      // Delete null persistence pairs and persistence thresholding 
      persistenceThresholding<dataType>(tree);
      
      // Merge saddle points according epsilon
      if(epsilon != 0)
        mergeSaddle<dataType>(tree, epsilon, treeNodeMerged);
    }
    
    template <class dataType>
    ftm::FTMTree_MT* computeBranchDecomposition(ftm::FTMTree_MT *tree,
                                                std::vector<std::vector<ftm::idNode>> &treeNodeMerged){
      ftm::FTMTree_MT *treeNew = copyTree<dataType>(tree);
      
      ftm::idNode root = getRoot(treeNew);          
      
      // Manage when there is only one pair
      if(isThereOnlyOnePersistencePair(tree)){
        ftm::idNode rootOrigin = treeNew->getNode(root)->getOrigin();
        treeNew->getNode(rootOrigin)->setOrigin(rootOrigin);
        return treeNew;
      }
      
      // Manage multi persistence pairing
      std::vector<std::vector<ftm::idNode>> treeMultiPers(tree->getNumberOfNodes());
      for(int i = 0; i < tree->getNumberOfNodes(); ++i)
        if(isLeaf(tree, i) and isNodeOriginDefined(tree, i) and not isNodeAlone(tree, i) and
            tree->getNode(tree->getNode(i)->getOrigin())->getOrigin() != i)
          treeMultiPers[tree->getNode(i)->getOrigin()].push_back(i);
        
      // General case
      std::vector<bool> nodeDone(tree->getNumberOfNodes(), false);
      std::queue<ftm::idNode> queueNodes;
      queueNodes.push(root);
      while(!queueNodes.empty()){
        ftm::idNode node = queueNodes.front();
        queueNodes.pop();
        ftm::idNode nodeOrigin = treeNew->getNode(node)->getOrigin();
        if(node == nodeOrigin or node < nodeOrigin)
          continue;
        
        // Init vector with all origins
        std::vector<std::tuple<ftm::idNode, int>> vecOrigins;
        for(auto nodeMergedOrigin : treeNodeMerged[node]){
          vecOrigins.push_back(std::make_tuple(nodeMergedOrigin, 0));
          for(auto multiPersOrigin : treeMultiPers[tree->getNode(nodeMergedOrigin)->getOrigin()])
            vecOrigins.push_back(std::make_tuple(multiPersOrigin, 1));
        }
        if(not isNodeMerged(tree, node))
          for(auto multiPersOrigin : treeMultiPers[node])
            vecOrigins.push_back(std::make_tuple(multiPersOrigin, 1));
        vecOrigins.push_back(std::make_tuple(nodeOrigin, 2));
        
        bool splitRoot = false and (vecOrigins.size() != 1 and isRoot(treeNew, node));
        
        // Process each origin
        for(auto stackTuple : vecOrigins){
          ftm::idNode nodeOriginT = std::get<0>(stackTuple);
          int nodeOriginTID = std::get<1>(stackTuple);
          if(nodeDone[nodeOriginT] and nodeDone[tree->getNode(nodeOriginT)->getOrigin()])
            continue;
          nodeDone[nodeOriginT] = true;
          nodeDone[tree->getNode(nodeOriginT)->getOrigin()] = true;
          
          // Manage new parent
          ftm::idNode newParent = node;
          // - if merged node
          if(nodeOriginTID == 0){
            newParent = treeNew->getNode(nodeOriginT)->getOrigin();
            setParent(treeNew, newParent, getParent(treeNew, node));
          // - if multi pers node or nodeOrigin and splitRoot
          }else if(nodeOriginTID == 1 or (nodeOriginTID == 2 and splitRoot)){ 
            newParent = nodeOriginT;
          }
          /*if(nodeOriginT != nodeOrigin){ //if(nodeOriginTID != 2)
            newParent = (nodeOriginTID == 0) ? treeNew->getNode(nodeOriginT)->getOrigin() : nodeOriginT;
            if(nodeOriginTID != 1) // if not a multi pers node
              setParent(treeNew, newParent, getParent(treeNew, node));
          }*/
          
          // Set nodes in the branch as childrens of the node
          ftm::idNode parentNodeOrigin = getParent(treeNew, nodeOriginT);
          while(parentNodeOrigin != node){
            ftm::idNode oldParentNodeOrigin = getParent(treeNew, parentNodeOrigin);
            setParent(treeNew, parentNodeOrigin, newParent);
            parentNodeOrigin = oldParentNodeOrigin;
          } 
          
          if(nodeOriginTID == 1 or (nodeOriginTID == 2 and splitRoot))
            setParent(treeNew, newParent, getParent(treeNew, node));
          else //if(nodeOriginTID != 1) // if not a multi pers node
            // Delete the other node of the pair
            deleteNode(treeNew, nodeOriginT);
          if(nodeOriginTID == 2 and splitRoot)
            tree->getNode(node)->setOrigin(node);
          
          // Push childrens of the node to the stack to process them
          std::vector<ftm::idNode> childrenNode = getChildren(treeNew, newParent);
          for(ftm::idNode children : childrenNode)
            if(! isLeaf(treeNew, children))
              queueNodes.push(children);
        }
      }
      
      // Verify inconsistency
      //verifyBranchDecompositionInconsistency<dataType>(treeNew);
      
      return treeNew;
    }
    
    template <class dataType>
    void dontUseMinMaxPair(ftm::FTMTree_MT *tree){
      ftm::idNode treeRoot = getRoot(tree);
      if(tree->getNode(treeRoot)->getOrigin() == treeRoot){
        ftm::idNode treeRootChildren0 = getChildren(tree, treeRoot)[0];
        bool isJt = (tree->getValue<dataType>(treeRootChildren0) < tree->getValue<dataType>(treeRoot));
        ftm::idNode nodeIdToDelete = treeRoot;
        for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i)
          if(tree->getNode(i)->getOrigin() == (int)treeRoot and i != treeRoot)
            if((isJt and tree->getValue<dataType>(i) < tree->getValue<dataType>(nodeIdToDelete)) or 
              (not isJt and tree->getValue<dataType>(i) > tree->getValue<dataType>(nodeIdToDelete)))
              nodeIdToDelete = i;
        if(nodeIdToDelete != treeRoot){
          if(isThereOnlyOnePersistencePair(tree))
            tree->getNode(nodeIdToDelete)->setOrigin(nodeIdToDelete);
          else
            deleteNode(tree, nodeIdToDelete);
        }
      }else{
        ftm::idNode rootOrigin = tree->getNode(treeRoot)->getOrigin();
        if(isThereOnlyOnePersistencePair(tree))
          tree->getNode(rootOrigin)->setOrigin(rootOrigin);
        else
          deleteNode(tree, rootOrigin);
      }
      
      tree->getNode(treeRoot)->setOrigin(treeRoot);
    }
    
    void verifyPairsTree(ftm::FTMTree_MT *tree){
      int cptBug = 0;
      for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i){
        //if((unsigned int)tree->getNode(i)->getOrigin() >= tree->getNumberOfNodes()){
        if(not isNodeOriginDefined(tree, i)){
          std::cout << i << " _ " << tree->getNode(i)->getOrigin() << " / " << tree->getNumberOfNodes() << std::endl;
          if(isNodeAlone(tree, i))
            std::cout << "alone" << std::endl;
          cptBug++;
        }
      }
      std::cout << cptBug << std::endl;
      myPause();
    }

    template<class dataType>
    ftm::FTMTree_MT* preprocessingPipeline(ftm::FTMTree_MT *&tree, double epsilonTree, double epsilon2Tree,
                                           int verboseT=1){
      Timer t_proc;
      std::vector<std::vector<ftm::idNode>> treeNodeMerged(tree->getNumberOfNodes());
      
      //tree->printTree2(); 
      preprocessTree<dataType>(tree, epsilonTree, treeNodeMerged);
      //tree->printTree2(); 
      //printPairsFromTree<dataType>(tree, BranchDecomposition); 

      //verifyPairsTree(tree);
      ftm::FTMTree_MT *treeOld = tree;
      if(BranchDecomposition)
        tree = computeBranchDecomposition<dataType>(tree, treeNodeMerged);
      //tree->printTree2(); 
      
      //verifyPairsTree(tree);
      if(not UseMinMaxPair)
        dontUseMinMaxPair<dataType>(tree);
      
      persistenceMerging<dataType>(tree, epsilon2Tree);
      
      //verifyPairsTree(tree);
      auto t_preproc_time = t_proc.getElapsedTime();
      if(verboseT > 0)
        std::cout << "TIME PREPROC.   = " << t_preproc_time << std::endl;
      
      //tree->printTree2(); 
      //printTreeScalars<dataType>(tree); 
      //printPairsFromTree<dataType>(tree, BranchDecomposition); 

      return treeOld;
    }
    
    // ----------------------------------------
    // Tree Postprocessing
    // ----------------------------------------
    // TODO manage full merge and dontUseMinMaxPair
    template <class dataType>
    void branchDecompositionToTree(ftm::FTMTree_MT *tree){
      bool isJT = isJoinTree<dataType>(tree);
      std::vector<std::tuple<ftm::idNode, ftm::idNode>> nodeParent;
      std::queue<ftm::idNode> queue;
      queue.push(getRoot(tree));
      while(!queue.empty()){
        ftm::idNode node = queue.front();
        queue.pop();
        if(isLeaf(tree, node)){
          if(isNodeAlone(tree, tree->getNode(node)->getOrigin()))
            nodeParent.push_back(std::make_tuple(tree->getNode(node)->getOrigin(), node));
          continue;
        }
        // Get children and sort them by scalar values
        std::vector<ftm::idNode> children = getChildren(tree, node);
        std::vector<std::tuple<ftm::idNode, dataType>> childrenScalars(children.size());
        for(int i = 0; i < children.size(); ++i)
          childrenScalars.push_back(std::make_tuple(children[i], tree->getValue<dataType>(children[i])));
        auto comp = [&](const std::tuple<ftm::idNode, dataType>& a, 
                        const std::tuple<ftm::idNode, dataType>& b){ 
                        return isJT ? std::get<1>(a) < std::get<1>(b) : std::get<1>(a) > std::get<1>(b); };
        std::sort(std::begin(childrenScalars), std::end(childrenScalars), comp);
        // Get new parent of children
        bool multiPersPair = tree->getNode(tree->getNode(node)->getOrigin())->getOrigin() != node;
        auto nodeOrigin = tree->getNode(node)->getOrigin();
        for(int i = 1; i < children.size(); ++i)
          nodeParent.push_back(std::make_tuple(children[i], children[i-1]));
        if(not multiPersPair)
          nodeParent.push_back(std::make_tuple(nodeOrigin, children[children.size()-1]));
        else{
          nodeParent.push_back(std::make_tuple(children[0], nodeOrigin));
          nodeParent.push_back(std::make_tuple(node, children[children.size()-1]));
        }
        // Push children to the queue
        for(auto child : children)
          queue.push(child);
      }
      // Set new parents for each node
      for(auto nodeParentT : nodeParent)
        setParent(tree, std::get<0>(nodeParentT), std::get<1>(nodeParentT));
      
      // TODO Manage merge of multi pers pairs
      /*std::queue<ftm::idNode> queue;
      queue.push(getRoot(tree));
      while(!queue.empty()){
        ftm::idNode node = queue.front();
        queue.pop();
        if(!isRoot(tree, node) and 
           tree->getValue<dataType>(node) == tree->getValue<dataType>(getParent(tree, node))){
          
        }
        // Push children to the queue
        std::vector<ftm::idNode> children = getChildren(node);
        for(child : children)
          queue.push(child);
      }*/
    }
    
    template<class dataType>
    void postprocessingPipeline(ftm::FTMTree_MT *&tree){
      if(BranchDecomposition) 
        branchDecompositionToTree<dataType>(tree);
    }
    
    // ----------------------------------------
    // Output Matching
    // ----------------------------------------
    void computeMatching(std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
                         std::vector<std::vector<std::vector<std::tuple<int, int>>>> &forestBackTable,
                         std::vector<std::tuple<ftm::idNode, ftm::idNode>> &outputMatching,
                         int startR, int startC){
      std::queue<std::tuple<int, int, bool>> backQueue;
      backQueue.push(std::make_tuple(startR, startC, true));
      while(! backQueue.empty()){
        std::tuple<int, int, bool> elem = backQueue.front();
        backQueue.pop();
        bool useTreeTable = std::get<2>(elem);
        int i = std::get<0>(elem);
        int j = std::get<1>(elem);
        
        if(useTreeTable){
          int tupleI = std::get<0>(treeBackTable[i][j]);
          int tupleJ = std::get<1>(treeBackTable[i][j]);                
          if(tupleI != 0 && tupleJ != 0){
            useTreeTable = (tupleI != i || tupleJ != j);
            backQueue.push(std::make_tuple(tupleI, tupleJ, useTreeTable));
            if(not useTreeTable){ // We have matched i and j
              outputMatching.push_back(std::make_tuple(tupleI-1, tupleJ-1));
              //std::cout << "T(" << i << ", " << j << ") " << tupleI << " - " << tupleJ << std::endl;
            }
          }
        }else{
          for(std::tuple<int, int> forestBackElem : forestBackTable[i][j]){
            int tupleI = std::get<0>(forestBackElem);
            int tupleJ = std::get<1>(forestBackElem);
            if(tupleI != 0 && tupleJ != 0){
              useTreeTable = (tupleI != i && tupleJ != j);
              backQueue.push(std::make_tuple(tupleI, tupleJ, useTreeTable));
              //std::cout << "F(" << i << ", " << j << ") " << tupleI << " - " << tupleJ << std::endl;
            }
          }
        }
      }
    }
    
    void convertBranchDecompositionMatching(ftm::FTMTree_MT *tree1, ftm::FTMTree_MT *tree2,
                         std::vector<std::tuple<ftm::idNode, ftm::idNode>> &outputMatching){
      std::vector<std::tuple<ftm::idNode, ftm::idNode>> toAdd;
      for(std::tuple<ftm::idNode, ftm::idNode> mTuple : outputMatching){
        ftm::idNode node1 = std::get<0>(mTuple);
        ftm::idNode node2 = std::get<1>(mTuple);
        ftm::idNode node1Origin = tree1->getNode(node1)->getOrigin();
        ftm::idNode node2Origin = tree2->getNode(node2)->getOrigin();
        ftm::idNode node1OriginOrigin = tree1->getNode(node1Origin)->getOrigin();
        ftm::idNode node2OriginOrigin = tree2->getNode(node2Origin)->getOrigin();
        if(! isNodeAlone(tree1, node1Origin) and ! isNodeAlone(tree2, node2Origin) and
           (node1OriginOrigin == node1 or node2OriginOrigin == node2) )
          toAdd.push_back(std::make_tuple(node1Origin, node2Origin));
      }
      outputMatching.insert(outputMatching.end(), toAdd.begin(), toAdd.end());
    }
    
    template<class dataType>
    void identifyRealMatching(ftm::FTMTree_MT *tree1, ftm::FTMTree_MT *tree2,
                         std::vector<std::tuple<ftm::idNode, ftm::idNode>> &outputMatching,
                         std::vector<std::tuple<ftm::idNode, ftm::idNode, bool>> &realMatching){
      for(std::tuple<ftm::idNode, ftm::idNode> mTuple : outputMatching){
        ftm::idNode tree1Node = std::get<0>(mTuple);
        ftm::idNode tree2Node = std::get<1>(mTuple);
        dataType relabelCostVal = relabelCostOnly<dataType>(tree1, tree1Node, tree2, tree2Node);
        dataType deleteInsertCostVal = deleteCost<dataType>(tree1, tree1Node) + 
                                       insertCost<dataType>(tree2, tree2Node);
        //std::cout << tree1Node << " _ " << tree2Node << " _ " << relabelCostVal << " _ " << deleteInsertCostVal << std::endl;
        bool isRealMatching = (relabelCostVal <= deleteInsertCostVal);
        realMatching.push_back(std::make_tuple(tree1Node, tree2Node, isRealMatching));
      }
    }
    
    // ----------------------------------------
    // Edit Costs
    // ----------------------------------------    
    template <class dataType>
    dataType computeDistance(dataType x1, dataType x2, dataType y1, dataType y2, double power=2){
      if(power <= 0)
        return std::max(abs(x1 - y1), abs(x2 - y2));
      else
        return std::pow(abs(x1 - y1), power) + std::pow(abs(x2 - y2), power);
    }

    template <class dataType>
    dataType regularizationCost(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
      dataType shiftMin1 = getMinMaxLocal<dataType>(tree, nodeId);
      dataType shiftMax1 = getMinMaxLocal<dataType>(tree, nodeId, false);
      dataType projec = (shiftMin1+shiftMax1)/2;
      return computeDistance<dataType>(shiftMin1, 0, shiftMax1, 0, WassersteinPower);
    }
    
    template <class dataType>
    dataType deleteCost(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
      dataType cost = 0;
      dataType newMin=0.0, newMax=1.0;
      if(NormalizedWasserstein and RescaledWasserstein){
        auto newMinMax = getNewMinMax<dataType>(tree, nodeId, tree, nodeId);
        newMin = std::get<0>(newMinMax);
        newMax = std::get<1>(newMinMax);
      }
      // Get birth/death
      auto birthDeath = NormalizedWasserstein ? getNormalizedBirthDeath<dataType>(tree, nodeId, newMin, 
                                                                                  newMax) : 
                                                getBirthDeath<dataType>(tree, nodeId);
      dataType birth = std::get<0>(birthDeath);
      dataType death = std::get<1>(birthDeath);
      dataType projec = (birth+death)/2;
      // Compute delete cost
      cost = computeDistance<dataType>(birth, death, projec, projec, WassersteinPower);
      // Divide cost by two if not branch decomposition and not merged
      /*if(! BranchDecomposition and ! isNodeMerged(tree, nodeId))
        cost /= 2;*/
      // Regularize
      if(NormalizedWasserstein and NormalizedWassersteinReg != 0)
        cost += NormalizedWassersteinReg * regularizationCost<dataType>(tree, nodeId);
      
      return cost;
    }

    template <class dataType>
    dataType insertCost(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
      return deleteCost<dataType>(tree, nodeId);
    }
    
    template <class dataType>
    dataType regularizationCost2(ftm::FTMTree_MT *tree1, ftm::idNode nodeId1,
                                 ftm::FTMTree_MT *tree2, ftm::idNode nodeId2){
      dataType shiftMin1 = getMinMaxLocal<dataType>(tree1, nodeId1);
      dataType shiftMax1 = getMinMaxLocal<dataType>(tree1, nodeId1, false);
      dataType shiftMin2 = getMinMaxLocal<dataType>(tree2, nodeId2);
      dataType shiftMax2 = getMinMaxLocal<dataType>(tree2, nodeId2, false);
      return computeDistance<dataType>(shiftMin1, shiftMax1, shiftMin2, shiftMax2, WassersteinPower);
    }

    template <class dataType>
    dataType relabelCostOnly(ftm::FTMTree_MT *tree1, ftm::idNode nodeId1,
                             ftm::FTMTree_MT *tree2, ftm::idNode nodeId2){
      dataType cost = 0;
      dataType newMin=0.0, newMax=1.0;
      if(NormalizedWasserstein and RescaledWasserstein){
        auto newMinMax = getNewMinMax<dataType>(tree1, nodeId1, tree2, nodeId2);
        newMin = std::get<0>(newMinMax);
        newMax = std::get<1>(newMinMax);
      }
      // Get birth/death of the first tree
      auto birthDeath1 = NormalizedWasserstein ? getNormalizedBirthDeath<dataType>(tree1, nodeId1, newMin, 
                                                                                   newMax) : 
                                                 getBirthDeath<dataType>(tree1, nodeId1);
      dataType birth1 = std::get<0>(birthDeath1);
      dataType death1 = std::get<1>(birthDeath1);
      // Get birth/death of the second tree
      auto birthDeath2 = NormalizedWasserstein ? getNormalizedBirthDeath<dataType>(tree2, nodeId2, newMin, 
                                                                                   newMax) : 
                                                 getBirthDeath<dataType>(tree2, nodeId2);
      dataType birth2 = std::get<0>(birthDeath2);
      dataType death2 = std::get<1>(birthDeath2);
      // Compute relabel cost
      cost = computeDistance<dataType>(birth1, death1, birth2, death2, WassersteinPower);
      // Divide cost by two if not branch decomposition and not merged
      /*bool merged = isNodeMerged(tree1, nodeId1) or isNodeMerged(tree2, nodeId2);
      if(! BranchDecomposition and ! merged)
        cost /= 2;*/
      // Regularize
      if(NormalizedWasserstein and NormalizedWassersteinReg != 0)
        cost += NormalizedWassersteinReg * regularizationCost2<dataType>(tree1, nodeId1, tree2, nodeId2);
        
      return cost;
    }
    
    template <class dataType>
    dataType relabelCost(ftm::FTMTree_MT *tree1, ftm::idNode nodeId1, 
                         ftm::FTMTree_MT *tree2, ftm::idNode nodeId2){
      // Compute relabel cost
      dataType cost = relabelCostOnly<dataType>(tree1, nodeId1, tree2, nodeId2);
      // Compute deleteInsert cost
      dataType deleteInsertCost = deleteCost<dataType>(tree1, nodeId1) 
                                  + insertCost<dataType>(tree2, nodeId2);

      if(KeepSubtree and deleteInsertCost < cost)
        cost = deleteInsertCost;
      
      return cost;
    }
    
    // ----------------------------------------
    // Edit Distance Dynamic Programming Equations
    // ----------------------------------------
    template<class dataType>
    void computeEquation8(ftm::FTMTree_MT *tree1, ftm::idNode nodeI, int i,
                          std::vector<std::vector<dataType>> &treeTable,
                          std::vector<std::vector<dataType>> &forestTable){
      std::vector<ftm::idNode> childrens = getChildren(tree1, nodeI);
      forestTable[i][0] = 0;
      for(ftm::idNode children : childrens)
        forestTable[i][0] += treeTable[children+1][0]; 
    }

    template<class dataType>
    void computeEquation9(ftm::FTMTree_MT *tree1, ftm::idNode nodeI, int i,
                          std::vector<std::vector<dataType>> &treeTable,
                          std::vector<std::vector<dataType>> &forestTable){
      treeTable[i][0] = forestTable[i][0] + deleteCost<dataType>(tree1, nodeI);
    }

    template<class dataType>
    void computeEquation10(ftm::FTMTree_MT *tree2, ftm::idNode nodeJ, int j,
                          std::vector<std::vector<dataType>> &treeTable,
                          std::vector<std::vector<dataType>> &forestTable){
      std::vector<ftm::idNode> childrens = getChildren(tree2, nodeJ);
      forestTable[0][j] = 0;
      for(ftm::idNode children : childrens)
        forestTable[0][j] += treeTable[0][children+1];
    }

    template<class dataType>
    void computeEquation11(ftm::FTMTree_MT *tree2, ftm::idNode nodeJ, int j,
                          std::vector<std::vector<dataType>> &treeTable,
                          std::vector<std::vector<dataType>> &forestTable){
      treeTable[0][j] = forestTable[0][j] + insertCost<dataType>(tree2, nodeJ);
    }  


    // Compute first or second term of equation 12 or 13
    template<class dataType>
    std::tuple<dataType, ftm::idNode> computeTerm1_2(std::vector<ftm::idNode> &childrens, int ind,  
                                                      std::vector<std::vector<dataType>> &table,
                                                      bool computeTerm1){
      dataType tempMin = (childrens.size() == 0) ? 
                              ((computeTerm1) ? table[ind][0] : table[0][ind])
                              : std::numeric_limits<dataType>::max();
      ftm::idNode bestIdNode = 0;
      for(ftm::idNode children : childrens){
        children += 1;
        dataType temp;
        if(computeTerm1){
          temp = table[ind][children] - table[0][children];
        }else{
          temp = table[children][ind] - table[children][0];
        }
        if(temp < tempMin){
          tempMin = temp;
          bestIdNode = children;
        }
      }
      return std::make_tuple(tempMin, bestIdNode);
    }

    template<class dataType>
    void computeEquation12(ftm::FTMTree_MT *tree1, ftm::FTMTree_MT *tree2, int i, int j,
                          ftm::idNode nodeI, ftm::idNode nodeJ,
                          std::vector<std::vector<dataType>> &treeTable,
                          std::vector<std::vector<dataType>> &forestTable,
                          std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
                          std::vector<ftm::idNode> &children1, std::vector<ftm::idNode> &children2){
      dataType treeTerm1, treeTerm2, treeTerm3;
      std::tuple<dataType, ftm::idNode> treeCoTerm1, treeCoTerm2;
      // Term 1
      treeCoTerm1 = computeTerm1_2<dataType>(children2, i, treeTable, true);
      treeTerm1 = treeTable[0][j] + std::get<0>(treeCoTerm1);
      
      // Term 2
      treeCoTerm2 = computeTerm1_2<dataType>(children1, j, treeTable, false);
      treeTerm2 = treeTable[i][0] + std::get<0>(treeCoTerm2);
      
      // Term 3
      treeTerm3 = forestTable[i][j] + relabelCost<dataType>(tree1, nodeI, tree2, nodeJ);
      
      // Compute table value
      treeTable[i][j] = KeepSubtree ? std::min(std::min(treeTerm1, treeTerm2), treeTerm3) :
                                      treeTerm3;
      
      // Add backtracking information
      if(treeTable[i][j] == treeTerm3){
        treeBackTable[i][j] = std::make_tuple(i, j);
      }else if(treeTable[i][j] == treeTerm2){
        treeBackTable[i][j] = std::make_tuple(std::get<1>(treeCoTerm2), j);
      }else{
        treeBackTable[i][j] = std::make_tuple(i, std::get<1>(treeCoTerm1));
      }
    }

    // ----------------------------------------
    // Assignment Problem Functions
    // ----------------------------------------
    template<class dataType>
    void createCostMatrix(std::vector<std::vector<dataType>> &treeTable,
                          std::vector<ftm::idNode> &children1, std::vector<ftm::idNode> &children2,
                          std::vector<std::vector<dataType>> &costMatrix){
      int nRows = children1.size(), nCols = children2.size();
      for(int i = 0; i < nRows; ++i){
        int forestTableI = children1[i] + 1;
        for(int j = 0; j < nCols; ++j){
          int forestTableJ = children2[j] + 1;
          // Cost of assigning i and j
          costMatrix[i][j] = treeTable[forestTableI][forestTableJ];
        }
        // Cost of not assigning i
        costMatrix[i][nCols] = treeTable[forestTableI][0];
      }
      for(int j = 0; j < nCols; ++j){
          int forestTableJ = children2[j] + 1;
          // Cost of not assigning j
          costMatrix[nRows][j] = treeTable[0][forestTableJ];
      }
      //costMatrix[nRows][nCols] = std::numeric_limits<dataType>::max();
      costMatrix[nRows][nCols] = 0;
      //costMatrix[nRows][nCols] = maxValue*10; 
    }

    template<class dataType>
    dataType postprocessAssignment(std::vector<matchingTuple> &matchings, 
                                   std::vector<ftm::idNode> &children1, std::vector<ftm::idNode> &children2,
                                   std::vector<std::tuple<int, int>> &forestAssignment){
      dataType cost = 0;
      for(matchingTuple mTuple : matchings){
        cost += std::get<2>(mTuple);            
        if(std::get<0>(mTuple) >= children1.size() || std::get<1>(mTuple) >= children2.size())
          continue;
        int tableId1 = children1[std::get<0>(mTuple)] + 1;
        int tableId2 = children2[std::get<1>(mTuple)] + 1;
        forestAssignment.push_back(std::make_tuple(tableId1, tableId2));
      }
      return cost;
    }

    // ----------------------------------------
    // Utils
    // ----------------------------------------
    template<class dataType>
    void printTable(dataType *table, int nRows, int nCols){
      std::streamsize ss = std::cout.precision();
      std::cout << "      ";
      for(int j = 0; j < nCols; ++j)
        std::cout << j-1 << "    ";
      std::cout << std::endl;
      for(int i = 0; i < nRows; ++i){
        std::cout << std::setw(3) << std::setfill('0') << std::internal << i-1 << " | ";
        for(int j = 0; j < nCols; ++j){
          std::cout << std::fixed << std::setprecision(2) << table[i * nCols + j] << " ";
        }
        std::cout << std::endl << std::endl;
      }
      std::cout.precision(ss);
      std::cout << " ------------------------------ " << std::endl;
    }

    template<class dataType>
    void printTableVector(std::vector<std::vector<dataType>> &table){
      std::streamsize ss = std::cout.precision();
      std::cout << "      ";
      for(int j = 0; j < table[0].size(); ++j)
        std::cout << j-1 << "    ";
      std::cout << std::endl;
      for(int i = 0; i < table.size(); ++i){
        std::cout << std::setw(3) << std::setfill('0') << std::internal << i-1 << " | ";
        for(int j = 0; j < table[0].size(); ++j){
          std::cout << std::fixed << std::setprecision(2) << table[i][j] << " ";
        }
        std::cout << std::endl << std::endl;
      }
      std::cout.precision(ss);
      std::cout << " ------------------------------ " << std::endl;
    }
    
    /*void print3DMatrixTuple(std::vector<std::vector<std::vector<std::tuple<int, int>>>> mat){
      for(auto vec1 : mat){
        for(auto vec2 : vec1){
          for(auto vec3 : vec2){
          }
        }
      }
    }*/

    void printMatching(std::vector<matchingTuple> &matchings){
      std::cout << " ------------------------------ " << std::endl;
      for(matchingTuple mTuple : matchings)
        std::cout << std::get<0>(mTuple) << " - " << std::get<1>(mTuple) << " - " << std::get<2>(mTuple) << std::endl;
      std::cout << " ------------------------------ " << std::endl;
    }
    
    void printMatching(std::vector<std::tuple<ftm::idNode, ftm::idNode>> &matchings){
      std::cout << " ------------------------------ " << std::endl;
      for(auto mTuple : matchings)
        std::cout << std::get<0>(mTuple) << " - " << std::get<1>(mTuple) << std::endl;
      std::cout << " ------------------------------ " << std::endl;
    }
    
    template<class dataType>
    void printPairs(std::vector<std::tuple<SimplexId, SimplexId, dataType>> &treePairs){
      for(auto pair : treePairs)
        std::cout << std::get<0>(pair) << " _ " << std::get<1>(pair) << " _ " << std::get<2>(pair) << std::endl;
      std::cout << " ------------------------------ " << std::endl;
    }
    
    template<class dataType>
    void printOutputMatching(std::vector<std::tuple<ftm::idNode, ftm::idNode>> &outputMatching,
                             ftm::FTMTree_MT *tree1, ftm::FTMTree_MT *tree2, bool computeCosts=true){
      dataType cost = 0;
      bool tree1Done[tree1->getNumberOfNodes()]{false};
      bool tree2Done[tree2->getNumberOfNodes()]{false};
      for(std::tuple<ftm::idNode, ftm::idNode> matching : outputMatching){
        ftm::idNode node0 = std::get<0>(matching);
        ftm::idNode node1 = std::get<1>(matching);
        ftm::idNode node0Origin = tree1->getNode(node0)->getOrigin();
        ftm::idNode node1Origin = tree2->getNode(node1)->getOrigin();
        std::cout << node0 << " - " << node1 << " _ [ ";
        std::cout << "f(" << node0 << ")=" << tree1->getValue<dataType>(node0) << " _ ";
        std::cout << "g(" << node1 << ")=" << tree2->getValue<dataType>(node1)<< " ]" << " _ [ ";
        std::cout << "f(" << node0Origin << ")=" << tree1->getValue<dataType>(node0Origin) << " _ ";
        std::cout << "g(" << node1Origin  << ")=" << tree2->getValue<dataType>(node1Origin) << " ] ";
        
        if(computeCosts){
          dataType tempCost = relabelCost<dataType>(tree1, node0, tree2, node1);
          std::cout << "cost = " << tempCost << std::endl;
          cost += tempCost;
        }else
          std::cout << std::endl;
        tree1Done[node0] = true;
        tree2Done[node1] = true;
      }
      
      for(int i = 0 ; i < tree1->getNumberOfNodes(); ++i)
        if(not tree1Done[i] and not isNodeAlone(tree1, i)){
          dataType tempCost = deleteCost<dataType>(tree1, i);
          ftm::idNode nodeOrigin = tree1->getNode(i)->getOrigin();
          std::cout << "T1 " << i << " _ [ f(" << i << ") = " << tree1->getValue<dataType>(i);
          std::cout << "_ f(" << nodeOrigin << ") = " << tree1->getValue<dataType>(nodeOrigin);
          std::cout << "]";
          if(computeCosts){
            std::cout << " _ cost = " << tempCost << std::endl;
            cost += tempCost;
          }else 
            std::cout << std::endl;
        }
      for(int i = 0 ; i < tree2->getNumberOfNodes(); ++i)
        if(not tree2Done[i] and not isNodeAlone(tree2, i)){
          dataType tempCost = deleteCost<dataType>(tree2, i);
          ftm::idNode nodeOrigin = tree2->getNode(i)->getOrigin();
          std::cout << "T2 " << i << " _ [ g(" << i << ") = " << tree2->getValue<dataType>(i);
          std::cout << "_ g(" << nodeOrigin << ") = " << tree2->getValue<dataType>(nodeOrigin);
          std::cout << "] _ cost = " << tempCost << std::endl;
          if(computeCosts){
            std::cout << " _ cost = " << tempCost << std::endl;
            cost += tempCost;
          }else 
            std::cout << std::endl;
        }
      if(computeCosts)
        std::cout << "total cost = " << cost << " (" << std::sqrt(cost) << ")" << std::endl;
      
      std::cout << " ------------------------------ " << std::endl;
    }
    
    template<class dataType>
    void printTreeScalars(ftm::FTMTree_MT *tree){
      std::streamsize ss = std::cout.precision();
      for(int i = 0; i < tree->getNumberOfNodes(); ++i)
        std::cout << i << " _ " << std::setprecision(12) << tree->getValue<dataType>(i) << std::endl;
      std::cout << " ------------------------------ " << std::endl;
      std::cout.precision(ss);
    }
}; // FTMTreeEditDistanceBase class

#endif
