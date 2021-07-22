/// \ingroup base
/// \class ttk::FTMTreeUtils
/// \author XXXXX
///

#ifndef _FTMTREEUTILS_H
#define _FTMTREEUTILS_H

#include "FTMTreePPCompute.h"
#include <stack>
#include <ttkAlgorithm.h>

using namespace ttk;
using namespace ftm;

// ----------------------------------------
// Tree Structures
// ----------------------------------------

struct MergeTree {
  ftm::FTMTree_MT *tree;
  ftm::Scalars *scalars;
  ftm::Params *params;
};

// ----------------------------------------
// Tree Functions
// ----------------------------------------

// --------------------
// Is / Get / Set / Delete
// --------------------

// --- Is

bool isRoot(ftm::FTMTree_MT *tree, ftm::idNode nodeId);

bool isLeaf(ftm::FTMTree_MT *tree, ftm::idNode nodeId);

bool isNodeAlone(ftm::FTMTree_MT *tree, ftm::idNode nodeId);

bool isNodeMerged(ftm::FTMTree_MT *tree, ftm::idNode nodeId);

bool isNodeOriginDefined(ftm::FTMTree_MT *tree, ftm::idNode nodeId);

bool isNodeIdInconsistent(ftm::FTMTree_MT *tree, ftm::idNode nodeId);

bool isThereOnlyOnePersistencePair(ftm::FTMTree_MT *tree);

// Do not normalize node is if root or son of a merged root
bool notNeedToNormalize(ftm::FTMTree_MT *tree, ftm::idNode nodeId);

bool isFullMerge(ftm::FTMTree_MT *tree);

bool isMultiPersPair(ftm::FTMTree_MT *tree, ftm::idNode nodeId);

bool isBranchOrigin(ftm::FTMTree_MT *tree, ftm::idNode nodeId);

// --- Get

int getNumberOfNodeAlone(ftm::FTMTree_MT *tree);

int getRealNumberOfNodes(ftm::FTMTree_MT *tree);

ftm::idNode getParent(ftm::FTMTree_MT *tree, ftm::idNode nodeId);

ftm::idNode getRoot(ftm::FTMTree_MT *tree);

std::vector<ftm::idNode> getAllRoots(ftm::FTMTree_MT *tree);

int getNumberOfRoot(ftm::FTMTree_MT *tree);

// Get childrens (idNode) of a node in the tree
std::vector<ftm::idNode> getChildren(ftm::FTMTree_MT *tree, ftm::idNode nodeId);

int getNumberOfChildren(ftm::FTMTree_MT *tree, ftm::idNode nodeId);

std::vector<ttk::ftm::idNode> getLeaves(ttk::ftm::FTMTree_MT *tree);

int getNumberOfLeaves(ftm::FTMTree_MT *tree);

int getTreeDepth(ftm::FTMTree_MT *tree);

int getNodeLevel(ftm::FTMTree_MT *tree, ftm::idNode nodeId);

std::vector<int> getAllNodeLevel(ftm::FTMTree_MT *tree);

void getTreeBranching(FTMTree_MT *tree, std::vector<idNode> &branching, std::vector<int> &branchingID,
                      std::vector<std::vector<ftm::idNode>> &nodeBranching);

void getTreeBranching(FTMTree_MT *tree, std::vector<idNode> &branching, std::vector<int> &branchingID);

// Node must be a non-leaf node
std::tuple<std::vector<ftm::idNode>, std::vector<ftm::idNode>> 
getBranchOriginsFromThisBranch(FTMTree_MT *tree, ftm::idNode node);

std::vector<ftm::idNode> getBranchSubtree(std::vector<ftm::idNode> &branching, FTMTree_MT *tree, 
                                          ftm::idNode branchRoot);

double getValue(FTMTree_MT *tree, ftm::idNode node, int dataType);

double getPersistence(FTMTree_MT *tree, ftm::idNode node, int dataType);

// --- Set

void setParent(ftm::FTMTree_MT *tree, ftm::idNode nodeId, ftm::idNode newParentNodeId);

// --- Delete

void deleteIthArc(ftm::FTMTree_MT *tree, ftm::idNode nodeId, int arcIth);

// Delete arc of the node to its parent
void deleteParent(ftm::FTMTree_MT *tree, ftm::idNode nodeId);

// Delete node by keeping subtree 
void deleteNode(ftm::FTMTree_MT *tree, ftm::idNode nodeId);

// Delete node without keeping subtree 
void deleteSubtree(ftm::FTMTree_MT *tree, ftm::idNode nodeId);

// --------------------

template <class dataType>
bool isJoinTree(ftm::FTMTree_MT *tree){
  auto root = getRoot(tree);
  ftm::idNode child;
  if(isFullMerge(tree)){
    dataType min = std::numeric_limits<dataType>::max();
    for(int i = 0; i < tree->getNumberOfNodes(); ++i){
      dataType value = tree->getValue<dataType>(i);
      if(not isNodeAlone(tree, i) and value < min){
        min = value;
        child = i;
      }
    }
  }else
    child = getChildren(tree, root)[0];
  return tree->getValue<dataType>(root) > tree->getValue<dataType>(child);
}

// --------------------
// Utils
// --------------------

void myPause();

template <typename type>
static type myAbs(const type var) {
  return (var >= 0) ? var : -var;
}

std::stringstream printTree(ftm::FTMTree_MT *tree, bool doPrint=true);

void printTree(MergeTree *tree, bool doPrint=true);

void printTreeStats(ftm::FTMTree_MT * tree);

void printTreeStats(MergeTree * tree);

void printTreeStats(std::vector<ftm::FTMTree_MT *> &trees);

void printTreeStats(std::vector<MergeTree *> &trees);

template<class dataType>
void printTreeScalars(ftm::FTMTree_MT *tree, bool printNodeAlone=true){
  std::streamsize ss = std::cout.precision();
  for(int i = 0; i < tree->getNumberOfNodes(); ++i){
    ftm::idNode iOrigin = isNodeOriginDefined(tree, i) ? tree->getNode(i)->getOrigin() : i;
    if(printNodeAlone or 
       (not printNodeAlone and (not isNodeAlone(tree, i) or not isNodeAlone(tree, iOrigin)) ))
      std::cout << i << " _ " << std::setprecision(12) << tree->getValue<dataType>(i) << std::endl;
  }
  std::cout << " ------------------------------ " << std::endl;
  std::cout.precision(ss);
}

// --------------------
// Create/Delete/Modify Tree
// --------------------

void freeTree(ftm::FTMTree_MT *tree);

template <class dataType>
void freeMergeTree(MergeTree *mt, bool freeScalars=true){
  if(freeScalars)
    delete [] (dataType*)mt->scalars->values;
  delete mt->scalars;
  delete mt->params;
  delete mt->tree;
  delete mt;
}

// --------------------
// Persistence
// --------------------

template <class dataType>
std::tuple<dataType, dataType> getBirthDeath(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
  ftm::idNode originId = tree->getNode(nodeId)->getOrigin();
  if(isNodeOriginDefined(tree, nodeId)){ // Avoid error if origin is not defined
    dataType pers1 = tree->getValue<dataType>(nodeId);
    dataType pers2 = tree->getValue<dataType>(originId);
    dataType birth = std::min(pers1, pers2);
    dataType death = std::max(pers1, pers2); 
    return std::make_tuple(birth, death);
  }
  if(not isNodeAlone(tree, nodeId)){
    std::cout << "return bad birthDeath for " << nodeId << std::endl;
    myPause();
  }
  return std::make_tuple(0.0, 0.0);
}

template <class dataType>
dataType getBirth(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
  return std::get<0>(getBirthDeath<dataType>(tree, nodeId));
}

template<class dataType>
void printNode(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
  auto birthDeath = getBirthDeath<dataType>(tree, nodeId);
  std::cout << "nodeId = " << nodeId << " (" << std::get<0>(birthDeath) << ") _ originId = " << tree->getNode(nodeId)->getOrigin() << " (" << std::get<1>(birthDeath) << ")" << std::endl;
}

template <class dataType>
bool isEqual(dataType first, dataType two, double eps=1e-6){
  return myAbs<dataType>(first - two) < eps*std::max(myAbs<dataType>(first), myAbs<dataType>(two));
}

template <class dataType>
bool isInferior(dataType first, dataType two){
  return first < two and not isEqual<dataType>(first, two);
}

template <class dataType>
bool isSuperior(dataType first, dataType two){
  return first > two and not isEqual<dataType>(first, two);
}

template <class dataType>
bool isParentInconsistent(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
  auto parentBirthDeath = getBirthDeath<dataType>(tree, getParent(tree, nodeId));
  dataType parentBirth = std::get<0>(parentBirthDeath);
  dataType parentDeath = std::get<1>(parentBirthDeath);
  auto birthDeath = getBirthDeath<dataType>(tree, nodeId);
  dataType birth = std::get<0>(birthDeath);
  dataType death = std::get<1>(birthDeath);
  bool parentInconsistent = isInferior<dataType>(parentDeath, death) or 
                            isSuperior<dataType>(parentBirth, birth);
  return parentInconsistent;
}

template <class dataType>
bool verifyBranchDecompositionInconsistency(ftm::FTMTree_MT *treeNew, bool doPause=false){
  bool inconsistency = false;
  std::queue<ftm::idNode> queue;
  queue.push(getRoot(treeNew));
  while(!queue.empty()){
    ftm::idNode node = queue.front();
    queue.pop();
    if(!isRoot(treeNew, node) and isParentInconsistent<dataType>(treeNew, node)){
      std::cout << "inconsistency" << std::endl;
      printNode<dataType>(treeNew, node);
      printNode<dataType>(treeNew, getParent(treeNew, node));
      if(doPause) myPause();
      inconsistency = true;
    }
    for(ftm::idNode child : getChildren(treeNew, node))
      queue.push(child);
  }
  return inconsistency;
}

/*template <class dataType>
dataType getMinMaxGlobal(ftm::FTMTree_MT *tree, bool getMin=true){
  dataType minMax = tree->getValue<dataType>(0);
  for(int i = 0 ; i < tree->getNumberOfNodes() ; ++i)
    if( (getMin and tree->getValue<dataType>(i) < minMax) or
        (not getMin and tree->getValue<dataType>(i) > minMax) )
      minMax = tree->getValue<dataType>(i);
  
  return minMax;
}*/

template <class dataType>
dataType getMinMaxLocal(ftm::FTMTree_MT *tree, ftm::idNode nodeId, bool getMin=true){
  auto nodeIdParent = getParent(tree, nodeId);
  
  if(notNeedToNormalize(tree, nodeId))
    return getMin ? 0.0 : 1.0;
  
  // General case
  std::tuple<dataType, dataType> birthDeath = getBirthDeath<dataType>(tree, nodeIdParent);
  
  // Verify inconsistency
  if(isParentInconsistent<dataType>(tree, nodeId)){
    std::cout << "inconsistency" << std::endl;
    //if(tree->getNumberOfNodes() < 1000) 
      printTree(tree);
      //printPairsFromTree<dataType>(tree, true);
    printNode<dataType>(tree, nodeId);
    printNode<dataType>(tree, nodeIdParent);
    myPause();
  }

  return getMin ? std::get<0>(birthDeath) : std::get<1>(birthDeath); 
}

template <class dataType>
std::tuple<dataType, dataType> getMinMaxLocalT(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
  dataType min = getMinMaxLocal<dataType>(tree, nodeId);
  dataType max = getMinMaxLocal<dataType>(tree, nodeId, false);
  return std::make_tuple(min, max);
}

template <class dataType>
std::tuple<double, double> getNormalizedBirthDeathDouble(ftm::FTMTree_MT *tree, ftm::idNode nodeId, 
                                                       dataType newMin=0.0, dataType newMax=1.0){
  auto birthDeath = getBirthDeath<dataType>(tree, nodeId);
  double birth = std::get<0>(birthDeath);
  double death = std::get<1>(birthDeath);
  dataType shiftMin = getMinMaxLocal<dataType>(tree, nodeId);
  dataType shiftMax = getMinMaxLocal<dataType>(tree, nodeId, false);
  birth = (newMax - newMin) * (birth - shiftMin) / (shiftMax - shiftMin);// + newMin;
  death = (newMax - newMin) * (death - shiftMin) / (shiftMax - shiftMin);// + newMin;
  return std::make_tuple(birth, death);
}

template <class dataType>
std::tuple<dataType, dataType> getNormalizedBirthDeath(ftm::FTMTree_MT *tree, ftm::idNode nodeId, 
                                                       dataType newMin=0.0, dataType newMax=1.0){
  auto birthDeath = getBirthDeath<dataType>(tree, nodeId);
  dataType birth = std::get<0>(birthDeath);
  dataType death = std::get<1>(birthDeath);
  dataType shiftMin = getMinMaxLocal<dataType>(tree, nodeId);
  dataType shiftMax = getMinMaxLocal<dataType>(tree, nodeId, false);
  birth = (newMax - newMin) * (birth - shiftMin) / (shiftMax - shiftMin);// + newMin;
  death = (newMax - newMin) * (death - shiftMin) / (shiftMax - shiftMin);// + newMin;
  return std::make_tuple(birth, death);
}

template <class dataType>
std::tuple<dataType, dataType> getNewMinMax(ftm::FTMTree_MT *tree1, ftm::idNode nodeId1, 
                                            ftm::FTMTree_MT *tree2, ftm::idNode nodeId2){
  dataType shiftMin1 = getMinMaxLocal<dataType>(tree1, nodeId1);
  dataType shiftMax1 = getMinMaxLocal<dataType>(tree1, nodeId1, false);
  dataType shiftMin2 = getMinMaxLocal<dataType>(tree2, nodeId2);
  dataType shiftMax2 = getMinMaxLocal<dataType>(tree2, nodeId2, false);
  return std::make_tuple((shiftMin1+shiftMin2)/2.0, (shiftMax1+shiftMax2)/2.0);
}

template <class dataType>
std::tuple<dataType, dataType> getRescaledBirthDeath(ftm::FTMTree_MT *tree1, ftm::idNode nodeId1,
                                                     ftm::FTMTree_MT *tree2, ftm::idNode nodeId2){
  auto newMinMax = getNewMinMax<dataType>(tree1, nodeId1, tree2, nodeId2);
  return getNormalizedBirthDeath<dataType>(tree1, nodeId1, std::get<0>(newMinMax), std::get<1>(newMinMax));
}

template <class dataType>
dataType getMinMaxLocalFromVector(ftm::FTMTree_MT *tree, ftm::idNode nodeId, 
                                  std::vector<dataType> &scalarsVector, bool getMin=true){
  auto nodeIdParent = getParent(tree, nodeId);
  
  if(notNeedToNormalize(tree, nodeId))
    return getMin ? 0.0 : 1.0;
  
  // General case
  dataType death = scalarsVector[nodeIdParent];
  dataType birth = scalarsVector[tree->getNode(nodeIdParent)->getOrigin()];
  if(death < birth){
    dataType temp = death;
    death = birth;
    birth = temp;
  }

  return getMin ? birth : death; 
}

template <class dataType>
std::tuple<dataType, dataType> getMinMaxLocalFromVectorT(ftm::FTMTree_MT *tree, ftm::idNode nodeId, 
                                                         std::vector<dataType> &scalarsVector){
  dataType min = getMinMaxLocalFromVector<dataType>(tree, nodeId, scalarsVector);
  dataType max = getMinMaxLocalFromVector<dataType>(tree, nodeId, scalarsVector, false);
  return std::make_tuple(min, max);
}

template <class dataType>
std::tuple<dataType, dataType> getNewMinMaxFromVector(ftm::FTMTree_MT *tree1, ftm::idNode nodeId1, 
                                                      ftm::FTMTree_MT *tree2, ftm::idNode nodeId2,
                                                      std::vector<dataType> &scalarsVector){
  dataType shiftMin1 = getMinMaxLocal<dataType>(tree1, nodeId1);
  dataType shiftMax1 = getMinMaxLocal<dataType>(tree1, nodeId1, false);
  dataType shiftMin2 = getMinMaxLocalFromVector<dataType>(tree2, nodeId2, scalarsVector);
  dataType shiftMax2 = getMinMaxLocalFromVector<dataType>(tree2, nodeId2, scalarsVector, false);
  return std::make_tuple((shiftMin1+shiftMin2)/2.0, (shiftMax1+shiftMax2)/2.0);
}

template <class dataType>
std::tuple<dataType, dataType> getRescaledBirthDeathFromVector(ftm::FTMTree_MT *tree1, ftm::idNode nodeId1,
                                                               ftm::FTMTree_MT *tree2, ftm::idNode nodeId2,
                                                               std::vector<dataType> &scalarsVector){
  auto newMinMax = getNewMinMaxFromVector<dataType>(tree1, nodeId1, tree2, nodeId2, scalarsVector);
  return getNormalizedBirthDeath<dataType>(tree1, nodeId1, std::get<0>(newMinMax), std::get<1>(newMinMax));
}

template<class dataType>
dataType getNodePersistence(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
  std::tuple<dataType, dataType> birthDeath = getBirthDeath<dataType>(tree, nodeId);
  return std::get<1>(birthDeath) - std::get<0>(birthDeath);
}

template<class dataType>
dataType getMaximumPersistence(ftm::FTMTree_MT *tree){
  ftm::idNode root = getRoot(tree);
  bool fullMerge = (tree->getNode(root)->getOrigin() == root);

  // Classic case
  if(not fullMerge)
    return getNodePersistence<dataType>(tree, getRoot(tree));

  // Full merge case
  dataType maxPers = std::numeric_limits<dataType>::lowest();
  for(int i = 0; i < tree->getNumberOfNodes(); ++i)
    if(not isNodeAlone(tree, i) and isNodeOriginDefined(tree, i) and
       tree->getNode(i)->getOrigin() == root)
      maxPers = std::max(maxPers, getNodePersistence<dataType>(tree, i));

  return maxPers;
}

template<class dataType>
dataType getSecondMaximumPersistence(ftm::FTMTree_MT *tree){
  ftm::idNode root = getRoot(tree);
  dataType pers = std::numeric_limits<dataType>::lowest();
  for(int i = 0; i < tree->getNumberOfNodes(); ++i){
    if(isLeaf(tree, i) and not isNodeAlone(tree, i) and isNodeOriginDefined(tree, i)){
      ftm::idNode nodeOrigin = tree->getNode(i)->getOrigin();
      if(not (nodeOrigin == root and tree->getNode(nodeOrigin)->getOrigin() == i))
        pers = std::max(pers, getNodePersistence<dataType>(tree, i));
    }
  }

  return pers;
}

template<class dataType>
void getPersistencePairs(ftm::FTMTree_MT *tree, 
                        //std::vector<std::tuple<SimplexId, SimplexId>> &pairs){
                        std::vector<std::tuple<SimplexId, SimplexId, dataType>> &pairs){
  ttk::ftm::FTMTreePPCompute pairsCompute;
  pairsCompute.computePersistencePairs<dataType>(tree, pairs);
}

template<class dataType>
std::vector<std::tuple<SimplexId, SimplexId, dataType>> computePersistencePairs(ftm::FTMTree_MT *tree){
  //std::vector<std::tuple<SimplexId, SimplexId>> pairs;
  std::vector<std::tuple<SimplexId, SimplexId, dataType>> pairs;
  getPersistencePairs<dataType>(tree, pairs);
  for(auto pair: pairs){
    // TODO verify this, we should remove saddle-saddle pairs
    if(tree->getNode(std::get<0>(pair))->getOrigin() < std::get<0>(pair) and
       tree->getNode(std::get<0>(pair))->getOrigin() >= 0)
      tree->getNode(tree->getNode(std::get<0>(pair))->getOrigin())->setOrigin(std::get<1>(pair));

    tree->getNode(std::get<0>(pair))->setOrigin(std::get<1>(pair));
    tree->getNode(std::get<1>(pair))->setOrigin(std::get<0>(pair));
  }
  return pairs;
}

template<class dataType>
void getPersistencePairsFromTree(ftm::FTMTree_MT *tree, 
                                 std::vector<std::tuple<ftm::idNode, ftm::idNode, dataType>> &pairs,
                                 bool useBD){
  std::vector<ftm::idNode> nodes;
  if(useBD){
    for(int i = 0; i < tree->getNumberOfNodes(); ++i)
      if(!isNodeAlone(tree, i) and tree->getNode(i)->getOrigin() != i)
        nodes.push_back(i);
  }else
    nodes = getLeaves(tree);
  for(auto node : nodes){
    auto pers = getNodePersistence<dataType>(tree, node);
    pairs.push_back(std::make_tuple(node, tree->getNode(node)->getOrigin(), pers));
  }
  auto comp = [&](const std::tuple<ftm::idNode, ftm::idNode, dataType> a,
              const std::tuple<ftm::idNode, ftm::idNode, dataType> b) {
    return std::get<2>(a) < std::get<2>(b);
  };
  sort(pairs.begin(), pairs.end(), comp);
}

template<class dataType>
std::stringstream printPairsFromTree(ftm::FTMTree_MT *tree, bool useBD=false, bool printPairs=true, bool doPrint=true){
  std::stringstream ss;
  std::vector<std::tuple<ftm::idNode, ftm::idNode, dataType>> pairs;
  getPersistencePairsFromTree(tree, pairs, useBD);
  ss << "size=" << pairs.size() << std::endl;
  if(printPairs)
    for(auto pair : pairs){
      ss << std::get<0>(pair) << " (" << tree->getValue<dataType>(std::get<0>(pair)) << ") _ ";
      ss << std::get<1>(pair) << " (" << tree->getValue<dataType>(std::get<1>(pair)) << ") _ ";
      ss << std::get<2>(pair) << std::endl;
    }
  ss << " ------------------------------ " << std::endl;
  
  if(doPrint)
    std::cout << ss.str();
  return ss;
}

template <class dataType>
std::vector<ftm::idNode> getMultiPersOrigins(ftm::FTMTree_MT *tree, bool useBD){
  std::vector<ftm::idNode> multiPersOrigins;
  
  std::vector<std::tuple<ftm::idNode, ftm::idNode, dataType>> pairs;
  getPersistencePairsFromTree(tree, pairs, useBD);
  //std::vector<ftm::idNode> origins(tree->getNumberOfNodes(), -1);
  std::vector<std::vector<ftm::idNode>> origins(tree->getNumberOfNodes());
  std::vector<bool> birthFound(tree->getNumberOfNodes(), false);
  for(auto pair : pairs){
    ftm::idNode nodeBirth = std::get<0>(pair);
    ftm::idNode nodeDeath = std::get<1>(pair);
    
    origins[nodeDeath].push_back(nodeBirth);
    birthFound[nodeBirth] = true;
  }
  
  for(int i = 0; i < origins.size(); ++i)
    if(birthFound[i])
      for(auto node : origins[i])
        multiPersOrigins.push_back(node);
  
  return multiPersOrigins;
}

std::vector<std::vector<ftm::idNode>> getMultiPersOriginsVectorFromTree(ftm::FTMTree_MT *tree);

std::stringstream printMultiPersOriginsVectorFromTree(ftm::FTMTree_MT *tree, bool doPrint=true);

template<class dataType>
void printMultiPersPairsFromTree(ftm::FTMTree_MT *tree, bool useBD=false, bool printPairs=true){
  std::vector<std::tuple<ftm::idNode, ftm::idNode, dataType>> pairs;
  getPersistencePairsFromTree(tree, pairs, useBD);
  std::vector<int> noOrigin(tree->getNumberOfNodes(), 0);
  int noMultiPers = 0;
  for(auto pair : pairs){
    noOrigin[std::get<0>(pair)]++;
    noMultiPers += (noOrigin[std::get<0>(pair)] > 1 ) ? 1 : 0;
    noOrigin[std::get<1>(pair)]++;
    noMultiPers += (noOrigin[std::get<1>(pair)] > 1 ) ? 1 : 0;
  }
  std::cout << "Number of multi pers pairs : " << noMultiPers << std::endl;
  if(printPairs){
    auto multiPers = getMultiPersOrigins<dataType>(tree, useBD);
    for(auto node : multiPers)
      std::cout << node << std::endl;
  }
  std::cout << " ------------------------------ " << std::endl;
}

template<class dataType>
bool isDeathHigher(ftm::FTMTree_MT *tree1, ftm::idNode tree1Node, 
                   ftm::FTMTree_MT *tree2, ftm::idNode tree2Node){
  dataType tree1NodeDeath = std::get<1>(getBirthDeath<dataType>(tree1, tree1Node));
  dataType tree2NodeDeath = std::get<1>(getBirthDeath<dataType>(tree2, tree2Node));
  return tree1NodeDeath > tree2NodeDeath;
}

template<class dataType>
int getNoPairs(ftm::FTMTree_MT *tree, bool useBD=false){
  std::vector<std::tuple<ftm::idNode, ftm::idNode, dataType>> pairs;
  getPersistencePairsFromTree(tree, pairs, useBD);
  return pairs.size();
}

template<class dataType>
bool isImportantPair(ftm::FTMTree_MT *tree, ftm::idNode nodeId, double threshold){
  dataType rootPers = getNodePersistence<dataType>(tree, getRoot(tree));
  if(threshold > 1)
    threshold /= 100.0;
  threshold = rootPers * threshold;
  return getNodePersistence<dataType>(tree, nodeId) > threshold;
}

// --------------------
// Create/Delete/Modify Tree
// --------------------

template <class dataType>
MergeTree* createEmptyMergeTree(int scalarSize){
  ftm::Scalars *scalars = new ftm::Scalars();
  scalars->size = scalarSize;
  dataType* scalarsValues = nullptr;
  scalars->values = (void*)scalarsValues;
  
  // Init Tree
  ftm::Params *params = new ftm::Params();
  //FTMTree_MT treeTemp(params, nullptr, scalars, Join_Split);
  ftm::FTMTree_MT *treeNew = new ftm::FTMTree_MT(params, nullptr, scalars, ftm::Join_Split);
  treeNew->makeAlloc();
  
  MergeTree *mergeTree = new MergeTree();
  mergeTree->tree = treeNew;
  mergeTree->scalars = scalars;
  mergeTree->params = params;
  
  return mergeTree;
}

// TODO Will be deprecated
template <class dataType>
ftm::FTMTree_MT* createEmptyTree(int scalarSize){
  //std::cout << "ftm::FTMTree_MT* createEmptyTree(int) will be deprecated" << std::endl;
  MergeTree *mergeTree = createEmptyMergeTree<dataType>(scalarSize);
  return mergeTree->tree;
}

// TODO Will be deprecated
template <class dataType>
void setTreeScalars(ftm::FTMTree_MT *tree, std::vector<dataType> &scalarsVector){
  //std::cout << "setTreeScalars(ftm::FTMTree_MT *, std::vector<dataType>) will be deprecated" << std::endl;
  dataType* scalarsValues = new dataType[scalarsVector.size()];
  for(int i = 0; i < scalarsVector.size(); ++i)
    scalarsValues[i] = scalarsVector[i];
  tree->setScalars((void*)scalarsValues);
}

template <class dataType>
void setTreeScalars(ftm::Scalars *scalars, std::vector<dataType> &scalarsVector){
  delete[] static_cast<dataType*>(scalars->values);
  scalars->values = nullptr;
  dataType* scalarsValues = new dataType[scalarsVector.size()];
  for(int i = 0; i < scalarsVector.size(); ++i)
    scalarsValues[i] = scalarsVector[i];
  scalars->values = (void*)(scalarsValues);
  scalars->size = scalarsVector.size();
}

template <class dataType>
void setTreeScalars(MergeTree *mergeTree, std::vector<dataType> &scalarsVector){
  setTreeScalars<dataType>(mergeTree->scalars, scalarsVector);
}

template <class dataType>
std::vector<dataType> getTreeScalars(ftm::FTMTree_MT *tree){
  std::vector<dataType> scalarsVector;
  for(int i = 0; i < tree->getNumberOfNodes(); ++i)
    scalarsVector.push_back(tree->getValue<dataType>(i));
  return scalarsVector;
}

template <class dataType>
std::vector<dataType> getTreeScalars(MergeTree *mergeTree){
  return getTreeScalars<dataType>(mergeTree->tree);
}

template <class dataType>
void appendTreeScalars(MergeTree *mergeTree, std::vector<dataType> &newScalarsVector){
  std::vector<dataType> scalarsVector = getTreeScalars<dataType>(mergeTree);
  scalarsVector.insert(scalarsVector.end(), newScalarsVector.begin(), newScalarsVector.end());
  setTreeScalars(mergeTree, scalarsVector);
}

void copyMergeTreeStructure(MergeTree *mergeTree, ftm::FTMTree_MT *tree);

template <class dataType>
MergeTree* copyMergeTree(ftm::FTMTree_MT *tree, bool doSplitMultiPersPairs=false){
  std::vector<dataType> scalarsVector = getTreeScalars<dataType>(tree);
  
  // Get multi persistence pairs
  std::vector<ftm::idNode> multiPersOrigins;
  if(doSplitMultiPersPairs){
    multiPersOrigins = getMultiPersOrigins<dataType>(tree, true);
    for(ftm::idNode nodeOrigin : multiPersOrigins){
      scalarsVector[nodeOrigin] = tree->getValue<dataType>(tree->getNode(nodeOrigin)->getOrigin());
      scalarsVector.push_back(tree->getValue<dataType>(nodeOrigin));
    }
  }
  
  // Create new tree  
  MergeTree *mergeTree = createEmptyMergeTree<dataType>(scalarsVector.size());
  ftm::FTMTree_MT *treeNew = mergeTree->tree;
  setTreeScalars<dataType>(mergeTree, scalarsVector);
  
  // Copy tree structure
  copyMergeTreeStructure(mergeTree, tree);
  
  // Add multi persistence nodes origins
  if(doSplitMultiPersPairs){
    for(ftm::idNode nodeOrigin : multiPersOrigins){
      int nodeCpt = treeNew->getNumberOfNodes();
      treeNew->makeNode(nodeCpt);
      treeNew->getNode(nodeCpt)->setOrigin(nodeOrigin);
      treeNew->getNode(nodeOrigin)->setOrigin(nodeCpt);
    }
  }

  return mergeTree;
}

template <class dataType>
MergeTree* copyMergeTree(MergeTree *mergeTree, bool doSplitMultiPersPairs=false){
  return copyMergeTree<dataType>(mergeTree->tree, doSplitMultiPersPairs);
}

// TODO  Will be deprecated
template <class dataType>
ftm::FTMTree_MT* copyTree(ftm::FTMTree_MT *tree){
  //std::cout << "ftm::FTMTree_MT* copyTree(ftm::FTMTree_MT *) will be deprecated" << std::endl;
  MergeTree *mergeTree = copyMergeTree<dataType>(tree);
  return mergeTree->tree;
}

template <class dataType>
ftm::FTMTree_MT* fusionTrees(std::vector<ftm::FTMTree_MT *> trees, bool copyScalars=true){  
  int scalarsSize = 0;
  for(int i = 0; i < trees.size(); ++i)
    scalarsSize += trees[i]->getNumberOfNodes();
  ftm::FTMTree_MT *treeNew = createEmptyTree<dataType>(scalarsSize);
  if(copyScalars){
    std::vector<dataType> scalarsVector;
    int actualTree = 0;
    int offset = 0;
    for(int i = 0; i < scalarsSize; ++i){
      scalarsVector.push_back(trees[actualTree]->getValue<dataType>(i-offset));
      if(i >= trees[actualTree]->getNumberOfNodes() + offset){
        offset += trees[actualTree]->getNumberOfNodes();
        actualTree++;
      }
    }
    setTreeScalars<dataType>(treeNew, scalarsVector);
  }
  
  // Add Nodes
  for(int i = 0; i < scalarsSize; ++i)
    treeNew->makeNode(i);
  int actualTree = 0;
  int offset = 0;
  for(int i = 0; i < scalarsSize; ++i){
    auto newOri = trees[actualTree]->getNode(i-offset)->getOrigin()+offset;
    treeNew->getNode(i)->setOrigin(newOri);
    if(i >= trees[actualTree]->getNumberOfNodes() + offset){
      offset += trees[actualTree]->getNumberOfNodes();
      actualTree++;
    }
  }
  
  return treeNew;
}

// Remove unused nodes
template <class dataType>
MergeTree* cleanMergeTree(ftm::FTMTree_MT *tree, std::vector<int> &nodeCorr, bool useBD=true){  
  // Create new tree
  int newNoNodes = getRealNumberOfNodes(tree)*2;
  MergeTree *mTreeNew = createEmptyMergeTree<dataType>(newNoNodes);
  ftm::FTMTree_MT *treeNew = mTreeNew->tree;
  std::vector<dataType> newScalarsVector(newNoNodes, 0);
  
  // Copy the old tree structure
  std::vector<int> nodeDone(tree->getNumberOfNodes(), 0);
  nodeCorr = std::vector<int>(tree->getNumberOfNodes(), -1);
  std::vector<std::vector<ftm::idNode>> treeMultiPers;
  if(not useBD)
    treeMultiPers = getMultiPersOriginsVectorFromTree(tree);
  std::queue<ftm::idNode> queue;
  for(auto leaf : getLeaves(tree))
    queue.push(leaf);
  while(!queue.empty()){
    ftm::idNode node = queue.front();
    queue.pop();
    ftm::idNode nodeOrigin = tree->getNode(node)->getOrigin();
    
    if(useBD){
      int nodeOriginIndex = treeNew->getNumberOfNodes();
      if(nodeCorr[nodeOrigin] == -1)
        treeNew->makeNode(nodeOriginIndex);
      else
        nodeOriginIndex = nodeCorr[nodeOrigin];
      int nodeIndex = treeNew->getNumberOfNodes();
      if(nodeCorr[node] == -1)
        treeNew->makeNode(nodeIndex);
      else
        nodeIndex = nodeCorr[node];
      
      if(nodeCorr[nodeOrigin] == -1)
        treeNew->getNode(nodeOriginIndex)->setOrigin(nodeIndex);
      treeNew->getNode(nodeIndex)->setOrigin(nodeOriginIndex);
      
      newScalarsVector[nodeOriginIndex] = tree->getValue<dataType>(nodeOrigin);
      newScalarsVector[nodeIndex] = tree->getValue<dataType>(node);
      nodeCorr[nodeOrigin] = nodeOriginIndex;
      nodeCorr[node] = nodeIndex;
    }else{
      int nodeCpt = treeNew->getNumberOfNodes();
      treeNew->makeNode(nodeCpt);
      if(!isLeaf(tree, node)){
        treeNew->getNode(nodeCpt)->setOrigin(nodeCorr[nodeOrigin]);
        if(not (isRoot(tree, node) and node == nodeOrigin))
          treeNew->getNode(nodeCorr[nodeOrigin])->setOrigin(nodeCpt);
        for(auto nodeMultiPers : treeMultiPers[node])
          treeNew->getNode(nodeCorr[nodeMultiPers])->setOrigin(nodeCpt);
      }else if(isNodeAlone(tree, nodeOrigin)){ // saddle merged
        treeNew->makeNode(nodeCpt+1);
        newScalarsVector[nodeCpt+1] = tree->getValue<dataType>(nodeOrigin);
        nodeCorr[nodeOrigin] = nodeCpt+1;
        treeNew->getNode(nodeCpt)->setOrigin(nodeCorr[nodeOrigin]);
        treeNew->getNode(nodeCorr[nodeOrigin])->setOrigin(nodeCpt);
      }
      newScalarsVector[nodeCpt] = tree->getValue<dataType>(node);
      nodeCorr[node] = nodeCpt;
    }
    
    for(auto child : getChildren(tree, node))
      treeNew->makeSuperArc(nodeCorr[child], nodeCorr[node]);
    
    if(!isRoot(tree, node)){
      ftm::idNode parent = getParent(tree, node);
      nodeDone[parent] += 1;
      if(nodeDone[parent] == getNumberOfChildren(tree, parent))
        queue.push(parent);
    }
  }
  
  // Set new scalars
  setTreeScalars<dataType>(mTreeNew, newScalarsVector);
  
  // Manage full merge
  auto treeRoot = getRoot(tree);
  if(tree->getNode(treeRoot)->getOrigin() == treeRoot){
    auto treeNewRoot = getRoot(treeNew);
    //auto treeNewRootOrigin = tree->getNode(treeNewRoot)->getOrigin();
    //treeNew->getNode(treeNewRootOrigin)->setOrigin(treeNewRootOrigin);
    treeNew->getNode(treeNewRoot)->setOrigin(treeNewRoot);
  } 
  
  // Return new tree
  return mTreeNew;
}

// TODO will be deprecated
template <class dataType> 
ftm::FTMTree_MT* cleanMergeTree2(ftm::FTMTree_MT *tree, std::vector<int> &nodeCorr, bool useBD=true){
  auto mT = cleanMergeTree<dataType>(tree, nodeCorr, useBD);
  return mT->tree;
}

// TODO will be deprecated
template <class dataType> 
ftm::FTMTree_MT* cleanMergeTree(ftm::FTMTree_MT *tree, bool useBD=true){
  std::vector<int> nodeCorr;
  auto mT = cleanMergeTree<dataType>(tree, nodeCorr, useBD);
  return mT->tree;
}

template <class dataType>
void cleanMergeTree(MergeTree *&mTree, std::vector<int> &nodeCorr, bool useBD=true){
  auto mTreeNew = cleanMergeTree<dataType>(mTree->tree, nodeCorr, useBD);
  
  // Free old tree and replace it with new tree
  freeMergeTree<dataType>(mTree, useBD);
  mTree = mTreeNew;
}

template <class dataType>
void cleanMergeTree(MergeTree *&mTree, bool useBD=true){
  std::vector<int> nodeCorr;
  cleanMergeTree<dataType>(mTree, nodeCorr, useBD);
}

template <class dataType>
void updateNodesAndScalars(MergeTree *&mTree1, std::vector<dataType> &newScalarsVector,
                           std::vector<std::tuple<ftm::idNode, bool>> &origins){
  ftm::FTMTree_MT *tree1 = mTree1->tree;

  // Create new tree
  MergeTree *mTreeNew = createEmptyMergeTree<dataType>(newScalarsVector.size());
  setTreeScalars<dataType>(mTreeNew, newScalarsVector);

  // Copy the old tree structure
  copyMergeTreeStructure(mTreeNew, tree1);

  //
  int noNodes = mTreeNew->tree->getNumberOfNodes();
  for(int i = noNodes; i < newScalarsVector.size(); ++i){
    mTreeNew->tree->makeNode(i);
    auto tup = origins[i-noNodes];
    ftm::idNode origin = std::get<0>(tup);
    bool reverse = std::get<1>(tup);
    mTreeNew->tree->getNode(i)->setOrigin(origin);
    if(reverse)
      mTreeNew->tree->getNode(origin)->setOrigin(i);
  }

  // Free old tree and replace it with new tree
  freeMergeTree<dataType>(mTree1);
  mTree1 = mTreeNew;
}

std::vector<ftm::FTMTree_MT*> mergeTreeToFTMTree(std::vector<MergeTree*> &trees);

// --------------------
// Barycenter Utils
// --------------------

template <class dataType>
ftm::idNode getNodesAndScalarsToAdd(MergeTree *mTree1, ftm::idNode nodeId1, ftm::FTMTree_MT *tree2, 
                                    ftm::idNode nodeId2, std::vector<dataType> &newScalarsVector, 
                                    std::vector<std::tuple<ftm::idNode, ftm::idNode, int>> &nodesToProcess,
                                    ftm::idNode nodeCpt, int i){
  ftm::FTMTree_MT *tree1 = mTree1->tree;
  // Get nodes and scalars to add
  std::queue<std::tuple<ftm::idNode, ftm::idNode>> queue;
  queue.push(std::make_tuple(nodeId2, nodeId1));
  nodesToProcess.push_back(std::make_tuple(nodeId2, nodeId1, i));
  while(!queue.empty()){
    auto queueTuple = queue.front();
    queue.pop();
    ftm::idNode node = std::get<0>(queueTuple);
    // Get scalars
    newScalarsVector.push_back(tree2->getValue<dataType>(tree2->getNode(node)->getOrigin()));
    newScalarsVector.push_back(tree2->getValue<dataType>(node));
    // Process children
    auto children = getChildren(tree2, node);
    for(auto child : children){
      queue.push(std::make_tuple(child, nodeCpt+1));
      nodesToProcess.push_back(std::make_tuple(child, nodeCpt+1, i));
    }
    nodeCpt += 2; // we will add two nodes (birth and death)
  }
  
  return nodeCpt;
}

template <class dataType>
std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>> addNodes(MergeTree *mTree1, int noTrees, 
                                    std::vector<std::tuple<ftm::idNode, ftm::idNode, int>> &nodesToProcess){
  ftm::FTMTree_MT *tree1 = mTree1->tree;
  
  // Add nodes
  std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>> nodesProcessed(noTrees);
  for(auto processTuple : nodesToProcess){
    ftm::idNode parent = std::get<1>(processTuple);
    ftm::idNode nodeTree1 = tree1->getNumberOfNodes();
    int index = std::get<2>(processTuple);
    nodesProcessed[index].push_back(std::make_tuple(nodeTree1+1, std::get<0>(processTuple)));
    // Make node and its origin
//std::cout << nodeTree1+1 << " _ " << parent << std::endl;
    tree1->makeNode(nodeTree1);
    tree1->makeNode(nodeTree1+1);    
    setParent(tree1, nodeTree1+1, parent);
    tree1->getNode(nodeTree1)->setOrigin(nodeTree1+1);
    tree1->getNode(nodeTree1+1)->setOrigin(nodeTree1);
    //tree1->printTree2();
  }
  
  return nodesProcessed;
}

template <class dataType>
std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>> updateNodesAndScalars(MergeTree *&mTree1, 
                                                                                     int noTrees, 
                                    std::vector<std::tuple<ftm::idNode, ftm::idNode, int>> &nodesToProcess,
                                    std::vector<dataType> &newScalarsVector){
  ftm::FTMTree_MT *tree1 = mTree1->tree;

  // Create new tree
  MergeTree *mTreeNew = createEmptyMergeTree<dataType>(newScalarsVector.size());
  ftm::FTMTree_MT *treeNew = mTreeNew->tree;
  setTreeScalars<dataType>(mTreeNew, newScalarsVector);

  // Copy the old tree structure
  copyMergeTreeStructure(mTreeNew, tree1);
//std::cout << "street 1" << std::endl;
//mTreeNew->tree->printTree2();
  // Add nodes in the other trees
  auto nodesProcessed = addNodes<dataType>(mTreeNew, noTrees, nodesToProcess);
//std::cout << "street 2" << std::endl;
//mTreeNew->tree->printTree2();
  // Free old tree and replace it with new tree
  freeMergeTree<dataType>(mTree1);
  mTree1 = mTreeNew;

  return nodesProcessed;
}

// --------------------
// Utils
// --------------------
template <class dataType>
ftm::FTMTree_MT * makeFakeTree(dataType *nodesScalar, std::vector<SimplexId> &nodes,
                               std::vector<std::tuple<ftm::idNode, ftm::idNode>> &arcs){
  // Init Scalars
  ftm::Scalars *scalars = new ftm::Scalars();
  scalars->size = nodes.size();
  scalars->values = (void*)nodesScalar;
  
  // Init Tree
  ftm::Params *params = new ftm::Params();
  ftm::FTMTree_MT *tree = new ftm::FTMTree_MT(params, nullptr, scalars, ftm::Join_Split);
  tree->makeAlloc();
  
  // Add Nodes
  for(auto i : nodes)
    tree->makeNode(i);
  
  // Add Arcs
  for(std::tuple<ftm::idNode, ftm::idNode> arc : arcs)
    tree->makeSuperArc(std::get<0>(arc), std::get<1>(arc)); // (down, Up)
  
  return tree;
}

template<class dataType>
void verifyMergeTreeStructure(FTMTree_MT *tree){
  bool problem = false;
  
  bool isJT = isJoinTree<dataType>(tree);
  std::vector<std::tuple<ftm::idNode, ftm::idNode>> problemNodes;
  std::queue<ftm::idNode> queue;
  queue.push(getRoot(tree));
  while(!queue.empty()){
    ftm::idNode node = queue.front();
    queue.pop();
    
    if(!isRoot(tree, node)){
      bool thisProblem;
      if(isJT)
        thisProblem = tree->getValue<dataType>(node) > tree->getValue<dataType>(getParent(tree, node));
      else
        thisProblem = tree->getValue<dataType>(node) < tree->getValue<dataType>(getParent(tree, node));
      
      if(thisProblem)
        problemNodes.push_back(std::make_tuple(node, getParent(tree, node)));
      
      problem |= thisProblem;
    }
    
    for(auto c : getChildren(tree, node))
      queue.push(c);
  }
  
  if(problem){
    std::cout << "merge tree in input is not valid" << std::endl;
    for(auto tup : problemNodes)
      std::cout << std::get<0>(tup) << " _ " << std::get<1>(tup) << std::endl;
    printTree(tree);
    printTreeScalars<dataType>(tree);
    myPause();
  }
}

#endif
