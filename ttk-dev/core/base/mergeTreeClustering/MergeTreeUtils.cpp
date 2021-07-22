#include "MergeTreeUtils.h"

// ----------------------------------------
// Tree Functions
// ----------------------------------------

// --------------------
// Is / Get / Set / Delete
// --------------------

// --- Is

bool isRoot(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
  if(nodeId >= tree->getNumberOfNodes())
    std::cout << nodeId << " _____ " << tree->getNumberOfNodes() << std::endl;
  return tree->getNode(nodeId)->getNumberOfUpSuperArcs() == 0;
}

bool isLeaf(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
  return tree->getNode(nodeId)->getNumberOfDownSuperArcs() == 0;
}

bool isNodeAlone(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
  return isRoot(tree, nodeId) and isLeaf(tree, nodeId);
}

bool isNodeMerged(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
  bool merged = isNodeAlone(tree, nodeId) or isNodeAlone(tree, tree->getNode(nodeId)->getOrigin());
  auto nodeIdOrigin = tree->getNode(nodeId)->getOrigin();
  merged = merged or nodeIdOrigin == tree->getNode(nodeIdOrigin)->getOrigin();
  return merged;
}

bool isNodeOriginDefined(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
  unsigned int origin = (unsigned int)tree->getNode(nodeId)->getOrigin();
  return origin != ftm::nullNodes and origin < tree->getNumberOfNodes() and origin >= 0;;
}

bool isNodeIdInconsistent(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
  return nodeId >= tree->getNumberOfNodes() or nodeId < 0;
}

bool isThereOnlyOnePersistencePair(ftm::FTMTree_MT *tree){
  ftm::idNode treeRoot = getRoot(tree);
  unsigned int cptNodeAlone = 0;
  ftm::idNode otherNode = treeRoot;
  for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i)
    if(isNodeAlone(tree, i))
      cptNodeAlone++;
    else if(i != treeRoot)
      otherNode = i;
  //unsigned int origin = (unsigned int)tree->getNode(otherNode)->getOrigin();  
  ftm::idNode treeRootOrigin = tree->getNode(treeRoot)->getOrigin();
  return (otherNode != treeRoot and tree->getNumberOfNodes()-cptNodeAlone == 2 and
          (treeRootOrigin == otherNode or treeRootOrigin == treeRoot));
  /*return (otherNode != treeRoot and cptNodeAlone == tree->getNumberOfNodes()-2 and 
         (origin == treeRoot or otherNode == origin));*/
}

// Do not normalize node is if root or son of a merged root
bool notNeedToNormalize(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
  auto nodeIdParent = getParent(tree, nodeId);
  return isRoot(tree, nodeId) or ( isRoot(tree, nodeIdParent) and 
         nodeIdParent == (unsigned int)tree->getNode(nodeIdParent)->getOrigin() ); 
  // and nodeIdOrigin == nodeIdParent) ) 
}

bool isFullMerge(ftm::FTMTree_MT *tree){
  ftm::idNode treeRoot = getRoot(tree);
  return (unsigned int)tree->getNode(treeRoot)->getOrigin() == treeRoot;
}

bool isMultiPersPair(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
  auto nodeOriginOrigin = (unsigned int)tree->getNode(tree->getNode(nodeId)->getOrigin())->getOrigin();
  return nodeOriginOrigin != nodeId;
}

// nodeId must be a saddle
bool isBranchOrigin(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
  return getParent(tree, tree->getNode(nodeId)->getOrigin()) != nodeId;
}

// --- Get

int getNumberOfNodeAlone(ftm::FTMTree_MT *tree){
  int cpt = 0;
  for(ftm::idNode i = 0; i < tree->getNumberOfNodes(); ++i)
    cpt += isNodeAlone(tree, i) ? 1 : 0;
  return cpt;
}

int getRealNumberOfNodes(ftm::FTMTree_MT *tree){
  return tree->getNumberOfNodes() - getNumberOfNodeAlone(tree);
}

ftm::idNode getParent(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
  if(! isRoot(tree, nodeId)){
    // _ Nodes in merge trees should have only one parent
    ftm::idSuperArc arcId = tree->getNode(nodeId)->getUpSuperArcId(0);
    ftm::idNode parentNodeId = tree->getSuperArc(arcId)->getUpNodeId();
    return parentNodeId;
  }
  return nodeId;
}

ftm::idNode getRoot(ftm::FTMTree_MT *tree){
  for(ftm::idNode node = 0; node < tree->getNumberOfNodes(); ++node)
    if(isRoot(tree, node) and !isLeaf(tree, node))
      return node;
  return ftm::nullNodes;
}

std::vector<ftm::idNode> getAllRoots(ftm::FTMTree_MT *tree){
  std::vector<ftm::idNode> roots;
  for(ftm::idNode node = 0; node < tree->getNumberOfNodes(); ++node)
    if(isRoot(tree, node) and !isLeaf(tree, node))
      roots.push_back(node);
  return roots;
}

int getNumberOfRoot(ftm::FTMTree_MT *tree){
  int noRoot = 0;
  for(ftm::idNode node = 0; node < tree->getNumberOfNodes(); ++node)
    if(isRoot(tree, node) and !isLeaf(tree, node))
      ++noRoot;
  return noRoot;
}

// Get childrens (idNode) of a node in the tree
std::vector<ftm::idNode> getChildren(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
  std::vector<ftm::idNode> childrens;
  for(ftm::idSuperArc i = 0; i < tree->getNode(nodeId)->getNumberOfDownSuperArcs(); ++i){
    ftm::idSuperArc arcId = tree->getNode(nodeId)->getDownSuperArcId(i);
    childrens.push_back(tree->getSuperArc(arcId)->getDownNodeId()); 
  }
  return childrens;
}

int getNumberOfChildren(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
  return tree->getNode(nodeId)->getNumberOfDownSuperArcs();
}

std::vector<ttk::ftm::idNode> getLeaves(ttk::ftm::FTMTree_MT *tree){
  std::vector<ttk::ftm::idNode> treeLeaves;
  for(ttk::ftm::idNode i = 0; i < tree->getNumberOfNodes(); ++i){
    if(isLeaf(tree, i) and !isRoot(tree, i))
      treeLeaves.push_back(i);
  }
  return treeLeaves;
}

int getNumberOfLeaves(ftm::FTMTree_MT *tree){
  auto leaves = getLeaves(tree);
  return leaves.size();
}

int getTreeDepth(ftm::FTMTree_MT *tree){
  int maxDepth=0;
  std::queue<std::tuple<ftm::idNode, int>> queue;
  queue.push(std::make_tuple(getRoot(tree), 0));
  while(!queue.empty()){
    auto tup = queue.front();
    queue.pop();
    ftm::idNode node = std::get<0>(tup);
    int depth = std::get<1>(tup);
    maxDepth = std::max(maxDepth, depth);
    for(ftm::idNode child : getChildren(tree, node))
      queue.push(std::make_tuple(child, depth+1));
  }
  return maxDepth;
}

int getNodeLevel(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
  int level = 0;
  auto root = getRoot(tree);
  int noRoot = getNumberOfRoot(tree);
  if(noRoot != 1){
    std::cout << "problem, there is " << noRoot << " root(s)" << std::endl;
    tree->printTree2();
    printTree(tree);
    myPause();
  }
  if(isNodeAlone(tree, nodeId))
    return 0;
  while(nodeId != root){
    nodeId = getParent(tree, nodeId);
    ++level;
  }
  return level;
}

std::vector<int> getAllNodeLevel(ftm::FTMTree_MT *tree){
  std::vector<int> allNodeLevel(tree->getNumberOfNodes());
  std::queue<std::tuple<ftm::idNode, int>> queue;
  queue.push(std::make_tuple(getRoot(tree), 0));
  while(!queue.empty()){
    auto tup = queue.front();
    queue.pop();
    ftm::idNode node = std::get<0>(tup);
    int level = std::get<1>(tup);
    allNodeLevel[node] = level;
    for(ftm::idNode child : getChildren(tree, node))
      queue.push(std::make_tuple(child, level+1));
  }
  return allNodeLevel;
}

void getTreeBranching(FTMTree_MT *tree, std::vector<idNode> &branching, std::vector<int> &branchingID,
                      std::vector<std::vector<ftm::idNode>> &nodeBranching){
  branching = std::vector<idNode>(tree->getNumberOfNodes());
  branchingID = std::vector<int>(tree->getNumberOfNodes(), -1);
  nodeBranching = std::vector<std::vector<ftm::idNode>>(tree->getNumberOfNodes());
  int branchID = 0;
  std::queue<idNode> queue;
  queue.push(getRoot(tree));
  while(!queue.empty()){
    idNode node = queue.front();
    queue.pop();
    if(isLeaf(tree, node))
      continue;
    auto nodeOrigin = tree->getNode(node)->getOrigin();
    idNode parentNodeOrigin = nodeOrigin;
    idNode oldParentNodeOrigin;
    while(parentNodeOrigin != node){
      branching[parentNodeOrigin] = node;
      branchingID[parentNodeOrigin] = branchID;
      nodeBranching[node].push_back(parentNodeOrigin);
      oldParentNodeOrigin = parentNodeOrigin;
      parentNodeOrigin = getParent(tree, parentNodeOrigin);
      if(oldParentNodeOrigin == parentNodeOrigin){
        printTree(tree);
        std::cout << getRoot(tree) << " _ " << parentNodeOrigin << std::endl;
        std::cout << node << " _ " << nodeOrigin << std::endl;
        std::cout << "getTreeBranching" << std::endl;
        myPause();
      }
    }
    if(isRoot(tree, node)){
      branching[node] = node;
      branchingID[node] = branchID;
    }
    ++branchID;
    for(auto child : getChildren(tree, node))
      queue.push(child);
  }
}

void getTreeBranching(FTMTree_MT *tree, std::vector<idNode> &branching, std::vector<int> &branchingID){
  std::vector<std::vector<ftm::idNode>> nodeBranching;
  getTreeBranching(tree, branching, branchingID, nodeBranching);
}

// Node must be a non-leaf node
std::tuple<std::vector<ftm::idNode>, std::vector<ftm::idNode>> 
getBranchOriginsFromThisBranch(FTMTree_MT *tree, ftm::idNode node){
  std::vector<ftm::idNode> branchOrigins, nonBranchOrigins;
  
  ftm::idNode nodeOrigin = tree->getNode(node)->getOrigin();
  ftm::idNode nodeParent = getParent(tree, nodeOrigin);
  while(nodeParent != node){
    if(isBranchOrigin(tree, nodeParent))
      branchOrigins.push_back(nodeParent);
    else
      nonBranchOrigins.push_back(nodeParent);
    nodeParent = getParent(tree, nodeParent);
  }
  
  return std::make_tuple(branchOrigins, nonBranchOrigins);
}

std::vector<ftm::idNode> getBranchSubtree(std::vector<ftm::idNode> &branching, FTMTree_MT *tree, 
                                          ftm::idNode branchRoot){
  std::vector<ftm::idNode> branchSubtree;
  std::queue<ftm::idNode> queue;
  queue.push(branchRoot);
  while(!queue.empty()){
    ftm::idNode node = queue.front();
    queue.pop();
    
    if(branching[node] != branchRoot and getParent(tree, node) == branchRoot and node != branchRoot)
      continue;
    
    branchSubtree.push_back(node);
    for(auto child : getChildren(tree, node))
      queue.push(child);
  }
  return branchSubtree;
}

double getValue(FTMTree_MT *tree, ftm::idNode node, int dataType){
  double res = 0;
  switch(dataType){
    vtkTemplateMacro(res = tree->getValue<VTK_TT>(node));
  }
  return res;
}

double getPersistence(FTMTree_MT *tree, ftm::idNode node, int dataType){
  double res = 0;
  switch(dataType){
    vtkTemplateMacro(res = getNodePersistence<VTK_TT>(tree, node));
  }
  return res;
}

// --- Set

void setParent(ftm::FTMTree_MT *tree, ftm::idNode nodeId, ftm::idNode newParentNodeId){
  deleteParent(tree, nodeId);
  tree->makeSuperArc(nodeId, newParentNodeId);
}

// --- Delete

void deleteIthUpArc(ftm::FTMTree_MT *tree, ftm::idNode nodeId, int arcIth){
  ftm::idSuperArc nodeArcId = tree->getNode(nodeId)->getUpSuperArcId(arcIth);
  // Delete down arc from old parent
  ftm::idNode parentNodeId = tree->getSuperArc(nodeArcId)->getUpNodeId();
  tree->getNode(parentNodeId)->removeDownSuperArc(nodeArcId);
  // Delete up arc from node
  tree->getNode(nodeId)->removeUpSuperArc(nodeArcId);
}

void deleteParent(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
  if(! isRoot(tree, nodeId)){
    // _ Nodes in trees should have only one parent
    deleteIthUpArc(tree, nodeId, 0);
  }
}

// Delete node by keeping subtree 
void deleteNode(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
  if(isRoot(tree, nodeId) and !isLeaf(tree, nodeId))
    std::cout << "deletion of root!" << std::endl;
  
  ftm::idNode parentNodeId = getParent(tree, nodeId);
  if(! isRoot(tree, nodeId)){
    // Delete down arc from parent node
    // _ Nodes in trees should have only one parent
    ftm::idSuperArc nodeArcId = tree->getNode(nodeId)->getUpSuperArcId(0);
    tree->getNode(parentNodeId)->removeDownSuperArc(nodeArcId);
  }
  // Delete up arc from child nodes
  for(ftm::idSuperArc i = 0; i < tree->getNode(nodeId)->getNumberOfDownSuperArcs(); ++i){
    ftm::idSuperArc arcId = tree->getNode(nodeId)->getDownSuperArcId(i);
    ftm::idNode childNodeId = tree->getSuperArc(arcId)->getDownNodeId();
    tree->getNode(childNodeId)->removeUpSuperArc(arcId);
    if(! isRoot(tree, nodeId))
      tree->makeSuperArc(childNodeId, parentNodeId);
  }
  // Reset deleted node
  tree->getNode(nodeId)->clearDownSuperArcs();
  tree->getNode(nodeId)->clearUpSuperArcs();
}

void deleteSubtree(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
  std::queue<ftm::idNode> queue;
  queue.push(nodeId);
  while(!queue.empty()){
    ftm::idNode node = queue.front();
    queue.pop();
    auto children = getChildren(tree, node);
    for(auto child : children)
      queue.push(child);
    deleteNode(tree, node);
  }
}

// --------------------
// Persistence
// --------------------
std::vector<std::vector<ftm::idNode>> getMultiPersOriginsVectorFromTree(ftm::FTMTree_MT *tree){
  std::vector<std::vector<ftm::idNode>> treeMultiPers(tree->getNumberOfNodes());
  for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i)
    if(isLeaf(tree, i) and isNodeOriginDefined(tree, i) and not isNodeAlone(tree, i) and
        tree->getNode(tree->getNode(i)->getOrigin())->getOrigin() != (int)i)
      treeMultiPers[tree->getNode(i)->getOrigin()].push_back(i);
  return treeMultiPers;
}

std::stringstream printMultiPersOriginsVectorFromTree(ftm::FTMTree_MT *tree, bool doPrint){
  std::stringstream ss;
  auto vec = getMultiPersOriginsVectorFromTree(tree);
  for(unsigned int i = 0; i < vec.size(); ++i)
    if(vec[i].size() != 0){
      ss << i << " : ";
      for(auto t : vec[i])
        ss << t << " ";
      ss << std::endl;
    }
    
  if(doPrint)
    std::cout << ss.str();
  return ss;
}

// --------------------
// Create/Delete/Modify Tree
// --------------------

void copyMergeTreeStructure(MergeTree *mergeTree, ftm::FTMTree_MT *tree){
  ftm::FTMTree_MT *treeNew = mergeTree->tree;
  
  // Add Nodes
  for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i)
    treeNew->makeNode(i);
  for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i)
    treeNew->getNode(i)->setOrigin(tree->getNode(i)->getOrigin());

  // Add Arcs
  ftm::idNode root = getRoot(tree);
  std::stack<ftm::idNode> stackNodes;
  stackNodes.push(root);
  while(!stackNodes.empty()){
    ftm::idNode node = stackNodes.top();
    stackNodes.pop();
    for(unsigned int i = 0; i < tree->getNode(node)->getNumberOfDownSuperArcs(); ++i){
      ftm::idSuperArc arcId = tree->getNode(node)->getDownSuperArcId(i);
      ftm::idNode childNodeId = tree->getSuperArc(arcId)->getDownNodeId();
      treeNew->makeSuperArc(childNodeId, node); // (down, Up)
      stackNodes.push(childNodeId);
    }
  }
}

void freeTree(ftm::FTMTree_MT *tree){
  delete tree;
}

std::vector<ftm::FTMTree_MT*> mergeTreeToFTMTree(std::vector<MergeTree*> &trees){
  std::vector<ftm::FTMTree_MT*> treesT;
  for(auto t : trees)
    treesT.push_back(t->tree);
  return treesT;
}

// --------------------
// Utils
// --------------------

void myPause(){ 
  std::cout << "pause" << std::endl; 
  int temp; 
  std::cin >> temp; 
}

void printNode(ftm::FTMTree_MT *tree, ftm::idNode node, std::stringstream &ss){
  ss << "(" << node << ") \\ ";
  for(auto child : getChildren(tree, node))
    ss << "+" << child << " ";

  if(!isRoot(tree, node))
    ss << " / +" << getParent(tree, node);
  ss << std::endl;
}

std::stringstream printTree(ftm::FTMTree_MT *tree, bool doPrint){
  std::stringstream ss;
  std::vector<ftm::idNode> allRoots = getAllRoots(tree);
  if(allRoots.size() != 1)
    ss << allRoots.size() << " roots" << std::endl;
  for(unsigned int i = 0; i < allRoots.size(); ++i){
    if(allRoots.size() != 1)
      ss << i << " _ ";
    ss << "Nodes----------" << std::endl;
    std::queue<ftm::idNode> queue;
    queue.push(allRoots[i]);
    while(!queue.empty()){
      ftm::idNode node = queue.front();
      queue.pop();

      printNode(tree, node, ss);

      for(auto child : getChildren(tree, node))
        queue.push(child);
    }
  }
  if(allRoots.size() == 0)
    for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i)
      //if(not isNodeAlone(tree, i))
        printNode(tree, i, ss);

  if(doPrint)
    std::cout << ss.str();
  return ss;
}

void printTree(MergeTree *tree, bool doPrint){
  printTree(tree->tree, doPrint);
}

void printTreeStats(ftm::FTMTree_MT * tree){
  auto noNodesT = tree->getNumberOfNodes();
  auto noNodes = getRealNumberOfNodes(tree);
  std::cout << "tree [node: " << noNodes << " / "<< noNodesT;
  std::cout << ", depth: " << getTreeDepth(tree) << "]" << std::endl;
}

void printTreeStats(MergeTree * tree){
  printTreeStats(tree->tree);
}

void printTreeStats(std::vector<ftm::FTMTree_MT *> &trees){
  int avgNodes = 0, avgNodesT = 0;
  double avgDepth=0;
  for(unsigned int i = 0; i < trees.size(); ++i){
    auto noNodesT = trees[i]->getNumberOfNodes();
    auto noNodes = getRealNumberOfNodes(trees[i]);
    avgNodes += noNodes;
    avgNodesT += noNodesT;
    avgDepth += getTreeDepth(trees[i]);
    /*printMultiPersPairsFromTree<dataType>(trees[i], BranchDecomposition);
    printPairsFromTree<dataType>(trees[i], BranchDecomposition);*/
  }
  avgNodes /= trees.size();
  avgNodesT /= trees.size();
  avgDepth /= trees.size();
  std::cout << trees.size() << " trees average [node: " << avgNodes << " / "<< avgNodesT;
  std::cout << ", depth: " << avgDepth << "]" << std::endl;
}

void printTreeStats(std::vector<MergeTree *> &trees){
  auto treesT = mergeTreeToFTMTree(trees);
  printTreeStats(treesT);
}
