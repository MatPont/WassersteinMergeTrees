#include "FTMTreeUtils.h"

// ----------------------------------------
// Tree Functions
// ----------------------------------------

// --------------------
// Is / Get / Set / Delete
// --------------------

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

void setParent(ftm::FTMTree_MT *tree, ftm::idNode nodeId, ftm::idNode newParentNodeId){
  deleteParent(tree, nodeId);
  tree->makeSuperArc(nodeId, newParentNodeId);
}

ftm::idNode getRoot(ftm::FTMTree_MT *tree){
  for(ftm::idNode node = 0; node < tree->getNumberOfNodes(); ++node){
    if(isRoot(tree, node) and !isLeaf(tree, node))
      return node;
  }
  return ftm::nullNodes;
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

// Delete node by keeping subtree 
void deleteNode(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
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

bool isNodeIdInconsistent(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
  return nodeId > tree->getNumberOfNodes() or nodeId < 0;
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
  unsigned int origin = (unsigned int)tree->getNode(otherNode)->getOrigin();
  return (otherNode != treeRoot and cptNodeAlone == tree->getNumberOfNodes()-2 and 
         (origin == treeRoot or otherNode == origin));
}

// Do not normalize node is if root or son of a merged root
bool notNeedToNormalize(ftm::FTMTree_MT *tree, ftm::idNode nodeId){
  auto nodeIdParent = getParent(tree, nodeId);
  return isRoot(tree, nodeId) or ( isRoot(tree, nodeIdParent) and 
         nodeIdParent == (unsigned int)tree->getNode(nodeIdParent)->getOrigin() ); 
  // and nodeIdOrigin == nodeIdParent) ) 
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

// --------------------
// Utils
// --------------------

void myPause(){ std::cout<<"pause"<<std::endl; int temp; std::cin>>temp; }
