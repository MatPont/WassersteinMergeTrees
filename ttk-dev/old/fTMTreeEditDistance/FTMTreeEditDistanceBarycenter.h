///
/// \ingroup base
/// \class ttk::FTMTreeEditDistanceBarycenter
/// \author Mathieu Pont <mathieu.pont@outlook.com>
/// \date 2020.
///
/// This module defines the %FTMTreeEditDistanceBarycenter class that computes the 
/// the barycenter of merge trees according the edit distance.
///

#ifndef _FTMTREEEDITDISTANCEBARYCENTER_H
#define _FTMTREEEDITDISTANCEBARYCENTER_H

#pragma once

#include <random>

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>

#include "FTMTreeEditDistanceBase.h"

namespace ttk {
  
    /**
    * The FTMTreeEditDistanceBarycenter class that computes the 
    * the barycenter merge trees according the edit distance.
    */
    class FTMTreeEditDistanceBarycenter : virtual public Debug, public FTMTreeEditDistanceBase {
      
    protected:
      double Tol = 0.0;
      bool AddNodes = true;
      bool Deterministic = true;
      bool IsCalled = false;

    public:
        FTMTreeEditDistanceBarycenter() {
          this->setDebugMsgPrefix(
              "FTMTreeEditDistanceBarycenter"); // inherited from Debug: prefix will be printed at the
          // beginning of every msg
          omp_set_nested(1);
        };
        ~FTMTreeEditDistanceBarycenter(){};

        void setTol(double tolT){
          Tol = tolT;
        }
        
        void setAddNodes(bool addNodesT){
          AddNodes = addNodesT;
        }
        
        void setDeterministic(bool deterministicT){
          Deterministic = deterministicT;
        }
        
        /**
        * Implementation of the algorithm.
        */
        int getBestInitTreeIndex(std::vector<ftm::FTMTree_MT *> trees){
          int bestIndex = -1;
          unsigned int bestValue = 0;
          std::vector<int> sizes(trees.size());
          for(int i = 0; i < trees.size(); ++i){
            if(getRealNumberOfNodes(trees[i]) > bestValue){
              bestIndex = i;
              bestValue = getRealNumberOfNodes(trees[i]);
            }
            sizes[i] = getRealNumberOfNodes(trees[i]);
          }
          if(not Deterministic){
            std::random_device rd;
            std::default_random_engine generator(rd());
            std::discrete_distribution<int> distribution(sizes.begin(), sizes.end());
            bestIndex = distribution(generator);
          }
          return bestIndex;
        }
        
        template <class dataType>
        MergeTree* initBarycenterTree(std::vector<ftm::FTMTree_MT *> trees){
          int bestIndex = getBestInitTreeIndex(trees);
          //bestIndex = 0;
          MergeTree *baryTree = copyMergeTree<dataType>(trees[bestIndex], true);

          return baryTree;
        }
        
        template<class dataType>
        void updateBarycenterTreeStructure(std::vector<ftm::FTMTree_MT *> &trees, MergeTree *&baryMergeTree,
                                  std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>> &matchings){
          ftm::FTMTree_MT *baryTree = baryMergeTree->tree;
          
          // Init matching matrix 
          // m[i][j] contains the node in the barycenter matched to the jth node of the ith tree
          std::vector<std::vector<ftm::idNode>> matrixMatchings(trees.size());
          std::vector<bool> baryMatched(baryTree->getNumberOfNodes(), false);
          for(int i = 0; i < matchings.size(); ++i){
            auto matching = matchings[i];
            std::vector<ftm::idNode> matchingT(trees[i]->getNumberOfNodes(), -1);
            for(auto match : matching){
              matchingT[std::get<1>(match)] = std::get<0>(match);
              baryMatched[std::get<0>(match)] = true;
            }
            matrixMatchings[i].insert(matrixMatchings[i].end(), matchingT.begin(), matchingT.end());
          }
          
          // Iterate through trees to get the nodes to add in the barycenter
          std::vector<std::vector<ftm::idNode>> nodesToAdd(trees.size());
          for(int i = 0; i < trees.size(); ++i){
            ftm::idNode root = getRoot(trees[i]);
            std::queue<ftm::idNode> queue;
            queue.push(root);
            while(!queue.empty()){
              ftm::idNode node = queue.front();
              queue.pop();
              bool processChildren = true;
              if(matrixMatchings[i][node] == -1){ // if node in trees[i] is not matched
                if(not KeepSubtree){
                  processChildren = false;
                  nodesToAdd[i].push_back(node);
                }else{
                  // TODO manage if keepSubtree=true
                  std::cout << "ERROR barycenter with KeepSubtree=true is not implemented yet" << std::endl;
                }
              }
              if(processChildren){
                auto children = getChildren(trees[i], node);
                for(auto child : children)
                  if(not (isThereOnlyOnePersistencePair(trees[i]) and isLeaf(trees[i], child)) )
                    queue.push(child);
              }
            }
          }
          
          // Delete nodes that are not matched in the barycenter
/*printNode<dataType>(baryTree, getRoot(baryTree));
std::cout << "matching :" << std::endl;
for(int i = 0; i < matchings.size(); ++i){
  std::cout << matchings[i].size() << " / " << getRealNumberOfNodes(trees[i]) << std::endl;
  printNode<dataType>(trees[i], getRoot(trees[i]));
}
std::cout << "matching end." << std::endl;
std::cout << getRealNumberOfNodes(baryTree) << " / " << baryTree->getNumberOfNodes() << " (" << trees.size() << ")" << std::endl;*/
          for(int i = 0; i < baryTree->getNumberOfNodes(); ++i)
            if(not baryMatched[i])
              deleteNode(baryTree, i);
//std::cout << getRealNumberOfNodes(baryTree) << " / " << baryTree->getNumberOfNodes() << " (" << trees.size() << ")" << std::endl;
          if(not KeepSubtree){
            // Add scalars and nodes not present in the barycenter
            ftm::idNode nodeCpt = baryTree->getNumberOfNodes();
            std::vector<std::tuple<ftm::idNode, ftm::idNode, int>> nodesToProcess;
            std::vector<dataType> newScalarsVector = getTreeScalars<dataType>(baryMergeTree);
            for(int i = 0; i < nodesToAdd.size(); ++i){
//std::cout << "tree " << i << std::endl;
              for(auto node : nodesToAdd[i]){
                ftm::idNode parent = matrixMatchings[i][getParent(trees[i], node)];
if(isNodeAlone(baryTree, parent)){
  std::cout << "parent " << parent << std::endl;
  myPause();
}
                std::vector<dataType> addedScalars;
                nodeCpt = getNodesAndScalarsToAdd<dataType>(baryMergeTree, parent, trees[i], node, 
                                                            addedScalars, nodesToProcess, nodeCpt, i);
                newScalarsVector.insert(newScalarsVector.end(), addedScalars.begin(), addedScalars.end());
              }
            }
            if(AddNodes){
              auto nodesProcessed = updateNodesAndScalars<dataType>(baryMergeTree, trees.size(), 
                                                                    nodesToProcess, newScalarsVector);
              for(int i = 0; i < matchings.size(); ++i)
                matchings[i].insert(matchings[i].end(), nodesProcessed[i].begin(), nodesProcessed[i].end());
            }
          }else{
            // TODO manage if keepSubtree=true
            std::cout << "ERROR barycenter with KeepSubtree=true is not implemented yet" << std::endl;
          }          
        }
        
        template<class dataType>
        std::tuple<double, double> getParametrizedBirthDeath(ftm::FTMTree_MT *tree1, ftm::idNode nodeId1,
                                                                 ftm::FTMTree_MT *tree2=nullptr, 
                                                                 ftm::idNode nodeId2=nullptr){
          std::tuple<double, double> birthDeath;
          // Normalized Wasserstein
          if(NormalizedWasserstein and not RescaledWasserstein)
            birthDeath = getNormalizedBirthDeathDouble<dataType>(tree1, nodeId1);
          // Rescaled Wasserstein
          else if(NormalizedWasserstein and RescaledWasserstein) 
            birthDeath = getRescaledBirthDeath<dataType>(tree1, nodeId1, tree2, nodeId2);
          // Classical Wasserstein
          else 
            birthDeath = getBirthDeath<dataType>(tree1, nodeId1);
          return birthDeath;
        }

        template<class dataType>
        std::tuple<dataType, dataType> getParametrizedBirthDeathFromVector(ftm::FTMTree_MT *tree1, 
                                              ftm::idNode nodeId1, ftm::FTMTree_MT *tree2, 
                                              ftm::idNode nodeId2, std::vector<dataType> &newScalarsVector){
          if(NormalizedWasserstein and RescaledWasserstein) 
            return getRescaledBirthDeathFromVector<dataType>(tree1, nodeId1, tree2, nodeId2, 
                                                             newScalarsVector);
          return getParametrizedBirthDeath<dataType>(tree1, nodeId1, tree2, nodeId2);
        }
        
        template<class dataType>
        std::tuple<dataType, dataType> interpolation(MergeTree *baryMergeTree, ftm::idNode nodeId,
                                                     std::vector<dataType> &newScalarsVector,
                                                     std::vector<ftm::FTMTree_MT *> &trees, 
                                                     std::vector<ftm::idNode> &nodes, 
                                                     std::vector<double> &alphas){
          ftm::FTMTree_MT *baryTree = baryMergeTree->tree;
          dataType mu_max = getMinMaxLocalFromVector<dataType>(baryTree, nodeId, newScalarsVector, false);
          dataType mu_min = getMinMaxLocalFromVector<dataType>(baryTree, nodeId, newScalarsVector);
          double newBirth=0, newDeath=0;
          
          // Compute projection 
          // TODO manage when alpha_i != 1/n
          double tempBirth=0, tempDeath=0;
          int offDiagonal = 0;
          for(int i = 0; i < trees.size(); ++i){
            if(nodes[i] != -1){ // if node is matched in trees[i]
              auto iBirthDeath = getParametrizedBirthDeathFromVector<dataType>(trees[i], nodes[i], baryTree, 
                                                                               nodeId, newScalarsVector);
              auto newMinMax = getNewMinMaxFromVector<dataType>(trees[i], nodes[i], baryTree, nodeId,
                                                                newScalarsVector);
              tempBirth += std::get<0>(iBirthDeath);
              tempDeath += std::get<1>(iBirthDeath);
              if(NormalizedWasserstein and RescaledWasserstein){
                auto newMinMax = getNewMinMaxFromVector<dataType>(trees[i], nodes[i], baryTree, nodeId,
                                                                newScalarsVector);
                tempBirth /= (std::get<1>(newMinMax) - std::get<0>(newMinMax));
                tempDeath /= (std::get<1>(newMinMax) - std::get<0>(newMinMax));
              }
              ++offDiagonal;
            }
          }
          tempBirth /= offDiagonal;
          tempDeath /= offDiagonal;
          double projec = (tempBirth+tempDeath)/2;
          //auto birthDeath = getParametrizedBirthDeath<dataType>(baryTree, nodeId, baryTree, nodeId);
          //dataType projec = (std::get<0>(birthDeath)+std::get<1>(birthDeath))/2;
//std::cout << "==============" << std::endl;
//std::cout << std::get<0>(birthDeath) << " _ " << std::get<1>(birthDeath) << " _ " << projec << std::endl;
          
          // Compute newBirth and newDeath
          dataType divisor = 0;
          for(int i = 0; i < trees.size(); ++i){
            double iBirth = projec, iDeath = projec;
            if(nodes[i] != -1){ // if node is matched in trees[i]
              auto iBirthDeath = getParametrizedBirthDeathFromVector<dataType>(trees[i], nodes[i], baryTree, 
                                                                               nodeId, newScalarsVector);
              iBirth = std::get<0>(iBirthDeath);
              iDeath = std::get<1>(iBirthDeath);
            }
//std::cout << "i" << i << " __ " << iBirth << " _ " << iDeath << std::endl;
            if(NormalizedWasserstein and RescaledWasserstein){
              dataType beta_max=1, beta_min=0;
              auto newMinMax = (nodes[i]==-1) ? getNewMinMax<dataType>(baryTree, nodeId, baryTree, nodeId):
                                               getNewMinMaxFromVector<dataType>(trees[i], nodes[i], baryTree,
                                                                                nodeId, newScalarsVector);
              if(nodes[i] == -1){
                beta_max = mu_max;
                beta_min = mu_min;
                iBirth *= (beta_max - beta_min);// / (std::get<1>(newMinMax) - std::get<0>(newMinMax));
                iDeath *= (beta_max - beta_min);// / (std::get<1>(newMinMax) - std::get<0>(newMinMax));
              }else{
                beta_min = std::get<0>(newMinMax);
                beta_max = std::get<1>(newMinMax);
              }
              iBirth *= (beta_max - beta_min);
              iDeath *= (beta_max - beta_min);
              divisor += alphas[i] * (beta_max - beta_min) * (beta_max - beta_min);
//std::cout << "i" << i << " __ " << beta_max << " _ " << beta_min << " _ " << alphas[i] << std::endl;
            }
            newBirth += alphas[i] * iBirth;
            newDeath += alphas[i] * iDeath;
          }          
          if(NormalizedWasserstein and RescaledWasserstein){
            newBirth /= divisor;
            newDeath /= divisor;
//std::cout << "divisor : " << divisor << std::endl;
          }
          if(NormalizedWasserstein or RescaledWasserstein){
            newBirth = newBirth * (mu_max - mu_min) + mu_min;
            newDeath = newDeath * (mu_max - mu_min) + mu_min;
//std::cout << "mu_max mu_min : " << mu_max << " _ " << mu_min << std::endl;
          }
          
//std::cout << newBirth << " _ " << newDeath << std::endl;
//int temp; std::cin >> temp;

          dataType newBirthT = newBirth;
          dataType newDeathT = newDeath;
          return std::make_tuple(newBirth, newDeath);
        }
        
        template<class dataType>
        std::tuple<dataType, dataType> interpolationAdded(ftm::FTMTree_MT * tree, ftm::idNode nodeId, 
                                                          double alpha, MergeTree *baryMergeTree, 
                                                          ftm::idNode nodeB, 
                                                          std::vector<dataType> &newScalarsVector){
          ftm::FTMTree_MT *baryTree = baryMergeTree->tree;
          dataType mu_max = getMinMaxLocalFromVector<dataType>(baryTree, nodeB, newScalarsVector, false);
          dataType mu_min = getMinMaxLocalFromVector<dataType>(baryTree, nodeB, newScalarsVector);
          
          auto birthDeath = getParametrizedBirthDeathFromVector<dataType>(tree, nodeId, baryTree, nodeB, 
                                                                          newScalarsVector);
          double newBirth = std::get<0>(birthDeath);
          double newDeath = std::get<1>(birthDeath);
          double projec = (newBirth+newDeath)/2;
//std::cout << "=_=_=_=_=_=_=_" << std::endl;
//std::cout << newBirth << " _ " << newDeath << " _ " << projec << " _ " << alpha << std::endl;
          
          dataType beta_min=0, beta_max=0, divisor=1;
          if(NormalizedWasserstein and RescaledWasserstein){
            auto newMinMax = getNewMinMaxFromVector<dataType>(tree, nodeId, baryTree, nodeB, 
                                                              newScalarsVector);
            beta_min = std::get<0>(newMinMax);
            beta_max = std::get<1>(newMinMax);
            newBirth *= (beta_max - beta_min);
            newDeath *= (beta_max - beta_min);
            /*projec = projec * (beta_max - beta_min);
            divisor = (beta_max - beta_min) * (beta_max - beta_min);*/
            projec = projec * (mu_max - mu_min) * (mu_max - mu_min) / (beta_max - beta_min);
            divisor = alpha * (beta_max - beta_min) * (beta_max - beta_min) + 
                      (1-alpha) * (mu_max - mu_min) * (mu_max - mu_min);
//std::cout << "beta_max beta_min : " << beta_max << " _ " << beta_min << std::endl;
          }
          
          newBirth = alpha * newBirth + (1-alpha) * projec;
          newDeath = alpha * newDeath + (1-alpha) * projec;
//std::cout << newBirth << " _ " << newDeath << std::endl;

          if(NormalizedWasserstein and RescaledWasserstein){
            newBirth /= divisor;
            newDeath /= divisor;
//std::cout << "divisor : " << divisor << std::endl;
          }
          
          if(NormalizedWasserstein or RescaledWasserstein){
            newBirth = newBirth * (mu_max - mu_min) + mu_min;
            newDeath = newDeath * (mu_max - mu_min) + mu_min;
//std::cout << "mu_max mu_min : " << mu_max << " _ " << mu_min << std::endl;
          }
//std::cout << newBirth << " _ " << newDeath << std::endl;
          
          dataType newBirthT = newBirth;
          dataType newDeathT = newDeath;
          return std::make_tuple(newBirthT, newDeathT);
        }
        
        template <class dataType>
        void purgeBarycenter(MergeTree *&baryMergeTree, std::vector<std::vector<ftm::idNode>> &baryMatching,
                             std::vector<ftm::FTMTree_MT *> &trees, std::vector<double> &alphas){
          ftm::FTMTree_MT *baryTree = baryMergeTree->tree;
          std::vector<bool> nodesProcessed(baryTree->getNumberOfNodes(), false);
          std::vector<dataType> nodesMatchingCost(baryTree->getNumberOfNodes(), 0);
          std::vector<dataType> nodesDestructCost(baryTree->getNumberOfNodes(), 0);
          auto leaves = getLeaves(baryTree);
          std::queue<ftm::idNode> queue;
          for(auto leaf : leaves)
            queue.push(leaf);
          while(!queue.empty()){
            ftm::idNode node = queue.front();
            queue.pop();
            for(int i = 0; i < trees.size(); ++i){
              dataType newMatchingCost = alphas[i];
              dataType newDestructCost = alphas[i];
              if(baryMatching[node][i] != -1){
                newMatchingCost *= relabelCost<dataType>(baryTree, node, trees[i], baryMatching[node][i]);
                newDestructCost *= deleteCost<dataType>(trees[i], baryMatching[node][i]);
              }else{
                newMatchingCost *= deleteCost<dataType>(baryTree, node);
                newDestructCost *= 0;
              }
              nodesMatchingCost[node] += newMatchingCost;
              nodesDestructCost[node] += newDestructCost;
            }
            auto children = getChildren(baryTree, node);
            for(auto child : children){
              nodesMatchingCost[node] += nodesMatchingCost[child];
              nodesDestructCost[node] += nodesDestructCost[child];
            }
            nodesProcessed[node] = true;
            if(not nodesProcessed[getParent(baryTree, node)])
              queue.push(getParent(baryTree, node));
            
            // Destruct subtree if better
            //std::cout << nodesDestructCost[node] << " _ " << 
            if(nodesDestructCost[node] < nodesMatchingCost[node]){
              deleteSubtree(baryTree, node);
              nodesDestructCost[node] = 0;
              nodesMatchingCost[node] = 0;
            }
          }
        }
        
        template<class dataType>
        void updateBarycenterTreeScalars(std::vector<ftm::FTMTree_MT *> &trees, MergeTree *&baryMergeTree,
                                  std::vector<double> &alphas, int indexAddedNodes,
                                  std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>> &matchings){
//std::cout << "updateBarycenterTreeScalars" << std::endl;
          ftm::FTMTree_MT *baryTree = baryMergeTree->tree;
          
          // Init matching matrix
          // m[i][j] contains the node in trees[j] matched to the node i in the barycenter
          //auto maxi = std::numeric_limits<ftm::idNode>::max();          
          std::vector<std::vector<ftm::idNode>> baryMatching(baryTree->getNumberOfNodes(), 
                                                             std::vector<ftm::idNode>(trees.size(), -1));
          std::vector<int> nodesAddedTree(baryTree->getNumberOfNodes(), -1);
          for(int i = 0; i < matchings.size(); ++i){
            auto matching = matchings[i];
            for(auto match : matching){
              baryMatching[std::get<0>(match)][i] = std::get<1>(match);
              if(std::get<0>(match) >= indexAddedNodes) // get the tree of this added node
                nodesAddedTree[std::get<0>(match)] = i;
            }
          }
          
//int cptBug = 0;
          // Interpolate scalars
          std::vector<dataType> newScalarsVector(baryTree->getNumberOfNodes());
          ftm::idNode root = getRoot(baryTree);
          std::queue<ftm::idNode> queue;
          queue.push(root);
          while(!queue.empty()){
            ftm::idNode node = queue.front();
            queue.pop();
            std::tuple<dataType, dataType> newBirthDeath;
            if(node < indexAddedNodes){
//std::cout << "not added" << std::endl;
              newBirthDeath = interpolation<dataType>(baryMergeTree, node, newScalarsVector, trees, 
                                                      baryMatching[node], alphas);
            }else{
//std::cout << "added" << std::endl;
              int i = nodesAddedTree[node];
              ftm::idNode nodeT = baryMatching[node][i];
              newBirthDeath = interpolationAdded<dataType>(trees[i], nodeT, alphas[i], baryMergeTree, node, 
                                                           newScalarsVector);
            }
            newScalarsVector[node] = std::get<1>(newBirthDeath);
            newScalarsVector[baryTree->getNode(node)->getOrigin()] = std::get<0>(newBirthDeath);
//////////////////////////////////////////
//if(newScalarsVector[baryTree->getNode(node)->getOrigin()] == newScalarsVector[node] and !isLeaf(baryTree, node)){
  //cptBug++;
  //std::cout << "bug" << std::endl;
  //std::cout << node << " " << baryTree->getNode(node)->getOrigin() << " (" << newScalarsVector[node] << ")" << std::endl;
  /*int cptNotMatched=0;
  for(int i = 0; i < baryMatching[node].size(); ++i)
    //std::cout << baryMatching[node][i] << std::endl;
    if(baryMatching[node][i] == -1)
      ++cptNotMatched;
  std::cout << cptNotMatched << " / " << baryMatching[node].size() << std::endl;*/
  //myPause();
//}
//////////////////////////////////////////
            auto children = getChildren(baryTree, node);
            for(auto child : children)
              queue.push(child);
          }
//std::cout << cptBug << " _/_ " << getRealNumberOfNodes(baryTree) << std::endl;
          
          setTreeScalars(baryMergeTree, newScalarsVector);
          if(NormalizedWasserstein and RescaledWasserstein)
            purgeBarycenter<dataType>(baryMergeTree, baryMatching, trees, alphas);
          persistenceThresholding<dataType>(baryMergeTree->tree, 0);
          //cleanMergeTree<dataType>(baryMergeTree);
        }
        
        int getNumberOfRoots(ftm::FTMTree_MT *tree){
          int noRoots = 0;
          for(int i = 0; i < tree->getNumberOfNodes(); ++i)
            noRoots += (isRoot(tree, i) and not isLeaf(tree, i)) ? 1 : 0;
          return noRoots;
        }

        template <class dataType>
        void updateBarycenterTree(std::vector<ftm::FTMTree_MT *> &trees, MergeTree *&baryMergeTree, 
                                  std::vector<double> &alphas, 
                                  std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>> &matchings){
          int indexAddedNodes = baryMergeTree->tree->getNumberOfNodes();
          updateBarycenterTreeStructure<dataType>(trees, baryMergeTree, matchings);
          updateBarycenterTreeScalars<dataType>(trees, baryMergeTree, alphas, indexAddedNodes, matchings);
        }
        
        template <class dataType>
        void computeOneDistance(ftm::FTMTree_MT *tree, ftm::FTMTree_MT *baryTree, 
                                std::vector<std::tuple<ftm::idNode, ftm::idNode>> &matching,
                                dataType &distance, int verboseT=0){
          FTMTreeEditDistance ftmTreeEditDistance;
          ftmTreeEditDistance.setVerbose(0);
          ftmTreeEditDistance.setProgressiveComputation(false);
          ftmTreeEditDistance.setPreprocess(false);
          ftmTreeEditDistance.setPostprocess(false);
          ftmTreeEditDistance.setBranchDecomposition(true);
          ftmTreeEditDistance.setNormalizedWasserstein(NormalizedWasserstein);
          ftmTreeEditDistance.setNormalizedWassersteinReg(NormalizedWassersteinReg);
          ftmTreeEditDistance.setRescaledWasserstein(RescaledWasserstein);
          ftmTreeEditDistance.setKeepSubtree(KeepSubtree);
          ftmTreeEditDistance.setAssignmentSolver(AssignmentSolver);
          //ftmTreeEditDistance.setParallelize(true);
          ftmTreeEditDistance.setParallelize(Parallelize);
          ftmTreeEditDistance.setNumberOfThreads(NumberOfThreads);
          ftmTreeEditDistance.setIsCalled(true);
          ftmTreeEditDistance.setDistanceSquared(true); // squared root
          distance = ftmTreeEditDistance.execute<dataType>(baryTree, tree, matching);
          if(verboseT >= 2){
            std::cout << "distance tree : " << distance << std::endl;
            std::cout << "distance²tree : " << distance*distance << std::endl;
            //printOutputMatching<dataType>(matching, baryTree, tree);
            //printMatching(matching);
          }
        }
        
        template <class dataType>
        void computeOneDistance(ftm::FTMTree_MT *tree, MergeTree *baryMergeTree, 
                                std::vector<std::tuple<ftm::idNode, ftm::idNode>> &matching,
                                dataType &distance, int verboseT=0){
          computeOneDistance<dataType>(tree, baryMergeTree->tree, matching, distance, verboseT);
        }
        
        template <class dataType>
        void computeOneDistance(MergeTree *baryMergeTree, MergeTree *baryMergeTree2, 
                                std::vector<std::tuple<ftm::idNode, ftm::idNode>> &matching,
                                dataType &distance, int verboseT=0){
          computeOneDistance<dataType>(baryMergeTree->tree, baryMergeTree2, matching, distance, verboseT);
        }
        
        template <class dataType>
        void assignment(std::vector<ftm::FTMTree_MT *> &trees, MergeTree *baryMergeTree,
                        std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>> &matchings,
                        std::vector<dataType> &distances, int verboseT=0){
#ifdef TTK_ENABLE_OPENMP
//#pragma omp parallel for schedule(dynamic) num_threads(NumberOfThreads) if(Parallelize)
#pragma omp parallel num_threads(NumberOfThreads) if(Parallelize and not IsCalled)
{
#pragma omp single nowait
{
//std::cout << "level 1 " << omp_get_num_threads() << std::endl;
#endif
          for(int i = 0; i < trees.size(); ++i)
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(i) untied
#endif
            computeOneDistance<dataType>(trees[i], baryMergeTree, matchings[i], distances[i], verboseT);
#ifdef TTK_ENABLE_OPENMP
} // pragma omp single nowait
#pragma omp taskwait
} // pragma omp parallel
#endif
        }
        
        // ----------------------------------------
        // Main Functions
        // ----------------------------------------
        template <class dataType>
        MergeTree* computeBarycenter(std::vector<ftm::FTMTree_MT *> &trees, MergeTree *&baryMergeTree,
                            std::vector<double> &alphas,
                            std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>> &finalMatchings,
                            int verboseT=0){
          Timer t_bary;

          ftm::FTMTree_MT *baryTree = baryMergeTree->tree;
          auto noNodesT = baryTree->getNumberOfNodes();
          auto noNodes = getRealNumberOfNodes(baryTree);
          if(verboseT >= 1)
            std::cout << "Barycenter number of nodes : " << noNodes << " / " << noNodesT << std::endl;
          
          /*baryTree->printTree2();
          printTreeScalars<dataType>(baryTree); 
          printPairsFromTree<dataType>(baryTree, BranchDecomposition); */
          
          bool converged = false;
          dataType frechetEnergy = -1;
          dataType minFrechet = std::numeric_limits<dataType>::max();
          int cptBlocked = 0;
          while(not converged){
            if(verboseT >= 1)
              std::cout << "------------------" << std::endl;
            // --- Assignment
            std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>> matchings(trees.size());
            std::vector<dataType> distances(trees.size(), -1);
            Timer t_assignment;
            assignment<dataType>(trees, baryMergeTree, matchings, distances);
            auto t_assignment_time = t_assignment.getElapsedTime();
            if(verboseT >= 1)
              std::cout << "assignment : " << t_assignment_time << std::endl;
            
            // --- Update
            Timer t_update;
            updateBarycenterTree<dataType>(trees, baryMergeTree, alphas, matchings);
            auto t_update_time = t_update.getElapsedTime();
            if(verboseT >= 1)
              std::cout << "update     : " << t_update_time << std::endl;
            baryTree = baryMergeTree->tree;
            auto noNodesT = baryTree->getNumberOfNodes();
            auto noNodes = getRealNumberOfNodes(baryTree);
            if(verboseT >= 1)
              std::cout << "Barycenter number of nodes : " << noNodes << " / " << noNodesT << std::endl;
            
            // --- Check convergence
            dataType currentFrechetEnergy = 0;
            for(int i = 0; i < trees.size(); ++i)
              currentFrechetEnergy += alphas[i] * distances[i]*distances[i];
            converged = (myAbs<dataType>(frechetEnergy - currentFrechetEnergy) <= Tol);
            frechetEnergy = currentFrechetEnergy;
            if(verboseT >= 1)
              std::cout << "Frechet energy : " << frechetEnergy << std::endl;
            //if(Tol == 0.0)
              Tol = frechetEnergy / 10000.0;
            
            minFrechet = std::min(minFrechet, frechetEnergy);
            if(not converged){
              cptBlocked += (minFrechet < frechetEnergy) ? 1 : 0;
              converged = (cptBlocked >= 10);
            }
            
            /*baryTree->printTree2();
            printTreeScalars<dataType>(baryTree); 
            printPairsFromTree<dataType>(baryTree, BranchDecomposition);*/
          }
          
          /*baryTree->printTree2();
          printTreeScalars<dataType>(baryTree); 
          printPairsFromTree<dataType>(baryTree, BranchDecomposition);*/
          
          std::vector<dataType> distances(trees.size(), -1);
          assignment<dataType>(trees, baryMergeTree, finalMatchings, distances, verboseT);
          dataType currentFrechetEnergy = 0;
          for(int i = 0; i < trees.size(); ++i)
            currentFrechetEnergy += alphas[i] * distances[i]*distances[i];
          if(verboseT >= 1)
            std::cout << "Frechet energy : " << currentFrechetEnergy << std::endl;
          
          if(verboseT >= 1)
            std::cout << "TIME BARYCENTER = " << t_bary.getElapsedTime() << std::endl;

          if(trees.size() == 2)
            verifyBarycenterTwoTrees<dataType>(trees, baryMergeTree, finalMatchings, distances);

/*std::cout << "matching :" << std::endl;
int avgMatch = 0, avgNoNodes=0;
for(int i = 0; i < finalMatchings.size(); ++i){
  //std::cout << finalMatchings[i].size() << " / " << getRealNumberOfNodes(trees[i]) << std::endl;
  avgMatch += finalMatchings[i].size();
  avgNoNodes += getRealNumberOfNodes(trees[i]);
}
avgMatch /= finalMatchings.size();
avgNoNodes /= finalMatchings.size();
std::cout << avgMatch << " / " << avgNoNodes << " (" << finalMatchings.size() << ")" << std::endl;*/
          
          return baryMergeTree;
        }
        
        template <class dataType>
        MergeTree* execute(std::vector<ftm::FTMTree_MT *> &trees, std::vector<double> &alphas,
                           std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>> &finalMatchings,
                           MergeTree *baryMergeTree=nullptr){
          // --- Preprocessing
          int avgNodes = 0, avgNodesT = 0;
          for(int i = 0; i < trees.size(); ++i){
            preprocessingPipeline<dataType>(trees[i], EpsilonTree2, Epsilon2Tree2, 0);
            auto noNodesT = trees[i]->getNumberOfNodes();
            auto noNodes = getRealNumberOfNodes(trees[i]);
            //std::cout << "tree " << i << " number of nodes : " << noNodes << " / " << noNodesT << std::endl;
            avgNodes += noNodes;
            avgNodesT += noNodesT;
            /*printMultiPersPairsFromTree<dataType>(trees[i], BranchDecomposition);
            printPairsFromTree<dataType>(trees[i], BranchDecomposition);*/
          }
          avgNodes /= trees.size();
          avgNodesT /= trees.size();
          std::cout << trees.size() << " trees average " << avgNodes << " / "<< avgNodesT << std::endl;
          
          if(baryMergeTree == nullptr)
            baryMergeTree = initBarycenterTree<dataType>(trees);
          
          // --- Execute
          baryMergeTree = computeBarycenter<dataType>(trees, baryMergeTree, alphas, finalMatchings, 2);
          
          // --- Postprocessing
          for(int i = 0; i < trees.size(); ++i)
            postprocessingPipeline<dataType>(trees[i]);
          postprocessingPipeline<dataType>(baryMergeTree->tree);
          for(int i = 0; i < trees.size(); ++i){
            convertBranchDecompositionMatching(baryMergeTree->tree, trees[i], finalMatchings[i]);
            //printOutputMatching<dataType>(finalMatchings[i], baryMergeTree->tree, trees[i], false);
          }
          
          return baryMergeTree;
        }
        
        template <class dataType>
        MergeTree* execute(std::vector<ftm::FTMTree_MT *> &trees, 
                           std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>> &finalMatchings,
                           MergeTree *baryMergeTree=nullptr){          
          //trees = makeMyTestBary<dataType>();
          
          std::vector<double> alphas;
          for(int i = 0; i < trees.size(); ++i)
            alphas.push_back(1.0/trees.size());
            
          return execute<dataType>(trees, alphas, finalMatchings, baryMergeTree);
        }
        
        // ----------------------------------------
        // Testing
        // ----------------------------------------
        template<class dataType>
        std::vector<ftm::FTMTree_MT *> makeMyTestBary(){
          std::vector<ftm::FTMTree_MT *> trees;
          
          /*float *nodesScalar1 = new float[4]{0, 1, 3, 8};
          std::vector<SimplexId> nodes1{0, 1, 2, 3};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs1{std::make_tuple(0, 2),
                                                                  std::make_tuple(1, 2),
                                                                  std::make_tuple(2, 3)};
          ftm::FTMTree_MT *treeTemp1 = makeFakeTree(nodesScalar1, nodes1, arcs1);
          trees.push_back(treeTemp1);
          
          float *nodesScalar2 = new float[2]{4, 12};
          std::vector<SimplexId> nodes2{0, 1};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs2{std::make_tuple(0, 1)};
          ftm::FTMTree_MT *treeTemp2 = makeFakeTree(nodesScalar2, nodes2, arcs2);
          trees.push_back(treeTemp2);
          
          float *nodesScalar3 = new float[2]{3, 11};
          std::vector<SimplexId> nodes3{0, 1};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs3{std::make_tuple(0, 1)};
          ftm::FTMTree_MT *treeTemp3 = makeFakeTree(nodesScalar3, nodes3, arcs3);
          trees.push_back(treeTemp3);*/
          float *nodesScalar1 = new float[2]{3, 12};
          std::vector<SimplexId> nodes1{0, 1};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs1{std::make_tuple(0, 1)};
          ftm::FTMTree_MT *treeTemp1 = makeFakeTree(nodesScalar1, nodes1, arcs1);
          trees.push_back(treeTemp1);
          
          float *nodesScalar2 = new float[2]{3, 12};
          std::vector<SimplexId> nodes2{0, 1};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs2{std::make_tuple(0, 1)};
          ftm::FTMTree_MT *treeTemp2 = makeFakeTree(nodesScalar2, nodes2, arcs2);
          trees.push_back(treeTemp2);
          
          float *nodesScalar3 = new float[2]{5, 10};
          std::vector<SimplexId> nodes3{0, 1};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs3{std::make_tuple(0, 1)};
          ftm::FTMTree_MT *treeTemp3 = makeFakeTree(nodesScalar3, nodes3, arcs3);
          trees.push_back(treeTemp3);
          
          float *nodesScalar4 = new float[4]{5, 6, 8, 10};
          std::vector<SimplexId> nodes4{0, 1, 2, 3};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs4{std::make_tuple(0, 2),
                                                                  std::make_tuple(1, 2),
                                                                  std::make_tuple(2, 3)};
          ftm::FTMTree_MT *treeTemp4 = makeFakeTree(nodesScalar4, nodes4, arcs4);
          trees.push_back(treeTemp4);
          
          return trees;
        }
        
        template <class dataType>
        void verifyBarycenterTwoTrees(std::vector<ftm::FTMTree_MT *> &trees, MergeTree *&baryMergeTree,
                            std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>> &finalMatchings,
                            std::vector<dataType> distances){
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> matching;
          dataType distance;
          computeOneDistance(trees[0], trees[1], matching, distance);
          if(distance != (distances[0] + distances[1])){
            std::cout << "distance T1 T2    : " << distance << std::endl;
            std::cout << "distance T1 T' T2 : " << distances[0] + distances[1] << std::endl;
          }
          /*printMatching(matching);
          
          for(int i = 0; i < finalMatchings[0].size(); ++i)
            for(int j = 0; j < finalMatchings[1].size(); ++j)
              if(std::get<0>(finalMatchings[0][i]) == std::get<0>(finalMatchings[1][j]))
                std::cout << std::get<1>(finalMatchings[1][j]) << " - " << std::get<1>(finalMatchings[0][i]) << " (" << std::get<0>(finalMatchings[0][i]) << ")" << std::endl;*/
          return;
          
          auto baryTree = baryMergeTree->tree;
          std::vector<std::vector<ftm::idNode>> baryMatched(baryTree->getNumberOfNodes(), 
                                                             std::vector<ftm::idNode>(trees.size(), -1));
          for(int i = 0; i < finalMatchings.size(); ++i)
            for(auto match : finalMatchings[i])
              baryMatched[std::get<0>(match)][i] = std::get<1>(match);
          
          std::queue<ftm::idNode> queue;
          queue.push(getRoot(baryTree));
          while(!queue.empty()){
            auto node = queue.front();
            queue.pop();
            std::vector<dataType> costs(trees.size());
            for(int i = 0; i < trees.size(); ++i)
              if(baryMatched[node][i] != -1)
                costs[i] = relabelCost<dataType>(baryTree, node, trees[i], baryMatched[node][i]);
              else
                costs[i] = deleteCost<dataType>(baryTree, node);
            dataType cost = 0;
            if(baryMatched[node][0] != -1 and baryMatched[node][1] != -1)
              cost = relabelCost<dataType>(trees[0], baryMatched[node][0], trees[1], baryMatched[node][1]);
            else if(baryMatched[node][0] == -1)
              cost = deleteCost<dataType>(trees[1], baryMatched[node][1]);
            else if(baryMatched[node][1] == -1)
              cost = deleteCost<dataType>(trees[0], baryMatched[node][0]);
            else
              std::cout << "problem" << std::endl;
            costs[0] = std::sqrt(costs[0]);
            costs[1] = std::sqrt(costs[1]);
            cost = std::sqrt(cost);
            if( myAbs<dataType>(costs[0] - costs[1]) > 1e-7 ){
              std::cout << "===================" << std::endl;
              std::cout << "cost T' T0    : " << costs[0] << std::endl;
              std::cout << "cost T' T1    : " << costs[1] << std::endl;
              std::cout << "cost T0 T1    : " << cost << std::endl;
              std::cout << "cost T0 T' T1 : " << costs[0] + costs[1] << std::endl;
              if(myAbs<dataType>((costs[0] + costs[1]) - cost) > 1e-7)
                std::cout << "diff          : " << myAbs<dataType>((costs[0] + costs[1]) - cost)<<std::endl;
              std::cout << "diff2         : " << myAbs<dataType>(costs[0] - costs[1]) << std::endl;
              printNode<dataType>(baryTree, node);
              printNode<dataType>(baryTree, getParent(baryTree, node));
              for(int i = 0; i < 2; ++i)
                if(baryMatched[node][i] != -1){
                  printNode<dataType>(trees[i], baryMatched[node][i]);
                  printNode<dataType>(trees[i], getParent(trees[i], baryMatched[node][i]));
                }
            }
            auto children = getChildren(baryTree, node);
            for(auto child : children)
              queue.push(child);
          }
        }
        
    }; // FTMTreeEditDistanceBarycenter class
    
} // namespace ttk

#endif
