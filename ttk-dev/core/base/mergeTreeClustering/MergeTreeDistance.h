/// Provide your information
///
/// \ingroup base
/// \class ttk::MergeTreeClustering
/// \author XXXXX
/// \date 2020.
///
/// This module defines the %MergeTreeClustering class that computes the edit distance
/// between two merge trees.
///

#ifndef _MERGETREEDISTANCE_H
#define _MERGETREEDISTANCE_H

#pragma once

#include <stack>
#include <thread>

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>

#include <BottleneckDistance.h>
#include <Munkres.h>
#include "AssignmentExhaustive.h"
#include "AssignmentAuction.h"
#include "MergeTreeBase.h"

#include "MergeTreeDistanceProgressive.h"

namespace ttk {
  
    /**
    * The MergeTreeDistance class provides methods to compute distance
    * between two merge trees.
    */
    class MergeTreeDistance : virtual public Debug, public MergeTreeBase {
      
    private:
        AssignmentExhaustive assignmentSolverExhaustive;
        double t_assignment_time = 0;
        
        bool Preprocess = true;
        bool Postprocess = true;
        bool SaveTree = false;
        bool OnlyEmptyTreeDistance = false;
        
        bool IsCalled = false;
        
        double AuctionEpsilon = -1;
        double AuctionEpsilonDiviser = 0;
        int AuctionNoRound = -1;

        // Just to get some stats about run
        //std::map<int, int> assignmentProblemSize, assignmentProblemIter;

        bool Testing = true; 
        
    public:
        MergeTreeDistance() {
          this->setDebugMsgPrefix(
              "MergeTreeDistance"); // inherited from Debug: prefix will be printed at the
          // beginning of every msg
          omp_set_nested(1);
        };
        ~MergeTreeDistance(){};
        
        void setIsCalled(bool ic){
          IsCalled = ic;
        }

        void setPreprocess(bool preproc){
          if(not ProgressiveComputation)
            Preprocess = preproc;
        }

        void setPostprocess(bool postproc){
          Postprocess = postproc;
        }
        
        void setTesting(bool test){
          Testing = test;
        }
        
        void setSaveTree(bool save){
          SaveTree = save;
        }
        
        void setAuctionEpsilon(double aucEpsilon){
          AuctionEpsilon = aucEpsilon;
        }

        void setAuctionEpsilonDiviser(double aucEpsilonDiviser){
          AuctionEpsilonDiviser = aucEpsilonDiviser;
        }
        
        void setAuctionNoRounds(double aucNoRounds){
          AuctionNoRound = aucNoRounds;
        }
        
        void setOnlyEmptyTreeDistance(double only){
          OnlyEmptyTreeDistance = only;
        }

        /**
        * Implementation of the algorithm.
        */
        
        // ----------------------------------------
        // Assignment Problem
        // ----------------------------------------
        template<class dataType>
        void forestAssignmentProblemMunkres(std::vector<std::vector<dataType>> &costMatrix, 
                                            std::vector<matchingTuple> &matchings){
          int nRows = costMatrix.size();
          int nCols = costMatrix[0].size();
          Munkres assignmentSolverMunkres;
          assignmentSolverMunkres.setInput(nRows, nCols, (void *)&costMatrix);
          assignmentSolverMunkres.run<dataType>(matchings);
        }
        
        template<class dataType>
        void forestAssignmentProblemExhaustiveSearch(std::vector<std::vector<dataType>> &costMatrix, 
                                                     std::vector<matchingTuple> &matchings){
          int nRows = costMatrix.size();
          int nCols = costMatrix[0].size();
          assignmentSolverExhaustive.setInput(nRows, nCols, (void *)&costMatrix);
          assignmentSolverExhaustive.run<dataType>(matchings);
        }
        
        template<class dataType>
        void forestAssignmentProblemAuction(std::vector<std::vector<dataType>> &costMatrix, 
                                            std::vector<matchingTuple> &matchings){
          int nRows = costMatrix.size();
          int nCols = costMatrix[0].size();
          AssignmentAuction<dataType> assignmentSolverAuction;
          assignmentSolverAuction.setInput(nRows, nCols, (void *)&costMatrix);
          assignmentSolverAuction.setBalanced(false);
          assignmentSolverAuction.run(matchings);
          assignmentSolverAuction.setEpsilon(AuctionEpsilon);
          assignmentSolverAuction.setEpsilonDiviserMultiplier(AuctionEpsilonDiviser);
          assignmentSolverAuction.setNumberOfRounds(AuctionNoRound);
          
          //assignmentProblemIter[assignmentSolverAuction.getIter()]++;
        }
        
        template<class dataType>
        void runAssignmentProblemSolver(std::vector<std::vector<dataType>> &costMatrix,
                                        std::vector<matchingTuple> &matchings){
          //AssignmentSolver = 2;
          switch(AssignmentSolver){
            case 1:
              forestAssignmentProblemExhaustiveSearch<dataType>(costMatrix, matchings);
              break;
            case 2:
              forestAssignmentProblemMunkres<dataType>(costMatrix, matchings);
              break;
            case 0:
            default:
              forestAssignmentProblemAuction<dataType>(costMatrix, matchings);
          }
        }
        
        template<class dataType>
        dataType forestAssignmentProblem(std::vector<std::vector<dataType>> &treeTable,
                                         std::vector<ftm::idNode> &children1,
                                         std::vector<ftm::idNode> &children2,
                                         std::vector<std::tuple<int, int>> &forestAssignment){
          // --- Create cost matrix
          int nRows = children1.size(), nCols = children2.size();
          std::vector<std::vector<dataType>> costMatrix(nRows+1, std::vector<dataType>(nCols+1));
          createCostMatrix(treeTable, children1, children2, costMatrix);
          //printTableVector(costMatrix);
          //std::cout << costMatrix.size() << " _ " << costMatrix[0].size() << std::endl;
          
          //assignmentProblemSize[costMatrix.size()*costMatrix[0].size()]++;
          
          // --- Solve assignment problem
          std::vector<matchingTuple> matchings;
          runAssignmentProblemSolver(costMatrix, matchings);
          //printMatching(matchings);
          
          // --- Postprocess matching to create output assignment          
          dataType cost = postprocessAssignment<dataType>(matchings, children1, children2, forestAssignment);
          
          return cost;
        }
        
        template<class dataType>
        void computeEquation13(int i, int j,
                        std::vector<std::vector<dataType>> &treeTable,
                        std::vector<std::vector<dataType>> &forestTable,
                        std::vector<std::vector<std::vector<std::tuple<int, int>>>> &forestBackTable,
                        std::vector<ftm::idNode> &children1, std::vector<ftm::idNode> &children2){
          //std::cout << i << " _ " << j << std::endl;
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
            forestTerm3 = forestAssignmentProblem<dataType>(treeTable, children1, 
                                                            children2, forestAssignment);
            if(not Parallelize)
              t_assignment_time += t_assignment.getElapsedTime();
            
            // Compute table value
            forestTable[i][j] = KeepSubtree ? std::min(std::min(forestTerm1, forestTerm2), forestTerm3) :
                                              forestTerm3;
            
            // Add backtracking information
            if(forestTable[i][j] == forestTerm3){
              forestBackTable[i][j] = forestAssignment;
              //forestBackTable[i][j].insert(forestBackTable[i][j].begin(), forestAssignment.begin(), forestAssignment.end());
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
        
        // ----------------------------------------
        // Main Functions
        // ----------------------------------------
        template <class dataType>
        dataType execute(ftm::FTMTree_MT *&tree1, ftm::FTMTree_MT *&tree2,
                         std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &outputMatching){
            Memory m;
            Timer t_total;
            
            // ---------------------
            // ----- Testing
            // --------------------
            Testing = false;
            if(Testing){
              //ProgressiveComputation = true;
              //BranchDecomposition = true;
              //BarycenterMergeTree = true;
              BranchDecomposition = false;
              NormalizedWasserstein = false;
              KeepSubtree = true;
              //Parallelize = true;
              EpsilonTree1 = 0; EpsilonTree2 = 0; Epsilon2Tree1 = 0; Epsilon2Tree2 = 0;
              makeMyTest<dataType>(tree1, tree2);
              //printTree(tree1); printTree(tree2); 
              //printTreeScalars<dataType>(tree1); printTreeScalars<dataType>(tree2);
              //printPairsFromTree<dataType>(tree1); printPairsFromTree<dataType>(tree2);
              
              //MergeTreeClusteringBarycenter bary;
              //bary.execute<dataType>(makeMyTest2<dataType>());
            }
            
            // ---------------------
            // ----- Preprocessing
            // --------------------
            //std::cout << "Preprocessing" << std::endl;
            ftm::FTMTree_MT *tree1Ori = tree1;
            ftm::FTMTree_MT *tree2Ori = tree2;
            //printTree(tree2);
            if(SaveTree){
              tree1 = copyTree<dataType>(tree1);
              tree2 = copyTree<dataType>(tree2);
            }
            if(not IsCalled){
              verifyMergeTreeStructure<dataType>(tree1); 
              verifyMergeTreeStructure<dataType>(tree2);
            }
            ftm::FTMTree_MT *tree1Old;
            ftm::FTMTree_MT *tree2Old;
            //printTree(tree2); printTreeScalars<dataType>(tree2);
            if(Preprocess){
              treesNodeCorr = std::vector<std::vector<int>>(2);
              tree1Old = preprocessingPipeline<dataType>(tree1, EpsilonTree1, Epsilon2Tree1, Epsilon3Tree1,
                                                         BranchDecomposition, UseMinMaxPair, CleanTree, 
                                                         treesNodeCorr[0], verbose);
              tree2Old = preprocessingPipeline<dataType>(tree2, EpsilonTree2, Epsilon2Tree2, Epsilon3Tree2,
                                                         BranchDecomposition, UseMinMaxPair, CleanTree, 
                                                         treesNodeCorr[1], verbose);
            }

            //printTree(tree1); printPairsFromTree<dataType>(tree1, BranchDecomposition); 
            //printTree(tree2); printPairsFromTree<dataType>(tree2, BranchDecomposition); 
            //printTreeStats(tree1); printTreeStats(tree2);
            
            // ---------------------
            // ----- Init dynamic progamming tables
            // --------------------
            //std::cout << "Init dynamic progamming tables" << std::endl;
            size_t nRows = tree1->getNumberOfNodes() + 1;
            size_t nCols = tree2->getNumberOfNodes() + 1;
            std::vector<std::vector<dataType>> treeTable(nRows, std::vector<dataType>(nCols));
            std::vector<std::vector<dataType>> forestTable(nRows, std::vector<dataType>(nCols));

            // Backtracking tables (output matching)
            std::vector<std::vector<std::tuple<int, int>>> 
                            treeBackTable(nRows, std::vector<std::tuple<int, int>>(nCols));
            std::vector<std::vector<std::vector<std::tuple<int, int>>>> 
                            forestBackTable(nRows, std::vector<std::vector<std::tuple<int, int>>>(nCols));
            //std::cout << "init dynamic programming tables done." << std::endl;
                            
            int indR = getRoot(tree1)+1;
            int indC = getRoot(tree2)+1;
            
            /*tree1Level = std::vector<int>(tree1->getNumberOfNodes(), 0);
            tree2Level = std::vector<int>(tree2->getNumberOfNodes(), 0);*/
            tree1Level = getAllNodeLevel(tree1);
            tree2Level = getAllNodeLevel(tree2);
            
            // ---------------------
            // ----- Compute edit distance
            // --------------------
            //std::cout << "Compute edit distance" << std::endl;
            if(ProgressiveComputation){
              MergeTreeDistanceProgressive<dataType> editDistanceProgressive;
              editDistanceProgressive.computeEditDistanceProgressive(tree1, tree2, 
                                  treeTable, forestTable, treeBackTable, forestBackTable);
            }else{
              computeEditDistance(tree1, tree2, treeTable, forestTable, 
                                  treeBackTable, forestBackTable, nRows, nCols);
            }
            dataType distance = treeTable[indR][indC];
            if(OnlyEmptyTreeDistance)
              distance = treeTable[indR][0];
            if(DistanceSquared)
              distance = std::sqrt(distance);
            //std::cout << distance << std::endl; 
            //printTableVector<dataType>(treeTable); 
            
            // ---------------------
            // ----- Compute matching
            // --------------------
            //std::cout << "Compute matching" << std::endl;
            Timer t_match;
            computeMatching<dataType>(tree1, tree2, treeBackTable, forestBackTable, outputMatching, 
                                      indR, indC);
            auto t_match_time = t_match.getElapsedTime();
            //printOutputMatching<dataType>(outputMatching, tree1, tree2);
            //printMatching(outputMatching);
            
            // ---------------------
            // ----- Postprocessing
            // --------------------
            //std::cout << "Postprocessing" << std::endl;
            if(Postprocess){
              // TODO see fixMergedRootOrigin in MergeTreeClusteringBase
              fixMergedRootOrigin<dataType>(tree1); 
              fixMergedRootOrigin<dataType>(tree2);
              postprocessingPipeline<dataType>(tree1); 
              //printTree(tree1); printPairsFromTree<dataType>(tree1);
              postprocessingPipeline<dataType>(tree2); 
              //printTree(tree2); printPairsFromTree<dataType>(tree2);
              if(BranchDecomposition)
                convertBranchDecompositionMatching<dataType>(tree1, tree2, outputMatching);
              
              /*printMatching(outputMatching);
              printTree(tree1); printPairsFromTree<dataType>(tree1, BranchDecomposition);
              printTree(tree2); printPairsFromTree<dataType>(tree2, BranchDecomposition); */
            }
            
            //printMatching(outputMatching);
            //printOutputMatching<dataType>(outputMatching, tree1, tree2);
            //printOutputMatching<dataType>(outputMatching, tree1, tree2, false);
            //printTableVector<dataType>(treeTable); printTableVector<dataType>(forestTable);
            //printMapIntInt(assignmentProblemSize);
            //printMapIntInt(assignmentProblemIter);
            
            if(verbose > 0){
              std::cout << "TIME COMP.MATC. = " << t_match_time << std::endl;
              std::cout << "TIME TOTAL      = " << t_total.getElapsedTime() << std::endl;
              std::cout << "--------------------------------" << std::endl;
              std::cout << "DISTANCE²       = " << treeTable[indR][indC] << std::endl;
              std::cout << "DISTANCE        = " << std::sqrt(treeTable[indR][indC]) << std::endl;
              /*std::cout << "DISTANCE X2     = " << 2*treeTable[indR][indC] << std::endl;
              std::cout << "DISTANCE X2SQRT = " << std::sqrt(2*treeTable[indR][indC]) << std::endl;*/
              std::cout << "--------------------------------" << std::endl;
              std::cout << "MEMORY          = " << m.getElapsedUsage() << std::endl;
              std::cout << "--------------------------------" << std::endl;
            }
            
            if(SaveTree){
              freeTree(tree1);
              freeTree(tree2);
              tree1 = tree1Ori;
              tree2 = tree2Ori;
            }
            
            return distance;
        }
        
        template <class dataType>
        dataType execute(ftm::FTMTree_MT *&tree1, ftm::FTMTree_MT *&tree2,
                         std::vector<std::tuple<ftm::idNode, ftm::idNode>> &outputMatching){
          std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> realOutputMatching;
          dataType res = execute<dataType>(tree1, tree2, realOutputMatching);
          for(auto tup : realOutputMatching)
            outputMatching.push_back(std::make_tuple(std::get<0>(tup), std::get<1>(tup)));
          return res;
        }
        
        template<class dataType>
        void computeEditDistance(ftm::FTMTree_MT *tree1, ftm::FTMTree_MT *tree2,
                               std::vector<std::vector<dataType>> &treeTable,
                               std::vector<std::vector<dataType>> &forestTable,
                               std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
                               std::vector<std::vector<std::vector<std::tuple<int, int>>>> &forestBackTable,
                               int nRows, int nCols){
          Timer t_dyn;
          t_assignment_time = 0;
          
          if(Parallelize){
            parallelEditDistance(tree1, tree2, treeTable, forestTable, 
                                 treeBackTable, forestBackTable, nRows, nCols);
          }else{
            /*classicEditDistance_v1(tree1, tree2, treeTable, forestTable, 
                                treeBackTable, forestBackTable, nRows, nCols);*/
            // Distance T1 to empty tree
            classicEditDistance(tree1, tree2, true, true, getRoot(tree1), getRoot(tree2), 
                                treeTable, forestTable, treeBackTable, forestBackTable, nRows, nCols);
            if(OnlyEmptyTreeDistance)
              return;
            // Distance T2 to empty tree
            classicEditDistance(tree1, tree2, false, true, getRoot(tree1), getRoot(tree2), 
                                treeTable, forestTable, treeBackTable, forestBackTable, nRows, nCols);
            // Distance T1 to T2
            classicEditDistance(tree1, tree2, true, false, getRoot(tree1), getRoot(tree2), 
                                treeTable, forestTable, treeBackTable, forestBackTable, nRows, nCols);
          }
          
          if(verbose > 0){
            std::cout << "TIME DYNA.PROG. = " << t_dyn.getElapsedTime() << std::endl;
            if(not Parallelize)
              std::cout << " - TIME ASSGNMT = " << t_assignment_time << std::endl;          
          }
        }
        
        /*template<class dataType>
        void classicEditDistance_v1(ftm::FTMTree_MT *tree1, ftm::FTMTree_MT *tree2, 
                              std::vector<std::vector<dataType>> &treeTable,
                              std::vector<std::vector<dataType>> &forestTable,
                              std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
                              std::vector<std::vector<std::vector<std::tuple<int, int>>>> &forestBackTable,
                              int nRows, int nCols){   
          for(int i = 1; i < nRows; ++i){
            ftm::idNode nodeI = i-1;
            if(isNodeAlone(tree1, nodeI)) continue;
            // --- Equation 8
            computeEquation8(tree1, nodeI, i, treeTable, forestTable);
            
            // --- Equation 9
            computeEquation9(tree1, nodeI, i, treeTable, forestTable);
          }
          
          for(int j = 1; j < nCols; ++j){
            ftm::idNode nodeJ = j-1;
            if(isNodeAlone(tree2, nodeJ)) continue;
            // --- Equation 10
            computeEquation10(tree2, nodeJ, j, treeTable, forestTable);
              
            // --- Equation 11
            computeEquation11(tree2, nodeJ, j, treeTable, forestTable);
          }
          
          for(int i = 1; i < nRows; ++i){
            ftm::idNode nodeI = i-1;
            if(isNodeAlone(tree1, nodeI)) continue;
            std::vector<ftm::idNode> children1 = getChildren(tree1, nodeI);
            for(int j = 1; j < nCols; ++j){
              ftm::idNode nodeJ = j-1;
              if(isNodeAlone(tree2, nodeJ)) continue;
              std::vector<ftm::idNode> children2 = getChildren(tree2, nodeJ);                
              // --- Equation 13
              computeEquation13(i, j, treeTable, forestTable, forestBackTable, children1, children2);
              
              // --- Equation 12
              computeEquation12(tree1, tree2, i, j, nodeI, nodeJ, treeTable, forestTable, treeBackTable,
                                children1, children2);
            }
          }
        }*/
        
        template<class dataType>
        void classicEditDistance(ftm::FTMTree_MT *tree1, ftm::FTMTree_MT *tree2, bool processTree1,
                              bool computeEmptyTree, ftm::idNode nodeI, ftm::idNode nodeJ,
                              std::vector<std::vector<dataType>> &treeTable,
                              std::vector<std::vector<dataType>> &forestTable,
                              std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
                              std::vector<std::vector<std::vector<std::tuple<int, int>>>> &forestBackTable,
                              int nRows, int nCols){       
          if(processTree1){
            std::vector<ftm::idNode> childrens = getChildren(tree1, nodeI);
            for(auto children : childrens)
              classicEditDistance(tree1, tree2, processTree1, computeEmptyTree, children, nodeJ, 
                                  treeTable, forestTable, treeBackTable, forestBackTable, nRows, nCols);
          }else{
            std::vector<ftm::idNode> childrens = getChildren(tree2, nodeJ);
            for(auto children : childrens)
              classicEditDistance(tree1, tree2, processTree1, computeEmptyTree, nodeI, children, 
                                  treeTable, forestTable, treeBackTable, forestBackTable, nRows, nCols);
          }
          
          if(processTree1){
            if(computeEmptyTree){
              int i = nodeI+1;
              // --- Equation 8
              computeEquation8(tree1, nodeI, i, treeTable, forestTable);
              
              // --- Equation 9
              computeEquation9(tree1, nodeI, i, treeTable, forestTable);
            }else
              classicEditDistance(tree1, tree2, false, false, nodeI, getRoot(tree2), treeTable, forestTable, 
                                  treeBackTable, forestBackTable, nRows, nCols);  
          }else{
            int j = nodeJ+1;
            if(computeEmptyTree){
              // --- Equation 10
              computeEquation10(tree2, nodeJ, j, treeTable, forestTable);
                
              // --- Equation 11
              computeEquation11(tree2, nodeJ, j, treeTable, forestTable);
            //}else{
            }else if(KeepSubtree or tree1Level[nodeI] == tree2Level[nodeJ]){
              int i = nodeI+1;
              std::vector<ftm::idNode> children1 = getChildren(tree1, nodeI);
              std::vector<ftm::idNode> children2 = getChildren(tree2, nodeJ);  
              // --- Equation 13
              computeEquation13(i, j, treeTable, forestTable, forestBackTable, children1, children2);
              
              // --- Equation 12
              computeEquation12(tree1, tree2, i, j, nodeI, nodeJ, treeTable, forestTable, treeBackTable,
                                children1, children2);
            }
          }
        }
                
        // ----------------------------------------
        // Parallel version
        // ----------------------------------------  
        template<class dataType>
        void parallelEditDistance(ftm::FTMTree_MT *tree1, ftm::FTMTree_MT *tree2, 
                              std::vector<std::vector<dataType>> &treeTable,
                              std::vector<std::vector<dataType>> &forestTable,
                              std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
                              std::vector<std::vector<std::vector<std::tuple<int, int>>>> &forestBackTable,
                              int nRows, int nCols){
          std::vector<int> tree1NodeChildSize, tree2NodeChildSize;
          for(int i = 0; i < tree1->getNumberOfNodes(); ++i)
            tree1NodeChildSize.push_back(getChildren(tree1, i).size());
          for(int j = 0; j < tree2->getNumberOfNodes(); ++j)
            tree2NodeChildSize.push_back(getChildren(tree2, j).size());
          
          std::vector<ftm::idNode> tree1Leaves = getLeaves(tree1);
          std::vector<ftm::idNode> tree2Leaves = getLeaves(tree2);

          // Distance T1 to empty tree
          parallelEmptyTreeDistance_v2(tree1, true, tree1Leaves, tree1NodeChildSize,
                                    treeTable, forestTable, treeBackTable, forestBackTable);
          if(OnlyEmptyTreeDistance)
            return;
          // Distance T2 to empty tree
          parallelEmptyTreeDistance_v2(tree2, false, tree2Leaves, tree2NodeChildSize,
                                    treeTable, forestTable, treeBackTable, forestBackTable);
          // Distance T1 to T2
          /*parallelDependTreeDistance(tree1, tree2, true, 0,
                               tree1Leaves, tree1NodeChildSize, tree2Leaves, tree2NodeChildSize,
                               treeTable, forestTable, treeBackTable, forestBackTable);*/
          parallelTreeDistance_v2(tree1, tree2, true, 0,
                               tree1Leaves, tree1NodeChildSize, tree2Leaves, tree2NodeChildSize,
                               treeTable, forestTable, treeBackTable, forestBackTable, true);
        }
      
        // Equation 12, 13
        template<class dataType>
        void parallelTreeDistance_v2(ftm::FTMTree_MT *tree1, ftm::FTMTree_MT *tree2, bool isTree1, int i,
                              std::vector<ftm::idNode> &tree1Leaves, std::vector<int> &tree1NodeChildSize,
                              std::vector<ftm::idNode> &tree2Leaves, std::vector<int> &tree2NodeChildSize,
                              std::vector<std::vector<dataType>> &treeTable,
                              std::vector<std::vector<dataType>> &forestTable,
                              std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
                              std::vector<std::vector<std::vector<std::tuple<int, int>>>> &forestBackTable,
                              bool firstCall=false){
          ftm::idNode nodeT = -1;
          ftm::FTMTree_MT *treeT = (isTree1) ? tree1 : tree2;
          std::vector<int> treeChildDone(treeT->getNumberOfNodes(), 0);
          std::vector<bool> treeNodeDone(treeT->getNumberOfNodes(), false);
          std::queue<ftm::idNode> treeQueue;
          if(isTree1)
            for(ftm::idNode leaf : tree1Leaves)
              treeQueue.push(leaf);
          else
            for(ftm::idNode leaf : tree2Leaves)
              treeQueue.push(leaf);
          if(not IsCalled) //and firstCall)
            parallelTreeDistancePara(tree1, tree2, isTree1, i, tree1Leaves, tree1NodeChildSize, 
                                     tree2Leaves, tree2NodeChildSize, treeTable, forestTable, 
                                     treeBackTable, forestBackTable, firstCall,
                                     nodeT, treeChildDone, treeNodeDone, treeQueue);
          else
            parallelTreeDistanceTask(tree1, tree2, isTree1, i, tree1Leaves, tree1NodeChildSize, 
                                     tree2Leaves, tree2NodeChildSize, treeTable, forestTable, 
                                     treeBackTable, forestBackTable, 
                                     nodeT, treeChildDone, treeNodeDone, treeQueue);
        }
        
        template<class dataType>
        void parallelTreeDistancePara(ftm::FTMTree_MT *tree1, ftm::FTMTree_MT *tree2, bool isTree1, int i,
                              std::vector<ftm::idNode> &tree1Leaves, std::vector<int> &tree1NodeChildSize,
                              std::vector<ftm::idNode> &tree2Leaves, std::vector<int> &tree2NodeChildSize,
                              std::vector<std::vector<dataType>> &treeTable,
                              std::vector<std::vector<dataType>> &forestTable,
                              std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
                              std::vector<std::vector<std::vector<std::tuple<int, int>>>> &forestBackTable,
                              bool firstCall, ftm::idNode nodeT, std::vector<int> &treeChildDone, 
                              std::vector<bool> &treeNodeDone, std::queue<ftm::idNode> &treeQueue){
#ifdef TTK_ENABLE_OPENMP
unsigned int nthreads = std::thread::hardware_concurrency();
#pragma omp parallel num_threads(NumberOfThreads) if(firstCall or nthreads == NumberOfThreads)
{
#pragma omp single nowait 
#endif
          parallelTreeDistanceTask(tree1, tree2, isTree1, i, tree1Leaves, tree1NodeChildSize, 
                                   tree2Leaves, tree2NodeChildSize, treeTable, forestTable, 
                                   treeBackTable, forestBackTable, 
                                   nodeT, treeChildDone, treeNodeDone, treeQueue);
#ifdef TTK_ENABLE_OPENMP
} // pragma omp parallel
#endif
        }
        
        template<class dataType>
        void parallelTreeDistanceTask(ftm::FTMTree_MT *tree1, ftm::FTMTree_MT *tree2, bool isTree1, int i,
                              std::vector<ftm::idNode> &tree1Leaves, std::vector<int> &tree1NodeChildSize,
                              std::vector<ftm::idNode> &tree2Leaves, std::vector<int> &tree2NodeChildSize,
                              std::vector<std::vector<dataType>> &treeTable,
                              std::vector<std::vector<dataType>> &forestTable,
                              std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
                              std::vector<std::vector<std::vector<std::tuple<int, int>>>> &forestBackTable,
                              ftm::idNode nodeT, std::vector<int> &treeChildDone, 
                              std::vector<bool> &treeNodeDone, std::queue<ftm::idNode> &treeQueue){
          while(!treeQueue.empty()){
            nodeT = treeQueue.front();
            treeQueue.pop();
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(nodeT) untied shared(treeTable, forestTable, treeBackTable, forestBackTable, treeChildDone, treeNodeDone)
{
#endif
            ftm::FTMTree_MT *treeT = (isTree1) ? tree1 : tree2;
            while(nodeT != -1){
              int t = nodeT+1;
              ftm::idNode nodeI = i-1;
              
              if(isTree1){
                parallelTreeDistance_v2(tree1, tree2, false, t,
                                     tree1Leaves, tree1NodeChildSize, tree2Leaves, tree2NodeChildSize, 
                                     treeTable, forestTable, treeBackTable, forestBackTable, false);
              //}else{
              }else if(KeepSubtree or tree1Level[nodeI] == tree2Level[nodeT]){
                int j = nodeT+1;
                std::vector<ftm::idNode> children1 = getChildren(tree1, nodeI);
                std::vector<ftm::idNode> children2 = getChildren(tree2, nodeT);  
                // --- Equation 13
                computeEquation13(i, j, treeTable, forestTable, forestBackTable, children1, children2);
                
                // --- Equation 12
                computeEquation12(tree1, tree2, i, j, nodeI, nodeT, treeTable, forestTable, treeBackTable,
                                  children1, children2);
              }
              
              // Manage parent
              ftm::idNode nodeTParent = getParent(treeT, nodeT);
              int childSize = (isTree1) ? tree1NodeChildSize[nodeTParent] : tree2NodeChildSize[nodeTParent];
              int oldTreeChildDone;
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
{
#endif
              oldTreeChildDone = treeChildDone[nodeTParent];
              treeChildDone[nodeTParent]++;
#ifdef TTK_ENABLE_OPENMP
} // pragma omp atomic capture
#endif
              if(not treeNodeDone[nodeTParent] and oldTreeChildDone+1 == childSize){
                nodeT = nodeTParent;
                treeNodeDone[nodeTParent] = true;
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskyield
#endif
                //std::cout << "notify " << nodeIParent << std::endl;
              }else
                nodeT = -1;
              
            } // while nodeI loop
#ifdef TTK_ENABLE_OPENMP
} // pragma omp task
#endif
          } // while treeQueue loop
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
        }
        
        // Equation 8, 9, 10, 11
        template<class dataType>
        void parallelEmptyTreeDistance_v2(ftm::FTMTree_MT *tree, bool isTree1, 
                              std::vector<ftm::idNode> &treeLeaves, std::vector<int> &treeNodeChildSize,
                              std::vector<std::vector<dataType>> &treeTable,
                              std::vector<std::vector<dataType>> &forestTable,
                              std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
                              std::vector<std::vector<std::vector<std::tuple<int, int>>>> &forestBackTable){          
          ftm::idNode nodeT = -1;
          std::vector<int> treeChildDone(tree->getNumberOfNodes(), 0);
          std::vector<bool> treeNodeDone(tree->getNumberOfNodes(), false);
          std::queue<ftm::idNode> treeQueue;
          for(ftm::idNode leaf : treeLeaves)
            treeQueue.push(leaf);
          if(not IsCalled)
            parallelEmptyTreeDistancePara(tree, isTree1, treeLeaves, treeNodeChildSize, treeTable, 
                                          forestTable, treeBackTable, forestBackTable,
                                          nodeT, treeChildDone, treeNodeDone, treeQueue);
          else
            parallelEmptyTreeDistanceTask(tree, isTree1, treeLeaves, treeNodeChildSize, treeTable, 
                                          forestTable, treeBackTable, forestBackTable,
                                          nodeT, treeChildDone, treeNodeDone, treeQueue);
        }
        
        template<class dataType>
        void parallelEmptyTreeDistancePara(ftm::FTMTree_MT *tree, bool isTree1, 
                              std::vector<ftm::idNode> &treeLeaves, std::vector<int> &treeNodeChildSize,
                              std::vector<std::vector<dataType>> &treeTable,
                              std::vector<std::vector<dataType>> &forestTable,
                              std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
                              std::vector<std::vector<std::vector<std::tuple<int, int>>>> &forestBackTable,
                              ftm::idNode nodeT, std::vector<int> &treeChildDone, 
                              std::vector<bool> &treeNodeDone, std::queue<ftm::idNode> &treeQueue){
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(NumberOfThreads)
{
#pragma omp single nowait
#endif
          parallelEmptyTreeDistanceTask(tree, isTree1, treeLeaves, treeNodeChildSize, treeTable, 
                                        forestTable, treeBackTable, forestBackTable,
                                        nodeT, treeChildDone, treeNodeDone, treeQueue);
#ifdef TTK_ENABLE_OPENMP
} // pragma omp parallel
#endif
        }
        
        template<class dataType>
        void parallelEmptyTreeDistanceTask(ftm::FTMTree_MT *tree, bool isTree1, 
                              std::vector<ftm::idNode> &treeLeaves, std::vector<int> &treeNodeChildSize,
                              std::vector<std::vector<dataType>> &treeTable,
                              std::vector<std::vector<dataType>> &forestTable,
                              std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
                              std::vector<std::vector<std::vector<std::tuple<int, int>>>> &forestBackTable,
                              ftm::idNode nodeT, std::vector<int> &treeChildDone, 
                              std::vector<bool> &treeNodeDone, std::queue<ftm::idNode> &treeQueue){
          while(!treeQueue.empty()){
            nodeT = treeQueue.front();
            treeQueue.pop();
            
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(nodeT) untied shared(treeTable, forestTable, treeBackTable, forestBackTable, treeChildDone, treeNodeDone)
{
#endif
            while(nodeT != -1){              
              if(isTree1){
                int i = nodeT+1;
                // --- Equation 8
                computeEquation8(tree, nodeT, i, treeTable, forestTable);
                
                // --- Equation 9
                computeEquation9(tree, nodeT, i, treeTable, forestTable);
              }else{
                int j = nodeT+1;
                // --- Equation 10
                computeEquation10(tree, nodeT, j, treeTable, forestTable);
                  
                // --- Equation 11
                computeEquation11(tree, nodeT, j, treeTable, forestTable);
              }
              
              // Manage parent
              ftm::idNode nodeTParent = getParent(tree, nodeT);
              int oldTreeChildDone;
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
{
#endif
              oldTreeChildDone = treeChildDone[nodeTParent];
              treeChildDone[nodeTParent]++;
#ifdef TTK_ENABLE_OPENMP
} // pragma omp atomic capture
#endif
              if(not treeNodeDone[nodeTParent] and oldTreeChildDone+1 == treeNodeChildSize[nodeTParent]){
                nodeT = nodeTParent;
                treeNodeDone[nodeTParent] = true;
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskyield
#endif
                //std::cout << "notify " << nodeIParent << std::endl;
              }else
                nodeT = -1;

            } // while nodeI loop
#ifdef TTK_ENABLE_OPENMP
} // pragma omp task
#endif
          } // while treeQueue loop
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
        }
        
        // ----- OLD v1
        // Equation 12, 13
        template<class dataType>
        void parallelTreeDistance(ftm::FTMTree_MT *tree1, ftm::FTMTree_MT *tree2, bool isTree1, int i,
                              std::vector<ftm::idNode> &tree1Leaves, std::vector<int> &tree1NodeChildSize,
                              std::vector<ftm::idNode> &tree2Leaves, std::vector<int> &tree2NodeChildSize,
                              std::vector<std::vector<dataType>> &treeTable,
                              std::vector<std::vector<dataType>> &forestTable,
                              std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
                              std::vector<std::vector<std::vector<std::tuple<int, int>>>> &forestBackTable,
                              bool firstCall=false){
          ftm::idNode nodeT = -1;
          ftm::FTMTree_MT *treeT = (isTree1) ? tree1 : tree2; 
          std::vector<int> treeChildDone(treeT->getNumberOfNodes(), 0);
          std::vector<bool> treeNodeDone(treeT->getNumberOfNodes(), false);
          std::queue<ftm::idNode> treeQueue;
          if(isTree1)
            for(ftm::idNode leaf : tree1Leaves)
              treeQueue.push(leaf);
          else
            for(ftm::idNode leaf : tree2Leaves)
              treeQueue.push(leaf);

#ifdef TTK_ENABLE_OPENMP
unsigned int nthreads = std::thread::hardware_concurrency();
#pragma omp parallel num_threads(NumberOfThreads) if((firstCall or nthreads == NumberOfThreads) and not IsCalled)
{
#pragma omp single nowait
#endif
          while(!treeQueue.empty()){
            nodeT = treeQueue.front();
            treeQueue.pop();
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(nodeT) untied //shared(treeTable, forestTable)//, treeBackTable, forestBackTable)
{
#endif
std::stringstream ss;
ss << &treeTable[0] << " _ " << tree1 << " _ " << tree2 << std::endl;
//std::cout << ss.str();
            while(nodeT != -1){
              int t = nodeT+1;
              
              if(isTree1){
                parallelTreeDistance(tree1, tree2, false, t,
                                     tree1Leaves, tree1NodeChildSize, tree2Leaves, tree2NodeChildSize, 
                                     treeTable, forestTable, treeBackTable, forestBackTable);
              }else{
                int j = nodeT+1;
                ftm::idNode nodeI = i-1;
                std::vector<ftm::idNode> children1 = getChildren(tree1, nodeI);
                std::vector<ftm::idNode> children2 = getChildren(tree2, nodeT);  
                // --- Equation 13
                computeEquation13(i, j, treeTable, forestTable, forestBackTable, children1, children2);
                
                // --- Equation 12
                computeEquation12(tree1, tree2, i, j, nodeI, nodeT, treeTable, forestTable, treeBackTable,
                                  children1, children2);
              }
              
              // Manage parent
              ftm::idNode nodeTParent = getParent(treeT, nodeT);
              int childSize = (isTree1) ? tree1NodeChildSize[nodeTParent] : tree2NodeChildSize[nodeTParent];
              int oldTreeChildDone;
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
{
#endif
              oldTreeChildDone = treeChildDone[nodeTParent];
              treeChildDone[nodeTParent]++;
#ifdef TTK_ENABLE_OPENMP
} // pragma omp atomic capture
#endif
              if(not treeNodeDone[nodeTParent] and oldTreeChildDone+1 == childSize){
                nodeT = nodeTParent;
                treeNodeDone[nodeTParent] = true;
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskyield
#endif
                //std::cout << "notify " << nodeIParent << std::endl;
              }else
                nodeT = -1;
              
            } // while nodeI loop
#ifdef TTK_ENABLE_OPENMP
} // pragma omp task
#endif
          } // while loop
#ifdef TTK_ENABLE_OPENMP
} // pragma omp parallel
#endif
        }
        
        // Equation 8, 9, 10, 11
        template<class dataType>
        void parallelEmptyTreeDistance(ftm::FTMTree_MT *tree, bool isTree1, 
                              std::vector<ftm::idNode> &treeLeaves, std::vector<int> &treeNodeChildSize,
                              std::vector<std::vector<dataType>> &treeTable,
                              std::vector<std::vector<dataType>> &forestTable,
                              std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
                              std::vector<std::vector<std::vector<std::tuple<int, int>>>> &forestBackTable){
          ftm::idNode nodeT = -1;
          std::vector<int> treeChildDone(tree->getNumberOfNodes(), 0);
          std::vector<bool> treeNodeDone(tree->getNumberOfNodes(), false);
          std::queue<ftm::idNode> treeQueue;
          for(ftm::idNode leaf : treeLeaves)
            treeQueue.push(leaf);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(NumberOfThreads) if(not IsCalled)
{
#pragma omp single nowait
#endif
          while(!treeQueue.empty()){
            nodeT = treeQueue.front();
            treeQueue.pop();
            
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(nodeT) untied //shared(treeTable, forestTable, treeBackTable, forestBackTable)
{
#endif
            while(nodeT != -1){              
              if(isTree1){
                int i = nodeT+1;
                // --- Equation 8
                computeEquation8(tree, nodeT, i, treeTable, forestTable);
                
                // --- Equation 9
                computeEquation9(tree, nodeT, i, treeTable, forestTable);
              }else{
                int j = nodeT+1;
                // --- Equation 10
                computeEquation10(tree, nodeT, j, treeTable, forestTable);
                  
                // --- Equation 11
                computeEquation11(tree, nodeT, j, treeTable, forestTable);
              }
              
              // Manage parent
              ftm::idNode nodeTParent = getParent(tree, nodeT);
              int oldTreeChildDone;
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
{
#endif
              oldTreeChildDone = treeChildDone[nodeTParent];
              treeChildDone[nodeTParent]++;
#ifdef TTK_ENABLE_OPENMP
} // pragma omp atomic capture
#endif
              if(not treeNodeDone[nodeTParent] and oldTreeChildDone+1 == treeNodeChildSize[nodeTParent]){
                nodeT = nodeTParent;
                treeNodeDone[nodeTParent] = true;
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskyield
#endif
                //std::cout << "notify " << nodeIParent << std::endl;
              }else
                nodeT = -1;

            } // while nodeI loop
#ifdef TTK_ENABLE_OPENMP
} // pragma omp task
#endif
          } // while loop
#ifdef TTK_ENABLE_OPENMP
} // pragma omp parallel
#endif
        }
        
        // ----------------------------------------
        // Other parallel version using dependency between tasks (seems wost than the above version)
        // ----------------------------------------
        // Equation 12, 13
        template<class dataType>
        void parallelDependTreeDistance(ftm::FTMTree_MT *tree1, ftm::FTMTree_MT *tree2, bool isTree1, int i,
                              std::vector<ftm::idNode> &tree1Leaves, std::vector<int> &tree1NodeChildSize,
                              std::vector<ftm::idNode> &tree2Leaves, std::vector<int> &tree2NodeChildSize,
                              std::vector<std::vector<dataType>> &treeTable,
                              std::vector<std::vector<dataType>> &forestTable,
                              std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
                              std::vector<std::vector<std::vector<std::tuple<int, int>>>> &forestBackTable){
          ftm::idNode nodeT = -1;
          ftm::FTMTree_MT *treeT = (isTree1) ? tree1 : tree2; 
          int treeChildPushed[treeT->getNumberOfNodes()] = {0};
          std::vector<bool> treeNodesPushed(treeT->getNumberOfNodes(), false);
          std::vector<ftm::idNode> *treeTLeaves = (isTree1) ? &tree1Leaves : &tree2Leaves; 
          std::vector<int> *treeTNodeChildSize = (isTree1) ? &tree1NodeChildSize : &tree2NodeChildSize; 
          std::queue<ftm::idNode> treeQueue;
          for(auto leaf : *treeTLeaves){
            treeNodesPushed[leaf] = true;
            treeChildPushed[getParent(treeT, leaf)] += 1;
            treeQueue.push(leaf);
          }
            
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(NumberOfThreads)
{
#pragma omp single nowait
#endif
          while(!treeQueue.empty()){
            nodeT = treeQueue.front();
            treeQueue.pop();
            ftm::idNode nodeTParent = getParent(treeT, nodeT);
            if(isLeaf(treeT, nodeT))
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(nodeT) depend(out:treeChildPushed[nodeTParent])
#endif
              parallelDependTreeTask(tree1, tree2, isTree1, i,
                                    tree1Leaves, tree1NodeChildSize, tree2Leaves, tree2NodeChildSize, 
                                    treeTable, forestTable, treeBackTable, forestBackTable, nodeT);
            else
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(nodeT) depend(in:treeChildPushed[nodeT]) depend(out:treeChildPushed[nodeTParent])
#endif
              parallelDependTreeTask(tree1, tree2, isTree1, i,
                                    tree1Leaves, tree1NodeChildSize, tree2Leaves, tree2NodeChildSize, 
                                    treeTable, forestTable, treeBackTable, forestBackTable, nodeT);
            // Manage parent
            if(not treeNodesPushed[nodeTParent] and 
               treeChildPushed[nodeTParent] == (*treeTNodeChildSize)[nodeTParent]){
              treeNodesPushed[nodeTParent] = true;
              treeChildPushed[getParent(treeT, nodeTParent)] += 1;
              treeQueue.push(nodeTParent);
            }
          } // while loop
#ifdef TTK_ENABLE_OPENMP
} // pragma omp parallel
#endif
        }
        
        template<class dataType>
        void parallelDependTreeTask(ftm::FTMTree_MT *tree1, ftm::FTMTree_MT *tree2, bool isTree1, int i,
                              std::vector<ftm::idNode> &tree1Leaves, std::vector<int> &tree1NodeChildSize,
                              std::vector<ftm::idNode> &tree2Leaves, std::vector<int> &tree2NodeChildSize,
                              std::vector<std::vector<dataType>> &treeTable,
                              std::vector<std::vector<dataType>> &forestTable,
                              std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
                              std::vector<std::vector<std::vector<std::tuple<int, int>>>> &forestBackTable,
                              ftm::idNode nodeT){
            int t = nodeT+1;
            
            if(isTree1){
              parallelDependTreeDistance(tree1, tree2, false, t,
                                  tree1Leaves, tree1NodeChildSize, tree2Leaves, tree2NodeChildSize, 
                                  treeTable, forestTable, treeBackTable, forestBackTable);
            }else{
              int j = nodeT+1;
              ftm::idNode nodeI = i-1;
              std::vector<ftm::idNode> children1 = getChildren(tree1, nodeI);
              std::vector<ftm::idNode> children2 = getChildren(tree2, nodeT);  
              // --- Equation 13
              computeEquation13(i, j, treeTable, forestTable, forestBackTable, children1, children2);
              
              // --- Equation 12
              computeEquation12(tree1, tree2, i, j, nodeI, nodeT, treeTable, forestTable, treeBackTable,
                                children1, children2);
            }
        }

        // ----------------------------------------
        // Utils
        // ----------------------------------------      
        void printMapIntInt(std::map<int, int> theMap){
          for (auto itr = theMap.begin(); itr != theMap.end(); ++itr) { 
            std::cout << '\t' << itr->first << '\t' << itr->second << '\n'; 
          } 
          std::cout << std::endl; 
        }
        
        // ----------------------------------------
        // Testing
        // ----------------------------------------
        template <class dataType>
        void classicalPersistenceAssignmentProblem(ftm::FTMTree_MT *tree1, ftm::FTMTree_MT *tree2){
          std::vector<std::tuple<ftm::idNode, ftm::idNode, dataType>> pairs1, pairs2;
          getPersistencePairsFromTree(tree1, pairs1);
          getPersistencePairsFromTree(tree2, pairs2);
          std::vector<std::vector<dataType>> costMatrix(pairs1.size()+1,
                                                        std::vector<dataType>(pairs2.size()+1));
          std::cout << costMatrix.size() << " _ " << costMatrix[0].size() << std::endl;
          for(int i = 0; i < costMatrix.size()-1; ++i){
            dataType nodeIValue = tree1->getValue<dataType>(std::get<0>(pairs1[i]));
            dataType nodeIOriginValue = tree1->getValue<dataType>(std::get<1>(pairs1[i]));
            for(int j = 0; j < costMatrix[0].size()-1; ++j){
              dataType nodeJValue = tree2->getValue<dataType>(std::get<0>(pairs2[j]));
              dataType nodeJOriginValue = tree2->getValue<dataType>(std::get<1>(pairs2[j]));
              costMatrix[i][j] = std::pow(nodeIValue - nodeJValue, 2) + 
                                 std::pow(nodeIOriginValue - nodeJOriginValue, 2);
            }
            costMatrix[i][costMatrix[0].size()-1] = 2*std::pow(std::get<2>(pairs1[i]), 2)/(std::pow(2, 2));
          }
          for(int j = 0; j < costMatrix[0].size()-1; ++j)
            costMatrix[costMatrix.size()-1][j] = 2*std::pow(std::get<2>(pairs2[j]), 2)/(std::pow(2, 2));
          std::vector<matchingTuple> matchings;
          forestAssignmentProblemMunkres(costMatrix, matchings);
          dataType cost = 0;
          for(auto tuple : matchings)
            cost += std::get<2>(tuple);
          std::cout << "cost      = " << cost << std::endl;
          std::cout << "cost sqrt = " << std::sqrt(cost) << std::endl;
        }
        
        template<class dataType>
        std::vector<ftm::FTMTree_MT*> makeMyTest2(){
          ftm::FTMTree_MT *tree1, *tree2;
          makeMyTest<dataType>(tree1, tree2);
          std::vector<ftm::FTMTree_MT*> vec;
          vec.push_back(tree1);
          vec.push_back(tree2);
          return vec;
        }
        
        template<class dataType>
        void makeMyTest(ftm::FTMTree_MT *&tree1, ftm::FTMTree_MT *&tree2){
          /*// Edge shifting
          float *nodesScalar1 = new float[6]{0, 1, 2, 3, 4, 5};
          std::vector<SimplexId> nodes1{0, 1, 2, 3, 4, 5};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs1{std::make_tuple(0, 4),
                                                                  std::make_tuple(1, 3),
                                                                  std::make_tuple(2, 3),
                                                                  std::make_tuple(3, 4),
                                                                  std::make_tuple(4, 5)};
          
          float *nodesScalar2 = new float[6]{0, 1, 2, 3, 4, 5};
          std::vector<SimplexId> nodes2{0, 1, 2, 3, 4, 5};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs2{std::make_tuple(0, 3),
                                                                  std::make_tuple(1, 4),
                                                                  std::make_tuple(2, 3),
                                                                  std::make_tuple(3, 4),
                                                                  std::make_tuple(4, 5)};*/
          
          ///////////////////////////////////////////////////////////////////////////////
          
          // Simple matching
          /*float *nodesScalar1 = new float[2]{0, 8};
          std::vector<SimplexId> nodes1{0, 1};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs1{std::make_tuple(0, 1)};
          
          float *nodesScalar2 = new float[2]{2, 5};
          std::vector<SimplexId> nodes2{0, 1};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs2{std::make_tuple(0, 1)};*/
          
          ///////////////////////////////////////////////////////////////////////////////
          
          /*float *nodesScalar1 = new float[8]{0, 1, 2, 8, 9, 10, 11, 20};
          std::vector<SimplexId>      nodes1{0, 1, 2, 3, 4, 5, 6, 7};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs1{std::make_tuple(0, 3),
                                                                  std::make_tuple(1, 3),
                                                                  std::make_tuple(3, 5),
                                                                  std::make_tuple(4, 5),
                                                                  std::make_tuple(5, 6),
                                                                  std::make_tuple(2, 6),
                                                                  std::make_tuple(6, 7)};
          
          float *nodesScalar2 = new float[4]{8, 9, 10, 11};
          std::vector<SimplexId>      nodes2{0, 1, 2, 3};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs2{std::make_tuple(0, 2), 
                                                                  std::make_tuple(1, 2),
                                                                  std::make_tuple(2, 3)};*/
          
          ///////////////////////////////////////////////////////////////////////////////
          
          // One node of a pair is not matched
          /*float *nodesScalar1 = new float[8]{0, 1, 2, 3, 4, 5, 6, 7};
          std::vector<SimplexId>      nodes1{0, 1, 2, 3, 4, 5, 6, 7};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs1{std::make_tuple(0, 6),
                                                                  std::make_tuple(1, 6),
                                                                  std::make_tuple(6, 7)};
          
          float *nodesScalar2 = new float[8]{0, 1, 2, 3, 4, 5, 6, 7};
          std::vector<SimplexId>      nodes2{0, 1, 2, 3, 4, 5, 6, 7};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs2{std::make_tuple(0, 6),
                                                                  std::make_tuple(5, 6),
                                                                  std::make_tuple(6, 7)};*/
          
          ///////////////////////////////////////////////////////////////////////////////
          
          /*float *nodesScalar1 = new float[8]{0, 1, 2, 3, 4, 5, 6, 7};
          std::vector<SimplexId>      nodes1{0, 1, 2, 3, 4, 5, 6, 7};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs1{std::make_tuple(0, 2),
                                                                  std::make_tuple(1, 2),
                                                                  std::make_tuple(3, 5),
                                                                  std::make_tuple(4, 5),
                                                                  std::make_tuple(2, 6),
                                                                  std::make_tuple(5, 6),
                                                                  std::make_tuple(6, 7)};
          
          float *nodesScalar2 = new float[8]{0, 1, 2, 3, 4, 5, 6, 7};
          std::vector<SimplexId>      nodes2{0, 1, 2, 3, 4, 5, 6, 7};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs2{std::make_tuple(0, 2),
                                                                  std::make_tuple(1, 2),
                                                                  std::make_tuple(2, 4),
                                                                  std::make_tuple(3, 4),
                                                                  std::make_tuple(4, 6),
                                                                  std::make_tuple(5, 6),
                                                                  std::make_tuple(6, 7)};*/

          ///////////////////////////////////////////////////////////////////////////////
          
          // BD interior destroyed node
          /*float *nodesScalar1 = new float[8]{0, 0.1, 4.75, 5.75, 6, 8};
          std::vector<SimplexId>      nodes1{0, 1, 2,    3,    4, 5};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs1{std::make_tuple(0, 4), 
                                                                  std::make_tuple(1, 3),
                                                                  std::make_tuple(2, 3),
                                                                  std::make_tuple(3, 4),
                                                                  std::make_tuple(4, 5)};

          float *nodesScalar2 = new float[8]{0, 4.75, 5.75, 8};
          std::vector<SimplexId>      nodes2{0, 1,    2,    3};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs2{std::make_tuple(0, 2),
                                                                  std::make_tuple(1, 2),
                                                                  std::make_tuple(2, 3)};*/
          
          ///////////////////////////////////////////////////////////////////////////////
          
          // Nodes of a pair map to different pairs
          /*float *nodesScalar1 = new float[8]{0, 1, 2, 3, 4, 5, 6, 7};
          std::vector<SimplexId>      nodes1{0, 1, 2, 3, 4, 5, 6, 7};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs1{std::make_tuple(0, 6), 
                                                                  std::make_tuple(5, 6),
                                                                  std::make_tuple(6, 7)};

          //float *nodesScalar2 = new float[8]{0, 1, 2, 3, 4.75, 5.75, 6, 7};
          //float *nodesScalar2 = new float[8]{0, 3, 2, 3, 4.75, 5.25, 6, 7};
          //float *nodesScalar2 = new float[8]{0, 2, 2, 3, 3.75, 5.75, 6, 7};
          float *nodesScalar2 = new float[8]{0, 2.5, 2.5, 3, 3.75, 5.75, 6, 7}; // without delete/insert
          //float *nodesScalar2 = new float[8]{0, 2, 2, 3, 3, 5, 6, 7};
          std::vector<SimplexId>      nodes2{0, 1, 2, 3, 4, 5, 6, 7};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs2{std::make_tuple(0, 6),
                                                                  std::make_tuple(1, 5),
                                                                  std::make_tuple(4, 5),                                                                  
                                                                  std::make_tuple(5, 6),
                                                                  std::make_tuple(6, 7)};
                                                                  
          // barycenter
          int mode = 0;
          if(mode == 1){
            nodesScalar1 = new float[8]{0, 3.375, 3.375, 3, 4.25, 5.375, 6, 7};
            nodes1 = std::vector<SimplexId>{0, 1, 2, 3, 4, 5, 6, 7};
            arcs1 = std::vector<std::tuple<ftm::idNode, ftm::idNode>>{std::make_tuple(0, 6),
                                                                      std::make_tuple(1, 5),
                                                                      std::make_tuple(4, 5),                                                                  
                                                                      std::make_tuple(5, 6),
                                                                      std::make_tuple(6, 7)};
          }
          if(mode == 2){
            nodesScalar2 = new float[8]{0, 3.375, 3.375, 3, 4.25, 5.375, 6, 7};
            nodes2 = std::vector<SimplexId>{0, 1, 2, 3, 4, 5, 6, 7};
            arcs2 = std::vector<std::tuple<ftm::idNode, ftm::idNode>>{std::make_tuple(0, 6),
                                                                      std::make_tuple(1, 5),
                                                                      std::make_tuple(4, 5),                                                                  
                                                                      std::make_tuple(5, 6),
                                                                      std::make_tuple(6, 7)};
          }*/

          ///////////////////////////////////////////////////////////////////////////////
          
          // Nodes of a pair map to different pairs 2
          /*float m = 40;
          float *nodesScalar1 = new float[8]{0*40, 1*40, 2*40, 3*40, 4*40, 5*40, 6*40, 7*40};
          std::vector<SimplexId>      nodes1{0, 1, 2, 3, 4, 5, 6, 7};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs1{std::make_tuple(0, 6), 
                                                                  std::make_tuple(5, 6),
                                                                  std::make_tuple(6, 7)};

          float *nodesScalar2 = new float[8]{0*40, 2.5*40, 2.5*40, 3*40, 3.75*40, 5.75*40, 6*40, 7*40}; 
          std::vector<SimplexId>      nodes2{0, 1, 2, 3, 4, 5, 6, 7};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs2{std::make_tuple(0, 6),
                                                                  std::make_tuple(1, 5),
                                                                  std::make_tuple(4, 5),                                                                  
                                                                  std::make_tuple(5, 6),
                                                                  std::make_tuple(6, 7)};
                                                                  
          // barycenter
          int mode = 2;
          if(mode == 1){
            //nodesScalar1 = new float[8]{0, 135, 135, 3*40, 170, 5.375*40, 6*40, 7*40};
            nodesScalar1 = new float[8]{0, 150, 150, 3*40, 175, 5.375*40, 6*40, 7*40};
            nodes1 = std::vector<SimplexId>{0, 1, 2, 3, 4, 5, 6, 7};
            arcs1 = std::vector<std::tuple<ftm::idNode, ftm::idNode>>{std::make_tuple(0, 6),
                                                                      std::make_tuple(1, 5),
                                                                      std::make_tuple(4, 5),                                                                  
                                                                      std::make_tuple(5, 6),
                                                                      std::make_tuple(6, 7)};
          }
          if(mode == 2){
            //nodesScalar2 = new float[8]{0, 135, 135, 3*40, 170, 5.375*40, 6*40, 7*40};
            nodesScalar2 = new float[8]{0, 150, 150, 3*40, 175, 5.375*40, 6*40, 7*40};
            nodes2 = std::vector<SimplexId>{0, 1, 2, 3, 4, 5, 6, 7};
            arcs2 = std::vector<std::tuple<ftm::idNode, ftm::idNode>>{std::make_tuple(0, 6),
                                                                      std::make_tuple(1, 5),
                                                                      std::make_tuple(4, 5),                                                                  
                                                                      std::make_tuple(5, 6),
                                                                      std::make_tuple(6, 7)};
          }*/
          
          ///////////////////////////////////////////////////////////////////////////////
          
          // Nodes of a pair map to different pairs 2
          /*float *nodesScalar1 = new float[8]{0, 1, 2, 3, 4, 5, 6, 7};
          std::vector<SimplexId>      nodes1{0, 1, 2, 3, 4, 5, 6, 7};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs1{std::make_tuple(0, 6), 
                                                                  std::make_tuple(5, 6),
                                                                  std::make_tuple(6, 7)};

          float *nodesScalar2 = new float[8]{0, 2, 3, 4, 3.75, 5.75, 6, 7};
          std::vector<SimplexId>      nodes2{0, 1, 2, 3, 4, 5, 6, 7};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs2{std::make_tuple(0, 6),
                                                                  std::make_tuple(1, 3),
                                                                  std::make_tuple(2, 3),
                                                                  std::make_tuple(3, 5),
                                                                  std::make_tuple(4, 5),                                                                  
                                                                  std::make_tuple(5, 6),
                                                                  std::make_tuple(6, 7)};*/
          
          ///////////////////////////////////////////////////////////////////////////////
          
          // Nodes of a pair map to different pairs 3
          // COUNTER EXAMPLE edit distance
          auto c = [&](int a){ return float(-a+20+1); };
          float *nodesScalar1 = new float[4]{0, c(20), c(13), 25};
          std::vector<SimplexId>      nodes1{0, 1, 2, 3};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs1{std::make_tuple(1, 0), 
                                                                  std::make_tuple(2, 1),
                                                                  std::make_tuple(3, 1)};

          float *nodesScalar2 = new float[6]{0, c(17), c(0), c(15), c(13), 25};
          std::vector<SimplexId>      nodes2{0, 1, 2, 3, 4, 5};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs2{std::make_tuple(1, 0),
                                                                  std::make_tuple(3, 1),
                                                                  std::make_tuple(2, 3),
                                                                  std::make_tuple(4, 3),
                                                                  std::make_tuple(5, 1)};
          
          ///////////////////////////////////////////////////////////////////////////////
          
          // Changement of persistence pairing during morphing
          /*float *nodesScalar1 = new float[4]{0, 1, 3, 8};
          std::vector<SimplexId> nodes1{0, 1, 2, 3};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs1{std::make_tuple(0, 2),
                                                                  std::make_tuple(1, 2),
                                                                  std::make_tuple(2, 3)};
          
          float *nodesScalar2 = new float[2]{4, 12};
          std::vector<SimplexId> nodes2{0, 1};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs2{std::make_tuple(0, 1)};*/
          
          ///////////////////////////////////////////////////////////////////////////////
          
          // Changement of persistence pairing during morphing 2
          /*float *nodesScalar1 = new float[6]{-2, 0, 1, 3, 8, 14};
          std::vector<SimplexId> nodes1{0, 1, 2, 3, 4, 5};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs1{std::make_tuple(0, 4),
                                                                  std::make_tuple(1, 3),
                                                                  std::make_tuple(2, 3),
                                                                  std::make_tuple(3, 4),
                                                                  std::make_tuple(4, 5)};
          
          float *nodesScalar2 = new float[4]{-2, 4, 12, 14};
          std::vector<SimplexId> nodes2{0, 1, 2, 3};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs2{std::make_tuple(0, 2),
                                                                  std::make_tuple(1, 2),
                                                                  std::make_tuple(2, 3)};*/
          
          ///////////////////////////////////////////////////////////////////////////////
          
          // Intermediary Merge Tree not the same matching
          /*float *nodesScalar1 = new float[6]{0, 0.01, 0.26, 0.33, 0.50, 1};
          std::vector<SimplexId>      nodes1{0, 1, 2, 3, 4, 5};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs1{std::make_tuple(0, 3),
                                                                  std::make_tuple(1, 3),
                                                                  std::make_tuple(3, 4),
                                                                  std::make_tuple(2, 4),
                                                                  std::make_tuple(4, 5)};
          
          float *nodesScalar2 = new float[6]{0, 0.35, 0.77, 0.39, 0.81, 1};
          std::vector<SimplexId>      nodes2{0, 1, 2, 3, 4, 5};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs2{std::make_tuple(0, 3), 
                                                                  std::make_tuple(1, 3),
                                                                  std::make_tuple(2, 4),
                                                                  std::make_tuple(3, 4),
                                                                  std::make_tuple(4, 5)};*/
          
          ///////////////////////////////////////////////////////////////////////////////
          
          // Multi pers pair (test: tree -> BD -> tree)
          /*float *nodesScalar1 = new float[9]{0, 1.5, 2, 3, 1, 5, 6, 3, 4};
          std::vector<SimplexId>      nodes1{0, 1, 2, 3, 4, 5, 6, 7, 8};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs1{std::make_tuple(0, 5),
                                                                  std::make_tuple(1, 3),
                                                                  std::make_tuple(2, 3),
                                                                  std::make_tuple(3, 5),
                                                                  std::make_tuple(7, 8),
                                                                  std::make_tuple(4, 8),
                                                                  std::make_tuple(8, 5),
                                                                  std::make_tuple(5, 6)};
          
          float *nodesScalar2 = new float[6]{0, 0.35, 0.77, 0.39, 0.81, 1};
          std::vector<SimplexId>      nodes2{0, 1, 2, 3, 4, 5};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs2{std::make_tuple(0, 3), 
                                                                  std::make_tuple(1, 3),
                                                                  std::make_tuple(2, 4),
                                                                  std::make_tuple(3, 4),
                                                                  std::make_tuple(4, 5)};*/
          
          ///////////////////////////////////////////////////////////////////////////////
          
          ftm::FTMTree_MT *treeTemp1 = makeFakeTree(nodesScalar1, nodes1, arcs1);
          tree1 = treeTemp1;
          ftm::FTMTree_MT *treeTemp2 = makeFakeTree(nodesScalar2, nodes2, arcs2);
          tree2 = treeTemp2;
        }

    }; // MergeTreeDistance class
    
} // namespace ttk

#endif
