/// Provide your information
///
/// \ingroup base
/// \class ttk::FTMTreeEditDistanceClustering
/// \author Mathieu Pont <mathieu.pont@outlook.com>
/// \date 2020.
///
/// This module defines the %FTMTreeEditDistanceClustering class that computes 
/// a clustering of an ensemble of merge trees in k clusters
/// It is actually a version of KMeans (dynamic clustering algorithm) 
/// where centroids are merge trees 
///

#ifndef _FTMTREEEDITDISTANCECLUSTERING_H
#define _FTMTREEEDITDISTANCECLUSTERING_H

#define treesMatchingVector std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>>
#define matchingVector std::vector<treesMatchingVector>

#pragma once

#include <random>

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>

#include "FTMTreeEditDistanceBarycenter.h"

namespace ttk {
  
    /**
    * The %FTMTreeEditDistanceClustering class that computes 
    * a clustering of an ensemble of merge trees in k clusters
    */
    class FTMTreeEditDistanceClustering : virtual public Debug, public FTMTreeEditDistanceBarycenter {
      
    private:
      int noCentroids = 2;
      bool ParallelizeUpdate = false;
      
      // Accelerated KMeans
      bool acceleratedInitialized = false;
      std::vector<std::vector<double>> lowerBound;
      std::vector<double> upperBound;
      std::vector<int> bestCentroid, oldBestCentroid;
      std::vector<double> bestDistance;
      std::vector<bool> recompute;
      std::vector<MergeTree*> oldCentroids, oldCentroids2;

    public:
        FTMTreeEditDistanceClustering() {
          this->setDebugMsgPrefix(
              "FTMTreeEditDistanceClustering"); // inherited from Debug: prefix will be printed at the
          // beginning of every msg
          //IsCalled = true; // for the call of barycenter functions
        };
        ~FTMTreeEditDistanceClustering(){};
        
        void setNoCentroids(int noCentroidsT){
          noCentroids = noCentroidsT;
        }
        
        /**
        * Implementation of the algorithm.
        */
        
        // KMeans++ init
        template <class dataType>
        std::vector<std::vector<MergeTree*>> initCentroids(std::vector<ftm::FTMTree_MT*> &trees,
                                                           std::vector<ftm::FTMTree_MT*> &trees2){
          std::vector<MergeTree*> centroids(noCentroids), centroids2(noCentroids);
          std::vector<dataType> distances(trees.size(), std::numeric_limits<dataType>::max() );
          for(int i = 0; i < noCentroids; ++i){
            if(i == 0){
              int bestIndex = getBestInitTreeIndex(trees);
//std::cout << i << " _ " << bestIndex << std::endl;
              centroids[i] = copyMergeTree<dataType>(trees[bestIndex], true);
              cleanMergeTree<dataType>(centroids[i]);
              if(trees2.size() != 0){
                centroids2[i] = copyMergeTree<dataType>(trees2[bestIndex], true);
                cleanMergeTree<dataType>(centroids2[i]);
              }              
            }else{
              // Create vector of probabilities 
              dataType sum = 0;
              for(auto val : distances)
                sum += val;
              int bestIndex = -1;
              dataType bestValue = std::numeric_limits<dataType>::min();
              std::vector<dataType> probabilities(trees.size());
              for(int j = 0; j < distances.size(); ++j){
                probabilities[j] = distances[j] / sum;
                if(probabilities[j] > bestValue){
                  bestValue = probabilities[j];
                  bestIndex = j;
                }
              }
              if(not Deterministic){
                std::random_device rd;
                std::default_random_engine generator(rd());
                std::discrete_distribution<int> distribution(probabilities.begin(), probabilities.end());
                bestIndex = distribution(generator);
              }
//std::cout << i << " _ " << bestIndex << std::endl;
              // Create new centroid
              centroids[i] = copyMergeTree<dataType>(trees[bestIndex], true);
              cleanMergeTree<dataType>(centroids[i]);
              if(trees2.size() != 0){
                centroids2[i] = copyMergeTree<dataType>(trees2[bestIndex], true);
                cleanMergeTree<dataType>(centroids2[i]);
              }
            }
            
            if(i == noCentroids-1)
              continue;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(NumberOfThreads) if(Parallelize)
#endif
            for(int j = 0; j < trees.size(); ++j){
              std::vector<std::tuple<ftm::idNode, ftm::idNode>> matching, matching2;
              dataType distanceT, distanceT2;
              computeOneDistance<dataType>(trees[j], centroids[i], matching, distanceT);
              if(trees2.size() != 0){
                computeOneDistance<dataType>(trees2[j], centroids2[i], matching2, distanceT2);
                distanceT += distanceT2;
              }
              distances[j] = std::min(distances[j], distanceT);
            }
          }

          std::vector<std::vector<MergeTree*>> allCentroids;
          allCentroids.push_back(centroids);
          allCentroids.push_back(centroids2);
          return allCentroids;
        }
        
        template <class dataType>
        MergeTree* initNewCentroid(std::vector<ftm::FTMTree_MT*> &trees){
          MergeTree *centroid;
          int bestIndex = -1;
          dataType bestValue = std::numeric_limits<dataType>::min();
          for(int i = 0; i < bestDistance.size(); ++i){
            if(bestValue < bestDistance[i]){
              bestValue = bestDistance[i];
              bestIndex = bestCentroid[i];
            }
          }
          centroid = copyMergeTree<dataType>(trees[bestIndex], true);
          cleanMergeTree<dataType>(centroid);
          
          return centroid;
        }

        template <class dataType>
        void assignmentCentroidsNaive(std::vector<ftm::FTMTree_MT*> &trees, 
                                      std::vector<MergeTree*> &centroids,
                                      std::vector<std::tuple<int, int>> &assignmentC, 
                                      std::vector<dataType> &bestDistanceT,
                                      std::vector<ftm::FTMTree_MT*> &trees2,
                                      std::vector<MergeTree*> &centroids2){
          std::vector<int> bestCentroidT(trees.size(), -1);
          
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(NumberOfThreads) if(Parallelize)
#endif
          for(int i = 0; i < trees.size(); ++i){
            for(int j = 0; j < centroids.size(); ++j){
              std::vector<std::tuple<ftm::idNode, ftm::idNode>> matching;
              std::vector<std::tuple<ftm::idNode, ftm::idNode>> matching2;
              dataType distance, distance2;
              computeOneDistance<dataType>(trees[i], centroids[j], matching, distance);
              if(trees2.size() != 0){
                computeOneDistance<dataType>(trees2[i], centroids2[j], matching2, distance2);
                distance += distance2;
              }
              if(distance < bestDistanceT[i]){
                bestDistanceT[i] = distance;
                bestDistance[i] = distance;
                bestCentroidT[i] = j;
                bestCentroid[i] = j;
              }
            }
          }
          
          for(int i = 0; i < bestCentroidT.size(); ++i)
            assignmentC.push_back(std::make_tuple(bestCentroidT[i], i));
            
        }
        
        template <class dataType>
        void getCentroidsDistanceMatrix(std::vector<MergeTree*> &centroids, 
                                        std::vector<std::vector<double>> &distanceMatrix){
          distanceMatrix = std::vector<std::vector<double>>(centroids.size(), 
                                                            std::vector<double>(centroids.size(), 0));
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(NumberOfThreads) if(Parallelize)
#endif
          for(int i = 0; i < centroids.size(); ++i)
            for(int j = i+1; j < centroids.size(); ++j){
              std::vector<std::tuple<ftm::idNode, ftm::idNode>> matching;
              dataType distance;
              computeOneDistance<dataType>(centroids[i], centroids[j], matching, distance);
              distanceMatrix[i][j] = distance;
              distanceMatrix[j][i] = distance;
            }
        }
        
        void hadamardProduct(std::vector<std::vector<double>> &distanceMatrix, 
                           std::vector<std::vector<double>> &distanceMatrix2){
          for(int i = 0; i < distanceMatrix.size(); ++i)
            for(int j = 0; j < distanceMatrix[i].size(); ++j)
              distanceMatrix[i][j] += distanceMatrix2[i][j];
        }
        
        template <class dataType>
        void initAcceleratedKMeansVectors(std::vector<ftm::FTMTree_MT*> &trees, 
                                          std::vector<MergeTree*> &centroids,
                                          std::vector<ftm::FTMTree_MT*> &trees2){
          lowerBound = std::vector<std::vector<double>>(trees.size(),
                                                        std::vector<double>(centroids.size(), 0));
          upperBound = std::vector<double>(trees.size(), std::numeric_limits<double>::max());
          bestCentroid = std::vector<int>(trees.size(), -1);
          oldBestCentroid = std::vector<int>(trees.size(), -1);
          bestDistance = std::vector<double>(trees.size(), std::numeric_limits<double>::max());
          recompute = std::vector<bool>(trees.size(), true);
        }
        
        template <class dataType>
        void initAcceleratedKMeans(std::vector<ftm::FTMTree_MT*> &trees, std::vector<MergeTree*> &centroids,
                                std::vector<ftm::FTMTree_MT*> &trees2, std::vector<MergeTree*> &centroids2){
          acceleratedInitialized = true;          
          std::vector<std::tuple<int, int>> assignmentC;
          std::vector<dataType> bestDistanceT(trees.size(), std::numeric_limits<dataType>::max());
          assignmentCentroidsNaive<dataType>(trees, centroids, assignmentC, bestDistanceT, 
                                             trees2, centroids2);
          for(int i = 0; i < bestDistanceT.size(); ++i)
            bestDistance[i] = bestDistanceT[i];
          for(auto asgn : assignmentC)
            bestCentroid[std::get<1>(asgn)] = std::get<0>(asgn);
          for(int i = 0; i < bestDistance.size(); ++i)
            upperBound[i] = bestDistance[i];
        }
        
        template <class dataType>
        void copyCentroids(std::vector<MergeTree*> &centroids, std::vector<MergeTree*> &oldCentroids){
          for(int i = 0; i < oldCentroids.size(); ++i)
            freeMergeTree<dataType>(oldCentroids[i]);
          oldCentroids.clear();
          for(int i = 0; i < centroids.size(); ++i)
            oldCentroids.push_back(copyMergeTree<dataType>(centroids[i]));
        }
        
        template <class dataType>
        void assignmentCentroidsAccelerated(std::vector<ftm::FTMTree_MT*> &trees, 
                                            std::vector<MergeTree*> &centroids, 
                                            std::vector<std::tuple<int, int>> &assignmentC, 
                                            std::vector<dataType> &bestDistanceT, 
                                            std::vector<ftm::FTMTree_MT*> &trees2, 
                                            std::vector<MergeTree*> &centroids2){
          if(not acceleratedInitialized){
            initAcceleratedKMeans<dataType>(trees, centroids, trees2, centroids2);
          }else{
            // Compute distance between old and new corresponding centroids
            std::vector<dataType> distanceShift(centroids.size()), distanceShift2(centroids2.size());
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(NumberOfThreads) if(Parallelize)
#endif
            for(int i = 0; i < centroids.size(); ++i){
              std::vector<std::tuple<ftm::idNode, ftm::idNode>> matching, matching2;
              computeOneDistance<dataType>(centroids[i], oldCentroids[i], matching, distanceShift[i]);
              if(trees2.size() != 0){
                computeOneDistance<dataType>(centroids2[i], oldCentroids2[i], matching2, distanceShift2[i]);
                distanceShift[i] += distanceShift2[i];
              }
            }
            
            // Step 5
            for(int i = 0; i < trees.size(); ++i)
              for(int c = 0; c < centroids.size(); ++c)
                lowerBound[i][c] = std::max(lowerBound[i][c] - distanceShift[c], 0.0);
              
            // Step 6
            for(int i = 0; i < trees.size(); ++i){
              upperBound[i] = upperBound[i] + distanceShift[bestCentroid[i]];
              recompute[i] = true;
            }
          }
          
          // Step 1
          std::vector<std::vector<double>> centroidsDistance, centroidsDistance2;
          getCentroidsDistanceMatrix<dataType>(centroids, centroidsDistance);
          if(trees2.size() != 0){
            getCentroidsDistanceMatrix<dataType>(centroids2, centroidsDistance2);
            hadamardProduct(centroidsDistance, centroidsDistance2);
          }
          std::vector<double> centroidScore(centroids.size(), std::numeric_limits<double>::max());
          for(int i = 0; i < centroids.size(); ++i)
            for(int j = i+1; j < centroids.size(); ++j){
              if(0.5 * centroidsDistance[i][j] < centroidScore[i])
                centroidScore[i] = 0.5 * centroidsDistance[i][j];
              if(0.5 * centroidsDistance[i][j] < centroidScore[j])
                centroidScore[j] = 0.5 * centroidsDistance[i][j];
            }
              
          // Step 2
          std::vector<bool> identified(trees.size());
          for(int i = 0; i < trees.size(); ++i)
            identified[i] = (upperBound[i] <= centroidScore[bestCentroid[i]]);

          // Step 3
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(NumberOfThreads) if(Parallelize)
#endif
          for(int i = 0; i < trees.size(); ++i)
            for(int c = 0; c < centroids.size(); ++c){
              if(not identified[i] and c != bestCentroid[i] and upperBound[i] > lowerBound[i][c] and 
                                       upperBound[i] > 0.5 * centroidsDistance[bestCentroid[i]][c]){
                // Step 3a
                if(recompute[i]){
                  std::vector<std::tuple<ftm::idNode, ftm::idNode>> matching, matching2;
                  dataType distance, distance2;
                  computeOneDistance<dataType>(trees[i], centroids[bestCentroid[i]], matching, distance);
                  if(trees2.size() != 0){
                    computeOneDistance<dataType>(trees2[i],centroids2[bestCentroid[i]],matching2,distance2);
                    distance += distance2;
                  }
                  recompute[i] = false;
                  lowerBound[i][bestCentroid[i]] = distance;
                  upperBound[i] = distance;
                  bestDistance[i] = distance;
                }else{
                  bestDistance[i] = upperBound[i];
                }
                // Step 3b
                if(bestDistance[i] > lowerBound[i][c] and 
                   bestDistance[i] > 0.5 * centroidsDistance[bestCentroid[i]][c]){
                  std::vector<std::tuple<ftm::idNode, ftm::idNode>> matching, matching2;
                  dataType distance, distance2;
                  computeOneDistance<dataType>(trees[i], centroids[c], matching, distance);
                  if(trees2.size() != 0){
                    computeOneDistance<dataType>(trees2[i], centroids2[c], matching2, distance2);
                    distance += distance2;
                  }
                  lowerBound[i][c] = distance;
                  if(distance < bestDistance[i]){
                    bestCentroid[i] = c;
                    upperBound[i] = distance;
                    bestDistance[i] = distance;
                  }
                }
              }
            }
          
          // Copy centroids for next step
          copyCentroids<dataType>(centroids, oldCentroids);
          if(trees2.size() != 0)
            copyCentroids<dataType>(centroids2, oldCentroids2);
          
          // Manage output
          for(int i = 0; i < bestDistance.size(); ++i)
            bestDistanceT[i] = bestDistance[i];
          for(int i = 0; i < bestCentroid.size(); ++i)
            assignmentC.push_back(std::make_tuple(bestCentroid[i], i));
        }
        
        template <class dataType>
        void assignmentCentroids(std::vector<ftm::FTMTree_MT*> &trees, std::vector<MergeTree*> &centroids, 
                        std::vector<std::tuple<int, int>> &assignmentC, std::vector<dataType> &bestDistanceT, 
                        std::vector<ftm::FTMTree_MT*> &trees2, std::vector<MergeTree*> &centroids2){
          oldBestCentroid = bestCentroid;
          if(NormalizedWasserstein and RescaledWasserstein)
            assignmentCentroidsNaive<dataType>(trees, centroids, assignmentC, bestDistanceT, 
                                               trees2, centroids2);
          else
            assignmentCentroidsAccelerated<dataType>(trees, centroids, assignmentC, bestDistanceT, 
                                                     trees2, centroids2);
        }
        
        template <class dataType>
        void finalAssignmentCentroids(std::vector<ftm::FTMTree_MT*> &trees, 
                                 std::vector<MergeTree*> &centroids, matchingVector &matchingsC, 
                                 std::vector<std::tuple<int, int>> &assignmentC,
                                 std::vector<dataType> &bestDistanceT, std::vector<ftm::FTMTree_MT*> &trees2, 
                                 std::vector<MergeTree*> &centroids2, matchingVector &matchingsC2){
          int noC = centroids.size();
          std::vector<std::vector<ftm::FTMTree_MT*>> assignedTrees(noC), assignedTrees2(noC);
          std::vector<std::vector<int>> assignedTreesIndex(noC);
          
          for(auto asgn : assignmentC){
            assignedTreesIndex[std::get<0>(asgn)].push_back(std::get<1>(asgn));
            assignedTrees[std::get<0>(asgn)].push_back(trees[std::get<1>(asgn)]);
            if(trees2.size() != 0)
              assignedTrees2[std::get<0>(asgn)].push_back(trees2[std::get<1>(asgn)]);
          }
          
          for(int i = 0; i < centroids.size(); ++i){
            std::vector<dataType> distances(assignedTrees[i].size(), 0); 
            std::vector<dataType> distances2(assignedTrees[i].size(), 0);
            treesMatchingVector matching(trees.size()), matching2(trees2.size());
            assignment<dataType>(assignedTrees[i], centroids[i], matching, distances);
            matchingsC[i] = matching;
            if(trees2.size() != 0){
              assignment<dataType>(assignedTrees2[i], centroids2[i], matching2, distances2);
              matchingsC2[i] = matching2;
            }
            for(int j = 0; j < assignedTreesIndex[i].size(); ++j){
              int index = assignedTreesIndex[i][j];
              bestDistanceT[index] = distances[j] + distances2[j];
            }
              
          }
        }
        
        void matchingCorrespondence(treesMatchingVector &matchingT, std::vector<int> &nodeCorr, 
                                    std::vector<int> assignedTreesIndex){
          for(int i : assignedTreesIndex){
            std::vector<std::tuple<ftm::idNode, ftm::idNode>> newMatching;
            for(auto tup : matchingT[i])
              newMatching.push_back(std::make_tuple(nodeCorr[std::get<0>(tup)], std::get<1>(tup)));
            matchingT[i] = newMatching;
          }
        }
        
        bool samePreviousAssignment(int clusterId){
          for(int i = 0; i < bestCentroid.size(); ++i)
            if(bestCentroid[i] == clusterId and bestCentroid[i] != oldBestCentroid[i])
              return false;
          return true;
        }

        template <class dataType>
        void updateCentroids(std::vector<ftm::FTMTree_MT*> &trees, std::vector<MergeTree*> &centroids, 
                             std::vector<double> &alphas, std::vector<std::tuple<int, int>> &assignmentC){
          int noC = centroids.size();
          std::vector<std::vector<ftm::FTMTree_MT*>> assignedTrees(noC);
          std::vector<std::vector<int>> assignedTreesIndex(noC);
          std::vector<std::vector<double>> assignedAlphas(noC);
          
          for(auto asgn : assignmentC){
            assignedTrees[std::get<0>(asgn)].push_back(trees[std::get<1>(asgn)]);
            assignedTreesIndex[std::get<0>(asgn)].push_back(std::get<1>(asgn));
            assignedAlphas[std::get<0>(asgn)].push_back(alphas[std::get<1>(asgn)]);
          }
#ifdef TTK_ENABLE_OPENMP
//#pragma omp parallel for schedule(dynamic) num_threads(NumberOfThreads) if(Parallelize and ParallelizeUpdate)
#pragma omp parallel num_threads(NumberOfThreads) if(Parallelize and ParallelizeUpdate)
{
#pragma omp single nowait
{
#endif
          for(int i = 0; i < centroids.size(); ++i){
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(i)
{
#endif
/*std::stringstream ss;
ss << "update centroid " << i << std::endl;
std::cout << ss.str();*/

            // Init new centroid if no trees are assigned to it
            if(assignedTrees[i].size() == 0){
              centroids[i] = initNewCentroid<dataType>(trees);
              for(int t = 0; t < trees.size(); ++t)
                lowerBound[t][i] = 0;
            }else
            // Do not update if same previous assignment
            if(not samePreviousAssignment(i)){
              // Otherwise compute barycenter of the assigned trees
              double alphasSum = 0;
              for(int j = 0; j < assignedAlphas[i].size(); ++j)
                alphasSum += assignedAlphas[i][j];
              for(int j = 0; j < assignedAlphas[i].size(); ++j)
                assignedAlphas[i][j] /= alphasSum;
              treesMatchingVector matching(assignedTrees[i].size());
              centroids[i] = computeBarycenter<dataType>(assignedTrees[i], centroids[i], 
                                                         assignedAlphas[i], matching);
              persistenceThresholding<dataType>(centroids[i]->tree, 0);
              cleanMergeTree<dataType>(centroids[i]);
            }
#ifdef TTK_ENABLE_OPENMP
} // pragma omp task
#endif
          }
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
        std::vector<MergeTree*> computeCentroids(std::vector<ftm::FTMTree_MT*> &trees, 
                                    std::vector<MergeTree*> &centroids, matchingVector &outputMatching, 
                                    std::vector<double> &alphas, std::vector<int> &clusteringAssignment,
                                    std::vector<ftm::FTMTree_MT*> &trees2,
                                    std::vector<MergeTree*> &centroids2, matchingVector &outputMatching2){
          Timer t_clust;
          
          int noCentroids = centroids.size();
          bool converged = false;
          dataType inertia = -1;
          dataType minInertia = std::numeric_limits<dataType>::max();
          int cptBlocked = 0;
          std::vector<std::tuple<int, int>> assignmentC;
          std::vector<dataType> bestDistanceT(trees.size(), std::numeric_limits<dataType>::max());
          while(not converged){
            std::cout << "==================" << std::endl;
            // --- Assignment
            Timer t_assignment;
            assignmentCentroids<dataType>(trees, centroids, assignmentC, bestDistanceT, trees2, centroids2);
            auto t_assignment_time = t_assignment.getElapsedTime();
            std::cout << "assignment : " << t_assignment_time << std::endl;
            
            // --- Update
            Timer t_update;
            updateCentroids<dataType>(trees, centroids, alphas, assignmentC);
            if(trees2.size() != 0)
              updateCentroids<dataType>(trees2, centroids2, alphas, assignmentC);
            auto t_update_time = t_update.getElapsedTime();
            std::cout << "update     : " << t_update_time << std::endl;
            
            // --- Check convergence
            dataType currentInertia = 0;
            for(auto distance : bestDistanceT)
              currentInertia += distance*distance;
            converged = myAbs<dataType>(inertia - currentInertia) < 0.01;
            inertia = currentInertia;
            std::cout << "Inertia : " << inertia << std::endl;
            
            minInertia = std::min(minInertia, inertia);
            if(not converged){
              cptBlocked += (minInertia < inertia) ? 1 : 0;
              converged = (cptBlocked >= 10);
            }
            
            // --- Reset vectors
            if(not converged){
              assignmentC.clear();
              bestDistanceT = std::vector<dataType>(trees.size(), std::numeric_limits<dataType>::max());
            }
          }
          
          /*for(int i = 0; i < centroids.size(); ++i){
            centroids[i]->tree->printTree2(); 
            printTreeScalars<dataType>(centroids[i]->tree); 
            printPairsFromTree<dataType>(centroids[i]->tree, BranchDecomposition);
          }*/
          
          matchingVector matchingsC(noCentroids);
          matchingVector matchingsC2(noCentroids);
          finalAssignmentCentroids<dataType>(trees, centroids, matchingsC, assignmentC, bestDistanceT,
                                             trees2, centroids2, matchingsC2);
          dataType currentInertia = 0;
          for(auto distance : bestDistanceT)
            currentInertia += distance*distance;
          std::cout << "Inertia : " << currentInertia << std::endl;
          
          // Manage output
          std::vector<int> cptCentroid(centroids.size(), 0);
          for(auto asgn : assignmentC){
            int centroid = std::get<0>(asgn);
            int tree = std::get<1>(asgn);
//std::cout << centroid << " " << tree << std::endl;
            clusteringAssignment[tree] = centroid;
            outputMatching[centroid][tree] = matchingsC[centroid][cptCentroid[centroid]];
            if(trees2.size() != 0)
              outputMatching2[centroid][tree] = matchingsC2[centroid][cptCentroid[centroid]];
            ++cptCentroid[centroid];
          }
          
          std::cout << "TIME CLUSTERING = " << t_clust.getElapsedTime() << std::endl;
          
          return centroids;
        }
        
        template <class dataType>
        void preprocessingClustering(std::vector<ftm::FTMTree_MT*> &trees){
          int avgNodes = 0, avgNodesT = 0;
          double avgDepth=0;
          for(int i = 0; i < trees.size(); ++i){
            preprocessingPipeline<dataType>(trees[i], EpsilonTree2, Epsilon2Tree2, 0);
            trees[i] = cleanMergeTree<dataType>(trees[i]);
            auto noNodesT = trees[i]->getNumberOfNodes();
            auto noNodes = getRealNumberOfNodes(trees[i]);
            //std::cout << "tree " << i << " number of nodes : " << noNodes << " / " << noNodesT << std::endl;
            avgNodes += noNodes;
            avgNodesT += noNodesT;
            avgDepth += getTreeDepth(trees[i]);
          }
          avgNodes /= trees.size();
          avgNodesT /= trees.size();
          avgDepth /= trees.size();
          std::cout << trees.size() << " trees average [node: " << avgNodes << " / "<< avgNodesT;
          std::cout << ", depth: " << avgDepth << "]" << std::endl;
        }

        template <class dataType>
        void postprocessingClustering(std::vector<ftm::FTMTree_MT*> &trees, 
                                      std::vector<MergeTree*> centroids, matchingVector &outputMatching,
                                      std::vector<int> &clusteringAssignment){
          for(int i = 0; i < trees.size(); ++i)
            postprocessingPipeline<dataType>(trees[i]);
          for(int i = 0; i < centroids.size(); ++i)
            postprocessingPipeline<dataType>(centroids[i]->tree);
          for(int c = 0; c < centroids.size(); ++c)
            for(int i = 0; i < trees.size(); ++i)
              if(clusteringAssignment[i] == c)
                convertBranchDecompositionMatching(centroids[c]->tree, trees[i], outputMatching[c][i]);
        }
        
        template <class dataType>
        std::vector<MergeTree*> execute(std::vector<ftm::FTMTree_MT*> &trees, matchingVector &outputMatching,
                                    std::vector<double> &alphas, std::vector<int> &clusteringAssignment,
                                    std::vector<ftm::FTMTree_MT*> &trees2, matchingVector &outputMatching2){
          // --- Preprocessing
          preprocessingClustering<dataType>(trees);
          if(trees2.size() != 0)
            preprocessingClustering<dataType>(trees2);
          
          //std::cout << "compute centroids..." << std::endl;
          auto allCentroids = initCentroids<dataType>(trees, trees2);
          //std::cout << "compute centroids done." << std::endl;
          std::vector<MergeTree*> centroids = allCentroids[0];
          std::vector<MergeTree*> centroids2;
          if(trees2.size() != 0)
            centroids2 = allCentroids[1];
          /*for(int i = 0; i < centroids.size(); ++i){
            verifyBranchDecompositionInconsistency<dataType>(centroids[i]->tree);
            if(trees2.size() != 0)
              verifyBranchDecompositionInconsistency<dataType>(centroids2[i]->tree);
          }*/
          initAcceleratedKMeansVectors<dataType>(trees, centroids, trees2);
          
          // --- Execute
          centroids = computeCentroids<dataType>(trees, centroids, outputMatching, alphas, 
                                                 clusteringAssignment, trees2, centroids2, outputMatching2);
          
          // --- Postprocessing
          postprocessingClustering<dataType>(trees, centroids, outputMatching, clusteringAssignment);
          if(trees2.size() != 0)
            postprocessingClustering<dataType>(trees2, centroids2, outputMatching2, clusteringAssignment);
              
          return centroids;
        }
        
        template <class dataType>
        std::vector<MergeTree*> execute(std::vector<ftm::FTMTree_MT*> &trees, matchingVector &outputMatching,
                                    std::vector<int> &clusteringAssignment, 
                                    std::vector<ftm::FTMTree_MT*> &trees2, matchingVector &outputMatching2){
          //trees = makeMyTestClust<dataType>();
          if(trees2.size() != 0)
            std::cout << "Use join and split trees" << std::endl;
          
          std::vector<double> alphas;
          for(int i = 0; i < trees.size(); ++i)
            alphas.push_back(1.0/trees.size());
          
          return execute<dataType>(trees, outputMatching, alphas, clusteringAssignment, trees2, 
                                   outputMatching2);
        }
        
        template <class dataType>
        std::vector<MergeTree*> execute(std::vector<ftm::FTMTree_MT*> &trees, matchingVector &outputMatching,
                                        std::vector<int> &clusteringAssignment){
          std::vector<ftm::FTMTree_MT*> trees2 = std::vector<ftm::FTMTree_MT*>();
          matchingVector outputMatching2 = matchingVector();
          return execute<dataType>(trees, outputMatching, clusteringAssignment, trees2, outputMatching2);
        }
        
        // ----------------------------------------
        // Testing
        // ----------------------------------------
        template<class dataType>
        std::vector<ftm::FTMTree_MT *> makeMyTestClust(){
          std::vector<ftm::FTMTree_MT *> trees;
          
          float *nodesScalar1 = new float[4]{0, 1, 3, 8};
          std::vector<SimplexId> nodes1{0, 1, 2, 3};
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs1{std::make_tuple(0, 2),
                                                                  std::make_tuple(1, 2),
                                                                  std::make_tuple(2, 3)};
          ftm::FTMTree_MT *treeTemp1 = makeFakeTree(nodesScalar1, nodes1, arcs1);
          trees.push_back(treeTemp1);
          
          float *nodesScalar2 = new float[4]{1, 2, 4, 9};
          std::vector<SimplexId> nodes2(nodes1);
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs2(arcs1);
          ftm::FTMTree_MT *treeTemp2 = makeFakeTree(nodesScalar2, nodes2, arcs2);
          trees.push_back(treeTemp2);
          
          float *nodesScalar3 = new float[4]{3, 4, 6, 11};
          std::vector<SimplexId> nodes3(nodes1);
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs3(arcs1);
          ftm::FTMTree_MT *treeTemp3 = makeFakeTree(nodesScalar3, nodes3, arcs3);
          trees.push_back(treeTemp3);
          
          float *nodesScalar4 = new float[4]{30, 40, 60, 110};
          std::vector<SimplexId> nodes4(nodes1);
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs4(arcs1);
          ftm::FTMTree_MT *treeTemp4 = makeFakeTree(nodesScalar4, nodes4, arcs4);
          trees.push_back(treeTemp4);
          
          float *nodesScalar5 = new float[4]{10, 20, 40, 90};
          std::vector<SimplexId> nodes5(nodes1);
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs5(arcs1);
          ftm::FTMTree_MT *treeTemp5 = makeFakeTree(nodesScalar5, nodes5, arcs5);
          trees.push_back(treeTemp5);
          
          float *nodesScalar6 = new float[4]{0, 10, 30, 80};
          std::vector<SimplexId> nodes6(nodes1);
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> arcs6(arcs1);
          ftm::FTMTree_MT *treeTemp6 = makeFakeTree(nodesScalar6, nodes6, arcs6);
          trees.push_back(treeTemp6);
          
          return trees;
        }
        
    }; // FTMTreeEditDistanceBarycenter class
    
} // namespace ttk

#endif
