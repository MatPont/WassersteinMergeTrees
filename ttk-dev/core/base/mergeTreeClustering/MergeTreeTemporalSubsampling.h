/// Provide your information
///
/// \ingroup base
/// \class ttk::MergeTreeTemporalSubsampling
/// \author XXXXX
/// \date 2021.
///
/// This module defines the %MergeTreeTemporalSubsampling class that computes 
///

#ifndef _MERGETREETEMPORALSUBSAMPLING_H
#define _MERGETREETEMPORALSUBSAMPLING_H

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>
#include <LDistance.h>

#include <vtkUnstructuredGrid.h>
#include <vtkImageData.h>
#include <vtkResampleToImage.h>
#include <vtkAppendFilter.h>
#include <vtkArrayCalculator.h>
#include <vtkPointData.h>

#include <vtkCharArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkDataArray.h>

namespace ttk {
  
    /**
    * The %MergeTreeTemporalSubsampling class that computes 
    * 
    */
    class MergeTreeTemporalSubsampling : virtual public Debug, public MergeTreeBase {
      
    private:
      double RemovalPercentage = 50.;
      bool UseL2Distance = false;
      std::vector<vtkDataSet *> treesSegmentation;
      std::vector<vtkSmartPointer<vtkDataSet>> images;
      bool DoResampleToImage = false;
      
    public:
        MergeTreeTemporalSubsampling() {
          this->setDebugMsgPrefix(
              "MergeTreeTemporalSubsampling"); // inherited from Debug: prefix will be printed at 
          // the beginning of every msg
        };
        ~MergeTreeTemporalSubsampling(){};
        
        void setRemovalPercentage(double rs){
          RemovalPercentage = rs;
        }
        
        void setUseL2Distance(bool useL2){
          UseL2Distance = useL2;
        }
        
        void setTreesSegmentation(std::vector<vtkDataSet *> &segmentation){
          treesSegmentation = segmentation;
        }
        
        void preprocessSegmentation(){
          for(int i = 0; i < treesSegmentation.size(); ++i){
            if(DoResampleToImage){
              auto seg = treesSegmentation[i];            
              //auto gridSeg = vtkUnstructuredGrid::SafeDownCast(seg);
              auto resampleFilter = vtkSmartPointer<vtkResampleToImage>::New();
              resampleFilter->SetInputDataObject(seg);
              resampleFilter->SetUseInputBounds(true);
              //resampleFilter->SetSamplingDimensions(100, 100, 100);
              resampleFilter->Update();
              auto resampleFilterOut = resampleFilter->GetOutput();
              //images.push_back(resampleFilterOut);
              vtkSmartPointer<vtkDataSet> image = vtkSmartPointer<vtkImageData>::New();
              image->DeepCopy(resampleFilterOut);
              images.push_back(image);
            }else{
              vtkSmartPointer<vtkDataSet> image;
              image.TakeReference(treesSegmentation[i]);
              images.push_back(image);
            }
          }
        }
        
        /**
        * Implementation of the algorithm.
        */
        vtkSmartPointer<vtkDataArray> createVtkDataArray(int dataType){
          vtkSmartPointer<vtkDataArray> outputScalarField;
          switch(dataType) {
            case VTK_CHAR:
              outputScalarField = vtkSmartPointer<vtkCharArray>::New();
              break;
            case VTK_DOUBLE:
              outputScalarField = vtkSmartPointer<vtkDoubleArray>::New();
              break;
            case VTK_FLOAT:
              outputScalarField = vtkSmartPointer<vtkFloatArray>::New();
              break;
            case VTK_INT:
              outputScalarField = vtkSmartPointer<vtkIntArray>::New();
              break;
            case VTK_ID_TYPE:
              outputScalarField = vtkSmartPointer<vtkIdTypeArray>::New();
              break;
            default: 
              outputScalarField = vtkSmartPointer<vtkDoubleArray>::New();
          }
          return outputScalarField;
        }
        
        vtkSmartPointer<vtkDataArray> createL2OutputScalarField(int dataType, int numberOfPoints){
          vtkSmartPointer<vtkDataArray> outputScalarField = createVtkDataArray(dataType);
          outputScalarField->SetNumberOfTuples(numberOfPoints);
          outputScalarField->SetName("L2-distance");
          return outputScalarField;
        }
        
        template <class dataType>
        dataType computeL2Distance(vtkSmartPointer<vtkDataSet> img1, vtkSmartPointer<vtkDataSet> img2,
                                   bool emptyFieldDistance=false){
          vtkDataArray *inputScalarField1 = img1->GetPointData()->GetArray("Scalars");
          vtkDataArray *inputScalarField2 = img2->GetPointData()->GetArray("Scalars");
          auto noPoints = inputScalarField1->GetNumberOfTuples();
          auto outputScalarField = createL2OutputScalarField(inputScalarField1->GetDataType(), noPoints);

          void *inputDataPointer2 = inputScalarField2->GetVoidPointer(0);
          if(emptyFieldDistance){
            auto temp = createVtkDataArray(inputScalarField1->GetDataType());
            temp->DeepCopy(inputScalarField1);
            for(int i = 0; i < inputScalarField1->GetNumberOfTuples(); ++i)
              temp->SetTuple1(i, 0);
            inputDataPointer2 = temp->GetVoidPointer(0);
          }
          
          LDistance lDistance;
          lDistance.setNumberOfPoints(noPoints);          
          lDistance.setOutputDataPointer(outputScalarField->GetVoidPointer(0));
          lDistance.setInputDataPointer1(inputScalarField1->GetVoidPointer(0));
          lDistance.setInputDataPointer2(inputDataPointer2);
          int errorCode;
          switch(inputScalarField1->GetDataType()) {
            vtkTemplateMacro(errorCode = lDistance.execute<VTK_TT>("2"));
          }
          if(errorCode != 0)
            std::cout << "LDistance errorCode=" << errorCode << std::endl;

          auto distance = lDistance.getResult();

          return distance;
        }
        
        template <class dataType>
        vtkSmartPointer<vtkDataSet> computeL2Barycenter(vtkSmartPointer<vtkDataSet> img1, 
                                                          vtkSmartPointer<vtkDataSet> img2, double alpha){
          // Rename and share scalars among inputs
          img1->GetPointData()->GetArray("Scalars")->SetName("Scalars1");
          img2->GetPointData()->GetArray("Scalars")->SetName("Scalars2");
          img1->GetPointData()->AddArray(img2->GetPointData()->GetArray("Scalars2"));
          img2->GetPointData()->AddArray(img1->GetPointData()->GetArray("Scalars1"));
          
          // Append inputs
          vtkSmartPointer<vtkAppendFilter> appendFilter = vtkSmartPointer<vtkAppendFilter>::New();
          appendFilter->AddInputData(img1);
          appendFilter->AddInputData(img2);
          appendFilter->MergePointsOn();
          appendFilter->Update();
          auto appendFilterOut = appendFilter->GetOutput();

          // Compute interpolation
          vtkSmartPointer<vtkArrayCalculator> arrayCalculator = vtkSmartPointer<vtkArrayCalculator>::New();
          arrayCalculator->SetInputData(appendFilterOut);
          arrayCalculator->AddScalarArrayName("Scalars1");
          arrayCalculator->AddScalarArrayName("Scalars2");
          std::stringstream ss;
          ss << alpha << "*Scalars1 + " << (1-alpha) << "*Scalars2";
          arrayCalculator->SetFunction(ss.str().data());
          arrayCalculator->SetResultArrayName("Scalars");
          arrayCalculator->Update();
          auto arrayCalculatorOut = arrayCalculator->GetOutput();
          
          // Rename and remove shared scalars
          img1->GetPointData()->GetArray("Scalars1")->SetName("Scalars");
          img2->GetPointData()->GetArray("Scalars2")->SetName("Scalars");
          img1->GetPointData()->RemoveArray("Scalars2");
          img2->GetPointData()->RemoveArray("Scalars1");
          
          // Create output
          vtkSmartPointer<vtkDataSet> barycenter;
          if(DoResampleToImage)
            barycenter = vtkSmartPointer<vtkImageData>::New();
          else
            barycenter = vtkSmartPointer<vtkUnstructuredGrid>::New();
          barycenter->DeepCopy(arrayCalculatorOut);
          barycenter->GetPointData()->RemoveArray("Scalars2");
          barycenter->GetPointData()->RemoveArray("Scalars1");
          
          return barycenter;
        }
        
        template <class dataType>
        dataType computeDistance(MergeTree *mTree1, MergeTree *mTree2, bool emptyTreeDistance=false){
          MergeTreeDistance ftmTreeEditDistance;
          ftmTreeEditDistance.setAssignmentSolver(AssignmentSolver);
          ftmTreeEditDistance.setEpsilonTree1(EpsilonTree1);
          ftmTreeEditDistance.setEpsilonTree2(EpsilonTree2);
          ftmTreeEditDistance.setEpsilon2Tree1(Epsilon2Tree1);
          ftmTreeEditDistance.setEpsilon2Tree2(Epsilon2Tree2);
          ftmTreeEditDistance.setEpsilon3Tree1(Epsilon3Tree1);
          ftmTreeEditDistance.setEpsilon3Tree2(Epsilon3Tree2);
          ftmTreeEditDistance.setProgressiveComputation(ProgressiveComputation);
          ftmTreeEditDistance.setBranchDecomposition(BranchDecomposition);
          ftmTreeEditDistance.setParallelize(Parallelize);
          ftmTreeEditDistance.setPersistenceThreshold(PersistenceThreshold);
          ftmTreeEditDistance.setNormalizedWasserstein(NormalizedWasserstein);
          ftmTreeEditDistance.setNormalizedWassersteinReg(NormalizedWassersteinReg);
          ftmTreeEditDistance.setRescaledWasserstein(RescaledWasserstein);
          ftmTreeEditDistance.setKeepSubtree(KeepSubtree);
          ftmTreeEditDistance.setUseMinMaxPair(UseMinMaxPair);
          ftmTreeEditDistance.setNumberOfThreads(NumberOfThreads);
          ftmTreeEditDistance.setDistanceSquared(true); // squared root
          ftmTreeEditDistance.setVerbose(0);
          ftmTreeEditDistance.setPreprocess(false);
          ftmTreeEditDistance.setPostprocess(false);
          //ftmTreeEditDistance.setIsCalled(true);
          ftmTreeEditDistance.setOnlyEmptyTreeDistance(emptyTreeDistance);
          
          std::vector<std::tuple<ftm::idNode, ftm::idNode>> matching;
          dataType distance = ftmTreeEditDistance.execute<dataType>(mTree1->tree, mTree2->tree, matching);
          
          return distance;
        }
        
        template<class dataType> 
        MergeTree* computeBarycenter(MergeTree *mTree1, MergeTree *mTree2, double alpha){
          MergeTreeBarycenter ftmTreeEditDistanceBary;
          ftmTreeEditDistanceBary.setAssignmentSolver(AssignmentSolver);
          ftmTreeEditDistanceBary.setEpsilonTree1(EpsilonTree1);
          ftmTreeEditDistanceBary.setEpsilonTree2(EpsilonTree2);
          ftmTreeEditDistanceBary.setEpsilon2Tree1(Epsilon2Tree1);
          ftmTreeEditDistanceBary.setEpsilon2Tree2(Epsilon2Tree2);
          ftmTreeEditDistanceBary.setEpsilon3Tree1(Epsilon3Tree1);
          ftmTreeEditDistanceBary.setEpsilon3Tree2(Epsilon3Tree2);
          ftmTreeEditDistanceBary.setProgressiveComputation(ProgressiveComputation);
          ftmTreeEditDistanceBary.setBranchDecomposition(BranchDecomposition);
          ftmTreeEditDistanceBary.setParallelize(Parallelize);
          ftmTreeEditDistanceBary.setPersistenceThreshold(PersistenceThreshold);
          ftmTreeEditDistanceBary.setNormalizedWasserstein(NormalizedWasserstein);
          ftmTreeEditDistanceBary.setNormalizedWassersteinReg(NormalizedWassersteinReg);
          ftmTreeEditDistanceBary.setRescaledWasserstein(RescaledWasserstein);
          ftmTreeEditDistanceBary.setKeepSubtree(KeepSubtree);
          ftmTreeEditDistanceBary.setUseMinMaxPair(UseMinMaxPair);
          ftmTreeEditDistanceBary.setNumberOfThreads(NumberOfThreads);
          ftmTreeEditDistanceBary.setAlpha(alpha);
          ftmTreeEditDistanceBary.setVerbose(0);
          ftmTreeEditDistanceBary.setPreprocess(false);
          ftmTreeEditDistanceBary.setPostprocess(false);
          //ftmTreeEditDistanceBary.setIsCalled(true);
          
          std::vector<FTMTree_MT *> intermediateTrees;
          intermediateTrees.push_back(mTree1->tree);
          intermediateTrees.push_back(mTree2->tree);
          std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>> outputMatchingBarycenter(2);
          MergeTree * barycenter = ftmTreeEditDistanceBary.execute<dataType>(intermediateTrees, 
                                                                             outputMatchingBarycenter);
          return barycenter;
        }
        
        double computeAlpha(int index1, int middleIndex, int index2){
          return 1-((double)middleIndex - index1) / (index2 - index1);
        }
        
        template<class dataType>
        void temporalSubsampling(std::vector<MergeTree*> &mTrees, std::vector<int> &removed, 
                                 std::vector<MergeTree*> &barycenters, 
                                 std::vector<vtkSmartPointer<vtkDataSet>> &barycentersL2){    
          int verbose = 0;
          
          std::vector<bool> treeRemoved(mTrees.size(), false);
          
          int toRemoved = mTrees.size() * RemovalPercentage / 100.;
          toRemoved = std::min(toRemoved, (int)(mTrees.size()-3));
          //std::cout << "will removed " << toRemoved << " of " << mTrees.size() << std::endl;
          for(int iter = 0; iter < toRemoved; ++iter){
            dataType bestCost = std::numeric_limits<dataType>::max();
            int bestMiddleIndex = -1, bestIndex1 = -1, bestIndex2 = -1;
            MergeTree* bestBarycenter;
            std::vector<std::tuple<MergeTree*, int>> bestBarycentersOnPath;
            vtkSmartPointer<vtkDataSet> bestBarycenterL2;
            std::vector<std::tuple<vtkSmartPointer<vtkDataSet>, int>> bestBarycentersL2OnPath;
            
            // Compute barycenter for each pair of trees
            if(verbose > 0) std::cout << "Compute barycenter for each pair of trees" << std::endl;
            int index1 = 0, index2 = 0;
            while(index2 != mTrees.size()-1){
              
              // Get index in the middle
              int middleIndex = index1+1;
              while(treeRemoved[middleIndex])
                ++middleIndex;
              
              // Get second index
              index2 = middleIndex+1;
              while(treeRemoved[index2])
                ++index2;
              
              // Compute barycenter
              if(verbose > 0) std::cout << "Compute barycenter" << std::endl;
              //std::cout << "compute barycenter of " << index1 << " and " << index2 << std::endl;
              double alpha = computeAlpha(index1, middleIndex, index2);
              MergeTree* barycenter; 
              vtkSmartPointer<vtkDataSet> barycenterL2;
              if(not UseL2Distance)
                barycenter = computeBarycenter<dataType>(mTrees[index1], mTrees[index2], alpha);
              else
                barycenterL2 = computeL2Barycenter<dataType>(images[index1], images[index2], alpha);
              
              // - Compute cost
              // Compute distance with middleIndex
              if(verbose > 0) std::cout << "Compute distance with middleIndex" << std::endl;
              dataType cost;
              if(not UseL2Distance)
                cost = computeDistance<dataType>(barycenter, mTrees[middleIndex]);
              else
                cost = computeL2Distance<dataType>(barycenterL2, images[middleIndex]);
              //std::cout << "distance barycenter and " << middleIndex << " is " << cost << std::endl;
              
              // Compute distances of previously removed trees on the path
              if(verbose > 0) std::cout << "Compute distances of previously removed trees" << std::endl;
              std::vector<std::tuple<MergeTree*, int>> barycentersOnPath;
              std::vector<std::tuple<vtkSmartPointer<vtkDataSet>, int>> barycentersL2OnPath;
              for(int i = 0; i < 2; ++i){
                int toReach = (i == 0 ? index1 : index2);
                int offset = (i == 0 ? -1 : 1);
                int tIndex = middleIndex + offset;
                while(tIndex != toReach){
                  
                  // Compute barycenter
                  double alpha = computeAlpha(index1, tIndex, index2);
                  MergeTree* barycenterP; 
                  vtkSmartPointer<vtkDataSet> barycenterPL2; 
                  if(not UseL2Distance)
                    barycenterP = computeBarycenter<dataType>(mTrees[index1], mTrees[index2], alpha);
                  else  
                    barycenterPL2 = computeL2Barycenter<dataType>(images[index1], images[index2], alpha);
                  
                  // Compute distance
                  dataType costP;
                  if(not UseL2Distance)
                    costP = computeDistance<dataType>(barycenterP, mTrees[tIndex]);
                  else
                    costP = computeL2Distance<dataType>(barycenterPL2, images[tIndex]);
                  
                  // Save results
                  if(not UseL2Distance)
                    barycentersOnPath.push_back(std::make_tuple(barycenterP, tIndex));
                  else
                    barycentersL2OnPath.push_back(std::make_tuple(barycenterPL2, tIndex));
                  cost += costP;
                  tIndex += offset;
                }
              }
              
              if(cost < bestCost){
                bestCost = cost;
                bestMiddleIndex = middleIndex;
                bestIndex1 = index1;
                bestIndex2 = index2;
                if(not UseL2Distance){
                  bestBarycenter = barycenter;
                  bestBarycentersOnPath = barycentersOnPath;
                }else{
                  bestBarycenterL2 = barycenterL2;
                  bestBarycentersL2OnPath = barycentersL2OnPath;
                }
              }
              
              // Go to the next index
              index1 = middleIndex;
            }
            
            // Removed the tree with the lowest cost
            if(verbose > 0) std::cout << "Removed the tree with the lowest cost" << std::endl;
            //std::cout << "remove " << bestMiddleIndex << " with cost = " << bestCost << std::endl;
            removed.push_back(bestMiddleIndex);
            treeRemoved[bestMiddleIndex] = true;
            if(not UseL2Distance){
              barycenters[bestMiddleIndex] = bestBarycenter;
              for(auto tup : bestBarycentersOnPath)
                barycenters[std::get<1>(tup)] = std::get<0>(tup);
            }else{
              barycentersL2[bestMiddleIndex] = bestBarycenterL2;
              for(auto tup : bestBarycentersL2OnPath)
                barycentersL2[std::get<1>(tup)] = std::get<0>(tup);
            }
          }
        }
        
        template <class dataType>
        std::vector<int> execute(std::vector<MergeTree*> &mTrees, 
                                 std::vector<std::vector<double>> &distanceMatrix, 
                                 std::vector<double> &emptyTreeDistances, std::vector<MergeTree*> &allMT){
          Timer t_tempSub;
          
          // --- Preprocessing
          if(UseL2Distance){
            if(DoResampleToImage) 
              std::cout << "resample to image..." << std::endl;
            preprocessSegmentation();
            if(DoResampleToImage) 
              std::cout << "resample to image done." << std::endl;
          }else{
            treesNodeCorr = std::vector<std::vector<int>>(mTrees.size());
            for(int i = 0; i < mTrees.size(); ++i){
              //std::cout << "preprocess tree " << i << " / " << mTrees.size() << std::endl;
              preprocessingPipeline<dataType>(mTrees[i]->tree, EpsilonTree2, Epsilon2Tree2, Epsilon3Tree2, 
                                              BranchDecomposition, UseMinMaxPair, CleanTree, 
                                              treesNodeCorr[i], 0);
              printTreeStats(mTrees[i]->tree);
            }
          }
          
          // --- Execute
          std::vector<MergeTree*> barycenters(mTrees.size());
          std::vector<vtkSmartPointer<vtkDataSet>> barycentersL2(mTrees.size());
          std::vector<int> removed;
          temporalSubsampling<dataType>(mTrees, removed, barycenters, barycentersL2);
          
          // --- Concatenate all trees/L2Images
          for(auto mt : mTrees)
            allMT.push_back(mt);
          std::vector<bool> removedB(mTrees.size(), false);
          for(auto r : removed)
            removedB[r] = true;
          for(int i = 0; i < barycenters.size(); ++i)
            if(removedB[i]){
              if(not UseL2Distance)
                allMT.push_back(barycenters[i]);
              else
                images.push_back(barycentersL2[i]);
            }
          
          // --- Compute empty tree distances
          int distMatSize = (not UseL2Distance ? allMT.size() : images.size());
          for(int i = 0; i < distMatSize; ++i){
            dataType distance;
            if(not UseL2Distance)
              distance = computeDistance<dataType>(allMT[i], allMT[i], true);
            else
              distance = computeL2Distance<dataType>(images[i], images[i], true);
            emptyTreeDistances.push_back(distance);
          }
          
          // --- Create distance matrix
          distanceMatrix = std::vector<std::vector<double>>(distMatSize, 
                                                            std::vector<double>(distMatSize, 0));
/*#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(NumberOfThreads) if(Parallelize)
#endif*/
          for(int i = 0; i < distMatSize; ++i)
            for(int j = i+1; j < distMatSize; ++j){
              std::vector<std::tuple<ftm::idNode, ftm::idNode>> matching;
              dataType distance;
              if(not UseL2Distance)
                distance = computeDistance<dataType>(allMT[i], allMT[j]);
              else
                distance = computeL2Distance<dataType>(images[i], images[j]);
              distanceMatrix[i][j] = distance;
              distanceMatrix[j][i] = distance;
            }
          
          // --- Postprocessing
          if(not UseL2Distance)
            for(int i = 0; i < allMT.size(); ++i)
              // TODO fix merged root origin when full merge
              postprocessingPipeline<dataType>(allMT[i]->tree);
          
          // --- Print results
          std::cout << "input   = " << mTrees.size() << std::endl;
          std::cout << "output  = " << distMatSize << std::endl;
          std::cout << "removed : ";
          for(int i = 0; i < removed.size(); ++i){
            auto r = removed[i];
            std::cout << r;
            if(i < removed.size()-1)
              std::cout << ", ";
          }
          std::cout << std::endl;
          
          sort(removed.begin(), removed.end());
          
          std::cout << "TIME TEMP.SUB   = " << t_tempSub.getElapsedTime() << std::endl;
          
          return removed;
        }
        
        // --------------------------------------------------------------------------------
        // Utils
        // --------------------------------------------------------------------------------
        void printVtkDataArray(vtkSmartPointer<vtkDataArray> array){
          for(int i = 0; i < array->GetNumberOfTuples(); ++i)
            std::cout << array->GetTuple1(i) << ", ";
          std::cout << std::endl;
        }
        
    }; // MergeTreeTemporalSubsampling class
    
} // namespace ttk

#endif
