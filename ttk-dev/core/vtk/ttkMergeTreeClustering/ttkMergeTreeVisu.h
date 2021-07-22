#ifndef _MERGETREEVISU_H
#define _MERGETREEVISU_H

#pragma once

#include <MergeTreeUtils.h>
#include <ttkMergeTreeUtilsVisu.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <ttkTriangulation.h>

#include <vtkFloatArray.h>
#include <vtkAppendFilter.h>
#include <vtkCellData.h>
#include <vtkPointData.h>

class ttkMergeTreeVisu {
private:  
  bool PlanarLayout = false;
  bool BranchDecompositionPlanarLayout = false;
  double BranchSpacing = 1.;
  bool RescaleTreesIndividually = false;
  bool OutputSegmentation = false;
  double DimensionSpacing = 1.;
  int DimensionToShift = 0;
  double ImportantPairs = 50.; // important pairs threshold
  double ImportantPairsSpacing = 1.;
  double NonImportantPairsSpacing = 1.;
  double NonImportantPairsProximity = 0.05;
  
  int ShiftMode = 0; // 0: Star ; 1: Barycenter ; 2: Line ; 3: Double Line
  
  int iSampleOffset = 0;
  int noSampleOffset = 0;
  
  std::vector<std::vector<SimplexId>> nodeCorr;
  std::vector<double> clusterShift;
  
public:
  ttkMergeTreeVisu(){};
  ~ttkMergeTreeVisu(){};
  
  // ==========================================================================
  // Getter / Setter
  // ==========================================================================
  void setPlanarLayout(bool b){
    PlanarLayout = b;
  }
  void setBranchDecompositionPlanarLayout(bool b){
    BranchDecompositionPlanarLayout = b;
  }
  void setBranchSpacing(double d){
    BranchSpacing = d;
  }
  void setRescaleTreesIndividually(bool b){
    RescaleTreesIndividually = b;
  }
  void setOutputSegmentation(bool b){
    OutputSegmentation = b;
  }
  void setDimensionSpacing(double d){
    DimensionSpacing = d;
  }
  void setDimensionToShift(int i){
    DimensionToShift = i;
  }
  void setImportantPairs(double d){
    ImportantPairs = d;
  }
  void setImportantPairsSpacing(double d){
    ImportantPairsSpacing = d;
  }
  void setNonImportantPairsSpacing(double d){
    NonImportantPairsSpacing = d;
  }
  void setNonImportantPairsProximity(double d){
    NonImportantPairsProximity = d;
  }
  
  void setShiftMode(int i){
    ShiftMode = i;
  }
  
  void setISampleOffset(int i){
    iSampleOffset = i;
  }
  void setNoSampleOffset(int i){
    noSampleOffset = i;
  }
  
  std::vector<std::vector<SimplexId>> getNodeCorr(){
    return nodeCorr;
  }
  
  std::vector<double> getClusterShift(){
    return clusterShift;
  }
  
  // ==========================================================================
  // Visu Functions
  // ==========================================================================
  void makeTreesOutput(FTMTree_MT *tree1, 
       vtkUnstructuredGrid *treesNodes, std::vector<int> &trees1NodeCorrMesh,
       vtkDataSet *treesSegmentation,
       std::tuple<double,double,double,double,double,double> allBounds,
       vtkSmartPointer<vtkUnstructuredGrid> vtkOutputNode1, 
       vtkSmartPointer<vtkUnstructuredGrid> vtkOutputArc1,
       vtkSmartPointer<vtkUnstructuredGrid> vtkOutputSegmentation1,
       int dataType, int verbose=0){
    std::vector<FTMTree_MT *> intermediateTrees{tree1};
    std::vector<vtkUnstructuredGrid *> treesNodesT{treesNodes};
    std::vector<std::vector<int>> trees1NodeCorrMeshT{trees1NodeCorrMesh};
    std::vector<vtkDataSet *> treesSegmentationT{treesSegmentation};
    std::vector<std::tuple<double,double,double,double,double,double>> allBoundsT{allBounds};
    
    makeTreesOutput(intermediateTrees, treesNodesT, trees1NodeCorrMeshT, treesSegmentationT, allBoundsT, 
                    vtkOutputNode1, vtkOutputArc1, vtkOutputSegmentation1, dataType, verbose);
  }
  
  void makeTreesOutput(FTMTree_MT *tree1, FTMTree_MT *tree2,
       std::vector<vtkUnstructuredGrid *> &treesNodes, std::vector<std::vector<int>> &trees1NodeCorrMesh,
       std::vector<vtkDataSet *> &treesSegmentation,
       std::vector<std::tuple<double,double,double,double,double,double>> &allBounds,
       vtkSmartPointer<vtkUnstructuredGrid> vtkOutputNode1, 
       vtkSmartPointer<vtkUnstructuredGrid> vtkOutputArc1,
       vtkSmartPointer<vtkUnstructuredGrid> vtkOutputSegmentation1,
       int dataType, int verbose=0){
    std::vector<FTMTree_MT *> intermediateTrees{tree1, tree2};
    
    makeTreesOutput(intermediateTrees, treesNodes, trees1NodeCorrMesh, treesSegmentation, allBounds, 
                    vtkOutputNode1, vtkOutputArc1, vtkOutputSegmentation1, dataType, verbose);
  }
  
  void makeTreesOutput(std::vector<FTMTree_MT *> &intermediateTrees,
       std::vector<vtkUnstructuredGrid *> &treesNodes, std::vector<std::vector<int>> &trees1NodeCorrMesh,
       std::vector<vtkDataSet *> &treesSegmentation,
       std::vector<std::tuple<double,double,double,double,double,double>> &allBounds,
       vtkSmartPointer<vtkUnstructuredGrid> vtkOutputNode1, 
       vtkSmartPointer<vtkUnstructuredGrid> vtkOutputArc1,
       vtkSmartPointer<vtkUnstructuredGrid> vtkOutputSegmentation1,
       int dataType, int verbose=0){
    std::vector<MergeTree*> barycenters;
    std::vector<int> clusteringAssignment(intermediateTrees.size(), 0);
    std::vector<std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>>> outputMatchingBarycenter;
    std::vector<std::tuple<double,double,double,double,double,double>> allBaryBounds;
    allBaryBounds.push_back(getMaximalBounds(allBounds, clusteringAssignment, 0));
    std::vector<std::vector<int>> allBaryBranchingID;
    
    makeTreesOutput(intermediateTrees, barycenters, clusteringAssignment, outputMatchingBarycenter,
                    treesNodes, trees1NodeCorrMesh, treesSegmentation, allBounds, 
                    allBaryBounds, allBaryBranchingID, vtkOutputNode1, vtkOutputArc1, 
                    vtkOutputSegmentation1, dataType, verbose);
  }
  
  void makeTreesOutput(std::vector<FTMTree_MT *> &trees, std::vector<MergeTree*> &barycenters, 
      std::vector<int> &clusteringAssignment,
      std::vector<std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>>> &outputMatchingBarycenter,
      std::vector<vtkUnstructuredGrid *> &treesNodes, std::vector<std::vector<int>> &trees1NodeCorrMesh,
      std::vector<vtkDataSet *> &treesSegmentation, 
      std::vector<std::tuple<double,double,double,double,double,double>> &allBounds,
      std::vector<std::tuple<double,double,double,double,double,double>> &allBaryBounds,
      std::vector<std::vector<int>> &allBaryBranchingID,
      vtkSmartPointer<vtkUnstructuredGrid> vtkOutputNode1, 
      vtkSmartPointer<vtkUnstructuredGrid> vtkOutputArc1,
      vtkSmartPointer<vtkUnstructuredGrid> vtkOutputSegmentation1,
      int dataType, int verbose=0){
    //verbose = 5;
    
    int numInputs = trees.size();
    int NumberOfBarycenters = barycenters.size();
    bool barycenterOutput = (NumberOfBarycenters != 0);
    NumberOfBarycenters = std::max(NumberOfBarycenters, 1); // to enter the outer loop if not barycenterOutput
    
    nodeCorr = std::vector<std::vector<SimplexId>>(numInputs);
    clusterShift = std::vector<double>(NumberOfBarycenters, 0);
    
    if(verbose > 0) std::cout << "// --- Input Trees" << std::endl;
    // - Declare VTK arrays
    vtkSmartPointer<vtkUnstructuredGrid> vtkArcs = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    
    // Node fields
    vtkNew<vtkIntArray> criticalType{};
    criticalType->SetName("CriticalType");
    vtkNew<vtkFloatArray> persistenceNode{};
    persistenceNode->SetName("Persistence");
    vtkNew<vtkIntArray> clusterIDNode{};
    clusterIDNode->SetName("ClusterID");
    vtkNew<vtkIntArray> treeIDNode{};
    treeIDNode->SetName("TreeID");
    vtkNew<vtkIntArray> isDummyNode{};
    isDummyNode->SetName("isDummyNode");
    vtkNew<vtkIntArray> branchNodeID{};
    branchNodeID->SetName("BranchNodeID");
    vtkNew<vtkFloatArray> scalar{};
    scalar->SetName("Scalar");
    vtkNew<vtkIntArray> branchBaryNodeID{};
    branchBaryNodeID->SetName("BranchBaryNodeID");
    vtkNew<vtkIntArray> isInterpolatedTreeNode{};
    isInterpolatedTreeNode->SetName("isInterpolatedTree");
    vtkNew<vtkIntArray> isImportantPairsNode{};
    isImportantPairsNode->SetName("isImportantPair");
    
    vtkNew<vtkIntArray> nodeID{};
    nodeID->SetName("NodeId");
    
    // Arc fields
    vtkNew<vtkFloatArray> persistenceArc{};
    persistenceArc->SetName("Persistence");
    vtkNew<vtkIntArray> clusterIDArc{};
    clusterIDArc->SetName("ClusterID");
    vtkNew<vtkIntArray> treeIDArc{};
    treeIDArc->SetName("TreeID");
    vtkNew<vtkIntArray> isImportantPairsArc{};
    isImportantPairsArc->SetName("isImportantPair");
    vtkNew<vtkIntArray> isDummyArc{};
    isDummyArc->SetName("isDummyArc");
    vtkNew<vtkIntArray> branchID{};
    branchID->SetName("BranchID");
    vtkNew<vtkIntArray> branchBaryID{};
    branchBaryID->SetName("BranchBaryID");
    vtkNew<vtkIntArray> isInterpolatedTreeArc{};
    isInterpolatedTreeArc->SetName("isInterpolatedTree");
    
    vtkNew<vtkIntArray> upNodeId{};
    upNodeId->SetName("upNodeId");
    vtkNew<vtkIntArray> downNodeId{};
    downNodeId->SetName("downNodeId");
    
    // Segmentation
    //vtkSmartPointer<vtkPoints> segPoints = vtkSmartPointer<vtkPoints>::New();
    //vtkSmartPointer<vtkUnstructuredGrid> segCells= vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkAppendFilter> appendFilter = vtkSmartPointer<vtkAppendFilter>::New();

    int cellCount = 0;
    bool foundOneInterpolatedTree = false;

    // Create trees output
    if(verbose > 0) std::cout << "// Create trees output" << std::endl;
    for(int c = 0; c < NumberOfBarycenters; ++c){
      double delta_max = std::numeric_limits<double>::lowest();
      int noSample = 0 + noSampleOffset;
      for(int i = 0; i < numInputs; ++i){
        delta_max = std::max((std::get<3>(allBounds[i]) - std::get<2>(allBounds[i])), delta_max);
        delta_max = std::max((std::get<1>(allBounds[i]) - std::get<0>(allBounds[i])), delta_max);
        if(clusteringAssignment[i] != c)
          continue;
        noSample += 1;
      }
      double radius = delta_max*2 * DimensionSpacing;
      double pi = 3.14159265359;
      int iSample = 0 + iSampleOffset;
      double prevXMax = 0;
      std::vector<double> allPrevXMax;
      double prevYMax = std::numeric_limits<double>::lowest();
      
      // Iterate through all trees
      for(int i = 0; i < numInputs; ++i){
        if(clusteringAssignment[i] != c)
          continue;
        bool isInterpolatedTree = !treesNodes[i] or trees1NodeCorrMesh[i].size() == 0 or 
                                  !treesSegmentation[i];
        foundOneInterpolatedTree |= isInterpolatedTree;
        std::vector<ftm::idNode> treeBranching;
        std::vector<int> treeBranchingID;
        getTreeBranching(trees[i], treeBranching, treeBranchingID);
        
        nodeCorr[i] = std::vector<SimplexId>(trees[i]->getNumberOfNodes());
        double angle = 360.0 / noSample * iSample;
        iSample += 1;
        double diff_x=0, diff_y=0;
        switch(ShiftMode){
          case 0: // Star
            diff_x = radius * std::cos(-1*angle * pi / 180) + clusterShift[c];
            diff_y = radius * std::sin(-1*angle * pi / 180);
            break;
          case 1: // Barycenter TODO
            break;
          case 2: // Line
            diff_x = prevXMax + radius;
            break;
          case 3: // Double Line 
            diff_x = prevXMax + radius;
            if(i >= numInputs / 2){
              diff_y = - (prevYMax + radius/2);
              diff_x = allPrevXMax[i - int(numInputs / 2)] + radius;
            }else
              allPrevXMax.push_back(prevXMax);
            break;
          default:
            break;
        }
        double diff_z = - std::get<4>(allBounds[i]);
        if(/*not barycenterOutput and*/ not PlanarLayout){
          if(DimensionToShift != 0){ // is not X
            if(DimensionToShift == 2) // is Z
              diff_z = diff_x;
            else if(DimensionToShift == 1) // is Y
              diff_y = diff_x;
            diff_x = - std::get<0>(allBounds[i]);
          }
        }
        
        std::vector<float> layout;
        if(PlanarLayout){
          double refPersistence;
          if(barycenterOutput)
            refPersistence = getPersistence(barycenters[0]->tree, getRoot(barycenters[0]->tree), dataType);
          else
            refPersistence = getPersistence(trees[0], getRoot(trees[0]), dataType);
          layout = treePlanarLayout(trees[i], dataType, allBaryBounds[c], refPersistence);
        }
        int cptNode = 0;
        std::vector<SimplexId> treeSimplexId(trees[i]->getNumberOfNodes());
        std::vector<SimplexId> treeDummySimplexId(trees[i]->getNumberOfNodes());
        std::vector<SimplexId> layoutCorr(trees[i]->getNumberOfNodes());
        std::vector<ftm::idNode> treeMatching(trees[i]->getNumberOfNodes(), -1);
        if(barycenterOutput)
          for(auto match : outputMatchingBarycenter[c][i])
            treeMatching[std::get<1>(match)] = std::get<0>(match);
        
        // Tree traversal
        if(verbose > 0) std::cout << "// Tree traversal" << std::endl;
        std::queue<idNode> queue;
        queue.push(getRoot(trees[i]));
        while(!queue.empty()){
          idNode node = queue.front();
          queue.pop();
          idNode nodeOrigin = trees[i]->getNode(node)->getOrigin();
          
          // Push children to the queue
          if(verbose > 2) std::cout << "// Push children to the queue" << std::endl;
          for(auto child : getChildren(trees[i], node))
            queue.push(child);
          
          // Get and insert point
          if(verbose > 2) std::cout << "// Get and insert point" << std::endl;
          int nodeMesh;
          double point[3] = {0, 0, 0};
          if(not isInterpolatedTree){
            nodeMesh = trees1NodeCorrMesh[i][node];
            double* pointP = treesNodes[i]->GetPoints()->GetPoint(nodeMesh);
            for(int i = 0; i < 3; ++i)
              point[i] = pointP[i];
          }
          if(PlanarLayout){
            layoutCorr[node] = cptNode;
            point[0] = layout[cptNode];
            point[1] = layout[cptNode+1];
            point[2] = 0;
            cptNode+=2;
          }else
            point[2] += diff_z;
          point[0] += diff_x;
          point[1] += diff_y;
          
          // Get x Max and y Min for next iteration if needed
          prevXMax = std::max(prevXMax, point[0]);
          if(ShiftMode == 3){ // Double line
            if(i < numInputs / 2) 
              prevYMax = std::max(prevYMax, point[1]);
            if(i == int(numInputs / 2) - 1) 
              prevXMax = 0;
          }

          // TODO to many dummy nodes are created
          bool dummyNode = PlanarLayout and not BranchDecompositionPlanarLayout
                                        and !isRoot(trees[i], node)
                                        /*and !isLeaf(trees[i], node) 
                                        and isBranchOrigin(trees[i], node)*/;
          if(dummyNode)
            treeDummySimplexId[node] = points->InsertNextPoint(point); // will be modified when processing son
          SimplexId nextPointId = points->InsertNextPoint(point);
          treeSimplexId[node] = nextPointId;
          nodeCorr[i][node] = nextPointId;
          if(dummyNode)
            nodeCorr[i][node] = treeDummySimplexId[node];
          
          // Insert cell connecting parent
          if(verbose > 2) std::cout << "// Add cell connecting parent" << std::endl;
          if(!isRoot(trees[i], node)){
            vtkIdType pointIds[2];
            pointIds[0] = treeSimplexId[node];
            
            ftm::idNode nodeParent = getParent(trees[i], node);
            // TODO to many dummy cells are created
            bool dummyCell = PlanarLayout and not BranchDecompositionPlanarLayout
                                          and treeBranching[node] == nodeParent 
                                          and !isRoot(trees[i], nodeParent);
            if(PlanarLayout and BranchDecompositionPlanarLayout){
              pointIds[1] = treeSimplexId[treeBranching[node]];
            }else if(dummyCell){
              double dummyPoint[3] = {point[0], layout[layoutCorr[nodeParent]+1]+diff_y, 0.};
              SimplexId dummyPointId = treeDummySimplexId[nodeParent];
              points->SetPoint(dummyPointId, dummyPoint);
              vtkIdType dummyPointIds[2];
              dummyPointIds[0] = dummyPointId;
              dummyPointIds[1] = treeSimplexId[nodeParent];
              vtkArcs->InsertNextCell(VTK_LINE, 2, dummyPointIds);
              pointIds[1] = dummyPointId;
            }else
              pointIds[1] = treeSimplexId[nodeParent];
            
            vtkArcs->InsertNextCell(VTK_LINE, 2, pointIds);
            
            // --------------
            // Arc field
            // --------------
            int toAdd = (dummyCell ? 2 : 1);
            for(int toAddT = 0; toAddT < toAdd; ++toAddT){              
              // Add branch bary ID
              if(verbose > 2) std::cout << "// Push arc bary branch id" << std::endl;
              if(barycenterOutput){
                int tBranchID = -1;
                if(treeMatching[node] >= 0 and treeMatching[node] < allBaryBranchingID[c].size()){
                  tBranchID = allBaryBranchingID[c][treeMatching[node]];
                  if(not isLeaf(trees[i], node) and treeMatching[nodeOrigin] >= 0 and 
                    treeMatching[nodeOrigin] < allBaryBranchingID[c].size())
                    tBranchID = allBaryBranchingID[c][treeMatching[nodeOrigin]];
                }
                branchBaryID->InsertNextTuple1(tBranchID);
              }
              
              // Add branch ID
              int tBranchID = treeBranchingID[node];
              branchID->InsertNextTuple1(tBranchID);
              
              // Add up and down nodeId
              upNodeId->InsertNextTuple1(treeSimplexId[nodeParent]);
              downNodeId->InsertNextTuple1(treeSimplexId[node]);
              
              // Add arc persistence
              if(verbose > 2) std::cout << "// Push arc persistence" << std::endl;
              idNode nodeToGetPers = treeBranching[node];
              if(PlanarLayout and BranchDecompositionPlanarLayout)
                nodeToGetPers = node;
              double persToAdd = getPersistence(trees[i], nodeToGetPers , dataType);
              persistenceArc->InsertNextTuple1(persToAdd);
              
              // Add arc cluster ID
              clusterIDArc->InsertNextTuple1(clusteringAssignment[i]);
              
              // Add arc tree ID
              treeIDArc->InsertNextTuple1(i);
              
              // Add isImportantPair
              bool isImportant;
              idNode nodeToGetImportance = treeBranching[node];
              if(PlanarLayout and BranchDecompositionPlanarLayout)
                nodeToGetImportance = node;
              switch(dataType){
                vtkTemplateMacro(isImportant = isImportantPair<VTK_TT>(trees[i], 
                                                                      nodeToGetImportance, ImportantPairs));
              }
              isImportantPairsArc->InsertNextTuple1(isImportant);
              
              // Add isDummyArc
              bool isDummy = toAdd == 2 and toAddT == 0;
              isDummyArc->InsertNextTuple1(isDummy);
              
              // Add isInterpolatedTree
              isInterpolatedTreeArc->InsertNextTuple1(isInterpolatedTree);
              
              cellCount++;
            }
          }
          
          // --------------
          // Node field
          // --------------
          int toAdd = (dummyNode ? 2 : 1);
          for(int toAddT = 0; toAddT < toAdd; ++toAddT){
            // Add node id
            nodeID->InsertNextTuple1(treeSimplexId[node]);
            
            // Add node scalar
            scalar->InsertNextTuple1(getValue(trees[i], node, dataType));
          
            // Add criticalType
            if(verbose > 2) std::cout << "// Add criticalType" << std::endl;
            int criticalTypeT = 0;
            if(not isInterpolatedTree)
              criticalTypeT = treesNodes[i]->GetPointData()->GetArray("CriticalType")->GetTuple1(nodeMesh);
            else{
              
            }
            criticalType->InsertNextTuple1(criticalTypeT);
            
            // Add node branch bary id
            if(verbose > 2) std::cout << "// Add node bary branch id" << std::endl;
            if(barycenterOutput){
              int tBranchID = -1;
              if(treeMatching[node] >= 0 and treeMatching[node] < allBaryBranchingID[c].size()){
                tBranchID = allBaryBranchingID[c][treeMatching[node]];
                if(not isLeaf(trees[i], node) and treeMatching[nodeOrigin] >= 0 and 
                  treeMatching[nodeOrigin] < allBaryBranchingID[c].size())
                  tBranchID = allBaryBranchingID[c][treeMatching[nodeOrigin]];
              }
              branchBaryNodeID->InsertNextTuple1(tBranchID);
            }
            
            // Add node branch bary id
            int tBranchID = treeBranchingID[node];
            if(not isLeaf(trees[i], node))
              tBranchID = treeBranchingID[nodeOrigin];
            branchNodeID->InsertNextTuple1(tBranchID);

            // Add node persistence
            if(verbose > 2) std::cout << "// Push node persistence" << std::endl;
            persistenceNode->InsertNextTuple1(getPersistence(trees[i], node, dataType));
            
            // Add node clusterID
            clusterIDNode->InsertNextTuple1(clusteringAssignment[i]);
            
            // Add arc tree ID
            treeIDNode->InsertNextTuple1(i);
            
            // Add isDummyNode
            bool isDummy = toAdd == 2 and toAddT == 1 and !isRoot(trees[i], node);
            isDummyNode->InsertNextTuple1(isDummy);
            
            // Add isInterpolatedTree
            isInterpolatedTreeNode->InsertNextTuple1(isInterpolatedTree);
            
            // Add isImportantPair
            bool isImportant;
            switch(dataType){
              vtkTemplateMacro(isImportant = isImportantPair<VTK_TT>(trees[i], node, ImportantPairs));
            }
            isImportantPairsNode->InsertNextTuple1(isImportant);
          }
          
          if(verbose > 2) std::cout << "end loop" << std::endl;
        }
        
        // Manage segmentation
        if(verbose > 1) std::cout << "// Shift segmentation" << std::endl;
        if(OutputSegmentation and not PlanarLayout){
          auto iTreesSegmentationCopy = vtkSmartPointer<vtkUnstructuredGrid>::New();
          iTreesSegmentationCopy->DeepCopy(treesSegmentation[i]);
          auto iVkOutputSegmentationTemp = vtkUnstructuredGrid::SafeDownCast(iTreesSegmentationCopy);
          for(int i = 0; i < iVkOutputSegmentationTemp->GetPoints()->GetNumberOfPoints(); ++i){
            double* point = iVkOutputSegmentationTemp->GetPoints()->GetPoint(i);
            point[0] += diff_x;
            point[1] += diff_y;
            point[2] += diff_z;
            iVkOutputSegmentationTemp->GetPoints()->SetPoint(i, point);
          }
          appendFilter->AddInputData(iVkOutputSegmentationTemp);
        }
      }
      if(c < NumberOfBarycenters-1)
        clusterShift[c+1] = radius*4 + clusterShift[c];
    }

    // --- Add VTK arrays to output
    if(verbose > 1) std::cout << "// Add VTK arrays to output" << std::endl;
    // Manage node output
    vtkOutputNode1->SetPoints(points);
    vtkOutputNode1->GetPointData()->AddArray(criticalType);
    vtkOutputNode1->GetPointData()->AddArray(persistenceNode);
    vtkOutputNode1->GetPointData()->AddArray(clusterIDNode);
    vtkOutputNode1->GetPointData()->AddArray(treeIDNode);
    vtkOutputNode1->GetPointData()->AddArray(isDummyNode);
    vtkOutputNode1->GetPointData()->AddArray(branchNodeID);
    vtkOutputNode1->GetPointData()->AddArray(nodeID);
    vtkOutputNode1->GetPointData()->AddArray(isImportantPairsNode);
    if(not BranchDecompositionPlanarLayout)
      vtkOutputNode1->GetPointData()->AddArray(scalar);
    if(barycenterOutput)
      vtkOutputNode1->GetPointData()->AddArray(branchBaryNodeID);
    if(foundOneInterpolatedTree)
      vtkOutputNode1->GetPointData()->AddArray(isInterpolatedTreeNode);
    
    // Manage arc output
    vtkArcs->SetPoints(points);
    vtkArcs->GetCellData()->AddArray(persistenceArc);
    vtkArcs->GetCellData()->AddArray(clusterIDArc);
    vtkArcs->GetCellData()->AddArray(treeIDArc);
    vtkArcs->GetCellData()->AddArray(isImportantPairsArc);
    vtkArcs->GetCellData()->AddArray(isDummyArc);
    vtkArcs->GetCellData()->AddArray(branchID);
    vtkArcs->GetCellData()->AddArray(upNodeId);
    vtkArcs->GetCellData()->AddArray(downNodeId);
    if(barycenterOutput)
      vtkArcs->GetCellData()->AddArray(branchBaryID);
    if(foundOneInterpolatedTree)
      vtkArcs->GetCellData()->AddArray(isInterpolatedTreeArc);
    vtkOutputArc1->ShallowCopy(vtkArcs);
    
    // Manage segmentation output
    if(OutputSegmentation and not PlanarLayout){
      appendFilter->Update();
      vtkOutputSegmentation1->ShallowCopy(appendFilter->GetOutput());
    }
  }
  
  // ==========================================================================
  // Planar Layout
  // ==========================================================================
  // TODO manage multi pers pairs
  template<class dataType>
  void treePlanarLayoutBDImpl(FTMTree_MT *tree, std::vector<float> &retVec, 
                              std::vector<LongSimplexId> &treeSimplexId, std::vector<idNode> &branching,
                              std::vector<std::vector<ftm::idNode>> &nodeBranching){
    idNode treeRoot = getRoot(tree);
    idNode treeRootOrigin = tree->getNode(treeRoot)->getOrigin();
    float rootY = retVec[treeSimplexId[treeRoot]*2+1];
    float rootOriginY = retVec[treeSimplexId[treeRootOrigin]*2+1];
    float rootYmin = std::min(rootY, rootOriginY);
    float rootYmax = std::max(rootY, rootOriginY);
    dataType rootPers = getNodePersistence<dataType>(tree, treeRoot);
    std::vector<std::tuple<float, float>> allNodeSpanX(tree->getNumberOfNodes());
    std::vector<std::tuple<float, float>> allNodeImportantSpanX(tree->getNumberOfNodes());
    
    // Compute gap
    float nonImportantPairsGap = (rootYmax - rootYmin) * 0.05 * NonImportantPairsSpacing;
    float importantPairsGap = std::max(nonImportantPairsGap, (float)ImportantPairsSpacing);
    
    // Some functions
    auto comp = [&](const ftm::idNode a, const ftm::idNode b) {
      return getNodePersistence<dataType>(tree, a) < getNodePersistence<dataType>(tree, b);
    };
    auto compSup = [&](const ftm::idNode a, const ftm::idNode b) { return not comp(a, b); };
    
    // Go
    std::vector<ftm::idNode> leaves = getLeaves(tree);
    sort(leaves.begin(), leaves.end(), comp);
    std::queue<ftm::idNode> queue;
    for(auto node : leaves)
      queue.push(node);
    while(!queue.empty()){
      idNode node = queue.front();
      queue.pop();
      idNode nodeOrigin = tree->getNode(node)->getOrigin();
      
      dataType nodePers = getNodePersistence<dataType>(tree, node);
      retVec[treeSimplexId[nodeOrigin]*2] = 0;
      retVec[treeSimplexId[nodeOrigin]*2+1] = nodePers / rootPers * (rootYmax - rootYmin) + rootYmin;
      
      // Positioning nodes in the branch
      if(isBranchOrigin(tree, nodeOrigin)){
        float prevX = 0;
        std::vector<idNode> nodeBranchingVector;
        for(int i = 1; i < nodeBranching[nodeOrigin].size(); ++i)
          nodeBranchingVector.push_back(nodeBranching[nodeOrigin][i]);
        sort(nodeBranchingVector.begin(), nodeBranchingVector.end(), compSup);
        
        // Iterate through each node of the branch
        int lastIndexImportant = -1;
        for(int i = 0; i < nodeBranchingVector.size(); ++i){
          idNode nodeBranchingI = nodeBranchingVector[i];
          
          // Get old node span X
          float oldMin = std::get<0>(allNodeSpanX[nodeBranchingI]);
          float oldMax = std::get<1>(allNodeSpanX[nodeBranchingI]);
          
          // Get x spacing
          float nodeSpacing = 0;
          if(i > 0){
            if(isImportantPair<dataType>(tree, nodeBranchingVector[i], ImportantPairs)){
              nodeSpacing = importantPairsGap;
            }else if(isImportantPair<dataType>(tree, nodeBranchingVector[i-1], ImportantPairs)){
              nodeSpacing = NonImportantPairsProximity;
              //prevX = 
            }else{
              nodeSpacing = nonImportantPairsGap;
            }
          }else if(not isImportantPair<dataType>(tree, nodeBranchingVector[i], ImportantPairs)
                  and isImportantPair<dataType>(tree, nodeOrigin, ImportantPairs))
            nodeSpacing = NonImportantPairsProximity;
          float newMin = prevX + nodeSpacing;
          float shiftX = newMin - oldMin;
          
          // Set y coordinate according difference in persistence
          dataType nodeBranchingIPers = getNodePersistence<dataType>(tree, nodeBranchingI);
          float shiftY = nodeBranchingIPers * BranchSpacing;
          float diffY = retVec[treeSimplexId[nodeBranchingI]*2+1];
          retVec[treeSimplexId[nodeBranchingI]*2+1] = retVec[treeSimplexId[nodeOrigin]*2+1] - shiftY;
          diffY = retVec[treeSimplexId[nodeBranchingI]*2+1] - diffY;
          
          // Shift this branch
          std::queue<idNode> queueBranching;
          queueBranching.push(nodeBranchingI);
          while(!queueBranching.empty()){
            idNode nodeBranchOrigin = queueBranching.front();
            queueBranching.pop();
            retVec[treeSimplexId[nodeBranchOrigin]*2] += shiftX;
            if(nodeBranchOrigin != nodeBranchingI)
              retVec[treeSimplexId[nodeBranchOrigin]*2+1] += diffY;
            if(isBranchOrigin(tree, nodeBranchOrigin))
              for(auto nodeB : nodeBranching[nodeBranchOrigin])
                queueBranching.push(nodeB);
          }
          
          // Update node span X
          allNodeSpanX[nodeBranchingI] = std::make_tuple(oldMin+shiftX, oldMax+shiftX);
          float oldMinImp = std::get<0>(allNodeImportantSpanX[nodeBranchingI]);
          float oldMaxImp = std::get<1>(allNodeImportantSpanX[nodeBranchingI]);
          allNodeImportantSpanX[nodeBranchingI] = std::make_tuple(oldMinImp+shiftX, oldMaxImp+shiftX);
          
          // Update x base for next iteration
          prevX = std::get<1>(allNodeSpanX[nodeBranchingI]);
          if(isImportantPair<dataType>(tree, nodeBranchingVector[i], ImportantPairs)){
            lastIndexImportant = i;
            prevX = std::get<1>(allNodeImportantSpanX[nodeBranchingI]);
            if(i < nodeBranchingVector.size()-1 and 
              not isImportantPair<dataType>(tree, nodeBranchingVector[i+1], ImportantPairs)){
              float spanMin = std::get<0>(allNodeSpanX[nodeBranchingVector[0]]);
              float spanMax = std::get<1>(allNodeImportantSpanX[nodeBranchingVector[lastIndexImportant]]);
              prevX = (spanMin + spanMax)/2;
            }
          }
        } // end for nodeBranching
        
        // Update node span X
        float spanMin = std::get<0>(allNodeSpanX[nodeBranchingVector[0]]);
        if(lastIndexImportant != -1){
          float spanMaxImp = std::get<1>(allNodeImportantSpanX[nodeBranchingVector[lastIndexImportant]]);
          allNodeImportantSpanX[nodeOrigin] = std::make_tuple(spanMin, spanMaxImp);
        }else{
          allNodeImportantSpanX[nodeOrigin] = std::make_tuple(0, 0);
          spanMin = 0;
        }
        float spanMax = std::get<1>(allNodeSpanX[nodeBranchingVector[nodeBranchingVector.size()-1]]);
        allNodeSpanX[nodeOrigin] = std::make_tuple(spanMin, spanMax);
      }else{
        allNodeSpanX[nodeOrigin] = std::make_tuple(0, 0);
        allNodeImportantSpanX[nodeOrigin] = std::make_tuple(0, 0);
      }
      
      // Positioning of this node x coordinate
      float spanMin = std::get<0>(allNodeImportantSpanX[nodeOrigin]);
      float spanMax = std::get<1>(allNodeImportantSpanX[nodeOrigin]);
      retVec[treeSimplexId[nodeOrigin]*2] = (spanMin + spanMax)/2;
    }
    
    // Copy coordinates of nodeOrigin
    for(auto node : leaves)
      queue.push(node);
    while(!queue.empty()){
      idNode node = queue.front();
      queue.pop();
      idNode nodeOrigin = tree->getNode(node)->getOrigin();
      retVec[treeSimplexId[node]*2] = retVec[treeSimplexId[nodeOrigin]*2];
      retVec[treeSimplexId[node]*2+1] = retVec[treeSimplexId[nodeOrigin]*2+1];
    }
  }

  // TODO manage multi pers pairs
  template<class dataType>
  std::vector<float> treePlanarLayoutImpl(FTMTree_MT *tree, int dataTypeInt,
                                       std::tuple<double, double, double, double, double, double> oldBounds,
                                       double refPersistence){
    int verbose = 0;
    //printTree(tree); printPairsFromTree<dataType>(tree);
    
    if(verbose > 0)
      std::cout << "================================" << std::endl << "Planar Layout" << std::endl;
    
    // ----------------------------------------------------
    // Init PlanarGraphLayout parameters
    // ----------------------------------------------------
    Timer t_init;
    if(verbose > 1) std::cout << "Init PlanarGraphLayout parameters" << std::endl;
    auto nPoints = getRealNumberOfNodes(tree);
    auto nEdges = nPoints - 1;
    auto outputArray = vtkSmartPointer<vtkFloatArray>::New();
    outputArray->SetName("Layout");
    outputArray->SetNumberOfComponents(2); // (x,y) position
    outputArray->SetNumberOfValues(nPoints * 2);
    std::vector<float> retVec(outputArray->GetNumberOfValues());
    
    //vtkSmartPointer<vtkUnstructuredGrid> connectivity = vtkSmartPointer<vtkUnstructuredGrid>::New();
    LongSimplexId connectivity[nEdges*3];
    
    dataType scalars[nPoints];
    dataType persistence[nPoints];
    int levels[nPoints];
    int branches[nPoints];
    float sizes[nPoints];
    
    // ----------------------------------------------------
    // Init internal parameters
    // ----------------------------------------------------
    if(verbose > 1) std::cout << "Init internal parameters" << std::endl;
    
    int cptNode = 0, cptEdge=0;
    std::vector<LongSimplexId> treeSimplexId(tree->getNumberOfNodes());
    std::vector<idNode> branching;
    std::vector<int> branchingID;
    std::vector<std::vector<ftm::idNode>> nodeBranching;
    getTreeBranching(tree, branching, branchingID, nodeBranching);
    
    if(verbose > 0) std::cout << "INIT            = " << t_init.getElapsedTime() << std::endl;
      
    // ----------------------------------------------------
    // Iterate through tree
    // ----------------------------------------------------
    Timer t_iterate;
    if(verbose > 1) std::cout << "Iterate through tree" << std::endl;
    
    std::queue<std::tuple<idNode, int, float>> queue;
    ftm::idNode treeRoot = getRoot(tree);
    ftm::idNode treeRootOrigin = tree->getNode(treeRoot)->getOrigin();
    queue.push(std::make_tuple(treeRoot, 0, getNumberOfLeaves(tree)));
    while(!queue.empty()){
      auto tup = queue.front();
      queue.pop();
      idNode node = std::get<0>(tup);
      int level = std::get<1>(tup);
      float pos = std::get<2>(tup);
      
      // Get and insert point
      treeSimplexId[node] = cptNode;
      scalars[cptNode] = tree->getValue<dataType>(node); //* (1/(level+1)*100);
      persistence[cptNode] = getNodePersistence<dataType>(tree, node);
      if(tree->getValue<dataType>(node) != scalars[cptNode])
        std::cout << tree->getValue<dataType>(node) << " _ " << scalars[cptNode] << std::endl;
      levels[cptNode] = level;
      branches[cptNode] = branchingID[node];
      sizes[cptNode] = pos;
      ++cptNode;
      
      // Add cell connecting parent
      if(!isRoot(tree, node)){
        vtkIdType pointIds[2];
        pointIds[0] = treeSimplexId[getParent(tree, node)];
        pointIds[1] = treeSimplexId[node];
        //connectivity->InsertNextCell(VTK_LINE, 2, pointIds);
        connectivity[cptEdge*3] = cptNode;
        connectivity[cptEdge*3+1] = treeSimplexId[getParent(tree, node)];
        connectivity[cptEdge*3+2] = treeSimplexId[node];
        ++cptEdge;
      }
      
      // Push children to the queue
      auto children = getChildren(tree, node);
      for(int i = 0; i < children.size(); ++i){
        auto child = children[i];
        float newPos = pos + i; //- children.size()/2.0; ///(16.0*(level+1)); //+1)/2 
        queue.push(std::make_tuple(child, level+1, newPos));
      }
    }
    
    if(verbose > 0) std::cout << "ITERATE TREE    = " << t_iterate.getElapsedTime() << std::endl;
    
    // ----------------------------------------------------
    // Call PlanarGraphLayout
    // ----------------------------------------------------
    /*Timer t_call;
    if(verbose > 1) std::cout << "Call PlanarGraphLayout" << std::endl;
    
    /auto idType = VTK_INT;
    PlanarGraphLayout layoutCompute;
    switch(vtkTemplate2PackMacro(idType, dataTypeInt)) {
      ttkTemplate2IdMacro(
        (//status =
          layoutCompute.execute<VTK_T1, VTK_T2>(
            // Output
            (float *)outputArray->GetVoidPointer(0),
            // Input
            //(LongSimplexId *)connectivity->GetCells()->GetPointer(), 
            (LongSimplexId *)connectivity, 
            nPoints, nEdges,
            (VTK_T2 *)scalars, nullptr, (VTK_T1 *)branches, nullptr //(float *)persistence
        )));
      // sequences (x-axis) ; sizes (y-axis) ; branches ; levels
      //(VTK_T1 *)levels (float *)sizes
    }
    auto ret = (float *)outputArray->GetVoidPointer(0);
    
    // Create output vector and swap x and y coordinates
    for(int i = 0; i < outputArray->GetNumberOfValues(); i+=2){
      //std::cout << ret[i] << std::endl;
      retVec[i] = ret[i+1];
      retVec[i+1] = ret[i];
    }
    
    if(verbose > 0) std::cout << "CALL            = " << t_call.getElapsedTime() << std::endl;*/
    
    // ----------------------------------------------------
    // Prepositioning coordinates (avoid call of PlanarGraphLayout)
    // ----------------------------------------------------
    std::queue<idNode> queue2;
    queue2.push(treeRoot);
    while(!queue2.empty()){
      idNode node= queue2.front();
      queue2.pop();
      
      retVec[treeSimplexId[node]*2] = tree->getValue<dataType>(branching[node]);
      retVec[treeSimplexId[node]*2+1] = tree->getValue<dataType>(node);
      
      // Push children to the queue
      for(auto child : getChildren(tree, node))
        queue2.push(child);
    }
    
    // ----------------------------------------------------
    // Rescale coordinates 
    // ----------------------------------------------------
    Timer t_rescale;
    if(verbose > 1) std::cout << "Rescale coordinates " << std::endl;
    
    float x_min, y_min, x_max, y_max;
    x_min = std::numeric_limits<float>::max();
    y_min = std::numeric_limits<float>::max();
    x_max = std::numeric_limits<float>::lowest();
    y_max = std::numeric_limits<float>::lowest();
    for(int i = 0; i < outputArray->GetNumberOfValues(); i+=2){
      x_min = std::min(x_min, retVec[i]);
      x_max = std::max(x_max, retVec[i]);
      y_min = std::min(y_min, retVec[i+1]);
      y_max = std::max(y_max, retVec[i+1]);
    }
    auto newBounds = std::make_tuple(x_min, x_max, y_min, y_max, 0, 0);
    
    // TODO correctly manage diff and offset if RescaleTreesIndividually
    double diff = std::max( (std::get<1>(oldBounds) - std::get<0>(oldBounds)), 
                                (std::get<3>(oldBounds) - std::get<2>(oldBounds)) );
    double offset = std::max(std::get<0>(oldBounds), std::get<2>(oldBounds));
    if(not RescaleTreesIndividually){
      //diff *= getNodePersistence<dataType>(tree, treeRoot) / refPersistence;
      diff = getNodePersistence<dataType>(tree, treeRoot);
      offset = getBirth<dataType>(tree, treeRoot);
    }
    
    for(int i = 0; i < outputArray->GetNumberOfValues(); i+=2){
      auto divisor1 = std::get<1>(newBounds) - std::get<0>(newBounds); // (x_max - x_min)
      divisor1 = (divisor1 == 0 ? 1 : divisor1);
      auto divisor2 = std::get<3>(newBounds) - std::get<2>(newBounds); // (y_max - y_min)
      divisor2 = (divisor2 == 0 ? 1 : divisor2);
      
      // x coordinate
      retVec[i] = (retVec[i] - std::get<0>(newBounds)) / divisor1;
      retVec[i] = retVec[i] * diff/2 + offset;
      
      // y coordinate
      retVec[i+1] = (retVec[i+1] - std::get<2>(newBounds)) / divisor2;
      retVec[i+1] = retVec[i+1] * diff + offset;
    }
    
    if(verbose > 0) std::cout << "RESCALE COORD.  = " << t_rescale.getElapsedTime() << std::endl;

    // ----------------------------------------------------
    // Call Branch Decomposition Planar Layout if asked
    // ----------------------------------------------------
    if(BranchDecompositionPlanarLayout){
      treePlanarLayoutBDImpl<dataType>(tree, retVec, treeSimplexId, branching, nodeBranching);
      return retVec;
    }
    
    // ----------------------------------------------------
    // Move nodes given scalars
    // ----------------------------------------------------
    Timer t_move;
    if(verbose > 1) std::cout << "Move nodes given scalars" << std::endl;
    
    float rootY = retVec[treeSimplexId[treeRoot]*2+1];
    float rootOriginY = retVec[treeSimplexId[treeRootOrigin]*2+1];
    float rootYmin = std::min(rootY, rootOriginY);
    float rootYmax = std::max(rootY, rootOriginY);
    auto rootBirthDeath = getBirthDeath<dataType>(tree, treeRoot);
    dataType rootBirth = std::get<0>(rootBirthDeath);
    dataType rootDeath = std::get<1>(rootBirthDeath);
    for(int i = 0; i < tree->getNumberOfNodes(); ++i){
      retVec[treeSimplexId[i]*2+1] = (tree->getValue<dataType>(i) - rootBirth) / (rootDeath - rootBirth);
      retVec[treeSimplexId[i]*2+1] = retVec[treeSimplexId[i]*2+1] * (rootYmax - rootYmin) + rootYmin;
    }
    
    if(verbose > 0) std::cout << "MOVE SCALAR     = " << t_move.getElapsedTime() << std::endl;
    
    // ----------------------------------------------------
    // Scale pairs given persistence
    // ----------------------------------------------------
    Timer t_scale;
    if(verbose > 1) std::cout << "Scale pairs given persistence" << std::endl;
    
    //printTree(tree);
    dataType rootPers = getNodePersistence<dataType>(tree, treeRoot);
    
    std::vector<ftm::idNode> leaves = getLeaves(tree);
    auto comp = [&](const ftm::idNode a, const ftm::idNode b) {
      return getNodePersistence<dataType>(tree, a) < getNodePersistence<dataType>(tree, b);
    };
    sort(leaves.begin(), leaves.end(), comp);
    std::stack<ftm::idNode> stack;
    for(auto node : leaves)
      stack.push(node);
    std::vector<bool> nodeDone(tree->getNumberOfNodes(), false);
    while(!stack.empty()){
      ftm::idNode node = stack.top();
      stack.pop();
      nodeDone[node] = true;    

      if(node == treeRoot or node == treeRootOrigin or isNodeAlone(tree, node))
        continue;
      
      dataType nodePers = getNodePersistence<dataType>(tree, node);
      ftm::idNode nodeOrigin = tree->getNode(node)->getOrigin();
      
      // Manage leaf
      if(isLeaf(tree, node)){
        float nodeDiff = (retVec[treeSimplexId[node]*2] - retVec[treeSimplexId[nodeOrigin]*2]);
        int sign = nodeDiff / std::abs(nodeDiff);
        auto inc = sign * nodePers / rootPers * (rootYmax - rootYmin) / 2;
        retVec[treeSimplexId[node]*2] = retVec[treeSimplexId[nodeOrigin]*2] + inc;
        
        // Push nodes in the branch to the stack
        ftm::idNode nodeParent = getParent(tree, node);
        ftm::idNode oldNodeParent = -1;
        while(nodeParent != nodeOrigin){
          if(not nodeDone[nodeParent])
            stack.push(nodeParent);
          else
            break;
          oldNodeParent = nodeParent;
          nodeParent = getParent(tree, nodeParent);
          if(oldNodeParent == nodeParent){
            std::cout << "treePlanarLayoutImpl oldNodeParent == nodeParent" << std::endl;
            break;
          }
        }
      }
      
      // Manage saddle
      if(not isLeaf(tree, node) and not isRoot(tree, node)){
        float branchY = retVec[treeSimplexId[tree->getNode(branching[node])->getOrigin()]*2];
        retVec[treeSimplexId[node]*2] = branchY;
      }
    }
    
    if(verbose > 0) std::cout << "SCALE PERS.     = " << t_scale.getElapsedTime() << std::endl;
    
    // ----------------------------------------------------
    // Branches positionning and avoid edges crossing
    // ----------------------------------------------------
    Timer t_avoid;
    if(verbose > 1) std::cout << "Avoid edges crossing" << std::endl;
    
    bool isJT = isJoinTree<dataType>(tree);
    auto compValue = [&](const ftm::idNode a, const ftm::idNode b) {
      return (isJT ? tree->getValue<dataType>(a) < tree->getValue<dataType>(b) : 
                    tree->getValue<dataType>(a) > tree->getValue<dataType>(b));
    };
    std::vector<std::tuple<float,float,float,float>> allBranchBounds(tree->getNumberOfNodes());
    std::vector<std::vector<ftm::idNode>> allBranchOrigins(tree->getNumberOfNodes());
    std::vector<int> allBranchOriginsSize(tree->getNumberOfNodes());
    std::queue<ftm::idNode> queueCrossing;
    
    // ----- Get important and non-important pairs gap and store saddle nodes of each branch
    int maxSize = std::numeric_limits<int>::lowest();
    for(auto leaf : leaves)
      queueCrossing.push(leaf);
    while(!queueCrossing.empty()){
      ftm::idNode node = queueCrossing.front();
      queueCrossing.pop();
      ftm::idNode nodeOrigin = tree->getNode(node)->getOrigin();
      
      // Get saddle nodes in the branch
      auto tupBranchOrigins = getBranchOriginsFromThisBranch(tree, nodeOrigin);
      allBranchOrigins[nodeOrigin] = std::get<0>(tupBranchOrigins);
      std::vector<ftm::idNode> nonBranchOrigins = std::get<1>(tupBranchOrigins);
      allBranchOrigins[nodeOrigin].insert(allBranchOrigins[nodeOrigin].end(), nonBranchOrigins.begin(), 
                                          nonBranchOrigins.end());
      sort(allBranchOrigins[nodeOrigin].begin(), allBranchOrigins[nodeOrigin].end(), compValue);
      allBranchOriginsSize[nodeOrigin] = allBranchOrigins[nodeOrigin].size();
      
      // Get sizes of sub-branches if they are non-important
      for(int i = 0; i < allBranchOrigins[nodeOrigin].size(); ++i){      
        ftm::idNode branchNodeOrigin = allBranchOrigins[nodeOrigin][i];
        bool isSubBranchImportant = isImportantPair<dataType>(tree, branchNodeOrigin, ImportantPairs);
        if(not isSubBranchImportant)
          allBranchOriginsSize[nodeOrigin] += allBranchOriginsSize[branchNodeOrigin];
      }

      if(isImportantPair<dataType>(tree, nodeOrigin, ImportantPairs))
        maxSize = std::max(maxSize, allBranchOriginsSize[nodeOrigin]);
    }
    double nonImportantPairsGap = (rootYmax - rootYmin) * 0.005 * NonImportantPairsSpacing;
    double importantPairsGap = (maxSize) * nonImportantPairsGap * 1.05;
    bool customImportantPairsSpacing = importantPairsGap < ImportantPairsSpacing;
    if(customImportantPairsSpacing)
      importantPairsGap = ImportantPairsSpacing;
    
    // ----- Positioning of branches and avoid conflict
    for(auto leaf : leaves)
      queueCrossing.push(leaf);
    while(!queueCrossing.empty()){
      ftm::idNode node = queueCrossing.front();
      queueCrossing.pop();
      ftm::idNode nodeOrigin = tree->getNode(node)->getOrigin();
      
      // Prepositioning of branches
      //auto restrictedBounds = getBranchBounds(retVec, treeSimplexId, branching, tree, nodeOrigin, true);
      auto restrictedBounds = std::make_tuple(retVec[treeSimplexId[node]*2], 
                                              retVec[treeSimplexId[node]*2], 0, 0);
      for(int i = 0; i < allBranchOrigins[nodeOrigin].size(); ++i){      
        ftm::idNode branchNodeOrigin = allBranchOrigins[nodeOrigin][i];
        ftm::idNode branchNode = tree->getNode(branchNodeOrigin)->getOrigin();
        
        bool isSubBranchImportant = isImportantPair<dataType>(tree, branchNodeOrigin, ImportantPairs);
        bool toLeft = not isSubBranchImportant;
        
        float branchNodeOriginXmin = std::get<0>(allBranchBounds[branchNodeOrigin]);
        float branchNodeOriginXmax = std::get<1>(allBranchBounds[branchNodeOrigin]);
        float shift = toLeft ? std::get<0>(restrictedBounds) - branchNodeOriginXmax :
                              //std::get<1>(restrictedBounds) - branchNodeOriginXmin;
                              std::get<1>(restrictedBounds) - retVec[treeSimplexId[branchNode]*2];
        shift += (toLeft ? -1 : 1) * (isSubBranchImportant ? importantPairsGap : nonImportantPairsGap);
        //shift += (toLeft ? -1 : 1) * nonImportantPairsGap;
        shiftBranchBounds(retVec, treeSimplexId, branching, allBranchBounds, 
                          allBranchOrigins[branchNodeOrigin], tree, branchNodeOrigin, shift);
      }

      // Shift a branch if conflict with another one
      for(int i = 1; i < allBranchOrigins[nodeOrigin].size(); ++i){
        ftm::idNode branchNodeOrigin = allBranchOrigins[nodeOrigin][i];
        ftm::idNode branchNode = tree->getNode(branchNodeOrigin)->getOrigin();
        for(int j = 0; j < i; ++j){
          auto first = allBranchBounds[branchNodeOrigin];
          ftm::idNode previousBranchNodeOrigin = allBranchOrigins[nodeOrigin][j];
          auto second = allBranchBounds[previousBranchNodeOrigin];
          
          bool branchConflict = isConflictingBranchAndBound(first, second, tree, previousBranchNodeOrigin, 
                                                            retVec, treeSimplexId);
          if(isConflictingBounds(first, second) or branchConflict){
            
            // Get left or right orientation given the branch
            int lastIndex = nodeBranching[branchNodeOrigin].size()-1;
            bool isLeft = (retVec[treeSimplexId[nodeBranching[branchNodeOrigin][lastIndex]]*2]
                          //< retVec[treeSimplexId[nodeOrigin]*2]);
                          < retVec[treeSimplexId[node]*2]);
            
            // Get shift
            float branchNodeOriginXmax = std::get<1>(first);
            float previousbranchNodeOriginXmin = std::get<0>(second);
            float branchNodeOriginXmin = std::get<0>(first);
            float previousbranchNodeOriginXmax = std::get<1>(second);
            float shift = isLeft ? previousbranchNodeOriginXmin - branchNodeOriginXmax: 
                                  //previousbranchNodeOriginXmax - branchNodeOriginXmin;
                                  previousbranchNodeOriginXmax - retVec[treeSimplexId[branchNode]*2];
            bool isSubBranchImportant = isImportantPair<dataType>(tree, branchNodeOrigin, 
                                                                  ImportantPairs);
            shift += (isLeft ? -1 : 1) * (isSubBranchImportant ? importantPairsGap : nonImportantPairsGap);
            //shift += (isLeft ? -1 : 1) * nonImportantPairsGap;
                                  
            // Shift bounds
            shiftBranchBounds(retVec, treeSimplexId, branching, allBranchBounds, 
                              allBranchOrigins[branchNodeOrigin], tree, branchNodeOrigin, shift);
          }
        }
      } // end for

      // TODO optimize get branch bounds by using results previously computed
      // Get branch x and y bounds
      allBranchBounds[nodeOrigin] = getBranchBounds(retVec, treeSimplexId, branching, tree, nodeOrigin);
      
      // Verify conflict (testing)
      /*std::cout << "====================================" << std::endl;
      for(int i = 0; i < allBranchOrigins[nodeOrigin].size(); ++i)
        printTuple(allBranchBounds[allBranchOrigins[nodeOrigin][i]]);*/
      /*for(int i = 1; i < allBranchOrigins[nodeOrigin].size(); ++i){
        ftm::idNode branchNodeOrigin = allBranchOrigins[nodeOrigin][i];
        for(int j = 0; j < i; ++j){
          ftm::idNode previousBranchNodeOrigin = allBranchOrigins[nodeOrigin][j];
          auto first = allBranchBounds[branchNodeOrigin];
          auto second = allBranchBounds[previousBranchNodeOrigin];
  std::cout << "-------------------------" << std::endl;
  std::cout << branchNodeOrigin << " _ " << previousBranchNodeOrigin << std::endl;
  printTuple(first); printTuple(second);
          if(isConflictingBounds(first, second))
            std::cout << i << " _ " << j << " conflict" << std::endl;
        }
      }*/
      
    } // end while  
    
    // ----- Correction of important/non-important pairs gap
    // Get real gap
    double realImportantPairsGap = std::numeric_limits<double>::lowest();
    /*if(customImportantPairsSpacing)
      realImportantPairsGap = importantPairsGap;
    else{
      for(auto leaf : leaves)
        queueCrossing.push(leaf);
      while(!queueCrossing.empty()){
        ftm::idNode node = queueCrossing.front();
        queueCrossing.pop();
        ftm::idNode nodeOrigin = tree->getNode(node)->getOrigin();
        //if(isRoot(tree, nodeOrigin)) continue; // TODO manage root gap and gap for the others
        
        bool isBranchImportant = isImportantPair<dataType>(tree, nodeOrigin, ImportantPairs);
        if(not isBranchImportant)
          continue;
        
        double gap = retVec[treeSimplexId[node]*2] - std::get<0>(allBranchBounds[nodeOrigin]);
        realImportantPairsGap = std::max(gap, realImportantPairsGap);
      }
      realImportantPairsGap *= 1.05;
    }*/
    realImportantPairsGap = importantPairsGap;
    
    // Shift given real gap
    for(auto leaf : leaves)
      queueCrossing.push(leaf);
    while(!queueCrossing.empty()){
      ftm::idNode node = queueCrossing.front();
      queueCrossing.pop();
      ftm::idNode nodeOrigin = tree->getNode(node)->getOrigin();
      
      bool isBranchImportant = isImportantPair<dataType>(tree, nodeOrigin, ImportantPairs);
      if(not isBranchImportant)
        continue;
      
      for(int i = 0; i < allBranchOrigins[nodeOrigin].size(); ++i){      
        ftm::idNode branchNodeOrigin = allBranchOrigins[nodeOrigin][i];
        bool isSubBranchImportant = isImportantPair<dataType>(tree, branchNodeOrigin, ImportantPairs);
        double shift = 0;
        if(not isSubBranchImportant){
          double gap = retVec[treeSimplexId[node]*2] - std::get<0>(allBranchBounds[nodeOrigin]);        
          double offset = (realImportantPairsGap-gap) * NonImportantPairsProximity;
          shift = -offset;
        }else{
          shift = -(importantPairsGap - realImportantPairsGap); //TODO multi shift depending on conflict
        }
        shiftBranchBounds(retVec, treeSimplexId, branching, allBranchBounds, 
                          allBranchOrigins[branchNodeOrigin], tree, branchNodeOrigin, shift);
      }
    }
    
    if(verbose > 0) std::cout << "AVOID CROSSING  = " << t_avoid.getElapsedTime() << std::endl;
    
    if(verbose > 0)
      std::cout << "================================" << std::endl;
    
    return retVec;
  }

  std::vector<float> treePlanarLayout(FTMTree_MT *tree, int dataType,
                                      std::tuple<double, double, double, double, double, double> oldBounds,
                                      double refPersistence){
    std::vector<float> res;
    switch(dataType){
      vtkTemplateMacro(res = treePlanarLayoutImpl<VTK_TT>(tree, dataType, oldBounds, refPersistence));
    }
    return res;
  }
  
  // ==========================================================================
  // Bounds Utils
  // ==========================================================================
  void printTuple(std::tuple<float,float,float,float> tup){
    std::cout << "-------------------------" << std::endl;
    std::cout << std::get<0>(tup) << " _ " << std::get<1>(tup) << " _ " << std::get<2>(tup) << " _ " << std::get<3>(tup) << " _ " << std::endl;
  }

  // Branchroot must be a non-leaf node
  std::tuple<float, float, float, float> 
  getBranchBounds(std::vector<float> &retVec, std::vector<LongSimplexId> &treeSimplexId,
                  std::vector<ftm::idNode> &branching,
                  FTMTree_MT *tree, ftm::idNode branchRoot, bool restricted=false){
    float x_min = std::numeric_limits<float>::max();
    float y_min = std::numeric_limits<float>::max();
    float x_max = std::numeric_limits<float>::lowest();
    float y_max = std::numeric_limits<float>::lowest();

    std::queue<ftm::idNode> queue;
    queue.push(branchRoot);
    while(!queue.empty()){
      ftm::idNode node = queue.front();
      queue.pop();
      
      // Skip if we go in the branch in which is branchRoot
      if(branching[node] != branchRoot and getParent(tree, node) == branchRoot and node != branchRoot)
        continue;

      // Skip if restricted
      if(restricted and !isLeaf(tree, node) and branching[node] != branchRoot and node != branchRoot)
        continue;
      
      y_min = std::min(y_min, retVec[treeSimplexId[node]*2+1]);
      y_max = std::max(y_max, retVec[treeSimplexId[node]*2+1]);
      if(node != branchRoot){
        x_min = std::min(x_min, retVec[treeSimplexId[node]*2]);
        x_max = std::max(x_max, retVec[treeSimplexId[node]*2]);
      }

      for(auto child : getChildren(tree, node))
        queue.push(child);
    }
    
    return std::make_tuple(x_min, x_max, y_min, y_max);
  }

  bool isConflictingBoundsXOneWay(std::tuple<float,float,float,float> first, 
                                  std::tuple<float,float,float,float> second){
    double eps = std::numeric_limits<float>::epsilon();
    return (std::get<0>(first) <= std::get<0>(second) and std::get<0>(second) <= std::get<1>(first)) or
          (std::get<0>(first) <= std::get<1>(second) and std::get<1>(second) <= std::get<1>(first)) or
          (isEqual<float>(std::get<0>(first), std::get<0>(second), eps) or 
            isEqual<float>(std::get<0>(second), std::get<1>(first), eps) or
            isEqual<float>(std::get<0>(first), std::get<1>(second), eps) or 
            isEqual<float>(std::get<1>(second), std::get<1>(first), eps) );
  }

  bool isConflictingBoundsX(std::tuple<float,float,float,float> first, 
                            std::tuple<float,float,float,float> second){
    return isConflictingBoundsXOneWay(first, second) or isConflictingBoundsXOneWay(second, first);
  }

  bool isConflictingBoundsYOneWay(std::tuple<float,float,float,float> first, 
                            std::tuple<float,float,float,float> second){
    double eps = std::numeric_limits<float>::epsilon();
    return (std::get<2>(first) <= std::get<2>(second) and std::get<2>(second) <= std::get<3>(first)) or
          (std::get<2>(first) <= std::get<3>(second) and std::get<3>(second) <= std::get<3>(first)) or
          (isEqual<float>(std::get<2>(first), std::get<2>(second), eps) or 
            isEqual<float>(std::get<2>(second), std::get<3>(first), eps) or
            isEqual<float>(std::get<2>(first), std::get<3>(second), eps) or 
            isEqual<float>(std::get<3>(second), std::get<3>(first), eps) );
  }

  bool isConflictingBoundsY(std::tuple<float,float,float,float> first, 
                            std::tuple<float,float,float,float> second){
    return isConflictingBoundsYOneWay(first, second) or isConflictingBoundsYOneWay(second, first);
  }

  bool isConflictingBounds(std::tuple<float,float,float,float> first, 
                          std::tuple<float,float,float,float> second){
    return isConflictingBoundsX(first, second) and isConflictingBoundsY(first, second);
  }

  bool isConflictingBranchAndBound(std::tuple<float,float,float,float> first, 
                                  std::tuple<float,float,float,float> second, 
                                  FTMTree_MT *tree, ftm::idNode branchNodeOrigin,
                                  std::vector<float> &retVec, std::vector<LongSimplexId> &treeSimplexId){
    float xBranchNodeOrigin = retVec[treeSimplexId[branchNodeOrigin]*2];
    float xBranchNode = retVec[treeSimplexId[tree->getNode(branchNodeOrigin)->getOrigin()]*2];
    float myMin = std::min(xBranchNode, xBranchNodeOrigin);
    float myMax = std::max(xBranchNode, xBranchNodeOrigin);
    auto branchBounds = std::make_tuple(myMin, myMax, 0, 0);
    return isConflictingBoundsX(first, branchBounds) and isConflictingBoundsY(first, second);
  }

  std::tuple<float,float,float,float>
  shiftBranchBoundsTuple(std::tuple<float,float,float,float> branchBound, float realShift){
    return std::make_tuple(std::get<0>(branchBound) + realShift, std::get<1>(branchBound) + realShift,
                          std::get<2>(branchBound), std::get<3>(branchBound));
  }

  void shiftBranchBounds(std::vector<float> &retVec, std::vector<LongSimplexId> &treeSimplexId,
                        std::vector<ftm::idNode> &branching,
                        std::vector<std::tuple<float,float,float,float>> &allBranchBounds, 
                        std::vector<ftm::idNode> &branchOrigins,
                        FTMTree_MT *tree, ftm::idNode branchRoot, float shift){
    std::queue<ftm::idNode> queue;
    queue.push(branchRoot);
    while(!queue.empty()){
      ftm::idNode node = queue.front();
      queue.pop();
      
      if(branching[node] != branchRoot and getParent(tree, node) == branchRoot and node != branchRoot)
        continue;
      
      if(node != branchRoot)
        retVec[treeSimplexId[node]*2] += shift;
      
      for(auto child : getChildren(tree, node))
        queue.push(child);
    }
    allBranchBounds[branchRoot] = shiftBranchBoundsTuple(allBranchBounds[branchRoot], shift);
    for(auto node : branchOrigins)
      allBranchBounds[node] = shiftBranchBoundsTuple(allBranchBounds[node], shift);
  }
};

#endif
