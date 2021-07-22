#include "ttkMergeTreeUtilsVisu.h"

// -------------------------------------------------------------------------------------
// Bounds Utils
// -------------------------------------------------------------------------------------

std::vector<double> tupleToVector(std::tuple<double,double,double,double,double,double> &tup){
  std::vector<double> vec{std::get<0>(tup), std::get<1>(tup), std::get<2>(tup), 
                          std::get<3>(tup), std::get<4>(tup), std::get<5>(tup)};
  return vec;
}

std::tuple<double,double,double,double,double,double> vectorToTuple(std::vector<double> &vec){
  return std::make_tuple(vec[0], vec[1], vec[2], vec[3], vec[4], vec[5]);
}

std::tuple<double, double, double, double, double, double> 
getMaximalBounds(std::vector<std::tuple<double, double, double, double, double, double>> &allBounds,
                 std::vector<int> &clusteringAssignment, int clusterID){
  double x_min = std::numeric_limits<double>::max();
  double y_min = std::numeric_limits<double>::max();
  double z_min = std::numeric_limits<double>::max();
  double x_max = std::numeric_limits<double>::lowest();
  double y_max = std::numeric_limits<double>::lowest();
  double z_max = std::numeric_limits<double>::lowest();
  for(int i = 0; i < allBounds.size(); ++i)
    if(clusteringAssignment[i] == clusterID){
      x_min = std::min(x_min, std::get<0>(allBounds[i]));
      x_max = std::max(x_max, std::get<1>(allBounds[i]));
      y_min = std::min(y_min, std::get<2>(allBounds[i]));
      y_max = std::max(y_max, std::get<3>(allBounds[i]));
      z_min = std::min(z_min, std::get<4>(allBounds[i]));
      z_max = std::max(z_max, std::get<5>(allBounds[i]));
    }
  return std::make_tuple(x_min, x_max, y_min, y_max, z_min, z_max);
}

std::tuple<double, double, double, double, double, double> 
getRealBounds(vtkUnstructuredGrid *treeNodes, FTMTree_MT *tree, std::vector<int> &nodeCorr){
  double x_min = std::numeric_limits<double>::max();
  double y_min = std::numeric_limits<double>::max();
  double z_min = std::numeric_limits<double>::max();
  double x_max = std::numeric_limits<double>::lowest();
  double y_max = std::numeric_limits<double>::lowest();
  double z_max = std::numeric_limits<double>::lowest();
  std::queue<idNode> queue;
  queue.push(getRoot(tree));
  while(!queue.empty()){
    idNode node = queue.front();
    queue.pop();
    double* point = treeNodes->GetPoints()->GetPoint(nodeCorr[node]);
    x_min = std::min(x_min, point[0]);
    x_max = std::max(x_max, point[0]);
    y_min = std::min(y_min, point[1]);
    y_max = std::max(y_max, point[1]);
    z_min = std::min(z_min, point[2]);
    z_max = std::max(z_max, point[2]);
    for(auto child : getChildren(tree, node))
      queue.push(child);
  }
  return std::make_tuple(x_min, x_max, y_min, y_max, z_min, z_max);
}

std::tuple<double, double, double, double, double, double> 
getRealBounds(vtkUnstructuredGrid *treeNodes, FTMTree_MT *tree){
  std::vector<int> nodeCorr(tree->getNumberOfNodes());
  for(int i = 0; i < nodeCorr.size(); ++i)
    nodeCorr[i] = i;
  
  return getRealBounds(treeNodes, tree, nodeCorr);
}

// -------------------------------------------------------------------------------------
// Temporal Subsampling MDS
// -------------------------------------------------------------------------------------
void makeTemporalSubsamplingOutput(std::vector<MergeTree*> &intermediateMTrees,
                                   std::vector<std::vector<double>> &embedding, 
                                   std::vector<MergeTree*> &allMT, std::vector<int> removed,
                                   vtkSmartPointer<vtkUnstructuredGrid> vtkOutputNode1, 
                                   vtkSmartPointer<vtkUnstructuredGrid> vtkOutputArc1,
                                   bool displayRemoved /*= false */){  
  auto fullSize = intermediateMTrees.size()+removed.size();
  std::vector<SimplexId> nodeSimplexId(fullSize);
  vtkNew<vtkIntArray> treeId{};
  treeId->SetName("TreeId");
  vtkNew<vtkIntArray> nodePathId{};
  nodePathId->SetName("NodePathId");
  vtkNew<vtkIntArray> pathId{};
  pathId->SetName("PathId");
  vtkNew<vtkIntArray> nodeRemoved{};
  nodeRemoved->SetName("isNodeRemoved");
  vtkNew<vtkIntArray> nodeId{};
  nodeId->SetName("NodeId");
  vtkNew<vtkIntArray> arcId{};
  arcId->SetName("ArcId");
  
  vtkSmartPointer<vtkUnstructuredGrid> vtkArcs = vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  for(int i = 0; i < embedding[0].size(); ++i){
    if(not displayRemoved and i >= intermediateMTrees.size())
      continue;
    
    // Next point
    double point[3] = {embedding[0][i], embedding[1][i], 0};
    nodeSimplexId[i] = points->InsertNextPoint(point);
    
    // Add TreeID field
    int index = i;
    if(i >= intermediateMTrees.size())
      index = removed[i-intermediateMTrees.size()];
    treeId->InsertNextTuple1(index);
    
    // Add NodePathId field
    int thisPathId = (i < intermediateMTrees.size() ? 0 : 1);
    nodePathId->InsertNextTuple1(thisPathId);
    
    // Add NodeId
    nodeId->InsertNextTuple1(i);
    
    // Add cell
    if(i > 0 and i < intermediateMTrees.size()){
      vtkIdType pointIds[2];
      pointIds[0] = nodeSimplexId[i-1];
      pointIds[1] = nodeSimplexId[i];
      vtkArcs->InsertNextCell(VTK_LINE, 2, pointIds);
      
      // Add PathId field
      pathId->InsertNextTuple1(0);
      
      // Add Arc Id
      arcId->InsertNextTuple1(i);
    }
  }
  
  if(displayRemoved){
    for(int i = 0; i < removed.size(); ++i){
      int index = removed[i];
      // First cell
      if(i == 0){
        vtkIdType pointIds[2];
        pointIds[0] = nodeSimplexId[index-1];
        pointIds[1] = nodeSimplexId[intermediateMTrees.size()+i];
        vtkArcs->InsertNextCell(VTK_LINE, 2, pointIds);
        pathId->InsertNextTuple1(1);
        
        if(removed[i+1] != index+1){
          vtkIdType pointIds[2];
          pointIds[0] = nodeSimplexId[intermediateMTrees.size()+i];
          pointIds[1] = nodeSimplexId[index+1];
          vtkArcs->InsertNextCell(VTK_LINE, 2, pointIds);
          pathId->InsertNextTuple1(1);
        }
      }
      // Between cell
      if(i > 0 and i < removed.size()-1){
        vtkIdType pointIds[2];
        int prevIndex = (removed[i-1] == index-1 ? intermediateMTrees.size()+i-1 : index-1);
        pointIds[0] = nodeSimplexId[prevIndex];
        pointIds[1] = nodeSimplexId[intermediateMTrees.size()+i];
        vtkArcs->InsertNextCell(VTK_LINE, 2, pointIds);
        pathId->InsertNextTuple1(1);
        
        vtkIdType pointIds2[2];
        int nextvIndex = (removed[i+1] == index+1 ? intermediateMTrees.size()+i+1 : index+1);
        pointIds2[0] = nodeSimplexId[intermediateMTrees.size()+i];
        pointIds2[1] = nodeSimplexId[nextvIndex];
        vtkArcs->InsertNextCell(VTK_LINE, 2, pointIds2);
        pathId->InsertNextTuple1(1);
      }
      // End cell
      if(i == removed.size()-1){
        vtkIdType pointIds[2];
        pointIds[0] = nodeSimplexId[intermediateMTrees.size()+i];
        pointIds[1] = nodeSimplexId[index+1];
        vtkArcs->InsertNextCell(VTK_LINE, 2, pointIds);
        pathId->InsertNextTuple1(1);
        
        if(removed[i-1] != index-1){
          vtkIdType pointIds[2];
          pointIds[0] = nodeSimplexId[index-1];
          pointIds[1] = nodeSimplexId[intermediateMTrees.size()+i];
          vtkArcs->InsertNextCell(VTK_LINE, 2, pointIds);
          pathId->InsertNextTuple1(1);
        }
      }
    }
  }else{
    for(int i = 0; i < removed.size(); ++i){
      int before = removed[i]-1;
      int index = i-1;
      if(i > 0)
        while(before == removed[index]){
          before = removed[index]-1;
          --index;
        }
      
      int after = removed[i]+1;
      index = i+1;
      if(i < removed.size()-1)
        while(after == removed[index]){
          after = removed[index]+1;
          ++index;
        }
      
      vtkIdType pointIds[2];
      pointIds[0] = nodeSimplexId[before];
      pointIds[1] = nodeSimplexId[after];
      vtkArcs->InsertNextCell(VTK_LINE, 2, pointIds);
      pathId->InsertNextTuple1(1);
      arcId->InsertNextTuple1(intermediateMTrees.size()+i);
    }
  }
  
  for(int i = 0; i < fullSize; ++i)
    nodeRemoved->InsertNextTuple1(0);
  for(int i = 0; i < removed.size(); ++i){
    nodeRemoved->SetTuple1(removed[i], 1);
    nodeRemoved->SetTuple1(intermediateMTrees.size()+i, 1);
  }
  
  vtkOutputNode1->SetPoints(points);
  vtkOutputNode1->GetPointData()->AddArray(treeId);
  vtkOutputNode1->GetPointData()->AddArray(nodePathId);
  vtkOutputNode1->GetPointData()->AddArray(nodeRemoved);
  vtkOutputNode1->GetPointData()->AddArray(nodeId);  
  vtkArcs->SetPoints(points);
  vtkArcs->GetCellData()->AddArray(pathId);
  vtkArcs->GetCellData()->AddArray(arcId);
  vtkOutputArc1->ShallowCopy(vtkArcs);
                                      
}
void makeTemporalSubsamplingMDSOutput(std::vector<MergeTree*> &intermediateMTrees, 
                                      std::vector<std::vector<double>> &distanceMatrix, 
                                      std::vector<MergeTree*> &allMT, std::vector<int> removed,
                                      vtkSmartPointer<vtkUnstructuredGrid> vtkOutputNode1, 
                                      vtkSmartPointer<vtkUnstructuredGrid> vtkOutputArc1, bool metricMDS){    
  // Call DimensionReduction filter
  std::vector<std::vector<double>> embedding{};
  std::vector<double> inputData;
  for(int i = 0; i < distanceMatrix.size(); ++i)
    for(int j = 0; j < distanceMatrix.size(); ++j)
      inputData.push_back(distanceMatrix[i][j]);

  /*DimensionReduction dimensionReduction;
  dimensionReduction.setInputModulePath("default");
  dimensionReduction.setInputModuleName("dimensionReduction");
  dimensionReduction.setInputFunctionName("doIt");
  dimensionReduction.setInputMatrixDimensions(distanceMatrix.size(), distanceMatrix.size());
  dimensionReduction.setInputMatrix(inputData.data());
  dimensionReduction.setInputMethod(2);
  dimensionReduction.setInputNumberOfComponents(2);
  dimensionReduction.setInputNumberOfNeighbors(5);
  dimensionReduction.setInputIsDeterministic(false);
  //mds_Metric, mds_Init, mds_MaxIteration, mds_Verbose, mds_Epsilon, InputIsADistanceMatrix
  dimensionReduction.setMDSParameters(metricMDS, 4, 300, 0, 0.001, true);
  dimensionReduction.setOutputComponents(&embedding);
  int errorCode = dimensionReduction.execute();
  if(errorCode != 0)
    std::cout << "errorCode = " << errorCode << std::endl;*/
  
  // Create output
  makeTemporalSubsamplingOutput(intermediateMTrees, embedding, allMT, removed, vtkOutputNode1, 
                                vtkOutputArc1);
}

void makeTemporalSubsamplingETDOutput(std::vector<MergeTree*> &intermediateMTrees, 
                                      std::vector<double> &emptyTreeDistances, 
                                      std::vector<MergeTree*> &allMT, std::vector<int> removed,
                                      vtkSmartPointer<vtkUnstructuredGrid> vtkOutputNode1, 
                                      vtkSmartPointer<vtkUnstructuredGrid> vtkOutputArc1, 
                                      double DistanceAxisStretch){
  std::vector<std::vector<double>> embedding(2);
  for(int i = 0; i < intermediateMTrees.size()+removed.size(); ++i){
    double y = emptyTreeDistances[i] * DistanceAxisStretch;
    double x = i;
    if(i >= intermediateMTrees.size())
      x = removed[i-intermediateMTrees.size()];
    
    embedding[0].push_back(x);
    embedding[1].push_back(y);
  }
  
  // Create output
  makeTemporalSubsamplingOutput(intermediateMTrees, embedding, allMT, removed, vtkOutputNode1, 
                                vtkOutputArc1);
}
