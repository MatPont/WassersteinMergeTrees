#include <ttkMergeTreeClustering.h>
#include <ttkMergeTreeUtils.h>
#include <ttkMergeTreeUtilsVisu.h>
//#include <ttkUtils.h>
#include <FTMTree.h>
#include <FTMStructures.h>
#include <MergeTreeUtils.h>
#include <ttkMergeTreeVisu.h>

#include <vtkDataObject.h> // For port information
#include <vtkObjectFactory.h> // for new macro

#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkCellData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkTable.h>
#include <vtkGeometryFilter.h>

#include <string>

using namespace ttk;
using namespace ftm;

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkMergeTreeClustering);

/**
 * Implement the filter constructor and destructor in the cpp file.
 *
 * The constructor has to specify the number of input and output ports
 * with the functions SetNumberOfInputPorts and SetNumberOfOutputPorts,
 * respectively. It should also set default values for all filter
 * parameters.
 *
 * The destructor is usually empty unless you want to manage memory
 * explicitly, by for example allocating memory on the heap that needs
 * to be freed when the filter is destroyed.
 */
ttkMergeTreeClustering::ttkMergeTreeClustering() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(8);
}

ttkMergeTreeClustering::~ttkMergeTreeClustering() {
}

/**
 * Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkMergeTreeClustering::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
  else if(port == 1){
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
  }else
    return 0;

  return 1;
}

/**
 * Specify the data object type of each output port
 *
 * This method specifies in the port information object the data type of the
 * corresponding output objects. It is possible to either explicitly
 * specify a type by adding a vtkDataObject::DATA_TYPE_NAME() key:
 *
 *      info->Set( ttkAlgorithm::DATA_TYPE_NAME(), "vtkUnstructuredGrid" );
 *
 * or to pass a type of an input port to an output port by adding the
 * ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT() key (see below).
 *
 * Note: prior to the execution of the RequestData method the pipeline will
 * initialize empty output data objects based on this information.
 */
int ttkMergeTreeClustering::FillOutputPortInformation(int port, vtkInformation *info) {
  switch(port){
    case 0:
    case 1:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      break;
    case 2:
      //info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      break;
    case 3:
    case 4:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      break;
    case 5:
      //info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      break;
    case 6: // Matching
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      break;
    case 7:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
      break;
    default:
      return 0;
  }

  return 1;
}

/**
 * Pass VTK data to the base code and convert base code output to VTK
 *
 * This method is called during the pipeline execution to update the
 * already initialized output data objects based on the given input
 * data objects and filter parameters.
 *
 * Note:
 *     1) The passed input data objects are validated based on the information
 *        provided by the FillInputPortInformation method.
 *     2) The output objects are already initialized based on the information
 *        provided by the FillOutputPortInformation method.
 */
void removePointsEpsilon(vtkSmartPointer<vtkUnstructuredGrid> vtkOutputNode, FTMTree_MT *tree){
  int noRemoved = 0;  
  vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
  for(int i = 0; i < vtkOutputNode->GetPoints()->GetNumberOfPoints(); ++i){
    if(not (tree->getNode(i)->getNumberOfUpSuperArcs() == 0 &&
        tree->getNode(i)->getNumberOfDownSuperArcs() == 0)){
      double* point = vtkOutputNode->GetPoints()->GetPoint(i);
      newPoints->InsertNextPoint(point);
    }else{
      // Remove position in arrays
      for(int j = 0; j < vtkOutputNode->GetPointData()->GetNumberOfArrays(); ++j){
        vtkOutputNode->GetPointData()->GetArray(j)->RemoveTuple(i-noRemoved);
      } 
      noRemoved++;      
    }
  }
  vtkOutputNode->SetPoints(newPoints);
}

void loadBlocks(std::vector<vtkMultiBlockDataSet *> &inputTrees, vtkMultiBlockDataSet* blocks){
  if(blocks != nullptr) {
    inputTrees.resize(blocks->GetNumberOfBlocks());
    for(size_t i = 0; i < inputTrees.size(); ++i) {
      inputTrees[i] = vtkMultiBlockDataSet::SafeDownCast(blocks->GetBlock(i));
    }
  }
}

void constructTrees(std::vector<vtkMultiBlockDataSet *> &inputTrees, 
                std::vector<MergeTree *> &intermediateTrees, std::vector<vtkUnstructuredGrid *> &treesNodes, 
                std::vector<vtkUnstructuredGrid *> &treesArcs, std::vector<vtkDataSet *> &treesSegmentation){
  const int numInputs = inputTrees.size();
  for(int i = 0; i < numInputs; i++) {
    //std::cout << "construct " << i << std::endl;
    treesNodes[i] = vtkUnstructuredGrid::SafeDownCast(inputTrees[i]->GetBlock(0));
    treesArcs[i] = vtkUnstructuredGrid::SafeDownCast(inputTrees[i]->GetBlock(1));
    treesSegmentation[i] = vtkDataSet::SafeDownCast(inputTrees[i]->GetBlock(2));
    intermediateTrees[i] = makeTree(treesNodes[i], treesArcs[i]);
  }
}

int ttkMergeTreeClustering::RequestData(vtkInformation *request,
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) { 
  int verbose = 0;
  // ------------------------------------------------------------------------------------
  // --- Get input object from input vector
  // ------------------------------------------------------------------------------------
  vtkUnstructuredGrid *tree1Nodes; //= vtkUnstructuredGrid::GetData(inputVector[0]);
  vtkUnstructuredGrid *tree1Arcs; //= vtkUnstructuredGrid::GetData(inputVector[1]);
  vtkDataSet *tree1Segmentation; //= vtkDataSet::GetData(inputVector[2]);
  
  vtkUnstructuredGrid *tree2Nodes; //= vtkUnstructuredGrid::GetData(inputVector[3]);
  vtkUnstructuredGrid *tree2Arcs; //= vtkUnstructuredGrid::GetData(inputVector[4]);
  vtkDataSet *tree2Segmentation; //= vtkDataSet::GetData(inputVector[5]);
  
  // ------------------------------------------------------------------------------------
  // --- Construct trees
  // ------------------------------------------------------------------------------------
  FTMTree_MT *tree1; //= makeTree(tree1Nodes, tree1Arcs);  
  //tree1->printTree2();
  
  FTMTree_MT *tree2; //= makeTree(tree2Nodes, tree2Arcs);
  //tree2->printTree2();

  // ------------------------------------------------------------------------------------
  // --- Get input object from input vector
  // ------------------------------------------------------------------------------------
  auto blocks = vtkMultiBlockDataSet::GetData(inputVector[0], 0);
  auto blocks2 = vtkMultiBlockDataSet::GetData(inputVector[1], 0);

  // ------------------------------------------------------------------------------------
  // --- Construct trees
  // ------------------------------------------------------------------------------------
  std::vector<vtkMultiBlockDataSet *> inputTrees, inputTrees2;
  loadBlocks(inputTrees, blocks);
  loadBlocks(inputTrees2, blocks2);
  const int numInputs = inputTrees.size();
  const int numInputs2 = inputTrees2.size();
  
  std::vector<MergeTree *> intermediateMTrees(numInputs), intermediateMTrees2(numInputs2);
  std::vector<FTMTree_MT *> intermediateTrees(numInputs), intermediateTrees2(numInputs2);
  std::vector<vtkUnstructuredGrid *> treesNodes(numInputs), treesNodes2(numInputs2);
  std::vector<vtkUnstructuredGrid *> treesArcs(numInputs), treesArcs2(numInputs2);
  std::vector<vtkDataSet *> treesSegmentation(numInputs), treesSegmentation2(numInputs2);
  constructTrees(inputTrees, intermediateMTrees, treesNodes, treesArcs, treesSegmentation);
  constructTrees(inputTrees2, intermediateMTrees2, treesNodes2, treesArcs2, treesSegmentation2);
  for(int i = 0; i < intermediateMTrees.size(); ++i)
    intermediateTrees[i] = intermediateMTrees[i]->tree;
  for(int i = 0; i < intermediateMTrees2.size(); ++i)
    intermediateTrees2[i] = intermediateMTrees2[i]->tree;
  
  int dataTypeInt = treesNodes[0]->GetPointData()->GetArray("Scalar")->GetDataType();
  
  tree1 = intermediateTrees[0];
  tree1Nodes = treesNodes[0];
  tree1Arcs = treesArcs[0];
  tree1Segmentation = treesSegmentation[0];
  tree2 = intermediateTrees[1];
  tree2Nodes = treesNodes[1];
  tree2Arcs = treesArcs[1];
  tree2Segmentation = treesSegmentation[1];
  
  // ------------------------------------------------------------------------------------
  // --- Call base
  // ------------------------------------------------------------------------------------
  std::vector<std::tuple<idNode, idNode, double>> outputMatching;
  std::vector<std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>>> 
         outputMatchingBarycenter(NumberOfBarycenters, 
                                  std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>>(numInputs)),
         outputMatchingBarycenter2(NumberOfBarycenters, 
                                  std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>>(numInputs));
  std::vector<MergeTree*> barycenters(NumberOfBarycenters);
  std::vector<int> clusteringAssignment(numInputs, 0);
  std::vector<std::vector<int>> trees1NodeCorrMesh, trees2NodeCorrMesh;
  bool AddNodes=true;
  // Classical distance
  double distance = 0;
  // Temporal subsampling
  std::vector<std::vector<double>> tsDistanceMatrix;
  std::vector<double> tsEmptyTreeDistances;
  std::vector<MergeTree*> tsAllMT;
  std::vector<int> tsRemoved;
  
  /*ComputeBarycenter=true;
  BranchDecomposition=true;
  NormalizedWasserstein=true;
  RescaledWasserstein=true;
  KeepSubtree=false;
  NumberOfBarycenters=2;*/
  //AddNodes=false;
  //RescaledWasserstein=false;
  //Deterministic = false;
  
  if(TemporalSubsampling){
    MergeTreeTemporalSubsampling mergeTreeTS;
    mergeTreeTS.setAssignmentSolver(AssignmentSolver);
    mergeTreeTS.setEpsilonTree1(EpsilonTree1);
    mergeTreeTS.setEpsilonTree2(EpsilonTree2);
    mergeTreeTS.setEpsilon2Tree1(Epsilon2Tree1);
    mergeTreeTS.setEpsilon2Tree2(Epsilon2Tree2);
    mergeTreeTS.setEpsilon3Tree1(Epsilon3Tree1);
    mergeTreeTS.setEpsilon3Tree2(Epsilon3Tree2);
    mergeTreeTS.setProgressiveComputation(ProgressiveComputation);
    mergeTreeTS.setBranchDecomposition(BranchDecomposition);
    mergeTreeTS.setParallelize(Parallelize);
    mergeTreeTS.setPersistenceThreshold(PersistenceThreshold);
    mergeTreeTS.setNormalizedWasserstein(NormalizedWasserstein);
    mergeTreeTS.setNormalizedWassersteinReg(NormalizedWassersteinReg);
    mergeTreeTS.setRescaledWasserstein(RescaledWasserstein);
    mergeTreeTS.setKeepSubtree(KeepSubtree);
    mergeTreeTS.setUseMinMaxPair(UseMinMaxPair);
    mergeTreeTS.setNumberOfThreads(NumberOfThreads);
    mergeTreeTS.setRemovalPercentage(RemovalPercentage);
    mergeTreeTS.setUseL2Distance(UseL2Distance);
    mergeTreeTS.setTreesSegmentation(treesSegmentation);
    mergeTreeTS.setDeleteMultiPersPairs(DeleteMultiPersPairs);
    mergeTreeTS.setEpsilon1UseFarthestSaddle(Epsilon1UseFarthestSaddle);
    switch(dataTypeInt){
      vtkTemplateMacro(
        tsRemoved = mergeTreeTS.execute<VTK_TT>(intermediateMTrees, tsDistanceMatrix, 
                                                          tsEmptyTreeDistances, tsAllMT);
      );
    }
    trees1NodeCorrMesh = mergeTreeTS.getTreesNodeCorr();
    for(int i = 0; i < intermediateMTrees.size(); ++i)
      intermediateTrees[i] = intermediateMTrees[i]->tree;
  }else if(not ComputeBarycenter){    
    MergeTreeDistance mergeTreeDistance;
    mergeTreeDistance.setVerbose(Verbose);
    mergeTreeDistance.setAssignmentSolver(AssignmentSolver);
    mergeTreeDistance.setEpsilonTree1(EpsilonTree1);
    mergeTreeDistance.setEpsilonTree2(EpsilonTree2);
    mergeTreeDistance.setEpsilon2Tree1(Epsilon2Tree1);
    mergeTreeDistance.setEpsilon2Tree2(Epsilon2Tree2);
    mergeTreeDistance.setEpsilon3Tree1(Epsilon3Tree1);
    mergeTreeDistance.setEpsilon3Tree2(Epsilon3Tree2);
    mergeTreeDistance.setProgressiveComputation(ProgressiveComputation);
    mergeTreeDistance.setBranchDecomposition(BranchDecomposition);
    mergeTreeDistance.setParallelize(Parallelize);
    mergeTreeDistance.setPersistenceThreshold(PersistenceThreshold);
    mergeTreeDistance.setNormalizedWasserstein(NormalizedWasserstein);
    mergeTreeDistance.setNormalizedWassersteinReg(NormalizedWassersteinReg);
    mergeTreeDistance.setRescaledWasserstein(RescaledWasserstein);
    mergeTreeDistance.setKeepSubtree(KeepSubtree);
    mergeTreeDistance.setUseMinMaxPair(UseMinMaxPair);
    mergeTreeDistance.setNumberOfThreads(NumberOfThreads);
    mergeTreeDistance.setCleanTree(true);
    mergeTreeDistance.setPostprocess(OutputTrees);
    mergeTreeDistance.setDeleteMultiPersPairs(DeleteMultiPersPairs);
    mergeTreeDistance.setEpsilon1UseFarthestSaddle(Epsilon1UseFarthestSaddle);
    switch(dataTypeInt){
      vtkTemplateMacro(
        distance = mergeTreeDistance.execute<VTK_TT>(tree1, tree2, outputMatching)
      );
    }
    trees1NodeCorrMesh = mergeTreeDistance.getTreesNodeCorr();
    intermediateTrees[0] = tree1;
    intermediateTrees[1] = tree2;
  }else{
    if(NumberOfBarycenters == 1){ //and numInputs2==0){
      MergeTreeBarycenter mergeTreeBarycenter;
      mergeTreeBarycenter.setAssignmentSolver(AssignmentSolver);
      mergeTreeBarycenter.setEpsilonTree1(EpsilonTree1);
      mergeTreeBarycenter.setEpsilonTree2(EpsilonTree2);
      mergeTreeBarycenter.setEpsilon2Tree1(Epsilon2Tree1);
      mergeTreeBarycenter.setEpsilon2Tree2(Epsilon2Tree2);
      mergeTreeBarycenter.setEpsilon3Tree1(Epsilon3Tree1);
      mergeTreeBarycenter.setEpsilon3Tree2(Epsilon3Tree2);
      mergeTreeBarycenter.setProgressiveComputation(ProgressiveComputation);
      mergeTreeBarycenter.setBranchDecomposition(BranchDecomposition);
      mergeTreeBarycenter.setParallelize(Parallelize);
      mergeTreeBarycenter.setPersistenceThreshold(PersistenceThreshold);
      mergeTreeBarycenter.setNormalizedWasserstein(NormalizedWasserstein);
      mergeTreeBarycenter.setNormalizedWassersteinReg(NormalizedWassersteinReg);
      mergeTreeBarycenter.setRescaledWasserstein(RescaledWasserstein);
      mergeTreeBarycenter.setKeepSubtree(KeepSubtree);
      mergeTreeBarycenter.setUseMinMaxPair(UseMinMaxPair);
      mergeTreeBarycenter.setTol(Tol);
      mergeTreeBarycenter.setAddNodes(AddNodes);
      mergeTreeBarycenter.setDeterministic(Deterministic);
      mergeTreeBarycenter.setNumberOfThreads(NumberOfThreads);
      mergeTreeBarycenter.setProgressiveBarycenter(ProgressiveBarycenter);
      mergeTreeBarycenter.setProgressiveSpeedDivisor(ProgressiveSpeedDivisor);
      mergeTreeBarycenter.setAlpha(Alpha);
      mergeTreeBarycenter.setPostprocess(OutputTrees);
      mergeTreeBarycenter.setDeleteMultiPersPairs(DeleteMultiPersPairs);
      mergeTreeBarycenter.setEpsilon1UseFarthestSaddle(Epsilon1UseFarthestSaddle);
      switch(dataTypeInt){
        vtkTemplateMacro(
          barycenters[0] = mergeTreeBarycenter.execute<VTK_TT>(intermediateTrees, 
                                                                   outputMatchingBarycenter[0])
        );
      }
      trees1NodeCorrMesh = mergeTreeBarycenter.getTreesNodeCorr();
    }else{
      MergeTreeClustering mergeTreeClustering;
      mergeTreeClustering.setAssignmentSolver(AssignmentSolver);
      mergeTreeClustering.setEpsilonTree1(EpsilonTree1);
      mergeTreeClustering.setEpsilonTree2(EpsilonTree2);
      mergeTreeClustering.setEpsilon2Tree1(Epsilon2Tree1);
      mergeTreeClustering.setEpsilon2Tree2(Epsilon2Tree2);
      mergeTreeClustering.setEpsilon3Tree1(Epsilon3Tree1);
      mergeTreeClustering.setEpsilon3Tree2(Epsilon3Tree2);
      mergeTreeClustering.setProgressiveComputation(ProgressiveComputation);
      mergeTreeClustering.setBranchDecomposition(BranchDecomposition);
      mergeTreeClustering.setParallelize(Parallelize);
      mergeTreeClustering.setPersistenceThreshold(PersistenceThreshold);
      mergeTreeClustering.setNormalizedWasserstein(NormalizedWasserstein);
      mergeTreeClustering.setNormalizedWassersteinReg(NormalizedWassersteinReg);
      mergeTreeClustering.setRescaledWasserstein(RescaledWasserstein);
      mergeTreeClustering.setKeepSubtree(KeepSubtree);
      mergeTreeClustering.setUseMinMaxPair(UseMinMaxPair);
      mergeTreeClustering.setTol(Tol);
      mergeTreeClustering.setAddNodes(AddNodes);
      mergeTreeClustering.setDeterministic(Deterministic);
      mergeTreeClustering.setNumberOfThreads(NumberOfThreads);
      mergeTreeClustering.setNoCentroids(NumberOfBarycenters);
      mergeTreeClustering.setProgressiveBarycenter(ProgressiveBarycenter);
      mergeTreeClustering.setProgressiveSpeedDivisor(ProgressiveSpeedDivisor);
      mergeTreeClustering.setPostprocess(OutputTrees);
      mergeTreeClustering.setDeleteMultiPersPairs(DeleteMultiPersPairs);
      mergeTreeClustering.setEpsilon1UseFarthestSaddle(Epsilon1UseFarthestSaddle);
      switch(dataTypeInt){
        vtkTemplateMacro(
          barycenters = mergeTreeClustering.execute<VTK_TT>(intermediateTrees, 
                                                              outputMatchingBarycenter, clusteringAssignment,
                                                              intermediateTrees2, outputMatchingBarycenter2)
        );
      }
      trees1NodeCorrMesh = mergeTreeClustering.getTreesNodeCorr();
      trees2NodeCorrMesh = mergeTreeClustering.getTrees2NodeCorr();
    }
  }
  
  /*for(int i = 0; i < intermediateMTrees.size(); ++i)
    intermediateTrees[i] = intermediateMTrees[i]->tree;
  for(int i = 0; i < intermediateMTrees2.size(); ++i)
    intermediateTrees2[i] = intermediateMTrees2[i]->tree;*/
  
  // Node correspondence if not got yet
  if(trees1NodeCorrMesh.size() == 0 or trees1NodeCorrMesh[0].size() == 0){ 
    trees1NodeCorrMesh = std::vector<std::vector<int>>(numInputs);
    trees2NodeCorrMesh = std::vector<std::vector<int>>(numInputs2);
    for(int i = 0; i < intermediateTrees.size(); ++i){
      trees1NodeCorrMesh[i] = std::vector<int>(intermediateTrees[i]->getNumberOfNodes());
      for(int j = 0; j < intermediateTrees[i]->getNumberOfNodes(); ++j)
        trees1NodeCorrMesh[i][j] = j;
    }
    for(int i = 0; i < intermediateTrees2.size(); ++i){
      trees2NodeCorrMesh[i] = std::vector<int>(intermediateTrees2[i]->getNumberOfNodes());
      for(int j = 0; j < intermediateTrees2[i]->getNumberOfNodes(); ++j)
        trees2NodeCorrMesh[i][j] = j;
    }
  }
  
  // ------------------------------------------------------------------------------------
  // --- Create output
  // ------------------------------------------------------------------------------------
  Timer t_makeTreesOutput;
  
  auto vtkOutputNode1 = vtkUnstructuredGrid::SafeDownCast(
    outputVector->GetInformationObject(0)->Get(vtkDataObject::DATA_OBJECT()));
  
  auto vtkOutputArc1 = vtkUnstructuredGrid::SafeDownCast(
    outputVector->GetInformationObject(1)->Get(vtkDataObject::DATA_OBJECT()));
  
  auto vtkOutputSegmentation1 = vtkUnstructuredGrid::SafeDownCast( //vtkDataSet::SafeDownCast(
    outputVector->GetInformationObject(2)->Get(vtkDataObject::DATA_OBJECT()));
  
  auto vtkOutputNode2 = vtkUnstructuredGrid::SafeDownCast(
    outputVector->GetInformationObject(3)->Get(vtkDataObject::DATA_OBJECT()));
  
  auto vtkOutputArc2 = vtkUnstructuredGrid::SafeDownCast(
    outputVector->GetInformationObject(4)->Get(vtkDataObject::DATA_OBJECT()));
  
  auto vtkOutputSegmentation2 = vtkUnstructuredGrid::SafeDownCast( //vtkDataSet::SafeDownCast(
    outputVector->GetInformationObject(5)->Get(vtkDataObject::DATA_OBJECT()));
  
  auto vtkOutputMatching = vtkUnstructuredGrid::SafeDownCast(
    outputVector->GetInformationObject(6)->Get(vtkDataObject::DATA_OBJECT()));
  
  auto vtkAssignment = vtkTable::SafeDownCast(
    outputVector->GetInformationObject(7)->Get(vtkDataObject::DATA_OBJECT()));
  
  // Get bounds
  std::vector<std::vector<SimplexId>> nodeCorr(numInputs);
  std::vector<std::tuple<double,double,double,double,double,double>> allBounds(numInputs);
  for(int i = 0; i < numInputs; ++i){
    if(OutputSegmentation){
      double* tBounds = treesSegmentation[i]->GetBounds();
      allBounds[i] = std::make_tuple(tBounds[0], tBounds[1], tBounds[2], tBounds[3], tBounds[4], tBounds[5]);
    }else
      allBounds[i] = getRealBounds(treesNodes[i], intermediateTrees[i], trees1NodeCorrMesh[i]);
    /*if(PlanarLayout){
      // TODO correctly manage bounds when planar layout
      std::vector<double> tBound = tupleToVector(allBounds[i]);
      int level = getTreeDepth(intermediateTrees[i]);
      for(int j = 0; j < tBound.size(); ++j)
        tBound[j] = tBound[j]*level/5;
      allBounds[i] = vectorToTuple(tBound);
    }*/
  }
  
  if(TemporalSubsampling){
    // ---------------------------------------------------------------
    // - Temporal Subsampling Output
    // ---------------------------------------------------------------
    if(TemporalSubsamplingMDS){
      makeTemporalSubsamplingMDSOutput(intermediateMTrees, tsDistanceMatrix, tsAllMT, tsRemoved, 
                                       vtkOutputNode1, vtkOutputArc1, MetricMDS);
    }else{
      makeTemporalSubsamplingETDOutput(intermediateMTrees, tsEmptyTreeDistances, tsAllMT, tsRemoved, 
                                       vtkOutputNode1, vtkOutputArc1, DistanceAxisStretch);
    }
    
    //
    ttkMergeTreeVisu visuMaker;
    visuMaker.setPlanarLayout(PlanarLayout);
    visuMaker.setBranchDecompositionPlanarLayout(BranchDecompositionPlanarLayout);
    visuMaker.setRescaleTreesIndividually(RescaleTreesIndividually);
    visuMaker.setOutputSegmentation(OutputSegmentation);
    visuMaker.setDimensionSpacing(DimensionSpacing);
    visuMaker.setDimensionToShift(DimensionToShift);
    visuMaker.setImportantPairs(ImportantPairs);
    visuMaker.setImportantPairsSpacing(ImportantPairsSpacing);
    visuMaker.setNonImportantPairsSpacing(NonImportantPairsSpacing);
    visuMaker.setNonImportantPairsProximity(NonImportantPairsProximity);
    visuMaker.setShiftMode(3);
    auto tsAllTree(intermediateTrees);
    for(int i = 0; i < intermediateTrees.size(); ++i){
      tsAllTree.push_back(intermediateTrees[i]);
      treesNodes.push_back(treesNodes[i]);
      trees1NodeCorrMesh.push_back(trees1NodeCorrMesh[i]);
      treesSegmentation.push_back(treesSegmentation[i]);
      allBounds.push_back(allBounds[i]);
    }
    for(int i = 0; i < tsRemoved.size(); ++i){
      int index = intermediateTrees.size()+tsRemoved[i];
      tsAllTree[index] = tsAllMT[intermediateTrees.size()+i]->tree;
      treesNodes[index] = nullptr;
      trees1NodeCorrMesh[index].clear();
      treesSegmentation[index] = nullptr;
    }
    visuMaker.makeTreesOutput(tsAllTree, treesNodes, trees1NodeCorrMesh, treesSegmentation, 
                              allBounds, vtkOutputNode2, vtkOutputArc2, vtkOutputSegmentation2, 
                              dataTypeInt, verbose);
    
    // - Add Distance Matrix Table
    // zero-padd column name to keep Row Data columns ordered
    const auto zeroPad
      = [](std::string &colName, const size_t numberCols, const size_t colIdx) {
        std::string max{std::to_string(numberCols - 1)};
        std::string cur{std::to_string(colIdx)};
        std::string zer(max.size() - cur.size(), '0');
        colName.append(zer).append(cur);
      };
    for(size_t i = 0; i < tsDistanceMatrix.size(); ++i) {
      std::string name{};
      zeroPad(name, tsDistanceMatrix.size(), i);

      vtkNew<vtkDoubleArray> col{};
      col->SetNumberOfTuples(tsDistanceMatrix.size());
      col->SetName(name.c_str());
      for(size_t j = 0; j < tsDistanceMatrix[i].size(); ++j) {
        col->SetTuple1(j, tsDistanceMatrix[i][j]);
      }
      vtkAssignment->AddColumn(col);
    }
    auto fullSize = intermediateMTrees.size()+tsRemoved.size();
    vtkNew<vtkDoubleArray> col{};
    col->SetNumberOfTuples(fullSize);
    col->SetName("FieldId");
    for(int i = 0; i < fullSize; ++i) {
      auto newTuple = (i < intermediateMTrees.size() ? i : tsRemoved[i-intermediateMTrees.size()]);
      col->SetTuple1(i, newTuple);
    }
    vtkAssignment->AddColumn(col);
  }else if(not ComputeBarycenter){
    if(OutputTrees){
      // ---------------------------------------------------------------
      // - Classical output
      // ---------------------------------------------------------------
      ttkMergeTreeVisu visuMaker;
      visuMaker.setPlanarLayout(PlanarLayout);
      visuMaker.setBranchDecompositionPlanarLayout(BranchDecompositionPlanarLayout);
      visuMaker.setBranchSpacing(BranchSpacing);
      visuMaker.setRescaleTreesIndividually(RescaleTreesIndividually);
      visuMaker.setOutputSegmentation(OutputSegmentation);
      visuMaker.setDimensionSpacing(DimensionSpacing);
      visuMaker.setDimensionToShift(DimensionToShift);
      visuMaker.setImportantPairs(ImportantPairs);
      visuMaker.setImportantPairsSpacing(ImportantPairsSpacing);
      visuMaker.setNonImportantPairsSpacing(NonImportantPairsSpacing);
      visuMaker.setNonImportantPairsProximity(NonImportantPairsProximity);
      if(not SepareOutputTrees){
        visuMaker.makeTreesOutput(tree1, tree2, treesNodes, trees1NodeCorrMesh, treesSegmentation, 
                                  allBounds, vtkOutputNode1, vtkOutputArc1, vtkOutputSegmentation1, 
                                  dataTypeInt, verbose);
        nodeCorr = visuMaker.getNodeCorr();
      }else{
        nodeCorr.clear();
        visuMaker.setNoSampleOffset(1);
        visuMaker.makeTreesOutput(tree1, treesNodes[0], trees1NodeCorrMesh[0], treesSegmentation[0], 
                                allBounds[0], vtkOutputNode1, vtkOutputArc1, vtkOutputSegmentation1, 
                                dataTypeInt, verbose);
        nodeCorr.push_back(visuMaker.getNodeCorr()[0]);
        visuMaker.setISampleOffset(1);
        visuMaker.makeTreesOutput(tree2, treesNodes[1], trees1NodeCorrMesh[1], treesSegmentation[1], 
                                allBounds[1], vtkOutputNode2, vtkOutputArc2, vtkOutputSegmentation2, 
                                dataTypeInt, verbose);
        nodeCorr.push_back(visuMaker.getNodeCorr()[0]);
      }
    
      // TODO merge with barycenter matching in ttkMergeTreeVisu class
      // Create matching
      vtkNew<vtkIntArray> matchingID{};
      matchingID->SetName("MatchingID");
      vtkNew<vtkIntArray> matchingType{};
      matchingType->SetName("MatchingType");
      vtkNew<vtkDoubleArray> matchPers{};
      matchPers->SetName("MeanMatchedPersistence");
      vtkNew<vtkDoubleArray> costArray{};
      costArray->SetName("Cost");
      vtkNew<vtkIntArray> tree1NodeIdField{};
      tree1NodeIdField->SetName("tree1NodeId");
      vtkNew<vtkIntArray> tree2NodeIdField{};
      tree2NodeIdField->SetName("tree2NodeId");
    
      vtkSmartPointer<vtkUnstructuredGrid> vtkMatching = vtkSmartPointer<vtkUnstructuredGrid>::New();
      vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
      int count = 0;
      for(auto match : outputMatching){
        vtkIdType pointIds[2];
        ftm::idNode tree1NodeId = std::get<0>(match);
        ftm::idNode tree2NodeId = std::get<1>(match);
        double cost = std::get<2>(match);
        
        // Get first point
        double* point1 = vtkOutputNode1->GetPoints()->GetPoint(nodeCorr[0][tree1NodeId]);
        const SimplexId nextPointId1 = points->InsertNextPoint(point1);
        pointIds[0] = nextPointId1;
        
        // Get second point
        double* point2;
        if(not SepareOutputTrees)
          point2 = vtkOutputNode1->GetPoints()->GetPoint(nodeCorr[1][tree2NodeId]);
        else
          point2 = vtkOutputNode2->GetPoints()->GetPoint(nodeCorr[1][tree2NodeId]);
        const SimplexId nextPointId2 = points->InsertNextPoint(point2);
        pointIds[1] = nextPointId2;
  
        // Add cell
        vtkMatching->InsertNextCell(VTK_LINE, 2, pointIds);
        
        // Add tree1 and tree2 node ids
        tree1NodeIdField->InsertNextTuple1(nodeCorr[0][tree1NodeId]);
        tree2NodeIdField->InsertNextTuple1(nodeCorr[1][tree2NodeId]);
        
        // Add matching ID
        matchingID->InsertNextTuple1(count);
        
        // Add matching type
        int thisType = 0;
        int tree1NodeDown = tree1->getNode(tree1NodeId)->getNumberOfDownSuperArcs();
        int tree1NodeUp = tree1->getNode(tree1NodeId)->getNumberOfUpSuperArcs();
        int tree2NodeDown = tree2->getNode(tree2NodeId)->getNumberOfDownSuperArcs();
        int tree2NodeUp = tree2->getNode(tree2NodeId)->getNumberOfUpSuperArcs();
        if(tree1NodeDown != 0 and tree1NodeUp != 0 and tree2NodeDown != 0 and tree2NodeUp != 0)
          thisType = 1; // Saddle to Saddle
        if(tree1NodeDown == 0 and tree1NodeUp != 0 and tree2NodeDown == 0 and tree2NodeUp != 0)
          thisType = 2; // Leaf to leaf
        if(tree1NodeDown != 0 and tree1NodeUp == 0 and tree2NodeDown != 0 and tree2NodeUp == 0)
          thisType = 3; // Root to root
        matchingType->InsertTuple1(count, thisType);
        
        // Add mean matched persistence
        double tree1Pers = getPersistence(tree1, tree1NodeId, dataTypeInt);
        double tree2Pers = getPersistence(tree2, tree2NodeId, dataTypeInt);
        double meanPersistence = (tree1Pers+tree2Pers)/2;
        matchPers->InsertTuple1(count, meanPersistence);
        
        // Add cost
        costArray->InsertNextTuple1(cost);
        
        count++;
      }
      vtkMatching->SetPoints(points);
      vtkMatching->GetCellData()->AddArray(matchingType);
      vtkMatching->GetCellData()->AddArray(matchPers);
      vtkMatching->GetCellData()->AddArray(matchingID);
      vtkMatching->GetCellData()->AddArray(costArray);
      vtkMatching->GetCellData()->AddArray(tree1NodeIdField);
      vtkMatching->GetCellData()->AddArray(tree2NodeIdField);
      vtkOutputMatching->ShallowCopy(vtkMatching);
      
      // Remove points of epsilon merging
      /*removePointsEpsilon(vtkOutputNode1, tree1);  
      removePointsEpsilon(vtkOutputNode2, tree2);*/
      
      // Add distance value to the table
      vtkNew<vtkDoubleArray> col{};
      col->SetName("Distance");
      col->SetNumberOfTuples(1);
      col->SetTuple1(0, distance);
      vtkAssignment->AddColumn(col);
    }
  }else {
    // ---------------------------------------------------------------
    // - Barycenter output
    // ---------------------------------------------------------------
    // - Assignment Table
    vtkNew<vtkDoubleArray> col{};
    col->SetNumberOfTuples(numInputs);
    col->SetName("ClusterAssignment");
    for(int i = 0; i < numInputs; ++i)
      col->SetTuple1(i, clusteringAssignment[i]);
    vtkAssignment->AddColumn(col);

    if(OutputTrees){
      // --- Declare internal arrays 
      std::vector<double> clusterShift(NumberOfBarycenters, 0);
      std::vector<std::tuple<double,double,double,double,double,double>> allBaryBounds(NumberOfBarycenters);
      std::vector<std::vector<ftm::idNode>> allBaryBranching(NumberOfBarycenters);
      std::vector<std::vector<int>> allBaryBranchingID(NumberOfBarycenters);
      for(int c = 0; c < NumberOfBarycenters; ++c){
        allBaryBounds[c] = getMaximalBounds(allBounds, clusteringAssignment, c);
        getTreeBranching(barycenters[c]->tree, allBaryBranching[c], allBaryBranchingID[c]);
      } 
     
      // ------------------------------------------
      // --- Input trees
      // ------------------------------------------
      ttkMergeTreeVisu visuMaker;
      visuMaker.setPlanarLayout(PlanarLayout);
      visuMaker.setBranchDecompositionPlanarLayout(BranchDecompositionPlanarLayout);
      visuMaker.setBranchSpacing(BranchSpacing);
      visuMaker.setRescaleTreesIndividually(RescaleTreesIndividually);
      visuMaker.setOutputSegmentation(OutputSegmentation);
      visuMaker.setDimensionSpacing(DimensionSpacing);
      visuMaker.setDimensionToShift(DimensionToShift);
      visuMaker.setImportantPairs(ImportantPairs);
      visuMaker.setImportantPairsSpacing(ImportantPairsSpacing);
      visuMaker.setNonImportantPairsSpacing(NonImportantPairsSpacing);
      visuMaker.setNonImportantPairsProximity(NonImportantPairsProximity);
      visuMaker.makeTreesOutput(intermediateTrees, barycenters, clusteringAssignment, 
                                outputMatchingBarycenter, treesNodes, trees1NodeCorrMesh, treesSegmentation, 
                                allBounds, allBaryBounds, allBaryBranchingID, 
                                vtkOutputNode1, vtkOutputArc1, vtkOutputSegmentation1, dataTypeInt, verbose);
      nodeCorr = visuMaker.getNodeCorr();
      clusterShift = visuMaker.getClusterShift();
    
      // ------------------------------------------
      // --- Barycenter
      // ------------------------------------------
      if(verbose > 0) std::cout << "// --- Barycenter" << std::endl;
      // TODO merge this in ttkMergeTreeVisu class
      // Declare VTK arrays
      vtkSmartPointer<vtkUnstructuredGrid> vtkArcs2 = vtkSmartPointer<vtkUnstructuredGrid>::New();
      vtkSmartPointer<vtkPoints> points2 = vtkSmartPointer<vtkPoints>::New();
      vtkNew<vtkIntArray> criticalType2{};
      criticalType2->SetName("CriticalType");
      vtkNew<vtkFloatArray> percentMatch2{};
      percentMatch2->SetName("PercentMatchNode");
      vtkNew<vtkFloatArray> percentMatchArc2{};
      percentMatchArc2->SetName("PercentMatchArc");
      vtkNew<vtkIntArray> branchID2{};    
      branchID2->SetName("BranchID");
      vtkNew<vtkIntArray> branchNodeID2{};
      branchNodeID2->SetName("BranchNodeID");
      vtkNew<vtkIntArray> nodeID2{};
      nodeID2->SetName("NodeId");
      vtkNew<vtkIntArray> upNodeId2{};
      upNodeId2->SetName("upNodeId");
      vtkNew<vtkIntArray> downNodeId2{};
      downNodeId2->SetName("downNodeId");
      vtkNew<vtkFloatArray> scalar2{};
      scalar2->SetName("Scalar");
      vtkNew<vtkFloatArray> persistenceNode2{};
      persistenceNode2->SetName("Persistence");
      vtkNew<vtkFloatArray> persistenceArc2{};
      persistenceArc2->SetName("Persistence");
      vtkNew<vtkIntArray> clusterIDNode2{};
      clusterIDNode2->SetName("ClusterID");
      vtkNew<vtkIntArray> clusterIDArc2{};
      clusterIDArc2->SetName("ClusterID");
      vtkNew<vtkIntArray> isImportantPairsArc{};
      isImportantPairsArc->SetName("isImportantPair");
      vtkNew<vtkIntArray> isDummyNode{};
      isDummyNode->SetName("isDummyNode");
      vtkNew<vtkIntArray> isDummyArc{};
      isDummyArc->SetName("isDummyArc");
      vtkNew<vtkIntArray> isImportantPairsNode{};
      isImportantPairsNode->SetName("isImportantPair");
      
      // Declare internal arrays
      std::vector<std::vector<SimplexId>> nodeCorrBary(NumberOfBarycenters);
      std::vector<std::vector<float>> allBaryPercentMatch(NumberOfBarycenters);
      double refPersistence = getPersistence(barycenters[0]->tree, getRoot(barycenters[0]->tree), dataTypeInt);
    
      // Create barycenters output
      if(verbose > 0) std::cout << "// Create barycenters output" << std::endl;
      for(int c = 0; c < NumberOfBarycenters; ++c){
        ftm::FTMTree_MT *baryTree = barycenters[c]->tree;
        allBaryPercentMatch[c] = std::vector<float>(baryTree->getNumberOfNodes(), 100);
        nodeCorrBary[c] = std::vector<SimplexId>(baryTree->getNumberOfNodes());
        
        std::vector<float> layout;
        if(PlanarLayout){
          ttkMergeTreeVisu visuMakerPlanarLayout;
          visuMakerPlanarLayout.setBranchDecompositionPlanarLayout(BranchDecompositionPlanarLayout);
          visuMakerPlanarLayout.setBranchSpacing(BranchSpacing);
          visuMakerPlanarLayout.setRescaleTreesIndividually(RescaleTreesIndividually);
          visuMakerPlanarLayout.setImportantPairs(ImportantPairs);
          visuMakerPlanarLayout.setImportantPairsSpacing(ImportantPairsSpacing);
          visuMakerPlanarLayout.setNonImportantPairsSpacing(NonImportantPairsSpacing);
          visuMakerPlanarLayout.setNonImportantPairsProximity(NonImportantPairsProximity);
          layout = visuMakerPlanarLayout.treePlanarLayout(baryTree, dataTypeInt, allBaryBounds[c], 
                                                          refPersistence);
        }
        int cptBadIndex = 0;
        int cptNode = 0;
        std::vector<SimplexId> treeSimplexId(baryTree->getNumberOfNodes());
        std::vector<SimplexId> treeDummySimplexId(baryTree->getNumberOfNodes());
        std::vector<SimplexId> layoutCorr(baryTree->getNumberOfNodes());
        
        // m[i][j] contains the node in trees[j] matched to the node i in the barycenter
        std::vector<std::vector<ftm::idNode>> baryMatching(baryTree->getNumberOfNodes(),
                                                           std::vector<ftm::idNode>(numInputs, -1));
        for(int i = 0; i < outputMatchingBarycenter[c].size(); ++i)
          for(auto match : outputMatchingBarycenter[c][i])
            baryMatching[std::get<0>(match)][i] = std::get<1>(match);
      
        // Tree traversal
        if(verbose > 0) std::cout << "// Tree traversal" << std::endl;
        std::queue<idNode> queue;
        queue.push(getRoot(baryTree));
        int realNoNodes = getRealNumberOfNodes(baryTree);
        while(!queue.empty()){
          idNode node = queue.front();
          queue.pop();
        
          // Push children to the queue
          if(verbose > 0) std::cout << "// Push children to the queue" << std::endl;
          for(auto child : getChildren(baryTree, node))
            queue.push(child);
        
          // Get and insert point
          int index = -1;
          double point[3] = {0, 0, 0};
          double noMatched = 0.0;
          for(int j = 0; j < numInputs; ++j){
            if(baryMatching[node][j] != -1){
              int nodeMesh = trees1NodeCorrMesh[j][baryMatching[node][j]];
              double* pointTemp = treesNodes[j]->GetPoints()->GetPoint(nodeMesh);
              for(int k = 0; k < 3; ++k)
                point[k] += pointTemp[k] - ( k!=2 ? 0 : std::get<4>(allBounds[j]) );
              noMatched += 1;
              index = j;
            }
          }
          cptBadIndex += (index == -1) ? 1 : 0;
          for(int k = 0; k < 3; ++k)
            point[k] /= noMatched; 
          if(PlanarLayout){
            layoutCorr[node] = cptNode;
            point[0] = layout[cptNode];
            point[1] = layout[cptNode+1];
            point[2] = 0;
            cptNode+=2;
          }
          point[0] += clusterShift[c];
          
          // TODO too many dummy nodes are created
          bool dummyNode = PlanarLayout and not BranchDecompositionPlanarLayout
                                        and !isRoot(baryTree, node)
                                        /*and !isLeaf(baryTree, node) 
                                        and isBranchOrigin(baryTree, node)*/;
          if(dummyNode)
            treeDummySimplexId[node] = points2->InsertNextPoint(point); // will be modified
          SimplexId nextPointId = points2->InsertNextPoint(point);
          treeSimplexId[node] = nextPointId;
          nodeCorrBary[c][node] = nextPointId;
          if(dummyNode)
            nodeCorrBary[c][node] = treeDummySimplexId[node];
          
          // Add cell connecting parent
          if(verbose > 0) std::cout << "// Add cell connecting parent" << std::endl;
          if(!isRoot(baryTree, node)){
            vtkIdType pointIds[2];
            pointIds[0] = treeSimplexId[node];
            
            ftm::idNode nodeParent = getParent(baryTree, node);
            // TODO too many dummy cells are created
            bool dummyCell = PlanarLayout and not BranchDecompositionPlanarLayout
                                          and allBaryBranching[c][node] == nodeParent 
                                          and !isRoot(baryTree, nodeParent);
            if(PlanarLayout and BranchDecompositionPlanarLayout){
              pointIds[1] = treeSimplexId[allBaryBranching[c][node]];
            }else if(dummyCell){
              double dummyPoint[3] = {point[0], layout[layoutCorr[nodeParent]+1], 0.};
              SimplexId dummyPointId = treeDummySimplexId[nodeParent];
              points2->SetPoint(dummyPointId, dummyPoint);
              vtkIdType dummyPointIds[2];
              dummyPointIds[0] = dummyPointId;
              dummyPointIds[1] = treeSimplexId[nodeParent];
              vtkArcs2->InsertNextCell(VTK_LINE, 2, dummyPointIds);
              pointIds[1] = dummyPointId;
            }else
              pointIds[1] = treeSimplexId[nodeParent];
            
            vtkArcs2->InsertNextCell(VTK_LINE, 2, pointIds);
          
            // --------------
            // Arc field
            // --------------
            int toAdd = (dummyCell ? 2 : 1);
            for(int toAddT = 0; toAddT < toAdd; ++toAddT){
              // Add arc matching percentage
              percentMatchArc2->InsertNextTuple1(allBaryPercentMatch[c][allBaryBranching[c][node]]);
              
              // Add branch id
              branchID2->InsertNextTuple1(allBaryBranchingID[c][node]);
            
              // Add up and down nodeId
              upNodeId2->InsertNextTuple1(treeSimplexId[nodeParent]);
              downNodeId2->InsertNextTuple1(treeSimplexId[node]);
              
              // Add arc persistence
              idNode nodeToGetPers = allBaryBranching[c][node];
              if(PlanarLayout and BranchDecompositionPlanarLayout)
                nodeToGetPers = node;
              double persToAdd = getPersistence(baryTree, nodeToGetPers , dataTypeInt);
              persistenceArc2->InsertNextTuple1(persToAdd);
              /*auto pers = getPersistence(baryTree, allBaryBranching[c][node], dataTypeInt);
              persistenceArc2->InsertNextTuple1(pers);*/
              
              // Add node clusterID
              clusterIDArc2->InsertNextTuple1(c);
              
              // Add isImportantPair
              bool isImportant;
              idNode nodeToGetImportance = allBaryBranching[c][node];
              if(PlanarLayout and BranchDecompositionPlanarLayout)
                nodeToGetImportance = node;
              switch(dataTypeInt){
                vtkTemplateMacro(isImportant = isImportantPair<VTK_TT>(baryTree, 
                                                                       nodeToGetImportance, ImportantPairs));
              }
              isImportantPairsArc->InsertNextTuple1(isImportant);
            
              // Add isDummyArc
              bool isDummy = toAdd == 2 and toAddT == 0;
              isDummyArc->InsertNextTuple1(isDummy);
            }
          }
          
          // --------------
          // Node field
          // --------------
          int toAdd = (dummyNode ? 2 : 1);
          for(int toAddT = 0; toAddT < toAdd; ++toAddT){
            // Add node id
            nodeID2->InsertNextTuple1(treeSimplexId[node]);
            
            // Add node scalar
            scalar2->InsertNextTuple1(getValue(baryTree, node, dataTypeInt));
            
            // Add criticalType
            int criticalTypeT = -1;
            if(index != -1){
              int nodeMesh = trees1NodeCorrMesh[index][baryMatching[node][index]];
              criticalTypeT = treesNodes[index]->GetPointData()->GetArray("CriticalType")->GetTuple1(nodeMesh);
            }
            criticalType2->InsertNextTuple1(criticalTypeT);
            
            // Add node matching percentage
            float percentMatchT = noMatched * 100 / numInputs;
            percentMatch2->InsertNextTuple1(percentMatchT);
            allBaryPercentMatch[c][node] = percentMatchT;
            
            // Add node branch id
            int tBranchID = allBaryBranchingID[c][node];
            if(not isLeaf(baryTree, node))
              tBranchID = allBaryBranchingID[c][baryTree->getNode(node)->getOrigin()];
            branchNodeID2->InsertNextTuple1(tBranchID);
            
            // Add node persistence
            persistenceNode2->InsertNextTuple1(getPersistence(baryTree, node, dataTypeInt));
            
            // Add node clusterID
            clusterIDNode2->InsertNextTuple1(c);
            
            // Add isDummyNode
            bool isDummy = toAdd == 2 and toAddT == 1 and !isRoot(baryTree, node);
            isDummyNode->InsertNextTuple1(isDummy);
            
            // Add isImportantPair
            bool isImportant;
            switch(dataTypeInt){
              vtkTemplateMacro(isImportant = isImportantPair<VTK_TT>(baryTree, node, ImportantPairs));
            }
            isImportantPairsNode->InsertNextTuple1(isImportant);
          }
        }
        if(cptBadIndex != 0)
          std::cout << "index -1 : " << cptBadIndex << std::endl;
      }
      vtkOutputNode2->SetPoints(points2);
      vtkOutputNode2->GetPointData()->AddArray(criticalType2);
      vtkOutputNode2->GetPointData()->AddArray(percentMatch2);
      vtkOutputNode2->GetPointData()->AddArray(branchNodeID2);
      vtkOutputNode2->GetPointData()->AddArray(nodeID2);
      vtkOutputNode2->GetPointData()->AddArray(scalar2);
      vtkOutputNode2->GetPointData()->AddArray(persistenceNode2);
      vtkOutputNode2->GetPointData()->AddArray(clusterIDNode2);
      vtkOutputNode2->GetPointData()->AddArray(isDummyNode);
      vtkOutputNode2->GetPointData()->AddArray(isImportantPairsNode);
      vtkArcs2->SetPoints(points2);
      vtkArcs2->GetCellData()->AddArray(percentMatchArc2);
      vtkArcs2->GetCellData()->AddArray(branchID2);
      vtkArcs2->GetCellData()->AddArray(upNodeId2);
      vtkArcs2->GetCellData()->AddArray(downNodeId2);
      vtkArcs2->GetCellData()->AddArray(persistenceArc2);
      vtkArcs2->GetCellData()->AddArray(clusterIDArc2);
      vtkArcs2->GetCellData()->AddArray(isDummyArc);
      vtkArcs2->GetCellData()->AddArray(isImportantPairsArc);
      vtkOutputArc2->ShallowCopy(vtkArcs2);
    
      // ------------------------------------------
      // --- Matching
      // ------------------------------------------
      vtkSmartPointer<vtkUnstructuredGrid> vtkMatching = vtkSmartPointer<vtkUnstructuredGrid>::New();
      vtkSmartPointer<vtkPoints> pointsM = vtkSmartPointer<vtkPoints>::New();
      vtkNew<vtkFloatArray> matchingPercentMatch{};
      matchingPercentMatch->SetName("MatchingPercentMatch");
      for(int c = 0; c < NumberOfBarycenters; ++c){
        for(int i = 0; i < numInputs; ++i){
          for(std::tuple<idNode, idNode> match : outputMatchingBarycenter[c][i]){
            vtkIdType pointIds[2];
          
            // Get first point
            double* point1 = vtkOutputNode2->GetPoints()->GetPoint(nodeCorrBary[c][std::get<0>(match)]);
            const SimplexId nextPointId1 = pointsM->InsertNextPoint(point1);
            pointIds[0] = nextPointId1;
          
            // Get second point
            double* point2 = vtkOutputNode1->GetPoints()->GetPoint(nodeCorr[i][std::get<1>(match)]);
            const SimplexId nextPointId2 = pointsM->InsertNextPoint(point2);
            pointIds[1] = nextPointId2;
          
            // Add cell
            vtkMatching->InsertNextCell(VTK_LINE, 2, pointIds);
          
            //Add arc matching percentage
            matchingPercentMatch->InsertNextTuple1(allBaryPercentMatch[c][std::get<0>(match)]);
          }
        }
      }
      vtkMatching->SetPoints(pointsM);
      vtkMatching->GetCellData()->AddArray(matchingPercentMatch);
      vtkOutputMatching->ShallowCopy(vtkMatching);
    }
  }
  
  if(Verbose > 0)
    std::cout << "TIME TREES OUT. = " << t_makeTreesOutput.getElapsedTime() << std::endl;
  
  // ------------------------------------------------------------------------------------
  // --- Free memory
  // ------------------------------------------------------------------------------------
  for(int i = 0; i < intermediateMTrees.size(); ++i){
    switch(dataTypeInt){
      vtkTemplateMacro(freeMergeTree<VTK_TT>(intermediateMTrees[i], false));
    }
  }
  for(int i = 0; i < intermediateMTrees2.size(); ++i)
    switch(dataTypeInt){
      vtkTemplateMacro(freeMergeTree<VTK_TT>(intermediateMTrees2[i], false));
    }
  
  return 1;
}
