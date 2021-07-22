#include <ttkFTMTreeEditDistance.h>
#include <ttkFTMTreeUtils.h>
#include <ttkUtils.h>
#include <FTMTree.h>
#include <FTMStructures.h>

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

#include <FTMTreeUtils.h>

using namespace ttk;
using namespace ftm;

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkFTMTreeEditDistance);

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
ttkFTMTreeEditDistance::ttkFTMTreeEditDistance() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(8);
}

ttkFTMTreeEditDistance::~ttkFTMTreeEditDistance() {
}

/**
 * Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkFTMTreeEditDistance::FillInputPortInformation(int port, vtkInformation *info) {
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
int ttkFTMTreeEditDistance::FillOutputPortInformation(int port, vtkInformation *info) {
  switch(port){
    case 0:
    case 1:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      break;
    case 2:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
      break;
    case 3:
    case 4:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      break;
    case 5:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
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
    treesNodes[i] = vtkUnstructuredGrid::SafeDownCast(inputTrees[i]->GetBlock(0));
    treesArcs[i] = vtkUnstructuredGrid::SafeDownCast(inputTrees[i]->GetBlock(1));
    treesSegmentation[i] = vtkDataSet::SafeDownCast(inputTrees[i]->GetBlock(2));
    intermediateTrees[i] = makeTree(treesNodes[i], treesArcs[i]);
  }
}

int ttkFTMTreeEditDistance::RequestData(vtkInformation *request,
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {  
  // ---------------------
  // --- Get input object from input vector
  // ---------------------
  vtkUnstructuredGrid *tree1Nodes; //= vtkUnstructuredGrid::GetData(inputVector[0]);
  vtkUnstructuredGrid *tree1Arcs; //= vtkUnstructuredGrid::GetData(inputVector[1]);
  vtkDataSet *tree1Segmentation; //= vtkDataSet::GetData(inputVector[2]);
  
  vtkUnstructuredGrid *tree2Nodes; //= vtkUnstructuredGrid::GetData(inputVector[3]);
  vtkUnstructuredGrid *tree2Arcs; //= vtkUnstructuredGrid::GetData(inputVector[4]);
  vtkDataSet *tree2Segmentation; //= vtkDataSet::GetData(inputVector[5]);
  
  // ---------------------
  // --- Construct trees
  // ---------------------
  FTMTree_MT *tree1; //= makeTree(tree1Nodes, tree1Arcs);  
  //tree1->printTree2();
  
  FTMTree_MT *tree2; //= makeTree(tree2Nodes, tree2Arcs);
  //tree2->printTree2();

  // ---------------------
  // --- Get input object from input vector
  // ---------------------
  auto blocks = vtkMultiBlockDataSet::GetData(inputVector[0], 0);
  auto blocks2 = vtkMultiBlockDataSet::GetData(inputVector[1], 0);

  // ---------------------
  // --- Construct trees
  // ---------------------
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
  
  int dataType = treesNodes[0]->GetPointData()->GetArray("Scalar")->GetDataType();
  
  tree1 = intermediateTrees[0];
  tree1Nodes = treesNodes[0];
  tree1Arcs = treesArcs[0];
  tree1Segmentation = treesSegmentation[0];
  tree2 = intermediateTrees[1];
  tree2Nodes = treesNodes[1];
  tree2Arcs = treesArcs[1];
  tree2Segmentation = treesSegmentation[1];
  
  // ---------------------
  // --- Call base
  // ---------------------
  std::vector<std::tuple<idNode, idNode>> outputMatching;
  std::vector<std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>>> 
         outputMatchingBarycenter(NumberOfBarycenters, 
                                  std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>>(numInputs)),
         outputMatchingBarycenter2(NumberOfBarycenters, 
                                  std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>>(numInputs));
  std::vector<MergeTree*> barycenters(NumberOfBarycenters);
  std::vector<int> clusteringAssignment(numInputs, 0);
  bool AddNodes=true;
  bool Deterministic=true;
  
  /*ComputeBarycenter=true;
  BranchDecomposition=true;
  NormalizedWasserstein=true;
  RescaledWasserstein=true;
  KeepSubtree=false;
  NumberOfBarycenters=2;*/
  //AddNodes=false;
  //RescaledWasserstein=false;
  Deterministic=false;
  
  //OutputTrees = not (ComputeBarycenter and NumberOfBarycenters > 1); // TODO fix clustering 
  
  if(not ComputeBarycenter){
    FTMTreeEditDistance ftmTreeEditDistance;
    ftmTreeEditDistance.setAssignmentSolver(AssignmentSolver);
    ftmTreeEditDistance.setEpsilonTree1(EpsilonTree1);
    ftmTreeEditDistance.setEpsilonTree2(EpsilonTree2);
    ftmTreeEditDistance.setEpsilon2Tree1(Epsilon2Tree1);
    ftmTreeEditDistance.setEpsilon2Tree2(Epsilon2Tree2);
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
    switch(dataType){
      vtkTemplateMacro(
        auto distance = ftmTreeEditDistance.execute<VTK_TT>(tree1, tree2, outputMatching)
      );
    }
  }else{
    if(NumberOfBarycenters == 1){ //and numInputs2==0){
      FTMTreeEditDistanceBarycenter ftmTreeEditDistanceBary;
      ftmTreeEditDistanceBary.setAssignmentSolver(AssignmentSolver);
      ftmTreeEditDistanceBary.setEpsilonTree1(EpsilonTree1);
      ftmTreeEditDistanceBary.setEpsilonTree2(EpsilonTree2);
      ftmTreeEditDistanceBary.setEpsilon2Tree1(Epsilon2Tree1);
      ftmTreeEditDistanceBary.setEpsilon2Tree2(Epsilon2Tree2);
      ftmTreeEditDistanceBary.setProgressiveComputation(ProgressiveComputation);
      ftmTreeEditDistanceBary.setBranchDecomposition(BranchDecomposition);
      ftmTreeEditDistanceBary.setParallelize(Parallelize);
      ftmTreeEditDistanceBary.setPersistenceThreshold(PersistenceThreshold);
      ftmTreeEditDistanceBary.setNormalizedWasserstein(NormalizedWasserstein);
      ftmTreeEditDistanceBary.setNormalizedWassersteinReg(NormalizedWassersteinReg);
      ftmTreeEditDistanceBary.setRescaledWasserstein(RescaledWasserstein);
      ftmTreeEditDistanceBary.setKeepSubtree(KeepSubtree);
      ftmTreeEditDistanceBary.setUseMinMaxPair(UseMinMaxPair);
      ftmTreeEditDistanceBary.setTol(Tol);
      ftmTreeEditDistanceBary.setAddNodes(AddNodes);
      ftmTreeEditDistanceBary.setNumberOfThreads(NumberOfThreads);
      switch(dataType){
        vtkTemplateMacro(
          barycenters[0] = ftmTreeEditDistanceBary.execute<VTK_TT>(intermediateTrees, 
                                                                   outputMatchingBarycenter[0])
        );
      }
    }else{
      FTMTreeEditDistanceClustering ftmTreeEditDistanceClust;
      ftmTreeEditDistanceClust.setAssignmentSolver(AssignmentSolver);
      ftmTreeEditDistanceClust.setEpsilonTree1(EpsilonTree1);
      ftmTreeEditDistanceClust.setEpsilonTree2(EpsilonTree2);
      ftmTreeEditDistanceClust.setEpsilon2Tree1(Epsilon2Tree1);
      ftmTreeEditDistanceClust.setEpsilon2Tree2(Epsilon2Tree2);
      ftmTreeEditDistanceClust.setProgressiveComputation(ProgressiveComputation);
      ftmTreeEditDistanceClust.setBranchDecomposition(BranchDecomposition);
      ftmTreeEditDistanceClust.setParallelize(Parallelize);
      ftmTreeEditDistanceClust.setPersistenceThreshold(PersistenceThreshold);
      ftmTreeEditDistanceClust.setNormalizedWasserstein(NormalizedWasserstein);
      ftmTreeEditDistanceClust.setNormalizedWassersteinReg(NormalizedWassersteinReg);
      ftmTreeEditDistanceClust.setRescaledWasserstein(RescaledWasserstein);
      ftmTreeEditDistanceClust.setKeepSubtree(KeepSubtree);
      ftmTreeEditDistanceClust.setUseMinMaxPair(UseMinMaxPair);
      ftmTreeEditDistanceClust.setTol(Tol);
      ftmTreeEditDistanceClust.setAddNodes(AddNodes);
      ftmTreeEditDistanceClust.setNoCentroids(NumberOfBarycenters);
      ftmTreeEditDistanceClust.setDeterministic(Deterministic);
      ftmTreeEditDistanceClust.setNumberOfThreads(NumberOfThreads);
      switch(dataType){
        vtkTemplateMacro(
          barycenters = ftmTreeEditDistanceClust.execute<VTK_TT>(intermediateTrees, outputMatchingBarycenter,
                                                               clusteringAssignment,
                                                               intermediateTrees2, outputMatchingBarycenter2)
        );
      }
    }
  }
  
  // ---------------------
  // --- Create output
  // ---------------------
  auto vtkOutputNode1 = vtkUnstructuredGrid::SafeDownCast(
    outputVector->GetInformationObject(0)->Get(vtkDataObject::DATA_OBJECT()));
  
  auto vtkOutputArc1 = vtkUnstructuredGrid::SafeDownCast(
    outputVector->GetInformationObject(1)->Get(vtkDataObject::DATA_OBJECT()));
  
  auto vtkOutputSegmentation1 = vtkDataSet::SafeDownCast(
    outputVector->GetInformationObject(2)->Get(vtkDataObject::DATA_OBJECT()));
  
  auto vtkOutputNode2 = vtkUnstructuredGrid::SafeDownCast(
      outputVector->GetInformationObject(3)->Get(vtkDataObject::DATA_OBJECT()));
  
  auto vtkOutputArc2 = vtkUnstructuredGrid::SafeDownCast(
    outputVector->GetInformationObject(4)->Get(vtkDataObject::DATA_OBJECT()));
  
  auto vtkOutputSegmentation2 = vtkDataSet::SafeDownCast(
    outputVector->GetInformationObject(5)->Get(vtkDataObject::DATA_OBJECT()));
  
  auto vtkOutputMatching = vtkUnstructuredGrid::SafeDownCast(
    outputVector->GetInformationObject(6)->Get(vtkDataObject::DATA_OBJECT()));
  
  auto vtkAssignment = vtkTable::SafeDownCast(
    outputVector->GetInformationObject(7)->Get(vtkDataObject::DATA_OBJECT()));
  
  if(not ComputeBarycenter){
    // TODO merge code with barycenter output
    // -------
    // - Classical output
    // -------    
    vtkOutputNode1->DeepCopy(tree1Nodes); // will (perhaps) be modified
    vtkOutputArc1->ShallowCopy(tree1Arcs);
    vtkOutputSegmentation1->ShallowCopy(tree1Segmentation);
    vtkOutputNode2->DeepCopy(tree2Nodes); // will be modified
    vtkOutputArc2->DeepCopy(tree2Arcs); // will be modified
    if(OutputSegmentation)
      vtkOutputSegmentation2->DeepCopy(tree2Segmentation); // will be modified

    // Shift
    double diffTree1_z = (vtkOutputNode1->GetBounds()[5] - vtkOutputNode1->GetBounds()[4])*2;
    double diff_z = (vtkOutputNode1->GetBounds()[4] - vtkOutputNode2->GetBounds()[4]);
    diff_z += diffTree1_z;
    
    // Shift nodes
    for(int i = 0; i < vtkOutputNode2->GetPoints()->GetNumberOfPoints(); ++i){
      double* point = vtkOutputNode2->GetPoints()->GetPoint(i);
      vtkOutputNode2->GetPoints()->SetPoint(i, point[0], point[1], point[2]+diff_z);
    }
    
    // Shift arcs
    for(int i = 0; i < vtkOutputArc2->GetPoints()->GetNumberOfPoints(); ++i){
      double* point = vtkOutputArc2->GetPoints()->GetPoint(i);
      vtkOutputArc2->GetPoints()->SetPoint(i, point[0], point[1], point[2]+diff_z);
    }
    
    // Shift segmentation
    if(OutputSegmentation){
      auto vtkOutputSegmentation2Temp = vtkUnstructuredGrid::SafeDownCast(vtkOutputSegmentation2);
      for(int i = 0; i < vtkOutputSegmentation2Temp->GetPoints()->GetNumberOfPoints(); ++i){
        double* point = vtkOutputSegmentation2Temp->GetPoints()->GetPoint(i);
        vtkOutputSegmentation2Temp->GetPoints()->SetPoint(i, point[0], point[1], point[2]+diff_z);
      }
      vtkOutputSegmentation2 = vtkDataSet::SafeDownCast(vtkOutputSegmentation2Temp);
    }
    
    // Create matching
    vtkNew<vtkIntArray> matchingType{};
    matchingType->SetName("MatchingType");
    
    vtkSmartPointer<vtkUnstructuredGrid> vtkMatching = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    int count = 0;
    for(std::tuple<idNode, idNode> match : outputMatching){
      vtkIdType pointIds[2];
      // Get first point
      double* point1 = vtkOutputNode1->GetPoints()->GetPoint(std::get<0>(match));
      const SimplexId nextPointId1 = points->InsertNextPoint(point1);
      pointIds[0] = nextPointId1;
      // Get second point
      double* point2 = vtkOutputNode2->GetPoints()->GetPoint(std::get<1>(match));
      const SimplexId nextPointId2 = points->InsertNextPoint(point2);
      pointIds[1] = nextPointId2;
      // Add cell
      vtkMatching->InsertNextCell(VTK_LINE, 2, pointIds);
      
      // Manage matching type
      int thisType = 0;
      int tree1NodeDown = tree1->getNode(std::get<0>(match))->getNumberOfDownSuperArcs();
      int tree1NodeUp = tree1->getNode(std::get<0>(match))->getNumberOfUpSuperArcs();
      int tree2NodeDown = tree2->getNode(std::get<1>(match))->getNumberOfDownSuperArcs();
      int tree2NodeUp = tree2->getNode(std::get<1>(match))->getNumberOfUpSuperArcs();
      if(tree1NodeDown != 0 and tree1NodeUp != 0 and tree2NodeDown != 0 and tree2NodeUp != 0)
        thisType = 1; // Saddle to Saddle
      if(tree1NodeDown == 0 and tree1NodeUp != 0 and tree2NodeDown == 0 and tree2NodeUp != 0)
        thisType = 2; // Leaf to leaf
      if(tree1NodeDown != 0 and tree1NodeUp == 0 and tree2NodeDown != 0 and tree2NodeUp == 0)
        thisType = 3; // Root to root
      matchingType->InsertTuple1(count, thisType);
      count++;
    }
    vtkMatching->SetPoints(points);
    vtkMatching->GetCellData()->AddArray(matchingType);
    vtkOutputMatching->ShallowCopy(vtkMatching);
    
    // Remove points of epsilon merging
    removePointsEpsilon(vtkOutputNode1, tree1);  
    removePointsEpsilon(vtkOutputNode2, tree2);    
  }else {
    // -------
    // - Barycenter output
    // -------
    // - Assignment Table
    vtkNew<vtkDoubleArray> col{};
    col->SetNumberOfTuples(numInputs);
    col->SetName("ClusterAssignment");
    for(int i = 0; i < numInputs; ++i)
      col->SetTuple1(i, clusteringAssignment[i]);
    vtkAssignment->AddColumn(col);

   if(OutputTrees){
    // --- Input Trees
    // Declare VTK arrays
    auto vtkOutputSegmentation1Temp = vtkUnstructuredGrid::SafeDownCast(vtkOutputSegmentation1);
    vtkSmartPointer<vtkUnstructuredGrid> vtkArcs = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkNew<vtkIntArray> criticalType{};
    criticalType->SetName("CriticalType");
    vtkNew<vtkIntArray> branchID{};
    branchID->SetName("BranchID");
    vtkNew<vtkIntArray> branchNodeID{};
    branchNodeID->SetName("BranchNodeID");
    vtkNew<vtkFloatArray> percentMatch{};
    percentMatch->SetName("PercentMatchNode");
    vtkNew<vtkFloatArray> percentMatchArc{};
    percentMatchArc->SetName("PercentMatchArc");
    
    // Declare internal arrays and initialize them
    std::vector<double> clusterShift(NumberOfBarycenters, 0);
    std::vector<std::vector<SimplexId>> nodeCorr(numInputs);
    std::vector<std::tuple<double,double,double,double,double,double>> allBounds(numInputs);
    std::vector<std::tuple<double,double,double,double,double,double>> allBaryBounds(NumberOfBarycenters);
    std::vector<std::vector<ftm::idNode>> allBaryBranching(NumberOfBarycenters);
    std::vector<std::vector<int>> allBaryBranchingID(NumberOfBarycenters);
    for(int i = 0; i < numInputs; ++i)
      allBounds[i] = getRealBounds(treesNodes[i], intermediateTrees[i]);
    for(int c = 0; c < NumberOfBarycenters; ++c){
      allBaryBounds[c] = getMaximalBounds(allBounds, clusteringAssignment, c);
      getTreeBranching(barycenters[c]->tree, allBaryBranching[c], allBaryBranchingID[c]);
    }
      
    // Create trees output
    for(int c = 0; c < NumberOfBarycenters; ++c){
      double delta_max = std::numeric_limits<double>::min();
      int noSample = 0;
      for(int i = 0; i < numInputs; ++i){
        if(clusteringAssignment[i] != c)
          continue;
        delta_max = std::max((std::get<3>(allBounds[i]) - std::get<2>(allBounds[i])), delta_max);
        delta_max = std::max((std::get<1>(allBounds[i]) - std::get<0>(allBounds[i])), delta_max);
        noSample += 1;
      }
      double radius = delta_max*2;
      double pi = 3.14159265359;
      int iSample = 0;
      for(int i = 0; i < numInputs; ++i){
        if(clusteringAssignment[i] != c)
          continue;
        nodeCorr[i] = std::vector<SimplexId>(intermediateTrees[i]->getNumberOfNodes());
        double angle = 360.0 / noSample * iSample;
        iSample += 1;
        double diff_x = radius * std::cos(-1*angle * pi / 180) + clusterShift[c];
        double diff_y = radius * std::sin(-1*angle * pi / 180);
        double diff_z = - std::get<4>(allBounds[i]);
        
        std::vector<float> layout;
        if(PlanarLayout)
          //layout = treePlanarLayout(intermediateTrees[i], dataType, allBounds[i]);
          layout = treePlanarLayout(intermediateTrees[i], dataType, allBaryBounds[c]);
        int cptNode = 0;
        std::vector<SimplexId> treeSimplexId(intermediateTrees[i]->getNumberOfNodes());
        std::vector<ftm::idNode> treeMatching(intermediateTrees[i]->getNumberOfNodes(), -1);
        for(auto match : outputMatchingBarycenter[c][i])
          treeMatching[std::get<1>(match)] = std::get<0>(match);
        
        // Tree traversal
        std::queue<idNode> queue;
        queue.push(getRoot(intermediateTrees[i]));
        while(!queue.empty()){
          idNode node = queue.front();
          queue.pop();
          
          // Get and insert point
          double* pointP = treesNodes[i]->GetPoints()->GetPoint(node);
          double point[3] = {pointP[0], pointP[1], pointP[2]};
          if(PlanarLayout){
            point[0] = layout[cptNode];
            point[1] = layout[cptNode+1];
            point[2] = 0;
            cptNode+=2;
          }else
            point[2] += diff_z;
          point[0] += diff_x;
          point[1] += diff_y;
          SimplexId nextPointId = points->InsertNextPoint(point);
          treeSimplexId[node] = nextPointId;
          nodeCorr[i][node] = nextPointId;
          
          // Add cell connecting parent
          if(!isRoot(intermediateTrees[i], node)){
            vtkIdType pointIds[2];
            pointIds[0] = treeSimplexId[node];
            pointIds[1] = treeSimplexId[getParent(intermediateTrees[i], node)];
            vtkArcs->InsertNextCell(VTK_LINE, 2, pointIds);
            
            // Add arc matching percentage
            //percentMatchArc->InsertNextTuple1(allBaryPercentMatch[c][allBaryBranching[c][treeMatching[node]]]);
            
            // Add branch ID
            branchID->InsertNextTuple1(allBaryBranchingID[c][treeMatching[node]]);
          }
          
          // Push children to the queue
          for(auto child : getChildren(intermediateTrees[i], node))
            queue.push(child);
          
          // Add criticalType
          int criticalTypeT = treesNodes[i]->GetPointData()->GetArray("CriticalType")->GetTuple1(node);
          criticalType->InsertNextTuple1(criticalTypeT);
          
          // Add node branch id
          int tBranchID = allBaryBranchingID[c][treeMatching[node]];
          if(not isLeaf(intermediateTrees[i], node))
            tBranchID=allBaryBranchingID[c][treeMatching[intermediateTrees[i]->getNode(node)->getOrigin()]];
          branchNodeID->InsertNextTuple1(tBranchID);
        }
        
        // Shift segmentation
        if(OutputSegmentation){
          auto iVkOutputSegmentationTemp = vtkUnstructuredGrid::SafeDownCast(treesSegmentation[i]);
          for(int j = 0; j < iVkOutputSegmentationTemp->GetPoints()->GetNumberOfPoints(); ++j){
            double* point = iVkOutputSegmentationTemp->GetPoints()->GetPoint(i);
            point[1] += diff_y;
            vtkOutputSegmentation1Temp->GetPoints()->InsertNextPoint(point);
          }
        }
      }
      if(c < NumberOfBarycenters-1)
        clusterShift[c+1] += radius*4;
    }
    
    // Add VTK arrays to output
    vtkOutputNode1->SetPoints(points);
    vtkOutputNode1->GetPointData()->AddArray(criticalType);
    vtkOutputNode1->GetPointData()->AddArray(branchNodeID);
    vtkArcs->SetPoints(points);
    vtkArcs->GetCellData()->AddArray(branchID);
    vtkOutputArc1->ShallowCopy(vtkArcs);
    vtkOutputSegmentation1 = vtkDataSet::SafeDownCast(vtkOutputSegmentation1Temp);
    
    // --- Barycenter
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
    
    // Declare internal arrays
    std::vector<std::vector<SimplexId>> nodeCorrBary(NumberOfBarycenters);
    std::vector<std::vector<float>> allBaryPercentMatch(NumberOfBarycenters);
    
    // Create barycenters output
    for(int c = 0; c < NumberOfBarycenters; ++c){
      ftm::FTMTree_MT *baryTree = barycenters[c]->tree;
      allBaryPercentMatch[c] = std::vector<float>(baryTree->getNumberOfNodes(), 100);
      nodeCorrBary[c] = std::vector<SimplexId>(baryTree->getNumberOfNodes());
      std::vector<SimplexId> treeSimplexId(baryTree->getNumberOfNodes());
      std::vector<float> layout;
      if(PlanarLayout)
        layout = treePlanarLayout(baryTree, dataType, allBaryBounds[c]);
      
      // m[i][j] contains the node in trees[j] matched to the node i in the barycenter
      std::vector<std::vector<ftm::idNode>> baryMatching(baryTree->getNumberOfNodes(),
                                                         std::vector<ftm::idNode>(numInputs, -1));
      for(int i = 0; i < outputMatchingBarycenter[c].size(); ++i)
        for(auto match : outputMatchingBarycenter[c][i])
          baryMatching[std::get<0>(match)][i] = std::get<1>(match);
      int cptBadIndex = 0;
      int cptNode = 0;
      
      // Tree traversal
      std::queue<idNode> queue;
      queue.push(getRoot(baryTree));
      while(!queue.empty()){
        idNode node = queue.front();
        queue.pop();
        
        // Get and insert point
        int index = -1;
        double point[3] = {0, 0, 0};
        double noMatched = 0.0;
        for(int j = 0; j < numInputs; ++j){
          if(baryMatching[node][j] != -1){
            double* pointTemp = treesNodes[j]->GetPoints()->GetPoint(baryMatching[node][j]);
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
          point[0] = layout[cptNode];
          point[1] = layout[cptNode+1];
          point[2] = 0;
          cptNode+=2;
        }
        point[0] += clusterShift[c];
        SimplexId nextPointId = points2->InsertNextPoint(point);
        treeSimplexId[node] = nextPointId;
        nodeCorrBary[c][node] = nextPointId;
        
        // Add cell connecting parent
        if(!isRoot(baryTree, node)){
          vtkIdType pointIds[2];
          pointIds[0] = treeSimplexId[node];
          pointIds[1] = treeSimplexId[getParent(baryTree, node)];
          vtkArcs2->InsertNextCell(VTK_LINE, 2, pointIds);
          
          // Add arc matching percentage
          percentMatchArc2->InsertNextTuple1(allBaryPercentMatch[c][allBaryBranching[c][node]]);
          
          // Add branch id
          branchID2->InsertNextTuple1(allBaryBranchingID[c][node]);
        }
        
        // Push children to the queue
        for(auto child : getChildren(baryTree, node))
          queue.push(child);
        
        // Add criticalType
        int criticalTypeT = (index == -1) ? -1 :
          treesNodes[index]->GetPointData()->GetArray("CriticalType")->GetTuple1(baryMatching[node][index]);
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
      }
      if(cptBadIndex != 0)
        std::cout << "index -1 : " << cptBadIndex << std::endl;
    }
    vtkOutputNode2->SetPoints(points2);
    vtkOutputNode2->GetPointData()->AddArray(criticalType2);
    vtkOutputNode2->GetPointData()->AddArray(percentMatch2);
    vtkOutputNode2->GetPointData()->AddArray(branchNodeID2);
    vtkArcs2->SetPoints(points2);
    vtkArcs2->GetCellData()->AddArray(percentMatchArc2);
    vtkArcs2->GetCellData()->AddArray(branchID2);
    vtkOutputArc2->ShallowCopy(vtkArcs2);

    // - Matching
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
  
  // ---------------------
  // --- Free memory
  // ---------------------
  for(int i = 0; i < intermediateMTrees.size(); ++i)
    switch(dataType){
      vtkTemplateMacro(freeMergeTree<VTK_TT>(intermediateMTrees[i], false));
    }
  for(int i = 0; i < intermediateMTrees2.size(); ++i)
    switch(dataType){
      vtkTemplateMacro(freeMergeTree<VTK_TT>(intermediateMTrees2[i], false));
    }
  
  return 1;
}
