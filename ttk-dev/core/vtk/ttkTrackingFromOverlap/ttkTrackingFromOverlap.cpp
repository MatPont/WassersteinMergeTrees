#include <ttkTrackingFromOverlap.h>

#include <vtkPointSet.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkUnstructuredGrid.h>

#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkLongLongArray.h>
#include <vtkPointData.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkTrackingFromOverlap)

  // Function to fetch data of a specificd block
  void getData(vtkMultiBlockDataSet *mb,
               size_t time,
               size_t level,
               string labelFieldName,
               vtkPointSet *&pointSet,
               vtkAbstractArray *&labels) {
  auto timesteps = vtkMultiBlockDataSet::SafeDownCast(mb->GetBlock(level));
  pointSet = vtkPointSet::SafeDownCast(timesteps->GetBlock(time));
  labels = pointSet->GetPointData()->GetAbstractArray(labelFieldName.data());
};

// Function to compute number of levels and timesteps contained in a
// vtkMultiBlockDataSet
void getNumberOfLevelsAndTimesteps(vtkMultiBlockDataSet *mb,
                                   size_t &nL,
                                   size_t &nT) {
  nL = mb->GetNumberOfBlocks();
  auto timesteps = vtkMultiBlockDataSet::SafeDownCast(mb->GetBlock(0));
  nT = timesteps->GetNumberOfBlocks();
};

// =============================================================================
// Finalize
// =============================================================================
template <typename labelType>
int finalize(vector<vector<Nodes>> &levelTimeNodesMap,
             vector<vector<Edges>> &levelTimeEdgesTMap,
             vector<vector<Edges>> &timeLevelEdgesNMap,
             int labelTypeId,
             string labelFieldName,

             vtkDataObject *trackingGraphObject) {
  auto trackingGraph = vtkUnstructuredGrid::SafeDownCast(trackingGraphObject);

  size_t nL = levelTimeNodesMap.size();
  size_t nT = levelTimeNodesMap[0].size();

  auto prepArray = [](vtkAbstractArray *array, string name, size_t nComponents,
                      size_t nValues) {
    array->SetName(name.data());
    array->SetNumberOfComponents(nComponents);
    array->SetNumberOfTuples(nValues);
  };

  // Add Points
  {
    size_t nNodes = 0;
    for(size_t t = 0; t < nT; t++)
      for(size_t l = 0; l < nL; l++)
        nNodes += levelTimeNodesMap[l][t].size();

    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(nNodes);
    auto pointCoords = (float *)points->GetVoidPointer(0);

    auto sequence = vtkSmartPointer<vtkLongLongArray>::New();
    prepArray(sequence, "SequenceIndex", 1, nNodes);
    auto sequenceData = (long long *)sequence->GetVoidPointer(0);

    auto level = vtkSmartPointer<vtkLongLongArray>::New();
    prepArray(level, "LevelIndex", 1, nNodes);
    auto levelData = (long long *)level->GetVoidPointer(0);

    auto size = vtkSmartPointer<vtkFloatArray>::New();
    prepArray(size, "Size", 1, nNodes);
    auto sizeData = (float *)size->GetVoidPointer(0);

    auto branch = vtkSmartPointer<vtkLongLongArray>::New();
    prepArray(branch, "BranchId", 1, nNodes);
    auto branchData = (long long *)branch->GetVoidPointer(0);

    auto label = vtkSmartPointer<vtkDataArray>::Take(
      vtkDataArray::CreateDataArray(labelTypeId));
    prepArray(label, labelFieldName, 1, nNodes);
    auto labelData = (labelType *)label->GetVoidPointer(0);

    size_t q1 = 0, q2 = 0;
    for(size_t t = 0; t < nT; t++) {
      for(size_t l = 0; l < nL; l++) {
        for(auto &node : levelTimeNodesMap[l][t]) {
          pointCoords[q1++] = node.x;
          pointCoords[q1++] = node.y;
          pointCoords[q1++] = node.z;

          sequenceData[q2] = t;
          levelData[q2] = l;
          sizeData[q2] = node.size;
          branchData[q2] = node.branchID;
          labelData[q2] = boost::get<labelType>(node.label);

          q2++;
        }
      }
    }

    trackingGraph->SetPoints(points);

    auto pointData = trackingGraph->GetPointData();
    pointData->AddArray(sequence);
    pointData->AddArray(level);
    pointData->AddArray(size);
    pointData->AddArray(label);
    pointData->AddArray(branch);
  }

  // Add Cells
  {
    // Build node index offset vector
    vector<size_t> timeLevelOffsetMap(nT * nL + 1);
    {
      timeLevelOffsetMap[0] = 0;
      size_t q = 1;
      for(size_t t = 0; t < nT; t++)
        for(size_t l = 0; l < nL; l++) {
          timeLevelOffsetMap[q]
            = timeLevelOffsetMap[q - 1] + levelTimeNodesMap[l][t].size();
          q++;
        }
    }

    size_t nEdgesT = 0;
    if(nT > 1)
      for(size_t t = 0; t < nT - 1; t++)
        for(size_t l = 0; l < nL; l++)
          nEdgesT += levelTimeEdgesTMap[l][t].size() / 4;

    size_t nEdgesN = 0;
    if(nL > 1)
      for(size_t l = 0; l < nL - 1; l++)
        for(size_t t = 0; t < nT; t++)
          nEdgesN += timeLevelEdgesNMap[t][l].size() / 4;

    auto cells = vtkSmartPointer<vtkIdTypeArray>::New();
    cells->SetNumberOfValues(3 * nEdgesT + 3 * nEdgesN);
    auto cellIds = (vtkIdType *)cells->GetVoidPointer(0);

    auto overlap = vtkSmartPointer<vtkFloatArray>::New();
    prepArray(overlap, "Overlap", 1, nEdgesT + nEdgesN);
    auto overlapData = (float *)overlap->GetVoidPointer(0);

    auto branch = vtkSmartPointer<vtkLongLongArray>::New();
    prepArray(branch, "BranchId", 1, nEdgesT + nEdgesN);
    auto branchData = (long long *)branch->GetVoidPointer(0);

    auto type = vtkSmartPointer<vtkCharArray>::New();
    prepArray(type, "Type", 1, nEdgesT + nEdgesN);
    auto typeData = (char *)type->GetVoidPointer(0);

    size_t q0 = 0, q1 = 0;

    // Tracking graphs
    if(nT > 1)
      for(size_t t = 1; t < nT; t++) {
        for(size_t l = 0; l < nL; l++) {
          auto &edges = levelTimeEdgesTMap[l][t - 1];
          for(size_t i = 0, j = edges.size(); i < j;) {
            cellIds[q0++] = 2;
            cellIds[q0++]
              = (vtkIdType)(timeLevelOffsetMap[(t - 1) * nL + l] + edges[i++]);
            cellIds[q0++]
              = (vtkIdType)(timeLevelOffsetMap[(t)*nL + l] + edges[i++]);
            typeData[q1] = 0;
            overlapData[q1] = edges[i++];
            branchData[q1] = edges[i++];
            q1++;
          }
        }
      }

    // Nesting trees
    if(nL > 1)
      for(size_t l = 1; l < nL; l++) {
        for(size_t t = 0; t < nT; t++) {
          auto &edges = timeLevelEdgesNMap[t][l - 1];
          size_t temp = t * nL;
          for(size_t i = 0, j = edges.size(); i < j;) {
            cellIds[q0++] = 2;
            cellIds[q0++]
              = (vtkIdType)(timeLevelOffsetMap[temp + (l - 1)] + edges[i++]);
            cellIds[q0++]
              = (vtkIdType)(timeLevelOffsetMap[temp + (l)] + edges[i++]);
            typeData[q1] = 1;
            overlapData[q1] = edges[i++];
            branchData[q1] = edges[i++];
            q1++;
          }
        }
      }

    auto cellArray = vtkSmartPointer<vtkCellArray>::New();
    cellArray->SetCells(nEdgesT + nEdgesN, cells);
    trackingGraph->SetCells(VTK_LINE, cellArray);

    auto cellData = trackingGraph->GetCellData();
    cellData->AddArray(type);
    cellData->AddArray(overlap);
    cellData->AddArray(branch);
  }

  return 1;
}

// =============================================================================
// Mesh Nested Tracking Graph
// =============================================================================
int ttkTrackingFromOverlap::meshNestedTrackingGraph(
  vtkDataObject *trackingGraph) {
  Timer t;

  dMsg(cout,
       "[ttkTrackingFromOverlap] "
       "=======================================================\n",
       infoMsg);
  dMsg(
    cout, "[ttkTrackingFromOverlap] Meshing nested tracking graph\n", infoMsg);

  switch(this->LabelDataType) {
    vtkTemplateMacro(
      finalize<VTK_TT>(this->levelTimeNodesMap, this->levelTimeEdgesTMap,
                       this->timeLevelEdgesNMap, this->LabelDataType,
                       this->GetLabelFieldName(), trackingGraph));
  }

  stringstream msg;
  msg << "[ttkTrackingFromOverlap] "
         "-------------------------------------------------------"
      << endl
      << "[ttkTrackingFromOverlap] Nested tracking graph meshed in "
      << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))."
      << endl;
  dMsg(cout, msg.str(), timeMsg);

  return 1;
}

// =============================================================================
// Reset
// =============================================================================
int ttkTrackingFromOverlap::reset() {
  this->LabelDataType = -1;

  this->levelTimeNodesMap.clear();
  this->levelTimeEdgesTMap.clear();
  this->timeLevelEdgesNMap.clear();

  this->previousIterationData = nullptr;

  return 1;
}

// =============================================================================
// Pack Data
// =============================================================================
int ttkTrackingFromOverlap::packInputData(
  vtkDataObject *inputDataObject, vtkMultiBlockDataSet *packedInput) const {
  /* Enforce following vtkMultiBlockDataSet structure:
      {
          level_0: {
              time_0: vtkPointSet,
              ...
              time_nT: vtkPointSet
          },
          ...
          level_nL: {
              time_0: vtkPointSet,
              ...
              time_nT: vtkPointSet
          }
      }
  */

  // Check inputDataObject depth: 2->(level->time), 1->(-,time), 0->(-,-)
  auto inputAsMB = vtkMultiBlockDataSet::SafeDownCast(inputDataObject);
  auto inputAsPS = vtkPointSet::SafeDownCast(inputDataObject);
  bool error = false;
  if(inputAsMB) {
    size_t n = inputAsMB->GetNumberOfBlocks();

    // Check if blocks are vtkPointSets or vtkMultiBlockDataSets ...
    size_t psCounter = 0;
    size_t mbCounter = 0;
    for(size_t i = 0; i < n; i++) {
      auto block = inputAsMB->GetBlock(i);
      auto blockAsMB = vtkMultiBlockDataSet::SafeDownCast(block);
      auto blockAsPS = vtkPointSet::SafeDownCast(block);
      if(blockAsMB)
        mbCounter++;
      if(blockAsPS)
        psCounter++;
    }

    if(mbCounter
       == n) { // input already contains timesteps per level -> nothing to do
      packedInput->ShallowCopy(inputAsMB);
    } else if(psCounter
              == n) { // if input is a single list of vtkPointSets over time
      auto level = vtkSmartPointer<vtkMultiBlockDataSet>::New();
      for(size_t i = 0; i < n; i++) {
        level->SetBlock(i, vtkPointSet::SafeDownCast(inputAsMB->GetBlock(i)));
      }
      packedInput->SetBlock(0, level);
    } else { // Unexpected input structure
      error = true;
    }
  } else if(inputAsPS) {
    auto level = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    level->SetBlock(0, inputAsPS);
    packedInput->SetBlock(0, level);
  } else { // Unexpected input structure
    error = true;
  }

  if(error) {
    dMsg(cerr,
         "[ttkTrackingFromOverlap] ERROR: Unable to convert input into "
         "'vtkPointSet' collection.\n",
         fatalMsg);
    return 0;
  }

  return 1;
}

// =============================================================================
// Check Data
// =============================================================================
int ttkTrackingFromOverlap::checkData(vtkMultiBlockDataSet *data) {
  size_t nL = data->GetNumberOfBlocks();
  size_t nT = 0;

  if(nL < 1) {
    dMsg(cerr,
         "[ttkTrackingFromOverlap] ERROR: Input must have at least one "
         "vtkPointSet.\n",
         fatalMsg);
    return 0;
  }

  for(size_t l = 0; l < nL; l++) {
    auto timesteps = vtkMultiBlockDataSet::SafeDownCast(data->GetBlock(l));
    size_t n = timesteps->GetNumberOfBlocks();
    if(n < 1) {
      dMsg(cerr,
           "[ttkTrackingFromOverlap] ERROR: Input must have at least one "
           "vtkPointSet.\n",
           fatalMsg);
      return 0;
    }
    if(nT == 0)
      nT = n;
    if(nT != n) {
      dMsg(cerr,
           "[ttkTrackingFromOverlap] ERROR: Timeseries have unequal length.\n",
           fatalMsg);
      return 0;
    }

    for(size_t t = 0; t < nT; t++) {
      auto pointSet = vtkPointSet::SafeDownCast(timesteps->GetBlock(t));
      size_t nPoints = pointSet->GetNumberOfPoints();
      auto labels = pointSet->GetPointData()->GetAbstractArray(
        this->GetLabelFieldName().data());

      if(nPoints > 0 && labels == nullptr) {
        dMsg(cout,
             "[ttkTrackingFromOverlap] ERROR: Point labels '"
               + this->GetLabelFieldName() + "' not found.\n",
             fatalMsg);
        return 0;
      }
      if(labels == nullptr)
        continue;

      int labelDataType = labels->GetDataType();
      if(this->LabelDataType < 0)
        this->LabelDataType = labelDataType;
      if(this->LabelDataType != labelDataType) {
        dMsg(cout,
             "[ttkTrackingFromOverlap] ERROR: Point labels do not have same "
             "type across point sets.\n",
             fatalMsg);
        return 0;
      }
    }
  }

  return 1;
}

// =============================================================================
// Pack Streamed Data
// =============================================================================
int ttkTrackingFromOverlap::packStreamedData(vtkMultiBlockDataSet *streamedData,
                                             vtkMultiBlockDataSet *data) const {
  size_t nL_PI = this->previousIterationData->GetNumberOfBlocks();
  size_t nL_CI = streamedData->GetNumberOfBlocks();
  if(nL_PI != nL_CI) {
    dMsg(cout,
         "[ttkTrackingFromOverlap] ERROR: Number of levels differ over time.\n",
         fatalMsg);
    return 0;
  }
  for(size_t l = 0; l < nL_PI; l++) {
    auto timestepsPI = vtkMultiBlockDataSet::SafeDownCast(
      this->previousIterationData->GetBlock(l));
    auto timestepsCI
      = vtkMultiBlockDataSet::SafeDownCast(streamedData->GetBlock(l));
    size_t nT_CI = timestepsCI->GetNumberOfBlocks();

    auto timesteps = vtkSmartPointer<vtkMultiBlockDataSet>::New();

    // First timestep is previous iteration
    timesteps->SetBlock(0, timestepsPI->GetBlock(0));

    // Remaining timesteps are from current iteration
    for(size_t t = 0; t < nT_CI; t++)
      timesteps->SetBlock(t + 1, timestepsCI->GetBlock(t));

    data->SetBlock(l, timesteps);
  }

  return 1;
}

// =============================================================================
// Store Streamed Data
// =============================================================================
int ttkTrackingFromOverlap::storeStreamedData(
  vtkMultiBlockDataSet *streamedData) {

  auto temp = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  size_t nL = streamedData->GetNumberOfBlocks();
  for(size_t l = 0; l < nL; l++) {
    auto timestepsSD
      = vtkMultiBlockDataSet::SafeDownCast(streamedData->GetBlock(l));
    size_t nT = timestepsSD->GetNumberOfBlocks();

    vtkSmartPointer<vtkMultiBlockDataSet> lastTimestep
      = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    lastTimestep->SetBlock(0, timestepsSD->GetBlock(nT - 1));

    temp->SetBlock(l, lastTimestep);
  }

  this->previousIterationData = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  this->previousIterationData->DeepCopy(temp);

  return 1;
}

// =============================================================================
// Compute Nodes
// =============================================================================
int ttkTrackingFromOverlap::computeNodes(vtkMultiBlockDataSet *data) {

  Timer timer;

  size_t nL, nT;
  getNumberOfLevelsAndTimesteps(data, nL, nT);

  // Reusable variables
  vtkPointSet *pointSet = nullptr;
  vtkAbstractArray *labels = nullptr;

  // Compute Nodes
  {
    dMsg(cout,
         "[ttkTrackingFromOverlap] "
         "=======================================================\n",
         infoMsg);
    dMsg(cout, "[ttkTrackingFromOverlap] Computing nodes\n", infoMsg);

    if(this->levelTimeNodesMap.size() != nL)
      this->levelTimeNodesMap.resize(nL);

    for(size_t l = 0; l < nL; l++) {
      {
        stringstream msg;
        msg << "[ttkTrackingFromOverlap] "
               "-------------------------------------------------------"
            << endl
            << "[ttkTrackingFromOverlap] Level Index: " << l << endl;
        dMsg(cout, msg.str(), infoMsg);
      }

      vector<Nodes> &timeNodesMap = this->levelTimeNodesMap[l];
      size_t timeOffset = timeNodesMap.size();
      timeNodesMap.resize(timeOffset + nT);

      for(size_t t = 0; t < nT; t++) {
        getData(data, t, l, this->GetLabelFieldName(), pointSet, labels);

        size_t nPoints = pointSet->GetNumberOfPoints();
        if(nPoints < 1)
          continue;

        switch(this->LabelDataType) {
          vtkTemplateMacro(this->trackingFromOverlap.computeNodes<VTK_TT>(
            (float *)pointSet->GetPoints()->GetVoidPointer(0),
            (VTK_TT *)labels->GetVoidPointer(0), nPoints,
            timeNodesMap[timeOffset + t]));
        }
      }
    }

    {
      stringstream msg;
      msg << "[ttkTrackingFromOverlap] "
             "-------------------------------------------------------"
          << endl
          << "[ttkTrackingFromOverlap] Nodes computed in "
          << timer.getElapsedTime() << " s. (" << threadNumber_
          << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return 1;
}

// =============================================================================
// Compute Tracking Graphs
// =============================================================================
int ttkTrackingFromOverlap::computeTrackingGraphs(vtkMultiBlockDataSet *data) {

  Timer timer;

  size_t nL, nT;
  getNumberOfLevelsAndTimesteps(data, nL, nT);

  if(nT < 2)
    return 1;

  // Reusable variables
  vtkPointSet *pointSet0 = nullptr;
  vtkPointSet *pointSet1 = nullptr;
  vtkAbstractArray *labels0 = nullptr;
  vtkAbstractArray *labels1 = nullptr;

  dMsg(cout,
       "[ttkTrackingFromOverlap] "
       "=======================================================\n",
       infoMsg);
  dMsg(cout, "[ttkTrackingFromOverlap] Computing tracking graphs\n", infoMsg);

  if(this->levelTimeEdgesTMap.size() != nL)
    this->levelTimeEdgesTMap.resize(nL);

  for(size_t l = 0; l < nL; l++) {
    {
      stringstream msg;
      msg << "[ttkTrackingFromOverlap] "
             "-------------------------------------------------------"
          << endl
          << "[ttkTrackingFromOverlap] Level Index: " << l << endl;
      dMsg(cout, msg.str(), infoMsg);
    }

    vector<Edges> &timeEdgesTMap = this->levelTimeEdgesTMap[l];
    size_t timeOffset = timeEdgesTMap.size();
    timeEdgesTMap.resize(timeOffset + nT - 1);

    for(size_t t = 1; t < nT; t++) {
      getData(data, t - 1, l, this->GetLabelFieldName(), pointSet0, labels0);
      getData(data, t, l, this->GetLabelFieldName(), pointSet1, labels1);

      size_t nPoints0 = pointSet0->GetNumberOfPoints();
      size_t nPoints1 = pointSet1->GetNumberOfPoints();
      if(nPoints0 < 1 || nPoints1 < 1)
        continue;

      switch(this->LabelDataType) {
        vtkTemplateMacro(this->trackingFromOverlap.computeOverlap<VTK_TT>(
          (float *)pointSet0->GetPoints()->GetVoidPointer(0),
          (float *)pointSet1->GetPoints()->GetVoidPointer(0),
          (VTK_TT *)labels0->GetVoidPointer(0),
          (VTK_TT *)labels1->GetVoidPointer(0), nPoints0, nPoints1,
          timeEdgesTMap[timeOffset + t - 1]));
      }
    }
  }

  {
    stringstream msg;
    msg << "[ttkTrackingFromOverlap] "
           "-------------------------------------------------------"
        << endl
        << "[ttkTrackingFromOverlap] Tracking graphs computed in "
        << timer.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))."
        << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  return 1;
}

// =============================================================================
// Compute Nesting Trees
// =============================================================================
int ttkTrackingFromOverlap::computeNestingTrees(vtkMultiBlockDataSet *data) {

  Timer timer;

  size_t nL, nT;
  getNumberOfLevelsAndTimesteps(data, nL, nT);

  if(nL < 2)
    return 1;

  // Reusable variables
  vtkPointSet *pointSet0 = nullptr;
  vtkPointSet *pointSet1 = nullptr;
  vtkAbstractArray *labels0 = nullptr;
  vtkAbstractArray *labels1 = nullptr;

  dMsg(cout,
       "[ttkTrackingFromOverlap] "
       "=======================================================\n",
       infoMsg);
  dMsg(cout, "[ttkTrackingFromOverlap] Computing nesting trees\n", infoMsg);

  size_t timeOffset = this->timeLevelEdgesNMap.size();
  this->timeLevelEdgesNMap.resize(timeOffset + nT);

  for(size_t t = 0; t < nT; t++) {
    {
      stringstream msg;
      msg << "[ttkTrackingFromOverlap] "
             "-------------------------------------------------------"
          << endl
          << "[ttkTrackingFromOverlap] Time Index: " << t << endl;
      dMsg(cout, msg.str(), infoMsg);
    }

    vector<Edges> &levelEdgesNMap = this->timeLevelEdgesNMap[timeOffset + t];
    levelEdgesNMap.resize(nL - 1);

    for(size_t l = 1; l < nL; l++) {
      getData(data, t, l - 1, this->GetLabelFieldName(), pointSet0, labels0);
      getData(data, t, l, this->GetLabelFieldName(), pointSet1, labels1);

      size_t nPoints0 = pointSet0->GetNumberOfPoints();
      size_t nPoints1 = pointSet1->GetNumberOfPoints();
      if(nPoints0 < 1 || nPoints1 < 1)
        continue;

      switch(this->LabelDataType) {
        vtkTemplateMacro(this->trackingFromOverlap.computeOverlap<VTK_TT>(
          (float *)pointSet0->GetPoints()->GetVoidPointer(0),
          (float *)pointSet1->GetPoints()->GetVoidPointer(0),
          (VTK_TT *)labels0->GetVoidPointer(0),
          (VTK_TT *)labels1->GetVoidPointer(0), nPoints0, nPoints1,
          levelEdgesNMap[l - 1]));
      }
    }
  }

  {
    stringstream msg;
    msg << "[ttkTrackingFromOverlap] "
           "-------------------------------------------------------"
        << endl
        << "[ttkTrackingFromOverlap] Nesting trees computed in "
        << timer.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))."
        << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 1;
}

// =============================================================================
// Compute Branches
// =============================================================================
int ttkTrackingFromOverlap::computeBranches() {

  size_t nL = this->levelTimeEdgesTMap.size();

  for(size_t l = 0; l < nL; l++)
    this->trackingFromOverlap.computeBranches(
      this->levelTimeEdgesTMap[l], this->levelTimeNodesMap[l]);

  return 1;
}

// =============================================================================
// Request Data
// =============================================================================
int ttkTrackingFromOverlap::RequestData(vtkInformation *request,
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {
  Timer timer;

  // Print Status
  {
    stringstream msg;
    msg << "==================================================================="
           "============="
        << endl
        << "[ttkTrackingFromOverlap] RequestData" << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  // Set Wrapper
  this->trackingFromOverlap.setWrapper(this);

  // Get Input Object
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  auto inputObject = inInfo->Get(vtkDataObject::DATA_OBJECT());

  // Get iteration information
  auto iterationInformation = vtkDoubleArray::SafeDownCast(
    inputObject->GetFieldData()->GetAbstractArray("_ttk_IterationInfo"));

  bool useStreamingOverTime = iterationInformation != nullptr;

  double iteration = 0;
  double nIterations = 0;
  if(useStreamingOverTime) {
    iteration = iterationInformation->GetValue(0);
    nIterations = iterationInformation->GetValue(1);
  }

  // On first iteration reset
  if(!useStreamingOverTime || iteration == 0)
    this->reset();

  // -------------------------------------------------------------------------
  // Prepare Input
  // -------------------------------------------------------------------------
  auto packedInput = vtkSmartPointer<vtkMultiBlockDataSet>::New();

  if(!this->packInputData(inputObject, packedInput))
    return 0;

  // Check input integrity
  if(!this->checkData(packedInput))
    return 0;

  // During streaming append data to previous timestep; Otherwise pass through
  auto data = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  if(useStreamingOverTime && this->previousIterationData != nullptr) {
    if(!this->packStreamedData(packedInput, data))
      return 0;
  } else
    data->ShallowCopy(packedInput);

  // -------------------------------------------------------------------------
  // Process Input
  // -------------------------------------------------------------------------
  // Compute nodes only for current input (previous iterations have already been
  // processed)
  if(!this->computeNodes(packedInput))
    return 0;

  // Compute tracking graphs
  if(!this->computeTrackingGraphs(data))
    return 0;

  // Compute nesting trees
  if(!this->computeNestingTrees(packedInput))
    return 0;

  // Store last timestep for next iteration
  if(useStreamingOverTime && !this->storeStreamedData(packedInput))
    return 0;

  // -------------------------------------------------------------------------
  // Generate Output
  // -------------------------------------------------------------------------
  if(!useStreamingOverTime || iteration == nIterations - 1) {
    // Get Output
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    auto trackingGraph = outInfo->Get(vtkDataObject::DATA_OBJECT());

    // Compute Branches
    if(!this->computeBranches())
      return 0;

    // Mesh Graph
    if(!this->meshNestedTrackingGraph(trackingGraph))
      return 0;
  }

  // -------------------------------------------------------------------------
  // Print total performance
  // -------------------------------------------------------------------------
  if(!useStreamingOverTime) {
    stringstream msg;
    msg << "[ttkTrackingFromOverlap] "
           "======================================================="
        << endl
        << "[ttkTrackingFromOverlap] Nested tracking graph generated in "
        << timer.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))."
        << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 1;
}