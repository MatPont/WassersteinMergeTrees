/// \ingroup base
/// \class vtk::ttkFTMTreeUtilsVisu
/// \author XXXXX
///

#ifndef _TTKFTMTREEUTILSVISU_H
#define _TTKFTMTREEUTILSVISU_H

//#include <ttkUtils.h>
#include <FTMTree.h>
#include <FTMStructures.h>
#include <MergeTreeUtils.h>
//#include <ttkFTMTreeUtils.h>

#include <vtkCellType.h>

#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkCellData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellArray.h>
#include <vtkAppendFilter.h>

#include <ttkMacros.h>
//#include <DimensionReduction.h>

using namespace ttk;
using namespace ftm;

// -------------------------------------------------------------------------------------
// Bounds Utils
// -------------------------------------------------------------------------------------

std::vector<double> tupleToVector(std::tuple<double,double,double,double,double,double> &tup);

std::tuple<double,double,double,double,double,double> vectorToTuple(std::vector<double> &vec);

std::tuple<double, double, double, double, double, double> 
getMaximalBounds(std::vector<std::tuple<double, double, double, double, double, double>> &allBounds,
                 std::vector<int> &clusteringAssignment, int clusterID);

std::tuple<double, double, double, double, double, double> 
getRealBounds(vtkUnstructuredGrid *treeNodes, FTMTree_MT *tree, std::vector<int> &nodeCorr);

std::tuple<double, double, double, double, double, double> 
getRealBounds(vtkUnstructuredGrid *treeNodes, FTMTree_MT *tree);

// -------------------------------------------------------------------------------------
// Temporal Subsampling MDS
// -------------------------------------------------------------------------------------
void makeTemporalSubsamplingOutput(std::vector<MergeTree*> &intermediateMTrees,
                                   std::vector<std::vector<double>> &embedding, 
                                   std::vector<MergeTree*> &allMT, std::vector<int> removed,
                                   vtkSmartPointer<vtkUnstructuredGrid> vtkOutputNode1, 
                                   vtkSmartPointer<vtkUnstructuredGrid> vtkOutputArc1,
                                   bool displayRemoved = false);

void makeTemporalSubsamplingMDSOutput(std::vector<MergeTree*> &intermediateMTrees, 
                                      std::vector<std::vector<double>> &distanceMatrix, 
                                      std::vector<MergeTree*> &allMT, std::vector<int> removed,
                                      vtkSmartPointer<vtkUnstructuredGrid> vtkOutputNode1, 
                                      vtkSmartPointer<vtkUnstructuredGrid> vtkOutputArc1, bool metricMDS);

void makeTemporalSubsamplingETDOutput(std::vector<MergeTree*> &intermediateMTrees, 
                                      std::vector<double> &emptyTreeDistances, 
                                      std::vector<MergeTree*> &allMT, std::vector<int> removed,
                                      vtkSmartPointer<vtkUnstructuredGrid> vtkOutputNode1, 
                                      vtkSmartPointer<vtkUnstructuredGrid> vtkOutputArc1, 
                                      double DistanceAxisStretch);

#endif
