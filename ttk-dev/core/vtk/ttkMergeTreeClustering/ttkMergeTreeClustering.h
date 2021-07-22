/// Provide your information
///
/// \ingroup vtk
/// \class ttkMergeTreeClustering
/// \author XXXXX
/// \date 2020.
///
/// \brief TTK VTK-filter that wraps the ttk::MergeTreeClustering module.
///
/// This VTK filter uses the ttk::MergeTreeClustering module to compute the edit distance
/// between two merge trees.
///
/// \param Input vtkUnstructuredGrid First Tree Skeleton Nodes
/// \param Input vtkUnstructuredGrid First Tree Skeleton Arcs
/// \param Input vtkUnstructuredGrid Second Tree Skeleton Nodes
/// \param Input vtkUnstructuredGrid Second Tree Skeleton Arcs
/// \param Output vtkUnstructuredGrid Matching
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::MergeTreeClustering
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkMergeTreeClusteringModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <ttkTriangulation.h>
//#include <ttkUtils.h>

// TTK Base Includes
#include <MergeTreeDistance.h>
#include <MergeTreeBarycenter.h>
#include <MergeTreeClustering.h>
#include <MergeTreeTemporalSubsampling.h>

class TTKMERGETREECLUSTERING_EXPORT ttkMergeTreeClustering
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  , protected ttk::MergeTreeDistance // and we inherit from the base class
{
private:
  //std::vector<matchingType> matchings_{};
  
  /**
   * Add all filter parameters only as private member variables and
   * initialize them here.
   */
  // Input Options
  bool Epsilon1UseFarthestSaddle = false;
  double EpsilonTree1 = 0.;
  double EpsilonTree2 = 0.;
  double Epsilon2Tree1 = 0.;
  double Epsilon2Tree2 = 0.;
  double Epsilon3Tree1 = 100.; 
  double Epsilon3Tree2 = 100.;
  double PersistenceThreshold = 0.;
  
  // Output Options
  double DistanceAxisStretch = 1.;
  bool OutputTrees = true;
  bool SepareOutputTrees = false;
  bool OutputSegmentation = false;
  bool PlanarLayout = false;
  bool BranchDecompositionPlanarLayout = false;
  double BranchSpacing = 1.;
  bool RescaleTreesIndividually = false;
  double DimensionSpacing = 1.;
  int DimensionToShift = 0;
  double ImportantPairs = 50.;
  double ImportantPairsSpacing = 1.;
  double NonImportantPairsSpacing = 1.;
  double NonImportantPairsProximity = 0.05;
  
  // Execution Options
  double Alpha = 0.5;
  int AssignmentSolver = 0;
  bool ProgressiveComputation = false;
  double ProgressiveSpeedDivisor = 4.0;
  bool BranchDecomposition = true;
  bool NormalizedWasserstein = true;
  double NormalizedWassersteinReg = 0.;
  bool RescaledWasserstein = false;
  bool KeepSubtree = false;
  bool UseMinMaxPair = true;
  double JoinSplitMixtureCoefficient = 0.5;
  bool DeleteMultiPersPairs = false;
  
  bool ComputeBarycenter = false;
  bool ProgressiveBarycenter = false;
  int NumberOfBarycenters = 1;
  double Tol = 0.0;
  
  bool Parallelize = true;
  int NumberOfThreads = 0;
  bool Deterministic = false;
  
  bool TemporalSubsampling = false;
  double RemovalPercentage = 50.;
  bool UseL2Distance = false;
  bool TemporalSubsamplingMDS = true;
  bool MetricMDS = true;
  
  bool Verbose = true;
  
public:
  /**
   * Automatically generate getters and setters of filter
   * parameters via vtkMacros.
   */
  vtkSetMacro(Epsilon1UseFarthestSaddle, bool);
  vtkGetMacro(Epsilon1UseFarthestSaddle, bool);
  
  vtkSetMacro(EpsilonTree1, double);
  vtkGetMacro(EpsilonTree1, double);
  
  vtkSetMacro(EpsilonTree2, double);
  vtkGetMacro(EpsilonTree2, double);
  
  vtkSetMacro(Epsilon2Tree1, double);
  vtkGetMacro(Epsilon2Tree1, double);
  
  vtkSetMacro(Epsilon2Tree2, double);
  vtkGetMacro(Epsilon2Tree2, double);
  
  vtkSetMacro(Epsilon3Tree1, double);
  vtkGetMacro(Epsilon3Tree1, double);
  
  vtkSetMacro(Epsilon3Tree2, double);
  vtkGetMacro(Epsilon3Tree2, double);
  
  vtkSetMacro(PersistenceThreshold, double);
  vtkGetMacro(PersistenceThreshold, double);
  
  vtkSetMacro(Alpha, double);
  vtkGetMacro(Alpha, double);
  
  vtkSetMacro(AssignmentSolver, int);
  vtkGetMacro(AssignmentSolver, int);
  
  vtkSetMacro(DistanceAxisStretch, double);
  vtkGetMacro(DistanceAxisStretch, double);
  
  vtkSetMacro(OutputTrees, bool);
  vtkGetMacro(OutputTrees, bool);
  
  vtkSetMacro(SepareOutputTrees, bool);
  vtkGetMacro(SepareOutputTrees, bool);
  
  vtkSetMacro(OutputSegmentation, bool);
  vtkGetMacro(OutputSegmentation, bool);
  
  vtkSetMacro(PlanarLayout, bool);
  vtkGetMacro(PlanarLayout, bool);
  
  vtkSetMacro(BranchDecompositionPlanarLayout, bool);
  vtkGetMacro(BranchDecompositionPlanarLayout, bool);
  
  vtkSetMacro(BranchSpacing, double);
  vtkGetMacro(BranchSpacing, double);  
  
  vtkSetMacro(RescaleTreesIndividually, bool);
  vtkGetMacro(RescaleTreesIndividually, bool);
  
  vtkSetMacro(DimensionSpacing, double);
  vtkGetMacro(DimensionSpacing, double);
  
  vtkSetMacro(DimensionToShift, int);
  vtkGetMacro(DimensionToShift, int);
  
  vtkSetMacro(ImportantPairs, double);
  vtkGetMacro(ImportantPairs, double);
  
  vtkSetMacro(ImportantPairsSpacing, double);
  vtkGetMacro(ImportantPairsSpacing, double);
  
  vtkSetMacro(NonImportantPairsSpacing, double);
  vtkGetMacro(NonImportantPairsSpacing, double);
  
  vtkSetMacro(NonImportantPairsProximity, double);
  vtkGetMacro(NonImportantPairsProximity, double);
  
  vtkSetMacro(ProgressiveComputation, bool);
  vtkGetMacro(ProgressiveComputation, bool);
  
  vtkSetMacro(ProgressiveSpeedDivisor, double);
  vtkGetMacro(ProgressiveSpeedDivisor, double);
  
  vtkSetMacro(BranchDecomposition, bool);
  vtkGetMacro(BranchDecomposition, bool);
  
  vtkSetMacro(Parallelize, bool);
  vtkGetMacro(Parallelize, bool);
  
  vtkSetMacro(NumberOfThreads, int);
  vtkGetMacro(NumberOfThreads, int);
  
  vtkSetMacro(Deterministic, bool);
  vtkGetMacro(Deterministic, bool);
  
  vtkSetMacro(NormalizedWasserstein, bool);
  vtkGetMacro(NormalizedWasserstein, bool);
  
  vtkSetMacro(NormalizedWassersteinReg, double);
  vtkGetMacro(NormalizedWassersteinReg, double);
  
  vtkSetMacro(RescaledWasserstein, bool);
  vtkGetMacro(RescaledWasserstein, bool);
  
  vtkSetMacro(KeepSubtree, bool);
  vtkGetMacro(KeepSubtree, bool);

  vtkSetMacro(UseMinMaxPair, bool);
  vtkGetMacro(UseMinMaxPair, bool);
  
  vtkSetMacro(DeleteMultiPersPairs, bool);
  vtkGetMacro(DeleteMultiPersPairs, bool);
  
  vtkSetMacro(JoinSplitMixtureCoefficient, double);
  vtkGetMacro(JoinSplitMixtureCoefficient, double);
  
  vtkSetMacro(ComputeBarycenter, bool);
  vtkGetMacro(ComputeBarycenter, bool);
  
  vtkSetMacro(ProgressiveBarycenter, bool);
  vtkGetMacro(ProgressiveBarycenter, bool);
  
  vtkSetMacro(NumberOfBarycenters, int);
  vtkGetMacro(NumberOfBarycenters, int);
  
  vtkSetMacro(Tol, double);
  vtkGetMacro(Tol, double);
  
  vtkSetMacro(TemporalSubsampling, bool);
  vtkGetMacro(TemporalSubsampling, bool);
  
  vtkSetMacro(RemovalPercentage, double);
  vtkGetMacro(RemovalPercentage, double);
  
  vtkSetMacro(UseL2Distance, double);
  vtkGetMacro(UseL2Distance, double);
  
  vtkSetMacro(TemporalSubsamplingMDS, bool);
  vtkGetMacro(TemporalSubsamplingMDS, bool);

  vtkSetMacro(MetricMDS, bool);
  vtkGetMacro(MetricMDS, bool);
  
  vtkSetMacro(Verbose, bool);
  vtkGetMacro(Verbose, bool);
  
  /**
   * This static method and the macro below are VTK conventions on how to
   * instantiate VTK objects. You don't have to modify this.
   */
  static ttkMergeTreeClustering *New();
  vtkTypeMacro(ttkMergeTreeClustering, ttkAlgorithm);

protected:
  /**
   * Implement the filter constructor and destructor
   * (see cpp file)
   */
  ttkMergeTreeClustering();
  ~ttkMergeTreeClustering() override;

  /**
   * Specify the input data type of each input port
   * (see cpp file)
   */
  int FillInputPortInformation(int port, vtkInformation *info) override;

  /**
   * Specify the data object type of each output port
   * (see cpp file)
   */
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  /**
   * Pass VTK data to the base code and convert base code output to VTK
   * (see cpp file)
   */
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
