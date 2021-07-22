/// \ingroup vtk
/// \class ttkMeshGraph
/// \author Jonas Lukasczyk (jl@jluk.de)
/// \date 01.12.2018
///
/// \brief TTK VTK-filter that generates a mesh for a graph.
///
/// This filter generates for each one dimensional cell (edge) of a
/// 'vtkUnstructuredGrid' a two dimensional cell by mapping a size value to the
/// width of the input cell. The output is a 'vtkUnstructuredGrid' consisting of
/// a set of either quadratic quads or linear polygons.
///
/// VTK wrapping code for the @MeshGraph package.
///
/// \param Input Graph (vtkUnstructuredGrid)
/// \param Output Graph (vtkUnstructuredGrid)
///
/// \sa ttk::MeshGraph

#pragma once

// VTK Module
#include <ttkMeshGraphModule.h>

// VTK includes
#include <ttkAlgorithm.h>

// TTK includes
#include <MeshGraph.h>

class TTKMESHGRAPH_EXPORT ttkMeshGraph : public ttkAlgorithm,
                                         protected ttk::MeshGraph {

private:
  bool UseVariableSize{false};
  int SizeAxis{0};
  float SizeScale{1};
  bool UseQuadraticCells{true};
  int Subdivisions{0};
  bool Tetrahedralize{false};

public:
  vtkSetMacro(UseVariableSize, bool);
  vtkGetMacro(UseVariableSize, bool);

  vtkSetMacro(SizeAxis, int);
  vtkGetMacro(SizeAxis, int);

  vtkSetMacro(SizeScale, float);
  vtkGetMacro(SizeScale, float);

  vtkSetMacro(UseQuadraticCells, bool);
  vtkGetMacro(UseQuadraticCells, bool);

  vtkSetMacro(Subdivisions, int);
  vtkGetMacro(Subdivisions, int);

  vtkSetMacro(Tetrahedralize, bool);
  vtkGetMacro(Tetrahedralize, bool);

  static ttkMeshGraph *New();
  vtkTypeMacro(ttkMeshGraph, ttkAlgorithm);

protected:
  ttkMeshGraph();
  ~ttkMeshGraph();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
