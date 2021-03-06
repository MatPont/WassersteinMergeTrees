/// \ingroup vtkWrappers
/// \class ttkTopologicalCompressionWriter
/// \author Maxime Soler <soler.maxime@total.com>
/// \date 21/04/2017
///
/// \brief VTK-filter that wraps the topologicalCompressionWriter processing
/// package.

#ifndef _VTK_TOPOLOGICALCOMPRESSIONWRITER_H
#define _VTK_TOPOLOGICALCOMPRESSIONWRITER_H

// TTK
#include <TopologicalCompression.h>
#include <ttkTriangulationAlgorithm.h>

// VTK
#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkWriter.h>

// VTK Module
#include <ttkTopologicalCompressionWriterModule.h>

// STD
#include <fstream>
#include <iostream>
#include <limits.h>
#include <string>

class TTKTOPOLOGICALCOMPRESSIONWRITER_EXPORT ttkTopologicalCompressionWriter
  : public vtkWriter {

public:
  static ttkTopologicalCompressionWriter *New();

  void SetThreads() {
    if(!UseAllCores)
      threadNumber_ = ThreadNumber;
    else {
      threadNumber_ = ttk::OsCall::getNumberOfCores();
    }
    Modified();
  }

  void SetThreadNumber(int threadNumber) {
    ThreadNumber = threadNumber;
    SetThreads();
  }

  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }

  vtkTypeMacro(ttkTopologicalCompressionWriter, vtkWriter);

  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  vtkSetMacro(ScalarField, std::string);
  vtkGetMacro(ScalarField, std::string);

  vtkSetMacro(ScalarFieldId, int);
  vtkGetMacro(ScalarFieldId, int);

  vtkGetMacro(Tolerance, double);
  vtkSetMacro(Tolerance, double);

  vtkGetMacro(MaximumError, double);
  vtkSetMacro(MaximumError, double);

  vtkGetMacro(ZFPBitBudget, double);
  vtkSetMacro(ZFPBitBudget, double);

  vtkGetMacro(ZFPOnly, bool);
  vtkSetMacro(ZFPOnly, bool);

  vtkGetMacro(CompressionType, int);
  vtkSetMacro(CompressionType, int);

  vtkGetMacro(NbSegments, int);
  vtkSetMacro(NbSegments, int);

  vtkGetMacro(NbVertices, int);
  vtkSetMacro(NbVertices, int);

  vtkGetMacro(SQMethod, std::string);
  vtkSetMacro(SQMethod, std::string);

  vtkSetMacro(Subdivide, bool);
  vtkGetMacro(Subdivide, bool);

  vtkSetMacro(UseTopologicalSimplification, bool);
  vtkGetMacro(UseTopologicalSimplification, bool);

  void SetDebugLevel(int debugLevel) {
    d.setDebugLevel(debugLevel);
  }

  inline void SetSQMethodPV(int c) {
    switch(c) {
      case 1:
        SetSQMethod("r");
        break;
      case 2:
        SetSQMethod("d");
        break;
      case 0:
      default:
        SetSQMethod("");
        break;
    }
  }

protected:
  // Regular writer management.
  ttkTopologicalCompressionWriter();
  ~ttkTopologicalCompressionWriter();
  virtual int FillInputPortInformation(int port, vtkInformation *info) override;
  void WriteData() override;

  // TTK management.
  vtkDataArray *GetInputScalarField(vtkImageData *vti);
  void ComputeTriangulation(vtkImageData *vti);
  int AllocateOutput(vtkDataArray *vda);
  void PerformCompression(vtkDataArray *vda);

protected:
  mutable int threadNumber_;

private:
  // Writer parameters.
  char *FileName;
  double ZFPBitBudget;
  bool ZFPOnly;
  int CompressionType;

  // Compression results.
  std::string ScalarField;
  int ScalarFieldId;
  int *Segmentation;
  double Tolerance;
  double MaximumError;
  int NbSegments;
  int NbVertices;
  std::string SQMethod;
  bool Subdivide;
  bool UseTopologicalSimplification;

  // TTK objects.
  vtkSmartPointer<vtkDataArray> outputScalarField;
  ttkTriangulation triangulation;
  ttk::TopologicalCompression topologicalCompression;
  int ThreadNumber;
  bool UseAllCores;

  ttk::Debug d;

  // Whatever.
  ttkTopologicalCompressionWriter(const ttkTopologicalCompressionWriter &);
  void operator=(const ttkTopologicalCompressionWriter &);

  // give ttkCinemaWriter access to ttkTopologicalCompressionWriter
  // protected member functions
  friend class ttkCinemaWriter;
};

#endif // _VTK_TOPOLOGICALCOMPRESSIONWRITER_H
