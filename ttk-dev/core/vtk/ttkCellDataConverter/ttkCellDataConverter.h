/// \ingroup vtk
/// \class ttkCellDataConverter
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date February 2016
///
/// \brief TTK VTK-filter that converts data types for cell-based scalar fields
/// (for instance, from double to float).
///
/// \param Input Input cell-based scalar field (vtkDataSet)
/// \param Output Output cell-based scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa vtkPointDataConverter

#ifndef _TTK_CELLDATACONVERTER_H
#define _TTK_CELLDATACONVERTER_H

// VTK includes -- to adapt
#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkShortArray.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedShortArray.h>

// VTK module

#include <ttkCellDataConverterModule.h>

// ttk code includes
#include <Wrapper.h>

#include <limits>

class TTKCELLDATACONVERTER_EXPORT ttkCellDataConverter
  : public vtkDataSetAlgorithm,
    protected ttk::Wrapper {

  enum SupportedType {
    Char = 0,
    Double,
    Float,
    Int,
    IdType,
    Short,
    UnsignedShort,
    UnsignedChar,
  };

public:
  static ttkCellDataConverter *New();

  vtkTypeMacro(ttkCellDataConverter, vtkDataSetAlgorithm);

  // default ttk setters
  void SetDebugLevel(int debugLevel) {
    setDebugLevel(debugLevel);
    Modified();
  }

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
  // end of default ttk setters

  vtkSetMacro(ScalarField, std::string);
  vtkGetMacro(ScalarField, std::string);

  void SetOutputType(int outputType) {
    OutputType = outputType;
    Modified();
  }
  vtkGetMacro(OutputType, int);

  void SetUseNormalization(bool onOff) {
    UseNormalization = onOff;
    Modified();
  }
  vtkGetMacro(UseNormalization, int);

protected:
  ttkCellDataConverter();
  ~ttkCellDataConverter() override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  template <typename A, typename B, typename C>
  int convert(vtkDataArray *inputData, vtkDataSet *output);

private:
  bool UseAllCores;
  int ThreadNumber;
  std::string ScalarField;
  int OutputType;
  bool UseNormalization;

  // base code features
  int doIt(vtkDataSet *input, vtkDataSet *output);
  bool needsToAbort() override;
  int updateProgress(const float &progress) override;
};

#endif // _TTK_CELLDATACONVERTER_H
