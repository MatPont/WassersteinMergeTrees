#include <ttkPlanarGraphLayout.h>

#include <ttkMacros.h>

#include <vtkAbstractArray.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(ttkPlanarGraphLayout);

ttkPlanarGraphLayout::ttkPlanarGraphLayout() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}
ttkPlanarGraphLayout::~ttkPlanarGraphLayout() {
}

int ttkPlanarGraphLayout::FillInputPortInformation(int port,
                                                   vtkInformation *info) {
  if(port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
  else
    return 0;
  return 1;
}

int ttkPlanarGraphLayout::FillOutputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0)
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
  else
    return 0;
  return 1;
}

int ttkPlanarGraphLayout::RequestData(vtkInformation *request,
                                      vtkInformationVector **inputVector,
                                      vtkInformationVector *outputVector) {
  // Get input and output objects
  auto input = vtkUnstructuredGrid::GetData(inputVector[0]);
  auto output = vtkUnstructuredGrid::GetData(outputVector);

  if(input == nullptr || output == nullptr) {
    return -1;
  }

  // Copy input to output
  output->ShallowCopy(input);

  size_t nPoints = output->GetNumberOfPoints();
  size_t nEdges = output->GetNumberOfCells();

  // Get input arrays
  auto sequenceArray = this->GetInputArrayToProcess(0, inputVector);
  if(this->GetUseSequences() && !sequenceArray) {
    this->printErr("Unable to retrieve sequence array.");
    return 0;
  }

  auto sizeArray = this->GetInputArrayToProcess(1, inputVector);
  if(this->GetUseSizes() && !sizeArray) {
    this->printErr("Unable to retrieve size array.");
    return 0;
  }

  auto branchArray = this->GetInputArrayToProcess(2, inputVector);
  if(this->GetUseBranches() && !branchArray) {
    this->printErr("Unable to retrieve branch array.");
    return 0;
  }

  auto levelArray = this->GetInputArrayToProcess(3, inputVector);
  if(this->GetUseLevels() && !levelArray) {
    this->printErr("Unable to retrieve level array.");
    return 0;
  }

  // Initialize output array
  auto outputArray = vtkSmartPointer<vtkFloatArray>::New();
  outputArray->SetName(this->GetOutputArrayName().data());
  outputArray->SetNumberOfComponents(2); // (x,y) position
  outputArray->SetNumberOfValues(nPoints * 2);

  auto dataType
    = this->GetUseSequences() ? sequenceArray->GetDataType() : VTK_CHAR;
  auto idType = this->GetUseBranches() ? branchArray->GetDataType()
                : this->GetUseLevels() ? levelArray->GetDataType()
                                       : VTK_CHAR;
  //   auto levelType
  //     = this->GetUseLevels()
  //         ? levels->GetDataType()
  //         : this->GetUseBranches() ? branches->GetDataType() : VTK_CHAR;

  //   if(branchType != levelType) {
  //     dMsg(cout,
  //          "[ttkPlanarGraphLayout] ERROR: Branch and Level array must have
  //          the " "same type.\n", fatalMsg);
  //     return 0;
  //   }

  int status = 1;

  switch(vtkTemplate2PackMacro(idType, dataType)) {
    ttkTemplate2IdMacro(
      (status = this->execute<VTK_T1, VTK_T2>(
         // Output
         (float *)outputArray->GetVoidPointer(0),
         // Input
         output->GetCells()->GetData()->GetPointer(0), nPoints, nEdges,
         !this->GetUseSequences() ? nullptr
                                  : (VTK_T2 *)sequenceArray->GetVoidPointer(0),
         !this->GetUseSizes() ? nullptr : (float *)sizeArray->GetVoidPointer(0),
         !this->GetUseBranches() ? nullptr
                                 : (VTK_T1 *)branchArray->GetVoidPointer(0),
         !this->GetUseLevels() ? nullptr
                               : (VTK_T1 *)levelArray->GetVoidPointer(0))));
  }

  //   // Compute layout with base code
  //   switch(vtkTemplate2PackMacro(branchType, sequenceType)) {
  //     vtkTemplate2Macro(
  //       (status = planarGraphLayout.execute<vtkIdType, VTK_T1, VTK_T2>(
  //          // Input
  //          !this->GetUseSequences() ? nullptr
  //                                   : (VTK_T2 *)sequences->GetVoidPointer(0),
  //          !this->GetUseSizes() ? nullptr : (float
  //          *)sizes->GetVoidPointer(0), !this->GetUseBranches() ? nullptr
  //                                  : (VTK_T1 *)branches->GetVoidPointer(0),
  //          !this->GetUseLevels() ? nullptr : (VTK_T1
  //          *)levels->GetVoidPointer(0), output->GetCells()->GetPointer(),
  //          nPoints, nEdges,
  //          // Output
  //          (float *)outputField->GetVoidPointer(0))));
  //   }
  if(status != 1)
    return 0;

  // Add output field to output
  output->GetPointData()->AddArray(outputArray);

  return 1;
}
