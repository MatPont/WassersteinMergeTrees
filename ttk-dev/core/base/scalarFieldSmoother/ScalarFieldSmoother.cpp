#include <ScalarFieldSmoother.h>

using namespace std;
using namespace ttk;

ScalarFieldSmoother::ScalarFieldSmoother() {
  inputData_ = nullptr;
  outputData_ = nullptr;
  dimensionNumber_ = 1;
  mask_ = nullptr;
  setDebugMsgPrefix("ScalarFieldSmoother");
}

ScalarFieldSmoother::~ScalarFieldSmoother() {
}
