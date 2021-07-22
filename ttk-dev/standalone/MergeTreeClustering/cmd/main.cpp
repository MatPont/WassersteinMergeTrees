/// \author Your Name Here <Your Email Address Here>.
/// \date The Date Here.
///
/// \brief dummy program example.

// include the local headers
#include <ttkMergeTreeClustering.h>
#include <ttkProgramBase.h>
#include <vtkDataObject.h>
#include <vtkCompositeDataSet.h>

#include "myVtkProgram.h"

using namespace std;
using namespace ttk;

int main(int argc, char **argv) {

    myVtkProgram<ttkMergeTreeClustering> program;
    //vtkProgram<ttkMergeTreeClustering> program;

    // specify local parameters to the TTK module with default values.
    double persistenceThreshold = 0.;
    double epsilon = 5.;
    double epsilon2 = 5.;
    double epsilon3 = 90.;
    double coef = 0.5;
    int noClusters = 1;
    int deterministic = 1;
    int noThreads = 0;

    // register these arguments to the command line parser
    program.parser_.setArgument("P", &persistenceThreshold, "Persistence Threshold", true);
    program.parser_.setArgument("E1", &epsilon, "Epsilon 1", true);
    program.parser_.setArgument("E2", &epsilon2, "Epsilon 2", true);
    program.parser_.setArgument("E3", &epsilon3, "Epsilon 3", true);
    //program.parser_.setArgument("C", &coef, "Coef", true);
    program.parser_.setArgument("K", &noClusters, "Number of Clusters", true);
    program.parser_.setArgument("D", &deterministic, "Deterministic", true);
    program.parser_.setArgument("T", &noThreads, "Number of Threads", true);

    int ret = 0;
    ret = program.init(argc, argv);

    if(ret != 0)
        return ret;

    // change here the arguments of the vtkWrapper that you want to update prior
    // to execution.
    program.ttkObject_->SetPersistenceThreshold(persistenceThreshold);
    program.ttkObject_->SetEpsilonTree1(epsilon);
    program.ttkObject_->SetEpsilonTree2(epsilon);
    program.ttkObject_->SetEpsilon2Tree1(epsilon);
    program.ttkObject_->SetEpsilon2Tree2(epsilon);
    program.ttkObject_->SetEpsilon3Tree1(epsilon);
    program.ttkObject_->SetEpsilon3Tree2(epsilon);
    program.ttkObject_->SetComputeBarycenter(true);
    program.ttkObject_->SetNumberOfBarycenters(noClusters);    
    program.ttkObject_->SetDeterministic(deterministic);
    program.ttkObject_->SetKeepSubtree(false);
    program.ttkObject_->SetPlanarLayout(true);
    program.ttkObject_->SetImportantPairs(25);
    program.ttkObject_->SetNonImportantPairsSpacing(0.1);
    program.ttkObject_->SetNumberOfThreads(noThreads);
    program.ttkObject_->SetInputDataObject(0, program.inputsMB_[0]);

    // execute data processing
    ret = program.run();

    if(ret != 0)
        return ret;

    // save the output
    ret = program.save();

    return ret;
}
