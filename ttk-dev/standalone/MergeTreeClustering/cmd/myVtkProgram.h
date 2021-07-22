#ifndef _MYVTKPROGRAM_H
#define _MYVTKPROGRAM_H

#include <ttkProgramBase.h>

#include <vtkXMLMultiBlockDataReader.h>

using namespace std;
using namespace ttk;

template <class ttkModule>
class myVtkProgram : public vtkProgram<ttkModule> {

public:
    std::vector<vtkCompositeDataSet *> inputsMB_;
    std::vector<vtkSmartPointer<vtkXMLMultiBlockDataReader>> multiBockDataReaders_;

    vtkCompositeDataSet *getInput(const int &inputId) {
        if((inputId < 0) || (inputId >= (int)inputsMB_.size()))
          return NULL;
        return inputsMB_[inputId];
    }

    int load(const vector<string> &inputPaths) {

        int ret = -1;

        for(int i = 0; i < (int)inputPaths.size(); i++) {

            string extension
            = inputPaths[i].substr(inputPaths[i].find_last_of('.') + 1);

            ret = load<vtkXMLMultiBlockDataReader>(inputPaths[i], multiBockDataReaders_);

            if(ret)
                return ret;
        }

        return 0;
    }
    
    template <class vtkReaderClass>
    int load(const std::string &fileName,
        std::vector<vtkSmartPointer<vtkReaderClass>> &readerList) {

        readerList.resize(readerList.size() + 1);
        readerList.back() = vtkSmartPointer<vtkReaderClass>::New();

        readerList.back()->SetFileName(fileName.data());

        // handle debug messages
        /*{
            std::stringstream msg;
            msg << "[ttkProgramBase] Reading input data..." << std::endl;
            dMsg(std::cout, msg.str(), 1);
        }*/

        readerList.back()->Update();
        inputsMB_.push_back(readerList.back()->GetOutput());

        /*{
            std::stringstream msg;
            msg << "[ttkProgramBase]   done! (read "
                << inputs_.back()->GetNumberOfPoints() << " vertices, "
                << inputs_.back()->GetNumberOfCells() << " cells)" << std::endl;
            dMsg(std::cout, msg.str(), Debug::infoMsg);
        }*/

        return 0;
    }
};

#endif
