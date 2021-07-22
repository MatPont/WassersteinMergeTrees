/// \ingroup base
/// \class vtk::ttkFTMTreeUtils
/// \author XXXXX
///

#ifndef _TTKFTMTREEUTILS_H
#define _TTKFTMTREEUTILS_H

#include <ttkUtils.h>
#include <FTMTree.h>
#include <FTMStructures.h>
#include <MergeTreeUtils.h>

#include <vtkCellType.h>

#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkCellData.h>
#include <vtkUnstructuredGrid.h>

#include <PlanarGraphLayout.h>
#include <ttkMacros.h>
#include <vtkCellArray.h>

using namespace ttk;
using namespace ftm;

// ------------------------------------------------------------------------------
// Tree Maker Utils
// ------------------------------------------------------------------------------
void manageInconsistentArcsMultiParent(FTMTree_MT *tree){
  ftm::idNode treeRoot = 0;
  for(int i = 0; i < tree->getNumberOfNodes(); ++i){
    if(tree->getNode(i)->getNumberOfDownSuperArcs() != 0 and tree->getNode(i)->getNumberOfUpSuperArcs() == 0)
      treeRoot = i;
  }
  
  for(int i = 0; i < tree->getNumberOfNodes(); ++i)
    if(tree->getNode(i)->getNumberOfUpSuperArcs() > 1){
      ftm::idNode lowestParent = std::numeric_limits<ftm::idNode>::max();
      for(int j = 0; j < tree->getNode(i)->getNumberOfUpSuperArcs(); ++j){
        auto tParent = tree->getSuperArc(tree->getNode(i)->getUpSuperArcId(j))->getUpNodeId();
        lowestParent = (lowestParent > tParent) ? tParent : lowestParent;
      }
      
      for(int j = 0; j < tree->getNode(i)->getNumberOfUpSuperArcs(); ++j){
        ftm::idSuperArc nodeArcId = tree->getNode(i)->getUpSuperArcId(j);
        auto tParent = tree->getSuperArc(nodeArcId)->getUpNodeId();
        if(tParent != lowestParent){
          
          if(tParent == treeRoot){
            for(int k = 0; k < tree->getNode(i)->getNumberOfDownSuperArcs(); ++k){
              ftm::idSuperArc nodeArcId2 = tree->getNode(i)->getDownSuperArcId(k);
              auto tChildren = tree->getSuperArc(nodeArcId2)->getDownNodeId();
              if(tChildren > i){
                tree->getNode(i)->removeDownSuperArc(nodeArcId2);
                tree->getNode(tChildren)->removeUpSuperArc(nodeArcId2);
                tree->makeSuperArc(tChildren, treeRoot);
                break;
              }
            }
          }
          
          // Delete down arc from old parent
          tree->getNode(tParent)->removeDownSuperArc(nodeArcId);
          // Delete up arc from node
          tree->getNode(i)->removeUpSuperArc(nodeArcId);
        }
      }
    }
}

void removeSelfLink(FTMTree_MT *tree){
  for(int i = 0; i < tree->getNumberOfNodes(); ++i){    
    for(int j = 0; j < tree->getNode(i)->getNumberOfUpSuperArcs(); ++j){
      ftm::idSuperArc nodeArcId = tree->getNode(i)->getUpSuperArcId(j);
      auto tParent = tree->getSuperArc(nodeArcId)->getUpNodeId();
      if(tParent == i){        
        // Delete down arc
        tree->getNode(i)->removeDownSuperArc(nodeArcId);
        // Delete up arc
        tree->getNode(i)->removeUpSuperArc(nodeArcId);
      }
    }
  }
}

// TODO how to detect multiple connected components to avoid some bugs?
MergeTree* makeTree(vtkUnstructuredGrid *treeNodes, vtkUnstructuredGrid *treeArcs){
  // Init Scalars
  //std::cout << "// Init Scalars" << std::endl;
  Scalars *scalars = new Scalars();
  vtkSmartPointer<vtkDataArray> nodesScalar = treeNodes->GetPointData()->GetArray("Scalar"); // 1: Scalar    
  scalars->size = nodesScalar->GetNumberOfTuples();
  scalars->values = ttkUtils::GetVoidPointer(nodesScalar);
  
  // Init Tree
  //std::cout << "// Init Tree" << std::endl;
  Params *params = new Params();
  //FTMTree_MT treeTemp(params, nullptr, scalars, Join_Split);
  FTMTree_MT *tree = new FTMTree_MT(params, nullptr, scalars, Join_Split);
  tree->makeAlloc();
  //tree->makeInit();
  //tree->initComp();
  
  // Add Nodes
  //std::cout << "// Add Nodes" << std::endl;
  vtkSmartPointer<vtkDataArray> nodesId = treeNodes->GetPointData()->GetArray("NodeId"); // 0: NodeId
  //vtkSmartPointer<vtkDataArray> vertexId = treeNodes->GetPointData()->GetArray("VertexId"); // 2: VertexId
  vtkIdType nodesNumTuples = nodesId->GetNumberOfTuples();
  for(vtkIdType i = 0; i < nodesNumTuples; ++i){
    tree->makeNode(i);      
    // perhaps todo: makeNode with vertexId and manage correctly everywhere (is it useful?)
    //auto tmpVertex = vertexId->GetTuple1(i);
    //tree->makeNode(tmpVertex);
  }
  
  // Add Arcs
  //std::cout << "// Add Arcs" << std::endl;
  vtkSmartPointer<vtkDataArray> arcsUp = treeArcs->GetCellData()->GetArray("upNodeId"); // 1: upNodeId
  vtkSmartPointer<vtkDataArray> arcsDown = treeArcs->GetCellData()->GetArray("downNodeId"); // 2: downNodeId
  vtkIdType arcsNumTuples = arcsUp->GetNumberOfTuples();
  std::set<std::tuple<double,double>> added_arcs; // Avoid duplicates
  for(vtkIdType i = 0; i < arcsNumTuples; ++i){
    double downId = arcsDown->GetTuple1(i);
    double upId = arcsUp->GetTuple1(i);
    auto it = added_arcs.find(std::make_tuple(downId, upId));
    if(it == added_arcs.end()){ // arc not added yet
      tree->makeSuperArc(downId, upId); // (down, Up)
      added_arcs.insert(std::make_tuple(downId, upId));
    }
  }
  
  /*treeNodes->GetPointData()->Print(std::cout);
  treeArcs->GetCellData()->Print(std::cout);
  int temp; std::cin >> temp;*/
  
  // Manage inconsistent arcs
  //std::cout << "// Manage inconsistent arcs" << std::endl;
  //tree->printTree2();
  manageInconsistentArcsMultiParent(tree);
  //tree->printTree2();
  
  // Remove self link
  removeSelfLink(tree);
  
  MergeTree *mergeTree = new MergeTree();
  mergeTree->tree = tree;
  mergeTree->scalars = scalars;
  mergeTree->params = params;
  
  return mergeTree;
}

void deleteTree(FTMTree_MT *tree){
  delete tree;
}

#endif
