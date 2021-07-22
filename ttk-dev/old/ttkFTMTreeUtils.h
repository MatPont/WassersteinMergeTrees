/// \ingroup base
/// \class vtk::ttkFTMTreeUtils
/// \author Mathieu Pont <mathieu.pont@outlook.com>
///

#ifndef _TTKFTMTREEUTILS_H
#define _TTKFTMTREEUTILS_H

#include <ttkUtils.h>
#include <FTMTree.h>
#include <FTMStructures.h>
#include <FTMTreeUtils.h>

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

template<class dataType>
std::vector<float> treePlanarLayoutImpl(FTMTree_MT *tree, int dataTypeInt,
                                      std::tuple<double, double, double, double, double, double> oldBounds){
  auto nPoints = getRealNumberOfNodes(tree);
  auto nEdges = tree->getNumberOfSuperArcs();
  auto outputArray = vtkSmartPointer<vtkFloatArray>::New();
  outputArray->SetName("Layout");
  outputArray->SetNumberOfComponents(2); // (x,y) position
  outputArray->SetNumberOfValues(nPoints * 2);
  
  vtkSmartPointer<vtkUnstructuredGrid> connectivity = vtkSmartPointer<vtkUnstructuredGrid>::New();
  
  dataType scalars[tree->getNumberOfNodes()];
  int levels[tree->getNumberOfNodes()];
  int branches[tree->getNumberOfNodes()];
  float sizes[tree->getNumberOfNodes()];
  
  int cptNode = 0;
  int noLeaves = getNumberOfLeaves(tree);
  std::vector<int> treeSimplexId(tree->getNumberOfNodes());
  std::queue<std::tuple<idNode, int, float>> queue;
  queue.push(std::make_tuple(getRoot(tree), 0, noLeaves));
  while(!queue.empty()){
    auto tup = queue.front();
    queue.pop();
    idNode node = std::get<0>(tup);
    int level = std::get<1>(tup);
    float pos = std::get<2>(tup);
    
    // Get and insert point
    treeSimplexId[node] = cptNode;
    scalars[cptNode] = -level; //tree->getValue<dataType>(node); //* (1/(level+1)*100);
    levels[cptNode] = level;
    branches[cptNode] = level;
    sizes[cptNode] = pos;
    ++cptNode;
    
    // Add cell connecting parent
    if(!isRoot(tree, node)){
      vtkIdType pointIds[2];
      pointIds[0] = treeSimplexId[getParent(tree, node)];
      pointIds[1] = treeSimplexId[node];
      connectivity->InsertNextCell(VTK_LINE, 2, pointIds);
    }
    
    // Push children to the queue
    auto children = getChildren(tree, node);
    for(int i = 0; i < children.size(); ++i){
      auto child = children[i];
      float newPos = pos + i; ///(2*(level+1)); //+1)/2 - children.size()/2;
      queue.push(std::make_tuple(child, level+1, newPos));
    }
  }
  
  // Call PlanarGraphLayout
  auto idType = VTK_INT;
  PlanarGraphLayout layoutCompute;
  switch(vtkTemplate2PackMacro(idType, dataTypeInt)) {
    ttkTemplate2IdMacro(
      (/*status =*/ layoutCompute.execute<VTK_T1, VTK_T2>(
         // Output
         (float *)outputArray->GetVoidPointer(0),
         // Input
         (VTK_T1 *)connectivity->GetCells()->GetPointer(), nPoints, nEdges,
         (VTK_T2 *)scalars, (float *)sizes, nullptr, (VTK_T1 *)levels // (VTK_T1 *)branches
      )));
  }
  
  auto ret = (float *)outputArray->GetVoidPointer(0);
  
  // Create output vector and swap x and y coordinates
  std::vector<float> retVec(outputArray->GetNumberOfValues());
  for(int i = 0; i < outputArray->GetNumberOfValues(); i+=2){
    //std::cout << ret[i] << std::endl;
    retVec[i] = ret[i+1];
    retVec[i+1] = ret[i];
  }
  
  // Rescale coordinates according original bounds
  float x_min, y_min, x_max, y_max;
  x_min = std::numeric_limits<float>::max();
  y_min = std::numeric_limits<float>::max();
  x_max = std::numeric_limits<float>::min();
  y_max = std::numeric_limits<float>::min();
  for(int i = 0; i < outputArray->GetNumberOfValues(); i+=2){
    x_min = std::min(x_min, retVec[i]);
    x_max = std::max(x_max, retVec[i]);
    y_min = std::min(y_min, retVec[i+1]);
    y_max = std::max(y_max, retVec[i+1]);
  }
  auto newBounds = std::make_tuple(x_min, x_max, y_min, y_max, 0, 0);
  
  /*std::cout << std::get<0>(newBounds) << " _ " << std::get<1>(newBounds) << " _ " << std::get<2>(newBounds) << " _ " << std::get<3>(newBounds) << " _ " << std::get<4>(newBounds) << " _ " << std::get<5>(newBounds) << std::endl;
  std::cout << std::get<0>(oldBounds) << " _ " << std::get<1>(oldBounds) << " _ " << std::get<2>(oldBounds) << " _ " << std::get<3>(oldBounds) << " _ " << std::get<4>(oldBounds) << " _ " << std::get<5>(oldBounds) << std::endl;*/
  
  for(int i = 0; i < outputArray->GetNumberOfValues(); i+=2){
    // x coordinate
    auto divisor1 = std::get<1>(newBounds)-std::get<0>(newBounds);
    divisor1 = (divisor1 == 0 ? 1 : divisor1);
    retVec[i] = (retVec[i]-std::get<0>(newBounds)) / divisor1;
    retVec[i] = retVec[i] * (std::get<1>(oldBounds) - std::get<0>(oldBounds)) + std::get<0>(oldBounds);
    
    // y coordinate
    auto divisor2 = std::get<3>(newBounds)-std::get<2>(newBounds);
    divisor2 = (divisor2 == 0 ? 1 : divisor2);
    retVec[i+1] = (retVec[i+1]-std::get<2>(newBounds)) / divisor2;
    retVec[i+1] = retVec[i+1] * (std::get<3>(oldBounds) - std::get<2>(oldBounds)) + std::get<2>(oldBounds);
  }
  
  /*tree->printTree2();
  outputArray->PrintSelf(std::cout,vtkIndent(2));*/
  
  return retVec;
}

std::vector<float> treePlanarLayout(FTMTree_MT *tree, int dataType,
                                    std::tuple<double, double, double, double, double, double> oldBounds){  
  std::vector<float> res;
  switch(dataType){
    vtkTemplateMacro(res = treePlanarLayoutImpl<VTK_TT>(tree, dataType, oldBounds));
  }
  return res;
}

void getTreeBranching(FTMTree_MT *tree, std::vector<idNode> &branching, std::vector<int> &branchingID){
  branching = std::vector<idNode>(tree->getNumberOfNodes());
  branchingID = std::vector<int>(tree->getNumberOfNodes(), -1);
  int branchID = 0;
  std::queue<idNode> queue;
  queue.push(getRoot(tree));
  while(!queue.empty()){
    idNode node = queue.front();
    queue.pop();
    if(isLeaf(tree, node))
      continue;
    auto nodeOrigin = tree->getNode(node)->getOrigin();
    idNode parentNodeOrigin = nodeOrigin;
    while(parentNodeOrigin != node){
      branching[parentNodeOrigin] = node;
      branchingID[parentNodeOrigin] = branchID;
      parentNodeOrigin = getParent(tree, parentNodeOrigin);
    }
    ++branchID;
    for(auto child : getChildren(tree, node))
      queue.push(child);
  }
}

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

MergeTree* makeTree(vtkUnstructuredGrid *treeNodes, vtkUnstructuredGrid *treeArcs){
    // Init Scalars
    Scalars *scalars = new Scalars();
    vtkSmartPointer<vtkDataArray> nodesScalar = treeNodes->GetPointData()->GetArray("Scalar"); // 1: Scalar    
    scalars->size = nodesScalar->GetNumberOfTuples();
    scalars->values = ttkUtils::GetVoidPointer(nodesScalar);
    
    // Init Tree
    Params *params = new Params();
    //FTMTree_MT treeTemp(params, nullptr, scalars, Join_Split);
    FTMTree_MT *tree = new FTMTree_MT(params, nullptr, scalars, Join_Split);
    tree->makeAlloc();
    //tree->makeInit();
    //tree->initComp();
    
    //FTMTree_MT tree_temp(params, nullptr, scalars, Join_Split);
    
    // Add Nodes
    vtkSmartPointer<vtkDataArray> nodesId = treeNodes->GetPointData()->GetArray("NodeId"); // 0: NodeId
    //vtkSmartPointer<vtkDataArray> vertexId = treeNodes->GetPointData()->GetArray("VertexId"); // 2: VertexId
    vtkIdType nodesNumTuples = nodesId->GetNumberOfTuples();    
    for (vtkIdType i = 0; i < nodesNumTuples; ++i){
        tree->makeNode(i);      
        // perhaps todo: makeNode with vertexId and manage correctly everywhere (is it useful?)
        //auto tmpVertex = vertexId->GetTuple1(i);
        //tree->makeNode(tmpVertex);
    }
    
    // Add Arcs
    vtkSmartPointer<vtkDataArray> arcsUp = treeArcs->GetCellData()->GetArray("upNodeId"); // 1: upNodeId
    vtkSmartPointer<vtkDataArray> arcsDown = treeArcs->GetCellData()->GetArray("downNodeId"); // 2: downNodeId
    vtkIdType arcsNumTuples = arcsUp->GetNumberOfTuples();
    std::set<std::tuple<double,double>> added_arcs; // Avoid duplicates
    for (vtkIdType i = 0; i < arcsNumTuples; ++i){
        double downId = arcsDown->GetTuple1(i);
        double upId = arcsUp->GetTuple1(i);
        auto it = added_arcs.find(std::make_tuple(downId, upId));
        if(it == added_arcs.end()){ // arc not added yet
            tree->makeSuperArc(downId, upId); // (down, Up)
            added_arcs.insert(std::make_tuple(downId, upId));
        }
    }
    
    // Manage inconsistent arcs
    //tree->printTree2();
    manageInconsistentArcsMultiParent(tree);
    //tree->printTree2();
    
    MergeTree *mergeTree = new MergeTree();
    mergeTree->tree = tree;
    mergeTree->scalars = scalars;
    mergeTree->params = params;
    
    return mergeTree;
}

void deleteTree(FTMTree_MT *tree){
    delete tree;
}

std::tuple<double, double, double, double, double, double> 
getMaximalBounds(std::vector<std::tuple<double, double, double, double, double, double>> &allBounds,
                 std::vector<int> &clusteringAssignment, int clusterID){
  double x_min = std::numeric_limits<double>::max();
  double y_min = std::numeric_limits<double>::max();
  double z_min = std::numeric_limits<double>::max();
  double x_max = std::numeric_limits<double>::min();
  double y_max = std::numeric_limits<double>::min();
  double z_max = std::numeric_limits<double>::min();
  for(int i = 0; i < allBounds.size(); ++i)
    if(clusteringAssignment[i] == clusterID){
      x_min = std::min(x_min, std::get<0>(allBounds[i]));
      x_max = std::max(x_max, std::get<1>(allBounds[i]));
      y_min = std::min(y_min, std::get<2>(allBounds[i]));
      y_max = std::max(y_max, std::get<3>(allBounds[i]));
      z_min = std::min(z_min, std::get<4>(allBounds[i]));
      z_max = std::max(z_max, std::get<5>(allBounds[i]));
    }
  return std::make_tuple(x_min, x_max, y_min, y_max, z_min, z_max);
}

std::tuple<double, double, double, double, double, double> 
getRealBounds(vtkUnstructuredGrid *treeNodes, FTMTree_MT *tree){
  double x_min = std::numeric_limits<double>::max();
  double y_min = std::numeric_limits<double>::max();
  double z_min = std::numeric_limits<double>::max();
  double x_max = std::numeric_limits<double>::min();
  double y_max = std::numeric_limits<double>::min();
  double z_max = std::numeric_limits<double>::min();
  std::queue<idNode> queue;
  queue.push(getRoot(tree));
  while(!queue.empty()){
    idNode node = queue.front();
    queue.pop();
    double* point = treeNodes->GetPoints()->GetPoint(node);
    x_min = std::min(x_min, point[0]);
    x_max = std::max(x_max, point[0]);
    y_min = std::min(y_min, point[1]);
    y_max = std::max(y_max, point[1]);
    z_min = std::min(z_min, point[2]);
    z_max = std::max(z_max, point[2]);
    for(auto child : getChildren(tree, node))
      queue.push(child);
  }
  return std::make_tuple(x_min, x_max, y_min, y_max, z_min, z_max);
}

#endif
