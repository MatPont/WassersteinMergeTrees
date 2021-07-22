/// \ingroup base
/// \class ttk::AssignmentAuction
/// \author XXXXX
///
/// Auction algorithm for Balanced and Unbalanced Assignement Problem
///
/// For the unbalanced problem:
///   The cost matrix in input has a size of (n + 1) x (m + 1)
///   - n is the number of jobs, m the number of workers
///   - the nth row contains the cost of not assigning workers
///   - the mth column is the same but with jobs
///   - the last cell (costMatrix[n][m]) is not used


#ifndef _ASSIGNMENTAUCTION_H
#define _ASSIGNMENTAUCTION_H

namespace ttk {

  template <typename dataType>
  class AssignmentAuction : public Debug {

  public:
    AssignmentAuction(){};

    ~AssignmentAuction(){};

    int run(std::vector<matchingTuple> &matchings);
    void runAuctionRound(std::vector<std::vector<dataType>> &costMatrix);
    
    void initFirstRound();
    void initBiddersAndGoods();
    void initEpsilon();
    void epsilonScaling();
    void makeBalancedMatrix(std::vector<std::vector<dataType>> &costMatrix);
    
    bool stoppingCriterion(std::vector<std::vector<dataType>> &costMatrix);
    dataType getRelativePrecision(std::vector<std::vector<dataType>> &costMatrix);    
    dataType getMatchingDistance(std::vector<std::vector<dataType>> &costMatrix);

    inline void clear() {
      rowSize = 0;
      colSize = 0;
    }

    inline void clearMatrix() {
      std::vector<std::vector<dataType>> C
        = *((std::vector<std::vector<dataType>> *)Cptr);
      for(int r = 0, rS0 = rowSize; r < rS0; ++r)
        for(int c = 0, cS0 = colSize; c < cS0; ++c)
          C[r][c] = 0.0;
    }

    inline int setInput(int rowSize_, int colSize_, void *C_) {
      rowSize = rowSize_;
      colSize = colSize_;

      Cptr = C_;
      
      setBalanced((rowSize == colSize));

      return 0;
    }
    
    inline void setBalanced(bool balanced){
      balancedAssignment = balanced;
      if(balancedAssignment)
        goodPrices.resize(colSize, 0);
      else
        goodPrices.resize((colSize-1)+(rowSize-1), 0);
    }
    
    inline void setNumberOfRounds(int noRounds){
      numberOfRounds = noRounds;
    }
    
    inline int getIter(){
      return iter;
    }
    
    inline void setEpsilon(double eps){
      epsilon = eps;
    }
    
    inline void setEpsilonDiviserMultiplier(double div){
      epsilonDiviserMultiplier = div;
    }
    
    inline void setPrices(std::vector<double> prices){
      goodPrices.resize(prices.size(), 0);
      for(int i = 0; i < prices.size(); ++i)
        goodPrices[i] = prices[i];
    }
    
    inline std::vector<double>& getPrices(){
      return goodPrices;
    }
    
    inline std::vector<std::vector<dataType>> getCostMatrix(){
      return *((std::vector<std::vector<dataType>> *)Cptr);
    }
    
    template<class vecType>
    void printVector(std::vector<vecType> &vec){
      for(auto valTemp : vec)
          std::cout << valTemp << " ";
      std::cout << std::endl;
      std::cout << " ------------------------------ " << std::endl;
    }

    void printTableVector(std::vector<std::vector<dataType>> &table){
      for(auto vecTemp : table){
        for(auto valTemp : vecTemp){
          std::cout << valTemp << " ";
        }
        std::cout << std::endl;
      }
      std::cout << " ------------------------------ " << std::endl;
    }
    
  private:
    void *Cptr;

    int rowSize = 0;
    int colSize = 0;
    
    bool balancedAssignment;
    int numberOfRounds = -1;
    int iter = 0;
    double epsilon = -1;
    double epsilonDiviserMultiplier = 0;
    dataType delta_lim = 0.01;
    
    std::vector<int> bidderAssignments{1, -1};
    std::vector<int> goodAssignments{};
    std::vector<double> goodPrices{};
  }; // AssignmentAuction Class

  template <typename type>
  static type abs(const type var) {
    return (var >= 0) ? var : -var;
  }
  
  template<class dataType>
  dataType getMaxValue(std::vector<std::vector<dataType>> &costMatrix, bool balancedAssignment){
    int nRows = costMatrix.size();
    int nCols = costMatrix[0].size();
    dataType maxValue = std::numeric_limits<dataType>::lowest();
    for(int i = 0; i < nRows; ++i){
      for(int j = 0; j < nCols; ++j){
        if(not balancedAssignment and i == nRows-1 and j == nCols-1)
          continue;
        maxValue = (costMatrix[i][j] > maxValue)? costMatrix[i][j] : maxValue;
      }
    }
    return maxValue;
  }
  
  template<class dataType>
  dataType getSecondMinValueVector(std::vector<dataType> &vec){
    dataType secondMin = std::numeric_limits<dataType>::max();
    dataType firstMin = std::numeric_limits<dataType>::max();
    for(auto elem : vec){
      if(elem < firstMin){
        secondMin = firstMin;
        firstMin = elem;
      }else if(elem < secondMin)
        secondMin = elem;
    }
    return secondMin;    
  }
  
  template <typename dataType>
  void AssignmentAuction<dataType>::initEpsilon(){
    if(epsilon == -1){
      std::vector<std::vector<dataType>> costMatrix = getCostMatrix();
      dataType maxValue = getMaxValue<dataType>(costMatrix, balancedAssignment);
      int tRowSize = balancedAssignment ? rowSize : (rowSize-1)+(colSize-1);
      int tColSize = balancedAssignment ? colSize : (rowSize-1)+(colSize-1);
      //epsilon = maxValue * std::min(tRowSize, tColSize)/2;
      epsilon = maxValue/4;
      //epsilon = std::pow(maxValue, 2)/4;
      //epsilon += *std::max_element(goodPrices.begin(), goodPrices.end());
      //epsilon += getSecondMinValueVector(goodPrices);
      if(epsilon == 0)
        epsilon = 1;
      epsilon /= ((epsilonDiviserMultiplier==0) ? 1 : epsilonDiviserMultiplier * 5);
    }
  }
  
  template <typename dataType>
  void AssignmentAuction<dataType>::epsilonScaling(){
    epsilon /= 5; 
    /*if(epsilon < 1e-6)
      epsilon = 1e-6;*/
  }
  
  template <typename dataType>
  void AssignmentAuction<dataType>::initBiddersAndGoods(){
    bidderAssignments.clear();
    goodAssignments.clear();
    if(balancedAssignment){
      bidderAssignments.resize(rowSize, -1);
      goodAssignments.resize(colSize, -1);
    }else{
      bidderAssignments.resize((rowSize-1)+(colSize-1), -1);
      goodAssignments.resize((colSize-1)+(rowSize-1), -1);
    }
    /*goodPrices.clear();    
    goodPrices.resize(colSize, 0);*/
  }
  
  template <typename dataType>
  void AssignmentAuction<dataType>::initFirstRound(){
    iter = 0;
    bidderAssignments[0] = -1;
    //epsilon /= ((epsilonDiviserMultiplier==0) ? 1 : epsilonDiviserMultiplier * 5);
  }
  
  template <typename dataType>
  void AssignmentAuction<dataType>::makeBalancedMatrix(std::vector<std::vector<dataType>> &costMatrix){
    int nRows = costMatrix.size();
    int nCols = costMatrix[0].size();
    costMatrix[nRows-1][nCols-1] = 0;
    
    // Add rows
    for(int i = 0; i < nCols-2; ++i){
      std::vector<dataType> newLine(costMatrix[nRows-1]);
      costMatrix.push_back(newLine);
    }
    // Add columns
    for(int i = 0; i < (nRows-1)+(nCols-1); ++i){
      for(int j = 0; j < nRows-2; ++j){
        costMatrix[i].push_back(costMatrix[i][nCols-1]);
      }
    }
  }
  
  // ----------------------------------------
  // Main Functions
  // ----------------------------------------
  template <typename dataType>
  void AssignmentAuction<dataType>::runAuctionRound(std::vector<std::vector<dataType>> &costMatrix){ 
    std::queue<int> unassignedBidders{};
    for(int i = 0; i < bidderAssignments.size(); ++i)
      unassignedBidders.push(i);
    
    while(! unassignedBidders.empty()){
      int bidderId = unassignedBidders.front();
      unassignedBidders.pop();
      
      // Get good with highest value
      dataType bestValue = std::numeric_limits<dataType>::lowest();
      dataType bestSecondValue = std::numeric_limits<dataType>::lowest();
      int bestGoodId = -1;
      for(int goodId = 0; goodId < goodPrices.size(); ++goodId){
        if(costMatrix[bidderId][goodId] == -1)
          continue;
        dataType goodPrice = goodPrices[goodId];
        dataType value = - costMatrix[bidderId][goodId] - goodPrice;
        if(value > bestValue){
          bestSecondValue = bestValue;
          bestValue = value;
          bestGoodId = goodId;
        }else if(value > bestSecondValue)
          bestSecondValue = value;
      }
      
      // Update assignments
      bidderAssignments[bidderId] = bestGoodId;
      if(goodAssignments[bestGoodId] != -1)
        unassignedBidders.push(goodAssignments[bestGoodId]);
      goodAssignments[bestGoodId] = bidderId;
      
      // Update price 
      double delta = abs<dataType>(bestValue - bestSecondValue) + epsilon;
      goodPrices[bestGoodId] = goodPrices[bestGoodId] + delta;
      
      //printVector(goodPrices);
    }
  }
   
  template <typename dataType>
  int AssignmentAuction<dataType>::run(std::vector<matchingTuple> &matchings){
    std::vector<std::vector<dataType>> costMatrix = getCostMatrix();
    
    initEpsilon();
    
    // Try to avoid price war
    double tempPrice = *std::max_element(goodPrices.begin(), goodPrices.end());
    std::vector<double> savedPrices;
    for(int i = 0; i < goodPrices.size(); ++i){
      auto old = goodPrices[i];
      goodPrices[i] = goodPrices[i] * epsilon / ((tempPrice == 0) ? 1 : tempPrice);
      auto t = old - goodPrices[i];
      savedPrices.push_back(t);
    }
    
    // Make balanced cost matrix
    if(not balancedAssignment)
      makeBalancedMatrix(costMatrix);
    //printTableVector(costMatrix);
    
    // Run acution
    initFirstRound();
    //printVector(goodPrices);
    while(not stoppingCriterion(costMatrix)){
      initBiddersAndGoods();
      runAuctionRound(costMatrix);
      //std::cout << epsilon << std::endl;
      //printVector(goodPrices);
      //printVector(bidderAssignments);
      epsilonScaling();
      iter++;
      if(numberOfRounds != -1 and iter >= numberOfRounds)
        break;
    }
    //printVector(goodPrices);
    
    // Create output matching
    for(int bidderId = 0; bidderId < bidderAssignments.size(); ++bidderId){
      /*int i = std::min(bidderId, rowSize-1);
      int j = std::min(bidderAssignments[bidderId], colSize-1);*/
      int i = bidderId;
      int j = bidderAssignments[bidderId];
      if(balancedAssignment or (not balancedAssignment and not ( i >= rowSize-1 and j >= colSize-1))){
        matchings.push_back(std::make_tuple(i, j, costMatrix[i][j]));
      }
    }
    
    // Set prices as before
    for(int i = 0; i < goodPrices.size(); ++i)
      goodPrices[i] += savedPrices[i];

    return 0;
  }
  
  // ----------------------------------------
  // Stopping Criterion Functions
  // ----------------------------------------
  // Adapted from Persistence Diagrams Auction
  template <typename dataType>
  bool AssignmentAuction<dataType>::stoppingCriterion(std::vector<std::vector<dataType>> &costMatrix) {
    if(bidderAssignments[0] == -1) // Auction not started
      return false;
    dataType delta = 5;
    delta = getRelativePrecision(costMatrix);
    //std::cout << "delta = " << delta << std::endl;
    return not (delta > delta_lim);
  }
  
  // Adapted from Persistence Diagrams Auction
  template <typename dataType>
  dataType AssignmentAuction<dataType>::getRelativePrecision(
                                        std::vector<std::vector<dataType>> &costMatrix) {
    dataType d = this->getMatchingDistance(costMatrix);
    if(d < 1e-12) {
      return 0;
    }
    dataType denominator = d - bidderAssignments.size() * epsilon;
    if(denominator <= 0) {
      return 1;
    } else {
      return d / denominator - 1;
    }
  }
  
  // Adapted from Persistence Diagrams Auction
  template <typename dataType>
  dataType AssignmentAuction<dataType>::getMatchingDistance(std::vector<std::vector<dataType>> &costMatrix) {
    dataType d = 0;
    for(int bidderId = 0; bidderId < bidderAssignments.size(); ++bidderId){
      int i = bidderId;
      int j = bidderAssignments[bidderId];
      d += costMatrix[i][j];
    }
    return d;
  }

} // namespace ttk

#endif
