#ifndef _BUILDWEIGHTS_H
#define _BUILDWEIGHTS_H

namespace ttk {
  template <typename dataType>
  std::vector<KDTree<dataType> *>
    KDTree<dataType>::build(dataType *data,
                            const int &ptNumber,
                            const int &dimension,
                            std::vector<std::vector<dataType>> &weights,
                            const int weight_number) {
    std::vector<KDTree<dataType> *> correspondance_map(ptNumber);
    // First, perform a argsort on the data
    // initialize original index locations
    for(int axis = 0; axis < dimension; axis++) {
      coords_min_.push_back(std::numeric_limits<dataType>::lowest());
      coords_max_.push_back(std::numeric_limits<dataType>::max());
    }
    std::vector<int> idx(ptNumber);
    for(int i = 0; i < ptNumber; i++) {
      idx[i] = i;
    }
    // sort indexes based on comparing values in coordinates
    sort(idx.begin(), idx.end(), [&](int i1, int i2) {
      return data[dimension * i1 + coords_number_]
             < data[dimension * i2 + coords_number_];
    });
    int median_loc = (int)(ptNumber - 1) / 2;
    int median_idx = idx[median_loc];
    correspondance_map[median_idx] = this;

    for(int axis = 0; axis < dimension; axis++) {
      coordinates_.push_back(data[dimension * median_idx + axis]);
    }
    for(int i = 0; i < weight_number; i++) {
      weight_.push_back(weights[i][median_idx]);
      min_subweights_.push_back(weights[i][median_idx]);
    }

    id_ = median_idx;
    parent_ = nullptr;
    level_ = 0;

    if(idx.size() > 2) {
      // Build left leaf
      std::vector<int> idx_left;
      for(int i = 0; i < median_loc; i++) {
        idx_left.push_back(idx[i]);
      }

      KDTree *left = new KDTree(this, (coords_number_ + 1) % dimension, true);
      left->buildRecursive(data, idx_left, ptNumber, dimension, this,
                           correspondance_map, weights, weight_number);
      left_ = left;
    }

    if(idx.size() > 1) {
      // Build right leaf
      std::vector<int> idx_right(ptNumber - median_loc - 1);
      for(int i = 0; i < ptNumber - median_loc - 1; i++) {
        idx_right[i] = idx[i + median_loc + 1];
      }
      KDTree *right = new KDTree(this, (coords_number_ + 1) % dimension, false);
      right->buildRecursive(data, idx_right, ptNumber, dimension, this,
                            correspondance_map, weights, weight_number);
      right_ = right;
    }

    return correspondance_map;
  }

  template <typename dataType>
  void KDTree<dataType>::buildRecursive(
    dataType *data,
    std::vector<int> idx_side,
    const int &ptNumber,
    const int &dimension,
    KDTree<dataType> *parent,
    std::vector<KDTree<dataType> *> &correspondance_map,
    std::vector<std::vector<dataType>> &weights,
    const int weight_number) {
    // First, perform a argsort on the data
    sort(idx_side.begin(), idx_side.end(), [&](int i1, int i2) {
      return data[dimension * i1 + coords_number_]
             < data[dimension * i2 + coords_number_];
    });
    int median_loc = (int)(idx_side.size() - 1) / 2;
    int median_idx = idx_side[median_loc];
    correspondance_map[median_idx] = this;

    for(int axis = 0; axis < dimension; axis++) {
      coordinates_.push_back(data[dimension * median_idx + axis]);
    }

    id_ = median_idx;
    parent_ = parent;
    level_ = parent->level_ + 1;

    for(int i = 0; i < weight_number; i++) {
      weight_.push_back(weights[i][median_idx]);
      min_subweights_.push_back(weights[i][median_idx]);
    }

    if(idx_side.size() > 1) {
      // Once we get to a leaf, update min_subweights of the parents
      for(int w = 0; w < weight_number; w++) {
        this->updateMinSubweight(w);
      }
    }

    // Create bounding box
    for(int axis = 0; axis < dimension; axis++) {
      coords_min_.push_back(parent_->coords_min_[axis]);
      coords_max_.push_back(parent_->coords_max_[axis]);
    }
    if(is_left_ && !this->isRoot()) {
      coords_max_[parent_->coords_number_]
        = parent_->coordinates_[parent_->coords_number_];
    } else if(!is_left_ && !this->isRoot()) {
      coords_min_[parent_->coords_number_]
        = parent_->coordinates_[parent_->coords_number_];
    }

    if(idx_side.size() > 2) {
      // Build left leaf
      std::vector<int> idx_left(median_loc);
      for(int i = 0; i < median_loc; i++) {
        idx_left[i] = idx_side[i];
      }

      KDTree *left = new KDTree(this, (coords_number_ + 1) % dimension, true);
      left->buildRecursive(data, idx_left, ptNumber, dimension, this,
                           correspondance_map, weights, weight_number);
      left_ = left;
    }

    if(idx_side.size() > 1) {
      // Build right leaf
      std::vector<int> idx_right(idx_side.size() - median_loc - 1);
      for(unsigned int i = 0; i < idx_side.size() - median_loc - 1; i++) {
        idx_right[i] = idx_side[i + median_loc + 1];
      }
      KDTree *right = new KDTree(this, (coords_number_ + 1) % dimension, false);
      right->buildRecursive(data, idx_right, ptNumber, dimension, this,
                            correspondance_map, weights, weight_number);
      right_ = right;
    }
    return;
  }
} // namespace ttk

#endif
