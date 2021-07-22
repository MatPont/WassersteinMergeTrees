/// \ingroup base
/// \class ttk::ftr::SuperArc
/// \author Gueunet Charles <charles.gueunet+ttk@gmail.com>
/// \date 2018-01-25
///
/// \brief TTK %FTRGraph graph arc
///
/// This class manage  arcs of the graph structure
/// along with their Segmentation
///
/// \sa ttk::FTRGraph

#ifndef FTR_SUPERARC_H
#define FTR_SUPERARC_H

// local includes
#include "FTRAtomicUF.h"
#include "FTRDataTypes.h"
#include "FTRScalars.h"
#include "FTRSegmentation.h"

// c++ includes
#ifndef TTK_ENABLE_KAMIKAZE
#include <iostream>
#endif

namespace ttk {
  namespace ftr {
    class Node;

    class SuperArc {
    private:
      idNode upNodeId_;
      idNode downNodeId_;
      AtomicUF *ufProp_;
      bool visible_;
      idVertex firstReg_, lastReg_, endV_;
      idSuperArc merged_;
      Segment segmentation_;
#ifndef NDEBUG
      bool fromUp_;
#endif

    public:
      SuperArc(const idNode down = nullNode, const idNode up = nullNode)
        : upNodeId_{up}, downNodeId_{down}, ufProp_{nullptr}, visible_{true},
          firstReg_{nullVertex}, lastReg_{nullVertex}, endV_{nullVertex},
          merged_{nullSuperArc}, segmentation_ {
      }
#ifndef NDEBUG
      , fromUp_ {
        false
      }
#endif
      {}

      idNode getUpNodeId(void) const {
        // Caution. can be nullNode
        return upNodeId_;
      }

      void setUpNodeId(const idNode id) {
        upNodeId_ = id;
      }

      idNode getDownNodeId(void) const {
#ifndef TTK_ENABLE_KAMIKAZE
        if(downNodeId_ == nullNode) {
          std::cerr << "[FTR Graph]: Arc have null down node" << std::endl;
        }
#endif
        return downNodeId_;
      }

      void setDownNodeId(const idNode id) {
        downNodeId_ = id;
      }

      Propagation *getPropagation(void) const {
#ifndef TTK_ENABLE_KAMIKAZE
        if(!ufProp_) {
          std::cerr << "[FTR Graph]: Arc have null UF propagation" << std::endl;
        }
#endif
        return ufProp_->find()->getPropagation();
      }

      void setUfProp(AtomicUF *const UFprop) {
        ufProp_ = UFprop;
      }

      bool hide(void) {
        bool old = visible_;
        visible_ = false;
        return old;
      }

      bool isVisible(void) const {
        return visible_;
      }

      // for regular vertices
      void visit(const idVertex v) {
        // firstV only set once
        if(firstReg_ == nullVertex)
          firstReg_ = v;

        lastReg_ = v;
      }

      // for all vertices:
      // used to close arc at the end if not close during construct
      void setEnd(const idVertex v) {
        // avoid an are being reclosed by a higher join
        if(endV_ == nullVertex) {
          endV_ = v;
        }
      }

      idVertex getEnd(void) const {
        return endV_;
      }

      bool isEmpty() const {
        return firstReg_ == nullVertex;
      }

      idVertex getFirstReg() const {
        return firstReg_;
      }

      idVertex getLastReg() const {
        return lastReg_;
      }

      void merge(const idSuperArc arc) {
        if(merged_ == nullSuperArc) {
          merged_ = arc;
          hide();
        }
      }

      bool merged(void) const {
        return merged_ != nullSuperArc;
      }

      idSuperArc mergedIn(void) const {
        return merged_;
      }

      void restore(void) {
        visible_ = true;
        merged_ = nullSuperArc;
      }

      const decltype(segmentation_) &segmentation() const {
        return segmentation_;
      }

      decltype(segmentation_) &segmentation() {
        return segmentation_;
      }

#ifndef NDEBUG
      void setFromUp(bool up) {
        fromUp_ = up;
      }

      bool getFromUp(void) const {
        return fromUp_;
      }
#endif
    };
  } // namespace ftr
} // namespace ttk

#endif /* end of include guard: FTR_SUPERARC_H */
