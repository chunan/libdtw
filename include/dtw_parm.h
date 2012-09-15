#ifndef __PARM_H__
#define __PARM_H__

#include <deque>
#include <cassert>
#include <vector>
#include "segtree.h"
#include "feature.h"
#include "utility.h"


namespace DtwUtil {

  class DtwParm {

    public:
      enum ParmType {FRAME = 0, SEGMENT = 1, SUPERSEG = 2};
      DtwParm() { Init(); }
      ~DtwParm() { ClearSeg(); }

      /* Load data */
      void LoadParm(const string feat_fname); // Load feature only
      void LoadParm(const string feat_fname,  // Load feature and segment
                    const double bseg_ratio,
                    const string tree_fname = "");
      void LoadParm(const string feat_fname,  // Load feature,segment,supersegment
                    const double bseg_ratio,
                    const double superseg_ratio,
                    const unsigned width,
                    const unsigned gran,
                    const string tree_fname = "");


      /* accessors */
      ParmType Type() const { return m_type; }
      unsigned NumFrame() const { return m_feat.LT(); }
      const DenseFeature& Feat() const { return m_feat; }
      // tree related
      const SegTree& Segtree() const { return m_segtree; }
      // basic segment related
      unsigned NumBasicSeg() const { return m_index.size(); }
      int BasicSegTreeIndex(unsigned int tidx) const {/*{{{*/
        assert(tidx < m_index.size());
        return m_index[tidx];
      }/*}}}*/
      int BasicSegStartT(unsigned int tidx) const { /*{{{*/
        assert(tidx < m_index.size());
        return m_segtree.StartT(m_index[tidx]);
      }/*}}}*/
      int BasicSegEndT(unsigned int tidx) const {/*{{{*/
        assert(tidx < m_index.size());
        return m_segtree.EndT(m_index[tidx]);
      }/*}}}*/
      int BasicSegLen(unsigned int sidx, unsigned int eidx) const {/*{{{*/
        assert(sidx <= eidx);
        assert(eidx < m_index.size());
        return m_segtree.EndT(m_index[eidx]) - m_segtree.StartT(m_index[sidx]) + 1;
      }/*}}}*/
      // supersegment related
      unsigned Width() const { return m_width; }
      unsigned Gran() const { return m_gran; }
      const deque<int>& Subseg(int l, int s) const { return m_subseg[l][s]; }
      float Subfrac(int l, int s, int g) const { return m_subfrac[l][s][g]; }
      const float* SegMean(int idx) const {/*{{{*/
        assert(idx >= 0 && idx < m_segtree.NumNode());
        if (idx < m_segtree.NumLeaf()) {
          return m_feat[idx];
        } else {
          return m_segmean[idx - m_segtree.NumLeaf()];
        }
      }/*}}}*/

      void DumpData() const;

    private:
      void Init();
      void GrowSub();
      void SplitSuperSeg(deque<int> *seg, deque<float> *frac,
                         TwoDimArray<float>& cum_feat);
      void ClearSeg();

      /* Data members */
      string m_featfname;
      string m_treefname;

      unsigned m_gran;
      unsigned m_width;

      //subseg[m_width][start_t][m_gran]
      deque<deque<deque<int> > > m_subseg;
      deque<deque<deque<float> > > m_subfrac;
      vector<int> m_index;
      vector<float *> m_segmean; // size() = m_segtree.NumIntNode()

      DenseFeature m_feat;
      SegTree m_segtree;
      float m_bseg_th;
      float m_superseg_th;

      ParmType m_type;
  };

} /* namespace DtwUtil */

#endif
