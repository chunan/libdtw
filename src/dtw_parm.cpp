#include "dtw_parm.h"
#include "utility.h"
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <algorithm>

namespace DtwUtil {

  void DtwParm::Init() {/*{{{*/

    m_gran = 0;
    m_width = 0;
    m_bseg_th = 0.0;
    m_superseg_th = 0.0;

    m_index.clear();
    m_segmean.clear();
    m_subseg.clear();
    m_subfrac.clear();

    m_type = FRAME;

  }/*}}}*/

  void DtwParm::LoadParm(const string feat_fname) {/*{{{*/

    /* Load feature file */
    m_featfname = feat_fname;
    m_feat.LoadFile(m_featfname);

  }/*}}}*/

  void DtwParm::LoadParm(const string feat_fname,/*{{{*/
                         const double bseg_ratio,
                         const string tree_fname) {

    m_type = SEGMENT;
    ClearSeg();

    /* Load feature file */
    m_featfname = feat_fname;
    m_feat.LoadFile(m_featfname);

    /* Load/Construct tree */
    if (tree_fname.empty()) {
      m_segtree.ConstructTree(m_feat);
    } else {
      m_treefname = tree_fname;
      m_segtree.Load_segtree(m_treefname);
    }

    /* Get basic segments */
    m_bseg_th = max(m_segtree.MergeMean() + bseg_ratio * m_segtree.MergeStd(),
                    0.0);
    m_segtree.GetBasicSeg(&m_index, m_bseg_th);
  }/*}}}*/

  void DtwParm::LoadParm(const string feat_fname,/*{{{*/
                         const double bseg_ratio,
                         const double superseg_ratio,
                         const unsigned width,
                         const unsigned gran,
                         const string tree_fname) {

    m_type = SUPERSEG;
    m_width = width;
    m_gran = gran;
    ClearSeg();

    if (m_width > m_gran) {
      fprintf(stderr, "Error: m_width(%d) > m_gran(%d).\n", m_width, m_gran);
      exit(1);
    }

    /* Load feature file */
    m_featfname = feat_fname;
    m_feat.LoadFile(m_featfname);

    if (m_feat.LT() < static_cast<int>(m_width)) {
      fprintf(stderr, "Length of feature (%d) < m_width (%d), set m_width to %d\n",
              m_feat.LT(), m_width, m_feat.LT());
    }

    /* Load/Construct tree */
    if (tree_fname.empty()) {
      m_segtree.ConstructTree(m_feat);
    } else {
      m_treefname = tree_fname;
      m_segtree.Load_segtree(m_treefname);
    }

    /* Compute thresholds */
    m_bseg_th = max(0.0, m_segtree.MergeMean() + bseg_ratio * m_segtree.MergeStd());
    m_superseg_th = max(0.0, m_segtree.MergeMean() + superseg_ratio * m_segtree.MergeStd());

    /* Get basic segments & supersegments */
    m_segtree.GetBasicSeg(&m_index, m_bseg_th);
    GrowSub();

  }/*}}}*/

  void DtwParm::DumpData() const {/*{{{*/
    string nameoftype[3] = {"FRAME", "SEGMENT", "SUPERSEG"};
    printf("Type:\n");
    printf("  %s\n", nameoftype[m_type].c_str());
    printf("FILENAMES:\n");
    if (!m_featfname.empty()) printf("  m_featfname = %s\n", m_featfname.c_str());
    if (!m_treefname.empty()) printf("  m_treefname = %s\n", m_treefname.c_str());
    printf("PARAMETER:\n");
    printf("  m_gran = %d\n m_width = %d\n", m_gran, m_width);
    //printf("DATA:m_feat\n");
    //m_feat.DumpData();
    printf("DATA:m_segtree\n");
    m_segtree.DumpData();
    printf("DATA:basic segment\n");
    for (unsigned i = 0; i < m_index.size(); i++) {
      printf("  basic_seg[%d](m_segtree idx = %d): %d ~ %d\n",
             i, m_index[i], BasicSegStartT(i), BasicSegEndT(i));
    }
    printf("DATA:m_segmean:\n  {");
    int segmean_count = 0;
    for (unsigned i = 0; i < m_segmean.size(); ++i) {
      if (m_segmean[i]) {
        printf(" %d", m_segtree.NumLeaf() + i);
        ++segmean_count;
      }
    }
    printf("} (%d/%d)\n", segmean_count, m_segtree.NumIntNode());
    printf("DATA:m_subseg: {\n");
    for (size_t w = 0; w < m_subseg.size(); ++w) {
      printf("  <width = %lu>:\n", w + 1);
      for (size_t s = 0; s < m_subseg[w].size(); ++s) {
        printf("   m_subseg[%lu][%lu] = (", w, s);
        for (size_t g = 0; g < m_subseg[w][s].size(); ++g) {
          printf(" %d", m_subseg[w][s][g]);
        }
        printf(")\n");
      }
    }
    printf("}\n");

  }/*}}}*/

  void DtwParm::SplitSuperSeg(deque<int> *seg,/*{{{*/
                              deque<float> *frac,
                              TwoDimArray<float>& cum_feat) {
    /* Find m_gran subsegments *//*{{{*/
    deque<int>::iterator itr_seg_tosplit;
    while (seg->size() < m_gran) {
      itr_seg_tosplit = max_element(seg->begin(), seg->end());
      /* Cannot split anymore. */
      if(*itr_seg_tosplit < m_feat.LT()) break;
      /* Can split. */
      int first_child = m_segtree.Child(*itr_seg_tosplit, 0);
      int second_child = m_segtree.Child(*itr_seg_tosplit, 1);
      *itr_seg_tosplit = second_child;
      seg->insert(itr_seg_tosplit, first_child);
    }/*}}}*/

    /* Not enough subsegments *//*{{{*/
    if(seg->size() < m_gran){
      int fr_s = m_segtree.StartT(seg->front());
      int fr_e = m_segtree.EndT(seg->back());

      seg->resize(m_gran);
      double per_share = double(fr_e - fr_s + 2 - 1e-3) / double(m_gran);
      for (unsigned k = 0; k < m_gran; k++)
        (*seg)[k] = fr_s + k * per_share;
    }
    /*}}}*/

    /* Calculate frac & m_segmean  */ /*{{{*/
    frac->resize(m_gran);
    double total_len = 0;
    /* For each subsegment in this supersegment */
    for (unsigned k = 0; k < m_gran; ++k) {
      int seg_idx = (*seg)[k]; // m_segtree index
      int segmean_idx = seg_idx - m_segtree.NumLeaf();
      (*frac)[k] = m_segtree.EndT(seg_idx) - m_segtree.StartT(seg_idx) + 1;
      total_len += (*frac)[k];
      /* Not a single frame segment and not allocated yet. */
      if (segmean_idx >= 0 && m_segmean[segmean_idx] == NULL) {
        /* Allocate space */
        m_segmean[segmean_idx] = new float [m_feat.LF()];
        /* Calculate mean */
        for (int f = 0; f < m_feat.LF(); ++f) {
          m_segmean[segmean_idx][f] = cum_feat(m_segtree.EndT(seg_idx), f);
        }
        if(m_segtree.StartT(seg_idx) != 0){
          int fr_s_1 = m_segtree.StartT(seg_idx) - 1;
          for (int f = 0; f < m_feat.LF(); ++f) {
            m_segmean[segmean_idx][f] -= cum_feat(fr_s_1, f);
          }
        }
        for (int f = 0; f < m_feat.LF(); f++) {
          m_segmean[segmean_idx][f] /= (*frac)[k];
        }
      }
    }
    for(unsigned k = 0; k < m_gran; k++){
      (*frac)[k] /= total_len;
    }
    /*}}}*/

  }/*}}}*/

  void DtwParm::GrowSub() {/*{{{*/

    /* Allocate m_segmean *//*{{{*/
    m_segmean.resize(m_segtree.NumIntNode(), NULL);
    /*}}}*/

    /* Allocate cum_feat *//*{{{*/
    TwoDimArray<float> cum_feat(m_feat.LT(), m_feat.LF());
    /*** Calculate cum_feat ***/
    for (int j = 0; j < m_feat.LF(); j++) {
      cum_feat(0, j) = m_feat(0, j);
    }
    for (int i = 1; i < m_feat.LT(); i++) {
      for (int j = 0; j < m_feat.LF(); j++) {
        cum_feat(i, j) = cum_feat(i - 1, j) + m_feat(i, j);
      }
    }
    /*}}}*/

    /* Allocate m_subseg & m_subfrac *//*{{{*/
    m_subseg.clear();
    m_subfrac.clear();
    m_subseg.resize(m_width);
    m_subfrac.resize(m_width);

    for (unsigned l = 0; l < m_width; l++) {
      m_subseg[l].resize(m_index.size() - l);
      m_subfrac[l].resize(m_index.size() - l);
    }
    /*}}}*/

    /* Find highest ancestor below threshold for every basic segment *//*{{{*/
    vector<int> highest_ancestor(m_index.size());
    for (unsigned t = 0; t < highest_ancestor.size(); ++t) {
      highest_ancestor[t] = m_index[t];
      int parent;
      while((parent = m_segtree.Parent(highest_ancestor[t])) > 0 &&
            m_segtree.MergeLoss(parent) < m_superseg_th) {
        highest_ancestor[t] = parent;
      }
      //printf("highest_ancestor[%d] = %d\n", t, highest_ancestor[t]);
    }/*}}}*/

    for (unsigned t = 0; t < m_index.size(); ++t) { // Starting segment time
      // Consists of l = i + 1 basic segments
      for (unsigned i = 0; i < m_width && t + i < m_index.size(); ++i) {
        int l = i + 1;

        // break if the newest basic segment does not belong to the same big tree
        if (highest_ancestor[t + i] != highest_ancestor[t]) break;

        // Init segment at hand /*{{{*/
        m_subseg[i][t].resize(l);
        for (int k = 0; k < l; k++) {
          m_subseg[i][t][k] = m_index[t + k];
        } /*}}}*/

        SplitSuperSeg(&m_subseg[i][t], &m_subfrac[i][t], cum_feat);

      }
    }
  }/*}}}*/

  void DtwParm::ClearSeg() {/*{{{*/
    if (m_type != SUPERSEG) return;
    /****** Free segmean ********/
    for (unsigned i = 0; i < m_segmean.size(); ++i) {
      if (m_segmean[i]) delete [] m_segmean[i];
    }
    m_segmean.clear();
    m_index.clear();
  }/*}}}*/

} /* namespace DtwUtil */
