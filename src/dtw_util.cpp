#include <climits>
#include <cmath>
#include <numeric>

#include "thread_util.h"
#include "dtw_util.h"
#include "utility.h"
#include "std_common.h"
#include "debug.h"


using StdCommonUtil::SnippetProfile;
using StdCommonUtil::SnippetProfileList;

namespace DtwUtil {

  Timer* timer = NULL;

  VectorDistFn DeterminePhiFn(string ext) {/*{{{*/

    VectorDistFn phi_fn = euclinorm;
    if (ext.compare("gp") == 0) {
      phi_fn = innernorm;
    } else if (ext.compare("mfc") == 0) {
      phi_fn = euclinorm;
    } else if (ext.compare("fbank") == 0) {
      phi_fn = euclinorm;
    } else if (ext.compare("plp") == 0) {
      phi_fn = euclinorm;
    } else {
      cerr << "unknown file extension: " << ext 
        << ", use euclinorm()\n";
    }

    return phi_fn;
  }/*}}}*/

  float innernorm(const float* a, const float* b, const int size) {/*{{{*/
    float ret = 0.0;
    for (int i = 0; i < size; i++) {
      ret += a[i] * b[i];
    }
    ret = -log(ret);
    return ret;
  }/*}}}*/

  float euclinorm(const float* a, const float* b, const int size) {/*{{{*/
    float ret = 0.0;
    for (int i = 0; i < size; i++) {
      float diff = a[i] - b[i];
      ret += diff * diff;
    }
    ret = sqrt(ret);
    return ret;
  }/*}}}*/


  void PrintPath(const vector<string>& Q_list,/*{{{*/
                 const vector<string>& D_list,
                 QDArray<vector<Path> >& path) {
    for (unsigned q = 0; q < Q_list.size(); ++q) {
      for (unsigned d = 0; d < D_list.size(); ++d) {
        if (path(q, d).empty()) continue;
        cerr << Q_list[q] << endl;
        cerr << D_list[d] << endl;
        cerr << path(q, d)[0].size() - 1 << endl;
        for (unsigned i = 1; i < path(q, d)[0].size(); ++i) {
          cerr << path(q, d)[0][i - 1].first + 1 << ' '
            << path(q, d)[0][i].first << ' '
            << path(q, d)[0][i - 1].second + 1 << ' '
            << path(q, d)[0][i].second << endl;
        }
      }
    }
  }/*}}}*/

  void PrintVecString(const vector<string>& list) {/*{{{*/
    cout << "index  name\n";
    for (unsigned i = 0; i < list.size(); ++i) {
      cout << setw(5) << i << "  " << list[i] << endl;
    }
  }/*}}}*/

  void PrintTwoDimFloat(TwoDimArray<float>& array) {/*{{{*/
    DUMP(__FUNCTION__);
    cout << fixed << noshowpos;
    cout << "   ";
    for (int c = 0; c < array.C(); ++c) {
      cout << setw(5) << c;
    }
    cout << endl;
    for (int r = 0; r < array.R(); ++r) {
      cout << noshowpos << setw(2) << r << ":";
      for (int c = 0; c < array.C(); ++c) {
        if (array(r, c) != float_inf) {
          cout << "  " << setw(3) << setprecision(0) << array(r, c);
        } else {
          cout << "  " << setw(3) << right << " .";
        }
      }
      cout << endl;
    }
    getc(stdin);
  }/*}}}*/

  void PrintTwoDimLast(TwoDimArray<IPair>& array) {/*{{{*/
    DUMP(__FUNCTION__);
    cout << "   ";
    for (int c = 0; c < array.C(); ++c) {
      cout << setw(5) << c;
    }
    cout << endl;
    for (int r = 0; r < array.R(); ++r) {
      cout << setw(2) << noshowpos << r << ":";
      for (int c = 0; c < array.C(); ++c) {
        cout << " " << showpos << array(r, c).first << array(r, c).second;
      }
      cout << endl;
    }
    getc(stdin);
  }/*}}}*/

  void PrintHeadDelBand(vector<float>& v) {/*{{{*/
    DUMP(__FUNCTION__);
    cout << fixed;
    for (unsigned q = 0; q < v.size(); ++q) {
      cout << setw(2) << q << ":" << setw(3) << setprecision(0) << v[q] << endl;
    }
    getc(stdin);
  }/*}}}*/



  void _LoadParm(const string &filename, DtwParm *parm,/*{{{*/
                 const DtwParm::ParmType type,
                 const double bseg_ratio,
                 const double superseg_ratio) {
    string base = filename;
    StripExtension(&base);

    string featfile = filename;
    string treefile = base + ".tb";

    if (type == DtwParm::FRAME) {
      parm->LoadParm(featfile);
    } else if (type == DtwParm::SEGMENT) {
      parm->LoadParm(featfile, bseg_ratio, treefile);
    } else { // SUPERSEG
      parm->LoadParm(featfile, bseg_ratio, superseg_ratio,
                     3, 3, treefile);
    }
    //parm->DumpData();
  }/*}}}*/

  void LoadParmList(const vector<string> &fnamelist, /*{{{*/
                    const FileDirExt &dirext,
                    vector<DtwParm> *parmlist,
                    const DtwParm::ParmType type,
                    const double bseg_ratio,
                    const double superseg_ratio) {
    DUMP(__FUNCTION__);
    /* Set number of Parms */
    parmlist->resize(fnamelist.size());
    for (unsigned i = 0; i < fnamelist.size(); ++i) {
      string fname;
      if (!dirext.dir.empty()) fname = dirext.dir + '/';
      fname += fnamelist[i];
      if (!dirext.ext.empty()) fname += '.' + dirext.ext;
      _LoadParm(fname, &(*parmlist)[i], type, bseg_ratio, superseg_ratio);
    }
  }/*}}}*/



  float SegPairwiseDist::measure_(int q, int d) {/*{{{*/
    const float* q_vec = qparm_->SegMean(q);
    const float* d_vec = dparm_->SegMean(d);
    assert(q_vec != NULL);
    assert(d_vec != NULL);
    return norm_(q_vec, d_vec, qparm_->Feat().LF());
  }/*}}}*/

  float FramePairwiseDist::measure_(int q, int d) {/*{{{*/
    const float* q_vec = (*q_feat_)[q_start_frm_ + q];
    const float* d_vec = (*d_feat_)[d_start_frm_ + d];
    assert(q_vec != NULL);
    assert(d_vec != NULL);
    return norm_(q_vec, d_vec, q_feat_->LF());
  }/*}}}*/




  /* Static variables for DtwRunner */
  float DtwRunner::del_ratio_ = 0.2;

  /* Static variables for FrameDtwRunner */
  unsigned FrameDtwRunner::nsnippet_ = 1;

  /* Static variables for SlopeConDtwRunner */
  int SlopeConDtwRunner::wq_ = 2;
  int SlopeConDtwRunner::wd_ = 2;

  /* Static variables for SegDtwRunner */
  int SegDtwRunner::wq_ = 3;
  int SegDtwRunner::wd_ = 3;
  float SegDtwRunner::rpen_ = 0.65;
  float SegDtwRunner::lratio_ = 3.0;
  unsigned SegDtwRunner::nsnippet_ = 1;

  /* Static variables for SegmentalDtwRunner */
  int SegmentalDtwRunner::adj_win = 15;
  int SegmentalDtwRunner::step_len = 15;
  unsigned SegmentalDtwRunner::nsnippet_ = 1;

  /* Static variables for FixFrameDtwRunner */
  int FixFrameDtwRunner::adj_win = 15;


  void DtwRunner::InitDtw(vector<float>* snippet_dist,/*{{{*/
                          vector<IPair>* snippet_bound,
                          vector<Path>* paths,
                          const DtwParm* q_parm,
                          const DtwParm* d_parm,
                          const IPair* qbound,
                          const IPair* dbound) {
    DUMP(__FUNCTION__);
    snippet_dist_ = snippet_dist;
    snippet_bound_ = snippet_bound;
    paths_ = paths;
    qparm_ = q_parm;
    dparm_ = d_parm;
    qbound_ = qbound;
    dbound_ = dbound;
  }/*}}}*/

  void DtwRunner::FindBestPaths(unsigned nsnippet,/*{{{*/
                                const float normalizer,
                                IPair drange) {

    if (drange.first < 0 || drange.first >= dL_) drange.first = 0;
    if (drange.second < 0 || drange.second  >= dL_) drange.second= dL_ - 1;

    vector<int> cand;
    _FindAllValley(&cand, drange);
    // cand sorted by increasing distance
#if 0
    if (cand.empty()) {
      cout << "No path found, score = \n";
      cout << score_;
    }
#endif

    unsigned num_fetched = 0;
    for (num_fetched = 0; num_fetched < cand.size(); ++num_fetched) {

      if (snippet_dist_->size() >= nsnippet) { // Full, eliminate the minimum

        float new_score = -score_(qL_ - 1, cand[num_fetched]) / normalizer;
        vector<float>::iterator itr_min =
          min_element(snippet_dist_->begin(), snippet_dist_->end());

        if (*itr_min > new_score) { // min-snippet > max-cand
          break;

        } else { // max-cand replace min-snippet
          int nth = itr_min - snippet_dist_->begin();
          *itr_min = new_score;
          if (snippet_bound_) {
            (*snippet_bound_)[nth] =
              IPair(dstart_ + root_(qL_ - 1, cand[num_fetched]),
                    dstart_ + cand[num_fetched]);
          }
          if (paths_) {
            _TracePath(&(*paths_)[nth], cand[num_fetched]);
          }
        }


      } else { // Not full, push back

        snippet_dist_->push_back(-score_(qL_ - 1, cand[num_fetched]) / normalizer);
        if (snippet_bound_) {
          snippet_bound_->push_back(
              IPair(dstart_ + root_(qL_ - 1, cand[num_fetched]),
                    dstart_ + cand[num_fetched]));
        }
        if (paths_) {
          paths_->push_back(Path());
          _TracePath(&paths_->back(), cand[num_fetched]);
        }

      } // if-else
    } // for

  }/*}}}*/

  struct lt_idx {/*{{{*/
    lt_idx(const TwoDimArray<float>& arr) : score_(arr) {}
    bool operator()(int i, int j) const {
      return score_.Entry(score_.R() - 1, i) < score_.Entry(score_.R() - 1, j);
    }
    const TwoDimArray<float>& score_;

  };/*}}}*/

  void DtwRunner::_FindAllValley(vector<int>* cand, IPair drange) {/*{{{*/

    if (drange.first < 0 || drange.first >= dL_) drange.first = 0;
    if (drange.second < 0 || drange.second  >= dL_) drange.second= dL_ - 1;

    cand->clear();
    const int q_last = qL_ - 1;

    int prev_root = -1;
    for (int d = drange.first; d <= drange.second; ++d) {

      if (score_(q_last, d) == float_inf) continue;

      if (prev_root != root_(q_last, d)) { // New root
        cand->push_back(d);
        prev_root = root_(q_last, d);
      } else { // In the same root tree
        if (score_(q_last, d) < score_(q_last, cand->back())) {
          cand->back() = d;
        }
      }

    }
    sort(cand->begin(), cand->end(), lt_idx(score_));
  }/*}}}*/

  void DtwRunner::_TracePath(Path* p, int d) {/*{{{*/
    p->clear();
    int len_q, len_d;
    int q;
    /* End pair */
    if ((len_d = lastL_(qL_ - 1, d).second) == 0) { // Deletion at end
      len_q = lastL_(qL_ - 1, d).first;
      q = qL_ - len_q - 1;
    } else { // No deletion at end
      q = qL_ - 1;
    }
    p->push_front(IPair(qstart_ + q, dstart_ + d));
    while (q >= 0 && d >= 0 && (len_q = lastL_(q, d).first) > 0) {
      len_d = lastL_(q, d).second;
      q -= len_q;
      d -= len_d;
      p->push_front(IPair(qstart_ + q, dstart_ + d));
    }
  }/*}}}*/



  void SegDtwRunner::DTW() {/*{{{*/
    DUMP(__FUNCTION__);

    qstart_ = qbound_ ? qbound_->first : 0;
    qend_ = qbound_ ? qbound_->second : -1;
    CheckStartEnd(*qparm_, &qstart_, &qend_, &qL_);
    dstart_ = dbound_ ? dbound_->first : 0;
    dend_ = dbound_ ? dbound_->second : -1;
    CheckStartEnd(*dparm_, &dstart_, &dend_, &dL_);

    /* Reallocate memory:
     * pdist.size = (parm_q.NumBasicSeg, parm_d.NumBasicSeg)
     * score_(paths_).size = (qL_, dL_)
     */
    pdist.Reset(qparm_, dparm_);



    score_.Resize(qL_, dL_);
    score_.Memfill(float_inf);
    root_.Resize(qL_, dL_);
    //root_.Memfill(-1);
    if (paths_) lastL_.Resize(qL_, dL_);

    MaxDelSeg();
    CalStartDel();
    CalScoreTable();
    CalEndDel();
    FindBestPaths(nsnippet_, qparm_->BasicSegLen(qstart_, qend_));
  }/*}}}*/

  /* Calculate the super segment distance that
   * - End at actual indices (q, d)
   * - Of length (len_q, len_d)
   */
  float SegDtwRunner::SegDist(const int q, const int len_q,/*{{{*/
                              const int d, const int len_d,
                              float upperbound) {
    assert(q >= 0 && q < qL_);
    assert(d >= 0 && d < dL_);
    if (upperbound < 0.0) {
      return float_inf;
    }
    const int lq = len_q - 1;
    const int ld = len_d - 1;
    if (q < lq || d < ld) {
      return float_inf;
    }
    if (qparm_->Subseg(lq, q - lq).empty() ||
        dparm_->Subseg(ld, d - ld).empty()) {
      return float_inf;
    }
    upperbound /= qparm_->BasicSegLen(q - lq, q);
    upperbound *= qparm_->Gran();
    float seg_dist = 0.0;
    for (unsigned g = 0; g < qparm_->Gran(); ++g) {
      int q_ssid = qparm_->Subseg(lq, q - lq)[g];
      int d_ssid = dparm_->Subseg(ld, d - ld)[g];
      float q_ssfrac = qparm_->Subfrac(lq, q - lq, g);
      float d_ssfrac = dparm_->Subfrac(ld, d - ld, g);

      float comp_pen = exp(rpen_ * abs(q_ssfrac - d_ssfrac));
      float vec_dist = pdist(q_ssid, d_ssid);
      seg_dist += comp_pen * vec_dist;
      if (seg_dist > upperbound) return float_inf;
    }
    seg_dist *= qparm_->BasicSegLen(q - lq, q);
    seg_dist /= qparm_->Gran();
    return seg_dist;
  }/*}}}*/

  void SegDtwRunner::CheckStartEnd(const DtwParm &parm,/*{{{*/
                                   int *s, int *e, int *l) {
    if (*s < 0) *s = 0;
    if (*e < 0 || *e >= static_cast<int>(parm.NumBasicSeg())) {
      *e = parm.NumBasicSeg() - 1;
    }
    *l = *e - *s + 1;
  }/*}}}*/

  void SegDtwRunner::MaxDelSeg() {/*{{{*/
    DUMP(__FUNCTION__);
    /* Calculate n_head_del and n_tail_del */
    n_head_del = n_tail_del = 0;
    if (del_ratio_ > 0.0) {
      float qL_T = qparm_->BasicSegLen(qstart_, qend_);
      for (int n = 0; n < qL_; ++n) {
        n_head_del = n;
        float ratio = qparm_->BasicSegLen(qstart_, qstart_ + n) / qL_T;
        if (ratio > del_ratio_) break;
      }
      for (int n = 0; n < qL_; ++n) {
        n_tail_del = n;
        float ratio = qparm_->BasicSegLen(qend_ - n, qend_) / qL_T;
        if (ratio > del_ratio_) break;
      }
      //printf("n_head_del = %d, n_tail_del = %d\n", n_head_del, n_tail_del);
    }
    head_del_band.resize(n_head_del);
    tail_del_band.Resize(n_tail_del, dL_);
  }/*}}}*/

  void SegDtwRunner::CalStartDel() {/*{{{*/
    /* When d = 0, deletion penalty stored in head_del_band */
    for (int q = 0; q < n_head_del; ++q) {
      head_del_band[q] = SegDist(qstart_ + q, 1, dstart_, 1);
    }
    for (int q = 0; q < n_head_del; ++q) {
      head_del_band[q] += head_del_band[q - 1];
    }
    /* d > 0, deletion penalty stored in score_.
     * Note that the penalty \sum SegDist (q1~qr, d) is
     * stored in score_(qr, d - 1), and that's why we
     * need head_del_band for d = 0.
     */
    for (int d = 1; d < dL_; ++d) {
      for (int q = 0; q < n_head_del; ++q) {
        score_(q, d - 1) = SegDist(qstart_ + q, 1, dstart_ + d, 1);
        root_(q, d - 1) = d;
        if (paths_) lastL_(q, d - 1) = IPair(-1, -1);
      }
      for (int q = 1; q < n_head_del; ++q) {
        score_(q, d - 1) +=  score_(q - 1, d - 1);
      }
    }
  }/*}}}*/

  void SegDtwRunner::CalScoreTable() {/*{{{*/
    DUMP(__FUNCTION__);
    /* q = 0 */
    for (int d = 0; d < dL_; ++d) {
      _BestScore(0, d);
    }
    /* d = 0 */
    for (int q = 1; q < qL_; ++q) {
      _BestScore(q, 0);
    }
    /* q > 1 && d > 1 */
    for (int q = 1; q < qL_; ++q) {
      for (int d = 1; d < dL_; ++d) {
        _BestScore(q, d);
      }
    }
  }/*}}}*/

  /* (q, d) are the relative index in search region */
  void SegDtwRunner::_BestScore(const int q, const int d) {/*{{{*/
    // Actual index in DtwParm
    const int q_act = qstart_ + q;
    const int d_act = dstart_ + d;
    for (int wq = 1; wq <= static_cast<int>(qparm_->Width()) && q - wq >= -1; ++wq) {
      int q_lenT = qparm_->BasicSegLen(q_act - wq + 1, q_act);
      int d_lb_lenT = ceil(q_lenT / lratio_);
      int d_ub_lenT = /* Unlimited d_ub_lenT at beginning and end */
        (q - wq == -1 || q == qL_ - 1) ? INT_MAX : floor(q_lenT * lratio_);
      for (int wd = 1; wd <= static_cast<int>(dparm_->Width()) && d - wd >= -1; ++wd) {
        int d_lenT = dparm_->BasicSegLen(d_act - wd + 1, d_act);
        if (d_lenT < d_lb_lenT) continue;
        if (d_lenT > d_ub_lenT) break;
        /* Determine the score of previous path */
        float prev_path_score;
        if (q - wq == -1) { // Q head
          prev_path_score = 0.0;
        } else if (d - wd != -1) { // Q not head, D not head
          prev_path_score = score_(q - wq, d - wd);
        } else { // Q not head, D head
          prev_path_score = q - wq < static_cast<int>(head_del_band.size())
            ? head_del_band[q - wq] : float_inf;
        }
        if (prev_path_score == float_inf) continue;
        float dist_upperbnd = score_(q, d) - prev_path_score;
        float new_seg_dist = SegDist(q_act, wq, d_act, wd, dist_upperbnd);
        if (new_seg_dist != float_inf) {
          score_(q, d) = new_seg_dist + prev_path_score;
          root_(q, d) = (q - wq == -1 || d - wd == -1) ?
            d - wd + 1 : root_(q - wq, d - wd);
          if (paths_) lastL_(q, d) = IPair(wq, wd);
        }
      }
    }
  }/*}}}*/

  void SegDtwRunner::CalEndDel() {/*{{{*/
    for (int d = 0; d < dL_; ++d) {
      // In case there's no need to calculate deletion
      int n;
      for (n = n_tail_del; n > 0; --n) {
        if (score_(qL_ - 1 - n, d) != float_inf) break;
      }
      int d_act = dstart_ + d;
      // Accumulated distance penalty from tail
      float accum_del_pen = 0.0;
      for (int len = 1; len <= n; ++len) {
        int q_act = qend_ - len;
        accum_del_pen += SegDist(q_act, 1, d_act, 1);
        float new_score = accum_del_pen + score_(qL_ - 1 - len, d);
        if (new_score != float_inf && new_score < score_(qL_ - 1, d)) {
          score_(qL_ - 1, d) = new_score;
          root_(qL_ - 1, d) = root_(qL_ - 1 - len, d);
          if (paths_) lastL_(qL_ - 1, d) = IPair(len, 0);
        }
      }
    }
  }/*}}}*/



  void FrameDtwRunner::DTW() {/*{{{*/
    DUMP(__FUNCTION__);

    qstart_ = qbound_ ? qbound_->first : 0;
    qend_ = qbound_ ? qbound_->second : -1;
    CheckStartEnd(*qparm_, &qstart_, &qend_, &qL_);

    dstart_ = dbound_ ? dbound_->first : 0;
    dend_ = dbound_ ? dbound_->second : -1;
    CheckStartEnd(*dparm_, &dstart_, &dend_, &dL_);

    /* Reallocate memory:
     * pdist.size = (parm_q.NumBasicSeg, parm_d.NumBasicSeg)
     * score_(paths_).size = (qL_, dL_)
     */
    pdist.Reset(&qparm_->Feat(), &dparm_->Feat(), qstart_, dstart_, qL_, dL_);
    score_.Resize(qL_, dL_);
    score_.Memfill(float_inf);
    root_.Resize(qL_, dL_);
    if (paths_) lastL_.Resize(qL_, dL_);

    MaxDelFrame();
    CalScoreTable();
    FindBestPaths(nsnippet_, qL_);
  }/*}}}*/

  void FrameDtwRunner::CheckStartEnd(const DtwParm &parm,/*{{{*/
                                     int *s, int *e, int *l) {
    if (*s < 0) *s = 0;
    if (*e < 0 || *e >= static_cast<int>(parm.Feat().LT())) {
      *e = parm.Feat().LT() - 1;
    }
    *l = *e - *s + 1;
  }/*}}}*/



  void SlopeConDtwRunner::CalScoreTable() {/*{{{*/
    DUMP(__FUNCTION__);

    /*
     *                    d    /wq\/wq\/1\
     *                    +-------+----+-+
     * max_d_when_q_is_0  |       +----+:+
     *                  \ |  +----+::::::|  min_d_when_q_is_end
     *                   `+--+:::::::::::| /
     *                    |:::::::::+----+'
     *                    |::::+----+    |
     *                    +----+---------+--> q
     *                     \wq/ \wq/
     */

    // max_d_when_q_is_0  = (dL_ - 1) - ceil((qL_ - 1) / wq_);
    const int qend = qL_ - 1;
    int d_hat = qend % wq_ ? qend / wq_ + 1 : qend / wq_;
    int max_d_when_q_is_0 = dL_ - 1 - d_hat;
    //int min_d_when_q_is_end = qL_ / wq_;

    if (max_d_when_q_is_0 >= -1) {

      // q == 0
      for (int d = 0; d <= max_d_when_q_is_0 + 1; ++d) {
        score_(0, d) = pdist(0, d);
        root_(0, d) = d;
        if (paths_) lastL_(0, d) = IPair(1, 1);
      }

      // d == 0
      for (int q = 1; q < wq_ && q < qL_; ++q) {
        score_(q, 0) = score_(q - 1, 0) + pdist(q, 0);
        root_(q, 0) = 0;
        if (paths_) lastL_(q, 0) = IPair(q + 1, 1);
      }


      for (int d = 1; d < dL_; ++d) {
        int d_reverse = dL_ - 1 - d;
        int begin = max(1, qend - wq_ * d_reverse);
        int end = min(qL_, wq_ * (d + 1));
        for (int q = begin; q < end; ++q) {
          _BestScoreInterior(q, d);
        }
      }

    }

  }/*}}}*/

  void SlopeConDtwRunner::_BestScoreInterior(const int q, const int d) {/*{{{*/

    assert(q > 0 && q < qL_);
    assert(d > 0 && d < dL_);

    float fixed_path_dist;

    /* 1 frame in Q : * frames in D */
    int max_wd = min(wd_, d); // score_(q - 1, d - wd) must be inbound

    // maximal possible wd
    // Inbound: q' < wq_ * (d' + 1)
    if (q - 1 >= wq_ * (d - max_wd + 1)) { // Lower-right out of bound
      max_wd = ((d + 1) * wq_ - q + 1) / wq_;
    }

    // Try last_fixed_path = 1 frame in Q : wd frames in D.
    int wd;
    for (wd = 1, fixed_path_dist = 0.0; wd <= max_wd ; ++wd) {
      fixed_path_dist += pdist(q, d - wd + 1); // sum_{j = d-wd+1}^{d} d(q,j)
      float new_score = fixed_path_dist + score_(q - 1, d - wd);
      // path coming from (q - 1, d - wd)
      if (new_score < score_(q, d)) {
        score_(q, d) = new_score;
        root_(q, d) = root_(q - 1, d - wd);
        if (paths_) lastL_(q, d) = IPair(1, wd);
      }
    }



    /* * frames in Q : 1 frame in D */
    int max_wq = min(wq_, q + 1); // possibly q - wq == -1

    // maximal possible wq
    if (max_wq < q + 1) { // Cannot reach first frame in Q
      for (; max_wq >= 0; --max_wq) {
        if (score_(q - max_wq, d - 1) != float_inf) break;
      }
    }

    // Try last_fixed_path = wq frames in Q : 1 frame in D.
    int wq;
    for (wq = 2, fixed_path_dist = pdist(q, d); wq <= max_wq; ++wq) {
      fixed_path_dist += pdist(q - wq + 1, d);
      float new_score;
      if (q - wq == -1) { // prev = 0;
        new_score = fixed_path_dist;
      } else { // has_prev
        new_score = fixed_path_dist + score_(q - wq, d - 1);
      }
      // coming from (q - wq, d - 1) (if not the first path)
      if (new_score < score_(q, d)) {
        score_(q, d) = new_score;
        root_(q, d) = (q - wq == -1) ? d : root_(q - wq, d - 1);
        if (paths_) lastL_(q, d) = IPair(wq, 1);
      }
    }

  }/*}}}*/



  void FreeFrameDtwRunner::CalScoreTable() {/*{{{*/
    // q == 0
    for (int d = 0; d < dL_; ++d) {
      score_(0, d) = pdist(0, d);
      root_(0, d) = d;
      if (paths_) lastL_(0, d) = IPair(1, 1);
    }

    // d == 0
    for (int q = 1; q < qL_; ++q) {
      score_(q, 0) = score_(q - 1, 0) + pdist(q, 0);
      root_(q, 0) = 0;
      if (paths_) lastL_(q, 0) = IPair(q + 1, 1);
    }

    // interior points
    for (int d = 1; d < dL_; ++d) {
      for (int q = 1; q < qL_; ++q) {
        // (-1, -1)
        score_(q, d) = score_(q - 1, d - 1);
        root_(q, d) = root_(q - 1, d - 1);
        if (paths_) lastL_(q, d) = IPair(1, 1);

        // (-1, 0)
        if (score_(q - 1, d) < score_(q, d)) {
          score_(q, d) = score_(q - 1, d);
          root_(q, d) = root_(q - 1, d);
          if (paths_) lastL_(q, d) = IPair(lastL_(q - 1, d).first + 1, 1);
        }

        // (0, -1)
        if (score_(q, d - 1) < score_(q, d)) {
          score_(q, d) = score_(q, d - 1);
          root_(q, d) = root_(q, d - 1);
          if (paths_) lastL_(q, d) = IPair(1, lastL_(q, d - 1).second + 1);
        }
        score_(q, d) += pdist(q, d);
      }
    }

  }/*}}}*/



  void FixFrameDtwRunner::CalScoreTable() {/*{{{*/

    float slope = static_cast<float>(dL_) / qL_;

    // (0, 0)
    score_(0, 0) = pdist(0, 0);
    root_(0, 0) = 0;

    // q == 0
    for (int d = 1; d <= adj_win && d < dL_; ++d) {
#ifdef FREEHEAD
      score_(0, d) = pdist(0, d);
      root_(0, d) = d;
      if (paths_) lastL_(0, d) = IPair(1, 1);
#else
      score_(0, d) = score_(0, d - 1) + pdist(0, d);
      root_(0, d) = 0;
      if (paths_) lastL_(0, d) = IPair(1, d + 1);
#endif
    }

    // d == 0
    for (int q = 1; q <= adj_win && q < qL_; ++q) {
      score_(q, 0) = score_(q - 1, 0) + pdist(q, 0);
      root_(q, 0) = 0;
      if (paths_) lastL_(q, 0) = IPair(q + 1, 1);
    }

    // interior points
    for (int q = 1; q < qL_; ++q) {
      int d_start = max<int>(1, slope * (q - adj_win));
      int d_last = min<int>(dL_ - 1, adj_win + q * slope);
      for (int d = d_start; d <= d_last; ++d) {

        // (-1, -1)
        score_(q, d) = score_(q - 1, d - 1);
        root_(q, d) = root_(q - 1, d - 1);
        if (paths_) lastL_(q, d) = IPair(1, 1);

        // (-1, 0)
        if (score_(q - 1, d) < score_(q, d)) {
          score_(q, d) = score_(q - 1, d);
          root_(q, d) = root_(q - 1, d);
          if (paths_) lastL_(q, d) = IPair(lastL_(q - 1, d).first + 1, 1);
        }

        // (0, -1)
        if (score_(q, d - 1) < score_(q, d)) {
          score_(q, d) = score_(q, d - 1);
          root_(q, d) = root_(q, d - 1);
          if (paths_) lastL_(q, d) = IPair(1, lastL_(q, d - 1).second + 1);
        }
        score_(q, d) += pdist(q, d);
      }
    }

  }/*}}}*/



  void SegmentalDtwRunner::DTW() {/*{{{*/
    DUMP(__FUNCTION__);

    qstart_ = qbound_ ? qbound_->first : 0;
    qend_ = qbound_ ? qbound_->second : -1;
    CheckStartEnd(*qparm_, &qstart_, &qend_, &qL_);

    dstart_ = dbound_ ? dbound_->first : 0;
    dend_ = dbound_ ? dbound_->second : -1;
    CheckStartEnd(*dparm_, &dstart_, &dend_, &dL_);

    /* Reallocate memory:
     * pdist.size = (parm_q.NumBasicSeg, parm_d.NumBasicSeg)
     * score_(paths_).size = (qL_, dL_)
     */
    pdist.Reset(&qparm_->Feat(), &dparm_->Feat(), qstart_, dstart_, qL_, dL_);
    score_.Resize(qL_, dL_);
    score_.Memfill(float_inf);
    root_.Resize(qL_, dL_);
    if (paths_) lastL_.Resize(qL_, dL_);

    MaxDelFrame();

    for (int d = 0; d + qL_ - 1 - adj_win < dL_ && d < dL_; d += step_len) {
      int d_targ_end = d + qL_ - 1;
      CalScoreTable(d);
      FindBestPaths(nsnippet_, qL_,
                    IPair(d_targ_end - adj_win, d_targ_end + adj_win));
    }
  }/*}}}*/

#define FREEHEAD
  void SegmentalDtwRunner::CalScoreTable(int d0) {/*{{{*/

    /* (0, d0) */
    score_(0, d0) = pdist(0, d0);
    root_(0, d0) = d0;
    if (paths_) lastL_(0, d0) = IPair(1, 1);

    /* (0, d0 + 1 ~ d0 + adj_win) */
    for (int d = d0 + 1; d <= d0 + adj_win && d < dL_; ++d) {
#ifdef FREEHEAD
      score_(0, d) = pdist(0, d);
      root_(0, d) = d;
      if (paths_) lastL_(0, d) = IPair(1, 1);
#else
      score_(0, d) = score_(0, d - 1) + pdist(0, d);
      root_(0, d) = d0;
      if (paths_) lastL_(0, d) = IPair(1, d - d0 + 1);
#endif
    }

    /* q \in [1, adj_win] */
    for (int q = 1; q <= adj_win && q < qL_; ++q) {

      int last_d = d0 + q + adj_win;

      // (q, d0)
      score_(q, d0) = score_(q - 1, d0) + pdist(q, d0);
      root_(q, d0) = d0;
      if (paths_) lastL_(q, d0) = IPair(q + 1, 1);

      // (q, d0+1 ~ last_d - 1)
      for (int d = d0 + 1; d < last_d && d < dL_; ++d) {
        ThreewayBestScore(q, d);
      }

      // (q, last_d)
      if (last_d < dL_) TwowayDallBestScore(q, last_d);

    }

    /* q \in [adj_win + 1, ... ] */
    for (int q = adj_win + 1; q < qL_; ++q) {

      int start_d = d0 + q - adj_win;
      int last_d = d0 + q + adj_win;

      // (q, start_d)
      TwowayQallBestScore(q, start_d);

      // (q, start_d + 1 ~ last_d - 1)
      for (int d = start_d + 1; d < last_d && d < dL_; ++d) {
        ThreewayBestScore(q, d);
      }

      // (q, last_d)
      if (last_d < dL_) TwowayDallBestScore(q, last_d);

    }

  }/*}}}*/

  void SegmentalDtwRunner::ThreewayBestScore(int q, int d) {/*{{{*/
    /* (-1, -1) */
    score_(q, d) = score_(q - 1, d - 1);
    root_(q, d) = root_(q - 1, d - 1);
    if (paths_) lastL_(q, d) = IPair(1, 1);

    /* (-1, 0) */
    if (score_(q - 1, d) < score_(q, d)) {
      score_(q, d) = score_(q - 1, d);
      root_(q, d) = root_(q - 1, d);
      if (paths_) lastL_(q, d) = IPair(lastL_(q - 1, d).first + 1, 1);
    }

    /* (0, -1) */
    if (score_(q, d - 1) < score_(q, d)) {
      score_(q, d) = score_(q, d - 1);
      root_(q, d) = root_(q, d - 1);
      if (paths_) lastL_(q, d) = IPair(1, lastL_(q, d - 1).second + 1);
    }
    score_(q, d) += pdist(q, d);
  }/*}}}*/

  void SegmentalDtwRunner::TwowayDallBestScore(int q, int d) {/*{{{*/

    /* (-1, -1) */
    if (score_(q - 1, d - 1) < score_(q, d - 1)) {
      score_(q, d) = score_(q - 1, d - 1);
      root_(q, d) = root_(q - 1, d - 1);
      if (paths_) lastL_(q, d) = IPair(1, 1);

    } else { /* (0, -1) */
      score_(q, d) = score_(q, d - 1);
      root_(q, d) = root_(q, d - 1);
      if (paths_) lastL_(q, d) = IPair(1, lastL_(q, d - 1).second + 1);
    }
    score_(q, d) += pdist(q, d);
  }/*}}}*/

  void SegmentalDtwRunner::TwowayQallBestScore(int q, int d) {/*{{{*/

    /* (-1, -1) */
    if (score_(q - 1, d - 1) < score_(q - 1, d)) {
      score_(q, d) = score_(q - 1, d - 1);
      root_(q, d) = root_(q - 1, d - 1);
      if (paths_) lastL_(q, d) = IPair(1, 1);

    } else { /* (-1, 0) */
      score_(q, d) = score_(q - 1, d);
      root_(q, d) = root_(q - 1, d);
      if (paths_) lastL_(q, d) = IPair(lastL_(q - 1, d).first + 1, 1);
    }
    score_(q, d) += pdist(q, d);
  }/*}}}*/





  bool QDDtwBatchParamSet::InstallDtwRunner(DtwRunner* runner) {/*{{{*/

    DUMP(__FUNCTION__);

    UPair* ticket = dispatcher_->GetObjPtr();
    if (ticket) {
      int qidx = ticket->first;
      int didx = ticket->second;

      vector<float>* snippet_dist = &(*snippet_dist_)(qidx, didx);
      vector<IPair>* snippet_bound = &snippet_bound_(qidx, didx);
      vector<Path>* paths = paths_ ? &(*paths_)(qidx, didx) : NULL;
      DtwParm* q_parm = &(*query_parm_list_)[qidx];
      DtwParm* d_parm = &(*doc_parm_list_)[didx];
      IPair* qboundary = vqboundary_ ? &(*vqboundary_)[qidx] : NULL;
      IPair* dboundary = vdboundary_ ? &(*vdboundary_)(qidx, didx) : NULL;

      runner->InitDtw(snippet_dist, snippet_bound , paths, q_parm, d_parm, qboundary, dboundary);
    }
    return ticket != NULL;
  }/*}}}*/



  bool QSDtwBatchParamSet::InstallDtwRunner(DtwRunner* runner) {/*{{{*/

    UPair* ticket = dispatcher_->GetObjPtr();
    if (ticket) {
      int qidx = ticket->first;
      int sidx = ticket->second;
      const SnippetProfile& pf = (*candidates_)[qidx].GetProfile(sidx);
      int didx = pf.Didx();

      vector<float>* snippet_dist = &snippet_dist_[qidx][sidx];
      vector<IPair>* snippet_bound = &snippet_bound_[qidx][sidx];
      DtwParm* q_parm = &(*qry_parm_list_)[qidx];
      DtwParm* d_parm = &(*doc_parm_list_)[didx];
      const IPair* qbound = vqboundary_ ? &(*vqboundary_)[qidx] : NULL;
      const IPair* dbound = &pf.Boundary();
      vector<Path>* paths = NULL;

      runner->InitDtw(snippet_dist, snippet_bound, paths,
                      q_parm, d_parm, qbound, dbound);
    }

    return ticket != NULL;

  }/*}}}*/

  void QSDtwBatchParamSet::CollectResult(bool do_sort) {/*{{{*/

    for (unsigned qidx = 0; qidx < candidates_->size(); ++qidx) {
      for (unsigned sidx = 0; sidx < (*candidates_)[qidx].size(); ++sidx) {

        SnippetProfile& snippet = (*candidates_)[qidx].ProfileRef(sidx);
        if (snippet_dist_[qidx][sidx].empty()) {
          snippet.ScoreRef() = -float_inf;
        } else {
          snippet.ScoreRef() = snippet_dist_[qidx][sidx][0];
          snippet.BoundaryRef() = snippet_bound_[qidx][sidx][0];
        }

      } /* for sidx */

      if (do_sort) (*candidates_)[qidx].Sort();

    } /* for qidx */

  }/*}}}*/



  bool PrfDtwBatchParamSet::InstallDtwRunner(DtwRunner* runner) {/*{{{*/

    UTriple* ticket = dispatcher_->GetObjPtr();
    if (ticket) {
      const int qidx = ticket->qidx;
      const int sidx = ticket->sidx;
      const int pidx = ticket->pidx;

      const SnippetProfile& pseu_region = (*candidates_)[qidx].GetProfile(pidx);
      const SnippetProfile& hypo_region = (*candidates_)[qidx].GetProfile(sidx);
      assert(pidx != sidx);

      vector<float>* snippet_dist = &output_dist_[qidx][sidx][pidx];
      vector<IPair>* snippet_bound = NULL;
      vector<Path>* paths = NULL;
      const DtwParm* q_parm = &(*doc_parm_list_)[pseu_region.Didx()];
      const DtwParm* d_parm = &(*doc_parm_list_)[hypo_region.Didx()];
      const IPair* qboundary = &pseu_region.Boundary();
      const IPair* dboundary = &hypo_region.Boundary();

      runner->InitDtw(snippet_dist, snippet_bound, paths,
                      q_parm, d_parm, qboundary, dboundary);
    }

    return ticket != NULL;
  }/*}}}*/

  void PrfDtwBatchParamSet::IntegratePrfDist(/*{{{*/
      const vector<float>& prf_weight) {

    for (unsigned qidx = 0; qidx < candidates_->size(); ++qidx) {

      int n_not_found = 0;

      /* push distance zero for each sidx < (*num_prf_)[qidx] */
      for (unsigned pidx = 0; pidx < (*num_prf_)[qidx]; ++pidx)
        output_dist_[qidx][pidx][pidx].push_back(0.0);

      /* for each snippet in qidx list */
      for (unsigned sidx = 0; sidx < (*candidates_)[qidx].size(); ++sidx) {

        SnippetProfile& snippet = (*candidates_)[qidx].ProfileRef(sidx);
        vector<vector<float> >& dist_vec = output_dist_[qidx][sidx];

        if (dist_vec[0].size() < (*num_prf_)[qidx]) {
          /* Colloct them in a single vector, namely, output_dist_[qidx][sidx][0] */
          if (dist_vec[0].empty()) {
            dist_vec[0].push_back(-float_inf);
            cout << "PRF not found: (pidx, sidx) = (0 , " << sidx << ")\n";
          }
          for (unsigned pidx = 1; pidx < (*num_prf_)[qidx]; ++pidx) {
            if (dist_vec[pidx].empty()) {// Cannot find path with PRF(pidx)
              dist_vec[0].push_back(-float_inf);
              cout << "PRF not found: (pidx, sidx) = ("
                << pidx << ", " << sidx << ")\n";
            } else {
              dist_vec[0].push_back(dist_vec[pidx][0]);
            }
          }
        }
        //cout << "DistVec(" << qidx << ", " << sidx << ") = " << dist_vec[0] << endl;

        // Calculate innerproduct: FIXME: += or =
        snippet.ScoreRef() = inner_product(
            dist_vec[0].begin(), dist_vec[0].end(), prf_weight.begin(), 0.0f);
        if (snippet.ScoreRef() == -float_inf)
          n_not_found++;

      } // end for sidx

      //(*candidates_)[qidx].Sort();
      //(*candidates_)[qidx].Normalize();
      if (n_not_found > 0) {
        cerr << "Warning: PrfDtwBatchParamSet::IntegratePrfDist():"
          << "qidx = " << qidx << ", not found = " << n_not_found << endl;
        //(*candidates_)[qidx].Resize(-n_not_found);
      }

    } // end for qidx

  }/*}}}*/


  void InitPrfDispatcher(Dispatcher<UTriple>* prf_dispatcher, /*{{{*/
                         const vector<SnippetProfileList>& candidates,
                         const vector<unsigned>& num_prf) {

    prf_dispatcher->Clear();

    for (unsigned qidx = 0; qidx < candidates.size(); ++qidx) {
      for (unsigned sidx = 0; sidx < candidates[qidx].size(); ++sidx) {
        for (unsigned pidx = 0; pidx < num_prf[qidx]; ++pidx) {
          if (sidx == pidx) continue;
          prf_dispatcher->Push(UTriple(qidx, sidx, pidx));
        }
      }
    }

  }/*}}}*/

  void* DtwManager::Run() {/*{{{*/
    // Try getting new ticket and install runner
    while (paramset_->InstallDtwRunner(runner_)) {
      runner_->DTW();
    }
    return NULL;
  }/*}}}*/

} /* namespace DtwUtil */
