#ifndef __DTW_UTIL_H__
#define __DTW_UTIL_H__

#include <limits>
#include <vector>
#include <cmath>
#include "thread_util.h"
#include "std_common.h"
#include "utility.h"
#include "dtw_parm.h"

#define float_inf std::numeric_limits<float>::infinity()

using StdCommonUtil::QDArray;
using StdCommonUtil::IPair;
using StdCommonUtil::UPair;
using StdCommonUtil::FileDirExt;
using StdCommonUtil::SnippetProfileList;


namespace DtwUtil {


  typedef float (*VectorDistFn)(const float*, const float*, const int);

  float innernorm(const float* a, const float* b, const int size);

  float euclinorm(const float* a, const float* b, const int size);

  typedef deque<IPair> Path;

  extern Timer* timer;


  VectorDistFn DeterminePhiFn(string ext);


  void PrintPath(const vector<string>& Q_list,
                 const vector<string>& D_list,
                 QDArray<vector<Path> >& path);

  void PrintVecString(const vector<string>& list);

  void PrintTwoDimFloat(TwoDimArray<float>& array);

  void PrintTwoDimLast(TwoDimArray<IPair>& array);

  void PrintHeadDelBand(vector<float>& v);


  void LoadParmList(const vector<string> &fnamelist,
                    const FileDirExt &dirext,
                    vector<DtwParm> *parmlist,
                    const DtwParm::ParmType type,
                    const double bseg_ratio = 0.15,
                    const double superseg_ratio = 4.0);



  template<class _Tp>
    void PrintCalCount(vector<_Tp>& runner) {/*{{{*/
      unsigned long tot_cal_count = 0;
      unsigned long n_digit_decimal = 1;
      unsigned  ndigits = 1;
      for (unsigned t = 0; t < runner.size(); ++t) {
        tot_cal_count += runner[t].PDist().CalCount();
      }
      while (tot_cal_count > n_digit_decimal) {
        ++ndigits;
        n_digit_decimal *= 10;
      }
      for (unsigned t = 0; t < runner.size(); ++t) {
        cout << "Total #call norm() [" << t << "] = "
          << setw(ndigits) << runner[t].PDist().CalCount() << endl;
      }
      cout << "Total #call norm() [*] = " << setw(ndigits) << tot_cal_count << endl;
    }/*}}}*/


  /* Interface: pairwise distance */
  class PairwiseDist {/*{{{*/
    public:
      PairwiseDist(VectorDistFn norm) {
        cal_count = 0;
        norm_ = norm;
      }
      virtual float operator()(const int q, const int d) {
        assert(q >= 0 && q < table_.R());
        assert(d >= 0 && d < table_.C());
        if (ready_[q][d] == false) {
          ++cal_count;
          table_(q, d) = measure_(q, d);
          ready_[q][d] = true;
        }
        return table_(q, d);
      }
      virtual unsigned long CalCount() const { return cal_count; }
    protected:
      virtual float measure_(int q, int d) = 0;
      /* data member */
      TwoDimArray<float> table_;
      QDArray<bool> ready_;
      unsigned long cal_count;
      VectorDistFn norm_;
  };/*}}}*/



  /* Pairwise distance of any two segments (two nodes in segtree) */
  class SegPairwiseDist : public PairwiseDist {/*{{{*/
    public:
      SegPairwiseDist(VectorDistFn norm) : PairwiseDist(norm) {
        qparm_ = dparm_ = NULL;
      }
      void Reset(const DtwParm* q_parm, const DtwParm* d_parm) {
        qparm_ = q_parm;
        dparm_ = d_parm;
        table_.Resize(qparm_->Segtree().NumNode(),
                      dparm_->Segtree().NumNode());
        ready_.Resize(qparm_->Segtree().NumNode(),
                      dparm_->Segtree().NumNode());
        ready_.memfill(false);
      }
    protected:
      virtual float measure_(int q, int d); /* Actual indices */
      const DtwParm *qparm_;
      const DtwParm *dparm_;
  };/*}}}*/


  /* Pairwise distance of any two frames */
  class FramePairwiseDist : public PairwiseDist {/*{{{*/
    public:
      FramePairwiseDist(VectorDistFn norm) : PairwiseDist(norm) {
        q_feat_ = d_feat_ = NULL;
      }
      void Reset(const DenseFeature* q_feat, const DenseFeature* d_feat,
                 int q_s, int d_s, int q_len, int d_len) {
        assert(q_s >= 0 && q_len > 0 && q_s + q_len <= q_feat->LT());
        assert(d_s >= 0 && d_len > 0 && d_s + d_len <= d_feat->LT());
        q_feat_ = q_feat;
        d_feat_ = d_feat;
        q_start_frm_ = q_s;
        d_start_frm_ = d_s;
        table_.Resize(q_len, d_len);
        ready_.Resize(q_len, d_len);
        ready_.memfill(false);
      }
    protected:
      virtual float measure_(int q, int d); /* Relative indices */
      const DenseFeature* q_feat_;
      const DenseFeature* d_feat_;
      int q_start_frm_;
      int d_start_frm_;
  };/*}}}*/



  /* Abstract class DtwRunner for common variables */
  class DtwRunner {/*{{{*/
    public:
      DtwRunner() {
        qparm_ = dparm_ = NULL;
        snippet_dist_ = NULL;
        paths_ = NULL;
        dbound_ = NULL;
      }
      virtual const PairwiseDist& PDist() const = 0;
      virtual void InitDtw(vector<float>* snippet_dist,
                           vector<IPair>* snippet_bound,
                           vector<Path>* paths,
                           const DtwParm* q_parm,
                           const DtwParm* d_parm,
                           const IPair* qbound,
                           const IPair* dbound);
      virtual void DTW() = 0;
    protected:
      virtual void CheckStartEnd(const DtwParm &parm, int *s, int *e, int *l) = 0;
      virtual void FindBestPaths(unsigned nsnippet, const float normalizer,
                                 IPair drange = IPair(0, -1));
      // FindBestPaths() helper
      virtual void _FindAllValley(vector<int>* cand, IPair drange);
      //virtual bool _IsValley(const int d);
      virtual void _TracePath(Path* p, int d);
      /* Input variables */
      const DtwParm* qparm_;
      const DtwParm* dparm_;
      const IPair* qbound_;
      const IPair* dbound_;
      /* Output variables */
      vector<float>* snippet_dist_;
      vector<IPair>* snippet_bound_;
      vector<Path>* paths_;
      /* Variables for DTW */
      int qstart_, qend_, qL_;
      int dstart_, dend_, dL_;
      TwoDimArray<IPair> lastL_;
      TwoDimArray<float> score_;
      TwoDimArray<int> root_;
      /* Tuned parameter */
      static float del_ratio_;
  };/*}}}*/

  class SegDtwRunner : public DtwRunner {/*{{{*/
    public:
      SegDtwRunner(VectorDistFn norm) : pdist(norm) {}
      const PairwiseDist& PDist() const { return pdist; }
      virtual void DTW();
      /* Tuned parameter */
      static int wq_;
      static int wd_;
      static float rpen_;
      static float lratio_;
      static unsigned nsnippet_;
      /* data */
    private:
      /* Supersegment distance */
      float SegDist(const int q, const int len_q,
                    const int d, const int len_d,
                    float upperbound = float_inf);
      /* DTW-related helpers */
      void CheckStartEnd(const DtwParm &parm, int *s, int *e, int *l);
      void MaxDelSeg();
      void CalStartDel();
      void CalScoreTable();
      void _BestScore(const int q, const int d);
      void CalEndDel();
      /* Variables for DTW */
      SegPairwiseDist pdist;
      int n_head_del, n_tail_del;
      vector<float> head_del_band;
      TwoDimArray<float> tail_del_band;
  };/*}}}*/

  class FrameDtwRunner: public DtwRunner {/*{{{*/
    public:
      FrameDtwRunner(VectorDistFn norm) : pdist(norm) {}
      const PairwiseDist& PDist() const { return pdist; }
      virtual void DTW();
      static unsigned nsnippet_;
    protected:
      virtual void CalScoreTable() = 0;
      virtual void CheckStartEnd(const DtwParm &parm, int *s, int *e, int *l);
      virtual void MaxDelFrame() { n_del = floor(del_ratio_ * qL_); }
      /* Variables for DTW */
      FramePairwiseDist pdist;
      int n_del;
  };/*}}}*/


  class SlopeConDtwRunner: public FrameDtwRunner {/*{{{*/
    public:
      SlopeConDtwRunner(VectorDistFn norm) : FrameDtwRunner(norm) {}
      /* Tuned parameter */
      static int wq_;
      static int wd_;
    protected:
      virtual void CalScoreTable();
      virtual void _BestScoreInterior(const int q, const int d);
  };/*}}}*/

  class FreeFrameDtwRunner : public FrameDtwRunner {/*{{{*/
    public:
      FreeFrameDtwRunner(VectorDistFn norm) : FrameDtwRunner(norm) {}
      static unsigned nsnippet_;
    protected:
      virtual void CalScoreTable();
  };/*}}}*/

  class FixFrameDtwRunner : public FrameDtwRunner {/*{{{*/
    public:
      FixFrameDtwRunner(VectorDistFn norm) : FrameDtwRunner(norm) {}
      /* Tuned parameter */
      static int adj_win;
    protected:
      virtual void CalScoreTable();
  };/*}}}*/

  class SegmentalDtwRunner : public FrameDtwRunner {/*{{{*/
    public:
      SegmentalDtwRunner(VectorDistFn norm) : FrameDtwRunner(norm) {}
      /* Overwrite FrameDtwRunner::DTW(), replaced with many iterations of DTWs */
      virtual void DTW();
      /* Tuned parameter */
      static unsigned nsnippet_;
      static int adj_win;
      static int step_len;
    protected:
      virtual void CalScoreTable(int d0);
      virtual void ThreewayBestScore(int q, int d);
      virtual void TwowayDallBestScore(int q, int d);
      virtual void TwowayQallBestScore(int q, int d);
    private:
      virtual void CalScoreTable() { assert(false); }
  };/*}}}*/



  /* Interface: Batch parameter sets */
  class DtwBatchParamSet {/*{{{*/
    public:
      virtual bool InstallDtwRunner(DtwRunner* runner) = 0;
  };/*}}}*/

  class QDDtwBatchParamSet : public DtwBatchParamSet {/*{{{*/
    public:
      virtual void Install(Dispatcher<UPair>* dispatcher,
                           vector<DtwParm>* query_parm_list,
                           vector<DtwParm>* doc_parm_list,
                           QDArray<vector<float> >* snippet_dist,
                           QDArray<vector<Path> >* paths,
                           vector<IPair>* vqboundary,
                           QDArray<IPair>* vdboundary) {
        dispatcher_ = dispatcher;
        query_parm_list_ = query_parm_list;
        doc_parm_list_ = doc_parm_list;
        snippet_dist_ = snippet_dist;
        paths_ = paths;
        vqboundary_ = vqboundary;
        vdboundary_ = vdboundary;

        snippet_dist_->Resize(query_parm_list_->size(), doc_parm_list_->size());
        snippet_bound_.Resize(query_parm_list_->size(), doc_parm_list_->size());
      }
      virtual bool InstallDtwRunner(DtwRunner* runner);
      QDArray<vector<IPair> > snippet_bound_;
    private:
      /* data */
      Dispatcher<UPair>* dispatcher_;
      vector<DtwParm>* query_parm_list_;
      vector<DtwParm>* doc_parm_list_;
      QDArray<vector<float> >* snippet_dist_;
      QDArray<vector<Path> >* paths_;
      vector<IPair>* vqboundary_;
      QDArray<IPair>* vdboundary_;
  };/*}}}*/

  class QSDtwBatchParamSet : public DtwBatchParamSet {/*{{{*/
    public:
      QSDtwBatchParamSet(Dispatcher<UPair>* dispatcher,
                         vector<DtwParm>* qry_parm_list,
                         vector<DtwParm>* doc_parm_list,
                         vector<IPair>* vqboundary,
                         vector<SnippetProfileList>* candidates) {
        dispatcher_ = dispatcher;
        qry_parm_list_ = qry_parm_list;
        doc_parm_list_ = doc_parm_list;
        vqboundary_ = vqboundary;
        candidates_ = candidates;

        /* Resize convenient variables */
        snippet_dist_.resize(candidates_->size());
        snippet_bound_.resize(candidates_->size());

        for (unsigned qidx = 0; qidx < candidates_->size(); ++qidx) {
          snippet_dist_[qidx].resize((*candidates_)[qidx].size());
          snippet_bound_[qidx].resize((*candidates_)[qidx].size());
        }

      }
      virtual bool InstallDtwRunner(DtwRunner* runner);
      virtual void CollectResult(bool do_sort = false);
    private:
      /* input */
      Dispatcher<UPair>* dispatcher_;
      vector<DtwParm>* qry_parm_list_;
      vector<DtwParm>* doc_parm_list_;
      vector<IPair>* vqboundary_;
      /* input & output */
      vector<SnippetProfileList>* candidates_;
      /* storage */
      vector<vector<vector<float> > > snippet_dist_;
      vector<vector<vector<IPair> > > snippet_bound_;
  };/*}}}*/

  struct UTriple {/*{{{*/
    UTriple (unsigned q, unsigned s, unsigned p) : qidx(q), sidx(s), pidx(p) {}
    unsigned qidx;
    unsigned sidx;
    unsigned pidx;
  };/*}}}*/

  class PrfDtwBatchParamSet : public DtwBatchParamSet {/*{{{*/
    public:
      PrfDtwBatchParamSet(Dispatcher<UTriple>* dispatcher,/*{{{*/
                          const vector<DtwParm>* doc_parm_list,
                          vector<SnippetProfileList>* candidates,
                          const vector<unsigned>* num_prf) {
        dispatcher_ = dispatcher;
        doc_parm_list_ = doc_parm_list;
        candidates_ = candidates;
        num_prf_ = num_prf;

        /* Resize output variables */
        output_dist_.resize(candidates_->size()); // Number of queries

        for (unsigned qidx = 0; qidx < candidates_->size(); ++qidx) {
          // Number of hypothesized regions
          output_dist_[qidx].resize((*candidates_)[qidx].size());

          for (unsigned sidx = 0; sidx < output_dist_[qidx].size(); ++sidx) {
            // Number of pseudo relevant regions
            output_dist_[qidx][sidx].resize((*num_prf_)[qidx]);
          }
        }

      }/*}}}*/
      virtual bool InstallDtwRunner(DtwRunner* runner);
      void IntegratePrfDist(const vector<float>& prf_weight);
    private:
      /* data */
      Dispatcher<UTriple>* dispatcher_;
      const vector<DtwParm>* doc_parm_list_;
      vector<SnippetProfileList>* candidates_;
      const vector<unsigned>* num_prf_;
      /* local */
      vector<vector<vector<vector<float> > > > output_dist_;
  };/*}}}*/

  void InitPrfDispatcher(Dispatcher<UTriple>* prf_dispatcher,
                         const vector<SnippetProfileList>& candidates,
                         const vector<unsigned>& num_prf);



  /* Glue DtwRunner and DtwBatchParamSet together */
  class DtwManager : public ThreadRunner {/*{{{*/
    public:
      DtwManager() {
        paramset_ = NULL;
        runner_ = NULL;
      }
      virtual void* Run();
      void Install(DtwBatchParamSet* paramset, DtwRunner* runner) {
        paramset_ = paramset;
        runner_ = runner;
      }

    private:
      DtwBatchParamSet* paramset_;
      DtwRunner* runner_;
  };/*}}}*/

} /* namespace DtwUtil */

#endif
