#include <iostream>
#include "dtw_parm.h"
#include "dtw_util.h"

using std::cout;
using DtwUtil::DtwParm;
using DtwUtil::SegDtwRunner;
using DtwUtil::FrameDtwRunner;
using DtwUtil::SlopeConDtwRunner;
using DtwUtil::DeterminePhiFn;

int main() {

  vector<string> q_fname;
  q_fname.push_back("data/iv1_1.fbank");
  q_fname.push_back("data/iv2_1.fbank");
  string d_fname = "data/N200108011200-09-14.fbank";
  DtwParm q_parm, d_parm;
  d_parm.LoadParm(d_fname);

  /* Use DeterminePhiFn() for local distance function selection */
  SlopeConDtwRunner scdtw_runner(DeterminePhiFn(GetExtension(q_fname[0])));
  // or you choose it your own like:
  // SlopeConDtwRunner scdtw_runner(DtwUtil::euclinorm);
  FrameDtwRunner::nsnippet_ = 5;
  vector<float> hypo_score;
  vector<pair<int, int> > hypo_bound;

  /* Frame-based DTW */
  //   DtwRunners can be reuse (SegDtwRunner too)
  //   Everytime you call InitDtw() then DTW(), it's like brand new!
  for (size_t i = 0; i < q_fname.size(); ++i) {
    q_parm.LoadParm(q_fname[i]);
    scdtw_runner.InitDtw(&hypo_score,
                         &hypo_bound, /* (start, end) frame */
                         NULL, /* do not backtracking */
                         &q_parm,
                         &d_parm,
                         NULL, /* full time span */
                         NULL); /* full time span */
    scdtw_runner.DTW();

    unsigned num_hypo = hypo_score.size();
    cout << "-- Slope-Constrained DTW \"" << q_fname[i] << "\"\n";
    for (unsigned i = 0; i < num_hypo; ++i) {
      cout << "hypothesized region[" << i << "]: score = " << hypo_score[i]
        << ", time span = (" << hypo_bound[i].first
        << ", " << hypo_bound[i].second << ")\n";
    }
  }


  /* Segment-based DTW */
  hypo_score.clear();
  hypo_bound.clear();
  float bseg_ratio = 0.5;
  float superseg_ratio = 4.0;
  int gran = 3, width = 3;
  d_parm.LoadParm(d_fname, bseg_ratio, superseg_ratio, width, gran, "");
  SegDtwRunner segdtw_runner(DtwUtil::euclinorm);
  SegDtwRunner::nsnippet_ = 5;

  for (size_t i = 0; i < q_fname.size(); ++i) {
    q_parm.LoadParm(q_fname[i], bseg_ratio, superseg_ratio, width, gran, "");
    //q_parm.DumpData();
    segdtw_runner.InitDtw(&hypo_score,
                          &hypo_bound, /* (start, end) basic segment */
                          NULL, /* do not backtracking */
                          &q_parm,
                          &d_parm,
                          NULL, /* full time span */
                          NULL); /* full time span */
    segdtw_runner.DTW();
    size_t num_hypo = hypo_score.size();
    cout << "-- Segment-based DTW \"" << q_fname[i] << "\"\n";
    for (unsigned i = 0; i < num_hypo; ++i) {
      cout << "hypothesized region[" << i << "]: score = " << hypo_score[i]
        << ", time span = (" << hypo_bound[i].first
        << ", " << hypo_bound[i].second << ")"
        << " in frame = (" << d_parm.BasicSegStartT(hypo_bound[i].first) << ", "
        << d_parm.BasicSegEndT(hypo_bound[i].second) << ")\n";

    }
  }


  return 0;
}
