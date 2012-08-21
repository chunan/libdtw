#include <iostream>
#include "dtw_parm.h"
#include "dtw_util.h"

using std::cout;
using DtwUtil::DtwParm;
using DtwUtil::SegDtwRunner;
using DtwUtil::FrameDtwRunner;
using DtwUtil::SlopeConDtwRunner;

int main() {

  string q_fname = "data/iv1_1.fbank";
  string d_fname = "data/N200108011200-09-14.fbank";
  DtwParm q_parm, d_parm;

  /* Frame-based DTW */
  q_parm.LoadParm(q_fname);
  d_parm.LoadParm(d_fname);
  SlopeConDtwRunner scdtw_runner(DtwUtil::euclinorm);
  FrameDtwRunner::nsnippet_ = 10;
  vector<float> hypo_score;
  vector<pair<int, int> > hypo_bound;
  scdtw_runner.InitDtw(&hypo_score,
                       &hypo_bound, /* (start, end) frame */
                       NULL, /* do not backtracking */
                       &q_parm,
                       &d_parm,
                       NULL, /* full time span */
                       NULL); /* full time span */
  scdtw_runner.DTW();

  unsigned num_hypo = hypo_score.size();
  cout << "-- Slope-Constrained DTW --\n";
  for (unsigned i = 0; i < num_hypo; ++i) {
    cout << "hypothesized region[" << i << "]: score = " << hypo_score[i]
      << ", time span = (" << hypo_bound[i].first
      << ", " << hypo_bound[i].second << ")\n";
  }


  /* Segment-based DTW */
  hypo_score.clear();
  hypo_bound.clear();
  float bseg_ratio = 0.5;
  float superseg_ratio = 4.0;
  int gran = 3, width = 3;
  q_parm.LoadParm(q_fname, bseg_ratio, superseg_ratio, width, gran, "");
  d_parm.LoadParm(d_fname, bseg_ratio, superseg_ratio, width, gran, "");
  SegDtwRunner segdtw_runner(DtwUtil::euclinorm);
  SegDtwRunner::nsnippet_ = 10;
  segdtw_runner.InitDtw(&hypo_score,
                       &hypo_bound, /* (start, end) basic segment */
                       NULL, /* do not backtracking */
                       &q_parm,
                       &d_parm,
                       NULL, /* full time span */
                       NULL); /* full time span */
  segdtw_runner.DTW();
  num_hypo = hypo_score.size();
  cout << "-- Segment-based DTW --\n";
  for (unsigned i = 0; i < num_hypo; ++i) {
    cout << "hypothesized region[" << i << "]: score = " << hypo_score[i]
      << ", time span = (" << hypo_bound[i].first
      << ", " << hypo_bound[i].second << ")\n";
  }


  return 0;
}
