# libdtw - Dynamic time warping library

## General information
- This library contains DTW related code.
- Two files: **dtw_parm.h/cpp**  <-depend on-- **dtw_util.h/cpp**.
- Support **multi-threading** for **batch processing** (multiple queries searching for multiple documents):
  - Create N threads.
  - Dispatcher holds all (query, doc) pairs.
  - Each thread requests a (query, doc) pair in a while loop.

## dtw_parm.h/cpp
- ``class DtwParm`` Contains parameters to run DTW.

### functions (partial)
- ``DtwParm::LoadParm`` is overloaded with three prototype.
  - All three types can be used in frame-based DTW.
  - The 2nd type has additional basic segment infomation (but not supersegment informations)
  - Only 3rd type can be used in segment-based DTW.
  See [section 3.3](https://drive.google.com/file/d/1TqvicnaZABsRSHIh3PjJWVRic495qcEP/view?usp=sharing)

```cpp
namespace DtwUtil {
  class DtwParm {
    public:
      enum ParmType {FRAME = 0, SEGMENT = 1, SUPERSEG = 2};
      DtwParm();

      /* Load data */
      void LoadParm(const string feat_fname); // 1st type: Load feature only
      void LoadParm(const string feat_fname,  // 2nd type: Load feature and segment
                    const double bseg_ratio,
                    const string tree_fname = "");
      void LoadParm(const string feat_fname,  // 3rd type: Load feature,segment,supersegment
                    const double bseg_ratio,
                    const double superseg_ratio,
                    const unsigned width,
                    const unsigned gran,
                    const string tree_fname = "");

      void DumpData() const;
  };
} /* namespace DtwUtil */

```

Examples
```cpp
#include "dtw_parm.h"

using DtwUtil::DtwParm;

int main() {
  string filename = "data/iv1_1.fbank";
  DtwParm parm;

  /* Load Parm for frame-based DTW */
  parm.LoadParm(filename);
  parm.DumpData();

  /* Load Parm for segment-based DTW but only generate basic segments.
   * An empty `tree_fname' invokes SegTree::ConstructTree(DenseFeature).
   */
  float bseg_ratio = 0.5;
  parm.LoadParm(filename, bseg_ratio, "");
  parm.DumpData();

  /* Load Parm for segment-based DTW with supersements loaded.
   * Each supersegment consist of no more then `width' basic segments.
   * When computing distance of supersegments, it is decomposed into `gran'
   * parts.
   * NOTE: width <= gran
   */
  float superseg_ratio = 4.0;
  int width = 3;
  int gran = 3;
  parm.LoadParm(filename, bseg_ratio, superseg_ratio, width, gran, "");
  parm.DumpData();

  return 0;
}
```
**NOTE**: DtwParm can be reuse. "parm" is used three times in the above code.
Whenever you need to reload it, call LoadParm().

## dtw_util.h/cpp
![DTW classes dependencies](img/dtw_class_diagram.jpg "DTW classes dependencies")

### Functions
* InitDtw(): setup input/output for DTW run
* DTW(): perform DTW on a query-document pair
```cpp
namespace DtwUtil {

  /* Frame-wise distance function */
  typedef float (*VectorDistFn)(const float*, const float*, const int);
  float innernorm(const float* a, const float* b, const int size);
  float euclinorm(const float* a, const float* b, const int size);

  class DtwRunner {
    public:
      /* ctor */
      DtwRunner();
      /* return reference to PairwiseDist in it, which has 
       * information of pairwise distances, number of
       * calculations performed, etc.
       */
      virtual const PairwiseDist& PDist() const = 0;
      /* Initialize DTW input/output
       * DTW will search for multiple hypothesized regions (snippets)
       * within document. The maximum number of snippets returned
       * can be specified by FrameDtwRunner::nsnippet_ or
       * SegDtwRunner::nsnippet_.
       */
      virtual void InitDtw(vector<float>* snippet_dist,
                           vector<IPair>* snippet_bound,
                           vector<Path>* paths,
                           const DtwParm* q_parm,
                           const DtwParm* d_parm,
                           const IPair* qbound,
                           const IPair* dbound);
      /* Run DTW (based on derived classes' implementation) */
      virtual void DTW() = 0;
  };

  /* !!!For concrete DtwRunner constructors, VectorDistFn must be specified.!!! */
  class SegDtwRunner {
    SegDtwRunner(VectorDistFn norm);
  };

  class SlopeConDtwRunner {
    SlopeConDtwRunner(VectorDistFn norm);
  };
  
  class SegmentalDtwRunner {
    SegmentalDtwRunner(VectorDistFn norm);
  };
} /* namespace DtwUtil */
```

### Example
Running DTW is simple:
```cpp
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
```
Output
```
-- Slope-Constrained DTW "data/iv1_1.fbank"
hypothesized region[0]: score = -2.69686, time span = (10, 71)
hypothesized region[1]: score = -3.54951, time span = (174, 245)
hypothesized region[2]: score = -4.04659, time span = (47, 120)
hypothesized region[3]: score = -4.08004, time span = (255, 310)
hypothesized region[4]: score = -4.09251, time span = (74, 155)
-- Slope-Constrained DTW "data/iv2_1.fbank"
hypothesized region[0]: score = -2.69686, time span = (10, 71)
hypothesized region[1]: score = -3.54951, time span = (174, 245)
hypothesized region[2]: score = -3.85096, time span = (255, 310)
hypothesized region[3]: score = -3.83785, time span = (11, 76)
hypothesized region[4]: score = -3.76874, time span = (314, 380)
-- Segment-based DTW "data/iv1_1.fbank"
hypothesized region[0]: score = -2.31454, time span = (0, 5) in frame = (0, 79)
hypothesized region[1]: score = -4.26458, time span = (11, 16) in frame = (157, 259)
hypothesized region[2]: score = -4.36012, time span = (3, 13) in frame = (40, 209)
hypothesized region[3]: score = -4.36789, time span = (15, 22) in frame = (219, 349)
hypothesized region[4]: score = -4.71572, time span = (25, 29) in frame = (381, 499)
-- Segment-based DTW "data/iv2_1.fbank"
hypothesized region[0]: score = -2.31454, time span = (0, 5) in frame = (0, 79)
hypothesized region[1]: score = -4.0183, time span = (20, 24) in frame = (304, 380)
hypothesized region[2]: score = -3.79656, time span = (13, 20) in frame = (202, 319)
hypothesized region[3]: score = -3.62125, time span = (0, 8) in frame = (0, 129)
hypothesized region[4]: score = -3.54052, time span = (25, 29) in frame = (381, 499)
```

## Batch: multithreading
Assuming you have a query, and you want to search for hypothesized regions in all documents. The DTWs performed on different (query, document) pairs are disjoint. Therefore we divide them into smaller tasks, namely a (query, document) pairs. All threads can perform DTW on each task simultaneously.

### DtwManager
```cpp
namespace DtwUtil {
  class DtwManager : public ThreadRunner {
    public:
      DtwManager();
      void Install(DtwBatchParamSet* paramset, DtwRunner* runner) {
      virtual void* Run() {
        while (paramset_->InstallDtwRunner(runner_)) runner_->DTW();
        return NULL;
      }
    private:
      DtwBatchParamSet* paramset_; // Make sure all DtwMangers share the same DtwBatchParamSet.
      DtwRunner* runner_; // Make sure all DtwManagers have *distinct* DtwRunners, otherwise it will fail.
  };
} /* namespace DtwUtil */
```
