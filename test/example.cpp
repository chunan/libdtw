#include "dtw_parm.h"

using DtwUtil::DtwParm;

int main() {

  string filename = "data/iv1_1.fbank";
  DtwParm parm;
  float bseg_ratio = 0.5;
  float superseg_ratio = 4.0;

  /* Load Parm for frame-based DTW */
  parm.LoadParm(filename);
  parm.DumpData();

  /* Load Parm for segment-based DTW but only generate basic segments.
   * An empty `tree_fname' invokes SegTree::ConstructTree(DenseFeature).
   */
  parm.LoadParm(filename, bseg_ratio, "");
  parm.DumpData();

  /* Load Parm for segment-based DTW with supersements loaded.
   * An empty `tree_fname' invokes SegTree::ConstructTree(DenseFeature).
   */
  parm.LoadParm(filename, bseg_ratio, superseg_ratio, 3, 3, "");
  parm.DumpData();

  return 0;
}
