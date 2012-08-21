#ifndef __DEBUG_H__
#define __DEBUG_H__

inline void DUMP(const char *str) {
#ifdef DEBUG
  cout << str << endl;
#endif
}

#endif
