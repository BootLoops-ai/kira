
#include "pyred/config.h"
// Needed to set the number PYRED_PP_NCOEFFCLASSES of coefficient classes
// to instantiate. The data types must be defined below.
#ifndef PYRED_PP_NCOEFFCLASSES
#define PYRED_PP_NCOEFFCLASSES 1
#endif

// Define the data types here
// (types beyond PYRED_PP_COEFFCLASS_<PYRED_PP_NCOEFFCLASSES-1>
//  will be ignored).
#define PYRED_PP_COEFFCLASS_0 Coeff_int
#define PYRED_PP_COEFFCLASS_1 Coeff_vec
#define PYRED_PP_COEFFCLASS_2 Coeff_arr<2>
#define PYRED_PP_COEFFCLASS_3 Coeff_arr<4>
#define PYRED_PP_COEFFCLASS_4 Coeff_arr<8>
#define PYRED_PP_COEFFCLASS_5 Coeff_arr<16>
#define PYRED_PP_COEFFCLASS_6 Coeff_arr<32>
#define PYRED_PP_COEFFCLASS_7 Coeff_arr<64>
#define PYRED_PP_COEFFCLASS_8 Coeff_arr<128>
#define PYRED_PP_COEFFCLASS_9 Coeff_arr<256>

// Use PYRED_PP_COEFFCLASSMEM(n) as names of members of a struct or union
// that contains one member per coefficient class.
#define PYRED_PP_COEFFCLASSMEM(n) cls_##n

// Nothing to edit beyond this line.

#define PYRED_PP_COEFFCLASS(n) PYRED_PP_COEFFCLASS_##n
#define PYRED_PP_NUMSYS(n) PYRED_PP_NUMSYS_##n
#define PYRED_PP_SOLMAP(n) PYRED_PP_SOLMAP_##n
#define PYRED_PP_COEFFMAP(n) PYRED_PP_COEFFMAP_##n

#define PYRED_PP_BOOL(n) PYRED_PP_BOOL_##n
#define PYRED_PP_BOOL_0 0
#define PYRED_PP_BOOL_1 1
#define PYRED_PP_BOOL_2 1
#define PYRED_PP_BOOL_3 1
#define PYRED_PP_BOOL_4 1
#define PYRED_PP_BOOL_5 1
#define PYRED_PP_BOOL_6 1
#define PYRED_PP_BOOL_7 1
#define PYRED_PP_BOOL_8 1
#define PYRED_PP_BOOL_9 1

// Repetition without further parameter
//   PYRED_PP_REPEAT0(n,m) --> m(0) m(1) m(2) ... m(n-1)
#define PYRED_PP_REPEAT0(n, m) PYRED_PP_REPEAT0_N(n, m)
// The step with PYRED_PP_REPEAT0_N is essential to expand n if it is a macro.
#define PYRED_PP_REPEAT0_N(n, m) PYRED_PP_REPEAT0_##n(m)
#define PYRED_PP_REPEAT0_0(m)
#define PYRED_PP_REPEAT0_1(m) m(0)
#define PYRED_PP_REPEAT0_2(m) PYRED_PP_REPEAT0_1(m) m(1)
#define PYRED_PP_REPEAT0_3(m) PYRED_PP_REPEAT0_2(m) m(2)
#define PYRED_PP_REPEAT0_4(m) PYRED_PP_REPEAT0_3(m) m(3)
#define PYRED_PP_REPEAT0_5(m) PYRED_PP_REPEAT0_4(m) m(4)
#define PYRED_PP_REPEAT0_6(m) PYRED_PP_REPEAT0_5(m) m(5)
#define PYRED_PP_REPEAT0_7(m) PYRED_PP_REPEAT0_6(m) m(6)
#define PYRED_PP_REPEAT0_8(m) PYRED_PP_REPEAT0_7(m) m(7)
#define PYRED_PP_REPEAT0_9(m) PYRED_PP_REPEAT0_8(m) m(8)

// The same as PYRED_PP_REPEAT0, but repeat n-2 times.
// This is to iterate over all instantiations of Coeff_arr.
//   PYRED_PP_REPEAT0(n,m) --> m(2) ... m(n-1)
#define PYRED_PP_REPEATCARR(n, m) PYRED_PP_REPEATCARR_N(n, m)
#define PYRED_PP_REPEATCARR_N(n, m) PYRED_PP_REPEATCARR_##n(m)
#define PYRED_PP_REPEATCARR_0(m)
#define PYRED_PP_REPEATCARR_1(m)
#define PYRED_PP_REPEATCARR_2(m)
#define PYRED_PP_REPEATCARR_3(m) m(2)
#define PYRED_PP_REPEATCARR_4(m) PYRED_PP_REPEATCARR_3(m) m(3)
#define PYRED_PP_REPEATCARR_5(m) PYRED_PP_REPEATCARR_4(m) m(4)
#define PYRED_PP_REPEATCARR_6(m) PYRED_PP_REPEATCARR_5(m) m(5)
#define PYRED_PP_REPEATCARR_7(m) PYRED_PP_REPEATCARR_6(m) m(6)
#define PYRED_PP_REPEATCARR_8(m) PYRED_PP_REPEATCARR_7(m) m(7)
#define PYRED_PP_REPEATCARR_9(m) PYRED_PP_REPEATCARR_8(m) m(8)

// If conditions
//   PYRED_PP_IF(n,t,f): for integer n>=0 return f if n==0, t otherwise.
//   PYRED_PP_IF_TRUE(n,t): if n>=0 return t, if n==0 return nothing.
//   PYRED_PP_IF_FALSE(n,f): if n==0 return f, if n>=0 return nothing.

#define PYRED_PP_IF_TRUE(n, t) PYRED_PP_IF(n, t, )
#define PYRED_PP_IF_FALSE(n, f) PYRED_PP_IF(n, , f)
#define PYRED_PP_IF(n, t, f) PYRED_PP_IF_01A(PYRED_PP_BOOL(n), t, f)
#define PYRED_PP_IF_01A(bl, t, f) PYRED_PP_IF_01B(bl, t, f)
// The step with PYRED_PP_IF_01B is essential to expand bl if it is a macro.
#define PYRED_PP_IF_01B(bl, t, f) PYRED_PP_IF_##bl(t, f)
#define PYRED_PP_IF_0(t, f) f
#define PYRED_PP_IF_1(t, f) t
