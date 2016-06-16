#ifndef MATH_STRUCT_H
#define MATH_STRUCT_H

typedef double math_float;

// 2d vector algebra
typedef struct VEC2D {
  math_float x, y, mod;
} VEC2D;

typedef struct VEC3D {
  math_float x, y, z, mod;
} VEC3D;

typedef struct MAT2D {
  math_float xx, xy, yx, yy;
} MAT2D;

typedef struct MAT3D {
  math_float xx, xy, xz, yx, yy, yz, zx, zy, zz;
} MAT3D;

typedef struct VEC6D {
  VEC3D a, g;
} VEC6D;

typedef struct MAT6D {
  MAT3D A[2][2];
} MAT6D;

typedef struct EigenSys {
  MAT2D A;
  math_float l1, l2;
  VEC2D v1, v2, u1, u2;
} EigenSys;

#endif
