#ifndef MATH_LIB_H
#define MATH_LIB_H

#include <stdio.h>
#include "math_struct.h"

#define EPSILON          1e-8

#ifdef __cplusplus
extern "C" {
#endif

  /* 2d vector algebra */
  VEC2D set_vec2d(math_float, math_float);
  math_float mod_vec2d(VEC2D *);
  void normalize_vec2d(VEC2D *);
  math_float prod_dot_2d(VEC2D, VEC2D);
  math_float prod_cross_2d(VEC2D, VEC2D);            /* returns z-component on a right-handed frame */
  void multiply_vec2d(math_float, VEC2D *);


  /* 3d vector algebra */
  VEC3D set_vec3d(math_float, math_float, math_float);
  math_float mod_vec3d(VEC3D *);
  void normalize_vec3d(VEC3D *);
  VEC3D _normalize_vec3d(VEC3D v);
  math_float prod_dot(VEC3D, VEC3D);
  VEC3D prod_cross(VEC3D, VEC3D);
  void multiply_vec3d(math_float, VEC3D *);


  /* 2d matrix algebra */
  MAT2D set_mat2d(math_float, math_float, math_float, math_float);
  MAT2D transpose_mat2d(MAT2D);
  MAT2D product_mat2d(MAT2D, MAT2D);
  void multiply_mat2d(math_float, MAT2D *);
  MAT2D make_rotation_2d_cs(math_float, math_float);
  MAT2D make_rotation_2d(math_float);
  VEC2D rotate_vec2d(MAT2D, VEC2D);
  MAT2D rotate_mat2d(MAT2D, MAT2D);


  /* 3d matrix algebra */
  MAT3D set_mat3d(math_float, math_float, math_float, math_float, math_float, math_float, math_float, math_float, math_float);
  MAT3D transpose_mat3d(MAT3D);
  MAT3D product_mat3d(MAT3D, MAT3D);
  void multiply_mat3d(math_float, MAT3D *);
  MAT3D make_rotation_cs(VEC3D, math_float, math_float);
  MAT3D make_rotation(VEC3D, math_float);
  VEC3D rotate_vec3d(MAT3D, VEC3D);
  MAT3D rotate_mat3d(MAT3D, MAT3D);


  /* 6d vector */
  VEC6D set_vec6d(math_float, math_float, math_float, math_float, math_float, math_float);
  VEC6D rotate_vec6d(MAT3D, VEC6D);


  /* 6d matrix */
  MAT6D set_mat6d(MAT3D, MAT3D, MAT3D, MAT3D);
  MAT6D rotate_mat6d(MAT3D, MAT6D);


  /* 2d eigenvalue problem */
  EigenSys eigs_2x2_sym_normalized(math_float, math_float);
  EigenSys eigs_2x2_sym(MAT2D);
  int check_eigs(EigenSys);

#ifdef _MSC_VER
  /* display function */
  void print_vec2d(VEC2D, const char *);
  void print_vec3d(VEC3D, const char *);
  void print_mat2d(MAT2D, const char *);
  void print_mat3d(MAT3D, const char *);
  void print_vec6d(VEC6D, const char *);
  void print_mat6d(MAT6D, const char *);
  void print_eigs(EigenSys, const char *);
#endif

#ifdef __cplusplus
}
#endif

#endif


