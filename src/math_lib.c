
#include "math_lib.h"
#include <math.h>

// 2d vector algebra
VEC2D set_vec2d(math_float x, math_float y) {
  VEC2D v;
  v.x = x;
  v.y = y;
  return v;
}

math_float mod_vec2d(VEC2D * v) {
  return sqrt(v->x*v->x + v->y*v->y);
}

void normalize_vec2d(VEC2D * v) {
  math_float mod = mod_vec2d(v);
  v->x /= mod;
  v->y /= mod;
}

math_float prod_dot_2d(VEC2D a, VEC2D b) {
  return a.x * b.x + a.y * b.y;
}

math_float prod_cross_2d(VEC2D a, VEC2D b) {
  return a.x * b.y - a.y * b.x;
}

void multiply_vec2d(math_float alpha, VEC2D * v) {
  v->x *= alpha;
  v->y *= alpha;
}


// 3d vector algebra
VEC3D set_vec3d(math_float x, math_float y, math_float z) {
  VEC3D v;
  v.x = x;
  v.y = y;
  v.z = z;
  return v;
}

math_float mod_vec3d(VEC3D * v) {
  return sqrt(v->x*v->x + v->y*v->y + v->z*v->z);
}

void normalize_vec3d(VEC3D * v) {
  math_float mod = mod_vec3d(v);
  v->x /= mod;
  v->y /= mod;
  v->z /= mod;
}

VEC3D _normalize_vec3d(VEC3D v) {
  VEC3D n;
  math_float mod = mod_vec3d(&v);
  n = set_vec3d(v.x / mod, v.y / mod, v.z / mod);
  return n;
}

math_float prod_dot(VEC3D a, VEC3D b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

VEC3D prod_cross(VEC3D a, VEC3D b) {
  return set_vec3d(a.y * b.z - a.z * b.y,
    a.z * b.x - a.x * b.z,
    a.x * b.y - a.y * b.x);
}

void multiply_vec3d(math_float alpha, VEC3D * v) {
  v->x *= alpha;
  v->y *= alpha;
  v->z *= alpha;
}


// 2d matrix algebra
MAT2D set_mat2d(math_float xx, math_float xy, math_float yx, math_float yy) {
  MAT2D A;
  A.xx = xx;
  A.xy = xy;
  A.yx = yx;
  A.yy = yy;
  return A;
}

MAT2D transpose_mat2d(MAT2D m) {
  MAT2D mt;
  mt.xx = m.xx;
  mt.xy = m.yx;
  mt.yx = m.xy;
  mt.yy = m.yy;
  return mt;
}

MAT2D product_mat2d(MAT2D A, MAT2D B) {
  MAT2D C;
  C.xx = A.xx * B.xx + A.xy * B.yx;
  C.xy = A.xx * B.xy + A.xy * B.yy;
  C.yx = A.yx * B.xx + A.yy * B.yx;
  C.yy = A.yx * B.xy + A.yy * B.yy;
  return C;
}

void multiply_mat2d(math_float alpha, MAT2D * A) {
  A->xx *= alpha;
  A->xy *= alpha;
  A->yx *= alpha;
  A->yy *= alpha;
}

MAT2D make_rotation_2d_cs(math_float c, math_float s) {
  MAT2D rot;
  rot.xx = c;
  rot.xy = -s;
  rot.yx = s;
  rot.yy = c;
  return rot;
}

MAT2D make_rotation_2d(math_float theta) {
  math_float c = cos(theta), s = sin(theta);
  MAT2D rot = make_rotation_2d_cs(c, s);
  return rot;
}

VEC2D rotate_vec2d(MAT2D R, VEC2D v) {
  VEC2D vr;
  vr.x = R.xx * v.x + R.xy * v.y;
  vr.y = R.yx * v.x + R.yy * v.y;
  return vr;
}

MAT2D rotate_mat2d(MAT2D R, MAT2D A) {
  MAT2D C;
  C = product_mat2d(product_mat2d(R, A), transpose_mat2d(R));
  return C;
}


// 3d matrix algebra
MAT3D set_mat3d(math_float xx, math_float xy, math_float xz, math_float yx, math_float yy, math_float yz, math_float zx, math_float zy, math_float zz) {
  MAT3D m;
  m.xx = xx;
  m.xy = xy;
  m.xz = xz;
  m.yx = yx;
  m.yy = yy;
  m.yz = yz;
  m.zx = zx;
  m.zy = zy;
  m.zz = zz;
  return m;
}

MAT3D transpose_mat3d(MAT3D A) {
  MAT3D C;
  C.xx = A.xx;
  C.xy = A.yx;
  C.xz = A.zx;
  C.yx = A.xy;
  C.yy = A.yy;
  C.yz = A.zy;
  C.zx = A.xz;
  C.zy = A.yz;
  C.zz = A.zz;
  return C;
}

MAT3D product_mat3d(MAT3D A, MAT3D B) {
  MAT3D C;
  C.xx = A.xx * B.xx + A.xy * B.yx + A.xz * B.zx;
  C.xy = A.xx * B.xy + A.xy * B.yy + A.xz * B.zy;
  C.xz = A.xx * B.xz + A.xy * B.yz + A.xz * B.zz;
  C.yx = A.yx * B.xx + A.yy * B.yx + A.yz * B.zx;
  C.yy = A.yx * B.xy + A.yy * B.yy + A.yz * B.zy;
  C.yz = A.yx * B.xz + A.yy * B.yz + A.yz * B.zz;
  C.zx = A.zx * B.xx + A.zy * B.yx + A.zz * B.zx;
  C.zy = A.zx * B.xy + A.zy * B.yy + A.zz * B.zy;
  C.zz = A.zx * B.xz + A.zy * B.yz + A.zz * B.zz;
  return C;
}

void multiply_mat3d(math_float alpha, MAT3D * m) {
  m->xx *= alpha;
  m->xy *= alpha;
  m->xz *= alpha;
  m->yx *= alpha;
  m->yy *= alpha;
  m->yz *= alpha;
  m->zx *= alpha;
  m->zy *= alpha;
  m->zz *= alpha;
}

MAT3D make_rotation_cs(VEC3D axis, math_float c, math_float s) {
  axis = _normalize_vec3d(axis);

  MAT3D rot;
  rot.xx = axis.x*axis.x*(1 - c) + c;   	    rot.xy = axis.x*axis.y*(1 - c) - axis.z*s;  rot.xz = axis.x*axis.z*(1 - c) + axis.y*s;
  rot.yx = axis.y*axis.x*(1 - c) + axis.z*s;	rot.yy = axis.y*axis.y*(1 - c) + c;	        rot.yz = axis.y*axis.z*(1 - c) - axis.x*s;
  rot.zx = axis.z*axis.x*(1 - c) - axis.y*s;	rot.zy = axis.z*axis.y*(1 - c) + axis.x*s;  rot.zz = axis.z*axis.z*(1 - c) + c;
  return rot;
}

MAT3D make_rotation(VEC3D axis, math_float angle) {
  math_float c, s;
  c = cos(angle);
  s = sin(angle);

  MAT3D rot = make_rotation_cs(axis, c, s);
  return rot;
}

VEC3D rotate_vec3d(MAT3D rot, VEC3D v) {
  VEC3D vr;
  vr.x = rot.xx * v.x + rot.xy * v.y + rot.xz * v.z;
  vr.y = rot.yx * v.x + rot.yy * v.y + rot.yz * v.z;
  vr.z = rot.zx * v.x + rot.zy * v.y + rot.zz * v.z;
  return vr;
}

MAT3D rotate_mat3d(MAT3D rot, MAT3D m) {
  MAT3D mr;
  mr = product_mat3d(rot, product_mat3d(m, transpose_mat3d(rot)));
  return mr;
}


// 6d vector
VEC6D set_vec6d(math_float ax, math_float ay, math_float az, math_float gx, math_float gy, math_float gz) {
  VEC6D v;
  v.a = set_vec3d(ax, ay, az);
  v.g = set_vec3d(gx, gy, gz);
  return v;
}

VEC6D rotate_vec6d(MAT3D rot, VEC6D v) {
  VEC6D vr;
  vr.a = rotate_vec3d(rot, v.a);
  vr.g = rotate_vec3d(rot, v.g);
  return vr;
}


// 6d matrix
MAT6D set_mat6d(MAT3D A00, MAT3D A01, MAT3D A10, MAT3D A11) {
  MAT6D m;
  m.A[0][0] = A00;
  m.A[0][1] = A01;
  m.A[1][0] = A10;
  m.A[1][1] = A11;
  return m;
}

MAT6D rotate_mat6d(MAT3D rot, MAT6D m) {
  MAT6D mr;
  mr.A[0][0] = rotate_mat3d(rot, m.A[0][0]);
  mr.A[0][1] = rotate_mat3d(rot, m.A[0][1]);
  mr.A[1][0] = rotate_mat3d(rot, m.A[1][0]);
  mr.A[1][1] = rotate_mat3d(rot, m.A[1][1]);
  return mr;
}


// 2d eigenvalue problem

// Returns the eigensystem of a 2x2 matrix of the form
// | 1 a |
// | a b |
// for real a,b
EigenSys eigs_2x2_sym_normalized(math_float a, math_float b) {
  EigenSys es;
  es.A.xx = 1.0;
  es.A.xy = a;
  es.A.yx = a;
  es.A.yy = b;
  es.l1 = (1 + b)*0.5 + 0.5*sqrt((b - 1)*(b - 1) + 4 * a*a);
  es.l2 = (1 + b)*0.5 - 0.5*sqrt((b - 1)*(b - 1) + 4 * a*a);
  es.v1 = set_vec2d(1.0, (es.l1 - 1) / a);
  es.v2 = set_vec2d(1.0, (es.l2 - 1) / a);
  if (prod_cross_2d(es.v1, es.v2) < 0) { // "right-handed", in the sense that z-axis is pointing "outward" the screen
    multiply_vec2d(-1.0, &es.v2);
  }
  es.u1 = es.v1; normalize_vec2d(&es.u1);
  es.u2 = es.v2; normalize_vec2d(&es.u2);
  return es;
}

// Returns the eigensystem of a 2x2 matrix of the form
// | a b |
// | b c |
// for real a,b,c
EigenSys eigs_2x2_sym(MAT2D A) {
  math_float trA, detA;
  EigenSys es;
  es.A = set_mat2d(A.xx, A.xy, A.yx, A.yy);
  trA = es.A.xx + es.A.yy;
  detA = es.A.xx * es.A.yy - es.A.xy * es.A.yx;
  es.l1 = (trA + sqrt(trA*trA - 4 * detA))*0.5;
  es.l2 = (trA - sqrt(trA*trA - 4 * detA))*0.5;
  es.v1 = set_vec2d(1.0, (es.l1 - es.A.xx) / es.A.xy);         // A.xx * 1 + A.xy * y = l1 * 1
  es.v2 = set_vec2d(1.0, (es.l2 - es.A.xx) / es.A.xy);         // A.xx * 1 + A.xy * y = l2 * 1
  if (prod_cross_2d(es.v1, es.v2) < 0) { // "right-handed", in the sense that z-axis is pointing "outward" the screen
    multiply_vec2d(-1.0, &es.v2);
  }
  es.u1 = es.v1; normalize_vec2d(&es.u1);
  es.u2 = es.v2; normalize_vec2d(&es.u2);
  return es;
}

int check_eigs(EigenSys es) {
  int check = 0;
  // check eigenvalues
  VEC2D lhs, rhs;
  lhs = rotate_vec2d(es.A, es.u1);
  rhs = es.u1; multiply_vec2d(es.l1, &rhs);
  if (fabs(lhs.x - rhs.x) < EPSILON) check += 1;
  if (fabs(lhs.y - rhs.y) < EPSILON) check += 10;
  lhs = rotate_vec2d(es.A, es.u2);
  rhs = es.u2; multiply_vec2d(es.l2, &rhs);
  if (fabs(lhs.x - rhs.x) < EPSILON) check += 100;
  if (fabs(lhs.y - rhs.y) < EPSILON) check += 1000;
  return check;
}

#ifdef _MSC_VER
// display function
void print_vec2d(VEC2D v, const char * name) {
  printf("%s = \n| %6.3f  %6.3f |\n", name, v.x, v.y);
}

void print_vec3d(VEC3D v, const char * name) {
  printf("%s = \n| %6.3f  %6.3f  %6.3f |\n", name, v.x, v.y, v.z);
}

void print_mat2d(MAT2D m, const char * name) {
  printf("%s = \n| %6.3f %6.3f |\n| %6.3f %6.3f |\n", name, m.xx, m.xy, m.yx, m.yy);
}

void print_mat3d(MAT3D m, const char * name) {
  printf("%s = \n| %6.3f %6.3f %6.3f |\n| %6.3f %6.3f %6.3f |\n| %6.3f %6.3f %6.3f |\n", name, m.xx, m.xy, m.xz, m.yx, m.yy, m.yz, m.zx, m.zy, m.zz);
}

void print_vec6d(VEC6D v, const char * name) {
  printf("%s = \n| %6.3f %6.3f %6.3f  ,  %6.3f %6.3f %6.3f |\n", name, v.a.x, v.a.y, v.a.z, v.g.x, v.g.y, v.g.z);
}

void print_mat6d(MAT6D m, const char * name) {
  printf("%s = \n", name);
  printf("| %6.3f %6.3f %6.3f \t %6.3f %6.3f %6.3f |\n", m.A[0][0].xx, m.A[0][0].xy, m.A[0][0].xz, m.A[0][1].xx, m.A[0][1].xy, m.A[0][1].xz);
  printf("| %6.3f %6.3f %6.3f \t %6.3f %6.3f %6.3f |\n", m.A[0][0].yx, m.A[0][0].yy, m.A[0][0].yz, m.A[0][1].yx, m.A[0][1].yy, m.A[0][1].yz);
  printf("| %6.3f %6.3f %6.3f \t %6.3f %6.3f %6.3f |\n", m.A[0][0].zx, m.A[0][0].zy, m.A[0][0].zz, m.A[0][1].zx, m.A[0][1].zy, m.A[0][1].zz);
  printf("| %8s                                    |\n", " ");
  printf("| %6.3f %6.3f %6.3f \t %6.3f %6.3f %6.3f |\n", m.A[1][0].xx, m.A[1][0].xy, m.A[1][0].xz, m.A[1][1].xx, m.A[1][1].xy, m.A[1][1].xz);
  printf("| %6.3f %6.3f %6.3f \t %6.3f %6.3f %6.3f |\n", m.A[1][0].yx, m.A[1][0].yy, m.A[1][0].yz, m.A[1][1].yx, m.A[1][1].yy, m.A[1][1].yz);
  printf("| %6.3f %6.3f %6.3f \t %6.3f %6.3f %6.3f |\n", m.A[1][0].zx, m.A[1][0].zy, m.A[1][0].zz, m.A[1][1].zx, m.A[1][1].zy, m.A[1][1].zz);
}

void print_eigs(EigenSys es, const char * description) {
  printf("%s\n", description);
  print_mat2d(es.A, "A ");
  printf("\nl1 = %6.3f\n", es.l1);
  print_vec2d(es.v1, "v1 ");
  print_vec2d(es.u1, "u1 ");
  printf("\nl2 = %6.3f\n", es.l2);
  print_vec2d(es.v2, "v2 ");
  print_vec2d(es.u2, "u2 ");
  printf("\ncheck (1111 is ok) = %04d\n", check_eigs(es));
}
#endif
