#pragma once
#include <cmath>

// 求两点之间的距离
float distance3D(float x1, float y1, float z1, float x2, float y2, float z2) {
  return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2));
}

//求三个空间点的重心
void centerOfMass(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3, float &x, float &y, float &z) {
  x = (x1 + x2 + x3) / 3;
  y = (y1 + y2 + y3) / 3;
  z = (z1 + z2 + z3) / 3;
}

// 求重心到三个空间点的距离
void distanceToCenterOfMass(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3, float &d1, float &d2, float &d3) {
  float x, y, z;
  centerOfMass(x1, y1, z1, x2, y2, z2, x3, y3, z3, x, y, z);
  d1 = distance3D(x1, y1, z1, x, y, z);
  d2 = distance3D(x2, y2, z2, x, y, z);
  d3 = distance3D(x3, y3, z3, x, y, z);
}

// 求重心到三个顶点的距离和
float sumDistanceToCenterOfMass(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3) {
  float d1, d2, d3;
  distanceToCenterOfMass(x1, y1, z1, x2, y2, z2, x3, y3, z3, d1, d2, d3);
  return d1 + d2 + d3;
}


// 求三角形的内心
void incenter(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3, float &x, float &y, float &z) {
  float d1, d2, d3;
  distanceToCenterOfMass(x1, y1, z1, x2, y2, z2, x3, y3, z3, d1, d2, d3);
  x = (d1 * x1 + d2 * x2 + d3 * x3) / (d1 + d2 + d3);
  y = (d1 * y1 + d2 * y2 + d3 * y3) / (d1 + d2 + d3);
  z = (d1 * z1 + d2 * z2 + d3 * z3) / (d1 + d2 + d3);
}

// 求三角形的内心到三个顶点的距离
void distanceToIncenter(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3, float &d1, float &d2, float &d3) {
  float x, y, z;
  incenter(x1, y1, z1, x2, y2, z2, x3, y3, z3, x, y, z);
  d1 = distance3D(x1, y1, z1, x, y, z);
  d2 = distance3D(x2, y2, z2, x, y, z);
  d3 = distance3D(x3, y3, z3, x, y, z);
}

// 求三角形的内心到三个顶点的距离和
float sumDistanceToIncenter(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3) {
  float d1, d2, d3;
  distanceToIncenter(x1, y1, z1, x2, y2, z2, x3, y3, z3, d1, d2, d3);
  return d1 + d2 + d3;
}

// 求费马点
void fermatPoint(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3, float &x, float &y, float &z) {
  float d1, d2, d3;
  distanceToCenterOfMass(x1, y1, z1, x2, y2, z2, x3, y3, z3, d1, d2, d3);
  float a = distance3D(x2, y2, z2, x3, y3, z3);
  float b = distance3D(x1, y1, z1, x3, y3, z3);
  float c = distance3D(x1, y1, z1, x2, y2, z2);
  float cosA = (b * b + c * c - a * a) / (2 * b * c);
  float cosB = (a * a + c * c - b * b) / (2 * a * c);
  float cosC = (a * a + b * b - c * c) / (2 * a * b);
  float xA, yA, zA, xB, yB, zB, xC, yC, zC;
  centerOfMass(x1, y1, z1, x2, y2, z2, x3, y3, z3, xA, yA, zA);
  xB = x1 + d1 * cosA;
  yB = y1 + d1 * sqrt(1 - cosA * cosA);
  zB = z1;
  xC = x1 + d1 * cosB;
  yC = y1 - d1 * sqrt(1 - cosB * cosB);
  zC = z1;
  float dA, dB, dC;
  distanceToCenterOfMass(xA, yA, zA, xB, yB, zB, xC, yC, zC, dA, dB, dC);
  if (dA < dB && dA < dC) {
    x = xA;
    y = yA;
    z = zA;
  } else if (dB < dA && dB < dC) {
    x = xB;
    y = yB;
    z = zB;
  } else {
    x = xC;
    y = yC;
    z = zC;
  }
}

// 求费马点到三个顶点的距离和
float sumDistanceToFermatPoint(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3) {
  float x, y, z;
  fermatPoint(x1, y1, z1, x2, y2, z2, x3, y3, z3, x, y, z);
  return distance3D(x1, y1, z1, x, y, z) + distance3D(x2, y2, z2, x, y, z) + distance3D(x3, y3, z3, x, y, z);
}


float minDistanceTLikeStructure(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3) {

  float d1, d2, d3;
  d1 = distance3D(x1, y1, z1, x2, y2, z2);
  d1 += distance3D((x1 + x2)/2, (y1 + y2)/2, (z1 + z2)/2, x3, y3, z3);
  d2 = distance3D(x1, y1, z1, x3, y3, z3);
  d2 += distance3D((x1 + x3)/2, (y1 + y3)/2, (z1 + z3)/2, x2, y2, z2);
  d3 = distance3D(x2, y2, z2, x3, y3, z3);
  d3 += distance3D((x2 + x3)/2, (y2 + y3)/2, (z2 + z3)/2, x1, y1, z1);

  if (d1 < d2 && d1 < d3) {
    return d1;
  } else if (d2 < d1 && d2 < d3) {
    return d2;
  } else {
    return d3;
  }
}
