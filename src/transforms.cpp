#include "transforms.h"

#include "CGL/matrix3x3.h"
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"
#include <cmath>

namespace CGL {

Vector2D operator*(const Matrix3x3 &m, const Vector2D &v) {
	Vector3D mv = m * Vector3D(v.x, v.y, 1);
	return Vector2D(mv.x / mv.z, mv.y / mv.z);
}

Matrix3x3 translate(float dx, float dy) {
	return Matrix3x3(
		1, 0, dx,
		0, 1, dy,
		0, 0, 1
	);
}

Matrix3x3 scale(float sx, float sy) {
	return Matrix3x3(
		sx, 0, 0,
		0, sy, 0,
		0, 0, 1
	);
}

// The input argument is in degrees counterclockwise
Matrix3x3 rotate(float deg) {
	float rad = radians(deg);
	return Matrix3x3(
		cos(rad), -sin(rad), 0,
		sin(rad), cos(rad),  0,
		0,           0,           1
	);
}

}
