#include <Novice.h>
#include <imgui.h>
#include <cmath>
#include "Matrix4x4.h"
#include "assert.h"
#include "Vector3.h"
#include <algorithm>

//加算
Vector3 Add(const Vector3& v1, const Vector3& v2) 
{
	Vector3 result;
	result.x = v1.x + v2.x;
	result.y = v1.y + v2.y;
	result.z = v1.z + v2.z;
	return result;
}

// 減算
Vector3 Subtract(const Vector3& v1, const Vector3& v2) 
{
	Vector3 result;
	result.x = v1.x - v2.x;
	result.y = v1.y - v2.y;
	result.z = v1.z - v2.z;
	return result;
}

// スカラー倍
Vector3 Multiply(float scalar, const Vector3& v)
{
	Vector3 result;
	result.x = scalar * v.x;
	result.y = scalar * v.y;
	result.z = scalar * v.z;
	return result;
}

// 行列の積
Matrix4x4 Multiply(const Matrix4x4& m1, const Matrix4x4& m2)
{
	Matrix4x4 result;

	for (int row = 0; row < 4; row++) {
		for (int column = 0; column < 4; column++) {
			result.m[row][column] = m1.m[row][0] * m2.m[0][column] + m1.m[row][1] * m2.m[1][column] + m1.m[row][2] * m2.m[2][column] + m1.m[row][3] * m2.m[3][column];
		}
	}
	return result;
}

// 内積
float Dot(const Vector3& v1, const Vector3& v2) 
{
	float result;
	result = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
	return result;
}

// X軸回転行列
Matrix4x4 MakeRotateXMatrix(float radian) 
{
	Matrix4x4 result;
	result = {
		1.0f, 0.0f, 0.0f, 0.0f,
		0.0f, std::cos(radian), std::sin(radian), 0.0f,
		0.0f, -std::sin(radian), std::cos(radian), 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f,
	};

	return result;
}

// Y軸回転行列
Matrix4x4 MakeRotateYMatrix(float radian) 
{
	Matrix4x4 result;
	result = {
		std::cos(radian), 0.0f, -std::sin(radian), 0.0f,
		0.0f, 1.0f, 0.0f, 0.0f, 
		std::sin(radian), 0.0f, std::cos(radian), 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f,
	};

	return result;
}

// Z軸回転行列
Matrix4x4 MakeRotateZMatrix(float radian)
{
	Matrix4x4 result;
	result = {
		std::cos(radian), std::sin(radian), 0.0f, 0.0f,
		-std::sin(radian), std::cos(radian), 0.0f, 0.0f,
		0.0f, 0.0f, 1.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f,
	};

	return result;
}

// 3次元のアフィン変換行列
Matrix4x4 MakeAffineMatrix(const Vector3& scale, const Vector3& rotate, const Vector3& translate)
{
	Matrix4x4 result;
	Matrix4x4 rotateXMatrix = MakeRotateXMatrix(rotate.x);
	Matrix4x4 rotateYMatrix = MakeRotateYMatrix(rotate.y);
	Matrix4x4 rotateZMatrix = MakeRotateZMatrix(rotate.z);
	Matrix4x4 rotateXYZMatrix = Multiply(rotateXMatrix, Multiply(rotateYMatrix, rotateZMatrix)); // R

	result = {
		scale.x * rotateXYZMatrix.m[0][0],
		scale.x * rotateXYZMatrix.m[0][1],
		scale.x * rotateXYZMatrix.m[0][2],
		0.0f,
		scale.y * rotateXYZMatrix.m[1][0],
		scale.y * rotateXYZMatrix.m[1][1],
		scale.y * rotateXYZMatrix.m[1][2],
		0.0f,
		scale.z * rotateXYZMatrix.m[2][0],
		scale.z * rotateXYZMatrix.m[2][1],
		scale.z * rotateXYZMatrix.m[2][2],
		0.0f,
		translate.x,
		translate.y,
		translate.z,
		1.0f };

	return result;
}

// 逆行列
Matrix4x4 Inverse(const Matrix4x4& m) 
{
	Matrix4x4 result;
	float a;
	Matrix4x4 b;

	a = 1 /
		(m.m[0][0] * m.m[1][1] * m.m[2][2] * m.m[3][3] + m.m[0][0] * m.m[1][2] * m.m[2][3] * m.m[3][1] + m.m[0][0] * m.m[1][3] * m.m[2][1] * m.m[3][2] - m.m[0][0] * m.m[1][3] * m.m[2][2] * m.m[3][1] -
			m.m[0][0] * m.m[1][2] * m.m[2][1] * m.m[3][3] - m.m[0][0] * m.m[1][1] * m.m[2][3] * m.m[3][2] - m.m[0][1] * m.m[1][0] * m.m[2][2] * m.m[3][3] - m.m[0][2] * m.m[1][0] * m.m[2][3] * m.m[3][1] -
			m.m[0][3] * m.m[1][0] * m.m[2][1] * m.m[3][2] + m.m[0][3] * m.m[1][0] * m.m[2][2] * m.m[3][1] + m.m[0][2] * m.m[1][0] * m.m[2][1] * m.m[3][3] + m.m[0][1] * m.m[1][0] * m.m[2][3] * m.m[3][2] +
			m.m[0][1] * m.m[1][2] * m.m[2][0] * m.m[3][3] + m.m[0][2] * m.m[1][3] * m.m[2][0] * m.m[3][1] + m.m[0][3] * m.m[1][1] * m.m[2][0] * m.m[3][2] - m.m[0][3] * m.m[1][2] * m.m[2][0] * m.m[3][1] -
			m.m[0][2] * m.m[1][1] * m.m[2][0] * m.m[3][3] - m.m[0][1] * m.m[1][3] * m.m[2][0] * m.m[3][2] - m.m[0][1] * m.m[1][2] * m.m[2][3] * m.m[3][0] - m.m[0][2] * m.m[1][3] * m.m[2][1] * m.m[3][0] -
			m.m[0][3] * m.m[1][1] * m.m[2][2] * m.m[3][0] + m.m[0][3] * m.m[1][2] * m.m[2][1] * m.m[3][0] + m.m[0][2] * m.m[1][1] * m.m[2][3] * m.m[3][0] + m.m[0][1] * m.m[1][3] * m.m[2][2] * m.m[3][0]);

	b.m[0][0] = m.m[1][1] * m.m[2][2] * m.m[3][3] + m.m[1][2] * m.m[2][3] * m.m[3][1] + m.m[1][3] * m.m[2][1] * m.m[3][2] - m.m[1][3] * m.m[2][2] * m.m[3][1] - m.m[1][2] * m.m[2][1] * m.m[3][3] -
		m.m[1][1] * m.m[2][3] * m.m[3][2];
	b.m[0][1] = -m.m[0][1] * m.m[2][2] * m.m[3][3] - m.m[0][2] * m.m[2][3] * m.m[3][1] - m.m[0][3] * m.m[2][1] * m.m[3][2] + m.m[0][3] * m.m[2][2] * m.m[3][1] + m.m[0][2] * m.m[2][1] * m.m[3][3] +
		m.m[0][1] * m.m[2][3] * m.m[3][2];
	b.m[0][2] = m.m[0][1] * m.m[1][2] * m.m[3][3] + m.m[0][2] * m.m[1][3] * m.m[3][1] + m.m[0][3] * m.m[1][1] * m.m[3][2] - m.m[0][3] * m.m[1][2] * m.m[3][1] - m.m[0][2] * m.m[1][1] * m.m[3][3] -
		m.m[0][1] * m.m[1][3] * m.m[3][2];
	b.m[0][3] = -m.m[0][1] * m.m[1][2] * m.m[2][3] - m.m[0][2] * m.m[1][3] * m.m[2][1] - m.m[0][3] * m.m[1][1] * m.m[2][2] + m.m[0][3] * m.m[1][2] * m.m[2][1] + m.m[0][2] * m.m[1][1] * m.m[2][3] +
		m.m[0][1] * m.m[1][3] * m.m[2][2];

	b.m[1][0] = -m.m[1][0] * m.m[2][2] * m.m[3][3] - m.m[1][2] * m.m[2][3] * m.m[3][0] - m.m[1][3] * m.m[2][0] * m.m[3][2] + m.m[1][3] * m.m[2][2] * m.m[3][0] + m.m[1][2] * m.m[2][0] * m.m[3][3] +
		m.m[1][0] * m.m[2][3] * m.m[3][2];
	b.m[1][1] = m.m[0][0] * m.m[2][2] * m.m[3][3] + m.m[0][2] * m.m[2][3] * m.m[3][0] + m.m[0][3] * m.m[2][0] * m.m[3][2] - m.m[0][3] * m.m[2][2] * m.m[3][0] - m.m[0][2] * m.m[2][0] * m.m[3][3] -
		m.m[0][0] * m.m[2][3] * m.m[3][2];
	b.m[1][2] = -m.m[0][0] * m.m[1][2] * m.m[3][3] - m.m[0][2] * m.m[1][3] * m.m[3][0] - m.m[0][3] * m.m[1][0] * m.m[3][2] + m.m[0][3] * m.m[1][2] * m.m[3][0] + m.m[0][2] * m.m[1][0] * m.m[3][3] +
		m.m[0][0] * m.m[1][3] * m.m[3][2];
	b.m[1][3] = +m.m[0][0] * m.m[1][2] * m.m[2][3] + m.m[0][2] * m.m[1][3] * m.m[2][0] + m.m[0][3] * m.m[1][0] * m.m[2][2] - m.m[0][3] * m.m[1][2] * m.m[2][0] - m.m[0][2] * m.m[1][0] * m.m[2][3] -
		m.m[0][0] * m.m[1][3] * m.m[2][2];

	b.m[2][0] = m.m[1][0] * m.m[2][1] * m.m[3][3] + m.m[1][1] * m.m[2][3] * m.m[3][0] + m.m[1][3] * m.m[2][0] * m.m[3][1] - m.m[1][3] * m.m[2][1] * m.m[3][0] - m.m[1][1] * m.m[2][0] * m.m[3][3] -
		m.m[1][0] * m.m[2][3] * m.m[3][1];
	b.m[2][1] = -m.m[0][0] * m.m[2][1] * m.m[3][3] - m.m[0][1] * m.m[2][3] * m.m[3][0] - m.m[0][3] * m.m[2][0] * m.m[3][1] + m.m[0][3] * m.m[2][1] * m.m[3][0] + m.m[0][1] * m.m[2][0] * m.m[3][3] +
		m.m[0][0] * m.m[2][3] * m.m[3][1];
	b.m[2][2] = m.m[0][0] * m.m[1][1] * m.m[3][3] + m.m[0][1] * m.m[1][3] * m.m[3][0] + m.m[0][3] * m.m[1][0] * m.m[3][1] - m.m[0][3] * m.m[1][1] * m.m[3][0] - m.m[0][1] * m.m[1][0] * m.m[3][3] -
		m.m[0][0] * m.m[1][3] * m.m[3][1];
	b.m[2][3] = -m.m[0][0] * m.m[1][1] * m.m[2][3] - m.m[0][1] * m.m[1][3] * m.m[2][0] - m.m[0][3] * m.m[1][0] * m.m[2][1] + m.m[0][3] * m.m[1][1] * m.m[2][0] + m.m[0][1] * m.m[1][0] * m.m[2][3] +
		m.m[0][0] * m.m[1][3] * m.m[2][1];

	b.m[3][0] = -m.m[1][0] * m.m[2][1] * m.m[3][2] - m.m[1][1] * m.m[2][2] * m.m[3][0] - m.m[1][2] * m.m[2][0] * m.m[3][1] + m.m[1][2] * m.m[2][1] * m.m[3][0] + m.m[1][1] * m.m[2][0] * m.m[3][2] +
		m.m[1][0] * m.m[2][2] * m.m[3][1];
	b.m[3][1] = +m.m[0][0] * m.m[2][1] * m.m[3][2] + m.m[0][1] * m.m[2][2] * m.m[3][0] + m.m[0][2] * m.m[2][0] * m.m[3][1] - m.m[0][2] * m.m[2][1] * m.m[3][0] - m.m[0][1] * m.m[2][0] * m.m[3][2] -
		m.m[0][0] * m.m[2][2] * m.m[3][1];
	b.m[3][2] = -m.m[0][0] * m.m[1][1] * m.m[3][2] - m.m[0][1] * m.m[1][2] * m.m[3][0] - m.m[0][2] * m.m[1][0] * m.m[3][1] + m.m[0][2] * m.m[1][1] * m.m[3][0] + m.m[0][1] * m.m[1][0] * m.m[3][2] +
		m.m[0][0] * m.m[1][2] * m.m[3][1];
	b.m[3][3] = m.m[0][0] * m.m[1][1] * m.m[2][2] + m.m[0][1] * m.m[1][2] * m.m[2][0] + m.m[0][2] * m.m[1][0] * m.m[2][1] - m.m[0][2] * m.m[1][1] * m.m[2][0] - m.m[0][1] * m.m[1][0] * m.m[2][2] -
		m.m[0][0] * m.m[1][2] * m.m[2][1];

	result.m[0][0] = a * b.m[0][0];
	result.m[0][1] = a * b.m[0][1];
	result.m[0][2] = a * b.m[0][2];
	result.m[0][3] = a * b.m[0][3];
	result.m[1][0] = a * b.m[1][0];
	result.m[1][1] = a * b.m[1][1];
	result.m[1][2] = a * b.m[1][2];
	result.m[1][3] = a * b.m[1][3];
	result.m[2][0] = a * b.m[2][0];
	result.m[2][1] = a * b.m[2][1];
	result.m[2][2] = a * b.m[2][2];
	result.m[2][3] = a * b.m[2][3];
	result.m[3][0] = a * b.m[3][0];
	result.m[3][1] = a * b.m[3][1];
	result.m[3][2] = a * b.m[3][2];
	result.m[3][3] = a * b.m[3][3];

	return result;
}

// 透視投影行列
Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio, float nearClip, float farClip)
{
	Matrix4x4 result;
	result = {
		1 / aspectRatio * (1 / std::tan(fovY / 2)), 0.0f, 0.0f, 0.0f, 0.0f, 1 / std::tan(fovY / 2), 0.0f, 0.0f, 0.0f, 0.0f, farClip / (farClip - nearClip), 1.0f, 0.0f, 0.0f,
		-nearClip * farClip / (farClip - nearClip), 0.0f,
	};

	return result;
}

Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height, float minDepth, float maxDepth)
{
	Matrix4x4 result;
	result = { width / 2, 0.0f, 0.0f, 0.0f, 0.0f, -height / 2, 0.0f, 0.0f, 0.0f, 0.0f, maxDepth - minDepth, 0.0f, left + (width / 2), top + (height / 2), minDepth, 1.0f };

	return result;
}

// 座標変換
Vector3 Transform(const Vector3& vector, const Matrix4x4& matrix)
{
	Vector3 result;
	result.x = vector.x * matrix.m[0][0] + vector.y * matrix.m[1][0] + vector.z * matrix.m[2][0] + 1.0f * matrix.m[3][0];
	result.y = vector.x * matrix.m[0][1] + vector.y * matrix.m[1][1] + vector.z * matrix.m[2][1] + 1.0f * matrix.m[3][1];
	result.z = vector.x * matrix.m[0][2] + vector.y * matrix.m[1][2] + vector.z * matrix.m[2][2] + 1.0f * matrix.m[3][2];
	float w = vector.x * matrix.m[0][3] + vector.y * matrix.m[1][3] + vector.z * matrix.m[2][3] + 1.0f * matrix.m[3][3];
	assert(w != 0.0f);
	result.x /= w;
	result.y /= w;
	result.z /= w;
	return result;
}

// 正規化|
Vector3 Normalize(const Vector3& v)
{
	Vector3 result;
	result.x = v.x / sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
	result.y = v.y / sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
	result.z = v.z / sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
	return result;
}

// クロス値
Vector3 Cross(const Vector3& v1, const Vector3& v2)
{
	Vector3 result;
	result = { v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x };

	return result;
}

// 距離
float Length(const Vector3& v1, const Vector3& v2)
{
	Vector3 distance;
	float d;

	distance = Subtract(v2, v1);
	d = sqrtf(powf(distance.x, 2) + powf(distance.y, 2) + powf(distance.z, 2));

	return d;
}


static const int kWindowWidth = 1280;
static const int kWindowHeight = 720;

//グリッド
void DrawGrid(const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix)
{
	const float kGridHalfwidth = 2.0f;                                     //Gridの半分
	const uint32_t kSubdivision = 10;                                      //分割数
	const float kGridEvery = (kGridHalfwidth * 2.0f) / float(kSubdivision);//1つ分の長さ

	//奥から手前への線を順々に引いていく
	for (uint32_t xIndex = 0; xIndex <= kSubdivision; ++xIndex)
	{
		//上の情報を使ってワールド座標系上の始点と終点を求める
		Vector3 zStart, zEnd;

		zStart = Vector3(xIndex * kGridEvery / 2 - kGridHalfwidth + 1, 0, 1);
		zEnd = Vector3(xIndex * kGridEvery / 2 - kGridHalfwidth + 1, 0, -3);

		//スクリーン座標系まで変換をかける
		Matrix4x4 startWorldMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, zStart);
		Matrix4x4 endWorldMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, zEnd);

		Matrix4x4 startwvpMatrix = Multiply(startWorldMatrix, viewProjectionMatrix);
		Matrix4x4 endwvpMatrix = Multiply(startWorldMatrix, viewProjectionMatrix);

		Vector3 startLocal = Transform(zStart, startwvpMatrix);
		Vector3 endLocal = Transform(zEnd, endwvpMatrix);

		Vector3 startScreen = Transform(startLocal, viewportMatrix);
		Vector3 endScreen = Transform(endLocal, viewportMatrix);

		//変換した座標を使って表示、色は薄い灰色(0xAAAAAAFF)。原点は黒
		if (xIndex == kSubdivision / 2)
		{
			Novice::DrawLine((int)startScreen.x, (int)startScreen.y, (int)endScreen.x, (int)endScreen.y, BLACK);
		}
		else
		{
			Novice::DrawLine((int)startScreen.x, (int)startScreen.y, (int)endScreen.x, (int)endScreen.y, 0xAAAAAAFF);
		}

	}

	//左から右も同じように順々に引いていく
	for (uint32_t zIndex = 0; zIndex <= kSubdivision; ++zIndex)
	{
		//奥から手前が左右に変わるだけ
		//上の情報を使ってワールド座標系上の始点と終点を求める
		Vector3 xStart, xEnd;

		xStart = Vector3(1, 0, zIndex * kGridEvery / 2 - kGridHalfwidth + 1);
		xEnd = Vector3(-1, 0, zIndex * kGridEvery / 2 - kGridHalfwidth + 1);

		//スクリーン座標系まで変換をかける
		Matrix4x4 startWorldMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, xStart);
		Matrix4x4 endWorldMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, xEnd);

		Matrix4x4 startwvpMatrix = Multiply(startWorldMatrix, viewProjectionMatrix);
		Matrix4x4 endwvpMatrix = Multiply(endWorldMatrix, viewProjectionMatrix);

		Vector3 startLocal = Transform(xStart, startwvpMatrix);
		Vector3 endLocal = Transform(xEnd, endwvpMatrix);

		Vector3 startScreen = Transform(startLocal, viewportMatrix);
		Vector3 endScreen = Transform(endLocal, viewportMatrix);

		//変換した座標を使って表示、色は薄い灰色(0xAAAAAAFF)。原点は黒
		if (zIndex == kSubdivision / 2)
		{
			Novice::DrawLine((int)startScreen.x, (int)startScreen.y, (int)endScreen.x, (int)endScreen.y, BLACK);
		}
		else
		{
			Novice::DrawLine((int)startScreen.x, (int)startScreen.y, (int)endScreen.x, (int)endScreen.y, 0xAAAAAAFF);
		}
	}
}

struct Line
{
	Vector3 origin;//始点
	Vector3 diff;  //終点への差分ベクトル
};

struct Ray
{
	Vector3 origin;//始点
	Vector3 diff;  //終点への差分ベクトル
};

struct Segment
{
	Vector3 origin;//始点
	Vector3 diff;  //終点への差分ベクトル
};

struct Sphere
{
	Vector3 center;//!<中心点
	float radius;  //!<半径
};

struct Plane
{
	Vector3 normal; // !< 法線
	float distance; // !< 距離
};

struct Triangle
{
	Vector3 vertices[3]; //!< 頂点
};

struct AABB
{
	Vector3 min;//!<最小点
	Vector3 max;//!< 最大点
};


///////////色々/////////////

Vector3 Project(const Vector3& v1, const Vector3& v2)
{
	Vector3 result;
	float t = Dot(v1, v2) / (sqrtf(powf(Dot(v2, v2), 2)));

	result = Multiply(t, v2);

	return result;
}

Vector3 ClosestPoint(const Vector3 point, const Segment& segment)
{
	Vector3 proja = Project(Subtract(point, segment.origin), Subtract(Add(segment.origin, segment.diff), segment.origin));
	Vector3 cp = Add(segment.origin, proja);
	// d = sqrtf((point.x - cp.x) + (point.y - cp.y) + (point.z - cp.z));

	return cp;
}

Vector3 Perpendicular(const Vector3& vector)
{
	if (vector.x != 0.0f || vector.y != 0.0f)
	{
		return { -vector.y,vector.x,0.0f };
	}
	return { 0.0f,-vector.z,vector.y };
}

//衝突点の媒介変数
float CreateParametricVariable(const float d, const Vector3& n, const Segment& segment)
{
	float t;
	
	t = d - Dot(segment.origin, n) / Dot(segment.diff, n);

	return t;
}

/////////////////////////////



////////////描画/////////////

void DrawSphere(const Sphere& sphere, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color)
{
	const float pi = 3.1415926535f;
	const uint32_t kSubdivision = 16;//分割数
	const float kLatEvery = pi / kSubdivision;   //経度分割1つ分の角度
	const float kLonEvery = (2 * pi) / kSubdivision;   //緯度分割1つ分の角度

	//緯度の方向に分割 -π/2 ～ π/2
	for (uint32_t latIndex = 0; latIndex < kSubdivision; ++latIndex)
	{
		sphere;
		float lat = -pi / 2.0f + kLatEvery * latIndex;//現在の緯度

		//経度の方向に分割 0 ～ 2π
		for (uint32_t lonIndex = 0; lonIndex < kSubdivision; ++lonIndex)
		{
			float lon = lonIndex * kLonEvery;//現在の経度

			//world座標系でのa,b,c,を求める
			Vector3 a, b, c;

			a = {
				sphere.radius * std::cos(lat) * std::cos(lon),
				sphere.radius * std::sin(lat),
				sphere.radius * std::cos(lat) * std::sin(lon)
			};

			b = {
				sphere.radius * std::cos(lat + (pi / kSubdivision)) * std::cos(lon),
				sphere.radius * std::sin(lat + (pi / kSubdivision)),
				sphere.radius * std::cos(lat + (pi / kSubdivision)) * std::sin(lon)
			};

			c = {
				sphere.radius * std::cos(lat) * std::cos(lon + ((pi * 2) / kSubdivision)),
				sphere.radius * std::sin(lat),
				sphere.radius * std::cos(lat) * std::sin(lon + ((pi * 2) / kSubdivision))
			};

			//a,b,cをScreen座標系まで変換...
			Matrix4x4 aWorldMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, sphere.center);
			Matrix4x4 bWorldMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, sphere.center);
			Matrix4x4 cWorldMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, sphere.center);

			Matrix4x4 awvpMatrix = Multiply(aWorldMatrix, viewProjectionMatrix);
			Matrix4x4 bwvpMatrix = Multiply(bWorldMatrix, viewProjectionMatrix);
			Matrix4x4 cwvpMatrix = Multiply(cWorldMatrix, viewProjectionMatrix);

			Vector3 aLocal = Transform(a, awvpMatrix);
			Vector3 bLocal = Transform(b, bwvpMatrix);
			Vector3 cLocal = Transform(c, cwvpMatrix);

			Vector3 aScreen = Transform(aLocal, viewportMatrix);
			Vector3 bScreen = Transform(bLocal, viewportMatrix);
			Vector3 cScreen = Transform(cLocal, viewportMatrix);

			//ab,bcで線を引く
			Novice::DrawLine((int)aScreen.x, (int)aScreen.y, (int)bScreen.x, (int)bScreen.y, color);
			Novice::DrawLine((int)aScreen.x, (int)aScreen.y, (int)cScreen.x, (int)cScreen.y, color);
		}

	}

}

void DrawPlane(const Plane& plane, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color)
{
	Vector3 center = Multiply(plane.distance, plane.normal); // 1
	Vector3 perpendiculars[4];
	perpendiculars[0] = Normalize(Perpendicular(plane.normal)); // 2
	perpendiculars[1] = { -perpendiculars[0].x,-perpendiculars[0].y,-perpendiculars[0].z }; // 3
	perpendiculars[2] = Cross(plane.normal, perpendiculars[0]); // 4
	perpendiculars[3] = { -perpendiculars[2].x,-perpendiculars[2].y,-perpendiculars[2].z }; // 5

	// 6
	Vector3 points[4];
	for (int32_t index = 0; index < 4; ++index)
	{
		Vector3 extend = Multiply(2.0f, perpendiculars[index]);
		Vector3 point = Add(center, extend);
		points[index] = Transform(Transform(point, viewProjectionMatrix), viewportMatrix);
	}
	//points をそれぞれ結んでDraw で矩形を描画する。DrawTringleを使って塗りつぶしても良いが、DepthがないのでMT3では分かりずらい
	Novice::DrawLine((int)points[0].x, (int)points[0].y, (int)points[2].x, (int)points[2].y, color);
	Novice::DrawLine((int)points[0].x, (int)points[0].y, (int)points[3].x, (int)points[3].y, color);
	Novice::DrawLine((int)points[1].x, (int)points[1].y, (int)points[2].x, (int)points[2].y, color);
	Novice::DrawLine((int)points[1].x, (int)points[1].y, (int)points[3].x, (int)points[3].y, color);
}

void DrawTriangle(const Triangle& triangle, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color)
{
	//スクリーン座標に変換
	Vector3 transform[3];
	Vector3 screen[3];

	for (int i = 0; i < 3; ++i) {
		transform[i] = Transform(triangle.vertices[i], viewProjectionMatrix);
		screen[i] = Transform(transform[i], viewportMatrix);
	}

	//三角形の描画
	Novice::DrawTriangle((int)screen[0].x, (int)screen[0].y,(int)screen[1].x, (int)screen[1].y,(int)screen[2].x, (int)screen[2].y,color, kFillModeWireFrame);
}

void DrawAABB(const AABB& aabb, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color)
{
	//8頂点求める
	Vector3 flontLeftUp{ aabb.min.x,aabb.max.y,aabb.min.z };
	Vector3 flontLeftUnder{ aabb.min.x,aabb.min.y,aabb.min.z };
	Vector3 flontRightUp{ aabb.max.x,aabb.max.y,aabb.min.z };
	Vector3 flontRightUnder{ aabb.max.x,aabb.min.y,aabb.min.z };

	Vector3 backLeftUp{ aabb.min.x,aabb.max.y,aabb.max.z };
	Vector3 backLeftUnder{ aabb.min.x,aabb.min.y,aabb.max.z };
	Vector3 backRightUp{ aabb.max.x,aabb.max.y,aabb.max.z };
	Vector3 backRightUnder{ aabb.max.x,aabb.min.y,aabb.max.z };

	//スクリーン座標に変換
	Vector3 transform[8];
	Vector3 screen[8];

	transform[0] = Transform(flontLeftUp, viewProjectionMatrix);
	transform[1] = Transform(flontLeftUnder, viewProjectionMatrix);
	transform[2] = Transform(flontRightUp, viewProjectionMatrix);
	transform[3] = Transform(flontRightUnder, viewProjectionMatrix);

	transform[4] = Transform(backLeftUp, viewProjectionMatrix);
	transform[5] = Transform(backLeftUnder, viewProjectionMatrix);
	transform[6] = Transform(backRightUp, viewProjectionMatrix);
	transform[7] = Transform(backRightUnder, viewProjectionMatrix);

	screen[0] = Transform(transform[0], viewportMatrix);
	screen[1] = Transform(transform[1], viewportMatrix);
	screen[2] = Transform(transform[2], viewportMatrix);
	screen[3] = Transform(transform[3], viewportMatrix);

	screen[4] = Transform(transform[4], viewportMatrix);
	screen[5] = Transform(transform[5], viewportMatrix);
	screen[6] = Transform(transform[6], viewportMatrix);
	screen[7] = Transform(transform[7], viewportMatrix);

	//8頂点線で結ぶ
	Novice::DrawLine((int)screen[0].x, (int)screen[0].y, (int)screen[1].x, (int)screen[1].y, color);
	Novice::DrawLine((int)screen[0].x, (int)screen[0].y, (int)screen[2].x, (int)screen[2].y, color);
	Novice::DrawLine((int)screen[0].x, (int)screen[0].y, (int)screen[4].x, (int)screen[4].y, color);

	Novice::DrawLine((int)screen[6].x, (int)screen[6].y, (int)screen[2].x, (int)screen[2].y, color);
	Novice::DrawLine((int)screen[6].x, (int)screen[6].y, (int)screen[4].x, (int)screen[4].y, color);
	Novice::DrawLine((int)screen[6].x, (int)screen[6].y, (int)screen[7].x, (int)screen[7].y, color);

	Novice::DrawLine((int)screen[5].x, (int)screen[5].y, (int)screen[1].x, (int)screen[1].y, color);
	Novice::DrawLine((int)screen[5].x, (int)screen[5].y, (int)screen[4].x, (int)screen[4].y, color);
	Novice::DrawLine((int)screen[5].x, (int)screen[5].y, (int)screen[7].x, (int)screen[7].y, color);

	Novice::DrawLine((int)screen[3].x, (int)screen[3].y, (int)screen[1].x, (int)screen[1].y, color);
	Novice::DrawLine((int)screen[3].x, (int)screen[3].y, (int)screen[2].x, (int)screen[2].y, color);
	Novice::DrawLine((int)screen[3].x, (int)screen[3].y, (int)screen[7].x, (int)screen[7].y, color);

}

/////////////////////////////



///////当たり判定////////////

//球・球
bool isCollision(const Sphere& s1, const Sphere& s2)
{
	float dis = s1.radius + s2.radius;

	Vector3 dist = { s1.center.x - s2.center.x ,s1.center.y - s2.center.y,s1.center.z - s2.center.z };

	float distance = sqrtf(powf(dist.x, 2) + powf(dist.y, 2) + powf(dist.z, 2));

	if (distance <= dis)
	{
		return true;
	}
	else
	{
		return false;
	}
}

//球・平面
bool isCollision(const Sphere& sphere, const Plane& plane)
{
	float k = sqrtf(powf(Dot(plane.normal, sphere.center) - plane.distance, 2));

	if (k <= sphere.radius)
	{
		return true;
	}
	else
	{
		return false;
	}
}

//線・平面
bool isCollision(const Segment& segment, const Plane& plane)
{
	//法線と線の内積
	float dot = Dot(plane.normal, segment.diff);

	//平行なので衝突しない
	if (dot == 0.0f)
	{
		return false;
	}

	//tを求める
	float t = (plane.distance - Dot(segment.origin, plane.normal)) / dot;

	//衝突判定
	if (t <= 1 && t >= 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

//線・三角形
bool isCollision(const Segment& segment, const Triangle& triangle)//平面求めてから判定を取るー＞三角形との積をなんやかんや
{

	Vector3 v01 = Subtract(triangle.vertices[1], triangle.vertices[0]);
	Vector3 v12 = Subtract(triangle.vertices[2], triangle.vertices[1]);
	Vector3 v20 = Subtract(triangle.vertices[0], triangle.vertices[2]);

	Vector3 n = Normalize(Cross(v01, v12));

	float d = Dot(triangle.vertices[0], n);

	//法線と線の内積
	float dot = Dot(n, segment.diff);

	float t = (d - Dot(segment.origin, n)) / dot;

	if (!(t <= 1 && t >= 0)) 
	{
		return false; 
	}

	Vector3 p = Add(segment.origin, Multiply(t, segment.diff));

	Vector3 v0p = Subtract(p, triangle.vertices[0]);
	Vector3 v1p = Subtract(p, triangle.vertices[1]);
	Vector3 v2p = Subtract(p, triangle.vertices[2]);

	//各辺を結んだベクトル、頂点、衝突点pを結んだベクトルのクロス積を取る
	Vector3 cross01 = Cross(v01, v1p);
	Vector3 cross12 = Cross(v12, v2p);
	Vector3 cross20 = Cross(v20, v0p);

	//衝突判定
	if (Dot(cross01, n) >= 0.0f &&
		Dot(cross12, n) >= 0.0f && 
		Dot(cross20, n) >= 0.0f) {
		return true;
	}
	else {
		return false;
	}
}

//AABB・AABB
bool isCollision(const AABB& aabb1, const AABB& aabb2)
{
	//衝突判定
	if ((aabb1.min.x <= aabb2.max.x && aabb1.max.x >= aabb2.min.x) &&
		(aabb1.min.y <= aabb2.max.y && aabb1.max.y >= aabb2.min.y) &&
		(aabb1.min.z <= aabb2.max.z && aabb1.max.z >= aabb2.min.z))
	{
		return true;
	}
	else
	{
		return false;
	}
}

//AABB・球
bool isCollision(const AABB& aabb, const Sphere& sphere)
{
	//最接点を求める
	Vector3 closestPoint{
	 std::clamp(sphere.center.x,aabb.min.x,aabb.max.x),
	 std::clamp(sphere.center.y,aabb.min.y,aabb.max.y),
	 std::clamp(sphere.center.z,aabb.min.z,aabb.max.z) };

	//最接点と球の中心との距離を求める
	float distance = Length(closestPoint, sphere.center);

	//距離が半径よりも小さければ衝突
	if (distance < sphere.radius)
	{
		return true;
	}
	else
	{
		return false;
	}

}

//AABB・線
bool isCollision(const AABB& aabb, const Segment segment)
{	
	Vector3 nX = { 1,0,0 };
	Vector3 nY = { 0,1,0 };
	Vector3 nZ = { 0,0,1 };

	float dotX = Dot(nX, segment.diff);
	float dotY = Dot(nY, segment.diff);
	float dotZ = Dot(nZ, segment.diff);

	//Vector3 dot = Add(segment.origin, segment.diff);

	float txMin = (aabb.min.x - segment.origin.x) / dotX;
	float txMax = (aabb.max.x - segment.origin.x) / dotX;

	float tyMin = (aabb.min.y - segment.origin.y) / dotY;
	float tyMax = (aabb.max.y - segment.origin.y) / dotY;

	float tzMin = (aabb.min.z - segment.origin.z) / dotZ;
	float tzMax = (aabb.max.z - segment.origin.z) / dotZ;


	float tNearX = min(txMin, txMax);
	float tFarX = max(txMin, txMax);

	float tNearY = min(tyMin, tyMax);
	float tFarY = max(tyMin, tyMax);

	float tNearZ = min(tzMin, tzMax);
	float tFarZ = max(tzMin, tzMax);

	

	//AABBとの衝突点（貫通点）のtが小さい方
	float tmin = max(max(tNearX, tNearY), tNearZ);

	//AABBとの衝突点（貫通点）のtが大きい方
	float tmax = min(min(tFarX, tFarY), tFarZ);
	


	if (txMax > INFINITY or txMin < -INFINITY or
		tyMax > INFINITY or tyMin < -INFINITY or
		tzMax > INFINITY or tzMin < -INFINITY)
	{
		return false;
	}

	/*if (std::isinf(txMax) or std::isinf(txMin) or
		std::isinf(tyMax) or std::isinf(tyMin) or
		std::isinf(tzMax) or std::isinf(tzMin))
	{
		return false;
	}*/
	

	if (std::isnan(txMax) or std::isnan(txMin) or
		std::isnan(tyMax) or std::isnan(tyMin) or
		std::isnan(tzMax) or std::isnan(tzMin))
	{
		return false;
	}

	//衝突判定
	if (tmin <= tmax)
	{
		if ((tmax <= 1 && tmax >= 0) or (tmin <= 1 && tmin >= 0))
		{
			return true;
		}
		else
		{
			return false;
		}
		
	}
	else
	{
		return false;
	}

}

/////////////////////////////

const char kWindowTitle[] = "LC1A_01_イイオカ_イサミ_MT3_02_07_確認課題";

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, 1280, 720);

	AABB aabb
	{
		.min{-0.5f,-0.5f,-0.5f},
		.max{0.5f,0.5f,0.5f},
	};

	Segment segment{
		.origin{-0.7f,0.3f,0.0f},
		.diff{2.0f,-0.5f,0.0f}
	};

	uint32_t color1 = WHITE;
	uint32_t color2 = WHITE;

	Vector3 cameraTranslate{ 0.0f,1.9f,-6.49f };
	Vector3 cameraRotate{ 0.26f,0.0f,0.0f };

	float cameraSpeed = 0.01f;

	// キー入力結果を受け取る箱
	char keys[256] = { 0 };
	char preKeys[256] = { 0 };

	// ウィンドウの×ボタンが押されるまでループ
	while (Novice::ProcessMessage() == 0) {
		// フレームの開始
		Novice::BeginFrame();

		// キー入力を受け取る
		memcpy(preKeys, keys, 256);
		Novice::GetHitKeyStateAll(keys);

		///
		/// ↓更新処理ここから
		///

		if (keys[DIK_W])
		{
			cameraTranslate.z += cameraSpeed;
		}
		if (keys[DIK_S])
		{
			cameraTranslate.z -= cameraSpeed;
		}
		if (keys[DIK_A])
		{
			cameraTranslate.x -= cameraSpeed;
		}
		if (keys[DIK_D])
		{
			cameraTranslate.x += cameraSpeed;
		}
		if (keys[DIK_Q])
		{
			cameraTranslate.y += cameraSpeed;
		}
		if (keys[DIK_E])
		{
			cameraTranslate.y -= cameraSpeed;
		}

		Matrix4x4 cameraMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, cameraRotate, cameraTranslate);
		Matrix4x4 viewMatrix = Inverse(cameraMatrix);
		Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(0.45f, (float)kWindowWidth / (float)kWindowHeight, 0.1f, 100.0f);
		Matrix4x4 viewProjectionMatrix = Multiply(viewMatrix, projectionMatrix);
		Matrix4x4 viewportMatrix = MakeViewportMatrix(0, 0, (float)kWindowWidth, (float)kWindowHeight, 0.0f, 1.0f);

		aabb.min.x = (std::min)(aabb.min.x, aabb.min.x);
		aabb.min.y = (std::min)(aabb.min.y, aabb.min.y);
		aabb.min.z = (std::min)(aabb.min.z, aabb.min.z);
		
		Vector3 start = Transform(Transform(segment.origin, viewProjectionMatrix), viewportMatrix);
		Vector3 end = Transform(Transform(Add(segment.origin, segment.diff), viewProjectionMatrix), viewportMatrix);


		if (isCollision(aabb,segment))
		{
			color1 = RED;
		}
		else
		{
			color1 = WHITE;
		}


		ImGui::Begin("window");
		ImGui::DragFloat3("CameraRotate", &cameraRotate.x, 0.01f);
		ImGui::DragFloat3("aabb.min", &aabb.min.x, 0.01f);
		ImGui::DragFloat3("aabb.max", &aabb.max.x, 0.01f);
		ImGui::DragFloat3("Segment.Origin", &segment.origin.x, 0.01f);
		ImGui::DragFloat3("Segment.Diff", &segment.diff.x, 0.01f);
		ImGui::End();

		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///

		DrawGrid(viewProjectionMatrix, viewportMatrix);
		DrawAABB(aabb, viewProjectionMatrix, viewportMatrix, color1);
		Novice::DrawLine((int)start.x, (int)start.y, (int)end.x, (int)end.y, color2);

		///
		/// ↑描画処理ここまで
		///

		// フレームの終了
		Novice::EndFrame();

		// ESCキーが押されたらループを抜ける
		if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
			break;
		}
	}

	// ライブラリの終了
	Novice::Finalize();
	return 0;
}
