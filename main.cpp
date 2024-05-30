#include <Novice.h>
#include <imgui.h>
#include <cmath>
#include "Matrix4x4.h"
#include "assert.h"
#include "Vector3.h"

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
		1.0f, 0.0f, 0.0f, 0.0f, 0.0f, std::cos(radian), std::sin(radian), 0.0f, 0.0f, -std::sin(radian), std::cos(radian), 0.0f, 0.0f, 0.0f, 0.0f, 1.0f,
	};

	return result;
}

// Y軸回転行列
Matrix4x4 MakeRotateYMatrix(float radian) 
{
	Matrix4x4 result;
	result = {
		std::cos(radian), 0.0f, -std::sin(radian), 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, std::sin(radian), 0.0f, std::cos(radian), 0.0f, 0.0f, 0.0f, 0.0f, 1.0f,
	};

	return result;
}

// Z軸回転行列
Matrix4x4 MakeRotateZMatrix(float radian)
{
	Matrix4x4 result;
	result = {
		std::cos(radian), std::sin(radian), 0.0f, 0.0f, -std::sin(radian), std::cos(radian), 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f,
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


static const int kWindowWidth = 1280;
static const int kWindowHeight = 720;

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
	uint32_t color;//色
};

struct Triangle
{
	Vector3 vertices[3]; //!< 頂点
};

void DrawTriangle(const Triangle& triangle, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color)
{
	Matrix4x4 WorldMatrix1 = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, triangle.vertices[0]);
	Matrix4x4 WorldMatrix2 = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, triangle.vertices[1]);
	Matrix4x4 WorldMatrix3 = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, triangle.vertices[2]);

	Matrix4x4 wvpMatrix1 = Multiply(WorldMatrix1, viewProjectionMatrix);
	Matrix4x4 wvpMatrix2 = Multiply(WorldMatrix2, viewProjectionMatrix);
	Matrix4x4 wvpMatrix3 = Multiply(WorldMatrix3, viewProjectionMatrix);

	Vector3 Local1 = Transform(triangle.vertices[0], wvpMatrix1);
	Vector3 Local2 = Transform(triangle.vertices[1], wvpMatrix2);
	Vector3 Local3 = Transform(triangle.vertices[2], wvpMatrix3);

	Vector3 Screen1 = Transform(Local1, viewportMatrix);
	Vector3 Screen2 = Transform(Local2, viewportMatrix);
	Vector3 Screen3 = Transform(Local3, viewportMatrix);

	//points をそれぞれ結んでDraw で矩形を描画する。DrawTringleを使って塗りつぶしても良いが、DepthがないのでMT3では分かりずらい
	Novice::DrawTriangle((int)Screen1.x, (int)Screen1.y, (int)Screen2.x, (int)Screen2.y, (int)Screen3.x, (int)Screen3.y, color, kFillModeWireFrame);
}

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

	//各辺を結んだベクトルと、頂点と衝突点pを結んだベクトルのクロス積を取る
	Vector3 cross01 = Cross(v01, v1p);
	Vector3 cross12 = Cross(v12, v2p);
	Vector3 cross20 = Cross(v20, v0p);

	//衝突判定
	if (Dot(cross01, n) >= 0.0f &&
		Dot(cross12, n) >= 0.0f &&
		Dot(cross20, n) >= 0.0f)
	{
		return true;
	}
	else
	{
		return false;
	}
}

const char kWindowTitle[] = "LC1A_01_イイオカ_イサミ_MT3_02_04_確認課題";

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, 1280, 720);

	Vector3 cameraTranslate{ 0.0f,1.9f,-6.49f };
	Vector3 cameraRotate{ 0.26f,0.0f,0.0f };

	Triangle triangle;
	triangle.vertices[0] = { -0.5f,0.0f,0.0f };
	triangle.vertices[1] = { 0.0f,0.5f,0.0f };
	triangle.vertices[2] = { 0.5f,0.0f,0.0f };

	Segment segment{ { -0.45f,0.41f,0.0f},{1.0f,0.58f,0.0f},WHITE };

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

		ImGui::Begin("window");
		ImGui::DragFloat3("CameraRotate", &cameraRotate.x, 0.01f);
		ImGui::DragFloat3("Triangle.v0", &triangle.vertices[0].x, 0.01f);
		ImGui::DragFloat3("Triangle.v1", &triangle.vertices[1].x, 0.01f);
		ImGui::DragFloat3("Triangle.v2", &triangle.vertices[2].x, 0.01f);
		ImGui::DragFloat3("Segment.Origin", &segment.origin.x, 0.01f);
		ImGui::DragFloat3("Segment.Diff", &segment.diff.x, 0.01f);
		ImGui::End();

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

		Matrix4x4 camelaMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, cameraRotate, cameraTranslate);
		Matrix4x4 viewMatriix = Inverse(camelaMatrix);
		Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(0.45f, (float)kWindowWidth / (float)kWindowHeight, 0.1f, 100.0f);
		Matrix4x4 viewProjectionMatrix = Multiply(viewMatriix, projectionMatrix);
		Matrix4x4 viewportMatrix = MakeViewportMatrix(0, 0, (float)kWindowWidth, (float)kWindowHeight, 0.0f, 1.0f);

		Vector3 start = Transform(Transform(segment.origin, viewProjectionMatrix), viewportMatrix);
		Vector3 end = Transform(Transform(Add(segment.origin, segment.diff), viewProjectionMatrix), viewportMatrix);


		if (isCollision(segment, triangle))
		{
			segment.color = RED;
		}
		else
		{
			segment.color = WHITE;
		}

		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///

		DrawGrid(viewProjectionMatrix, viewportMatrix);
		DrawTriangle(triangle, viewProjectionMatrix, viewportMatrix, WHITE);
		Novice::DrawLine((int)start.x, (int)start.y, (int)end.x, (int)end.y, segment.color);

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
