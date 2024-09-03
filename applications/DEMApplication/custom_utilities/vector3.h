#ifndef VECTOR3_H
#define  VECTOR3_H

#include <cmath>
#include <iostream>
#include <algorithm>

namespace Kratos
{
    class Vector3 {
    public:
        union {
            struct {
                double x;
                double y;
                double z;
            };
            double array[3];
        };
    public:
        constexpr Vector3(void) : x(0.0), y(0.0), z(0.0) {}

        constexpr Vector3(double xVal, double yVal, double zVal) : x(xVal), y(yVal), z(zVal) {}

        ~Vector3(void) {}

        Vector3 Normalised() const {
            Vector3 temp(x, y, z);
            temp.Normalise();
            return temp;
        }

        void			Normalise() {
            double length = Length();

            if (length != 0.0) {
                length = 1.0 / length;
                x = x * length;
                y = y * length;
                z = z * length;
            }
        }

        double	Length() const {
            return sqrt((x*x) + (y*y) + (z*z));
        }

        constexpr double	LengthSquared() const {
            return ((x*x) + (y*y) + (z*z));
        }

        constexpr double		GetMaxElement() const {
            double v = x;
            if (y > v) {
                v = y;
            }
            if (z > v) {
                v = z;
            }
            return v;
        }

        double		GetAbsMaxElement() const {
            double v = abs(x);
            if (abs(y) > v) {
                v = abs(y);
            }
            if (abs(z) > v) {
                v = abs(z);
            }
            return v;
        }

        static constexpr double	Dot(const Vector3 &a, const Vector3 &b) {
            return (a.x*b.x) + (a.y*b.y) + (a.z*b.z);
        }

        static Vector3	Cross(const Vector3 &a, const Vector3 &b) {
            return Vector3((a.y*b.z) - (a.z*b.y), (a.z*b.x) - (a.x*b.z), (a.x*b.y) - (a.y*b.x));
        }

        inline Vector3  operator+(const Vector3  &a) const {
            return Vector3(x + a.x, y + a.y, z + a.z);
        }

        inline Vector3  operator-(const Vector3  &a) const {
            return Vector3(x - a.x, y - a.y, z - a.z);
        }

        inline Vector3  operator-() const {
            return Vector3(-x, -y, -z);
        }

        inline Vector3  operator*(double a)	const {
            return Vector3(x * a, y * a, z * a);
        }

        inline Vector3  operator*(const Vector3  &a) const {
            return Vector3(x * a.x, y * a.y, z * a.z);
        }

        inline Vector3  operator/(const Vector3  &a) const {
            return Vector3(x / a.x, y / a.y, z / a.z);
        };

        inline Vector3  operator/(double v) const {
            return Vector3(x / v, y / v, z / v);
        };

        inline constexpr void operator+=(const Vector3  &a) {
            x += a.x;
            y += a.y;
            z += a.z;
        }

        inline void operator-=(const Vector3  &a) {
            x -= a.x;
            y -= a.y;
            z -= a.z;
        }


        inline void operator*=(const Vector3  &a) {
            x *= a.x;
            y *= a.y;
            z *= a.z;
        }

        inline void operator/=(const Vector3  &a) {
            x /= a.x;
            y /= a.y;
            z /= a.z;
        }

        inline void operator*=(double f) {
            x *= f;
            y *= f;
            z *= f;
        }

        inline void operator/=(double f) {
            x /= f;
            y /= f;
            z /= f;
        }

        inline double operator[](int i) const {
            return array[i];
        }

        inline double& operator[](int i) {
            return array[i];
        }

        inline bool	operator==(const Vector3 &A)const { return (A.x == x && A.y == y && A.z == z) ? true : false; };
        inline bool	operator!=(const Vector3 &A)const { return (A.x == x && A.y == y && A.z == z) ? false : true; };

        inline friend std::ostream& operator<<(std::ostream& o, const Vector3& v) {
            o << "Vector3(" << v.x << "," << v.y << "," << v.z << ")" << std::endl;
            return o;
        }
    };
}

#endif // VECTOR3_H