#ifndef VECTOR2_H
#define  VECTOR2_H

#include <cmath>
#include <iostream>
#include <algorithm>

namespace Kratos
{
    class Vector2 {
    public:
        union {
            struct {
                double x;
                double y;
            };
            double array[2];
        };
    public:
        constexpr Vector2(void) : x(0.0), y(0.0) {}

        constexpr Vector2(double xVal, double yVal) : x(xVal), y(yVal) {}

        ~Vector2(void) {}

        Vector2 Normalised() const {
            Vector2 temp(x, y);
            temp.Normalise();
            return temp;
        }

        void			Normalise() {
            double length = Length();

            if (length != 0.0) {
                length = 1.0 / length;
                x = x * length;
                y = y * length;
            }
        }

        double	Length() const {
            return sqrt((x*x) + (y*y));
        }

        constexpr double	LengthSquared() const {
            return ((x*x) + (y*y));
        }

        constexpr double		GetMaxElement() const {
            double v = x;
            if (y > v) {
                v = y;
            }
            return v;
        }

        double		GetAbsMaxElement() const {
            double v = abs(x);
            if (abs(y) > v) {
                v = abs(y);
            }
            return v;
        }

        static constexpr double	Dot(const Vector2 &a, const Vector2 &b) {
            return (a.x*b.x) + (a.y*b.y);
        }

        static double Cross(const Vector2 &a, const Vector2 &b) {
            return (a.x * b.y - a.y * b.x);
        }

        inline Vector2  operator+(const Vector2  &a) const {
            return Vector2(x + a.x, y + a.y);
        }

        inline Vector2  operator-(const Vector2  &a) const {
            return Vector2(x - a.x, y - a.y);
        }

        inline Vector2  operator-() const {
            return Vector2(-x, -y);
        }

        inline Vector2  operator*(double a)	const {
            return Vector2(x * a, y * a);
        }

        inline Vector2  operator*(const Vector2  &a) const {
            return Vector2(x * a.x, y * a.y);
        }

        inline Vector2  operator/(const Vector2  &a) const {
            return Vector2(x / a.x, y / a.y);
        };

        inline Vector2  operator/(double v) const {
            return Vector2(x / v, y / v);
        };

        inline constexpr void operator+=(const Vector2  &a) {
            x += a.x;
            y += a.y;
        }

        inline void operator-=(const Vector2  &a) {
            x -= a.x;
            y -= a.y;
        }


        inline void operator*=(const Vector2  &a) {
            x *= a.x;
            y *= a.y;
        }

        inline void operator/=(const Vector2  &a) {
            x /= a.x;
            y /= a.y;
        }

        inline void operator*=(double f) {
            x *= f;
            y *= f;
        }

        inline void operator/=(double f) {
            x /= f;
            y /= f;
        }

        inline double operator[](int i) const {
            return array[i];
        }

        inline double& operator[](int i) {
            return array[i];
        }

        inline bool	operator==(const Vector2 &A)const { return (A.x == x && A.y == y) ? true : false; };
        inline bool	operator!=(const Vector2 &A)const { return (A.x == x && A.y == y) ? false : true; };

        inline friend std::ostream& operator<<(std::ostream& o, const Vector2& v) {
            o << "Vector2(" << v.x << "," << v.y << ")" << std::endl;
            return o;
        }
    };
}

#endif // VECTOR2_H