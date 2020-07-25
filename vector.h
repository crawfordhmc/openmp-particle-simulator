#ifndef VECTOR_H
#define VECTOR_H

#include <cmath>

//Class is used for vector calculations
//vectors are used vor velocities and positions
class Vector{
public:
    Vector(double x, double y):x_(x), y_(y){}

    double getX() const { return x_; }
    double getY() const { return y_; }
    void setX(const double x) { x_ = x; }
    void setY(const double y) { y_ = y; }

    //returns unit vector
    Vector getUnitVec (){
        double length = getLength();
        if(length == 0.){
            return Vector(0,0);
        }
        double x = getX()/length;
        double y = getY()/length;
        return Vector(x, y);
    }

    //returns length sqrt(x^2+y^2)
    double getLength () {
        double v = std::sqrt(getX() * getX() + getY() * getY());
        return v;
    }

    //returns sum of a vector added to current vector
    Vector add(const Vector &helpVec){
        double x = getX()+helpVec.getX();
        double y = getY()+helpVec.getY();
        return Vector(x, y);
    }

    //substracts a vector from current vector
    Vector substract(const Vector &helpVec){
        double x = getX()-helpVec.getX();
        double y = getY()-helpVec.getY();
        return Vector(x, y);
    }

    //return dot product
    double dot(const Vector &helpVec){
        double prod = getX()*helpVec.getX()+getY()*helpVec.getY();
        return prod;
    }


private:
    double x_;
    double y_;

};

//multiplication by scalar
inline Vector operator*(const double scal, const Vector &helpVec) {
    double x = scal*helpVec.getX();
    double y = scal*helpVec.getY();
    return Vector(x, y);
}
#endif // VECTOR_H
