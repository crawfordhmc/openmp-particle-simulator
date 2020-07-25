#ifndef PARTICLE_H
#define PARTICLE_H
#include "vector.h"


//class represent particles in the simulation
class Particle{
public:
    Particle(int N, double x, double y, double vx, double vy):
        N_(N), pos_(Vector(x,y)), v_(Vector(vx, vy)) {}

	void setN(int n) { N_ = n; }
    int getN() {return N_;}

    double getX() {return pos_.getX();}
    double getY() {return pos_.getY();}
    Vector getPos() {return pos_;}

    double getVx() {return v_.getX();}
    double getVy() {return v_.getY();}
    Vector getV() {return v_;}

    void setX(const double x) {pos_.setX(x);}
    void setY(const double y) {pos_.setY(y);}
    void setPos(const Vector &helpVec) {pos_=helpVec;}

    void setVx(const double vx) {v_.setX(vx);}
    void setVy(const double vy) {v_.setY(vy);}
    void setV(const Vector &v) {v_=v;}

    //moves particle acording to it's velocity for a given time
    void moveParticle (double t){
        double x = getX() + getVx()*t;
        double y = getY() + getVy()*t;
        setX(x);
        setY(y);
    }

    void addCollisionWall (){collisionsWall_++;}
    int getCollisionsWall (){return collisionsWall_;}

    void addCollisionPart (){collisionsPart_++;}
    int getCollisionsPart (){return collisionsPart_;}

	void setCollided(bool value) { collided_ = value; }
    bool getCollided (){return collided_;}

	void setCollidedWith(int index) { collidedWith_ = index; }
	int getCollidedWith() { return collidedWith_; }

    void setStopPart (bool value){stopPart_ = value;}
    bool getStopPart() {return stopPart_;}

private:
    //amount of collisions with walls
    int collisionsWall_ = 0;
    //amount of collisions with other particles
    int collisionsPart_ = 0;
    //truth value whether particle has collided on this step
    //true= it has, false=it hasn't
    bool collided_ = false;
	//index of particle collided with if applicable
	int collidedWith_ = -1;
    //truthvalue that tells if particle is stopped at a wall true = stopped
    bool stopPart_ = false;

    //index of particle
    int N_;
    //position (x,y)
    Vector pos_;
    //velocity (vx,vy)
    Vector v_;
};

#endif // PARTICLE_H
