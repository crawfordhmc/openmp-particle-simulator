#include "engine.h"
//#include <omp.h>

Engine::Engine(vector<Particle> particles,
	int N, int L, int r, int S, bool printing) :
	particles_(particles), N_(N), L_(L), r_(r), S_(S), printing_(printing)
{

	//print initial positions
	this->print();

	//go through every step of simulation
	while (curS_ < S) {

		if (curS_ != 0) {

			//Set collided to false and release stopped particles
			for (auto& part : particles_) {
				part.setCollided(false);
				part.setStopPart(false);
				part.setCollidedWith(-1);
			}

			//we only print if it is told in input file
			if (printing_)
				this->print();
		}
		advanceStep();
		curS_++;
	}
	//print final positions
	for (auto& part : particles_) {
		printf("%d %d %10.8f %10.8f %10.8f %10.8f %d %d\n", curS_, part.getN(),
			part.getX(), part.getY(), part.getVx(), part.getVy(), part.getCollisionsPart(), part.getCollisionsWall());
	}
}

Engine::~Engine() {

}

//Prints info of particle
void Engine::print() {
	for (auto& part : particles_) {
		//step id x y vx vy
		printf("%d %d %10.8f %10.8f %10.8f %10.8f\n", curS_, part.getN(), part.getX(), part.getY(), part.getVx(), part.getVy());
	}
}

// function that returns true if collision 1 is earlier
bool Engine::earlier(Collision c1, Collision c2) {
	return c1.t < c2.t;
}

void Engine::advanceStep()
{
	collisions.clear();
	forward_list<Collision> collided;
	Collision stop; //empty collision for stopping particles at walls
	allCollisions(collisions, 1.); //inserts into a list of collisions
	collisions.sort(earlier);

	//sequential checking of valid collisions
	for (Collision c : collisions) { //assuming id corresponds with index like in testcases
		if (c.type != PARTCOLL) { //wall
			if (!particles_[c.partId1].getCollided()) {
				particles_[c.partId1].setCollided(true);
				collided.push_front(c);
			} // else ignore collision
		} //collisions must be particles
		else if (particles_[c.partId1].getCollided() || particles_[c.partId2].getCollided()) {
			continue; //ignore collision
		}
		else { //valid collision is with another particle
			particles_[c.partId1].setCollided(true);
			particles_[c.partId2].setCollided(true);
			collided.push_front(c);
		}
	}

	//#pragma omp parallel for
	for(Collision c: collided){
		particles_[c.partId1].moveParticle(c.t); //move to collision time
		if(c.type != PARTCOLL)
			collideWall(particles_[c.partId1], c.type); //change velocity & increment count
		if(c.type == PARTCOLL){
			particles_[c.partId2].moveParticle(c.t); //move to collision time
			collideParticles(particles_[c.partId1], particles_[c.partId2]); //change velocity & increment count
			stop = getWallCol(1. - c.t, particles_[c.partId2]); // get stopping wall
			if(stop.type != NOCOLL){
				particles_[c.partId2].moveParticle(stop.t); //move additional time to stop
			}
			else {
				particles_[c.partId2].moveParticle(1 - c.t); //move to end of step
			}
			particles_[c.partId2].setStopPart(true);
		}
		stop = getWallCol(1. - c.t, particles_[c.partId1]); // get stopping wall
		if (stop.type != NOCOLL) {
			particles_[c.partId1].moveParticle(stop.t); //move additional time to stop
		}
		else {
			particles_[c.partId1].moveParticle(1. - c.t);
		}
		particles_[c.partId1].setStopPart(true);
	}
	//move non-collided particles
	advanceParticles(1);
}

void Engine::allCollisions(forward_list<Collision> &collisions, double t)
{
	//#pragma omp parallel for
	for (auto& p : particles_) {
		Collision wall = getWallCol(t, p);
		if (wall.type != NOCOLL) {
			//#pragma omp critical (push2list)
			collisions.push_front(wall);
		}
		findPartCols(t, p, collisions);
	}
}

void Engine::advanceParticles(double t)
{
	//#pragma omp parallel for
	for (int p = 0; p < N_; p++) {
		//If particle isn't already collided
		if (!particles_[p].getStopPart()) {
			//move particle
			particles_[p].moveParticle(t);
		}
	}
}

void Engine::findPartCols(double t, Particle& part1, forward_list<Collision> &collisions)
{
	Collision helpCol;
	int helpN = NO_VALUE;
	int p2;

	//#pragma omp parallel for
	for (p2 = 0; p2 < N_; p2++) {
		//If particles are distinct and collision is not a repeat, get collision
		if (p2 != part1.getN() && particles_[p2].getCollidedWith() != part1.getN()) {
			auto& part2 = particles_[p2];
			helpCol = particleCollisionTime(part1, part2, t, r_);
			//Check if we had a collision
			if (helpCol.type == PARTCOLL) {
				part1.setCollidedWith(p2);
				//#pragma omp critical (push2list)
				collisions.push_front(helpCol);
			}
		}
	}
}


//for a particle check the time it takes to collide with a wall,
//return the earliest collision within t, untyped collision if none
Collision Engine::getWallCol(double t, Particle& part)
{
	Collision retCol;
	retCol.type = NOCOLL;
	retCol.t = NO_VALUE;
	retCol.partId1 = part.getN();
	//casting integers to doubles before subtracting other doubles for accuracy
	double Lr = L_ - r_;
	double L = L_;

		//if vx is negative we can have a collision w LEFT wall
	if (part.getVx() < 0) {
		double t_help = (Lr - (L - part.getX())) / fabs(part.getVx());
		//time needs to be positive, smaller than the given time to hit the other walls
		if (t_help >= 0. && t_help < t) {
			retCol.t = t_help;
			retCol.type = LEFTWALL;
		}
	}
	else if (part.getVx() > 0) { //if vx is positive nonzero we can have a collision w Right wall
		double t_help = (Lr - part.getX()) / part.getVx();
		//time needs to be positive, smaller than the given time to hit the other walls
		if (t_help >= 0. && t_help < t) {
			retCol.t = t_help;
			retCol.type = RIGHTWALL;
		}
	}
	if (part.getVy() > 0) {
		//if (comp(part.getVy(), 0) == 0) { //if vy is positive nonzero we can have a collision w top wall
		double t_help = (Lr - part.getY()) / part.getVy();
		//time needs to be positive, smaller than the given time and smaller
		//than lowest time to hit the other walls
		if (t_help >= 0. && t_help < t) {
			if (t_help < retCol.t) { //new time < old
				retCol.t = t_help;
				retCol.type = TOPWALL;
			}
			else if (t_help == retCol.t) { //collide with 2 walls at same time
				retCol.type = CORNER;
			}
		}
	}
	else if (part.getVy() < 0) {
		//else if (comp(part.getVy(), 0) == 1) { //if vy is negative we can have a collision w Bottom wall
		double t_help = (Lr - (L - part.getY())) / fabs(part.getVy());
		//time needs to be positive, smaller than the given time and smaller
		//than lowest time to hit the other walls
		if (t_help >= 0. && t_help < t) {
			if (t_help < retCol.t) { //new time < old
				retCol.t = t_help;
				retCol.type = BOTTOMWALL;
			}
			else if (t_help == retCol.t) { //collide with 2 walls at same time
				retCol.type = CORNER;
			}
		}
	}
	return retCol;
}

//collides given particles -> changes velocities to those after collision
void collideParticles(Particle& part1, Particle& part2)
{
	//projecting the velocity unit normal and unit tangent vectors onto the unit
	//normal and unit tangent vectors, which is done by taking the dot product
	Vector vn = part2.getPos().substract(part1.getPos());
	Vector vun = vn.getUnitVec();
	Vector vut(-vun.getY(), vun.getX());
	double v1n = vun.dot(part1.getV());
	double v2n = vun.dot(part2.getV());
	double v1t = vut.dot(part1.getV());
	double v2t = vut.dot(part2.getV());

	//find the new tangential velocities (after the collision).
	double v1tPrime = v1t;
	double v2tPrime = v2t;

	//find the new normal velocities.
	//this is where we use the one-dimensional collision formulas.
	//for m1 = m2 we get
	double v1nPrime = v2n;
	double v2nPrime = v1n;

	//convert the scalar normal and tangential velocities into vectors.
	//multiply the unit normal vector by the scalar normal velocity
	Vector v_v1nPrime = v1nPrime * vun;
	Vector v_v2nPrime = v2nPrime * vun;// Multiplication by a scalar
	Vector v_v1tPrime = v1tPrime * vut;
	Vector v_v2tPrime = v2tPrime * vut;

	//the final velocity vectors by adding the normal and tangential
	//components for each object
	part1.setVx(v_v1nPrime.getX() + v_v1tPrime.getX());
	part1.setVy(v_v1nPrime.getY() + v_v1tPrime.getY());
	part2.setVx(v_v2nPrime.getX() + v_v2tPrime.getX());
	part2.setVy(v_v2nPrime.getY() + v_v2tPrime.getY());

	//info that collided in cycle and up the collision count
	part1.addCollisionPart();
	part2.addCollisionPart();
}

//collides given particle with the given wall (changes the velocity)
void collideWall(Particle& part, CollType col)
{
	part.addCollisionWall();

	//depending on which wall is hit, either vx or vy will
	//change sign
	if (col == LEFTWALL || col == RIGHTWALL || col == CORNER) {
		part.setVx(-part.getVx());
	}
	if (col == TOPWALL || col == BOTTOMWALL || col == CORNER) {
		part.setVy(-part.getVy());
	}
}

//finds time it takes for given particles to collide with each other
//if it happens within t, collision info is returned
Collision particleCollisionTime(Particle& part1, Particle& part2, double t, int r)
{
	Collision retCol;
	retCol.type = NOCOLL;
	retCol.t = NO_VALUE;

	//Calculate t for collision of particles
	//a=(v2x-v1x)^2+(v2y-v1y)^2
	double a = pow(part2.getVx() - part1.getVx(), 2.0) + pow(part2.getVy() - part1.getVy(), 2.0);
	//b=2*((x20-x10)*(v2x-v1x)+(y20-y10)*(v2y-v1y))
	double b = 2.0 * ((part2.getX() - part1.getX()) * (part2.getVx() - part1.getVx())
		+ (part2.getY() - part1.getY()) * (part2.getVy() - part1.getVy()));
	//c=(x20-x10)^2+(y20-y10)^2-(r1+r2)^ 2
	double c = pow(part2.getX() - part1.getX(), 2.0)
		+ pow(part2.getY() - part1.getY(), 2.0)
		- pow(2 * r, 2.0);

	//Det=b^2-4ac
	double det = pow(b, 2.0) - 4 * a * c;

	// If a==0 here will be no collision
	if (a != 0.) {
		double t_help = (-b - sqrt(det)) / (2.0 * a);
		//collision happens within given t & time is positive
		if (t_help >= 0. && t_help < t) {
			retCol.type = PARTCOLL;
			retCol.t = t_help;
			retCol.partId1 = part1.getN();
			retCol.partId2 = part2.getN();
		}
	}
	return retCol;
}