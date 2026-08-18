#ifndef PARTICLES_H
#define PARTICLES_H

#include<vector>
#include <fstream>

//Particle class to keep track of every particle's position, velocity, acceleration, and energy
class Particles {
private:
	float vx = 0.0f;
	float vy;

	float x;
	float y;

	float a = -9.8f;

	float mass = 0.5f;

	int cellRow;
	int cellCol;

public:

	//setter and getter vx
	void setVx(float num) {
		vx = num;
	}

	float getVx() const {
		return vx;
	}

	//setter and getter vy
	void setVy(float num) {
		vy = num;
	}

	float getVy() const {
		return vy;
	}

	//setter and getter x and y
	void setX(float num) {
		x = num;
	}

	float getX() const {
		return x;
	}

	void setY(float num) {
		y = num;
	}

	float getY() const {
		return y;
	}

	//setter and getter acceleration
	void setA(float num) {
		a = num;
	}

	float getA() const {
		return a;
	}

	//setter and getter mass (for when user can dynamically customize mass)
	void setMass(float num) {
		mass = num;
	}

	float getMass() const {
		return mass;
	}

	//setter and getter for grid coordinates

	void setRow(int r) {
		cellRow = r;
	}

	int getRow() {
		return cellRow;
	}

	void setCol(int c) {
		cellCol = c;
	}

	int getCol() {
		return cellCol;
	}



};


class Spring {
private:
	float k;
	float h = 6;

public:
	float getK() const {
		return k;
	}

	void setK(float num) {
		k = num;
	}

	float getHeight() const {
		return h;
	}

};

//functions
void updatePos(std::vector<Particles>&, Spring&);
void wallCollis(std::vector<Particles>&);
void clearAndFix(std::vector<Particles>&, std::vector<std::vector<std::vector<Particles*>>>&, int);
void particleCollis(std::vector<Particles*>&);
void csvDump(std::vector<Particles>&, std::fstream&,float);
float checkAccel(float, Spring&, float);

#endif


