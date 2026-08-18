//cpp file for the program's functions

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <fstream>
#include<string>
#include "Particles.h"

using namespace std;

const float dt = 0.1f;

//updates position using euclidean approximation
void updatePos(vector<Particles>& part, Spring& s) {
	float currentY, currentVy, k1Y, k1Vy, testY2, testVy2, k2Y, k2Vy, testY3, testVy3, k3Y, k3Vy, testY4, testVy4, k4Y, k4Vy, finalY, finalVy;

	for (int index = 0; index < part.size(); index++) {
		//current Y and velocity
		currentY = part[index].getY();
		currentVy = part[index].getVy();

		//Slope of postition and velosity respectively
		k1Y = part[index].getVy();
		k1Vy = checkAccel(currentY,s,part[index].getMass());

		//predicted y and velocity in 0.05 seconds
		testY2 = currentY + k1Y * (dt / 2);
		testVy2 = currentVy + k1Vy * (dt / 2);

		//Slope of postition and velosity respectively
		k2Y = testVy2;
		k2Vy= checkAccel(testY2, s, part[index].getMass());

		//predicted y and velocity in 0.05 seconds
		testY3 = testY2 + k2Y * (dt/2);
		testVy3 = testVy2 + k2Vy * (dt/2);

		//Slope of postition and velosity respectively
		k3Y = testVy3;
		k3Vy= checkAccel(testY3, s, part[index].getMass());

		//predicted y and velocity in 0.1 seconds
		testY4 = testY3 + k3Y * dt;
		testVy4 = testVy3 + k3Vy * dt;

		k4Y = testVy4;
		k4Vy = checkAccel(testY4, s, part[index].getMass());

		//weighted avg of y and vy to find final values (the midpoint values are twice as accurate as the values at begining and end of time interval)
		finalY = currentY + (dt / 6) * (k1Y+(2*k2Y)+(2*k3Y)+k4Y);
		finalVy = currentVy + (dt / 6) * (k1Vy + (2 * k2Vy) + (2 * k3Vy) + k4Vy);

		//accealration is constant so euler approximation is accurate enough
		part[index].setX(part[index].getX() + part[index].getVx() * dt);

		//set height and vy to final values
		part[index].setY(finalY);
		part[index].setVy(finalVy);
		
	}
}

//checks if any of the particles collided with the wall
void  wallCollis(vector<Particles>& part) {
	for (int index = 0; index < part.size(); index++) {
		//if the particle is beyond or colliding with the left edge: reverse its x velocity and change its x position to 1
		if (part[index].getX() <= 0) {
			part[index].setVx(part[index].getVx() * -1);
			part[index].setX(1);

			//if the particle's vx is under 30: increase its x velocity 0.1
			if (part[index].getVx() < 30) {
				part[index].setVx(part[index].getVx() + 0.1);
			}
		}

		//if the particle is beyond or colliding with the right edge: reverse its x velocity and change its x position to 799
		else if (part[index].getX() >= 800) {
			part[index].setVx(part[index].getVx() * -1);
			part[index].setX(799);
		}

		//if the particle is beyond or colliding with the upper edge: reverse its y velocity and change its y position to 799
		if (part[index].getY() >= 800) {
			part[index].setVy(part[index].getVy() * -1);
			part[index].setY(799);
		}
	}
}

//clears and places the updated particles into their correct grid positions
void clearAndFix(vector<Particles>& part, vector<vector<vector<Particles*>>>& grid, int width) {
	int row, col;

	//assign every particle a grid cell based on its position
	for (int index = 0; index < part.size(); index++) {
		row = part[index].getY() / width;
		col = part[index].getX() / width;

		int nR = grid.size(), nC = grid[0].size();

		//safeguards to make sure that every particle gets assigned a valid cell 

		if (row < 0) {
			row = 0;
		}
		else if (row >= nR) {
			row = nR - 1;
		}
		if (col < 0) {
			col = 0;
		}
		else if (col >= nC) {
			col = nC - 1;
		}

		int r = part[index].getRow(), c = part[index].getCol();

		//if the particle's row or column is different from where it was before in the grid, then it needs to be moved to it's correct cell coordinates
		if (r != row || c != col) {

			//for readability and to loop though the grid cell faster use a reference for grid[r][c] instead of recalculating it multiple times 
			vector<Particles*>&oldCell = grid[r][c];

			for (int i = 0; i < oldCell.size(); i++) {
				if (oldCell[i] == &part[index]) {
					oldCell[i] = oldCell.back();
					oldCell.pop_back();
					break;
				}
			}

			part[index].setRow(row);
			part[index].setCol(col);

			grid[row][col].emplace_back(&part[index]);

		}

	}
}

//checks if particles collided with each other by passing an active grid cell and looping through its particles
void particleCollis(vector<Particles*>& grid) {
	float d, dx, dy;

	for (int start = 0; start < grid.size() - 1; start++) {
		for (int index = start + 1; index < grid.size(); index++) {

			//compute dy and dx to figure out if they collided on x or y axis
			dy = abs(grid[index]->getY() - grid[start]->getY());
			dx = abs(grid[index]->getX() - grid[start]->getX());

			//if distance^2 is less than the particle radius^2: then they collided. Using squared to budget cpu resources and be accurate at the same time
			if (dx * dx + dy * dy < 100) {

				//if they collided on y axis
				if (dx >= dy) {
					//if the second particle is higher then the first one
					if (grid[index]->getY() >= grid[start]->getY()) {

						//increase the lower particle's vy by 0.1 if its vy is less than 30, reverse vy of the higher one and slow it down on y axis by 0.1
						grid[start]->setVy(grid[start]->getVy() + (grid[start]->getVy() < 30) ? 0.1 : 0);
						grid[index]->setVy(grid[index]->getVy() * -1);
						grid[index]->setVy(grid[index]->getVy() - 0.1);
					}
					else {
						//increase the lower particle's vy by 0.1 if its vy is less than 30, reverse vy of the higher one and slow it down on y axis by 0.1
						grid[index]->setVy(grid[index]->getVy() + (grid[index]->getVy() < 30) ? 0.1 : 0);
						grid[start]->setVy(grid[start]->getVy() * -1);
						grid[start]->setVy(grid[start]->getVy() - 0.1);
					}
				}
				// if it collided on the x axis: reverse the vx of both particles
				else {
					grid[start]->setVx(grid[start]->getVx() * -1);
					grid[index]->setVx(grid[index]->getVx() * -1);
				}
			}
		}
	}
}

//saves particle state into csv file
void csvDump(vector<Particles>&part,std::fstream& file,float time) {
	for (int index = 0; index < part.size(); index++) {
		file << to_string(index + 1) <<"," << to_string(part[index].getX()) <<"," << to_string(part[index].getY()) <<"," << to_string(part[index].getVx())
			<< "," << to_string(part[index].getVy()) << "," << to_string(part[index].getA()) << "," << "," << to_string(part[index].getRow()) << "," << to_string(part[index].getCol()) <<", "<<to_string(time)<<"\n";
	}
}

//checks accelaration at different points
float checkAccel(float h, Spring& s, float m) {
	if (h <= s.getHeight()) {
		return ((s.getK() / m) * -1) * (h- s.getHeight()) - 9.8;
	}
	else {
		return -9.8f;
	}
}