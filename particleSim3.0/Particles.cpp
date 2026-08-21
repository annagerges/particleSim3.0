//cpp file for the program's functions

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <fstream>
#include<string>
#include "Particles.h"

using namespace std;

const float dt = 0.001f;

//updates position using rk4 approximation
void updatePos(vector<Particles>& particles, Spring& s) {
	for (auto& p : particles) {
		// Determine if particle is on spring ONCE at start of step
		bool onSpring = (p.getY() <= s.getHeight());

		// RK4 stages
		float k1y = p.getVy();

		//if it's on the spring change it's accelaration to k*dy/m-g and if not than keep accelaration to -9.8
		float k1v = onSpring ? (s.getK() / p.getMass()) * (s.getHeight() - p.getY()) - 9.8f : -9.8f;

		float k2y = p.getVy() + 0.5f * k1v * dt;
		float k2v = onSpring ? (s.getK() / p.getMass()) * (s.getHeight() - (p.getY() + 0.5f * k1y * dt)) - 9.8f : -9.8f;

		float k3y = p.getVy() + 0.5f * k2v * dt;
		float k3v = onSpring ? (s.getK() / p.getMass()) * (s.getHeight() - (p.getY() + 0.5f * k2y * dt)) - 9.8f : -9.8f;

		float k4y = p.getVy() + k3v * dt;
		float k4v = onSpring ? (s.getK() / p.getMass()) * (s.getHeight() - (p.getY() + k3y * dt)) - 9.8f : -9.8f;

		// Update
		p.setY(p.getY() + (dt / 6.0f) * (k1y + 2 * k2y + 2 * k3y + k4y));
		p.setVy(p.getVy() + (dt / 6.0f) * (k1v + 2 * k2v + 2 * k3v + k4v));

	}
}


//checks if any of the particles collided with the wall
void  wallCollis(vector<Particles>& part) {
	for (int index = 0; index < part.size(); index++) {
		//if the particle is beyond or colliding with the left edge: reverse its x velocity and change its x position to 1
		if (part[index].getX() <= 0) {
			part[index].setVx(abs(part[index].getVx()));
			part[index].setX(1);

		}

		//if the particle is beyond or colliding with the right edge: reverse its x velocity and change its x position to 799
		else if (part[index].getX() >= 800) {
			part[index].setVx(abs(part[index].getVx()) * -1);
			part[index].setX(799);
		}

		//if the particle is beyond or colliding with the upper edge: reverse its y velocity and change its y position to 799
		if (part[index].getY() >= 800) {
			part[index].setVy(abs(part[index].getVy()) * -1);
			part[index].setY(799);
		}
		// if the particle is beyond or colliding with the floor: reverse its y velocity and clamp it
		else if (part[index].getY() <= 1) {
			part[index].setVy(abs(part[index].getVy()));
			part[index].setY(1);
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
			dy = grid[index]->getY() - grid[start]->getY();
			dx = grid[index]->getX() - grid[start]->getX();

			//absolute distance
			float absDx = abs(dx);
			float absDy = abs(dy);

			//if distance^2 is less than the particle radius^2: then they collided. Using squared to budget cpu resources and be accurate at the same time
			if (absDx * absDx + absDy * absDy < 100) {

				// Calculate relative velocity
				float dvx = grid[index]->getVx() - grid[start]->getVx();
				float dvy = grid[index]->getVy() - grid[start]->getVy();

				//only swap velo if they are moving towards eachother. Determines if they are pointing to eachother and acts accordingly
				if (dx * dvx + dy * dvy < 0) {

					if (absDx <= absDy) {
						float tempVy = grid[start]->getVy();
						grid[start]->setVy(grid[index]->getVy());
						grid[index]->setVy(tempVy);
					}
					else {
						float tempVx = grid[start]->getVx();
						grid[start]->setVx(grid[index]->getVx());
						grid[index]->setVx(tempVx);
					}
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
//checks accelaration at different points
float checkAccel(float h, Spring& s, float m) {
	if (h <= s.getHeight()){
		return (s.getK() / m) * (s.getHeight() - h) - 9.8f;
	}
	else {
		return -9.8f;
	}
}
