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
	for (int index = 0; index < part.size(); index++) {
		//changes acceleration because the ball is on the spring
		if (part[index].getY() <= s.getHeight()) {
			part[index].setA((((s.getK() / part[index].getMass()) * -1) * (part[index].getY() - s.getHeight())) - 9.8);
		}
		else {
			part[index].setA(-9.8);
		}
		part[index].setVy(part[index].getVy() + part[index].getA() * 0.01);
		part[index].setY(part[index].getY() + part[index].getVy() * 0.01);
		part[index].setX(part[index].getX() + part[index].getVx() * 0.01);
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
			vector<Particles*>oldCell = grid[r][c];

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
void csvDump(vector<Particles>&part,std::fstream& file) {
	for (int index = 0; index < part.size(); index++) {
		file << to_string(index + 1) <<"," << to_string(part[index].getX()) <<"," << to_string(part[index].getY()) <<"," << to_string(part[index].getVx())
			<< "," << to_string(part[index].getVy()) << "," << to_string(part[index].getA()) << "," << to_string(part[index].getEp()) << ","
			<< to_string(part[index].getEk()) << "," << to_string(part[index].getRow()) << "," << to_string(part[index].getCol()) << "\n";
	}
}