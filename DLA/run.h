#ifndef RUN_H
#define RUN_H
#include <vector>
#include <fstream>
#include "rand.h"
#include "comn.h"
#include "particle.h"
#include "grid.h"

void launch(Disk &p, Ran *myran);
bool one_step(Disk &p0, std::vector<Disk> &cluster, Ran *myran, Grid &grid);
void run(std::vector<Disk> &cluster, int nPar, Ran *myran, Grid &grid);
void cal_fractal_dim(const std::vector<Disk> &cluster);

void launch(Rect &p, Ran *myran);
bool one_step(Rect &p, std::vector<Rect> &cluster, Ran *myran);
void run(std::vector<Rect> &cluster, int nPar, Ran *myran);
#endif