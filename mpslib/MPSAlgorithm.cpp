// (c) 2015-2016 I-GIS (www.i-gis.dk) and Solid Earth Geophysics, Niels Bohr Institute (http://imgp.nbi.ku.dk)
//
//    This file is part of MPSlib.
//
//    MPSlib is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Lesser General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    MPSlib is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public License
//    along with MPSlib (COPYING.LESSER).  If not, see <http://www.gnu.org/licenses/>.
//
#include <iomanip>      // std::setprecision
#include <iostream>
#include <cctype>       // isspace
#include <time.h>
#include <map>

#include "MPSAlgorithm.h"
#include "Coords3D.h"
#include "Utility.h"
#include "IO.h"

/**
* @brief Constructors
*/
MPS::MPSAlgorithm::MPSAlgorithm(void) {
	//Initialize parameters
	_hdSearchRadius = 1;
	_totalGridsLevel = 1;
	_sgDimX = 0;
	_sgDimY = 0;
	_sgDimZ = 0;
	_sgWorldMinX = 0;
	_sgWorldMinY = 0;
	_sgWorldMinZ = 0;
	_sgCellSizeX = 0;
	_sgCellSizeY = 0;
	_sgCellSizeZ = 0;
	_tiDimX = 0;
	_tiDimY = 0;
	_tiDimZ = 0;
	_hasMaskData = false;	
}

/**
* @brief Destructors
*/
MPS::MPSAlgorithm::~MPSAlgorithm(void) {

}

/**
* @brief Initialize the Hard Data Grid with a value, default is NaN
* @param hdg the simulation GRID
* @param sgDimX dimension X of the grid
* @param sgDimY dimension Y of the gri
* @param sgDimZ dimension Z of the grid
* @param value value of each grid node default is NAN
*/
void MPS::MPSAlgorithm::_initializeHDG(std::vector<std::vector<std::vector<float>>>& hdg, const int& sgDimX, const int& sgDimY, const int& sgDimZ, const float& value) {

	hdg.resize(sgDimZ);
	for (int z=0; z<sgDimZ; z++) {
		hdg[z].resize(sgDimY);
		for (int y=0; y<sgDimY; y++) {
			hdg[z][y].resize(sgDimX);
			for (int x=0; x<sgDimX; x++) {
				hdg[z][y][x] = value;
			}
		}
	}
}

/**
* @brief Initialize the Simulation Grid with a value, default is NaN
* @param sg the simulation GRID
* @param sgDimX dimension X of the grid
* @param sgDimY dimension Y of the gri
* @param sgDimZ dimension Z of the grid
* @param value value of each grid node default is NAN
*/
void MPS::MPSAlgorithm::_initializeSG(std::vector<std::vector<std::vector<float>>>& sg, const int& sgDimX, const int& sgDimY, const int& sgDimZ, const float& value) {
	sg.resize(sgDimZ);
	for (int z=0; z<sgDimZ; z++) {
		sg[z].resize(sgDimY);
		for (int y=0; y<sgDimY; y++) {
			sg[z][y].resize(sgDimX);
			for (int x=0; x<sgDimX; x++) {
				sg[z][y][x] = value;
			}
		}
	}
	//std::cout << "SG (X, Y, Z): " << SG[0][0].size() << " " << SG[0].size() << " " << SG.size() << " " << SG[0][0][0]<< std::endl;
}

/**
* @brief Initialize the Simulation Grid from a 3D grid with values
* The Copy and simulation grid must have the same dimension in X, Y and Z
* @param sg the simulation GRID
* @param sgDimX dimension X of the grid
* @param sgDimY dimension Y of the gri
* @param sgDimZ dimension Z of the grid
* @param grid the grid to copy data in
* @param value value in copy grid considered to be a NAN value default is -1
*/
void MPS::MPSAlgorithm::_initializeSG(std::vector<std::vector<std::vector<float>>>& sg, const int& sgDimX, const int& sgDimY, const int& sgDimZ, const std::vector<std::vector<std::vector<float>>>& grid, const float& nanValue) {

	sg.resize(sgDimZ);
	for (int z=0; z<sgDimZ; z++) {
		sg[z].resize(sgDimY);
		for (int y=0; y<sgDimY; y++) {
			sg[z][y].resize(sgDimX);
			for (int x=0; x<sgDimX; x++) {
				if(!MPS::utility::is_nan(grid[z][y][x])) sg[z][y][x] = grid[z][y][x];
				else sg[z][y][x]= std::numeric_limits<float>::quiet_NaN();
				// std::cout << "init:: "  << x << y << z << " VAL " <<  _sg[z][y][x] << " " << grid[z][y][x] << std::endl;
			}
		}
	}
}

/**
* @brief Initilize a sequential simulation path
* @param sgDimX dimension X of the path
* @param sgDimY dimension Y of the path
* @param sgDimZ dimension Z of the path
* @param path output simulation path
*/
void MPS::MPSAlgorithm::_initilizePath(const int& sgDimX, const int& sgDimY, const int& sgDimZ, std::vector<int>& path) {
	//Putting sequential indices
	int cnt = 0;
	for (int z=0; z<sgDimZ; z++) {
		for (int y=0; y<sgDimY; y++) {
			for (int x=0; x<sgDimX; x++) {
				path[cnt] = cnt;
				cnt++;
			}
		}
	}
}


/**
* @brief Generate a realization from a PDF defined as a map
* @param the pdf as a std::map
* @param a realization from the pdf
*/
float MPS::MPSAlgorithm::_sampleFromPdf(std::map<float, float>& Pdf) {
	float simulatedValue;
	float randomValue	;
	randomValue = ((float) rand() / (RAND_MAX));

	float cumsum_pdf=0; // integral conditional probability density (conditionalPdfFromTi)
	//std::cout << "   cPDF=" << Pdf[0] << "," << Pdf[1] << std::endl;
	for(std::map<float,float>::iterator iter = Pdf.begin(); iter != Pdf.end(); ++iter) {
		cumsum_pdf = cumsum_pdf + iter->second;
		if (cumsum_pdf >= randomValue) {
			//std::cout << "   randomValue=" << randomValue << "<=" <<  cumsum_pdf << ", THEN index=" << iter->first << std::endl;
			simulatedValue = iter->first;
			break;
		}
	}

	return simulatedValue;
}

/**
* @brief Check if the current node is closed to a node in a given grid
* The Given grid could be softdata or harddata grid
* @param x coordinate X of the current node
* @param y coordinate Y of the current node
* @param z coordinate Z of the current node
* @param level current grid level
* @param grid given grid
* @param searchRadius search radius for closed point search
* @param closestCoordinates closest coordinates found
* @return True if found a closed node
*/
bool MPS::MPSAlgorithm::_IsClosedToNodeInGrid(const int& x, const int& y, const int& z, const int& level, const std::vector<std::vector<std::vector<float>>>& grid, const float& searchRadius, MPS::Coords3D& closestCoordinates) {
	//Using circular search
	//TODO: Need to check this again, it runs slow with HD
	std::vector<MPS::Coords3D> L;
	std::vector<float> V;
	_circularSearch(x, y, z, grid, 1, searchRadius, L, V);
	//_circularSearch(x, y, z, grid, 1, std::pow(2, level), L, V);
	bool foundClosest = L.size() > 0;
	if (foundClosest) {
		closestCoordinates.setX(x + L[0].getX());
		closestCoordinates.setY(y + L[0].getY());
		closestCoordinates.setZ(z + L[0].getZ());
	}
	return foundClosest;

	////define search area
	//int fromX, toX, fromY, toY, fromZ, toZ;
	///*fromX = (int)((x - (searchRadius)) / _sgCellSizeX);
	//toX = (int)((x + (searchRadius)) / _sgCellSizeX);
	//fromY = (int)((y - (searchRadius)) / _sgCellSizeY);
	//toY = (int)((y + (searchRadius)) / _sgCellSizeY);
	//fromZ = (int)((z - (searchRadius)) / _sgCellSizeZ);
	//toZ = (int)((z + (searchRadius)) / _sgCellSizeZ);*/
	//float levelSearchRadius = std::pow(2, level) *  searchRadius;
	//fromX = (int)(x - (levelSearchRadius));
	//toX = (int)(x + (levelSearchRadius));
	//fromY = (int)(y - (levelSearchRadius));
	//toY = (int)(y + (levelSearchRadius));
	//fromZ = (int)(z - (levelSearchRadius));
	//toZ = (int)(z + (levelSearchRadius));
	////clamp value to avoid out of bounds
	//if (fromX < 0) fromX = 0;
	//if (toX > _sgDimX) toX = _sgDimX;
	//if (fromY < 0) fromY = 0;
	//if (toY > _sgDimY) toY = _sgDimY;
	//if (fromZ < 0) fromZ = 0;
	//if (toZ > _sgDimZ) toZ = _sgDimZ;
	////searching for the closest point TODO: use circular search for better performance
	//float distMin = pow(searchRadius, 2); //avoid the sqrt later for computing the 3d distance
	//float dist = 0;
	//bool foundClosest = false;
	//int closestX = fromX, closestY = fromY, closestZ = fromZ;
	//int step = std::pow(2, level);
	////std::cout << fromX << " " << toX << " " << fromY << " " << toY << " " << fromZ << " " << toZ << " " << step << std::endl;
	//for (int cZ = fromZ; cZ<toZ; cZ+=step) {
	//	for (int cY = fromY; cY<toY; cY+=step) {
	//		for (int cX = fromX; cX<toX; cX+=step) {
	//			//Compute distance and get the closest point
	//			if(!MPS::utility::is_nan(grid[cZ][cY][cX])) {
	//				dist = (pow((cX - x)*_sgCellSizeX, 2) + pow((cY - y)*_sgCellSizeY, 2) + pow((cZ - z)*_sgCellSizeZ, 2));
	//				if (dist <= distMin) {
	//					distMin = dist;
	//					closestX = cX;
	//					closestY = cY;
	//					closestZ = cZ;
	//					foundClosest = true;
	//				}
	//			}
	//		}
	//	}
	//}
	//closestCoordinates.setX(closestX);
	//closestCoordinates.setY(closestY);
	//closestCoordinates.setZ(closestZ);
	//return foundClosest;
}

/**
* @brief Compute cpdf from softdata as map
* @param x coordinate X of the current node
* @param y coordinate Y of the current node
* @param z coordinate Z of the current node
* @param level current grid level
* @param softPdf softPdf list
* @param closestCoords closest point used with relocation, if not then the current x, y, z is used
* @return computed value from softdata
*/
bool MPS::MPSAlgorithm::_getCpdfFromSoftData(const int& x, const int& y, const int& z, const int& level, std::map<float, float>& softPdf, MPS::Coords3D& closestCoords) {

	// clear soft data - just in case - needed for preferential path to work
	softPdf.clear();

	//Empty grids check
	if (_softDataGrids.empty()) return false;
	//Out of bound check
	int sdgDimX = (int)_softDataGrids[0][0][0].size();
	int sdgDimY = (int)_softDataGrids[0][0].size();
	int sdgDimZ = (int)_softDataGrids[0].size();
	if ((x >= sdgDimX) || (y >= sdgDimY) || (z >= sdgDimZ)) return false;
	//Define a default closest node at the current position
	closestCoords.setX(x);
	closestCoords.setY(y);
	closestCoords.setZ(z);
	if (level == 0) { //For coarse level, then check for the same node localtion in softdata grid
		//No value check
		if (MPS::utility::is_nan(_softDataGrids[0][z][y][x])) return false;
	} else {
		//Check if node is closed to a search radius if using multiple level, doing the relocation here
		if (MPS::utility::is_nan(_softDataGrids[0][closestCoords.getZ()][closestCoords.getY()][closestCoords.getX()])) {
			if (!_IsClosedToNodeInGrid(x, y, z, level, _softDataGrids[0], ceil(pow(2, level) / 2), closestCoords)) return false;
		}
		//if (!_IsClosedToNodeInGrid(x, y, z, _softDataGrids[0], _hdSearchRadius, closestCoords)) return false;
		//Check if the closest node found already in the current relocated node, if so then stop
		//if(std::find(addedNodes.begin(), addedNodes.end(), closestCoords) != addedNodes.end()) return false;
	}
	//Perform computation
	//Looping through all the softdata categories grids and fill the conditional points
	//std::multimap<float, float> softPdf;
	float sumProbability = 0;
	unsigned int lastIndex = (int)_softDataCategories.size() - 1;
	for (unsigned int i=0; i<lastIndex; i++) {
		sumProbability += _softDataGrids[i][closestCoords.getZ()][closestCoords.getY()][closestCoords.getX()];
		//std::cout << "i="<< i << " :: "<< _softDataGrids[i][closestCoords.getZ()][closestCoords.getY()][closestCoords.getX()] << std::endl;
		softPdf.insert(std::pair<float, float>(_softDataCategories[i], _softDataGrids[i][closestCoords.getZ()][closestCoords.getY()][closestCoords.getX()]));
	}
	//Last categorie
	softPdf.insert(std::pair<float, float>(_softDataCategories[lastIndex], 1 - sumProbability));

	if (_debugMode>1) {
		std::cout << "_getCpdfFromSoftData->[x,y,z]=" << x << "," << y << "," << z  << " ### closest #### ";
		std::cout << "" << closestCoords.getX() << "," << closestCoords.getY() << "," << closestCoords.getZ() << std::endl;
		std::cout << "_getCpdfFromSoftData->  -  [" << softPdf[0] << ","<< softPdf[1] << "]"<<std::endl;
	}

	return true;
}

/**
* @brief Show the SG in the terminal
*/
void MPS::MPSAlgorithm::_showSG(void) const {
	for (int z=0; z<_sgDimZ; z++) {
		std::cout << "Z: " << (z + 1) << "/" << _sgDimZ << std::endl;
		for (int y=0; y<_sgDimY; y++) {
			for (int x=0; x<_sgDimX; x++) {
				std::cout << MPS::io::onscreenChars[int(_sg[z][y][x]) % MPS::io::onscreenChars.size()];
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
}

/**
* @brief Read different data (TI, hard and softdata from files)
*/
void MPS::MPSAlgorithm::_readDataFromFiles(void) {
	//Reading TI file
	bool readSucessfull = false;
	std::string fileExtension = MPS::utility::getExtension(_tiFilename);

	if (_debugMode>0) {
		std::cout << "READING TI: " << _tiFilename << std::endl;
	}
	// Read TI
	_readTIFromFiles();

	if (_debugMode>0) {
		std::cout << "READING HARD DATA: " << _hardDataFileNames << std::endl;
	}
	//Reading Hard conditional data
	_readHardDataFromFiles();

	if (_debugMode>0) {
		std::cout << "READING SOFT DATA: " << _softDataFileNames[0] << std::endl;
	}
	//Reading Soft conditional data
	_readSoftDataFromFiles();

	if (_debugMode>0) {
		std::cout << "READING MASK: " << _maskDataFileName << std::endl;
	}
	//Reading Mask data
	_readMaskDataFromFile();
}

/**
* @brief Read Hard data
*/
void MPS::MPSAlgorithm::_readHardDataFromFiles(void) {
	bool readSucessfull = false;
	std::string fileExtension = MPS::utility::getExtension(_tiFilename);

	//Reading Hard conditional data
	readSucessfull = false;
	fileExtension = MPS::utility::getExtension(_hardDataFileNames);
	if (fileExtension == "csv" || fileExtension == "txt") readSucessfull = MPS::io::readTIFromGS3DCSVFile(_hardDataFileNames, _hdg);
	else if (fileExtension == "gslib" || fileExtension == "sgems" || fileExtension == "SGEMS") readSucessfull = MPS::io::readTIFromGSLIBFile(_hardDataFileNames, _hdg);
	else if (fileExtension == "dat") readSucessfull = MPS::io::readHardDataFromEASFile(_hardDataFileNames, -999, _sgDimX, _sgDimY, _sgDimZ, _sgWorldMinX, _sgWorldMinY, _sgWorldMinZ, _sgCellSizeX, _sgCellSizeY, _sgCellSizeZ, _hdg);
	else if (fileExtension == "grd3") readSucessfull = MPS::io::readTIFromGS3DGRD3File(_hardDataFileNames, _hdg);
	if ((!readSucessfull)&(_debugMode>-1)) {
		std::cout << "Error reading harddata " << _hardDataFileNames << std::endl;
	}
}

/**
* @brief Read soft data
*/
void MPS::MPSAlgorithm::_readSoftDataFromFiles(void) {
	bool readSucessfull = false;
	std::string fileExtension = MPS::utility::getExtension(_tiFilename);

	for (unsigned int i = 0; i<_softDataFileNames.size(); i++) {
		readSucessfull = false;
		fileExtension = MPS::utility::getExtension(_softDataFileNames[i]);
		if (fileExtension == "csv" || fileExtension == "txt") readSucessfull = MPS::io::readTIFromGS3DCSVFile(_softDataFileNames[i], _softDataGrids[i]);
		else if (fileExtension == "gslib" || fileExtension == "sgems" || fileExtension == "SGEMS") readSucessfull = MPS::io::readTIFromGSLIBFile(_softDataFileNames[i], _softDataGrids[i]);
		else if (fileExtension == "dat") readSucessfull = MPS::io::readSoftDataFromEASFile(_softDataFileNames[i], _softDataCategories, _sgDimX, _sgDimY, _sgDimZ, _sgWorldMinX, _sgWorldMinY, _sgWorldMinZ, _sgCellSizeX, _sgCellSizeY, _sgCellSizeZ, _softDataGrids); //EAS read only 1 file
		else if (fileExtension == "grd3") readSucessfull = MPS::io::readTIFromGS3DGRD3File(_softDataFileNames[i], _softDataGrids[i]);
		if (!readSucessfull) {
			_softDataGrids.clear();
			if (_debugMode>-1) {
				std::cout << "Error reading softdata " << _softDataFileNames[i] << std::endl;
			}
		}
	}

}


/**
* @brief Read soft data
*/
void MPS::MPSAlgorithm::_readTIFromFiles(void) {
	bool readSucessfull = false;
	std::string fileExtension = MPS::utility::getExtension(_tiFilename);

	if (fileExtension == "csv" || fileExtension == "txt") readSucessfull = MPS::io::readTIFromGS3DCSVFile(_tiFilename, _TI);
	else if (fileExtension == "dat" || fileExtension == "gslib" || fileExtension == "sgems" || fileExtension == "SGEMS") readSucessfull = MPS::io::readTIFromGSLIBFile(_tiFilename, _TI);
	else if (fileExtension == "grd3") readSucessfull = MPS::io::readTIFromGS3DGRD3File(_tiFilename, _TI);
	if (!readSucessfull) {
		std::cerr << "Error reading TI " << _tiFilename << std::endl;
		exit(-1);
	}

	//Checking the TI array dimensions
	_tiDimX = (int)_TI[0][0].size();
	_tiDimY = (int)_TI[0].size();
	_tiDimZ = (int)_TI.size();

	// _getCategories();

}


/**
* @brief Mask data
*/
void MPS::MPSAlgorithm::_readMaskDataFromFile(void) {
	bool readSucessfull = false;
	std::string fileExtension = MPS::utility::getExtension(_maskDataFileName);

	if (fileExtension == "csv" || fileExtension == "txt") readSucessfull = MPS::io::readTIFromGS3DCSVFile(_maskDataFileName, _maskDataGrid);
	else if (fileExtension == "dat" || fileExtension == "gslib" || fileExtension == "sgems" || fileExtension == "SGEMS") readSucessfull = MPS::io::readTIFromGSLIBFile(_maskDataFileName, _maskDataGrid);
	else if (fileExtension == "grd3") readSucessfull = MPS::io::readTIFromGS3DGRD3File(_maskDataFileName, _maskDataGrid);
	if (!readSucessfull) {
		if (_debugMode>-1) {
			std::cout << "Mask Data is missing" << _maskDataFileName << std::endl;
		}
	}

	//Initialize the flag for mask data availibility
	_hasMaskData = readSucessfull;

	// Out of bound check
	// The mask grid has to have the same dimension as the simulation grid, if there is any different then just ignore the mask grid and fire an error message
	if (readSucessfull) {
		int mdgDimX = (int)_maskDataGrid[0][0].size();
		int mdgDimY = (int)_maskDataGrid[0].size();
		int mdgDimZ = (int)_maskDataGrid.size();
		if ((_sgDimX != mdgDimX) || (_sgDimY != mdgDimY) || (_sgDimZ != mdgDimZ)) {
			std::cout << "The mask grid has to have the same dimension as the simulation grid " << _maskDataFileName << std::endl;
			_hasMaskData = false;
		}
	}
}

/**
* @brief Mask data
*/


void MPS::MPSAlgorithm::_getCategories(void) {

	_DataCategories.clear();
	
	// Loop through TI to find unique categories _Categegories
	for (int z=0; z<_tiDimZ; z+=1) {
		for (int y=0; y<_tiDimY; y+=1) {
			for (int x=0; x<_tiDimX; x+=1) {
				//For each pixel
				float val = _TI[z][y][z];
				if (std::find(_DataCategories.begin(), _DataCategories.end(), val) != _DataCategories.end()) {
					//
				} else {
					std::cout << val << "  ";
					_DataCategories.push_back(val);
				}
			}
		}
	}
	
	// sort Data Categories	
	std::sort (_DataCategories.begin(), _DataCategories.end() );
	// 

	if (_debugMode>-1) {
		std::cout << "Found " << _DataCategories.size() << " unique categories" << std::endl;
	}
	
	// Update the soft data categories read in par file
	_softDataCategories = _DataCategories;


  
  if (_DataCategories.size() < 20 ) {

		if (_debugMode>-1) {
			std::cout << "_DataCategories are: ";
			for (unsigned int i = 0; i < _DataCategories.size(); i++) {
				std::cout << _DataCategories[i] << " ";
			}
			std::cout << std::endl;
		

			std::cout << "_softDataCategories are: ";
			for (unsigned int i = 0; i < _softDataCategories.size(); i++) {
				std::cout << _softDataCategories[i] << " ";
			}
			std::cout << std::endl;
		}
	} else {
		if (_debugMode>-1) {
			std::cout << "-- probably a continious training image!";		
		}
	}
}

/**
* @brief Fill a simulation grid node from hard data and a search radius
* @param x coordinate X of the current node
* @param y coordinate Y of the current node
* @param z coordinate Z of the current node
* @param level current grid level
* @param addedNodes list of added nodes
* @param putbackNodes list of node to put back later
*/
void MPS::MPSAlgorithm::_fillSGfromHD(const int& x, const int& y, const int& z, const int& level, std::vector<MPS::Coords3D>& addedNodes, std::vector<MPS::Coords3D>& putbackNodes) {
	//Searching closest value in hard data and put that into the simulation grid
	//Only search for the node not NaN and within the radius
	//Do this only if have a hard data defined or
	//If current node already has value
	if(!_hdg.empty() && MPS::utility::is_nan(_sg[z][y][x])) {
		MPS::Coords3D closestCoords;
		//if (_IsClosedToNodeInGrid(x, y, z, level, _hdg, _hdSearchRadius, closestCoords)) {
		if (_IsClosedToNodeInGrid(x, y, z, level, _hdg, ceil(pow(2, level) / 2), closestCoords)) { //Limit within the direct neighbor
			//Adding the closest point to a list to desallocate after
			//addedNodes.push_back(closestCoords);
			putbackNodes.push_back(closestCoords);
			//Adding the current location to a list to desallocate after
			addedNodes.push_back(MPS::Coords3D(x, y, z));
			//Temporally put the closest node found to the sg cell
			_sg[z][y][x] = _hdg[closestCoords.getZ()][closestCoords.getY()][closestCoords.getX()];
			//Temporally put NaN value to the hdg value relocated so the same point will not be relocated more than 2 times
			_hdg[closestCoords.getZ()][closestCoords.getY()][closestCoords.getX()] = std::numeric_limits<float>::quiet_NaN();
			//_sg[z][y][x] = std::numeric_limits<float>::quiet_NaN();
			//std::cout << x << " " << y << " " << z << " " << pow(2, level) << " " << level << std::endl;
		}
	}
}

/**
* @brief Clear the SG nodes from the list of added nodes found by _fillSGfromHD
* @param addedNodes list of added nodes
* @param putbackNodes list of node to put back later
*/
void MPS::MPSAlgorithm::_clearSGFromHD(std::vector<MPS::Coords3D>& addedNodes, std::vector<MPS::Coords3D>& putbackNodes) {
	//Cleaning the allocated data from the SG
	for (int i=0; i<addedNodes.size(); i++) {
		_hdg[putbackNodes[i].getZ()][putbackNodes[i].getY()][putbackNodes[i].getX()] = _sg[addedNodes[i].getZ()][addedNodes[i].getY()][addedNodes[i].getX()];
		if ((addedNodes[i].getZ() != putbackNodes[i].getZ()) || (addedNodes[i].getY() != putbackNodes[i].getY()) || (addedNodes[i].getX() != putbackNodes[i].getX())) {
			_sg[addedNodes[i].getZ()][addedNodes[i].getY()][addedNodes[i].getX()] = std::numeric_limits<float>::quiet_NaN();
		}
	}
	//for (std::list<MPS::Coords3D>::iterator ptToBeRelocated = addedNodes.begin(); ptToBeRelocated != addedNodes.end(); ++ptToBeRelocated) {
	//	_sg[ptToBeRelocated->getZ()][ptToBeRelocated->getY()][ptToBeRelocated->getX()] = std::numeric_limits<float>::quiet_NaN();
	//}
	addedNodes.clear();
	putbackNodes.clear();
}



/**
* @brief Read a line of configuration file and put the result inside a vector data
* @param file filestream
* @param ss string stream
* @param data output data
* @param s string represents each data
* @param str string represents the line
* @return true if the line contains data
*/
bool MPS::MPSAlgorithm::_readLineConfiguration(std::ifstream& file, std::stringstream& ss, std::vector<std::string>& data, std::string& s, std::string& str) {
	ss.clear();
	data.clear();
	getline(file, str);
	ss << str;
	while (getline(ss, s, '#')) {
		s.erase(std::remove_if(s.begin(), s.end(), [](char x){return std::isspace(x);}), s.end()); //Removing spaces
		if(!s.empty()) {
			data.push_back(s);
		}
	}

	return (data.size() > 1);
}


/**
* @brief Read a line of configuration file and put the result inside a vector data
* @param file filestream
* @param ss string stream
* @param data output data
* @param s string represents each data
* @param str string represents the line
* @return true if the line contains data
*/
bool MPS::MPSAlgorithm::_readLineConfiguration_mul(std::ifstream& file, std::stringstream& ss, std::vector<std::string>& data, std::string& s, std::string& str) {

  std::string s2;
	ss.clear();
	data.clear();
	getline(file, str);
	ss << str;
	int count=0;
	while (getline(ss, s, '#')) {
    count++;

    if (count==1) {
			// BEFORE #
			//std::cout << "  count=" << count << std::endl;
			s.erase(std::remove_if(s.begin(), s.end(), [](char x){return std::isspace(x);}), s.end()); //Removing spaces
			if(!s.empty()) {
					data.push_back(s);
			}

		} else {
			// AFTER #: SPLIT AT #
			std::stringstream sss;
			sss << s;
			//std::cout << "s=" << s << std::endl << std::endl;
			int count2=0;
			while (getline(sss, s2, ' ')) {
				count2++;
				if (count2>1) {
					//std::cout << "count2=" << count2 << " s2=" << s2 << std::endl;
					data.push_back(s2);
				}
			}
		}

	}

	//for (unsigned int j=0; j<data.size(); j++) {
	//	std::cout << "data[" << j << "]=" <<  data[j] << std::endl;
	//}
	return (data.size() > 1);
}



/**
* @brief Shuffle the simulation grid path based preferential to soft data
* @param level current multi grid level
*/               
bool MPS::MPSAlgorithm::_shuffleSgPathPreferentialToSoftData(const int& level, std::list<MPS::Coords3D>& allocatedNodesFromSoftData) {
	// PATH PREFERENTIAL BY ENTROPY OF SOFT DATA
	// facEntropy=0 --> prefer soft data, but disregard entroopy
	// facEntropy>0 --> higher number means more deterministic path based increasingly on Entropy
	float facEntropy = _shuffleEntropyFactor;
	int offset = int(std::pow(2, level));
	float randomValue;				// ePath: forst col is random number, second col is an integer (index of node in SG grid)
	std::multimap<float, int> ePath;
	int node1DIdx;
	bool isRelocated = false;
	bool isAlreadyAllocated = false;
	//std::list<MPS::Coords3D> allocatedNodesFromSoftData; //Using to allocate the multiple grid with closest hd values
	std::map<float, float> softPdf;
	MPS::Coords3D closestCoords;

	int c=0;
	//Looping through each index of each multiple grid
	for (int z=0; z<_sgDimZ; z+= offset) {
		for (int y=0; y<_sgDimY; y+= offset) {
			for (int x=0; x<_sgDimX; x+= offset) {
				isRelocated = false;
				randomValue = ((float) rand() / (RAND_MAX));
				// Any Soft data??
				if (_getCpdfFromSoftData(x, y, z, level, softPdf, closestCoords)) {

					// We found a conditional data in the vicinity
					if (_debugMode > 2) {
						c = c + 1;
						std::cout << "nsoft=" << c << std::endl;
						std::cout << "  [x,y,z]   = " << x << "," << y << "," << z << " ## ";
						std::cout << softPdf[0] << "," << softPdf[1] << std::endl;
						std::cout << "  [x,y,z]C  = " << closestCoords.getX() << "," << closestCoords.getY() << "," << closestCoords.getZ() << "," << std::endl;
					}
					
					//}
					//Check if the closest node found already in the current relocated node
					isAlreadyAllocated = (std::find(allocatedNodesFromSoftData.begin(), allocatedNodesFromSoftData.end(), closestCoords) != allocatedNodesFromSoftData.end());
					float E; // total Entropy
					float Ei; // Partial Entropy
					float I; // Information content
					float p; // probability
					float q; // a priori 1D marginal
					if (!isAlreadyAllocated) {
						isRelocated = x != closestCoords.getX() || y != closestCoords.getY() || z != closestCoords.getZ();
						E=0;
						for(auto it = softPdf.cbegin(); it != softPdf.cend(); ++it) {
							
							//Put the relocated softdata into the softdata grid to continue the simulation
							// Remember to remove AFTER simulation is done
							if (isRelocated) {
								if (_debugMode > 2) {
									std::cout << "  RELOCATING CAT:" << it->first << " p:" << it->second << std::endl;
								}
								for (unsigned int i = 0; i<_softDataCategories.size(); i++) {
									if (it->first == _softDataCategories[i]) {
										_softDataGrids[i][z][y][x] = it->second;										
										break;
									}
								}
							}
							
							//// COMPUTE ENTROPY FOR SOFT PDF
							if (_shuffleSgPath==2) {
								p = it->second;
								if (p==0) {
									Ei=0;
								} else {
									Ei=p*(-1*log2(p));
									//Ei=p*(-1*log(p) / log(2.)); //msvc2012 doesnt have log2 yet
								}
								E=E+Ei;
							}
						}

						I=1-E; // Information content

						// Order Soft preference based on Entropy, and use a factor to control
						// the importance of Entropy
						randomValue = randomValue - 1 - facEntropy*I;

						if (_debugMode>1) {
							std::cout <<  "SOFT DATA -- ";
							std::cout << "cnt=" << node1DIdx << " x=" << x << " y=" << y << " z=" << z;
							std::cout << " -- E=" << E;
							std::cout << ", I=" << I;
							std::cout << ", randomV=" << randomValue << std::endl;
						}
					}
				}
				// assign random value to current node!
				MPS::utility::treeDto1D(x,y,z, _sgDimX, _sgDimY, node1DIdx);
				// assign random value to node in fine gird closest to current node!
				//MPS::utility::treeDto1D(closestCoords.getX(), closestCoords.getY(), closestCoords.getZ(), _sgDimX, _sgDimY, node1DIdx);
				ePath.insert ( std::pair<float,int>(randomValue, node1DIdx) );
				//ePath[randomValue] = cnt;
				//If closestCoords are different than current coordinate that mean there is a relocation so save the node to reinitialize
				if (isRelocated) {
					MPS::Coords3D nodeToBeReinitialized;
					nodeToBeReinitialized.setX(x);
					nodeToBeReinitialized.setY(y);
					nodeToBeReinitialized.setZ(z);
					allocatedNodesFromSoftData.push_back(nodeToBeReinitialized);
				}
			}
		}
	}


	if (_debugMode>-1) {
		std::cout << "Shuffling path, using preferential path with facEntropy=" << facEntropy << std::endl;
	}

	// Update simulation path
	int i=0;
	for(auto it = ePath.cbegin(); it != ePath.cend(); ++it)	{
		_simulationPath[i] = it->second;
		i++;
	}

	if (_debugMode>1) {
		std::cout << "PATH = " << std::endl;
		int tmpX, tmpY, tmpZ;

		for (auto i = _simulationPath.begin(); i != _simulationPath.end(); ++i) {
			MPS::utility::oneDTo3D(*i, _sgDimX, _sgDimY, tmpX, tmpY, tmpZ);
			std::cout << tmpX << "," << tmpY << "," << tmpZ << "  ";
		}
		std::cout << std::endl << "PATH END" << std::endl;
	}
	return 1;
}

/**
* @brief Start the simulation
* Virtual function implemented from MPSAlgorithm
*/
void MPS::MPSAlgorithm::startSimulation(void) {

	// Write license information to screen
	if (_debugMode>-1) {
		std::cout << "__________________________________________________________________________________" << std::endl;
	 	std::cout << "MPSlib: a C++ library for multiple point simulation" << std::endl;
	 	std::cout << "(c) 2015-2016 I-GIS (www.i-gis.dk) and" << std::endl;
		std::cout << "              Solid Earth Geophysics, Niels Bohr Institute (http://imgp.nbi.ku.dk)" << std::endl;
		std::cout << "This program comes with ABSOLUTELY NO WARRANTY;"  << std::endl;
    	std::cout << "This is free software, and you are welcome to redistribute it"  << std::endl;
    	std::cout << "under certain conditions. See 'COPYING.LESSER'for details."  << std::endl;
		std::cout << "__________________________________________________________________________________" << std::endl;
	}


	//Intitialize random seed or not
	if (_seed != 0) srand(_seed); //same seed
	else srand(time(NULL)); //random seed

	//Get the output filename
	size_t found;
	found = _tiFilename.find_last_of("/\\");
	std::string outputFilename = _outputDirectory + "/" + _tiFilename.substr(found+1);

	//Doing the simulation
	double totalSecs = 0;
	clock_t endNode, beginRealization, endRealization;
	double elapsedRealizationSecs, elapsedNodeSecs;
	int nodeEstimatedSeconds, seconds, hours, minutes, lastProgress;
	// HARD DATA RELOCATION
	std::vector<MPS::Coords3D> allocatedNodesFromHardData; //Using to allocate the multiple grid with closest hd values
	std::vector<MPS::Coords3D> nodeToPutBack; //Using to allocate the multiple grid with closest hd values
	// SOFT DATA RELOCATION
	std::list<MPS::Coords3D> allocatedNodesFromSoftData; //Using to allocate the multiple grid with closest hd values


	lastProgress = 0;
	int nodeCnt = 0, totalNodes = 0;
	int sg1DIdx, offset;

	for (int n=0; n<_realizationNumbers; n++) {

		if (_debugMode >= 2) {
			std::cout << "MPSLIB: simulation realization " << n + 1 << "/" << _realizationNumbers << std::endl;
		}

		beginRealization = clock();
		//Initialize the iteration count grid
		_initializeSG(_sgIterations, _sgDimX, _sgDimY, _sgDimZ, 0);
		//Initialize Simulation Grid from hard data or with NaN value
		_initializeSG(_sg, _sgDimX, _sgDimY, _sgDimZ);
		//Initialize temporary grids if debugMode is high
		if (_debugMode>1) {
			// Initialize some extra grids for extra information
			_initializeSG(_tg1, _sgDimX, _sgDimY, _sgDimZ);
			_initializeSG(_tg2, _sgDimX, _sgDimY, _sgDimZ);
			_initializeSG(_tg3, _sgDimX, _sgDimY, _sgDimZ);
			_initializeSG(_tg4, _sgDimX, _sgDimY, _sgDimZ);
			_initializeSG(_tg5, _sgDimX, _sgDimY, _sgDimZ);
		}

		/*if(!_hdg.empty()) {
		std::cout << "Initialize from hard data " << _hardDataFileNames << std::endl;
		_initializeSG(_sg, _sgDimX, _sgDimY, _sgDimZ, _hdg, std::numeric_limits<float>::quiet_NaN());
		} else {
		if (_debugMode>0) {
		std::cout << "Initialize with NaN value" << std::endl;
		}
		_initializeSG(_sg, _sgDimX, _sgDimY, _sgDimZ);
		}*/

		//Multi level grids
		for (int level=_totalGridsLevel; level>=0; level--) {

			if (_debugMode >= 2) {
				std::cout << "MPSLIB: starting multigrid  " << level <<  std::endl;
			}


			if (_debugMode > 2) {
				//Writting SG to file before simulation is done at current level
				MPS::io::writeToGSLIBFile(outputFilename + "_sg_A_before_" + std::to_string(n) + "_level_" + std::to_string(level) + ".gslib", _sg, _sgDimX, _sgDimY, _sgDimZ);
			}

			_InitStartSimulationEachMultipleGrid(level);

			//For each space level from coarse to fine
			offset = int(std::pow(2, level));

			//Define a simulation path for each level
			if (_debugMode > -1) {
				std::cout << "Define simulation path for level " << level << std::endl;
			}
			_simulationPath.clear();
			//std::cout << allocatedNodesFromHardData.size() << std::endl;

			nodeCnt = 0;
			totalNodes = (_sgDimX / offset) * (_sgDimY / offset) * (_sgDimZ / offset);
			for (int z=0; z<_sgDimZ; z+=offset) {
				for (int y=0; y<_sgDimY; y+=offset) {
					for (int x=0; x<_sgDimX; x+=offset) {
						//Fill the simulation node with the value from hard data grid
						if (!_hdg.empty() && MPS::utility::is_nan(_sg[z][y][x])) {
							_sg[z][y][x] = _hdg[z][y][x];
						}
						//Changing the simulation path based on a mask grid
						//Doing the simulation within the mask grid
						if (_hasMaskData) {
							if (_maskDataGrid[z][y][x] == 1) { //Doing the simulation only if mask data is 1
								MPS::utility::treeDto1D(x, y, z, _sgDimX, _sgDimY, sg1DIdx);
								_simulationPath.push_back(sg1DIdx);
							}
						}
						else { //No mask data, just do the simulation like before
							MPS::utility::treeDto1D(x, y, z, _sgDimX, _sgDimY, sg1DIdx);
							_simulationPath.push_back(sg1DIdx);
						}

						//The relocation process happens if the current simulation grid value is still NaN
						//Moving hard data to grid node only on coarsed level
						if (level != 0) _fillSGfromHD(x, y, z, level, allocatedNodesFromHardData, nodeToPutBack);

						//Progression
						if (_debugMode > -1 && !_hdg.empty()) {
							nodeCnt++;
							//Doing the progression
							//Print progression on screen
							int progress = (int)((nodeCnt / (float)totalNodes) * 100);
							if ((progress % 10) == 0 && progress != lastProgress) { //Report every 10%
								lastProgress = progress;
								std::cout << "Relocating hard data to the simulation grid at level: " << level << " Progression (%): " << progress << std::endl;
							}
						}
					}
				}
			}
			
			if (_debugMode > 2) {
				MPS::io::writeToGSLIBFile(outputFilename + "_sg_B_after_hard_relocation_before_simulation_" + std::to_string(n) + "_level_" + std::to_string(level) + ".gslib", _sg, _sgDimX, _sgDimY, _sgDimZ);
				std::cout << "After relocation" << std::endl;
				_showSG();
			}


			if (!_softDataGrids.empty() && _debugMode > 2) {
				MPS::io::writeToGSLIBFile(outputFilename + "_sdg_before_shuffling_" + std::to_string(n) + "_level_" + std::to_string(level) + ".gslib", _softDataGrids, _sgDimX, _sgDimY, _sgDimZ, _softDataCategories.size());
			}

			//std::cout << allocatedNodesFromHardData.size() << std::endl << std::endl;
			//Shuffle simulation path indices vector for a random path
			if (_debugMode > -1) {
				std::cout << "Shuffling simulation path using type " << _shuffleSgPath << std::endl;
			}

			//Back to random path if no soft data
			if (_softDataGrids.empty() && _shuffleSgPath == 2) {
				if (_debugMode > 2) {
					std::cout << "WARNING: no soft data found, switch to random path" << std::endl;
				}
				_shuffleSgPath = 1;				
			}
			//Shuffling
			if (_shuffleSgPath==1) {
				// random shuffling
				std::random_shuffle ( _simulationPath.begin(), _simulationPath.end() );
			} else if (_shuffleSgPath>1) {
				// shuffling preferential to soft data
				_shuffleSgPathPreferentialToSoftData(level, allocatedNodesFromSoftData);
			}

			if (!_softDataGrids.empty() && _debugMode > 2) {
				MPS::io::writeToGSLIBFile(outputFilename + "_sdg_after_shuffling_" + std::to_string(n) + "_level_" + std::to_string(level) + ".gslib", _softDataGrids, _sgDimX, _sgDimY, _sgDimZ, _softDataCategories.size());
			}

			//Performing the simulation
			//For each value of the path
			int progressionCnt = 0;
			int totalNodes = (int)_simulationPath.size();

			int SG_idxX, SG_idxY, SG_idxZ;
			if (_debugMode > 1) {
				std::cout << "________________________________________" << std::endl;
				std::cout << "_______ START SIMULATION _______________" << std::endl; 
			}

			////Cleaning the allocated data from the SG
			//_clearSGFromHD(allocatedNodesFromHardData);
			
			if (_debugMode > 1) {
				MPS::io::writeToGSLIBFile(outputFilename + "_path_" + std::to_string(n) + "_level_" + std::to_string(level) + ".gslib", _simulationPath, _simulationPath.size(), 1, 1);
			}
			for (unsigned int ii=0; ii<_simulationPath.size(); ii++) {
				//Get node coordinates
				MPS::utility::oneDTo3D(_simulationPath[ii], _sgDimX, _sgDimY, SG_idxX, SG_idxY, SG_idxZ);

				if (_debugMode > 2) {
					std::cout << " " << std::endl;
					std::cout << "at node = "<< ii <<"/"<< _simulationPath.size()<< std::endl;
				}

				//Performing simulation for non NaN value ...
				if (MPS::utility::is_nan(_sg[SG_idxZ][SG_idxY][SG_idxX]))
					_sg[SG_idxZ][SG_idxY][SG_idxX] = _simulate(SG_idxX, SG_idxY, SG_idxZ, level);

				if (_debugMode > -1) {
					//Doing the progression
					//Print progression on screen
					int progress = (int)((progressionCnt / (float)totalNodes) * 100);
					progressionCnt ++;
					if ((progress % 5) == 0 && progress != lastProgress) { //Report every 5%
						lastProgress = progress;
						endNode = clock();
						elapsedNodeSecs = double(endNode - beginRealization) / CLOCKS_PER_SEC;
						nodeEstimatedSeconds = (int)((elapsedNodeSecs/(float)(progressionCnt)) * (float)(totalNodes - progressionCnt));
						MPS::utility::secondsToHrMnSec(nodeEstimatedSeconds, hours, minutes, seconds);
						if (progress > 0) //Ignore first time that cant provide any time estimation
							std::cout << "Realization: " << n + 1 << "/" << _realizationNumbers << " Level: " << level << " Progression (%): " << progress << " finish in: " << hours << " h " << minutes << " mn " << seconds << " sec" << std::endl;
					}
				}
			}
			if (_debugMode > 2) {
				MPS::io::writeToGSLIBFile(outputFilename + "_sg_C_after_simulation_before_cleaning_relocation_" + std::to_string(n) + "_level_" + std::to_string(level) + ".gslib", _sg, _sgDimX, _sgDimY, _sgDimZ);
				std::cout << "After simulation" << std::endl;
				_showSG();
			}

			//Cleaning the allocated data from the SG
			if(level != 0) _clearSGFromHD(allocatedNodesFromHardData, nodeToPutBack);

			//Cleaning the allocated data from the Soft Data Grid
			if (level != 0 && allocatedNodesFromSoftData.size()>0) {
					for (std::list<MPS::Coords3D>::iterator ptToBeRelocated = allocatedNodesFromSoftData.begin(); ptToBeRelocated != allocatedNodesFromSoftData.end(); ++ptToBeRelocated) {
						// Clear relocated Soft data
						for (unsigned int i=0; i<_softDataCategories.size(); i++) {
							_softDataGrids[i][ptToBeRelocated->getZ()][ptToBeRelocated->getY()][ptToBeRelocated->getX()] = std::numeric_limits<float>::quiet_NaN();
						}
						// Clear simulated har data at relocated soft data location
						_sg[ptToBeRelocated->getZ()][ptToBeRelocated->getY()][ptToBeRelocated->getX()] = std::numeric_limits<float>::quiet_NaN();
				}
			}
			// clear list of allocated nodes on soft data
			allocatedNodesFromSoftData.clear();


			if (_debugMode > 2) {
				MPS::io::writeToGSLIBFile(outputFilename + "_sg_D_after_hard_cleaning_relocation_" + std::to_string(n) + "_level_" + std::to_string(level) + ".gslib", _sg, _sgDimX, _sgDimY, _sgDimZ);
				std::cout << "After cleaning relocation" << std::endl;
				_showSG();
			}

			
			
		}

		if (_debugMode > 0) {
			_showSG();
		}

		if (_debugMode > -1) {
			endRealization = clock();
			elapsedRealizationSecs = double(endRealization - beginRealization) / CLOCKS_PER_SEC;
			totalSecs += elapsedRealizationSecs;
			std::cout << "Elapsed time (sec): " << elapsedRealizationSecs << "\t" << " total: " << totalSecs << std::endl;
		}

		if (_debugMode > -2) {
			//Write result to file
			if (_debugMode > -1) {
				std::cout << "Write simulation grid to hard drive..." << std::endl;
			}
			MPS::io::writeToGSLIBFile(outputFilename + "_sg_" + std::to_string(n) + ".gslib", _sg, _sgDimX, _sgDimY, _sgDimZ);
			if (_saveGrd3 ) {
				MPS::io::writeToGRD3File(outputFilename + "_sg_gs3d_" + std::to_string(n) + ".grd3", _sg, _sgDimX, _sgDimY, _sgDimZ, _sgWorldMinX, _sgWorldMinY, _sgWorldMinZ, _sgCellSizeX, _sgCellSizeY, _sgCellSizeZ, 3);
			}
			//MPS::io::writeToGS3DCSVFile(outputFilename + "_sg_gs3d_" + std::to_string(n) + ".csv", _sg, _sgDimX, _sgDimY, _sgDimZ, _sgWorldMinX, _sgWorldMinY, _sgWorldMinZ, _sgCellSizeX, _sgCellSizeY, _sgCellSizeZ);
			//MPS::io::writeToASCIIFile(outputFilename + "_sg_ascii" + std::to_string(n) + ".txt", _sg, _sgDimX, _sgDimY, _sgDimZ, _sgWorldMinX, _sgWorldMinY, _sgWorldMinZ, _sgCellSizeX, _sgCellSizeY, _sgCellSizeZ);
			//MPS::io::writeToGS3DCSVFile(outputFilename + "_ti_gs3d_" + std::to_string(n) + ".csv", _TI, _tiDimX, _tiDimY, _tiDimZ, _sgWorldMinX, _sgWorldMinY, _sgWorldMinZ, _sgCellSizeX, _sgCellSizeY, _sgCellSizeZ);
		}


		if (_debugMode>1) {
			//Write temporary grids to  file
			MPS::io::writeToGSLIBFile(outputFilename + "_tg1_" + std::to_string(n) + ".gslib", _tg1, _sgDimX, _sgDimY, _sgDimZ);
			MPS::io::writeToGSLIBFile(outputFilename + "_tg2_" + std::to_string(n) + ".gslib", _tg2, _sgDimX, _sgDimY, _sgDimZ);
			MPS::io::writeToGSLIBFile(outputFilename + "_tg3_" + std::to_string(n) + ".gslib", _tg3, _sgDimX, _sgDimY, _sgDimZ);
			MPS::io::writeToGSLIBFile(outputFilename + "_tg4_" + std::to_string(n) + ".gslib", _tg4, _sgDimX, _sgDimY, _sgDimZ);
			MPS::io::writeToGSLIBFile(outputFilename + "_tg5_" + std::to_string(n) + ".gslib", _tg5, _sgDimX, _sgDimY, _sgDimZ);
		}


		if (_debugMode > 1) {
			//Write random path to file
			MPS::io::writeToGSLIBFile(outputFilename + "_path_" + std::to_string(n) + ".gslib", _simulationPath, _sgDimX, _sgDimY, _sgDimZ);
		}
	}

	if (_debugMode > -1) {
		MPS::utility::secondsToHrMnSec((int)(totalSecs/_realizationNumbers), hours, minutes, seconds);
		std::cout << "Total simulation time " << totalSecs << "s"<<std::endl;
		std::cout << "Average time for " << _realizationNumbers << " simulations (hours:minutes:seconds) : " << hours << ":" << minutes << ":" << seconds << std::endl;
	}

	if(_debugMode > -1 ) {
		std::cout << "Number of threads: " << _numberOfThreads << std::endl;
		// std::cout << "_maxNeighbours: " << _maxNeighbours << std::endl;
		// std::cout << " maxIterations: " << _maxIterations << std::endl;
		std::cout << "SG: " << _sgDimX << " " << _sgDimY << " " << _sgDimZ << std::endl;
		std::cout << "TI: " << _tiFilename << " " << _tiDimX << " " << _tiDimY << " " << _tiDimZ << " " << _TI[0][0][0]<< std::endl;
	}
}

/**
* @brief Filling L and V vectors with data
* @param grid the grid to search data
* @param idxX search index X
* @param idxY search index Y
* @param idxZ search index Z
* @param foundCnt how many node found
* @param maxNeighboursLimit maximum number of neigbor nodes count
* @param sgIdxX index X of the node in the 3D grid
* @param sgIdxY index Y of the node in the 3D grid
* @param sgIdxZ index Z of the node in the 3D grid
* @param L output vector distances between a found nodes and the currrent node
* @param V output vector values of the found nodes
* @return true if foundCnt is greater than max neighbors allowed
*/
bool MPS::MPSAlgorithm::_addingData(const std::vector<std::vector<std::vector<float>>>& grid, const int& idxX, const int& idxY, const int& idxZ, int& foundCnt, const int& maxNeighboursLimit, const int& sgIdxX, const int& sgIdxY, const int& sgIdxZ, std::vector<MPS::Coords3D>& L, std::vector<float>& V) {
	if (!MPS::utility::is_nan(grid[idxZ][idxY][idxX])) {
		foundCnt++;
		if (foundCnt > maxNeighboursLimit) return true;
		MPS::Coords3D aCoords;
		aCoords.setX(idxX - sgIdxX);
		aCoords.setY(idxY - sgIdxY);
		aCoords.setZ(idxZ - sgIdxZ);
		bool isPresent = (std::find(L.begin(), L.end(), aCoords) != L.end());
		if(!isPresent) {
			L.push_back(aCoords);
			V.push_back(grid[idxZ][idxY][idxX]);
		}
	}
	return false;
}

/**
* @brief Search data in a direction
* @param grid the grid to search data
* @param direction search direction (0: direction X, 1: direction Y, 2: direction Z)
* @param idxX search index X
* @param idxY search index Y
* @param idxZ search index Z
* @param foundCnt counter of found nodes
* @param maxNeighboursLimit maximum number of neigbor nodes count
* @param xOffset offset in X dimension of searching node
* @param yOffset offset in Y dimension of searching node
* @param zOffset offset in Z dimension of searching node
* @param sgIdxX index X of the node in the 3D grid
* @param sgIdxY index Y of the node in the 3D grid
* @param sgIdxZ index Z of the node in the 3D grid
* @param L output vector distances between a found nodes and the currrent node
* @param V output vector values of the found nodes
* @return true if foundCnt is greater than max neighbors allowed
*/
void MPS::MPSAlgorithm::_searchDataInDirection(const std::vector<std::vector<std::vector<float>>>& grid, const int& direction, int& idxX, int& idxY, int& idxZ, int& foundCnt, const int& maxNeighboursLimit, const int& xOffset, const int& yOffset, const int& zOffset, const int& sgIdxX, const int& sgIdxY, const int& sgIdxZ, std::vector<MPS::Coords3D>& L, std::vector<float>& V) {
	if(direction == 0) { //Direction X
		for(int k=-yOffset; k<=yOffset; k++) {
			idxY = sgIdxY + k;
			for(int j=-zOffset; j<=zOffset; j++) {
				idxZ = sgIdxZ + j;
				//Adding value inside viewport only
				if((idxX >= 0 && idxX < _sgDimX) && (idxY >= 0 && idxY < _sgDimY) && (idxZ >= 0 && idxZ < _sgDimZ)) {
					if (_addingData(grid, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, sgIdxX, sgIdxY, sgIdxZ, L, V)) break;
				}
			}
		}
	} else if(direction == 1) { //Direction Y
		for(int k=-xOffset+1; k<xOffset; k++) {
			idxX = sgIdxX + k;
			for(int j=-zOffset+1; j<zOffset; j++) {
				idxZ = sgIdxZ + j;
				//Adding value inside viewport only
				if((idxX >= 0 && idxX < _sgDimX) && (idxY >= 0 && idxY < _sgDimY) && (idxZ >= 0 && idxZ < _sgDimZ)) {
					if (_addingData(grid, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, sgIdxX, sgIdxY, sgIdxZ, L, V)) break;
				}
			}
		}
	} else if(direction == 2) { //Direction Z
		for(int k=-xOffset+1; k<xOffset; k++) {
			idxX = sgIdxX + k;
			for(int j=-yOffset+1; j<yOffset; j++) {
				idxY = sgIdxY + j;
				//Adding value inside viewport only
				if((idxX >= 0 && idxX < _sgDimX) && (idxY >= 0 && idxY < _sgDimY) && (idxZ >= 0 && idxZ < _sgDimZ)) {
					if (_addingData(grid, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, sgIdxX, sgIdxY, sgIdxZ, L, V)) break;
				}
			}
		}
	}
}


/**
* @brief Searching a neighbor node using circular search and return back a vector L (distance between a found node and current node) and a vector V (value of the found node)
* @param sgIdxX index X in the simulation grid
* @param sgIdxY index Y in the simulation grid
* @param sgIdxZ index Z in the simulation grid
* @param grid the grid to search data
* @param maxNeighboursLimit maximum number of neigbor nodes count
* @param maxRadiusLimit maximum search radius allowed (-1 if not used any search radius limit)
* @param L output vector distances between a found nodes and the currrent node
* @param V output vector values of the found nodes
*/
void MPS::MPSAlgorithm::_circularSearch(const int& sgIdxX, const int& sgIdxY, const int& sgIdxZ, const std::vector<std::vector<std::vector<float>>>& grid, const int& maxNeighboursLimit, const float& maxRadiusLimit, std::vector<MPS::Coords3D>& L, std::vector<float>& V) {
	int foundCnt = 0;
	int idxX, idxY, idxZ;
	int xOffset, yOffset, zOffset;
	int maxXOffset = _sgDimX - 1;
	int maxYOffset = _sgDimY - 1;
	int maxZOffset = _sgDimZ - 1;

	int maxDim = std::max(std::max(maxXOffset, maxYOffset), maxZOffset);

	//Check center point
	if (!MPS::utility::is_nan(grid[sgIdxZ][sgIdxY][sgIdxX])) {
		foundCnt++;
		MPS::Coords3D aCoords(0, 0, 0);
		L.push_back(aCoords);
		V.push_back(grid[sgIdxZ][sgIdxY][sgIdxX]);		
	}

	//random direction
	//int randomDirection;

	for(int i=1; i<maxDim; i++) {
  	//maximum neighbor count check
		if (foundCnt > maxNeighboursLimit) break;
		//if (L.size() > maxNeighboursLimit) break;

		//maximum search radius check
		if (i > maxRadiusLimit && maxRadiusLimit != -1) break;

		//Initialize offset
		xOffset = yOffset = zOffset = i;

		// //Get a random search direction
		// randomDirection = rand() % 6;
		// if (_debugMode > 2) {
			// std::cout << "_circularSearch: i=" << i << ", maxDim=" << maxDim << ", Random search direction = " << randomDirection;
			// std::cout << ", foundCnt=" << foundCnt << ", maxNeighboursLimit=" << maxNeighboursLimit;
			// std::cout << ", L.size=" << L.size() << std::endl;
		// }
		// switch (randomDirection) {
		// case 0 : //X Y Z
				//direction +X
				idxX = sgIdxX + xOffset;
				_searchDataInDirection(grid, 0, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				//direction -X
				idxX = sgIdxX - xOffset;
				_searchDataInDirection(grid, 0, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				//direction +Y
				idxY = sgIdxY + yOffset;
				_searchDataInDirection(grid, 1, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				//direction -Y
				idxY = sgIdxY - yOffset;
				_searchDataInDirection(grid, 1, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				//direction +Z
				idxZ = sgIdxZ + zOffset;
				_searchDataInDirection(grid, 2, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				//direction -Z
				idxZ = sgIdxZ - zOffset;
				_searchDataInDirection(grid, 2, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);
		// case 1 : //X Z Y
				// //direction +X
				// idxX = sgIdxX + xOffset;
				// _searchDataInDirection(grid, 0, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				// //direction -X
				// idxX = sgIdxX - xOffset;
				// _searchDataInDirection(grid, 0, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				// //direction +Z
				// idxZ = sgIdxZ + zOffset;
				// _searchDataInDirection(grid, 2, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				// //direction -Z
				// idxZ = sgIdxZ - zOffset;
				// _searchDataInDirection(grid, 2, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				// //direction +Y
				// idxY = sgIdxY + yOffset;
				// _searchDataInDirection(grid, 1, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				// //direction -Y
				// idxY = sgIdxY - yOffset;
				// _searchDataInDirection(grid, 1, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);
		// case 2 : //Y X Z
				// //direction +Y
				// idxY = sgIdxY + yOffset;
				// _searchDataInDirection(grid, 1, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				// //direction -Y
				// idxY = sgIdxY - yOffset;
				// _searchDataInDirection(grid, 1, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				// //direction +X
				// idxX = sgIdxX + xOffset;
				// _searchDataInDirection(grid, 0, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				// //direction -X
				// idxX = sgIdxX - xOffset;
				// _searchDataInDirection(grid, 0, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				// //direction +Z
				// idxZ = sgIdxZ + zOffset;
				// _searchDataInDirection(grid, 2, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				// //direction -Z
				// idxZ = sgIdxZ - zOffset;
				// _searchDataInDirection(grid, 2, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);
		// case 3 : //Y Z X
				// //direction +Y
				// idxY = sgIdxY + yOffset;
				// _searchDataInDirection(grid, 1, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				// //direction -Y
				// idxY = sgIdxY - yOffset;
				// _searchDataInDirection(grid, 1, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				// //direction +Z
				// idxZ = sgIdxZ + zOffset;
				// _searchDataInDirection(grid, 2, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				// //direction -Z
				// idxZ = sgIdxZ - zOffset;
				// _searchDataInDirection(grid, 2, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				// //direction +X
				// idxX = sgIdxX + xOffset;
				// _searchDataInDirection(grid, 0, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				// //direction -X
				// idxX = sgIdxX - xOffset;
				// _searchDataInDirection(grid, 0, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);
		// case 4 : //Z X Y
				// //direction +Z
				// idxZ = sgIdxZ + zOffset;
				// _searchDataInDirection(grid, 2, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				// //direction -Z
				// idxZ = sgIdxZ - zOffset;
				// _searchDataInDirection(grid, 2, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				// //direction +X
				// idxX = sgIdxX + xOffset;
				// _searchDataInDirection(grid, 0, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				// //direction -X
				// idxX = sgIdxX - xOffset;
				// _searchDataInDirection(grid, 0, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				// //direction +Y
				// idxY = sgIdxY + yOffset;
				// _searchDataInDirection(grid, 1, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				// //direction -Y
				// idxY = sgIdxY - yOffset;
				// _searchDataInDirection(grid, 1, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);
		// default : //Z Y X
				// //direction +Z
				// idxZ = sgIdxZ + zOffset;
				// _searchDataInDirection(grid, 2, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				// //direction -Z
				// idxZ = sgIdxZ - zOffset;
				// _searchDataInDirection(grid, 2, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				// //direction +Y
				// idxY = sgIdxY + yOffset;
				// _searchDataInDirection(grid, 1, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				// //direction -Y
				// idxY = sgIdxY - yOffset;
				// _searchDataInDirection(grid, 1, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				// //direction +X
				// idxX = sgIdxX + xOffset;
				// _searchDataInDirection(grid, 0, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);

				// //direction -X
				// idxX = sgIdxX - xOffset;
				// _searchDataInDirection(grid, 0, idxX, idxY, idxZ, foundCnt, maxNeighboursLimit, xOffset, yOffset, zOffset, sgIdxX, sgIdxY, sgIdxZ, L, V);
		// }
	}
	//cout << "After searching: " << L.size() << " " << V.size() << endl;
}
