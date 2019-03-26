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
#pragma once

#include <vector>
#include <list>
#include <string>
#include <cmath>
#include <algorithm>
#include <atomic>
#include <thread>
#include <fstream>
#include <sstream>
/**
* @brief Main simulation package
*/
namespace MPS {
	class MPSAlgorithm;
	class Coords3D;
}

/**
* @brief Abstract Multiple Points Statistic algorithm
*
* This class contains some shared procedures and functions for different implementation of MPS algorithms.
* This class cannot be used directly, view DirectSimulation for an example of how to implement an algorithm from this base.
*/
class MPS::MPSAlgorithm
{
private:
	/**
	* @brief Initialize the Hard Data Grid with a value, default is NaN
	* @param hdg the Hard Data GRID
	* @param sgDimX dimension X of the grid
	* @param sgDimY dimension Y of the gri
	* @param sgDimZ dimension Z of the grid
	* @param value value of each grid node default is NAN
	*/
	void _initializeHDG(std::vector<std::vector<std::vector<float>>>& sg, const int& sgDimX, const int& sgDimY, const int& sgDimZ, const float& value = std::numeric_limits<float>::quiet_NaN());
	/**
	* @brief Check if the current node is closed to a node in a given grid
	* The Given grid could be softdata or harddata grid
	* @param x coordinate X of the current node
	* @param y coordinate Y of the current node
	* @param z coordinate Z of the current node
	* @param level current grid level
	* @param grid given grid
	* @param searchRadius search radius for closed point search (in world unit)
	* @param closestCoordinates closest coordinates found
	* @return True if found a closed node
	*/
	bool _IsClosedToNodeInGrid(const int& x, const int& y, const int& z, const int& level, const std::vector<std::vector<std::vector<float>>>& grid, const float& searchRadius, MPS::Coords3D& closestCoordinates);
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
	bool _addingData(const std::vector<std::vector<std::vector<float>>>& grid, const int& idxX, const int& idxY, const int& idxZ, int& foundCnt, const int& maxNeighboursLimit, const int& sgIdxX, const int& sgIdxY, const int& sgIdxZ, std::vector<MPS::Coords3D>& L, std::vector<float>& V);
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
	void _searchDataInDirection(const std::vector<std::vector<std::vector<float>>>& grid, const int& direction, int& idxX, int& idxY, int& idxZ, int& foundCnt, const int& maxNeighboursLimit, const int& xOffset, const int& yOffset, const int& zOffset, const int& sgIdxX, const int& sgIdxY, const int& sgIdxZ, std::vector<MPS::Coords3D>& L, std::vector<float>& V);

protected:
	/**
	* @brief The simulation grid
	*/
	std::vector<std::vector<std::vector<float>>> _sg;
	/**
	* @brief The hard data grid (same size as simulation grid)
	*/
	std::vector<std::vector<std::vector<float>>> _hdg;
	/**
	* @brief Temporary grid 1 - meaning define by type of sim-algorithm (same size as simulation grid)
	*/
	std::vector<std::vector<std::vector<float>>> _tg1;
	/**
	* @brief Temporary grid 2 - meaning define by type of sim-algorithm (same size as simulation grid)
	*/
	std::vector<std::vector<std::vector<float>>> _tg2;
	/**
	* @brief Temporary grid 3 - meaning define by type of sim-algorithm (same size as simulation grid)
	*/
	std::vector<std::vector<std::vector<float>>> _tg3;
	/**
	* @brief Temporary grid 4 - meaning define by type of sim-algorithm (same size as simulation grid)
	*/
	std::vector<std::vector<std::vector<float>>> _tg4;
	/**
	* @brief Temporary grid 4 - meaning define by type of sim-algorithm (same size as simulation grid)
	*/
	std::vector<std::vector<std::vector<float>>> _tg5;
	/**
	* @brief hard data search radius for multiple grids
	*/
	float _hdSearchRadius;
	/**
	* @brief A copy of the simulation grid used for debugging which counts the number of iterations
	*/
	std::vector<std::vector<std::vector<float>>> _sgIterations;
	/**
	* @brief The simulation path
	*/
	std::vector<int> _simulationPath;
	/**
	* @brief multigrids levels
	*/
	int _totalGridsLevel;
	/**
	* @brief Dimension X of the simulation Grid
	*/
	int _sgDimX;
	/**
	* @brief Dimension Y of the simulation Grid
	*/
	int _sgDimY;
	/**
	* @brief Dimension Z of the simulation Grid
	*/
	int _sgDimZ;
	/**
	* @brief Coordinate X Min of the simulation grid in world coordinate
	*/
	float _sgWorldMinX;
	/**
	* @brief Coordinate Y Min of the simulation grid in world coordinate
	*/
	float _sgWorldMinY;
	/**
	* @brief Coordinate Z Min of the simulation grid in world coordinate
	*/
	float _sgWorldMinZ;
	/**
	* @brief Size of a cell in X direction in world coordinate
	*/
	float _sgCellSizeX;
	/**
	* @brief Size of a cell in Y direction in world coordinate
	*/
	float _sgCellSizeY;
	/**
	* @brief Size of a cell in Z direction in world coordinate
	*/
	float _sgCellSizeZ;
	/**
	* @brief Maximum conditional data allowed
	*/
	int _maxCondData;
	/**
	* @brief Define type of random simulation grid path
	*/
	int _shuffleSgPath;
	/**
	* @brief Define entropy factor for random simulation path
	*/
	int _shuffleEntropyFactor;
	/**
	* @brief The number of realization created
	*/
	int _realizationNumbers;
	/**
	* @brief If in debug mode, some extra files and informations will be created
	* Different levels of debug are available :
	* -1 : no information
	* 0 : with on console text about time elapsed
	* 1 : with grid preview on console
	* 2 : extra files are exported (iteration counter) to output folder
	*/
	int _debugMode = 0;
	/**
	* @brief Show the simulation grid result in the console
	*/
	bool _showPreview = 1;
	/**
	* @brief Save output as GRD3 files for GeoScene3D
	*/
	bool _saveGrd3 = 0;
	/**
	* @brief Initial value of the simulation
	*/
	float _seed = 0;
	/**
	* @brief Maximum number of iterations
	*/
	unsigned int _maxIterations = std::numeric_limits<unsigned int>::max();
	
	/**
	* @brief Dimension X of the training image
	*/
	int _tiDimX;
	/**
	* @brief Dimension Y of the training image
	*/
	int _tiDimY;
	/**
	* @brief Dimension Z of the training image
	*/
	int _tiDimZ;
	/**
	* @brief Maximum neighbour allowed when doing the neighbour search function
	*/
	int _maxNeighbours = std::numeric_limits<int>::max();
	/**
	* @brief Make a random training image path
	*/
	bool _shuffleTiPath;
	/**
	* @brief Training image search path
	*/
	std::vector<int> _tiPath;
	/**
	* @brief Maximum threads used for the simulation
	*/
	int _numberOfThreads;
	/**
	* @brief Traininng image's filename
	*/
	std::string _tiFilename;
	/**
	* @brief Output directory to store the result
	*/
	std::string _outputDirectory;
	/**
	* @brief Hard data filename used for the simulation
	*/
	std::string _hardDataFileNames;
	/**
	* @brief Soft data filenames used for the simulation
	*/
	std::vector<std::string> _softDataFileNames;
	/**
	* @brief Mask data filename used for the simulation
	* the mask is used to limit the area used for the simulation
	* when the grid node is initialize to 1 that mean the node will be simulated
	* the mask grid has the same size as the simulation grid
	*/
	std::string _maskDataFileName;
	/**
	* @brief Soft data categories
	*/
	std::vector<float> _softDataCategories;
	//std::vector<float> _DataCategories;
	/**
	* @brief Softdata grid
	*/
	std::vector<std::vector<std::vector<std::vector<float>>>> _softDataGrids;
	/**
	* @brief The training image
	*/
	std::vector<std::vector<std::vector<float>>> _TI;
	/**
	* @brief The mask grid
	*/
	std::vector<std::vector<std::vector<float>>> _maskDataGrid;
	/**
	* @brief Flag to know if the maskData is available
	*/
	bool _hasMaskData;
	/**
	* @brief Threads used for the simulation
	*/
	std::vector<std::thread> _threads;
	/**
	* @brief Atomic flag used to sychronize threads
	*/
	std::atomic<bool> _jobDone;
	/**
	* @brief Read a line of configuration file and put the result inside a vector data
	* @param file filestream
	* @param ss string stream
	* @param data output data
	* @param s string represents each data
	* @param str string represents the line
	* @return true if the line contains data
	*/
	bool _readLineConfiguration(std::ifstream& file, std::stringstream& ss, std::vector<std::string>& data, std::string& s, std::string& str);
	/**
	* @brief Read a line of configuration file and put the result inside a vector data
	* @param file filestream
	* @param ss string stream
	* @param data output data
	* @param s string represents each data
	* @param str string represents the line
	* @return true if the line contains data
	*/
	bool _readLineConfiguration_mul(std::ifstream& file, std::stringstream& ss, std::vector<std::string>& data, std::string& s, std::string& str);
	/**
	* @brief Read different data (TI, hard and softdata from files)
	*/
	void _readDataFromFiles(void);
	
	/**
	* @param Read hard data
	*/
	void _readHardDataFromFiles(void);
	
	/**
	* @param Read soft data
	*/
	void _readSoftDataFromFiles(void);

	/**
	* @brief Read TI data
	*/
	void _readTIFromFiles(void);

	/**
	* @brief Read Mask data
	*/
	void _readMaskDataFromFile(void);

	/**
	* @brief Get unique set of categories from TI
	*/
	//void _getCategories(void);

	/**
	* @brief Initialize a sequential simulation path
	* @param sgDimX dimension X of the path
	* @param sgDimY dimension Y of the path
	* @param sgDimZ dimension Z of the path
	* @param path output simulation path
	*/
		void _initilizePath(const int& sgDimX, const int& sgDimY, const int& sgDimZ, std::vector<int>& path);
	/**
	* @brief Initialize the Simulation Grid with a value, default is NaN
	* @param sg the simulation GRID
	* @param sgDimX dimension X of the grid
	* @param sgDimY dimension Y of the gri
	* @param sgDimZ dimension Z of the grid
	* @param value value of each grid node default is NAN
	*/
	void _initializeSG(std::vector<std::vector<std::vector<float>>>& sg, const int& sgDimX, const int& sgDimY, const int& sgDimZ, const float& value = std::numeric_limits<float>::quiet_NaN());
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
	void _initializeSG(std::vector<std::vector<std::vector<float>>>& sg, const int& sgDimX, const int& sgDimY, const int& sgDimZ, const std::vector<std::vector<std::vector<float>>>& grid, const float& nanValue = -1);
	/**
	* @brief Show the SG in the terminal
	*/
	void _showSG(void) const;
	/**
	* @brief Fill a simulation grid node from hard data and a search radius
	* @param x coordinate X of the current node
	* @param y coordinate Y of the current node
	* @param z coordinate Z of the current node
	* @param level current grid level
	* @param addedNodes list of added nodes
	* @param putbackNodes list of node to put back later
	*/
	void _fillSGfromHD(const int& x, const int& y, const int& z, const int& level, std::vector<MPS::Coords3D>& addedNodes, std::vector<MPS::Coords3D>& putbackNodes);
	/**
	* @brief Clear the SG nodes from the list of added nodes found by _fillSGfromHD
	* @param addedNodes list of added nodes
	* @param putbackNodes list of node to put back later
	*/
	void _clearSGFromHD(std::vector<MPS::Coords3D>& addedNodes, std::vector<MPS::Coords3D>& putbackNodes);
	/**
	* @brief Shuffle the simulation grid path based preferential to soft data
	* @param level current multi grid level
	*/
	bool _shuffleSgPathPreferentialToSoftData(const int& level, std::list<MPS::Coords3D>& allocatedNodesFromSoftData);
	/**
	* @brief Generate a realization from a PDF defined as a map
	* @param the pdf as a std::map
	* @param a realization from the pdf
	*/
	float _sampleFromPdf(std::map<float, float>& Pdf);
	/**
	* @brief Compute cpdf from softdata
	* @param x coordinate X of the current node
	* @param y coordinate Y of the current node
	* @param z coordinate Z of the current node
	* @param level current grid level
	* @param map with conditional pdf
	* @param closestCoords closest point used with relocation, if not then the current x, y, z is used
	* @return true if found a value
	*/
	bool _getCpdfFromSoftData(const int& x, const int& y, const int& z, const int& level, std::map<float, float>& conditionalPoints, MPS::Coords3D& closestCoords);
	/**
	* @brief simulation algorithm main function
	* @param sgIdxX index X of a node inside the simulation grind
	* @param sgIdxY index Y of a node inside the simulation grind
	* @param sgIdxZ index Z of a node inside the simulation grind
	* @param level level of the current grid
	* @return found node's value
	*/
	virtual float _simulate(const int& sgIdxX, const int& sgIdxY, const int& sgIdxZ, const int& level) = 0;
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
	void _circularSearch(const int& sgIdxX, const int& sgIdxY, const int& sgIdxZ, const std::vector<std::vector<std::vector<float>>>& grid, const int& maxNeighboursLimit, const float& maxRadiusLimit, std::vector<MPS::Coords3D>& L, std::vector<float>& V);

	/**
	* @brief Abstract function allow acces to the beginning of each simulation of each multiple grid
	* @param level the current grid level
	*/
	virtual void _InitStartSimulationEachMultipleGrid(const int& level) = 0;

public:
	/**
	* @brief Constructors
	*/
	MPSAlgorithm(void);
	/**
	* @brief Abstract function for starting the simulation
	*/
	virtual void startSimulation(void) = 0;
	/**
	* @brief Destructors
	*/
	virtual ~MPSAlgorithm(void);
	/**
	* @brief get the simulation grid
	* @return the simulation grid
	*/
	inline std::vector<std::vector<std::vector<float> > > sg() const {return _sg;}
	/**
	* @brief set the simulation grid
	* @param sg new simulation grid
	*/
	inline void setSg(const std::vector<std::vector<std::vector<float> > > &sg) {_sg = sg;}
	/**
	* @brief get the iterations simulation grid
	* @return iterations simulation grid
	*/
	inline std::vector<std::vector<std::vector<float> > > sgIterations() const {return _sgIterations;}
	/**
	* @brief set the iterations simulation grid
	* @param sgIterations new iterations simulation grid
	*/
	inline void setsgIterations(const std::vector<std::vector<std::vector<float> > > &sgIterations) {_sgIterations = sgIterations;}
	/**
	* @brief get the simulation path
	* @return the simulation path
	*/
	inline std::vector<int> simulationPath() const {return _simulationPath;}
	/**
	* @brief set the simulation path
	* @param simulationPath new simulation path
	*/
	inline void setSimulationPath(const std::vector<int> &simulationPath) {_simulationPath = simulationPath;}
	/**
	* @brief get the simulation grid dimension X
	* @return simulation grid dimension X
	*/
	inline int sgDimX() const {return _sgDimX;}
	/**
	* @brief set new dimension X to the simulation grid. The sg need to be reinitialized after calling this function
	* @param sgDimX new dimension X of the sg
	*/
	inline void setSgDimX(int sgDimX) {_sgDimX = sgDimX;}
	/**
	* @brief get the simulation grid dimension Y
	* @return simulation grid dimension Y
	*/
	inline int sgDimY() const {return _sgDimY;}
	/**
	* @brief set new dimension Y to the simulation grid. The sg need to be reinitialized after calling this function
	* @param sgDimY new dimension Y of the sg
	*/
	inline void setSgDimY(int sgDimY) {_sgDimY = sgDimY;}
	/**
	* @brief get the simulation grid dimension Z
	* @return simulation grid dimension Z
	*/
	inline int sgDimZ() const {return _sgDimZ;}
	/**
	* @brief set new dimension Z to the simulation grid. The sg need to be reinitialized after calling this function
	* @param sgDimZ new dimension Z of the sg
	*/
	inline void setSgDimZ(int sgDimZ) {_sgDimZ = sgDimZ;}
	/**
	* @brief Getter shuffleSgPath
	* @return shuffleSgPath
	*/
	inline int shuffleSgPath() const {return _shuffleSgPath;}
	/**
	* @brief Setter ShuffleSgPath
	* @param shuffleSgPath new value
	*/
	inline void setShuffleSgPath(int shuffleSgPath) {_shuffleSgPath = shuffleSgPath;}
	/**
	* @brief get the realization numbers
	* @return realization number
	*/
	inline int realizationNumbers() const {return _realizationNumbers;}
	/**
	* @brief set the realization numbers
	* @param realizationNumbers new realization numbers
	*/
	inline void setRealizationNumbers(int realizationNumbers) {_realizationNumbers = realizationNumbers;}
	/**
	* @brief get is debug mode
	* @return debug mode
	*/
	inline int debugMode() const {return _debugMode;}
	/**
	* @brief set debug mode
	* @param debugMode new debug mode
	*/
	inline void setDebugMode(int debugMode) {_debugMode = debugMode;}
	/**
	* @brief get is preview showed
	* @return is preview showed
	*/
	inline bool showPreview() const {return _showPreview;}
	/**
	* @brief set ShowPreview
	* @param showPreview new showpreview
	*/
	inline void setShowPreview(bool showPreview) {_showPreview = showPreview;}
	/**
	* @brief Getter maxIterations
	* @return maxIterations
	*/
	inline int maxIterations() const {return _maxIterations;}
	/**
	* @brief Setter MaxIterations
	* @param maxIterations new value
	*/
	inline void setMaxIterations(int maxIterations) {_maxIterations = maxIterations;}
	/**
	* @brief Getter numberOfThreads
	* @return numberOfThreads
	*/
	inline int numberOfThreads() const {return _numberOfThreads;}
	/**
	* @brief Setter NumberOfThreads
	* @param numberOfThreads new value
	*/
	inline void setNumberOfThreads(int numberOfThreads) {_numberOfThreads = numberOfThreads;}
	/**
	* @brief Getter tiFilename
	* @return tiFilename
	*/
	inline std::string tiFilename() const {return _tiFilename;}
	/**
	* @brief Setter TiFilename
	* @param tiFilename new value
	*/
	inline void setTiFilename(const std::string &tiFilename) {_tiFilename = tiFilename;}
	/**
	* @brief Getter outputDirectory
	* @return outputDirectory
	*/
	inline std::string outputDirectory() const {return _outputDirectory;}
	/**
	* @brief Setter OutputDirectory
	* @param outputDirectory new value
	*/
	inline void setOutputDirectory(const std::string &outputDirectory) {_outputDirectory = outputDirectory;}
	/**
	* @brief Getter hardDataFileNames
	* @return hardDataFileNames
	*/
	inline std::string hardDataFileNames() const {return _hardDataFileNames;}
	/**
	* @brief Setter HardDataFileNames
	* @param hardDataFileNames new value
	*/
	inline void setHardDataFileNames(const std::string &hardDataFileNames) {_hardDataFileNames = hardDataFileNames;}
	/**
	* @brief Getter TI
	* @return TI
	*/
	inline std::vector<std::vector<std::vector<float> > > TI() const {return _TI;}
	/**
	* @brief Setter TI
	* @param TI new value
	*/
	inline void setTI(const std::vector<std::vector<std::vector<float> > > &TI) {_TI = TI;}
	/**
	* @brief Getter tiDimX
	* @return tiDimX
	*/
	inline int tiDimX() const {return _tiDimX;}
	/**
	* @brief Setter TiDimX
	* @param tiDimX new value
	*/
	inline void setTiDimX(int tiDimX) {_tiDimX = tiDimX;}
	/**
	* @brief Getter tiDimY
	* @return tiDimY
	*/
	inline int tiDimY() const {return _tiDimY;}
	/**
	* @brief Setter TiDimY
	* @param tiDimY new value
	*/
	inline void setTiDimY(int tiDimY) {_tiDimY = tiDimY;}
	/**
	* @brief Getter tiDimZ
	* @return tiDimZ
	*/
	inline int tiDimZ() const {return _tiDimZ;}
	/**
	* @brief Setter TiDimZ
	* @param tiDimZ new value
	*/
	inline void setTiDimZ(int tiDimZ) {_tiDimZ = tiDimZ;}
	/**
	* @brief get the maximum neighbours number used for neighbour search
	* @return maximum neighbour number
	*/
	inline int maxNeighbours() const {return _maxNeighbours;}
	/**
	* @brief set maximum neighbours number
	* @param maxNeighbours new maximum neighbour number
	*/
	inline void setMaxNeighbours(int maxNeighbours) {_maxNeighbours = maxNeighbours;}
	/**
	* @brief Getter shuffleTiPath
	* @return shuffleTiPath
	*/
	inline bool shuffleTiPath() const {return _shuffleTiPath;}
	/**
	* @brief Setter ShuffleTiPath
	* @param shuffleTiPath new value
	*/
	inline void setShuffleTiPath(bool shuffleTiPath) {_shuffleTiPath = shuffleTiPath;}
	/**
	* @brief Getter tiPath
	* @return tiPath
	*/
	inline std::vector<int> tiPath() const {return _tiPath;}
	/**
	* @brief Setter TiPath
	* @param tiPath new value
	*/
	inline void setTiPath(const std::vector<int> &tiPath) {_tiPath = tiPath;}

};
