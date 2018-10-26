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
#include <cctype>       // isspace
#include <algorithm>    // std::previousMatchedCntuffle std::remove_if
#include <map>

#include "ENESIM.h"
#include "Utility.h"
#include "IO.h"
#include "Coords3D.h"

/**
* @brief Constructors
*/
MPS::ENESIM::ENESIM(void) : MPS::MPSAlgorithm(){
	// Multiple grids are not used in ENESIM type algorithms
	_totalGridsLevel=0;
}

/**
* @brief Destructors
*/
MPS::ENESIM::~ENESIM(void) {

}

/**
* @brief Read configuration file
* @param fileName configuration filename
*/
void MPS::ENESIM::_readConfigurations(const std::string& fileName) {

	// Check that parameter file exist
	if  (!(MPS::io::file_exist(fileName))) {
		std::cout << "Paremeter file '" << fileName << "' does not exist -> quitting" << std::endl;
		exit(0);
		return;
	}


	std::ifstream file(fileName);
	std::string str;
	std::stringstream ss;
	std::string s;
	std::vector<std::string> data;
	int read_int; // dummy integer

	// Process each line
	// Number of realizations
	_readLineConfiguration(file, ss, data, s, str);
	_realizationNumbers = stoi(data[1]);
	// Seed value
	_readLineConfiguration(file, ss, data, s, str);
	_seed = stof(data[1]);
	// Maximum number of counts for seeting up the conditional pdf
	_readLineConfiguration(file, ss, data, s, str);
	read_int  =  stoi(data[1]);
	if (read_int<0) {
		_nMaxCountCpdf = std::numeric_limits<int>::max();
	} else {
		_nMaxCountCpdf=read_int;
	}

	// Maximum number of conditional data (_maxNeighbours)
	_readLineConfiguration_mul(file, ss, data, s, str);
	_maxNeighbours = stoi(data[1]);
	if (data.size()>2) {
		_maxNeighbours_soft = stoi(data[2]);
	}
	else {
		_maxNeighbours_soft = 0;
	}
	//std::cout << "_maxNeighbours_soft= " << _maxNeighbours_soft << std::endl;

	// Maximum iterations
	_readLineConfiguration(file, ss, data, s, str);
	read_int  =  stoi(data[1]);
	if (read_int<0) {
		_maxIterations = std::numeric_limits<int>::max();
	} else {
		_maxIterations = read_int;
	}
	

	// Distance Measure and minimum distance
	_readLineConfiguration_mul(file, ss, data, s, str);
	_distance_measure = stoi(data[1]);
	// optional. Template size x - base when n_mulgrid=0
	if (data.size()>2) {
		_distance_threshold = stof(data[2]);
	} else {
		_distance_threshold = 0;
	}
	if (data.size()>3) {
		_distance_power_order = stof(data[3]);
	} else {
		_distance_power_order = 0;
	}
		
	// ColocateDimension
	_readLineConfiguration(file, ss, data, s, str);
	_ColocatedDimension = stof(data[1]);

	// Maximum search Radius
	_readLineConfiguration_mul(file, ss, data, s, str);
	_maxSearchRadius = stof(data[1]);
	if (data.size()>2) {
		_maxSearchRadius_soft = stof(data[2]);
	}
	else {
		_maxSearchRadius_soft = 0;
	}
	//std::cout << "_maxSearchRadius_soft= " << _maxSearchRadius_soft << std::endl;

	// Simulation Grid size X
	_readLineConfiguration(file, ss, data, s, str);
	_sgDimX = stoi(data[1]);
	// Simulation Grid size Y
	_readLineConfiguration(file, ss, data, s, str);
	_sgDimY = stoi(data[1]);
	// Simulation Grid size Z
	_readLineConfiguration(file, ss, data, s, str);
	_sgDimZ = stoi(data[1]);
	// Simulation Grid World min X
	_readLineConfiguration(file, ss, data, s, str);
	_sgWorldMinX = stof(data[1]);
	// Simulation Grid World min Y
	_readLineConfiguration(file, ss, data, s, str);
	_sgWorldMinY = stof(data[1]);
	// Simulation Grid World min Z
	_readLineConfiguration(file, ss, data, s, str);
	_sgWorldMinZ = stof(data[1]);
	// Simulation Grid Cell Size X
	_readLineConfiguration(file, ss, data, s, str);
	_sgCellSizeX = stof(data[1]);
	// Simulation Grid Cell Size Y
	_readLineConfiguration(file, ss, data, s, str);
	_sgCellSizeY = stof(data[1]);
	// Simulation Grid Cell Size Z
	_readLineConfiguration(file, ss, data, s, str);
	_sgCellSizeZ = stof(data[1]);
	// TI filename
	_readLineConfiguration(file, ss, data, s, str);
	data[1].erase(std::remove_if(data[1].begin(), data[1].end(), [] (char x){return std::isspace(x); }), data[1].end()); //Removing spaces
	_tiFilename = data[1];
	// Output directory
	_readLineConfiguration(file, ss, data, s, str);
	data[1].erase(std::remove_if(data[1].begin(), data[1].end(), [] (char x){return std::isspace(x); }), data[1].end()); //Removing spaces
	_outputDirectory = data[1];
	// Shuffle SGPATH
	_readLineConfiguration_mul(file, ss, data, s, str);
	_shuffleSgPath = (stoi(data[1]));
	if (data.size()>2) {
		// EntropyFactor is the second value, if it exists
		_shuffleEntropyFactor = stof(data[2]);
	} else {
		_shuffleEntropyFactor = 4;
	}
	// Shuffle TI Path
	_readLineConfiguration(file, ss, data, s, str);
	_shuffleTiPath = (stoi(data[1]) != 0);
	// Hard data
	_readLineConfiguration(file, ss, data, s, str);
	data[1].erase(std::remove_if(data[1].begin(), data[1].end(), [] (char x){return std::isspace(x); }), data[1].end()); //Removing spaces
	_hardDataFileNames = data[1];
	// Hard data search radius
	_readLineConfiguration(file, ss, data, s, str);
	_hdSearchRadius = stof(data[1]);
	// Softdata categories
	if(_readLineConfiguration(file, ss, data, s, str)) {
		// std::cout << data.size() << " " << data.empty() << std::endl;
		ss.clear();
		ss << data[1];
		while (getline(ss, s, ';')) {
			//Parsing categories
			if(!s.empty()) {
				_softDataCategories.push_back(stof(s));
			}
		}
	}
	// Softdata filenames
	if(_readLineConfiguration(file, ss, data, s, str)) {
		ss.clear();
		ss << data[1];
		while (getline(ss, s, ';')) {
			//Parsing categories
			if(!s.empty()) {
				s.erase(std::remove_if(s.begin(), s.end(), [] (char x){return std::isspace(x); }), s.end()); //Removing spaces
				_softDataFileNames.push_back(s);
			}
		}
		// Resize softdata grids
		_softDataGrids.resize(_softDataFileNames.size());
	}
	// Number of threads
	_readLineConfiguration(file, ss, data, s, str);
	_numberOfThreads = stoi(data[1]);
	// DEBUG MODE
	_readLineConfiguration(file, ss, data, s, str);
	_debugMode = stoi(data[1]);
}

/**
* @brief Compute cpdf from training image using ENESIM type
* @param x coordinate X of the current node
* @param y coordinate Y of the current node
* @param z coordinate Z of the current node
* @param map with conditional pdf
* @return true if found a value
*/
bool MPS::ENESIM::_getCpdfTiEnesim(const int& sgIdxX, const int& sgIdxY, const int& sgIdxZ, std::map<float, float>& cPdf) {
	int maxCpdfCount = _nMaxCountCpdf;

	// map containing the count of conditional data values
	std::map<float, int> conditionalCount;

	int CpdfCount = 0; // the number of matches found, hence the number of counts used to construct the condtitional
	float valueFromTI; 
	float LC_dist_current; // the 'current' distance'
	float LC_dist_min; // The minimal distance found
	float TI_x_min, TI_y_min, TI_z_min; // The position in the TI associated with the minimal distance
	int TI_idxX, TI_idxY, TI_idxZ;
	float V_center_ti = -1; // value of the central node in the TI
	int node1DIdx; // Integer to hold 1d index from 3D point
	std::vector<MPS::Coords3D> L_s; // relative location in soft grid
	std::vector<float> V_s; // value at realtive location in soft grid


	// Find (_maxNeighbours) the closest data in simulation grid (_sg)
	// and store the relive locatoin in L_h, and the values in V_h
	std::vector<MPS::Coords3D> L_c;
	std::vector<float> V_c;
	_circularSearch(sgIdxX, sgIdxY, sgIdxZ, _sg, _maxNeighbours, _maxSearchRadius, L_c, V_c);
	
	if (_debugMode>2) {
		std::cout << "_getCpdfTiEnesim: -- ENESIM TOP --"<< std::endl;
		std::cout << "SGxyx=(" << sgIdxX << "," << sgIdxY <<  "," << sgIdxZ <<")" << std::endl;
		std::cout << "nLxyz=" << L_c.size() << "(of max "<< _maxNeighbours << ")" << std::endl;
		std::cout << "circularSearch: _maxNeighbours=" << _maxNeighbours << " _maxSearchRadius=" << _maxSearchRadius << std::endl;

	}

	// Compute relative distance for each conditional data
	std::vector<float> L_dist(L_c.size());
	for (unsigned int i=0; i<L_c.size(); i++) {
		 if (_debugMode>3) {
			std::cout << "  Lxyz=("<<L_c[i].getX() << ","<<L_c[i].getY() <<  ","<<L_c[i].getZ() << ")"<< " V="<<V_c[i]<< std::endl;
		}
		if (_distance_power_order!=0) {
			L_dist[i] = sqrt(L_c[i].getX()*L_c[i].getX() + L_c[i].getY()*L_c[i].getY() + L_c[i].getZ()*L_c[i].getZ());
			L_dist[i]  = pow(L_dist[i],-1*_distance_power_order) ;

		} else {
			L_dist[i]=1;
		}
	}

	
	// The path scanning the training image is shifted such that a random start location is chosen
	int ti_shift;
	ti_shift = (std::rand() % (int)(_tiDimX*_tiDimY*_tiDimZ));
	std::rotate(_tiPath.begin(), _tiPath.begin() + ti_shift, _tiPath.end());

	// At this point V_c represents the conditional value(s), and L_c realtive position(s) oF of the contitional (hard) values(s).
	// to find the local conditional pdf f_TI(m_i | m_c)

	
	// For SOFT DATA
	std::map<float, std::vector<float>> conditionalSoftProb;

	
	// Loop over training image to find a close match to L_c and find ,V_c
	
	CpdfCount = 0; // Set number of count in Cpdf to zero
	LC_dist_min = std::numeric_limits<float>::max(); // set the current 'minimum' distance to a very high number
	float LC_dist; // distance of L,V to value in TI
	unsigned int i_ti_path;
	
	for (i_ti_path=0; i_ti_path<_tiPath.size(); i_ti_path++) {
		MPS::utility::oneDTo3D(_tiPath[i_ti_path], _tiDimX, _tiDimY, TI_idxX, TI_idxY, TI_idxZ);

		// Get the centered value in the TI
		V_center_ti = _TI[TI_idxZ][TI_idxY][TI_idxX];

		// Get the distance between the conditional data in TI and SIM grid
		LC_dist = _computeDistanceLV_TI(L_c, V_c, TI_idxX, TI_idxY, TI_idxZ, L_dist);
		


		// Check if current L,T in TI match conditional observations better
		if (LC_dist<LC_dist_min) {
			
			// We have a new MIN distance, update the found conditional values 
			//   as well as the minimum distance and the associated location in the TI
			valueFromTI = V_center_ti; 
			LC_dist_current = LC_dist; // needed?
			LC_dist_min = LC_dist;
			TI_x_min = TI_idxX;
			TI_y_min = TI_idxY;
			TI_z_min = TI_idxZ;
			  
			
			// Add a count to the Cpdf if the current node match L,V according to some threshold..
			if (LC_dist<=_distance_threshold) {
				CpdfCount++;
				// Update conditionalCount Counter (from which the local cPdf can be computed)
				if (conditionalCount.find(valueFromTI) == conditionalCount.end()) {
					// Then we've encountered the word for a first time.
					// Is this slow?
					conditionalCount[valueFromTI] = 1; // Initialize it to	 1.
				} else {
					// Then we've already seen it before..
					conditionalCount[valueFromTI]++; // Just increment it.
				}
				if (_debugMode > 2) {
					std::cout << "Matching event  i_ti_path=" << i_ti_path << ", V_center_ti=" << valueFromTI << std::endl;
				}
				// MAKE SURE TO RESET LC_dist_min
				LC_dist_min = 1e+9;
			

			}

		}


		if (_debugMode > 3) {
			std::cout << "  " << "i_ti_path=" << i_ti_path<< " ===";
			std::cout << "  TI	xyz=(" << TI_idxX << "," << TI_idxY << "," << TI_idxZ << ")";
			std::cout << " - LC_dist=" << LC_dist << " min(LC_dist)=" << LC_dist_current;
			std::cout << std::endl;
		}




		if (maxCpdfCount <= CpdfCount ) {
			if (_debugMode > 2) {
				std::cout << "MaxCpdfCount Reached i_ti_path=" << i_ti_path << " - CpdfCount=" << CpdfCount << std::endl;
			}
			break;
		}


		// Stop looking if we have reached the maximum number of allowed iterations in TI
		if (i_ti_path>_maxIterations) {
			if (_debugMode > 2) {
				std::cout << "Max Ite Reached i_ti_path=" << i_ti_path << " - nmax_ite=" << _maxIterations << std::endl;
			}
			break;
		}

		//std::cout << "-- DONE --" << std::endl;

	} // END SCAN OF TI FOR CPDF
  //std::cout << sgIdxX << "," << sgIdxY << "  :h_dist_cum = " << h_dist_cum << "  LC_dist = "<< LC_dist << std::endl;

	
	// At this point set the local soft data as NaN to avoid simulate it again
	_softDataGrids[0][sgIdxZ][sgIdxY][sgIdxX] = std::numeric_limits<float>::quiet_NaN();

	// assign minimum distance to temporary grid 1
	if (_debugMode>1) {
		_tg1[sgIdxZ][sgIdxY][sgIdxX] = LC_dist_current;
		_tg2[sgIdxZ][sgIdxY][sgIdxX] = CpdfCount;
		MPS::utility::treeDto1D(TI_x_min, TI_y_min, TI_z_min, _tiDimX, _tiDimY, node1DIdx);
		_tg3[sgIdxZ][sgIdxY][sgIdxX] = node1DIdx;
		//_tg4[sgIdxZ][sgIdxY][sgIdxX] = _tiPath[i_ti_path]; // POSITION ON TI
		_tg4[sgIdxZ][sgIdxY][sgIdxX] = i_ti_path;
		_tg5[sgIdxZ][sgIdxY][sgIdxX] = L_c.size(); // number of used conditional points
	}



	// CHECK THAT conditionalCount HAS AT LEAST ONE COUNT
	// This may not always be the case when no conditional event has been found!
	if (CpdfCount==0 ) {
		// NO MATCHES HAS BEEN FOUND. ADD THE CURRENT BEST MACTH TO THE conditionalCount
		// std::cout << " CpdfCount " << CpdfCount << std::endl;
		if (conditionalCount.find(valueFromTI) == conditionalCount.end()) {
			// Then we've encountered the word for a first time. Is this slow?
			conditionalCount[valueFromTI] = 1; // Initialize it to	 1.
		} else {
			// Then we've already seen it before..
			conditionalCount[valueFromTI]++; // Just increment it.
		}
	}

	// COMPUTE ConditionalPdf from conditionalCount
	//Loop through all the conditional points to calculate the pdf (probability distribution function)
	//Getting the sum of counter and sort the map using the counter
	float totalCounter = 0;
	for(std::map<float,int>::iterator iter = conditionalCount.begin(); iter != conditionalCount.end(); ++iter) {
		totalCounter += iter->second;
	}
	// std::cout << "totalCounter=" << totalCounter << std::endl;

	//std::map<float, float> conditionalPdfFromTi;
	//Looping again in conditionalCount, compute and add the probabilities from TI
	for(std::map<float, int>::iterator iter = conditionalCount.begin(); iter != conditionalCount.end(); ++iter) {
		cPdf.insert(std::pair<float, float>(iter->first, (float)(iter->second) / (float)totalCounter));
	}
	// Now the conditional PDF from the TI is available as conditionalPdfFromTi;


	if (_debugMode>2) {
		std::cout << "NsimInTI=" << i_ti_path;
		std::cout << " --CpdfCount = " << CpdfCount;
		std::cout << " valueFromTI = " << valueFromTI;
		std::cout << " LC_dist_current = " << LC_dist_current;

		std::cout << "SGxyx=(" << sgIdxX << "," << sgIdxY <<  "," << sgIdxZ <<")" << std::endl;
		std::cout << "_getCpdfTiEnesim: -- ENESIM END --"<< std::endl;

	}



	return true;
}



/**
* @brief Compute cpdf from training image using ENESIM type
* @param x coordinate X of the current node
* @param y coordinate Y of the current node
* @param z coordinate Z of the current node
* @param map with conditional pdf
* @return true if found a value
*/
bool MPS::ENESIM::_getCpdfTiEnesimNew(const int& sgIdxX, const int& sgIdxY, const int& sgIdxZ, std::map<float, float>& cPdf, float& SoftProbability) {

	// MAKE SURE TO COUNT THE NUMBER OF SOFT CONDTIONING DATA ONCE, SUCH THAT WE DO NOT NEED
	// TO LOOK FOR SOFT DATA ONCE THEY HAVE BEEN SIMULATED!
	// COULD SPEED UP SNESIM AS WELL

	int maxCpdfCount = _nMaxCountCpdf;

	// map containing the count of conditional hard data values
	std::map<float, int> conditionalCount;
	// map containing the probability of conditional conditional data values
	std::map<float, float> conditionalSoftProb;
	
	SoftProbability = 1;

	int CpdfCount = 0; // the number of matches found, hence the number of counts used to construct the condtitional
	float valueFromTI;
	float LC_dist_current; // the 'current' distance'
	float LC_dist_min; // The minimal distance found
	float TI_x_min, TI_y_min, TI_z_min; // The position in the TI associated with the minimal distance
	int TI_idxX, TI_idxY, TI_idxZ;
	float V_center_ti = -1; // value of the central node in the TI
	int node1DIdx; // Integer to hold 1d index from 3D point
	std::vector<MPS::Coords3D> L_s; // relative location in soft grid
	std::vector<float> V_s; // value at realtive location in soft grid

	std::vector<float> P_soft_i;
	std::vector<float> P_soft_max_i;

	// initialize counters of hard and soft data stats 
	for (unsigned int i = 0; i < _softDataCategories.size(); i++) {
		conditionalCount[_softDataCategories[i]] = 0;
		conditionalSoftProb[_softDataCategories[i]] = 0;
	}
						
	// Find (_maxNeighbours) the closest data in simulation grid (_sg)
	// and store the relive locatoin in L_h, and the values in V_h
	std::vector<MPS::Coords3D> L_c;
	std::vector<float> V_c;
	_circularSearch(sgIdxX, sgIdxY, sgIdxZ, _sg, _maxNeighbours, _maxSearchRadius, L_c, V_c);

	if (_debugMode>2) {
		std::cout << "_getCpdfTiEnesimNew: -- ENESIM TOP --"; 
		std::cout << " SGxyz=(" << sgIdxX << "," << sgIdxY << "," << sgIdxZ << ")";
		std::cout << " nLxyz=" << L_c.size() << "(of max " << _maxNeighbours << ")";
		std::cout << " CircularSearch: _maxNeighbours=" << _maxNeighbours << " _maxSearchRadius=" << _maxSearchRadius << std::endl;
	}

//	std::cout << "START" << std::endl;
//	std::cout << "f(m_i|m_c)=[" << conditionalCount[_softDataCategories[0]] << "," << conditionalCount[_softDataCategories[1]] << "]";
//	std::cout << "  f_soft(m_i)=[" << conditionalSoftProb[_softDataCategories[0]] << "," << conditionalSoftProb[_softDataCategories[1]] << "]" << std::endl;


	// Compute relative distance for each conditional data
	std::vector<float> L_dist(L_c.size());
	for (unsigned int i = 0; i<L_c.size(); i++) {
		if (_debugMode>3) {
			std::cout << "  Lxyz=(" << L_c[i].getX() << "," << L_c[i].getY() << "," << L_c[i].getZ() << ")" << " V=" << V_c[i] << std::endl;
		}
		if (_distance_power_order != 0) {
			L_dist[i] = sqrt(L_c[i].getX()*L_c[i].getX() + L_c[i].getY()*L_c[i].getY() + L_c[i].getZ()*L_c[i].getZ());
			L_dist[i] = pow(L_dist[i], -1 * _distance_power_order);

		}
		else {
			L_dist[i] = 1;
		}
	}


	// Setup the path through the training image

	// Tirst check if the last dimension should be treated as individual types of random varibles
	// Then the only the lower dimensional X-Y TI-grid shoudl be search for matches
	// _StationaryLastDim

	if (_ColocatedDimension > 0) {
	
		// Check if the current _tiPath is set of for full size, and if so, reduce the size
		int smallTiSize = _tiDimX * _tiDimY; // size of the TI-grid to scan through
		if (_tiPath.size() != smallTiSize) {
			if (_debugMode > 2) {
				std::cout << "Using Co-located dimension = " << _ColocatedDimension << std::endl;
				std::cout << _tiPath.size() << "<" << smallTiSize << std::endl;
			}
			_tiPath.resize(smallTiSize);
			if (_debugMode > 2) {
				std::cout << "New reduced Ti path size = " << _tiPath.size() << std::endl;
			}
		}
	
				int node1DIdx;
		int cnt = 0;
		int z = sgIdxZ; // The current position in the simulation grid
		for (int y = 0; y<_tiDimY; y++) {
			for (int x = 0; x<_tiDimX; x++) {			
				MPS::utility::treeDto1D(x, y, z, _tiDimX, _tiDimY, node1DIdx);
				//_tiPath[cnt] = cnt + 1; // sequential path
				_tiPath[cnt] = node1DIdx;
				cnt++;
			}
		}
		/*
		std::cout << "z=" << z << std::endl;
		std::cout << "  _tiPath[0]=" << _tiPath[0] << ", _tiPath[end]=" << _tiPath[_tiPath.size() - 1] <<  ", N=" << _tiPath.size()<< std::endl;
		std::cout << "  _tiPath[0]=" << _tiPath[0] << ", _tiPath[end]=" << _tiPath.back() << std::endl;
		*/

	}

	// The path scanning the training image is shifted such that a random start location is chosen
	int ti_shift;
	ti_shift = (std::rand() % (int)(_tiPath.size() ));
	std::rotate(_tiPath.begin(), _tiPath.begin() + ti_shift, _tiPath.end());
	

	// At this point V_c represents the conditional value(s), and L_c realtive position(s) oF of the contitional (hard) values(s).
	// to find the local conditional pdf f_TI(m_i | m_c)

	// Loop over training image to find a close match to L_c and find ,V_c

	CpdfCount = 0; // Set number of count in Cpdf to zero
	LC_dist_min = std::numeric_limits<float>::max(); // set the current 'minimum' distance to a very high number
	float LC_dist; // distance of L,V to value in TI
	unsigned int i_ti_path;

	for (i_ti_path = 0; i_ti_path<_tiPath.size(); i_ti_path++) {
		MPS::utility::oneDTo3D(_tiPath[i_ti_path], _tiDimX, _tiDimY, TI_idxX, TI_idxY, TI_idxZ);

		// Get the centered value in the TI
		V_center_ti = _TI[TI_idxZ][TI_idxY][TI_idxX];

		// Get the distance between the conditional data in TI and SIM grid
		LC_dist = _computeDistanceLV_TI(L_c, V_c, TI_idxX, TI_idxY, TI_idxZ, L_dist);



		// Check if current L,T in TI match conditional observations better
		if (LC_dist<LC_dist_min) {

			// We have a new MIN distance, update the found conditional values 
			//   as well as the minimum distance and the associated location in the TI
			valueFromTI = V_center_ti;
			LC_dist_current = LC_dist; // needed?
			LC_dist_min = LC_dist;
			TI_x_min = TI_idxX;
			TI_y_min = TI_idxY;
			TI_z_min = TI_idxZ;

			// Add a count to the Cpdf if the current node match L,V according to some threshold..
			if (LC_dist <= _distance_threshold) {
				// We have a match of the hard data,
				// Increment total counter, and counter for the specific found center node node (V_center_ti)
				CpdfCount++;
				conditionalCount[valueFromTI]++; // Just increment it.


				if (_debugMode > 2) {
					std::cout << "Matching event  i_ti_path=" << i_ti_path << ", V_center_ti=" << valueFromTI << std::endl;
				}
				// MAKE SURE TO RESET LC_dist_min
				LC_dist_min = 1e+9;


				// Check if soft data are available! 
				float P_soft = 1;
				float P_soft_max = 1;

				
				if (!_softDataGrids.empty()) {

					int sdX, sdY, sdZ;
					int NLs;

					// FIND THE SOFT PROBABILITY ASSOCIATED to the  THE CURRENT MATCH!
					L_s.clear();
					V_s.clear();
					_circularSearch(sgIdxX, sgIdxY, sgIdxZ, _softDataGrids[0], _maxNeighbours_soft, _maxSearchRadius_soft, L_s, V_s);
					if (_debugMode>2) {
						std::cout << "--> Found N_soft_cond data = " << L_s.size() << std::endl;
					}
					P_soft_i.clear();
					P_soft_max_i.clear();

					if (L_s.size() > 0) {
						//std::cout << CpdfCount;
						//std::cout << "  - Found L_s.size()=" << L_s.size() << " conditional data using _maxNeighbours=" << _maxNeighbours_soft << "  _maxSearchRadius=" << _maxSearchRadius_soft;
						//std::cout << " at _sgd=(" << sgIdxX << "," << sgIdxY << "," << sgIdxZ << ")" << std::endl;
						for (unsigned int i = 0; i < L_s.size(); i++) {

							// find location of soft data in simulation/soft-data grid
							sdX = sgIdxX + L_s[i].getX();
							sdY = sgIdxY + L_s[i].getY();
							sdZ = sgIdxZ + L_s[i].getZ();
							float dist = sqrt(L_s[i].getZ()*L_s[i].getZ() + L_s[i].getY()*L_s[i].getY() + L_s[i].getX()*L_s[i].getX());
							// find location of soft data in the training image, _TI
							int TI_x_soft = TI_x_min + +L_s[i].getX();
							int TI_y_soft = TI_y_min + +L_s[i].getY();
							int TI_z_soft = TI_z_min + +L_s[i].getZ();


							// If soft data location in _TI is located within the define grid then use it
							if (TI_x_soft >= 0 && TI_x_soft < _tiDimX && TI_y_soft >= 0 && TI_y_soft  < _tiDimY && TI_z_soft >= 0 && TI_z_soft  < _tiDimZ) {
								// get value at soft position in TI
								//float sdValTI = _sg[sdZ][sdY][sdX];
								if (_debugMode > 4) {
									std::cout << " Soft pos in TI=(" << TI_x_soft << "," << TI_y_soft << "," << TI_z_soft << ")";
									std::cout << " ValInTI=" << _TI[TI_z_soft][TI_y_soft][TI_x_soft];
									std::cout << " relpos=(" << L_s[i].getX() << "," << L_s[i].getY() << "," << L_s[i].getZ() << ")";
									std::cout << " dist=" << dist;
									std::cout << " center=(" << sgIdxX << "," << sgIdxY << "," << sgIdxZ << ")";

									std::cout << " P=[" << _softDataGrids[0][sdZ][sdY][sdX] << "," << _softDataGrids[1][sdZ][sdY][sdX] << "]";
									std::cout << " V_C=" << V_s[i] << std::endl;
								}
								// find the soft probability of one conditional soft data, as P_soft_i
								// find the MAX soft probability of the one conditional soft data, as P_soft_max_i
								float p_max = 0;
								for (unsigned int icat = 0; icat < _softDataGrids.size(); icat++) {
									if (_softDataGrids[icat][sdZ][sdY][sdX] > p_max) {
										p_max = _softDataGrids[icat][sdZ][sdY][sdX];
									}
								}

								// find category id of found value in the TI at locations of the soft data
								unsigned int icat;	
								for (unsigned int i = 0; i < _softDataCategories.size(); i++) {
									//std::cout << "_softDataCategories[i]=" << _softDataCategories[i] << " == _TI[TI_z_soft][TI_y_soft][TI_x_soft]=" << _TI[TI_z_soft][TI_y_soft][TI_x_soft] << std::endl;
									if (_softDataCategories[i] == _TI[TI_z_soft][TI_y_soft][TI_x_soft]) {
										icat = i;
									}									
								}
								//std::cout << "--- icat=" << icat << std::endl;
								P_soft_i.push_back(_softDataGrids[icat][sdZ][sdY][sdX]);
								P_soft_max_i.push_back(p_max);


								//std::cout << "P=" << P_soft_i.back() << " TI(soft)=" << _TI[TI_z_soft][TI_y_soft][TI_x_soft] << std::endl;


							}
							else {
								// relative soft data located OUTSIDE TI grid. 
								// In this case assign the lowest soft probability is the local soft probabilty
								// assign the lowest possible soft probability to assign th eleast weight to this data
								
								// find the MIN and MAX soft probability of the one conditional soft data, as P_soft_max_i
								float p_min = 1;
								float p_max = 0;
								for (unsigned int icat = 0; icat < _softDataGrids.size(); icat++) {
									if (_softDataGrids[icat][sdZ][sdY][sdX] < p_min) {
										p_min = _softDataGrids[icat][sdZ][sdY][sdX];
									}
									if (_softDataGrids[icat][sdZ][sdY][sdX] > p_max) {
										p_max = _softDataGrids[icat][sdZ][sdY][sdX];
									}
								}
								P_soft_i.push_back(p_min);
								P_soft_max_i.push_back(p_max);

							}
						}
					
					}
					
					// We now have the probability for the maxNeighbours_soft closest soft data, with probabilities stored in P_soft_i
					if (_debugMode > 2) {
						std::cout << "--> Using N_soft_cond data = " << P_soft_i.size() << std::endl;
					}
					// get the combined soft probability
					P_soft = 0;
					P_soft_max = 0;
					if (P_soft_i.size() > 0) {
						P_soft = P_soft_i[0];
						P_soft_max = P_soft_max_i[0];
						if (_debugMode>2) {
							std::cout << "  [P_soft,P_soft_max]=" << P_soft_i[0] << "," << P_soft_max_i[0] << std::endl;
						}
						if (P_soft_i.size() > 1) {
							for (int index = 1; index < P_soft_i.size(); ++index) {
								if (_debugMode > 2) {
									std::cout << "  [P_soft,P_soft_max]=" << P_soft_i[index] << "," << P_soft_max_i[index] << std::endl;
								}
								P_soft = P_soft * P_soft_i[index];
								P_soft_max = P_soft_max * P_soft_max_i[index];
							}
						}
						// std::cout << std::endl;

						// Add weighed count to conditionalSoftProb
						SoftProbability = P_soft / P_soft_max;
						conditionalSoftProb[valueFromTI] = conditionalSoftProb[valueFromTI] + SoftProbability; 
						
					}
					
					
				} else {
					P_soft = 1;
				}
				
			}

		}

		if (_debugMode > 3) {
			std::cout << "  " << "i_ti_path=" << i_ti_path << " ===";
			std::cout << "  TI	xyz=(" << TI_idxX << "," << TI_idxY << "," << TI_idxZ << ")";
			std::cout << " - LC_dist=" << LC_dist << " min(LC_dist)=" << LC_dist_current;
			std::cout << std::endl;
		}

		if (maxCpdfCount <= CpdfCount) {
			if (_debugMode > 2) {
				std::cout << "MaxCpdfCount Reached i_ti_path=" << i_ti_path << " - CpdfCount=" << CpdfCount << std::endl;
			}
			break;
		}

		// Stop looking if we have reached the maximum number of allowed iterations in TI
		if (i_ti_path>_maxIterations) {
			if (_debugMode > 2) {
				std::cout << "Max Ite Reached i_ti_path=" << i_ti_path << " - nmax_ite=" << _maxIterations << std::endl;
			}
			break;
		}

		//std::cout << "-- DONE --" << std::endl;

	} // END SCAN OF TI FOR CPDF
	
	// We now have counts to establish the conditional from hard data f(m_i|m_c)
	// and, if avaialble, information from soft probabilities f(m_i|m_c)*f_soft(m)


	if (_debugMode > 2) {
		std::cout << "  f(m_i|m_c)=[" << conditionalCount[_softDataCategories[0]] << "," << conditionalCount[_softDataCategories[1]] << "]";
		std::cout << "  f(m_i|m_c)*f_soft(m)=[" << conditionalSoftProb[_softDataCategories[0]] << "," << conditionalSoftProb[_softDataCategories[1]] << "]" << std::endl;
	}

	

	
	// assign minimum distance to temporary grid 1
	if (_debugMode>1) {
		_tg1[sgIdxZ][sgIdxY][sgIdxX] = LC_dist_current;
		_tg2[sgIdxZ][sgIdxY][sgIdxX] = CpdfCount;
		MPS::utility::treeDto1D(TI_x_min, TI_y_min, TI_z_min, _tiDimX, _tiDimY, node1DIdx);
		_tg3[sgIdxZ][sgIdxY][sgIdxX] = node1DIdx;
		//_tg4[sgIdxZ][sgIdxY][sgIdxX] = _tiPath[i_ti_path]; // POSITION ON TI
		_tg4[sgIdxZ][sgIdxY][sgIdxX] = i_ti_path;
		_tg5[sgIdxZ][sgIdxY][sgIdxX] = L_c.size(); // number of used conditional points
	}



	// CHECK THAT conditionalCount HAS AT LEAST ONE COUNT
	// This may not always be the case when no conditional event has been found!
	if (CpdfCount == 0) {
		// NO MATCHES HAS BEEN FOUND. ADD THE CURRENT BEST MACTH TO THE conditionalCount
		// std::cout << " CpdfCount " << CpdfCount << std::endl;
		if (conditionalCount.find(valueFromTI) == conditionalCount.end()) {
			// Then we've encountered the word for a first time. Is this slow?
			conditionalCount[valueFromTI] = 1; // Initialize it to	 1.
		}
		else {
			// Then we've already seen it before..
			conditionalCount[valueFromTI]++; // Just increment it.
		}
	}

	// Get the the sum of the unnormalized conditional probabilities
	float totalCounter_conditionalCount = 0;
	for (std::map<float, int>::iterator iter = conditionalCount.begin(); iter != conditionalCount.end(); ++iter) {
		totalCounter_conditionalCount += iter->second;
	}
	float totalCounter_conditionalSoftProb = 0;
	for (std::map<float, float>::iterator iter = conditionalSoftProb.begin(); iter != conditionalSoftProb.end(); ++iter) {
		totalCounter_conditionalSoftProb += iter->second;
	}
	//std::cout << "totalCounter_conditionalSoftProb=" << totalCounter_conditionalSoftProb << std::endl;
	//std::cout << "totalCounter_conditionalCount=" << totalCounter_conditionalCount << std::endl;

	// normalize output
	if (totalCounter_conditionalSoftProb==0) {
		// NO SOFT INFORMATION
		// f(m_i | m_c)
		for (std::map<float, int>::iterator iter = conditionalCount.begin(); iter != conditionalCount.end(); ++iter) {
			cPdf.insert(std::pair<float, float>(iter->first, (float)(iter->second) / (float)totalCounter_conditionalCount));
		}

	} else {
		// USE CPDF CONDITIONAL TO SOFT
		// f(m_i|m_c)* \prod f_soft(m_i)
		for (std::map<float, float>::iterator iter = conditionalSoftProb.begin(); iter != conditionalSoftProb.end(); ++iter) {
			cPdf.insert(std::pair<float, float>(iter->first, (float)(iter->second) / (float)totalCounter_conditionalSoftProb));
		}
	}
	
	
	// Compute the Conditional Entropy 
	if (_computeConditionalEntropy) {
		float H = 0;
		float P = 0;
		float I = 0;
		int N = 2; // number of outcomes.. global parameter?
		for (std::map<float, float>::iterator iter = cPdf.begin(); iter != cPdf.end(); ++iter) {
			if (_debugMode > 2) {
				std::cout << " P(" << iter->first << ")=" << iter->second;
			}
			cPdf.insert(std::pair<float, float>(iter->first, (float)(iter->second) / (float)totalCounter_conditionalSoftProb));
			P = iter->second;
			if (P > 0) {
				I = -P * log(P) / log(N);
				H = H + I;
			}
		}
		if (_debugMode > 2) {
			std::cout << " H= " << H << ".";
			std::cout << std::endl;
		}
		_ceg[sgIdxZ][sgIdxY][sgIdxX] = H;
		
	}


	if (_debugMode > 2) {
		std::cout << "NsimInTI=" << i_ti_path;
		std::cout << " --CpdfCount = " << CpdfCount;
		std::cout << " valueFromTI = " << valueFromTI;
		std::cout << " LC_dist_current = " << LC_dist_current;

		std::cout << "SGxyx=(" << sgIdxX << "," << sgIdxY << "," << sgIdxZ << ")" << std::endl;
		std::cout << "_getCpdfTiEnesim: -- ENESIM END --" << std::endl;
	}


	return true;
}


/**
* @brief Compute distance between conditional data in TI and template L
* @param TIi_dxX coordinate X of the current node in TI
* @param TIi_dxY coordinate Y of the current node in TI
* @param TIi_dxZ coordinate Z of the current node in TI
* @param L_dist Precompute Distance!
* @return distance
*/
float MPS::ENESIM::_computeDistanceLV_TI(std::vector<MPS::Coords3D>& L, std::vector<float>& V, const int& TI_idxX, const int& TI_idxY, const int& TI_idxZ, std::vector<float>& L_dist) {

	float h_dist=0;
	float h_dist_cum=0;
	float LC_dist;
	float V_ti;
	int TI_x, TI_y, TI_z;

	LC_dist=0;
	//h_dist_cum=0;

	if (L.size()>0) {
		for (unsigned int i=0; i<L.size(); i++) {
			//For each pixel relatively to the current pixel based on vector L
			TI_x = TI_idxX + L[i].getX();
			TI_y = TI_idxY + L[i].getY();
			TI_z = TI_idxZ + L[i].getZ();

		  h_dist = L_dist[i];
			h_dist_cum = h_dist_cum + h_dist;

			// Check wgether the current conditional point
			// is located within the training image limits
			if((TI_x >= 0 && TI_x < _tiDimX) && (TI_y >= 0 && TI_y < _tiDimY) && (TI_z >= 0 && TI_z < _tiDimZ)) {
				V_ti = _TI[TI_z][TI_y][TI_x];

				if (_distance_measure==1) {
					// Discrete measure: no matching pixel means added distance of 1
					if (V_ti!=V[i]) {
						// add a distance of 1, if case of no matching pixels
						LC_dist=LC_dist+1*h_dist;
					}
				}	else if (_distance_measure==2){
					LC_dist = LC_dist + ((V_ti-V[i])*(V_ti-V[i]))/(256*256)*h_dist;
				}
			} else {
				// The conditioning location of the point to compare to is located outside
				// the traning image, and will be treated not a match
				if (_distance_measure==1) {
					LC_dist = LC_dist+1*h_dist;
				} else if (_distance_measure==2) {
					LC_dist = 1e+9*h_dist;
				}
			}
		}
		if (_distance_measure==1) {
			// Normaize LC_dist
			//LC_dist = LC_dist / h_dist_cum;
		} else if (_distance_measure==2) {
			LC_dist = LC_dist / L.size();
		}
	} else {
		// No conditional points --> distance zero
		LC_dist=0;
	}

	return LC_dist;
}

/**
* @brief Compute distance between conditional data in TI and template L.
* @brief Much slower than using Precomputed distance
* @param TIi_dxX coordinate X of the current node in TI
* @param TIi_dxY coordinate Y of the current node in TI
* @param TIi_dxZ coordinate Z of the current node in TI
* @return distance
*/
/*
float MPS::ENESIM::_computeDistanceLV_TI(std::vector<MPS::Coords3D>& L, std::vector<float>& V, const int& TI_idxX, const int& TI_idxY, const int& TI_idxZ) {

	float h_dist;
	float h_dist_cum;
	float LC_dist;
	float V_ti;
	int TI_x, TI_y, TI_z;

	LC_dist=0;
	h_dist_cum=0;

	if (L.size()>0) {
		for (unsigned int i=0; i<L.size(); i++) {
			//For each pixel relatively to the current pixel based on vector L
			TI_x = TI_idxX + L[i].getX();
			TI_y = TI_idxY + L[i].getY();
			TI_z = TI_idxZ + L[i].getZ();
			if (_distance_power_order!=0) {
				h_dist = sqrt(L[i].getX()*L[i].getX() + L[i].getY()*L[i].getY() + L[i].getZ()*L[i].getZ());
				h_dist = pow(h_dist,-1*_distance_power_order) ;
			} else {
				h_dist=1;
			}
			h_dist_cum = h_dist_cum + h_dist;

			// Check wgether the current conditional point
			// is located within the training image limits
			if((TI_x >= 0 && TI_x < _tiDimX) && (TI_y >= 0 && TI_y < _tiDimY) && (TI_z >= 0 && TI_z < _tiDimZ)) {
				V_ti = _TI[TI_z][TI_y][TI_x];

				if (_distance_measure==1) {
					// Discrete measure: no matching picel means added distance of 1
					if (V_ti!=V[i]) {
						// add a distance of 1, if case of no matching pixels
						LC_dist=LC_dist+1*h_dist;
					}
				}	else if (_distance_measure==2){
					LC_dist = LC_dist + ((V_ti-V[i])*(V_ti-V[i]))*h_dist;
				}
			} else {
				// The conditioning location of the point to compare to is located outside
				// the traning image, and will be treated not a match
				if (_distance_measure==1) {
					LC_dist = LC_dist+1*h_dist;
				} else if (_distance_measure==2) {
					LC_dist = 1000000*h_dist;
				}
			}
		}
		if (_distance_measure==1) {
			// Normaize LC_dist
			LC_dist = LC_dist / h_dist_cum;
		}
	} else {
		// No conditional points --> distance zero
		LC_dist=0;
	}

	return LC_dist;
}

*/



/**
* @brief Compbine two PDFs assuming independence
* @param probability density 1, pdf1
* @param probability density 2, pdf2
* @return true if succesfull, and update pdf is available in pdf1
*/
bool MPS::ENESIM::_combinePdf(std::map<float, float>& cPdf, std::map<float, float>& softPdf) {

	float key;
	float val_softpdf;
	float val_cpdf;
	float sum_cpdf;
	sum_cpdf=0;
	for(std::map<float,float>::iterator iter = softPdf.begin(); iter != softPdf.end(); ++iter) {
		key=iter->first;
		val_softpdf=softPdf[key];

		if (cPdf.count(key)) {
			// the key is set in both the soft and coniditional pdf.
			val_cpdf=cPdf[key];
			//std::cout << "  key=" << key << " P SO="<< val_softpdf << std::endl;

		} else {
			// the key is not set for the conditional pdf
			val_cpdf=0.00001;
		}
		//std::cout << "sum pdf = " << sum_cpdf << std::endl;
		//std::cout << "CPDF_TI   INDEX=" << iter->first << "  COUNT =" << iter->second << " ->";

		val_cpdf = val_cpdf*val_softpdf;
		sum_cpdf = sum_cpdf + val_cpdf;

		// Update the coniditonal pdf assuming independence of soft data and conditional pdf from TI
		cPdf[key]=val_cpdf;

	}

	// normalize combined pdf
	for(std::map<float,float>::iterator iter = cPdf.begin(); iter != cPdf.end(); ++iter) {
		cPdf[iter->first] = cPdf[iter->first] / sum_cpdf;
		// std::cout << "key=" << iter->first << " P TI="<< conditionalPdfFromTi[iter->first] << std::endl;
	}

	return true;
}

/**
* @brief Compute cpdf from softdata and set the value to the current node
* @param sgIdxX coordinate X of the current node
* @param sgIdxY coordinate Y of the current node
* @param sgIdxZ coordinate Z of the current node
* @param iterationCnt Iterations counter
* @return simulated value
*/
float MPS::ENESIM::_getRealizationFromCpdfTiEnesim(const int& sgIdxX, const int& sgIdxY, const int& sgIdxZ, float& iterationCnt) {

	// Get cPdf from training image using ENESIM style
	std::map<float, float> conditionalPdfFromTi;
	_getCpdfTiEnesim(sgIdxX, sgIdxY,sgIdxZ, conditionalPdfFromTi);

	// Check if any SoftData are available?
	//std::multimap<float, float> softPdf;
	std::map<float, float> softPdf;
	MPS::Coords3D closestCoords;
	if (_getCpdfFromSoftData(sgIdxX, sgIdxY, sgIdxZ, 0, softPdf, closestCoords)) {
		// Combine with conditional pdf from TI
		_combinePdf(conditionalPdfFromTi,softPdf);
	}

	// sample from conditional pdf
	float simulatedValue;
	simulatedValue = _sampleFromPdf(conditionalPdfFromTi);

	// DONE SIMULATING
	if (_debugMode > 2) {
		// std::cout << "Simulated value="<< simulatedValue << ", LC_dist_min=" << LC_dist_min << " " << std::endl ;
		std::cout << "Simulated value="<< simulatedValue << " " << std::endl ;
	}
	return simulatedValue;
}

/**
* @brief get realization from Cpdf and Soft data using a metropolis sampler when soft data is available
* @param sgIdxX coordinate X of the current node
* @param sgIdxY coordinate Y of the current node
* @param sgIdxZ coordinate Z of the current node
* @param iterationCnt Iterations counter
* @return simulated value
*/                 
float MPS::ENESIM::_getRealizationFromCpdfTiEnesimRejection(const int& sgIdxX, const int& sgIdxY, const int& sgIdxZ, float& iterationCnt) {

	std::map<float, float> conditionalPdfFromTi;
	// Check if any SoftData are available?
	//std::multimap<float, float> softPdf;
	float simulatedValue;
	float randomValue	;
	float pAcc;

	std::map<float, float> softPdf;
	int i=0;;
	int maxIterations = 100;  // decide whether soft or ti cpdf takes preference..
	bool isAccepted = false;
	MPS::Coords3D closestCoords;
	if (_getCpdfFromSoftData(sgIdxX, sgIdxY, sgIdxZ, 0, softPdf, closestCoords)) {

		do {
			// MAKE SURE THAT conditionalPdfFromTi IS OBTAINED INDEPENDENTLTY EACH TIME!!!
			// ESPECIELLY WHEN BASED ON ONLY ONE COUNT

			conditionalPdfFromTi.clear();
			_getCpdfTiEnesim(sgIdxX, sgIdxY,sgIdxZ, conditionalPdfFromTi);
			simulatedValue = _sampleFromPdf(conditionalPdfFromTi);

			randomValue = ((float) rand() / (RAND_MAX));
			pAcc = softPdf[simulatedValue];

			if(_debugMode > 2 ) {
				std::cout << "i= " << i << " maxIterations=" << maxIterations;
				std::cout << "   p_ti = [" <<  conditionalPdfFromTi[0] << "," << conditionalPdfFromTi[1] << "]";
				std::cout << " simval = " << simulatedValue;
				std::cout << "   p_soft = [" <<  softPdf[0] << "," << softPdf[1] << "]";
				std::cout << " - pAcc = " << pAcc << std::endl;
			}

			// accept simulatedValue with probabilty from soft data
			if (randomValue<pAcc) {
				isAccepted=true;
			}
			i++;
		} while ((i<maxIterations)&(!isAccepted));
		if(_debugMode > 2 ) {
			std::cout << " simval = " << simulatedValue << " - randomValue = " << randomValue << std::endl;
			std::cout << "________________________________" << std::endl;
		}


	} else {
		// obtain conditional and generate a realization wihtout soft data
		_getCpdfTiEnesim(sgIdxX, sgIdxY,sgIdxZ, conditionalPdfFromTi);

		// sample from conditional pdf
		simulatedValue = _sampleFromPdf(conditionalPdfFromTi);
	}



	// DONE SIMULATING
	if (_debugMode > 2) {
		// std::cout << "Simulated value="<< simulatedValue << ", LC_dist_min=" << LC_dist_min << " " << std::endl ;
		std::cout << "Simulated value="<< simulatedValue << " " << std::endl ;
	}
	return simulatedValue;
}



/**
* @brief get realization from Cpdf and Soft data using a metropolis sampler when soft data is available
* @param sgIdxX coordinate X of the current node
* @param sgIdxY coordinate Y of the current node
* @param sgIdxZ coordinate Z of the current node
* @param iterationCnt Iterations counter
* @return simulated value
*/
float MPS::ENESIM::_getRealizationFromCpdfTiEnesimRejectionNonCo(const int& sgIdxX, const int& sgIdxY, const int& sgIdxZ, float& iterationCnt) {


	std::map<float, float> conditionalPdfFromTi;
	// Check if any SoftData are available?
	//std::multimap<float, float> softPdf;
	float simulatedValue;
	float randomValue;
	float pAcc;
	float SoftProbability;
	float pAccOpt = 1e-9;
	float simulatedValueOpt;


	if (_nMaxCountCpdf == 1) {
		// REJECTION TYPE SIMULATION
		
		int i = 0;;
		int maxIterations = 200;  // decide whether soft or ti cpdf takes preference..
		bool isAccepted = false;


		do {
			// MAKE SURE THAT conditionalPdfFromTi IS OBTAINED INDEPENDENTLTY EACH TIME!!!
			// ESPECIELLY WHEN BASED ON ONLY ONE COUNT
			
			// obtain conditional and generate a realization wihtout soft data
			conditionalPdfFromTi.clear();
			SoftProbability = 1;

			_getCpdfTiEnesimNew(sgIdxX, sgIdxY, sgIdxZ, conditionalPdfFromTi, SoftProbability);
			simulatedValue = _sampleFromPdf(conditionalPdfFromTi);
			randomValue = ((float)rand() / (RAND_MAX));
			pAcc = SoftProbability;
				
			// Make sure to initialize optimal values
			if (i == 0) {
				pAccOpt = pAcc;
				simulatedValueOpt = simulatedValue;
			}

			// Store the most probable so far
			if (pAcc > pAccOpt) {
				simulatedValueOpt = simulatedValue;
				pAccOpt = pAcc;
			}
			
			if (_debugMode > 2) {
				std::cout << "  i= " << i << " maxIterations=" << maxIterations;
				std::cout << "   p_ti = [" << conditionalPdfFromTi[0] << "," << conditionalPdfFromTi[1] << "]";
				std::cout << " simval = " << simulatedValue;
				std::cout << "   pAcc =p_soft=" << SoftProbability;
				std::cout << " simvalOpt = " << simulatedValueOpt;
				std::cout << "   pAccOpt=" << pAccOpt << std::endl;
			}

			// accept simulatedValue with probabilty from soft data
			if (randomValue<pAcc) {
				isAccepted = true;
				simulatedValueOpt = simulatedValue;
				pAccOpt = pAcc;
				//std::cout << "  current hard data ACCEPTED  as " << randomValue << "<" << pAcc << std::endl;
			}
			else {
				//std::cout << "  current hard data real rejected as " << randomValue << "!<" << pAcc << std::endl;
			}
			i++;

		} while ((i<maxIterations)&(!isAccepted));
		simulatedValue = simulatedValueOpt;

	}
	else {
		// SIMULTATE DIRECTLY FROM CONDITIONAL
		// obtain conditional and generate a realization wihtout soft data
		_getCpdfTiEnesimNew(sgIdxX, sgIdxY, sgIdxZ, conditionalPdfFromTi, SoftProbability);
		simulatedValue = _sampleFromPdf(conditionalPdfFromTi);

	}

	if (_debugMode >1) {
		std::cout << "at [ix,y,ix]=[" << sgIdxX << "," << sgIdxY << "," << sgIdxZ << "] real=" << simulatedValue << std::endl;
		// std::cout << "Simulated value="<< simulatedValue << ", LC_dist_min=" << LC_dist_min << " " << std::endl ;
		// std::cout << "Simulated value=" << simulatedValue << " " << std::endl;
	}

	// DONE SIMULATING
	
	// At this point set the local soft data as NaN to avoid simulate it again
	if (!_softDataGrids.empty()) {
		_softDataGrids[0][sgIdxZ][sgIdxY][sgIdxX] = std::numeric_limits<float>::quiet_NaN();
		// Check if there are any more soft data left. If so, empty the _softDataGrid!
	}


	return simulatedValue;
}


