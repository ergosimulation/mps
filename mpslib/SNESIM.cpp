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
#include <cctype>       // isspace
#include <algorithm>    // std::random_shuffle std::remove_if
#include <list>

#include "SNESIM.h"
#include "IO.h"
#include "Coords3D.h"

/**
* @brief Constructors
*/
MPS::SNESIM::SNESIM(void) : MPS::MPSAlgorithm(){

}

/**
* @brief Destructors
*/
MPS::SNESIM::~SNESIM(void) {

}

/**
* @brief Read configuration file
* @param fileName configuration filename
*/
void MPS::SNESIM::_readConfigurations(const std::string& fileName) {

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
	// Process each line
	// Number of realizations
	_readLineConfiguration(file, ss, data, s, str);
	_realizationNumbers = stoi(data[1]);
	// Initial value
	_readLineConfiguration(file, ss, data, s, str);
	_seed = stof(data[1]);
	// Multiple grid level
	_readLineConfiguration(file, ss, data, s, str);
	_totalGridsLevel = stoi(data[1]);
	//Min node count
	_readLineConfiguration(file, ss, data, s, str);
	_minNodeCount = stoi(data[1]);
	//Max conditional count
	_readLineConfiguration(file, ss, data, s, str);
	_maxCondData = stoi(data[1]);
	
	// Template size x
	_readLineConfiguration_mul(file, ss, data, s, str);
	_templateSizeX = stoi(data[1]);
	// optional. Template size x - base when n_mulgrid=0
	if (data.size()>2) {
		_templateSizeX_base = stoi(data[2]);
	} else {
		_templateSizeX_base = _templateSizeX;
	}
  // Template size y
	_readLineConfiguration_mul(file, ss, data, s, str);
	_templateSizeY = stoi(data[1]);
	// optional. Template size y - base when n_mulgrid=0
	if (data.size()>2) {
		_templateSizeY_base = stoi(data[2]);
	} else {
		_templateSizeY_base = _templateSizeY;
	}
  // Template size z
	_readLineConfiguration_mul(file, ss, data, s, str);
	_templateSizeZ = stoi(data[1]);
	// optional. Template size z - base when n_mulgrid=0
	if (data.size()>2) {
		_templateSizeZ_base = stoi(data[2]);
	} else {
		_templateSizeZ_base = _templateSizeZ;
	}
	
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
	data[1].erase(std::remove_if(data[1].begin(), data[1].end(), [](char x){return std::isspace(x);}), data[1].end()); //Removing spaces
	_tiFilename = data[1];
	// Output directory
	_readLineConfiguration(file, ss, data, s, str);
	data[1].erase(std::remove_if(data[1].begin(), data[1].end(), [](char x){return std::isspace(x);}), data[1].end()); //Removing spaces
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
	data[1].erase(std::remove_if(data[1].begin(), data[1].end(), [](char x){return std::isspace(x);}), data[1].end()); //Removing spaces
	_hardDataFileNames = data[1];
	// Hard data search radius
	_readLineConfiguration(file, ss, data, s, str);
	_hdSearchRadius = stof(data[1]);
	//Get the maximum value of template size into the search radius
	_hdSearchRadius = std::max(std::max(_templateSizeX, _templateSizeY), _templateSizeZ);
	// Softdata categories
	//std::cout << s << std::endl;
	if(_readLineConfiguration(file, ss, data, s, str)) {
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
				s.erase(std::remove_if(s.begin(), s.end(), [](char x){return std::isspace(x);}), s.end()); //Removing spaces
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
	// Mask data grid
	if (_readLineConfiguration(file, ss, data, s, str)) {
		data[1].erase(std::remove_if(data[1].begin(), data[1].end(), [](char x) {return std::isspace(x); }), data[1].end()); //Removing spaces
		_maskDataFileName = data[1];
	}


}

/**
* @brief Construct templatefaces and sort them around template center
* @param sizeX template size X
* @param sizeX template size Y
* @param sizeX template size Z
*/
void MPS::SNESIM::_constructTemplateFaces(const int& sizeX, const int& sizeY, const int& sizeZ) {
	int templateCenterX = (int)floor(sizeX / 2);
	int templateCenterY = (int)floor(sizeY / 2);
	int templateCenterZ = (int)floor(sizeZ / 2);
	int totalTemplateIndices = sizeX * sizeY * sizeZ;
	int totalTemplates = 0;
	//std::cout << "template center (x, y, z): " << templateCenterX << " " << templateCenterY << " " << templateCenterZ << std::endl;

	//Create a template path
	std::vector<int> templatePath;
	templatePath.resize(totalTemplateIndices);
	//Initialize a sequential order
	_initilizePath(sizeX, sizeY, sizeZ, templatePath);

	//Sort template path with the distance to the template center
	std::sort(templatePath.begin(), templatePath.end(), TemplateSorter(sizeX, sizeY, sizeZ));

	//Loop through all template indices
	//initialize faces
	_templateFaces.clear();
	//Initialize the faces with center point(0, 0, 0)
	_templateFaces.push_back(MPS::Coords3D(0, 0, 0));
	int offsetX, offsetY, offsetZ;
	int templateIdxX, templateIdxY, templateIdxZ;
	for (int i=0; i<totalTemplateIndices; i++) {
		//std::cout << templatePath[i] << " ";
		MPS::utility::oneDTo3D(templatePath[i], sizeX, sizeY, templateIdxX, templateIdxY, templateIdxZ);
		offsetX = templateIdxX - templateCenterX;
		offsetY = templateIdxY - templateCenterY;
		offsetZ = templateIdxZ - templateCenterZ;
		//Ignore center point
		if (offsetX != 0 || offsetY != 0 || offsetZ != 0) {
			_templateFaces.push_back(MPS::Coords3D(offsetX, offsetY, offsetZ));
		}
	}

	////Forcing a template like in article
	//_templateFaces.clear();
	//_templateFaces.push_back(MPS::Coords3D(0, 0, 0));
	//_templateFaces.push_back(MPS::Coords3D(0, 1, 0));
	//_templateFaces.push_back(MPS::Coords3D(1, 0, 0));
	//_templateFaces.push_back(MPS::Coords3D(0, -1, 0));
	//_templateFaces.push_back(MPS::Coords3D(-1, 0, 0));
	////Showing the template faces
	//std::cout << _templateFaces.size() << std::endl;
	//for (unsigned int i=0; i<_templateFaces.size(); i++) {
	//	std::cout << _templateFaces[i].getX() << " " << _templateFaces[i].getY() << " " << _templateFaces[i].getZ() << std::endl;
	//}

}

/**
* @brief Getting a node value by calculating the cummulatie probability distribution function
* @param conditionalPoints list of all found conditional points
* @param x coordinate X of the current node
* @param y coordinate Y of the current node
* @param z coordinate Z of the current node
*/
float MPS::SNESIM::_cpdf(std::map<float, int>& conditionalPoints, const int& x, const int& y, const int& z) {
	// TMH: RENAME TO _cpdfFromConditionalCount
	//      RENAME conditionalPoints to conditionalCount

	//Fill probability from TI
	float foundValue = std::numeric_limits<float>::quiet_NaN();

	// Loop through all the found conditional data values, and find the number of occurences
	// as the sum of the occurences of each conditional data value.
	float totalCounter = 0;
	for(std::map<float,int>::iterator iter = conditionalPoints.begin(); iter != conditionalPoints.end(); ++iter) {
		totalCounter += iter->second;
		//std::cout << "cpdf ti1: " << iter->first << " " << iter->second << std::endl;
	}
	//std::cout << "totalCounter " << totalCounter << std::endl;

	// obtain the conditional probability from the training image
	// by dividing the number of counts for each condtional data value, with the
	// number of count (totalCounter)
	std::map<float, float> probabilitiesFromTI;
	for(std::map<float, int>::iterator iter = conditionalPoints.begin(); iter != conditionalPoints.end(); ++iter) {
		probabilitiesFromTI.insert(std::pair<float, float>(iter->first, (float)(iter->second) / (float)totalCounter));
		// std::cout << "cpdf ti2: " << iter->first << " " << (float)(iter->second) / (float)totalCounter << std::endl;
	}

	// Perhaps some soft data needs to be combined with the
	// conditional probability found from the previously simulated data.
	std::multimap<float, float> probabilitiesCombined;

	//Fill probability from Softdata
	bool useSoftData = true;
	//Empty grids check
	if (_softDataGrids.empty()) useSoftData = false;
	//No value check
	else if (MPS::utility::is_nan(_softDataGrids[0][z][y][x])) useSoftData = false;
	//Not same size check
	//else if (_softDataCategories.size() != probabilitiesFromTI.size()) useSoftData = false;

	//Perform computation of probabilities from soft data
	if (useSoftData) {
		std::map<float, float> probabilitiesFromSoftData;
		float sumProbability = 0;
		unsigned int lastIndex = (int)_softDataCategories.size() - 1;
		//_softDataGrids[0][z][y][x] = 0.8;
		for (unsigned int i=0; i<lastIndex; i++) {
			sumProbability += _softDataGrids[i][z][y][x];
			probabilitiesFromSoftData.insert(std::pair<float, float>(_softDataCategories[i], _softDataGrids[i][z][y][x]));
			//std::cout << "cpdf sd: " << _softDataCategories[i] << " " << _softDataGrids[i][z][y][x] << std::endl;
		}
		//Last categorie
		probabilitiesFromSoftData.insert(std::pair<float, float>(_softDataCategories[lastIndex], 1 - sumProbability));
		//std::cout << "cpdf sd: " << _softDataCategories[lastIndex] << " " << 1 - sumProbability << std::endl;
		//Compute the combined probabilities from TI and Softdata
		
		// totalValue is the sum of the prodict fo P_TI and P_SOFT, IF THIS IS ZERO, THE the two are conflict
		// Then rge soft data wins..
		float totalValue = 0;
		std::map<float,float>::iterator searchIter;
		for(std::map<float,float>::iterator iter = probabilitiesFromSoftData.begin(); iter != probabilitiesFromSoftData.end(); ++iter) {
			//Size check
			searchIter = probabilitiesFromTI.find(iter->first);
			if (searchIter != probabilitiesFromTI.end()) totalValue += iter->second * probabilitiesFromTI[iter->first];
			else totalValue += 0; //iter->second; //If probability from TI has less elements than sd then ignore it
		}


		//Normalize and compute cummulated value
		float cumulateValue = 0;
		for(std::map<float,float>::iterator iter = probabilitiesFromSoftData.begin(); iter != probabilitiesFromSoftData.end(); ++iter) {
			//Size check
			searchIter = probabilitiesFromTI.find(iter->first);
			if (searchIter != probabilitiesFromTI.end()) {
				//std::cout << "P_soft=" << iter->second << " P_TI=" << probabilitiesFromTI[iter->first] << " TOT=" << totalValue << std::endl;
				cumulateValue += (iter->second * probabilitiesFromTI[iter->first]) / totalValue;
			} else {
				cumulateValue += 0; //(iter->second) / totalValue;  //If probability from TI has less elements than sd then ignore it
			}
			//std::cout << "  cumulateValue=" << cumulateValue << std::endl;

			probabilitiesCombined.insert(std::pair<float, float>(cumulateValue, iter->first));
		}

		if (totalValue == 0) {
			// This means there is inconsistentcy between probabilitiesFromTI and probabilitiesFromSoftData
			// In this case the 'data' takes precedence of the conditional pdf from TI
			probabilitiesCombined.clear();
			probabilitiesCombined.insert(probabilitiesFromSoftData.begin(), probabilitiesFromSoftData.end());
		}


		if (_debugMode > 2) {
			std::cout << "totalCounter " << totalCounter << std::endl;
			// PRINT SOFT CONDITIONAL AND LCOAL SOFT PROB
			std::cout << "TI  : ";
			for (std::map<float, float>::iterator iter = probabilitiesFromTI.begin(); iter != probabilitiesFromTI.end(); ++iter) {
				std::cout << " P[" << iter->first << "]=";
				std::cout << iter->second;
			}
			std::cout << std::endl;
			std::cout << "SOFT: ";
			for (std::map<float, float>::iterator iter = probabilitiesFromSoftData.begin(); iter != probabilitiesFromSoftData.end(); ++iter) {
				std::cout << " P[" << iter->first << "]=";
				std::cout << iter->second;
			}
			std::cout << std::endl;
			std::cout << "COMB: ";
			for (std::map<float, float>::iterator iter = probabilitiesCombined.begin(); iter != probabilitiesCombined.end(); ++iter) {
				std::cout << " P[" << iter->first << "]=";
				std::cout << iter->second;
			}
			std::cout << std::endl;
			std::cout << std::endl;
		}



	} else {
		//Only from TI
		float cumulateValue = 0;
		for(std::map<float,float>::iterator iter = probabilitiesFromTI.begin(); iter != probabilitiesFromTI.end(); ++iter) {
			cumulateValue += iter->second;
			probabilitiesCombined.insert(std::pair<float, float>(cumulateValue, iter->first));
		}
	}

	//Getting the probability distribution value
	//Random possible value between 0 and 1 then scale into the maximum probability
	float randomValue = ((float) rand() / (RAND_MAX));
	//randomValue = 0.8;
	for(std::map<float,float>::iterator iter = probabilitiesCombined.begin(); iter != probabilitiesCombined.end(); ++iter) {
		if (iter->first >= randomValue) {
			foundValue = iter->second;
			break;
		}
	}
	//std::cout << randomValue << " " << foundValue << std::endl;
	return foundValue;
}
