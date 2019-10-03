#include <cmath> 
#include <cstdio> 
#include <random> 
#include <functional> 
#include <iostream> 
#include <fstream> 
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

template<typename T = float>
inline T linearInterpolation(const T &lo, const T &hi, const T &t)
{
	return lo * (1 - t) + hi * t;
}

class PerlinNoise
{
	static const unsigned tableSize = 256;
	static const unsigned tableSizeMask = tableSize - 1;
	glm::vec3 gradientTable[tableSize];
	unsigned permutationTable[tableSize * 2];

public:
	PerlinNoise()
	{
		unsigned seed = 2019;
		std::mt19937 generator(seed);
		std::uniform_real_distribution<float> distrFloat;
		auto dice = std::bind(distrFloat, generator);
		float gradientLen2;
		for (unsigned i = 0; i < tableSize; ++i) {
			do {
				gradientTable[i] = glm::vec3(2 * dice() - 1, 2 * dice() - 1, 2 * dice() - 1);
				gradientLen2 = gradientTable[i].x * gradientTable[i].x + gradientTable[i].y * gradientTable[i].y + gradientTable[i].z * gradientTable[i].z;
			} while (gradientLen2 > 1);
			gradientTable[i] /= sqrtf(gradientLen2); // normalize gradient 
			permutationTable[i] = i;
		}

		std::uniform_int_distribution<unsigned> distributionInt;
		auto diceInt = std::bind(distributionInt, generator);
		for (unsigned i = 0; i < tableSize; ++i)
			std::swap(permutationTable[i], permutationTable[diceInt() & tableSizeMask]);
		
		// extend the permutation table in the index range [256:512]
		for (unsigned i = 0; i < tableSize; ++i) {
			permutationTable[tableSize + i] = permutationTable[i];
		}
	}
	virtual ~PerlinNoise() {}

	int hash(const int &x, const int &y, const int &z) const
	{
		return permutationTable[permutationTable[permutationTable[x] + y] + z];
	}

public:
	float eval(const glm::vec3 &p) const
	{
		int xi0 = ((int)std::floor(p.x)) & tableSizeMask;
		int yi0 = ((int)std::floor(p.y)) & tableSizeMask;
		int zi0 = ((int)std::floor(p.z)) & tableSizeMask;

		int xi1 = (xi0 + 1) & tableSizeMask;
		int yi1 = (yi0 + 1) & tableSizeMask;
		int zi1 = (zi0 + 1) & tableSizeMask;

		float tx = p.x - ((int)std::floor(p.x));
		float ty = p.y - ((int)std::floor(p.y));
		float tz = p.z - ((int)std::floor(p.z));

		float u = tx; 
		float v = ty;  
		float w = tz; 

		// gradientTable at the corner of the cell
		const glm::vec3 &c000 = gradientTable[hash(xi0, yi0, zi0)];
		const glm::vec3 &c100 = gradientTable[hash(xi1, yi0, zi0)];
		const glm::vec3 &c010 = gradientTable[hash(xi0, yi1, zi0)];
		const glm::vec3 &c110 = gradientTable[hash(xi1, yi1, zi0)];

		const glm::vec3 &c001 = gradientTable[hash(xi0, yi0, zi1)];
		const glm::vec3 &c101 = gradientTable[hash(xi1, yi0, zi1)];
		const glm::vec3 &c011 = gradientTable[hash(xi0, yi1, zi1)];
		const glm::vec3 &c111 = gradientTable[hash(xi1, yi1, zi1)];

		// generate vectors going from the grid points to p
		float x0 = tx, x1 = tx - 1;
		float y0 = ty, y1 = ty - 1;
		float z0 = tz, z1 = tz - 1;

		glm::vec3 p000 = glm::vec3(x0, y0, z0);
		glm::vec3 p100 = glm::vec3(x1, y0, z0);
		glm::vec3 p010 = glm::vec3(x0, y1, z0);
		glm::vec3 p110 = glm::vec3(x1, y1, z0);

		glm::vec3 p001 = glm::vec3(x0, y0, z1);
		glm::vec3 p101 = glm::vec3(x1, y0, z1);
		glm::vec3 p011 = glm::vec3(x0, y1, z1);
		glm::vec3 p111 = glm::vec3(x1, y1, z1);

		// linear interpolation
		float a = linearInterpolation(glm::dot(c000, p000), glm::dot(c100, p100), u);
		float b = linearInterpolation(glm::dot(c010, p010), glm::dot(c110, p110), u);
		float c = linearInterpolation(glm::dot(c001, p001), glm::dot(c101, p101), u);
		float d = linearInterpolation(glm::dot(c011, p011), glm::dot(c111, p111), u);

		float e = linearInterpolation(a, b, v);
		float f = linearInterpolation(c, d, v);

		return linearInterpolation(e, f, w); 
	}
};

int main(int argc, char **argv)
{
	unsigned imageWidth = 512;
	unsigned imageHeight = 512;
	float *noiseMap = new float[imageWidth * imageHeight]{ 0 };

	PerlinNoise noise;
	glm::vec3 vector = glm::vec3(0, 0, 0);
	float frequency = 0.01f;
	for (unsigned j = 0; j < imageHeight; ++j) {
		for (unsigned i = 0; i < imageWidth; ++i) {
			vector.x = i * frequency * 10;
			vector.y = j * frequency * 10;

			float g = noise.eval(vector) ;
			noiseMap[j * imageWidth + i] = g - (int)g;
		}
	}

	// output noise map to PPM
	std::ofstream ofs;
	ofs.open("./heightMap.ppm", std::ios::out | std::ios::binary);
	ofs << "P6\n" << imageWidth << " " << imageHeight << "\n255\n";
	for (unsigned k = 0; k < imageWidth * imageHeight; ++k) {
		unsigned char n = static_cast<unsigned char>(noiseMap[k] * 255);
		ofs << n << n << n;
	}
	ofs.close();

	delete[] noiseMap;

	return 0;
}