#pragma once
#include "glm/glm.hpp"
#include "glm/gtc/constants.hpp"
namespace SHFunctions {


	// order 6, by
	/*@article{Sloan2013SH,
	author = { Peter - Pike Sloan },
		title = { Efficient Spherical Harmonic Evaluation },
		year = { 2013 },
		month = { September },
		day = { 8 },
		journal = { Journal of Computer Graphics Techniques(JCGT) },
		volume = { 2 },
		number = { 2 },
		pages = { 84--83 },
		url = { http://jcgt.org/published/0002/02/06/},
		issn = {2331 - 7418}
	}*/

	void SHEval6(const float fX, const float fY, const float fZ, float *pSH)
	{
		float fC0, fC1, fS0, fS1, fTmpA, fTmpB, fTmpC;
		float fZ2 = fZ * fZ;

		pSH[0] = 0.2820947917738781f;
		pSH[2] = 0.4886025119029199f*fZ;
		pSH[6] = 0.9461746957575601f*fZ2 + -0.31539156525252f;
		pSH[12] = fZ * (1.865881662950577f*fZ2 + -1.119528997770346f);
		pSH[20] = 1.984313483298443f*fZ*pSH[12] + -1.006230589874905f*pSH[6];
		pSH[30] = 1.98997487421324f*fZ*pSH[20] + -1.002853072844814f*pSH[12];
		fC0 = fX;
		fS0 = fY;

		fTmpA = -0.48860251190292f;
		pSH[3] = fTmpA * fC0;
		pSH[1] = fTmpA * fS0;
		fTmpB = -1.092548430592079f*fZ;
		pSH[7] = fTmpB * fC0;
		pSH[5] = fTmpB * fS0;
		fTmpC = -2.285228997322329f*fZ2 + 0.4570457994644658f;
		pSH[13] = fTmpC * fC0;
		pSH[11] = fTmpC * fS0;
		fTmpA = fZ * (-4.683325804901025f*fZ2 + 2.007139630671868f);
		pSH[21] = fTmpA * fC0;
		pSH[19] = fTmpA * fS0;
		fTmpB = 2.03100960115899f*fZ*fTmpA + -0.991031208965115f*fTmpC;
		pSH[31] = fTmpB * fC0;
		pSH[29] = fTmpB * fS0;
		fC1 = fX * fC0 - fY * fS0;
		fS1 = fX * fS0 + fY * fC0;

		fTmpA = 0.5462742152960395f;
		pSH[8] = fTmpA * fC1;
		pSH[4] = fTmpA * fS1;
		fTmpB = 1.445305721320277f*fZ;
		pSH[14] = fTmpB * fC1;
		pSH[10] = fTmpB * fS1;
		fTmpC = 3.31161143515146f*fZ2 + -0.47308734787878f;
		pSH[22] = fTmpC * fC1;
		pSH[18] = fTmpC * fS1;
		fTmpA = fZ * (7.190305177459987f*fZ2 + -2.396768392486662f);
		pSH[32] = fTmpA * fC1;
		pSH[28] = fTmpA * fS1;
		fC0 = fX * fC1 - fY * fS1;
		fS0 = fX * fS1 + fY * fC1;

		fTmpA = -0.5900435899266435f;
		pSH[15] = fTmpA * fC0;
		pSH[9] = fTmpA * fS0;
		fTmpB = -1.770130769779931f*fZ;
		pSH[23] = fTmpB * fC0;
		pSH[17] = fTmpB * fS0;
		fTmpC = -4.403144694917254f*fZ2 + 0.4892382994352505f;
		pSH[33] = fTmpC * fC0;
		pSH[27] = fTmpC * fS0;
		fC1 = fX * fC0 - fY * fS0;
		fS1 = fX * fS0 + fY * fC0;

		fTmpA = 0.6258357354491763f;
		pSH[24] = fTmpA * fC1;
		pSH[16] = fTmpA * fS1;
		fTmpB = 2.075662314881041f*fZ;
		pSH[34] = fTmpB * fC1;
		pSH[26] = fTmpB * fS1;
		fC0 = fX * fC1 - fY * fS1;
		fS0 = fX * fS1 + fY * fC1;

		fTmpC = -0.6563820568401703f;
		pSH[35] = fTmpC * fC0;
		pSH[25] = fTmpC * fS0;
	}

}