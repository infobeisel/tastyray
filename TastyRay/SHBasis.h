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

		// fz d: (14.809976568*fZ2*fZ - 6.347132815*fZ) 
		pSH[20] = 1.984313483298443f*fZ*pSH[12] + -1.006230589874905f*pSH[6];
		
		 //fz d: (36.839351573*fZ2*fZ2 - 24.559567715*fZ2 + 1.754254837)  
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

	void SHEval6Derivative(const float fX, const float fY, const float fZ, float * pSH, float * pSHdX, float * pSHdY, float * pSHdZ) {

		float fC0, fC1, fS0, fS1, fTmpA, fTmpB, fTmpC, fTmpAdZ, fTmpBdZ, fTmpCdZ;
		float fC0dX, fC1dX, fS0dX, fS1dX, fC0dY, fC1dY, fS0dY, fS1dY;
		float fZ2 = fZ * fZ;
		float fX2 = fX * fX;
		float fY2 = fY * fY;
		for(int i = 0; i < 36; i++) {
			pSHdX[i] = pSHdY[i] = pSHdZ[i] = 0.0f;
		}
		//pSHdX[0] = pSHdX[2] = pSHdX[6] = pSHdX[12] = pSHdX[20] = pSHdX[30] = 0.0f;
		//pSHdY[0] = pSHdY[2] = pSHdY[6] = pSHdY[12] = pSHdY[20] = pSHdY[30] = 0.0f;
		//pSHdZ[0] = pSHdZ[3] = pSHdZ[1] = pSHdZ[4] = pSHdZ[8] = pSHdZ[15] = pSHdZ[9] = pSHdZ[24] = pSHdZ[16] = pSHdZ[35] = pSHdZ[25] = 0.0f; 
		pSH[0] = 0.2820947917738781f;
		pSH[2] = 0.4886025119029199f*fZ;
		pSHdZ[2] = 0.4886025119029199f;
		pSH[6] = 0.9461746957575601f*fZ2 + -0.31539156525252f;
		pSHdZ[6] = 2.0f * 0.9461746957575601f * fZ;
		pSH[12] = fZ * (1.865881662950577f*fZ2 + -1.119528997770346f);
		pSHdZ[12] = 3.0f * 1.865881662950577f*fZ2 - 1.119528997770346f;
		pSH[20] = 1.984313483298443f*fZ*pSH[12] + -1.006230589874905f*pSH[6];
		pSHdZ[20] = (14.809976568f*fZ2*fZ - 6.347132815f*fZ) ;
		pSH[30] = 1.98997487421324f*fZ*pSH[20] + -1.002853072844814f*pSH[12];
		pSHdZ[30] = (36.839351573f*fZ2*fZ2 - 24.559567715f*fZ2 + 1.754254837f);
		fC0 = fX;
		fS0 = fY;

		fTmpA = -0.48860251190292f;
		pSH[3] = fTmpA * fC0;
		pSHdX[3] = fTmpA;
		pSHdY[3] = 0.0f;
		pSH[1] = fTmpA * fS0;
		pSHdX[1] = 0.0f;
		pSHdY[1] = fTmpA;
		fTmpB = -1.092548430592079f*fZ;
		fTmpBdZ = -1.092548430592079f;
		pSH[7] = fTmpB * fC0;
		pSHdX[7] = fTmpB;
		pSHdY[7] = 0.0f;
		pSHdZ[7] = fTmpBdZ * fC0;
		pSH[5] = fTmpB * fS0;
		pSHdX[5] = 0.0f;
		pSHdY[5] = fTmpB;
		pSHdZ[5] = fTmpBdZ * fS0;
		fTmpC = -2.285228997322329f*fZ2 + 0.4570457994644658f;
		fTmpCdZ = -2.285228997322329f*2.0f*fZ;
		pSH[13] = fTmpC * fC0;
		pSHdX[13] = fTmpC;
		pSHdY[13] = 0.0f;
		pSHdZ[13] = fTmpCdZ * fC0;
		pSH[11] = fTmpC * fS0;
		pSHdX[11] = 0.0f;
		pSHdY[11] = fTmpB;
		pSHdZ[11] = fTmpCdZ * fS0;
		fTmpA = fZ * (-4.683325804901025f*fZ2 + 2.007139630671868f);
		
		fTmpAdZ = (-14.049977415f*fZ2 + 2.007139630671868f);
		pSH[21] = fTmpA * fC0;
		pSHdX[21] = fTmpA;
		pSHdY[21] = 0.0f;
		pSHdZ[21] = fTmpAdZ * fC0;
		pSH[19] = fTmpA * fS0;
		pSHdX[19] = 0.0f;
		pSHdY[19] = fTmpA;
		pSHdZ[19] = fTmpAdZ * fC0;
		fTmpB = 2.03100960115899f*fZ*fTmpA + -0.991031208965115f*fTmpC;

		fTmpBdZ = -38.0475187*fZ2*fZ + 12.682506233 * fZ;
		pSH[31] = fTmpB * fC0;
		pSHdX[31] = fTmpB;
		pSHdY[31] = 0.0f;
		pSHdZ[31] = fTmpBdZ * fC0;
		pSH[29] = fTmpB * fS0;
		pSHdX[29] = 0.0f;
		pSHdY[29] = fTmpB;
		pSHdZ[29] = fTmpBdZ * fS0;
		fC1 = fX * fC0 - fY * fS0;
		fS1 = fX * fS0 + fY * fC0;


		fTmpA = 0.5462742152960395f;
		fC1dX = 2.0f * fX;
		fC1dY = -2.0f * fY;
		fS1dX = 2.0f * fY;
		fS1dY = 2.0f * fX;
		pSH[8] = fTmpA * fC1;
		pSHdX[8] = fTmpA * fC1dX;
		pSHdY[8] = fTmpA * fC1dY;
		pSH[4] = fTmpA * fS1;
		pSHdX[4] = fTmpA * fS1dX;
		pSHdY[4] = fTmpA * fS1dY;
		fTmpB = 1.445305721320277f*fZ;
		fTmpBdZ = 1.445305721320277f;
		pSH[14] = fTmpB * fC1;
		pSHdX[14] = fTmpB * fC1dX;
		pSHdY[14] = fTmpB * fC1dY;
		pSHdZ[14] = fTmpBdZ * fC1;
		pSH[10] = fTmpB * fS1;
		pSHdX[10] = fTmpB * fS1dX;
		pSHdY[10] = fTmpB * fS1dY;
		pSHdZ[10] = fTmpBdZ * fS1;
		fTmpC = 3.31161143515146f*fZ2 + -0.47308734787878f;
		fTmpCdZ = 2.0f*3.31161143515146f*fZ;
		pSH[22] = fTmpC * fC1;
		pSHdX[22] = fTmpC * fC1dX;
		pSHdY[22] = fTmpC * fC1dY;
		pSHdZ[22] = fTmpBdZ * fC1;
		pSH[18] = fTmpC * fS1;
		pSHdX[18] = fTmpC * fS1dX;
		pSHdY[18] = fTmpC * fS1dY;
		pSHdZ[18] = fTmpBdZ * fS1;
		fTmpA = fZ * (7.190305177459987f*fZ2 + -2.396768392486662f);
		fTmpAdZ = 21.570915532f*fZ2-2.396768392f;
		pSH[32] = fTmpA * fC1;
		pSHdX[32] = fTmpA * fC1dX;
		pSHdY[32] = fTmpA * fC1dY;
		pSHdZ[32] = fTmpAdZ * fC1;
		pSH[28] = fTmpA * fS1;
		pSHdX[28] = fTmpA * fS1dX;
		pSHdY[28] = fTmpA * fS1dY;
		pSHdZ[28] = fTmpAdZ * fS1;
		fC0 = fX * fC1 - fY * fS1;
		fS0 = fX * fS1 + fY * fC1;

		
		

		fTmpA = -0.5900435899266435f;
		fC0dX = 3.0f * fX2 - fY2 - 2.0f * fY2; 
		fC0dY = -6.0f * fY * fX;
		fS0dX = 6.0f * fX * fY;
		fS0dY = 2.0f * fX2 + fX2 - 3.0f * fY2;
		pSH[15] = fTmpA * fC0;
		pSHdX[15] = fTmpA * fC0dX;
		pSHdY[15] = fTmpA * fC0dY;
		pSH[9] = fTmpA * fS0;
		pSHdX[9] = fTmpA * fS0dX;
		pSHdY[9] = fTmpA * fS0dY;
		fTmpB = -1.770130769779931f*fZ;
		fTmpBdZ = -1.770130769779931f;
		pSH[23] = fTmpB * fC0;
		pSHdX[23] = fTmpB * fC0dX;
		pSHdY[23] = fTmpB * fC0dY;
		pSHdZ[23] = fTmpBdZ * fC0;
		pSH[17] = fTmpB * fS0;
		pSHdX[17] = fTmpB * fS0dX;
		pSHdY[17] = fTmpB * fS0dY;
		pSHdZ[17] = fTmpBdZ * fS0;
		fTmpC = -4.403144694917254f*fZ2 + 0.4892382994352505f;
		fTmpCdZ = -4.403144694917254f*2.0f*fZ;
		pSH[33] = fTmpC * fC0;
		pSHdX[33] = fTmpC * fC0dX;
		pSHdY[33] = fTmpC * fC0dY;
		pSHdZ[33] = fTmpCdZ * fC0;
		pSH[27] = fTmpC * fS0;
		pSHdX[27] = fTmpC * fS0dX;
		pSHdY[27] = fTmpC * fS0dY;
		pSHdZ[27] = fTmpCdZ * fS0;
		fC1 = fX * fC0 - fY * fS0;
		fS1 = fX * fS0 + fY * fC0;

		fTmpA = 0.6258357354491763f;
		fC1dX =  4.0f * glm::pow(fX,3.0f) - 12.0f * fX * fY2;
		fC1dY = -12.0f * fX2 * fY + 4.0f * fY2 * fY;
		fS1dX = 12.0f * fX2 * fY - 4.0f * fY2 * fY;
		fS1dY = -12.0f * fY2 * fX + 4.0f * fX2 * fX;
		pSH[24] = fTmpA * fC1;
		pSHdX[24] = fTmpA * fC1dX;
		pSHdY[24] = fTmpA * fC1dY;
		pSH[16] = fTmpA * fS1;
		pSHdX[16] = fTmpA * fS1dX;
		pSHdY[16] = fTmpA * fS1dY;
		fTmpB = 2.075662314881041f*fZ;
		fTmpBdZ = 2.075662314881041f;
		pSH[34] = fTmpB * fC1;
		pSHdX[34] = fTmpB * fC1dX;
		pSHdY[34] = fTmpB * fC1dY;
		pSHdZ[34] = fTmpBdZ * fC1;
		pSH[26] = fTmpB * fS1;
		pSHdX[26] = fTmpB * fS1dX;
		pSHdY[26] = fTmpB * fS1dY;
		pSHdZ[26] = fTmpBdZ * fS1;
		fC0 = fX * fC1 - fY * fS1;
		fS0 = fX * fS1 + fY * fC1;

		fC0dX = 5.0f * fX2 * fX2 - 30.0f * fX2 * fY2 + 5.0f * fY2 * fY2 ; 
		fC0dY = 20.0f * fX * fY2 * fY - 20.0f * fX2 * fX * fY;
		fS0dX = 20.0f * fX2 * fX * fY - 20.0f * fX * fY2 * fY;
		fS0dY = 5.0f * fX2 * fX2 - 30.0f * fX2 * fY2 + 5.0f * fY2 * fY2;
		fTmpC = -0.6563820568401703f;
		pSH[35] = fTmpC * fC0;
		pSHdX[35] = fTmpC * fC0dX;
		pSHdY[35] = fTmpC * fC0dY;
		pSH[25] = fTmpC * fS0;
		pSHdX[25] = fTmpC * fS0dX;
		pSHdY[25] = fTmpC * fS0dY;
	}

}