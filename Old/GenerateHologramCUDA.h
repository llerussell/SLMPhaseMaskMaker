/*   
   Header file for GenerateHologramCUDA
   Hologram generating algorithms for CUDA Devices
   
   Copyright 2009, 2010, 2011, 2012 Martin Persson 
   martin.persson@physics.gu.se


   This file is part of GenerateHologramCUDA.

    GenerateHologramCUDA is free software: you can redistribute it and/or 
    modify it under the terms of the GNU Lesser General Public License as published 
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GenerateHologramCUDA is distributed in the hope that it will be 
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with GenerateHologramCUDA.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GENERATEHOLOGRAMCUDA_H
#define GENERATEHOLOGRAMCUDA_H
    
#ifdef __cplusplus
extern "C" {  // only need to export C interface if
              // used by C++ source code
#endif

__declspec(dllexport) int GenerateHologram(float *h_checkData, unsigned char *h_pSLM_uc, float *x_spots, float *y_spots, float *z_spots, float *I_spots, int N_spots, int N_iterations, float *h_Iobtained, int method);
__declspec(dllexport) int Corrections(int UseAberrationCorr, float *h_AberrationCorr, int UseLUTPol, int PolOrder, float *h_LUTPolCoeff, int saveAmplitudes, float alpha, int DCborderWidth, int UseLUT, unsigned char *h_LUT_uc);
__declspec(dllexport) int startCUDAandSLM(int EnableSLM, float *h_pSLMstart, char* LUTFile, unsigned short TrueFrames, int deviceId);
__declspec(dllexport) int stopCUDAandSLM();

#ifdef __cplusplus
}
#endif

#endif