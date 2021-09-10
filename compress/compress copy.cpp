/*
****************************

FasTensor (FT) Copyright (c) 2021, The Regents of the University of
California, through Lawrence Berkeley National Laboratory (subject to
receipt of any required approvals from the U.S. Dept. of Energy). 
All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department
of Energy and the U.S. Government consequently retains certain rights.  As
such, the U.S. Government has been granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative 
works, and perform publicly and display publicly, and to permit others to do so.


****************************


*** License Agreement ***

FasTensor (FT) Copyright (c) 2021, The Regents of the University of
California, through Lawrence Berkeley National Laboratory (subject to
receipt of any required approvals from the U.S. Dept. of Energy). 
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

(2) Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

(3) Neither the name of the University of California, Lawrence Berkeley
National Laboratory, U.S. Dept. of Energy nor the names of its contributors
may be used to endorse or promote products derived from this software
without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches,
or upgrades to the features, functionality or performance of the source
code ("Enhancements") to anyone; however, if you choose to make your
Enhancements available either publicly, or directly to Lawrence Berkeley
National Laboratory, without imposing a separate written license agreement
for such Enhancements, then you hereby grant the following license: a
non-exclusive, royalty-free perpetual license to install, use, modify,
prepare derivative works, incorporate into other computer software,
distribute, and sublicense such enhancements or derivative works thereof,
in binary and source code form.
*/

/**
 * @example ft_example_array.cpp
 * @author Bin Dong (dbin@lbl.gov)
 * @brief This example show how to use Array in FasTensor to access data.  
 * @version 0.1
 * @date 2021-05-28
 * 
 * @copyright Copyright (c) 2021
 * 
 */

/**
 *
 * Author: Bin Dong dbin@lbl.gov
 * Web: https://crd.lbl.gov/bin-dong
 * Scientific Data Management Research Group
 * Lawrence Berkeley National Laboratory
 *
 */

#include <iostream>
#include <stdarg.h>
#include <vector>
#include <stdlib.h>
#include "ft.h"

using namespace std;
using namespace FT;

int main(int argc, char *argv[])
{
	//Init the MPICH, etc.
	AU_Init(argc, argv);

	// set up the chunk size and the overlap size
	// 11648, 30000 for each dataset
	std::vector<int> chunk_size(2);
	std::vector<int> overlap_size = {0, 0};

	//Input data
	Array<short> *A = new Array<float>("EP_HDF5:./test-data/testf-16x16-array.h5:/testg/testd");
	std::vector<unsigned long long> array_size;
	A->GetArraySize(array_size);
	chunk_size[0] = array_size[0];
	chunk_size[1] = array_size[1];

	std::vector<unsigned long long> start = {0, 0};
	std::vector<unsigned long long> end = {5, 5};
	std::vector<short> data;

	data.resize(array_size[0] * array_size[1]);
	A->SetChunkSize(chunk_size);
	A->SetOverlapSize(overlap_size);

	end[0] = array_size[0];
	end[1] = array_size[1];
	A->ReadEndpoint(start, end, static_cast<void *>(data.data()));

	//Input data
	Array<short> *B = new Array<float>("EP_HDF5:./test-data/testf-16x16-array-sz.h5:/testg/testd");
	std::vector<std::string> filter_paramter;
	//#define H5Z_FILTER_SZ 32017
	int filter_id = 32017;

	filter_paramter.push_back(std::to_string(filter_id));
	filter_paramter.push_back(std::to_string(array_size[0] * array_size[1]));

	std::vector<unsigned int> filter_cd_values;
	filter_paramter.push_back(Vector2String(filter_cd_values)));

	std::vector<hsize_t> filter_chunk_size;
	filter_chunk_size.push_back(array_size[0]);
	filter_chunk_size.push_back(array_size[1]);
	filter_paramter.push_back(Vector2String(filter_chunk_size)));
	//  parameter_v[0]: filter_id
	// parameter_v[1]: filter_cd_nelmts
	//  parameter_v[2]: filter_cd_values
	//  parameter_v[3]: chunk_size
	B->ControlEndpoint(HDF5_ENABLE_FILTER, filter_paramter);
	B->WriteEndpoint(start, end, static_cast<void *>(data.data()));

	delete A;

	AU_Finalize();

	return 0;
}
