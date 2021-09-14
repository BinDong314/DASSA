/**
 *
 * Author: Bin Dong
 * Email questions to dbin@lbl.gov
 * Scientific Data Management Research Group
 * Lawrence Berkeley National Laboratory
 *
 */

//
// This  example code is for compressing the DAS data

#include <iostream>
#include <stdarg.h>
#include <vector>
#include <stdlib.h>
#include <string>

#include "ft.h"
#include "../Tools3rd/termcolor.hpp"
#include "../Tools3rd/INIReader.h"

using namespace std;
using namespace AU;
//using namespace DasLib;

int read_config_file(std::string file_name, int mpi_rank);
std::string config_file = "./decimate.config";

bool is_input_single_file = false;
std::string input_dir_file = "./test-data/dir";
std::string input_file_type = "EP_HDF5";
std::string input_h5_dataset = "/dat";

bool is_input_search_rgx = false;
std::string input_search_rgx = "^(.*)[1234]\\.tdms$";

int chs_per_file = 11648;
int lts_per_file = 30000;
int n_files_to_concatenate = 1;

bool is_filter_chunk_size = false;
std::vector<int> filter_chunk_size = {30000, 6912};

bool is_output_single_file = false;
std::string output_type = "EP_HDF5";
std::string output_file_dir = "./test-data/dir-decimate/";
std::string output_dataset = "/dat";

bool is_dir_output_match_replace_rgx = true;
std::string dir_output_match_rgx = "^(.*)\\.tdms$";
std::string dir_output_replace_rgx = "$1.h5";

void printf_help(char *cmd);

bool is_column_major = true;
bool is_column_major_from_config = false;

std::string MeasureLengthName = "MeasureLength[m]";
std::string SpatialResolutionName = "SpatialResolution[m]";
std::string SamplingFrequencyName = "SamplingFrequency[Hz]";

std::string compression_method_name = "H5Z_FILTER_DEFLATE";
std::string compression_method_parameter = "";

bool is_preprocessing = false;
std::string preprocessing_method_name = "H5Z_FILTER_SHUFFLE";
std::string preprocessing_method_parameter = "";

bool is_stencil_tag_once = false;

bool is_file_range = false;
int file_range_start_index = 0;
int file_range_end_index = 1;
std::string file_range_indexes_str;

bool is_tag_flag = true;

/**
 * @brief input data type
 *    0 : short (by default)
 *    1 : unsigned short 
 */
int input_data_type = 0;

bool isNumber(const string &str)
{
	for (char const &c : str)
	{
		if (std::isdigit(c) == 0)
			return false;
	}
	return true;
}

/**
 * @brief Split in_str into a vector string
 * 
 * @param in_str 
 * @param result_str_vector 
 */
void SplitStr2Vector(std::string in_str, char split_cha, std::vector<std::string> &result_str_vector)
{
	//std::vector<string> result;
	result_str_vector.clear();
	std::stringstream s_stream(in_str); //create string stream from the string
	while (s_stream.good())
	{
		std::string substr;
		std::getline(s_stream, substr, split_cha); //get first string delimited by comma
		result_str_vector.push_back(substr);
	}
}

// H5Z_FILTER_DEFLATE cd_values[0] = 9	Data compression filter, employing the gzip algorithm
// H5Z_FILTER_SHUFFLE	Data shuffling filter
// H5Z_FILTER_FLETCHER32  	Error detection filter, employing the Fletcher32 checksum algorithm
// H5Z_FILTER_SZIP	Data compression filter, employing the SZIP algorithm
// H5Z_FILTER_NBIT	Data compression filter, employing the N-Bit algorithm
// H5Z_FILTER_SCALEOFFSET	Data compression filter, employing the scale-offset algorithm

void FindCompressMethod(const std::string &compression_method_name, const std::string &compression_method_parameter, int &filter_id, std::vector<unsigned int> &filter_cd_values)
{
	String2Vector(compression_method_parameter, filter_cd_values);
	if (isNumber(compression_method_name))
	{
		filter_id = std::stoi(compression_method_name);
		return;
	}
	else
	{
		if (compression_method_name == "H5Z_FILTER_DEFLATE")
		{
			filter_id = H5Z_FILTER_DEFLATE;
		}
		else if (compression_method_name == "H5Z_FILTER_SHUFFLE")
		{
			filter_id = H5Z_FILTER_SHUFFLE;
		}
		else if (compression_method_name == "H5Z_FILTER_FLETCHER32")
		{
			filter_id = H5Z_FILTER_FLETCHER32;
		}
		else if (compression_method_name == "H5Z_FILTER_SZIP")
		{
			filter_id = H5Z_FILTER_SZIP;
		}
		else if (compression_method_name == "H5Z_FILTER_NBIT")
		{
			filter_id = H5Z_FILTER_NBIT;
		}
		else if (compression_method_name == "H5Z_FILTER_SCALEOFFSET")
		{
			filter_id = H5Z_FILTER_SCALEOFFSET;
		}
		else
		{
			AU_EXIT("Don't recognize filter from HDF5 library =" + compression_method_name);
		}
	}
}

//Here we just read and write data
//Compression happens during the writing
inline Stencil<std::vector<short>> udf_compress(const Stencil<short> &iStencil)
{
	std::vector<int> max_offset_upper;
	iStencil.GetOffsetUpper(max_offset_upper);
	PrintVector("max_offset_upper = ", max_offset_upper);

	//int chs_per_file_udf = max_offset_upper[0] + 1, lts_per_file_udf = max_offset_upper[1] + 1;
	int chs_per_file_udf, lts_per_file_udf;
	std::vector<int> start_offset = {0, 0};
	std::vector<int> end_offset = {max_offset_upper[0], max_offset_upper[1]};

	std::vector<short> ts_short;
	iStencil.ReadNeighbors(start_offset, end_offset, ts_short);

	Stencil<std::vector<short>> oStencil;
	//
	// Deal with tag
	//
	if (is_tag_flag && iStencil.HasTagMap() && !is_stencil_tag_once)
	{
		std::map<std::string, std::string> tag_map;
		iStencil.GetTagMap(tag_map);
		oStencil.SetTagMap(tag_map);
		if (is_output_single_file) //We only deal with meta once
			is_stencil_tag_once = true;
	}

	std::vector<size_t> vector_shape(2);
	vector_shape[0] = max_offset_upper[0] + 1;
	vector_shape[1] = max_offset_upper[1] + 1;
	PrintVector("vector_shape: ", vector_shape);
	oStencil.SetShape(vector_shape);

	// for (size_t i = 0; i < ts_short.size(); i++)
	// {
	// 	ts_short[i] = std::abs(ts_short[i]);
	// }

	// if (input_data_type == 1)
	// {
	// 	std::vector<unsigned short> ts_unsigned_short;
	// 	ts_unsigned_short.resize(ts_short.size());
	// 	for (size_t i = 0; i < ts_short.size(); i++)
	// 	{
	// 		ts_unsigned_short.push_back(ts_short[i]);
	// 	}
	// 	oStencil = ts_unsigned_short;
	// }
	// else
	// {
	oStencil = ts_short;
	//}
	return oStencil;
}

int main(int argc, char *argv[])
{
	int copt;
	bool has_config_file_flag = false;
	while ((copt = getopt(argc, argv, "c:h")) != -1)
		switch (copt)
		{
		case 'c':
			config_file.assign(optarg);
			has_config_file_flag = true;
			break;
		case 'h':
			printf_help(argv[0]);
			exit(0);
		default:
			printf("Wrong option [%c] for %s \n", copt, argv[0]);
			printf_help(argv[0]);
			exit(-1);
			break;
		}

	//Init the MPICH, etc.
	AU_Init(argc, argv);

	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int numprocs, rank, namelen;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Get_processor_name(processor_name, &namelen);
	printf("Process %d on %s out of %d\n", rank, processor_name, numprocs);

	if (has_config_file_flag)
		read_config_file(config_file, ft_rank);

	// set up the chunk size and the overlap size
	// 11648, 30000 for each dataset
	std::vector<int> chunk_size(2);
	std::vector<int> overlap_size = {0, 0};

	std::string A_endpoint_id;
	if (!is_input_single_file)
	{
		A_endpoint_id = "EP_DIR:" + input_file_type + ":" + input_dir_file;
	}
	else
	{
		A_endpoint_id = input_file_type + ":" + input_dir_file;
	}
	//+ ":" + input_h5_dataset

	if (input_file_type == "EP_HDF5")
	{
		A_endpoint_id += ":";
		A_endpoint_id += input_h5_dataset;
	}

	//Input data
	AU::Array<short> *A = new AU::Array<short>(A_endpoint_id);

	std::vector<std::string> aug_merge_index, aug_dir_sub_cmd, aug_input_search_rgx;

	//Set fhe search reges on file
	if (is_input_search_rgx && !is_input_single_file)
	{
		aug_input_search_rgx.push_back(input_search_rgx);
		A->ControlEndpoint(DIR_INPUT_SEARCH_RGX, aug_input_search_rgx);
	}

	if (is_file_range && !is_input_single_file)
	{
		std::vector<size_t> file_range_index;
		for (size_t i = file_range_start_index; i <= file_range_end_index; i++)
		{
			file_range_index.push_back(i);
		}
		file_range_indexes_str = Vector2String(file_range_index);
		std::cout << " file_range_indexes_str =" << file_range_indexes_str << "\n";
		std::vector<std::string> index_param;
		index_param.push_back(file_range_indexes_str);
		A->ControlEndpoint(DIR_FILE_SORT_INDEXES, index_param);
	}

	//Set the index to merge file
	if (!is_input_single_file)
	{
		if (is_column_major)
		{
			aug_merge_index.push_back("0");
			A->ControlEndpoint(DIR_MERGE_INDEX, aug_merge_index);
		}
		else
		{
			aug_merge_index.push_back("1");
			A->ControlEndpoint(DIR_MERGE_INDEX, aug_merge_index);
		}
	}

	if (!is_input_single_file)
	{
		std::vector<std::string> file_size_str;
		A->ControlEndpoint(DIR_GET_FILE_SIZE, file_size_str);
		String2Vector(file_size_str[0], chunk_size);
	}
	else
	{
		std::vector<unsigned long long> array_size;
		A->GetArraySize(array_size);
		chunk_size[0] = array_size[0];
		chunk_size[1] = array_size[1];
	}

	//By default, the TDMS file is column major
	if (input_file_type == "EP_TDMS")
	{
		is_column_major_from_config = true; //Set to skip the check from file
		is_column_major = true;
	}

	if (!is_column_major_from_config)
	{
		/*
         * Here we try to read  MeasureLengthName/SpatialResolutionName/SamplingFrequencyName or nTrace/nPoint
         *    to detect the layout of the file.
         * 
         */
		int meta_chs = -1, meta_time_series_points = -1;
		std::string MeasureLength, SpatialResolution, SamplingFrequency;
		A->GetTag(MeasureLengthName, MeasureLength);
		A->GetTag(SpatialResolutionName, SpatialResolution);
		if (MeasureLength.empty() || SpatialResolution.empty())
		{
			A->GetTag("nTrace", MeasureLength);
			if (MeasureLength.empty())
			{
				std::cout << "Error: Can not find out [" << MeasureLengthName << "/" << SpatialResolutionName << "] or [nTrace] in the file to detect organization. You can also set [is_column_vector = true/false] in config file to disable this detection ! \n";
				exit(-1);
			}
			meta_chs = std::stoi(MeasureLength);
		}
		else
		{
			meta_chs = std::stoi(MeasureLength) / std::stoi(SpatialResolution);
		}

		A->GetTag(SamplingFrequencyName, SamplingFrequency);
		if (SamplingFrequency.empty())
		{
			A->GetTag("nPoint", SamplingFrequency);
			if (SamplingFrequency.empty())
			{
				std::cout << "Error: Can not find out [" << SamplingFrequencyName << "] or [nPoint] in the file to detect organization. You can also set [is_column_vector = true/false] in config file to disable this detection ! \n";
				exit(-1);
			}
			meta_time_series_points = std::stoi(SamplingFrequency);
		}
		else
		{
			meta_time_series_points = 60 * std::stoi(SamplingFrequency);
		}

		//std::cout << "meta_time_series_points = " << meta_time_series_points << " , meta_chs =  " << meta_chs << " \n";
		if (chunk_size[0] == meta_time_series_points && chunk_size[1] == meta_chs)
		{
			is_column_major = true;
			if (!ft_rank)
			{
				std::cout << termcolor::reset << "\n";
				std::cout << termcolor::magenta << "Found data organization = " << termcolor::green << " column vector (time by channel) \n";
				std::cout << termcolor::reset << "\n";
			}
		}
		else if (chunk_size[0] == meta_chs && chunk_size[1] == meta_time_series_points)
		{
			is_column_major = false;
			if (!ft_rank)
			{
				std::cout << termcolor::reset << "\n";
				std::cout << termcolor::magenta << "Found data organization = " << termcolor::green << " row vector (channel by time)\n";
				std::cout << termcolor::reset << "\n";
			}
		}
		else
		{
			std::cout << "Metadata and data are inconsistent in the size ! " << std::endl;
			exit(-1);
		}
	}

	if (is_column_major)
	{
		chunk_size[0] = chunk_size[0] * n_files_to_concatenate;
		lts_per_file = chunk_size[0];
		chs_per_file = chunk_size[1];
	}
	else
	{
		chunk_size[1] = chunk_size[1] * n_files_to_concatenate;
		lts_per_file = chunk_size[1];
		chs_per_file = chunk_size[0];
	}

	PrintVector("chunk_size = ", chunk_size);
	if (!ft_rank)
		std::cout << "lts_per_file = " << lts_per_file << ",chs_per_file = " << chs_per_file << "\n";

	///Users/dbin/work/arrayudf-git-svn-test-on-bitbucket/examples/das/tdms-dir", chunk_size, overlap_size);
	A->SetChunkSize(chunk_size);
	A->SetOverlapSize(overlap_size);

	//Extract tag to be along with Input Stencil
	A->GetStencilTag();

	if (!is_column_major && input_file_type == "EP_TDMS")
	{
		aug_dir_sub_cmd.push_back("BINARY_ENABLE_TRANSPOSE_ON_READ");
		A->ControlEndpoint(DIR_SUB_CMD_ARG, aug_dir_sub_cmd);
	}

	//Result data
	//AU::Array<short> *B;

	AU::Array<short> *B;
	//AuEndpointDataType t;
	if (is_output_single_file)
	{
		//Store into a single file
		//B = new AU::Array<short>(output_type + ":" + output_file_dir + ":" + output_dataset);
		// if (input_data_type == 0)
		// {
		// 	t = AuEndpointDataType::AU_SHORT;
		B = new AU::Array<short>(output_type + ":" + output_file_dir + ":" + output_dataset);
		// }
		// else if (input_data_type == 1)
		// {
		// 	t = AuEndpointDataType::AU_USHORT;
		// 	B = new AU::Array<unsigned short>(output_type + ":" + output_file_dir + ":" + output_dataset);
		// }
	}
	else
	{
		//Store into multiple file
		//B = new AU::Array<short>("EP_DIR:" + output_type + ":" + output_file_dir + ":" + output_dataset);

		// if (input_data_type == 0)
		// {
		// 	t = AuEndpointDataType::AU_SHORT;
		B = new AU::Array<short>("EP_DIR:" + output_type + ":" + output_file_dir + ":" + output_dataset);
		// }
		// else if (input_data_type == 1)
		// {
		// 	t = AuEndpointDataType::AU_USHORT;
		// 	B = new AU::Array<unsigned short>("EP_DIR:" + output_type + ":" + output_file_dir + ":" + output_dataset);
		// }

		//Use the below rgx pattern to name the file
		std::vector<std::string> aug_output_replace_arg;
		aug_output_replace_arg.push_back(dir_output_match_rgx);
		aug_output_replace_arg.push_back(dir_output_replace_rgx);
		B->ControlEndpoint(DIR_MERGE_INDEX, aug_merge_index);

		if (is_dir_output_match_replace_rgx)
			B->ControlEndpoint(DIR_OUPUT_REPLACE_RGX, aug_output_replace_arg);
	}

	if (is_preprocessing)
	{
		std::vector<std::string> filter_preprocessing_parameter;
		int preprocessing_id;
		std::vector<unsigned int> preprocessing_cd_values;

		FindCompressMethod(preprocessing_method_name, preprocessing_method_parameter, preprocessing_id, preprocessing_cd_values);

		filter_preprocessing_parameter.push_back(std::to_string(preprocessing_id));
		filter_preprocessing_parameter.push_back(Vector2String(preprocessing_cd_values));

		if (is_output_single_file)
		{
			B->ControlEndpoint(HDF5_ENABLE_FILTER_PREPROCESSING, filter_preprocessing_parameter);
		}
		else
		{
			filter_preprocessing_parameter.insert(filter_preprocessing_parameter.begin(), std::to_string(HDF5_ENABLE_FILTER_PREPROCESSING));
			B->ControlEndpoint(DIR_SUB_CMD_ARG, filter_preprocessing_parameter);
		}
	}

	//Enable the filter ID
	std::vector<std::string> filter_paramter;
	std::vector<std::string> control_paramter_null;
	// H5Z_FILTER_DEFLATE cd_values[0] = 9	Data compression filter, employing the gzip algorithm
	// H5Z_FILTER_SHUFFLE	Data shuffling filter
	// H5Z_FILTER_FLETCHER32  	Error detection filter, employing the Fletcher32 checksum algorithm
	// H5Z_FILTER_SZIP	Data compression filter, employing the SZIP algorithm
	// H5Z_FILTER_NBIT	Data compression filter, employing the N-Bit algorithm
	// H5Z_FILTER_SCALEOFFSET	Data compression filter, employing the scale-offset algorithm
	//#define H5Z_FILTER_SZ 32017

	std::vector<std::string> compression_method_name_vector;
	std::vector<std::string> compression_method_parameter_vector;
	SplitStr2Vector(compression_method_name, ',', compression_method_name_vector);
	SplitStr2Vector(compression_method_parameter, ':', compression_method_parameter_vector);

	assert(compression_method_name_vector.size() == compression_method_parameter_vector.size());

	for (int i = 0; i < compression_method_name_vector.size(); i++)
	{
		compression_method_name = compression_method_name_vector[i];
		compression_method_parameter = compression_method_parameter_vector[i];
		std::cout << i << " id = : " << compression_method_name << ", parameter = " << compression_method_parameter << "\n";
		int filter_id = H5Z_FILTER_DEFLATE;
		std::vector<unsigned int> filter_cd_values;

		if (!is_filter_chunk_size)
		{
			filter_chunk_size = chunk_size;
		}
		filter_paramter.clear();
		FindCompressMethod(compression_method_name, compression_method_parameter, filter_id, filter_cd_values);
		filter_paramter.push_back(std::to_string(filter_id));
		filter_paramter.push_back(Vector2String(filter_cd_values));
		filter_paramter.push_back(Vector2String(filter_chunk_size));
		//  parameter_v[0]: filter_id
		//  parameter_v[2]: filter_cd_values
		//  parameter_v[3]: chunk_size
		if (is_output_single_file)
		{
			B->ControlEndpoint(HDF5_ENABLE_FILTER, filter_paramter);
			B->ControlEndpoint(HDF5_DISABLE_MPI_IO, control_paramter_null);
			B->ControlEndpoint(HDF5_DISABLE_COLLECTIVE_IO, control_paramter_null);
		}
		else
		{
			filter_paramter.insert(filter_paramter.begin(), std::to_string(HDF5_ENABLE_FILTER));
			B->ControlEndpoint(DIR_SUB_CMD_ARG, filter_paramter);

			control_paramter_null.clear();
			control_paramter_null.push_back(to_string(HDF5_DISABLE_MPI_IO));
			B->ControlEndpoint(DIR_SUB_CMD_ARG, control_paramter_null);

			control_paramter_null.clear();
			control_paramter_null.push_back(std::to_string(HDF5_DISABLE_COLLECTIVE_IO));
			B->ControlEndpoint(DIR_SUB_CMD_ARG, control_paramter_null);
		}
	}
	//
	//Stride on execution
	//Each chunk only runs the udf_compress once
	A->EnableApplyStride(chunk_size);

	A->SetVectorDirection(AU_FLAT_OUTPUT_ROW);

	//Run
	A->Transform(udf_compress, B);
	// if (t == AuEndpointDataType::AU_SHORT)
	// {
	// 	//AU::Array<short> *BB = static_cast<FT::Array<short> *>(B);
	// 	A->Transform2<std::vector<short>>(udf_compress, B);
	// }
	// else
	// {
	// 	//AU::Array<unsigned short> *BB = static_cast<FT::Array<unsigned short> *>(B);
	// 	A->Transform2<std::vector<unsigned short>>(udf_compress, B);
	// }

	A->ReportTime();
	//Clear
	delete A;
	delete B;

	AU_Finalize();

	return 0;
}

void printf_help(char *cmd)
{
	char *msg = (char *)"Usage: %s [OPTION]\n\
      	  -h help (--help)\n\
         -c config file for parameters  \n\
          Example: mpirun -n 1 %s -c  compress.config\n";

	fprintf(stdout, msg, cmd, cmd);
}

int read_config_file(std::string file_name, int mpi_rank)
{
	INIReader reader(file_name);

	if (reader.ParseError() < 0)
	{
		std::cout << "Can't load [" << file_name << " ]\n";
		return 1;
	}

	input_dir_file = reader.Get("parameter", "input_dir_file", "/Users/dbin/work/arrayudf-git-svn-test-on-bitbucket/examples/das/tdms-dir");

	input_file_type = reader.Get("parameter", "input_file_type", "EP_TDMS");

	std::string temp_str = reader.Get("parameter", "is_input_search_rgx", "false");

	is_input_search_rgx = (temp_str == "false") ? false : true;

	if (is_input_search_rgx)
	{
		input_search_rgx = reader.Get("parameter", "input_search_rgx", "^(.*)[1234]\\.tdms$");
	}

	std::string is_file_range_str = reader.Get("parameter", "is_file_range", "false");
	if (is_file_range_str == "false" || is_file_range_str == "0")
	{
		is_file_range = false;
	}
	else if (is_file_range_str == "true" || is_file_range_str == "1")
	{
		is_file_range = true;
	}
	else
	{
		AU_EXIT("Don't read the is_file_range's value " + is_file_range_str);
	}

	if (is_file_range)
	{
		file_range_start_index = reader.GetInteger("parameter", "file_range_start_index", 0);
		file_range_end_index = reader.GetInteger("parameter", "file_range_end_index", 1);
	}

	chs_per_file = reader.GetInteger("parameter", "chs_per_file", 11648);
	lts_per_file = reader.GetInteger("parameter", "lts_per_file", 30000);

	//bool is_filter_chunk_size = false;
	//std::vector<int> filter_chunk_size = {30000, 6912};
	temp_str = reader.Get("parameter", "is_compression_chunk_size", "false");
	is_filter_chunk_size = (temp_str == "false" || temp_str == "0") ? false : true;
	if (is_filter_chunk_size)
	{
		temp_str = reader.Get("parameter", "compression_chunk_size", "30000, 6912");
		String2Vector(temp_str, filter_chunk_size);
	}

	temp_str = reader.Get("parameter", "is_column_vector", "NULL");
	if (temp_str == "NULL")
	{
		is_column_major_from_config = false;
	}
	else
	{
		is_column_major = (temp_str == "false" || temp_str == "0") ? false : true;
		is_column_major_from_config = true;
	}

	input_h5_dataset = reader.Get("parameter", "input_dataset", "/dat");

	temp_str = reader.Get("parameter", "is_input_single_file", "false");
	is_input_single_file = (temp_str == "false") ? false : true;

	n_files_to_concatenate = reader.GetInteger("parameter", "n_files_to_concatenate", 2);

	MeasureLengthName = reader.Get("parameter", "attribute_name_measure_length", "MeasureLength[m]");
	SpatialResolutionName = reader.Get("parameter", "attribute_name_spatial_resolution", "SpatialResolution[m]");
	SamplingFrequencyName = reader.Get("parameter", "attribute_name_sampling_frequency", "SamplingFrequency[Hz]");

	compression_method_name = reader.Get("parameter", "compression_method_name", "H5Z_FILTER_DEFLATE");
	compression_method_parameter = reader.Get("parameter", "compression_method_parameter", "");

	temp_str = reader.Get("parameter", "is_preprocessing", "false");
	is_preprocessing = (temp_str == "false" || temp_str == "0") ? false : true;
	if (is_preprocessing)
	{
		preprocessing_method_name = reader.Get("parameter", "preprocessing_method_name", "H5Z_FILTER_SHUFFLE");
		preprocessing_method_parameter = reader.Get("parameter", "preprocessing_method_parameter", "");
	}
	temp_str = reader.Get("parameter", "is_output_single_file", "false");

	is_output_single_file = (temp_str == "false") ? false : true;

	output_type = reader.Get("parameter", "output_type", "EP_HDF5");

	output_file_dir = reader.Get("parameter", "output_file_dir", "./tdms-dir-dec/test.h5");

	output_dataset = reader.Get("parameter", "output_dataset", "/data");

	temp_str = reader.Get("parameter", "is_dir_output_match_replace_rgx", "false");

	is_dir_output_match_replace_rgx = (temp_str == "false") ? false : true;

	if (is_dir_output_match_replace_rgx)
	{
		dir_output_match_rgx = reader.Get("parameter", "output_file_regex_match", "^(.*)\\.tdms$");

		dir_output_replace_rgx = reader.Get("parameter", "output_file_regex_replace", "$1.h5");
	}

	if (!mpi_rank)
	{
		std::cout << "\n\n";
		std::cout << termcolor::red << "Parameters to run the Decimate: ";

		std::cout << termcolor::blue << "\n\n Input parameters: ";

		std::cout << termcolor::magenta << "\n        is_input_single_file = " << termcolor::green << is_input_single_file;

		std::cout << termcolor::magenta << "\n        input_dir_file = " << termcolor::green << input_dir_file;
		std::cout << termcolor::magenta << "\n        input_file_type = " << termcolor::green << input_file_type;
		if (is_filter_chunk_size)
			std::cout << termcolor::magenta << "\n        compression_chunk_size = " << termcolor::green << filter_chunk_size[0] << ", " << filter_chunk_size[1];

		std::cout << termcolor::magenta << "\n        compression_method_name = " << termcolor::green << compression_method_name << " \n";

		std::cout << termcolor::magenta << "\n        compression_method_parameter = " << termcolor::green << compression_method_parameter << " \n";

		if (is_preprocessing)
		{
			std::cout << termcolor::magenta << "\n        preprocessing_method_name = " << termcolor::green << preprocessing_method_name << " \n";

			std::cout << termcolor::magenta << "\n        preprocessing_method_parameter = " << termcolor::green << preprocessing_method_parameter << " \n";
		}

		if (is_input_search_rgx)
		{
			std::cout << termcolor::magenta << "\n        input_search_rgx = " << termcolor::green << input_search_rgx;
		}
		std::cout << termcolor::blue << "\n\n Runtime parameters: ";
		std::cout << termcolor::magenta << "\n\n        lts_per_file = " << termcolor::green << lts_per_file;
		std::cout << termcolor::magenta << "\n        chs_per_file = " << termcolor::green << chs_per_file;
		std::cout << termcolor::magenta << "\n        n_files_to_concatenate = " << termcolor::green << n_files_to_concatenate;

		std::cout << termcolor::blue << "\n\n Output parameters: ";

		std::cout << termcolor::magenta << "\n        is_output_single_file = " << termcolor::green << is_output_single_file;
		std::cout << termcolor::magenta << "\n        output_type = " << termcolor::green << output_type;
		std::cout << termcolor::magenta << "\n        output_file_dir = " << termcolor::green << output_file_dir;
		std::cout << termcolor::magenta << "\n        output_dataset = " << termcolor::green << output_dataset;

		if (is_dir_output_match_replace_rgx)
		{
			std::cout << termcolor::magenta << "\n        dir_output_match_rgx = " << termcolor::green << dir_output_match_rgx;
			std::cout << termcolor::magenta << "\n        dir_output_replace_rgx = " << termcolor::green << dir_output_replace_rgx;
		}

		std::cout << termcolor::reset << "\n\n";
	}
	fflush(stdout);

	return 0;
}
