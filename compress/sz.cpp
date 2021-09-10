

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>
#include "hdf5.h"
#include "H5Z_SZ.h"

#define DATASET "testdata_compressed"

//argc is like the dims of the data
/**
 * @brief 
 * 
 * @param dims the number of dimensions of the data
 * @param argv 
 *            argv[0]: data type
 * @param filename 
 * @param dataset 
 * @return int 
 */
int sz_compress(int dims, std::vector<std::string> argv, std::string filename, std::string dataset)
{
    size_t r5 = 0, r4 = 0, r3 = 0, r2 = 0, r1 = 0;
    char outDir[640], oriFilePath[640], outputFilePath[640];
    size_t cd_nelmts, nbEle;
    unsigned int *cd_values = NULL;
    //unsigned int cd_values[7];

    herr_t status;
    htri_t avail;
    unsigned filter_config;

    hid_t sid, idsid, cpid, fid;

    // if (argc < 4)
    // {
    //     printf("Test case: szToHDF5 [dataType] [config_file] [srcFilePath] [dimension sizes...]\n");
    //     printf("Example1 : szToHDF5 -f sz.config testdata/x86/testfloat_8_8_128.dat 8 8 128\n");
    //     printf("Example 2: szToHDF5 -i32 sz.config testdata/x86/testint32_8x8x8.dat 8 8 8\n");
    //     exit(0);
    // }

    printf("config file = %s\n", argv[2]);

    int dataType = 0;
    if (argv[0] == "-f")
        dataType = SZ_FLOAT;
    else if (argv[0] == "-d")
        dataType = SZ_DOUBLE;
    else if (argv[0] == "-i8")
        dataType = SZ_INT8;
    else if (argv[0] == "-i16")
        dataType = SZ_INT16;
    else if (argv[0] == "-i32")
        dataType = SZ_INT32;
    else if (argv[0] == "-i64")
        dataType = SZ_INT64;
    else if (argv[0] == "-u8")
        dataType = SZ_UINT8;
    else if (argv[0] == "-u16")
        dataType = SZ_UINT16;
    else if (argv[0] == "-u32")
        dataType = SZ_UINT32;
    else if (argv[0] == "-u64")
        dataType = SZ_UINT64;
    else
    {
        printf("Error: unknown data type in szToHDF5.c!\n");
        exit(0);
    }

    strcpy(cfgFile, argv[2]);
    sprintf(oriFilePath, "%s", argv[3]);
    if (argc >= 5)
        r1 = atoi(argv[4]); //8
    if (argc >= 6)
        r2 = atoi(argv[5]); //8
    if (argc >= 7)
        r3 = atoi(argv[6]); //128
    if (argc >= 8)
        r4 = atoi(argv[7]);
    if (argc >= 9)
        r5 = atoi(argv[8]);

    printf("cfgFile=%s\n", cfgFile);
    sprintf(outputFilePath, "%s.sz.h5", oriFilePath);
    printf("output file: %s \n", outputFilePath);

    //	printf("argv[1]=%s, dataType=%d\n", argv[1], dataType);
    nbEle = computeDataLength(r5, r4, r3, r2, r1);

    //	printf("nbEle=%u\n", nbEle);

    //Create cd_values
    printf("Dimension sizes: n5=%u, n4=%u, n3=%u, n2=%u, n1=%u\n", r5, r4, r3, r2, r1);
    SZ_metaDataToCdArray(&cd_nelmts, &cd_values, dataType, r5, r4, r3, r2, r1);
    //load_conffile_flag = 0;
    //SZ_metaDataErrToCdArray(&cd_nelmts, &cd_values, dataType, r5, r4, r3, r2, r1, REL, 0.01, 0.01, 0, 0); //SZ_FLOAT or SZ_DOUBLE or SZ_INT 100x500x500 : 0, 0, 100, 500, 500, ABS, REL (0.01, 0.01*(max-min), PW_REL (0.01, 5, 6, 7, 8, 9 --> 5*0.01, 6*0.01, ...), PSNR (mean squared error)).
    /*cd_nelmts = 5;
	cd_values[0] = 3;
	cd_values[1] = 0;
	cd_values[2] = 128;
	cd_values[3] = 8;
	cd_values[4] = 8;
	cd_values[5] = 0;				
	cd_values[6] = 0;*/

    int i = 0;
    //	for(i=0;i<cd_nelmts;i++)
    //		printf("cd_values[%d]=%u\n", i, cd_values[i]);

    //compute dimension
    int dim = computeDimension(r5, r4, r3, r2, r1);

    hsize_t dims[5] = {0, 0, 0, 0, 0}, chunk[5] = {0, 0, 0, 0, 0};
    init_dims_chunk(dim, dims, chunk, nbEle, r5, r4, r3, r2, r1);

    /* create HDF5 file */
    if (0 > (fid = H5Fcreate(outputFilePath, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)))
        ERROR(H5Fcreate);

    /*Create dataspace. Setting maximum size */
    if (0 > (sid = H5Screate_simple(dim, dims, NULL)))
        ERROR(H5Screate_simple);

    /* setup dataset creation properties */
    if (0 > (cpid = H5Pcreate(H5P_DATASET_CREATE)))
        ERROR(H5Pcreate);

    /* Add the SZ compression filter and set the chunk size */
    if (0 > H5Pset_filter(cpid, H5Z_FILTER_SZ, H5Z_FLAG_MANDATORY, cd_nelmts, cd_values))
        ERROR(H5Pset_filter);
    avail = H5Zfilter_avail(H5Z_FILTER_SZ);
    if (avail)
    {
        status = H5Zget_filter_info(H5Z_FILTER_SZ, &filter_config);
        if (filter_config & H5Z_FILTER_CONFIG_ENCODE_ENABLED)
            printf("sz filter is available for encoding and decoding.\n");
    }
    if (0 > H5Pset_chunk(cpid, dim, chunk))
        ERROR(H5Pset_chunk);

    //Initialize the configuration for SZ
    //You can also use the global variable conf_params to set the configuratoin for sz without cfgFile.
    //Example of setting an absolute error bound:
    //			sz_params* params = H5Z_SZ_Init_Default();
    //			params->errorBoundMode = ABS;
    //			params->absErrBound = 1E-4;

    //H5Z_SZ_Init(cfgFile);

    printf("....Writing SZ compressed data.............\n");

    if (dataType == SZ_FLOAT)
    {
        float *data = readFloatData(oriFilePath, &nbEle, &status);

        printf("original data = ");
        for (i = 0; i < 20; i++)
            printf("%f ", data[i]);
        printf("....\n");

        if (dataEndianType == LITTLE_ENDIAN_DATA)
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_IEEE_F32LE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT)))
                ERROR(H5Dcreate);
            if (0 > H5Dwrite(idsid, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data))
                ERROR(H5Dwrite);
        }
        else //BIG_ENDIAN_DATA
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_IEEE_F32BE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT)))
                ERROR(H5Dcreate);
            if (0 > H5Dwrite(idsid, H5T_IEEE_F32BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data))
                ERROR(H5Dwrite);
        }
        free(data);
        if (0 > H5Dclose(idsid))
            ERROR(H5Dclose);
    }
    else if (dataType == SZ_DOUBLE)
    {
        double *data = readDoubleData(oriFilePath, &nbEle, &status);

        printf("original data = ");
        for (i = 0; i < 20; i++)
            printf("%f ", data[i]);
        printf("....\n");

        if (dataEndianType == LITTLE_ENDIAN_DATA)
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_IEEE_F64LE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT)))
                ERROR(H5Dcreate);
            if (0 > H5Dwrite(idsid, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data))
                ERROR(H5Dwrite);
        }
        else //BIG_ENDIAN_DATA
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_IEEE_F64BE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT)))
                ERROR(H5Dcreate);
            if (0 > H5Dwrite(idsid, H5T_IEEE_F64BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data))
                ERROR(H5Dwrite);
        }
        free(data);
        if (0 > H5Dclose(idsid))
            ERROR(H5Dclose);
    }
    else if (dataType == SZ_INT8)
    {
        char *data = (char *)readInt8Data(oriFilePath, &nbEle, &status);

        printf("original data = ");
        for (i = 0; i < 20; i++)
            printf("%d ", data[i]);
        printf("....\n");

        if (dataEndianType == LITTLE_ENDIAN_DATA)
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_I8LE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT)))
                ERROR(H5Dcreate);
            if (0 > H5Dwrite(idsid, H5T_STD_I8LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data))
                ERROR(H5Dwrite);
        }
        else //BIG_ENDIAN_DATA
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_I8BE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT)))
                ERROR(H5Dcreate);
            if (0 > H5Dwrite(idsid, H5T_STD_I8BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data))
                ERROR(H5Dwrite);
        }
        free(data);
        if (0 > H5Dclose(idsid))
            ERROR(H5Dclose);
    }
    else if (dataType == SZ_UINT8)
    {
        unsigned char *data = readByteData(oriFilePath, &nbEle, &status);

        printf("original data = ");
        for (i = 0; i < 20; i++)
            printf("%d ", data[i]);
        printf("....\n");

        if (dataEndianType == LITTLE_ENDIAN_DATA)
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_U8LE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT)))
                ERROR(H5Dcreate);
            if (0 > H5Dwrite(idsid, H5T_STD_U8LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data))
                ERROR(H5Dwrite);
        }
        else //BIG_ENDIAN_DATA
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_U8BE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT)))
                ERROR(H5Dcreate);
            if (0 > H5Dwrite(idsid, H5T_STD_U8BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data))
                ERROR(H5Dwrite);
        }
        free(data);
        if (0 > H5Dclose(idsid))
            ERROR(H5Dclose);
    }
    else if (dataType == SZ_INT16)
    {
        short *data = readInt16Data(oriFilePath, &nbEle, &status);

        printf("original data = ");
        for (i = 0; i < 20; i++)
            printf("%d ", data[i]);
        printf("....\n");

        if (dataEndianType == LITTLE_ENDIAN_DATA)
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_I16LE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT)))
                ERROR(H5Dcreate);
            if (0 > H5Dwrite(idsid, H5T_STD_I16LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data))
                ERROR(H5Dwrite);
        }
        else //BIG_ENDIAN_DATA
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_I16BE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT)))
                ERROR(H5Dcreate);
            if (0 > H5Dwrite(idsid, H5T_STD_I16BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data))
                ERROR(H5Dwrite);
        }
        free(data);
        if (0 > H5Dclose(idsid))
            ERROR(H5Dclose);
    }
    else if (dataType == SZ_UINT16)
    {
        unsigned short *data = readUInt16Data(oriFilePath, &nbEle, &status);

        printf("original data = ");
        for (i = 0; i < 20; i++)
            printf("%d ", data[i]);
        printf("....\n");

        if (dataEndianType == LITTLE_ENDIAN_DATA)
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_U16LE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT)))
                ERROR(H5Dcreate);
            if (0 > H5Dwrite(idsid, H5T_STD_U16LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data))
                ERROR(H5Dwrite);
        }
        else //BIG_ENDIAN_DATA
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_U16BE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT)))
                ERROR(H5Dcreate);
            if (0 > H5Dwrite(idsid, H5T_STD_U16BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data))
                ERROR(H5Dwrite);
        }
        free(data);
        if (0 > H5Dclose(idsid))
            ERROR(H5Dclose);
    }
    else if (dataType == SZ_INT32)
    {
        int *data = readInt32Data(oriFilePath, &nbEle, &status);

        printf("original data = ");
        for (i = 0; i < 20; i++)
            printf("%d ", data[i]);
        printf("....\n");

        if (dataEndianType == LITTLE_ENDIAN_DATA)
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_I32LE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT)))
                ERROR(H5Dcreate);
            if (0 > H5Dwrite(idsid, H5T_STD_I32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data))
                ERROR(H5Dwrite);
        }
        else //BIG_ENDIAN_DATA
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_I32BE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT)))
                ERROR(H5Dcreate);
            if (0 > H5Dwrite(idsid, H5T_STD_I32BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data))
                ERROR(H5Dwrite);
        }
        free(data);
        if (0 > H5Dclose(idsid))
            ERROR(H5Dclose);
    }
    else if (dataType == SZ_UINT32)
    {
        unsigned int *data = readUInt32Data(oriFilePath, &nbEle, &status);

        printf("original data = ");
        for (i = 0; i < 20; i++)
            printf("%d ", data[i]);
        printf("....\n");

        if (dataEndianType == LITTLE_ENDIAN_DATA)
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_U32LE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT)))
                ERROR(H5Dcreate);
            if (0 > H5Dwrite(idsid, H5T_STD_U32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data))
                ERROR(H5Dwrite);
        }
        else //BIG_ENDIAN_DATA
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_U32BE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT)))
                ERROR(H5Dcreate);
            if (0 > H5Dwrite(idsid, H5T_STD_U32BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data))
                ERROR(H5Dwrite);
        }
        free(data);
        if (0 > H5Dclose(idsid))
            ERROR(H5Dclose);
    }
    else if (dataType == SZ_INT64)
    {
        long *data = (long *)readInt64Data(oriFilePath, &nbEle, &status);

        printf("original data = ");
        for (i = 0; i < 20; i++)
            printf("%ld ", data[i]);
        printf("....\n");

        if (dataEndianType == LITTLE_ENDIAN_DATA)
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_I64LE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT)))
                ERROR(H5Dcreate);
            if (0 > H5Dwrite(idsid, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data))
                ERROR(H5Dwrite);
        }
        else //BIG_ENDIAN_DATA
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_I64BE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT)))
                ERROR(H5Dcreate);
            if (0 > H5Dwrite(idsid, H5T_STD_I64BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data))
                ERROR(H5Dwrite);
        }
        free(data);
        if (0 > H5Dclose(idsid))
            ERROR(H5Dclose);
    }
    else if (dataType == SZ_UINT64)
    {
        unsigned long *data = (unsigned long *)readUInt64Data(oriFilePath, &nbEle, &status);

        printf("original data = ");
        for (i = 0; i < 20; i++)
            printf("%ld ", data[i]);
        printf("....\n");

        if (dataEndianType == LITTLE_ENDIAN_DATA)
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_U64LE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT)))
                ERROR(H5Dcreate);
            if (0 > H5Dwrite(idsid, H5T_STD_U64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data))
                ERROR(H5Dwrite);
        }
        else //BIG_ENDIAN_DATA
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_U64BE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT)))
                ERROR(H5Dcreate);
            if (0 > H5Dwrite(idsid, H5T_STD_U64BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data))
                ERROR(H5Dwrite);
        }
        free(data);
        if (0 > H5Dclose(idsid))
            ERROR(H5Dclose);
    }
    else
    {
        printf("Error: unknown data type in szToHDF5.c!\n");
        exit(0);
    }

    /*Close and release reosurces*/
    if (0 > H5Sclose(sid))
        ERROR(H5Sclose);
    if (0 > H5Pclose(cpid))
        ERROR(H5Pclose);
    if (0 > H5Fclose(fid))
        ERROR(H5Fclose);
    free(cd_values);
    printf("Output hdf5 file: %s\n", outputFilePath);
    H5Z_SZ_Finalize();
    H5close();
    return 0;
}

int sz_decompress(int argc, char *argv[])
{
    int dimSize = 0;
    size_t r5 = 0, r4 = 0, r3 = 0, r2 = 0, r1 = 0, nbEle = 0;
    char hdf5FilePath[640], outputFilePath[640];
    hid_t file, dset, dcpl, space_id, dtype; /*Handles*/
    H5Z_filter_t filter_id = 0;
    herr_t status;
    H5T_class_t type_class;
    H5T_sign_t dsign;
    H5T_order_t dorder;

    htri_t avail;
    char filter_name[80];
    unsigned int flags = 0;
    size_t nelmts = 0, dsize;
    unsigned int values_out[7] = {0, 0, 0, 0, 0, 0, 0}; //at most 7 parameters

    if (argc < 2)
    {
        printf("Test case: dszFromHDF5 [hdf5FilePath]\n");
        printf("Example 1: dszFromHDF5 testdata/x86/testfloat_8_8_128.dat.sz.hdf5\n");
        printf("Example 2: dszFromHDF5 testdata/x86/testint32_8x8x8.dat.sz.hdf5\n");
        exit(0);
    }

    sprintf(hdf5FilePath, "%s", argv[1]);
    sprintf(outputFilePath, "%s.out.h5", hdf5FilePath);

    /*Open the hdf5 file with SZ-compressed data*/
    file = H5Fopen(hdf5FilePath, H5F_ACC_RDONLY, H5P_DEFAULT);
    dset = H5Dopen(file, DATASET, H5P_DEFAULT);

    /*Retrieve dataset creation property list.*/
    dcpl = H5Dget_create_plist(dset);

    /*Check that filter is not registered with the library yet*/
    avail = H5Zfilter_avail(H5Z_FILTER_SZ);
    if (!avail)
        printf("sz filter is not yet available after the H5Pget_filter call.\n");
    else
        printf("sz filter is available.\n");

    space_id = H5Dget_space(dset);
    nbEle = H5Sget_simple_extent_npoints(space_id);

    if ((dtype = H5Dget_type(dset)) < 0)
        printf("Error: H5Dget_type(dset) < 0\n");

    /*Read the data using the default properties.*/
    printf("....Reading SZ compressed data .....................\n");

    if ((type_class = H5Tget_class(dtype)) < 0)
    {
        printf("Error: H5Tget_class<0\n");
        exit(0);
    }
    if (0 == (dsize = H5Tget_size(dtype)))
    {
        printf("Error: H5Tget_size==0\n");
        exit(0);
    }

    if ((dorder = H5Tget_order(dtype)) < 0)
        printf("Error: H5Tget_order<0\n");

    switch (type_class)
    {
    case H5T_FLOAT:
        if (H5Tequal(dtype, H5T_IEEE_F32BE) == 1 || H5Tequal(dtype, H5T_IEEE_F32LE) == 1 || H5Tequal(dtype, H5T_NATIVE_FLOAT) == 1)
        {
            printf("data type: float\n");
            float *data = (float *)malloc(sizeof(float) * nbEle);
            if (dorder == H5T_ORDER_LE)
                status = H5Dread(dset, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
            else //H5T_ORDER_BE
                status = H5Dread(dset, H5T_IEEE_F32BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
            /*Print the first 20 data values to check the correctness.*/
            int i;
            printf("reconstructed data = ");
            for (i = 0; i < 20; i++)
                printf("%f ", data[i]);
            printf("\n");
            free(data);
        }
        else //64bit: double
        {
            printf("data type: double\n");
            double *data = (double *)malloc(sizeof(double) * nbEle);
            if (dorder == H5T_ORDER_LE)
                status = H5Dread(dset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
            else
                status = H5Dread(dset, H5T_IEEE_F64BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
            /*Print the first 10 data values to check the correctness.*/
            int i;
            printf("reconstructed data = ");
            for (i = 0; i < 20; i++)
                printf("%f ", data[i]);
            printf("\n");
            free(data);
        }
        break;
    case H5T_INTEGER:
        if (0 > (dsign = H5Tget_sign(dtype)))
        {
            printf("Error in calling H5Tget_sign(type_id)....\n");
            exit(0);
        }
        if (dsign == H5T_SGN_NONE) //unsigned
        {
            if (dsize == 1)
            {
                printf("data type: unsigned char\n");
                unsigned char *data = (unsigned char *)malloc(sizeof(unsigned char) * nbEle);
                if (dorder == H5T_ORDER_LE)
                    status = H5Dread(dset, H5T_STD_U8LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                else
                    status = H5Dread(dset, H5T_STD_U8BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                int i;
                printf("reconstructed data = ");
                for (i = 0; i < 20; i++)
                    printf("%d ", data[i]);
                printf("\n");
                free(data);
            }
            else if (dsize == 2)
            {
                printf("data type: unsigned short\n");
                unsigned short *data = (unsigned short *)malloc(sizeof(unsigned short) * nbEle);
                if (dorder == H5T_ORDER_LE)
                    status = H5Dread(dset, H5T_STD_U16LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                else
                    status = H5Dread(dset, H5T_STD_U16BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                int i;
                printf("reconstructed data = ");
                for (i = 0; i < 20; i++)
                    printf("%d ", data[i]);
                printf("\n");
                free(data);
            }
            else if (dsize == 4)
            {
                printf("data type: unsigned int\n");
                unsigned int *data = (unsigned int *)malloc(sizeof(unsigned int) * nbEle);
                if (dorder == H5T_ORDER_LE)
                    status = H5Dread(dset, H5T_STD_U32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                else
                    status = H5Dread(dset, H5T_STD_U32BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                int i;
                printf("reconstructed data = ");
                for (i = 0; i < 20; i++)
                    printf("%d ", data[i]);
                printf("\n");
                free(data);
            }
            else if (dsize == 8)
            {
                printf("data type: unsigned long\n");
                unsigned long *data = (unsigned long *)malloc(sizeof(unsigned long) * nbEle);
                if (dorder == H5T_ORDER_LE)
                    status = H5Dread(dset, H5T_STD_U64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                else
                    status = H5Dread(dset, H5T_STD_U64BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                int i;
                printf("reconstructed data = ");
                for (i = 0; i < 20; i++)
                    printf("%ld ", data[i]);
                printf("\n");
                free(data);
            }
        }
        else
        {
            if (dsize == 1)
            {
                printf("data type: char\n");
                char *data = (char *)malloc(sizeof(char) * nbEle);
                if (dorder == H5T_ORDER_LE)
                    status = H5Dread(dset, H5T_STD_I8LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                else
                    status = H5Dread(dset, H5T_STD_I8BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                int i;
                printf("reconstructed data = ");
                for (i = 0; i < 20; i++)
                    printf("%d ", data[i]);
                printf("\n");
                free(data);
            }
            else if (dsize == 2)
            {
                printf("data type: short\n");
                short *data = (short *)malloc(sizeof(short) * nbEle);
                if (dorder == H5T_ORDER_LE)
                    status = H5Dread(dset, H5T_STD_I16LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                else
                    status = H5Dread(dset, H5T_STD_I16BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                int i;
                printf("reconstructed data = ");
                for (i = 0; i < 20; i++)
                    printf("%d ", data[i]);
                printf("\n");
                free(data);
            }
            else if (dsize == 4)
            {
                printf("data type: int\n");
                int *data = (int *)malloc(sizeof(int) * nbEle);
                if (dorder == H5T_ORDER_LE)
                    status = H5Dread(dset, H5T_STD_I32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                else
                    status = H5Dread(dset, H5T_STD_I32BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                int i;
                printf("reconstructed data = ");
                for (i = 0; i < 20; i++)
                    printf("%d ", data[i]);
                printf("\n");
                free(data);
            }
            else if (dsize == 8)
            {
                printf("data type: long\n");
                long *data = (long *)malloc(sizeof(long) * nbEle);
                if (dorder == H5T_ORDER_LE)
                    status = H5Dread(dset, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                else
                    status = H5Dread(dset, H5T_STD_I64BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
                int i;
                printf("reconstructed data = ");
                for (i = 0; i < 20; i++)
                    printf("%ld ", data[i]);
                printf("\n");
                free(data);
            }
        }

        break;
    default:
        printf("Error: H5Z-SZ supports only float, double or integers.\n");
        exit(0);
    }

    status = H5Pclose(dcpl);
    status = H5Dclose(dset);
    status = H5Fclose(file);

    return 0;
}
