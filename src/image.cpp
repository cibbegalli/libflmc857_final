/*
 *Created by Deangeli Gomes Neves
 *
 * This software may be freely redistributed under the terms
 * of the MIT license.
 *
 */
#include "image.h"

Image* createImage(int nx, int ny,int nchannels){
    Image* image = (Image*)calloc(1,sizeof(Image));
    image->nx = nx;
    image->ny = ny;
    image->nz = 1;
    image->numberPixels = image->nx*image->ny*image->nz;
    image->dx = 0;
    image->dy = 0;
    image->dz = 0;
    image->scalingFactor = 255;
    image->channelDepth = 8;
    image->unid[0] = 'm';
    image->unid[1] = 'm';
    image->channel = (float**)calloc(nchannels,sizeof(float*));
    image->nchannels = nchannels;
    findAppropriateColorSpace(image);
    image->dataTypeId = FLOAT;
    for (int i = 0; i < image->nchannels; ++i) {
        image->channel[i] = (float*)calloc(image->numberPixels,sizeof(float));
        //image->channel[i] = (float*)calloc(image->numberPixels,sizeof(float*));
        //image->channel->floatChannel[i] = (float*)calloc(image->numberPixels,sizeof(float));
    }
    return image;
}

Image* createImage(int nx, int ny,int nz, int nchannels){
    Image* image = (Image*)calloc(1,sizeof(Image));
    image->nx = nx;
    image->ny = ny;
    image->nz = nz;
    image->numberPixels = image->nx*image->ny*image->nz;
    image->dx = 0;
    image->dy = 0;
    image->dz = 0;
    image->scalingFactor = 255;
    image->channelDepth = 8;
    image->unid[0] = 'm';
    image->unid[1] = 'm';
    image->channel = (float**)calloc(nchannels,sizeof(float*));
    image->nchannels = nchannels;
    findAppropriateColorSpace(image);
    image->dataTypeId = FLOAT;
    for (int i = 0; i < image->nchannels; ++i) {
        image->channel[i] = (float*)calloc(image->numberPixels,sizeof(float));
        //image->channel[i] = (float*)calloc(image->numberPixels,sizeof(float*));
        //image->channel->floatChannel[i] = (float*)calloc(image->numberPixels,sizeof(float));
    }
    return image;
}

Image* createImage(int nx, int ny){
    Image* image = (Image*)calloc(1,sizeof(Image));
    image->nx = nx;
    image->ny = ny;
    image->numberPixels = nx*ny;
    image->dx = 0;
    image->dy = 0;
    image->scalingFactor = 255;
    image->channelDepth = 8;
    image->unid[0] = 'm';
    image->unid[1] = 'm';
    image->channel = (float**)calloc(1,sizeof(float*));
    image->nchannels = 1;
    image->colorSpace = GRAYSCALE;
    image->dataTypeId = FLOAT;
    for (int i = 0; i < image->nchannels; ++i) {
        image->channel[i] = (float*)calloc(image->numberPixels,sizeof(float*));
    }
    return image;
}

void findAppropriateColorSpace(Image* image){
    if (image->nchannels == 1){
        image->colorSpace = GRAYSCALE;
    }else if(image->nchannels == 3){
        image->colorSpace = RGB;
    }else if(image->nchannels == 4){
        image->colorSpace = RGBA;
    }else{
        image->colorSpace = UNKNOWN;
    }
}


void destroyImage(Image**image ){
    if((*image) == NULL){
        return;
    }
    for (int i = 0; i < (*image)->nchannels; ++i) {
        free((*image)->channel[i]);
    }

    free((*image)->channel);
    free((*image));
    (*image) = NULL;
}

void destroyImageVoidPointer(void* data){
    Image** aux = ((Image**)data);
    destroyImage(aux);
    *aux = NULL;
}

Image* readImage(char *filename){
    char ext[10],*pos;
    pos = strrchr(filename,'.') + 1;
    sscanf(pos,"%s",ext);
    Image* image = NULL;
    //gray image
    if (strcmp(ext,"pgm")==0){
        image = readImagePGM(filename);
    }else if(strcmp(ext,"ppm")==0){
        image = readImagePPM(filename);
    }else if(strcmp(ext,"png")==0){
        image = readImagePNG(filename);
    }else if( (strcmp(ext,"jpg")==0) || (strcmp(ext,"jpeg")==0) ){
        image = readImageJPEG(filename);
    }else if( (strcmp(ext,"tif")==0) || (strcmp(ext,"TIF")==0)){
        image = readImageTIFF(filename);
    }
    else{
        printf("image format unknown\n");
    }
    return image;
}

Image *readImagePGM(char *filename)
{
    FILE *fp=NULL;
    unsigned char *value=NULL;
    char type[10];
    int  maxValue,ncols,nrows,n;
    Image *img = NULL;
    char z[256];

    fp = fopen(filename,"rb");
    if (fp == NULL){
        fprintf(stderr,"Cannot open %s\n",filename);
        exit(-1);
    }
    fscanf(fp,"%s\n",type);
    NCFgets(z,255,fp);
    if((strcmp(type,"P5")==0)){
        int c = getc(fp);
        while (c == '#') {
            while (getc(fp) != '\n');
            c = getc(fp);
        }
        ungetc(c, fp);
        sscanf(z,"%d %d\n",&ncols,&nrows);
        n = ncols*nrows;
        NCFgets(z,255,fp);
        sscanf(z,"%d\n",&maxValue);
        fgetc(fp);
        value = (unsigned char *)calloc(n,sizeof(unsigned char));
        if (value != NULL){
            fread(value,sizeof(unsigned char),n,fp);
        }else{
            fprintf(stderr,"Insufficient memory in ReadImage\n");
            exit(-1);
        }
        fclose(fp);
        img = createImage(ncols,nrows,1);
        img->scalingFactor = maxValue;
        img->channelDepth = (log(maxValue+1)/log(2))+0.1;
        img->colorSpace = GRAYSCALE;
        for (int i=0; i < n; i++)
            img->channel[0][i]=(float)value[i];
        free(value);
    }else{
        if((strcmp(type,"P2")==0)){
            int c = getc(fp);
            while (c == '#') {
                while (getc(fp) != '\n');
                c = getc(fp);
            }
            ungetc(c, fp);
            sscanf(z,"%d %d\n",&ncols,&nrows);
            n = ncols*nrows;
            NCFgets(z,255,fp);
            sscanf(z,"%d\n",&maxValue);
            img = createImage(ncols,nrows,1);
            img->scalingFactor = maxValue;
            img->channelDepth = (log(maxValue+1)/log(1))+0.1;
            img->colorSpace = GRAYSCALE;

            for (int i=0; i < n; i++)
                fscanf(fp,"%f",&img->channel[0][i]);
            fclose(fp);
        }else{
            fprintf(stderr,"Input image must be P2 or P5\n");
            exit(-1);
        }
    }
    return(img);
}

Image *readImagePPM(char *filename){
    FILE *fp=NULL;
    char type[10];
    int  x,y,nx,ny,Imax;
    Image *I=NULL;
    ushort rgb16[3];
    char z[256];

    fp = fopen(filename,"r");
    if (fp != NULL) {
        fscanf(fp,"%s",type);
        NCFgets(z,255,fp);

        if (strcmp(type,"P6")==0){
            int c = getc(fp);
            while (c == '#') {
                while (getc(fp) != '\n');
                c = getc(fp);
            }
            ungetc(c, fp);

            if (fscanf(fp, "%d %d", &nx, &ny) != 2) {
                fprintf(stderr, "Invalid image size (error loading '%s')\n", filename);
                exit(1);
            }
            I = createImage(nx,ny,3);

            if (fscanf(fp, "%d", &Imax) != 1) {
                fprintf(stderr, "Invalid rgb component (error loading '%s')\n",
                        filename);
                exit(1);
            }

            I->scalingFactor = Imax;
            I->channelDepth = (log(Imax+1)/log(2))+0.1;
            I->colorSpace = RGB;
            while (fgetc(fp) != '\n');

            if (Imax <= 255){
                int iter = 0;
                for (y=0; y < ny; y++) {
                    for (x = 0; x < nx; x++) {
                        I->channel[0][iter] =  (float)fgetc(fp);//red
                        I->channel[1][iter] =  (float)fgetc(fp);//green
                        I->channel[2][iter] =  (float)fgetc(fp);//blue
                        iter++;
                    }
                }
            }else if (Imax > 255 && Imax <= 65535){

                int iter = 0;
                for (y=0; y < ny; y++) {
                    for (x=0; x < nx; x++){
                        if (fread(rgb16, 2, 3, fp) == 3) {
                            I->channel[0][iter]=(((rgb16[0] & 0xff) << 8) | ((ushort) rgb16[0] >> 8));
                            I->channel[1][iter]=(((rgb16[1] & 0xff) << 8) | ((ushort) rgb16[1] >> 8));
                            I->channel[2][iter]=(((rgb16[2] & 0xff) << 8) | ((ushort) rgb16[2] >> 8));
                        } else{
                            Error("Reading 16-bit error","ReadColorImage");
                        }
                        iter++;
                    }
                }
            }
        }else if(strcmp(type,"P3")==0){
            int c = getc(fp);
            while (c == '#') {
                while (getc(fp) != '\n');
                c = getc(fp);
            }
            ungetc(c, fp);
            if (fscanf(fp, "%d %d", &nx, &ny) != 2) {
                fprintf(stderr, "Invalid image size (error loading '%s')\n", filename);
                exit(1);
            }

            I = createImage(nx,ny,3);
            I->colorSpace = RGB;

            if (fscanf(fp, "%d", &Imax) != 1) {
                fprintf(stderr, "Invalid rgb component (error loading '%s')\n",
                        filename);
                exit(1);
            }

            I->scalingFactor = Imax;
            I->channelDepth = (log(Imax+1)/log(2))+0.1;

            int R=0,G=0,B=0;
            int iter=0;
            for (y=0; y < ny; y++) {
                for (x=0; x < nx; x++){
                    (void)fscanf(fp, "%d %d %d",&R,&G,&B);
                    I->channel[0][iter] = R;
                    I->channel[1][iter] = G;
                    I->channel[2][iter] = B;
                    iter++;
                }
            }
        }

        else{
            Error("ColorImage type invalid","ReadColormage");
        }

        fclose(fp);
    }else{
        Error(MSG2,"ReadColorImage");
    }

    return(I);

}

//TODO: suppport to 8-bit channel image
Image *readImageTIFF(char *filename){
    printf("Support to tiff images is unenabled\n");
    return NULL;
    // TIFF *tif=TIFFOpen(filename, "r");
    // uint32 width;
    // uint32 height;
    // uint32 npixels;
    // TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
    // TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);

    // npixels = width*height;
    // uint32 * raster=(uint32 *) _TIFFmalloc(npixels *sizeof(uint32));
    // int numberOfChannels = 4;

    // int status = TIFFReadRGBAImage(tif, width, height, raster, 0);
    // if(status==0){
    //     printf("[readImageTIFF] An error occured when reading image %s\n",filename);
    //     return NULL;
    // }

    // Image* image = createImage(width,height,numberOfChannels);
    // int k = 0;
    // for (int y = height-1; y >=0 ; --y) {
    //     for (int x = 0; x < ((int)width); ++x) {
    //         //printf("%u %u\n",x,y);
    //         int A = raster[k] >> 24;
    //         int B = (raster[k] >> 16) & 0x0000FF;
    //         int G = (raster[k] >> 8) & 0x0000FF;
    //         int R =  (raster[k]) & 0x0000FF;
    //         //printf("%d %d %d %d\n",R,G,B,A);
    //         imageValCh(image,x,y,0) = R;
    //         imageValCh(image,x,y,1) = G;
    //         imageValCh(image,x,y,2) = B;
    //         imageValCh(image,x,y,3) = A;
    //         k++;
    //     }
    // }

    // _TIFFfree(raster);
    // TIFFClose(tif);
    // return image;
}

Image *readImagePNG(char *filename){
/*
 * Copyright 2002-2010 Guillaume Cottenceau.
 *
 * This software may be freely redistributed under the terms
 * of the X11 license.
 *
 */
    Image* image = NULL;
    int width, height;
    png_byte color_type;
    png_byte bit_depth;

    png_structp png_ptr;
    png_infop info_ptr;

    //int number_of_passes;
    png_bytep * row_pointers;
    char header[8];    // 8 is the maximum size that can be checked

    /* open file and test for it being a png */
    FILE *fp = fopen(filename, "rb");
    if (!fp){
        printf("[read_png_file] File %s could not be opened for reading\n", filename);
        return image;
    }
    fread(header, 1, 8, fp);
    if (png_sig_cmp(((png_const_bytep)header), 0, 8)){
        printf("[read_png_file] File %s is not recognized as a PNG file\n",filename);
        return image;
    }


    /* initialize stuff */
    png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

    if (!png_ptr){
        printf("[read_png_file] png_create_read_struct failed\n");
        return image;
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr){
        printf("[read_png_file] png_create_info_struct failed\n");
        return image;
    }

    if (setjmp(png_jmpbuf(png_ptr))) {
        printf("[read_png_file] Error during init_io\n");
        return image;
    }

    png_init_io(png_ptr, fp);
    png_set_sig_bytes(png_ptr, 8);

    png_read_info(png_ptr, info_ptr);

    width = png_get_image_width(png_ptr, info_ptr);
    height = png_get_image_height(png_ptr, info_ptr);
    //color_type = png_get_color_type(png_ptr, info_ptr);
    bit_depth = png_get_bit_depth(png_ptr, info_ptr);

    //number_of_passes = png_set_interlace_handling(png_ptr);
    png_read_update_info(png_ptr, info_ptr);


    /* read file */
    if (setjmp(png_jmpbuf(png_ptr))){
        printf("[read_png_file] Error during read_image\n");
        return image;
    }

    row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);
    for (int y=0; y<height; y++) {
        row_pointers[y] = (png_byte *) malloc(png_get_rowbytes(png_ptr, info_ptr));
    }
    png_read_image(png_ptr, row_pointers);


    fclose(fp);

    ColorSpace colorSpace;
    unsigned int numberChannels = png_get_channels(png_ptr, info_ptr);
    color_type = png_get_color_type(png_ptr, info_ptr);
    if (color_type == PNG_COLOR_TYPE_GRAY){
        colorSpace = GRAYSCALE;
    }else if(color_type == PNG_COLOR_TYPE_GRAY_ALPHA){
        colorSpace = GRAYSCALE_ALPHA;
    }else if(color_type == PNG_COLOR_TYPE_RGB){
        colorSpace = RGB;
    }else if(color_type == PNG_COLOR_TYPE_RGB_ALPHA){
        colorSpace = RGBA;
    }else if(color_type == PNG_COLOR_TYPE_PALETTE){
        printf("[readImagePNG] PNG_COLOR_TYPE_PALETTE is not implemented yet\n");
        colorSpace = UNKNOWN;
    }else{
        printf("[readImagePNG] PNG_COLOR_TYPE_ unknown\n");
        colorSpace = UNKNOWN;
    }
    image = createImage(width,height,numberChannels);
    image->scalingFactor = pow(2,bit_depth)-1;
    image->channelDepth = bit_depth;
    image->colorSpace = colorSpace;
    for (int y=0; y<height; y++) {
        png_byte* row = row_pointers[y];
        for (int x=0; x<width; x++) {
            png_byte* ptr = &(row[x*numberChannels]);
            for (unsigned int c = 0; c < numberChannels; ++c) {
                imageValCh(image,x,y,c) = ptr[c];
            }
        }
    }

    for (int y=0; y<height; y++) {
        free(row_pointers[y]);
    }
    free(row_pointers);
    png_destroy_read_struct(&png_ptr, &info_ptr, NULL);

    return image;
}

Image *readImageJPEG(char *filename){
    //code from jpeg-9b/example.c
    Image* image = NULL;
    /* This struct contains the JPEG decompression parameters and pointers to
    * working space (which is allocated as needed by the JPEG library).
    */
    struct jpeg_decompress_struct cinfo;
    /* We use our private extension JPEG error handler.
    * Note that this struct must live as long as the main JPEG parameter
    * struct, to avoid dangling-pointer problems.
    */
    struct jpeg_error_mgr jerr;
    //struct my_error_mgr jerr;

    /* More stuff */
    FILE * infile;		/* source file */
    JSAMPARRAY buffer;		/* Output row buffer */
    int row_stride;		/* physical row width in output buffer */
    /* In this example we want to open the input file before doing anything else,
 * so that the setjmp() error recovery below can assume the file is open.
 * VERY IMPORTANT: use "b" option to fopen() if you are on a machine that
 * requires it in order to read binary files.
 */

    if ((infile = fopen(filename, "rb")) == NULL) {
        printf("[readImageJPEG] can't open %s\n",filename);
        return NULL;
    }

    /* Step 1: allocate and initialize JPEG decompression object */

    /* We set up the normal JPEG error routines, then override error_exit. */
    cinfo.err = jpeg_std_error(&jerr);
    //jerr.pub.error_exit = my_error_exit;

    /* Establish the setjmp return context for my_error_exit to use. */
    jmp_buf setjmp_buffer;
    if (setjmp(setjmp_buffer)) {
        /* If we get here, the JPEG code has signaled an error.
         * We need to clean up the JPEG object, close the input file, and return.
         */
        jpeg_destroy_decompress(&cinfo);
        printf("[readImageJPEG] code has signaled an error\n");
        fclose(infile);
        return NULL;
    }
    /* Now we can initialize the JPEG decompression object. */
    jpeg_create_decompress(&cinfo);
    /* Step 2: specify data source (eg, a file) */
    jpeg_stdio_src(&cinfo, infile);
    /* Step 3: read file parameters with jpeg_read_header() */
    (void) jpeg_read_header(&cinfo, TRUE);
    /* We can ignore the return value from jpeg_read_header since
     *   (a) suspension is not possible with the stdio data source, and
     *   (b) we passed TRUE to reject a tables-only JPEG file as an error.
     * See libjpeg.txt for more info.
     */

    /* Step 4: set parameters for decompression */

    /* In this example, we don't need to change any of the defaults set by
     * jpeg_read_header(), so we do nothing here.
     */

    /* Step 5: Start decompressor */
    (void) jpeg_start_decompress(&cinfo);
    /* We can ignore the return value since suspension is not possible
 * with the stdio data source.
 */

    /* We may need to do some setup of our own at this point before reading
     * the data.  After jpeg_start_decompress() we have the correct scaled
     * output image dimensions available, as well as the output colormap
     * if we asked for color quantization.
     * In this example, we need to make an output work buffer of the right size.
     */
    /* JSAMPLEs per row in output buffer */
    row_stride = cinfo.output_width * cinfo.output_components;
    /* Make a one-row-high sample array that will go away when done with image */
    buffer = (*cinfo.mem->alloc_sarray)
            ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);

    /* Step 6: while (scan lines remain to be read) */
    /*           jpeg_read_scanlines(...); */

    /* Here we use the library's state variable cinfo.output_scanline as the
     * loop counter, so that we don't have to keep track ourselves.
     */
    image = createImage(cinfo.output_width,cinfo.output_height,cinfo.output_components);
    image->channelDepth = cinfo.data_precision;
    image->scalingFactor = pow(2,image->channelDepth)-1;


    //printf("%d %d\n",cinfo.jpeg_color_space, cinfo.out_color_space);

    if(cinfo.out_color_space == JCS_GRAYSCALE){
        image->colorSpace = GRAYSCALE;
    }else if(cinfo.out_color_space == JCS_RGB){//maybe a lib bug?
        image->colorSpace = RGB;
    }else if(cinfo.out_color_space == JCS_YCbCr){//maybe a lib bug?
        image->colorSpace = YCbCr;
    }else if(cinfo.out_color_space == JCS_CMYK){
        image->colorSpace = CMYK;
    }else if(cinfo.out_color_space == JCS_YCCK){
        image->colorSpace = CMYK;
    } else if(cinfo.out_color_space == JCS_BG_RGB){
        printf("[readImageJPEG] big gamut R/G/B not implemented yet, image color space will be set to RGB\n");
        image->colorSpace = RGB;
        //return NULL;
    }else if(cinfo.out_color_space == JCS_BG_YCC){
        printf("[readImageJPEG] big gamut Y/Cb/Cr not implemented yet, image color space will be set to YCbCr\n");
        image->colorSpace = YCbCr;
        //return NULL;
    }else{
        printf("[readImageJPEG] unknown color space\n");
        image->colorSpace = UNKNOWN;
    }
    unsigned  int imageRow = 0;
    while (cinfo.output_scanline < cinfo.output_height){
        jpeg_read_scanlines(&cinfo, buffer, 1);
        unsigned int imageCol = 0;
        for (unsigned int i = 0; i < (unsigned int)cinfo.output_width; i++) {
            int shift = i*cinfo.num_components;
            for (int c = 0; c < cinfo.num_components; ++c) {
                imageValCh(image,imageCol,imageRow,c) = buffer[0][shift+c];
            }
            imageCol++;
        }
        imageRow++;
    }

    /* Step 7: Finish decompression */
    (void) jpeg_finish_decompress(&cinfo);

    /* This is an important step since it will release a good deal of memory. */
    jpeg_destroy_decompress(&cinfo);

    /* After finish_decompress, we can close the input file.
     * Here we postpone it until after no more JPEG errors are possible,
     * so as to simplify the setjmp error logic above.  (Actually, I don't
     * think that jpeg_destroy can do an error exit, but why assume anything...)
     */

    fclose(infile);
    /* At this point you may want to check to see whether any corrupt-data
     * warnings occurred (test whether jerr.pub.num_warnings is nonzero).
     */

    //jerr.num_warnings; //useful to know about corrupted data
    //printf("%ld\n",jerr.num_warnings);
    return image;

}

//TODO: suppport to 16-bit channel image
void writeImageTIFF(Image *image,char *filename){
    printf("Support to tiff images is unenabled\n");
    return;
//     TIFF *out= TIFFOpen(filename, "w");
//     int sampleperpixel = image->nchannels;    // or 3 if there is no alpha channel, you should get a understanding of alpha in class soon.

//     TIFFSetField (out, TIFFTAG_IMAGEWIDTH, image->nx);  // set the width of the image
//     TIFFSetField(out, TIFFTAG_IMAGELENGTH, image->ny);    // set the height of the image
//     TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, sampleperpixel);   // set number of channels per pixel
//     TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8);    // set the size of the channels
//     TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);    // set the origin of the image.
// //   Some other essential fields to set that you do not have to understand for now.
//     TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
//     TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);

//     tsize_t linebytes = sampleperpixel * image->nx;     // length in memory of one row of pixel in the image.
//     unsigned char *buf = NULL;        // buffer used to store the row of pixel information for writing to file
// //    Allocating memory to store the pixels of current row
//     buf =(unsigned char *)_TIFFmalloc(linebytes);
//     // We set the strip size of the file to be size of one row of pixels
//     TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(out, image->nx*sampleperpixel));

//     //Now writing image to the file one strip at a time
//     int k=0;
//     for (uint32 row = 0; row < (unsigned long)image->ny; row++)
//     {
//         k=0;
//         for (uint32 col = 0; col < (unsigned long)image->nx; ++col) {

//             for (int c = 0; c < sampleperpixel; ++c) {
//                 buf[k] = (unsigned char)imageValCh(image,col,row,c);
//                 k++;
//             }

//         }
//         //memcpy(buf, &image[(image->ny-row-1)*linebytes], linebytes);    // check the index here, and figure out why not using h*linebytes
//         if (TIFFWriteScanline(out, buf, row, 0) < 0){
//             break;
//         }
//     }
//     (void) TIFFClose(out);
//     if (buf)
//         _TIFFfree(buf);
}

void writeImageJPEG(Image *image,char *filename){
    //code from jpeg-9b/example.c
    /* This struct contains the JPEG compression parameters and pointers to
 * working space (which is allocated as needed by the JPEG library).
 * It is possible to have several such structures, representing multiple
 * compression/decompression processes, in existence at once.  We refer
 * to any one struct (and its associated working data) as a "JPEG object".
 */
    struct jpeg_compress_struct cinfo;

    /* This struct represents a JPEG error handler.  It is declared separately
 * because applications often want to supply a specialized error handler
 * (see the second half of this file for an example).  But here we just
 * take the easy way out and use the standard error handler, which will
 * print a message on stderr and call exit() if compression fails.
 * Note that this struct must live as long as the main JPEG parameter
 * struct, to avoid dangling-pointer problems.
 */
    struct jpeg_error_mgr jerr;

    /* More stuff */
    FILE * outfile;		/* target file */
    JSAMPARRAY buffer;
    /* Step 1: allocate and initialize JPEG compression object */

    /* We have to set up the error handler first, in case the initialization
     * step fails.  (Unlikely, but it could happen if you are out of memory.)
     * This routine fills in the contents of struct jerr, and returns jerr's
     * address which we place into the link field in cinfo.
     */
    cinfo.err = jpeg_std_error(&jerr);

    /* Now we can initialize the JPEG compression object. */
    jpeg_create_compress(&cinfo);
    /* Step 2: specify data destination (eg, a file) */
    /* Note: steps 2 and 3 can be done in either order. */

    /* Here we use the library-supplied code to send compressed data to a
     * stdio stream.  You can also write your own code to do something else.
     * VERY IMPORTANT: use "b" option to fopen() if you are on a machine that
     * requires it in order to write binary files.
     */
    if ((outfile = fopen(filename, "wb")) == NULL) {
        fprintf(stderr, "can't open %s\n", filename);
        exit(1);
    }
    jpeg_stdio_dest(&cinfo, outfile);

    /* First we supply a description of the input image.
* Four fields of the cinfo struct must be filled in:
*/
    cinfo.image_width = image->nx; 	/* image width and height, in pixels */
    cinfo.image_height = image->ny;
    cinfo.data_precision = image->channelDepth;

    if(image->colorSpace == RGB || image->colorSpace == RGBA){
        cinfo.input_components = 3;
        cinfo.in_color_space = JCS_RGB;
        cinfo.jpeg_color_space = JCS_RGB;
    }else if(image->colorSpace == GRAYSCALE || image->colorSpace == GRAYSCALE_ALPHA){
        cinfo.input_components = 1;
        cinfo.in_color_space = JCS_GRAYSCALE;
        cinfo.jpeg_color_space = JCS_GRAYSCALE;
    }else if(image->colorSpace == YCbCr){
        cinfo.input_components = 3;
        cinfo.in_color_space = JCS_YCbCr;
        cinfo.jpeg_color_space = JCS_YCbCr;
    }else if(image->colorSpace == CMYK){
        cinfo.input_components = 4;
        cinfo.in_color_space = JCS_CMYK;
        cinfo.jpeg_color_space = JCS_CMYK;
    }else if(image->colorSpace == YCbCrK){
        cinfo.input_components = 4;
        cinfo.in_color_space = JCS_YCCK;
        cinfo.jpeg_color_space = JCS_YCCK;
    }else{
        cinfo.input_components = image->nchannels;
        cinfo.in_color_space = JCS_UNKNOWN;
        cinfo.jpeg_color_space = JCS_UNKNOWN;
    }

    /* Now use the library's routine to set default compression parameters.
* (You must set at least cinfo.in_color_space before calling this,
* since the defaults depend on the source color space.)
*/
    jpeg_set_defaults(&cinfo);

    /* Now you can set any non-default parameters you wish to.
* Here we just illustrate the use of quality (quantization table) scaling:
*/
    int quality = 100;
    jpeg_set_quality(&cinfo, quality, TRUE /* limit to baseline-JPEG values */);

    /* Step 4: Start compressor */
    /* TRUE ensures that we will write a complete interchange-JPEG file.
     * Pass TRUE unless you are very sure of what you're doing.
     */
    jpeg_start_compress(&cinfo, TRUE);

    /* Step 5: while (scan lines remain to be written) */
    /*           jpeg_write_scanlines(...); */

    /* Here we use the library's state variable cinfo.next_scanline as the
     * loop counter, so that we don't have to keep track ourselves.
     * To keep things simple, we pass one scanline per call; you can pass
     * more if you wish, though.
     */
    int row_stride = cinfo.image_width * cinfo.num_components;
    buffer = (*cinfo.mem->alloc_sarray)
            ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);
    unsigned int imageRow = 0;
    while (cinfo.next_scanline < cinfo.image_height) {
        /* jpeg_write_scanlines expects an array of pointers to scanlines.
         * Here the array is only one element long, but you could pass
         * more than one scanline at a time if that's more convenient.
         */
        unsigned int imageCol = 0;
        for (unsigned int i = 0; i < (unsigned int)cinfo.image_width; i++) {
            int shift = i*cinfo.num_components;
            for (int c = 0; c < cinfo.num_components; ++c) {
                buffer[0][(shift+c)] = (unsigned char)imageValCh(image,imageCol,imageRow,c);
            }
            imageCol++;
        }
        imageRow++;
        (void) jpeg_write_scanlines(&cinfo, buffer, 1);
    }

    /* Step 6: Finish compression */

    jpeg_finish_compress(&cinfo);
    /* After finish_compress, we can close the output file. */
    fclose(outfile);

    /* Step 7: release JPEG compression object */

    /* This is an important step since it will release a good deal of memory. */
    jpeg_destroy_compress(&cinfo);
}

void writeImagePNG(Image *image,char *filename){
/*
 * Copyright 2002-2010 Guillaume Cottenceau.
 *
 * This software may be freely redistributed under the terms
 * of the X11 license.
 *
 */
    /* create file */
    FILE *fp = fopen(filename, "wb");
    if (!fp){
        printf("[write_png_file] File %s could not be opened for writing\n", filename);
    }

    /* initialize stuff */
    png_structp png_ptr;
    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

    if (!png_ptr){
        printf("[write_png_file] png_create_write_struct failed\n");
    }

    png_infop info_ptr;
    info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr){
        printf("[write_png_file] png_create_info_struct failed\n");
    }

    if (setjmp(png_jmpbuf(png_ptr))){
        printf("[write_png_file] Error during init_io\n");
    }

    png_init_io(png_ptr, fp);


    /* write header */
    if (setjmp(png_jmpbuf(png_ptr))){
        printf("[write_png_file] Error during writing header\n");
    }
    png_uint_32 width, height;
    width = image->nx;
    height = image->ny;
    png_byte bit_depth = image->channelDepth;
    png_byte color_type;
    if(image->colorSpace == RGB) {
        color_type = PNG_COLOR_TYPE_RGB;
    } else if(image->colorSpace == RGBA) {
        color_type = PNG_COLOR_TYPE_RGB_ALPHA;
    }else if(image->colorSpace == GRAYSCALE){
        color_type = PNG_COLOR_TYPE_GRAY;
    }else if(image->colorSpace == GRAYSCALE_ALPHA){
        color_type = PNG_COLOR_TYPE_GRAY_ALPHA;
    }else{
        if(image->nchannels == 1){
            color_type = PNG_COLOR_TYPE_GRAY;
        }
        else if(image->nchannels == 2){
            color_type = PNG_COLOR_TYPE_GRAY_ALPHA;
        }
        else if(image->nchannels == 3){
            color_type = PNG_COLOR_TYPE_RGB;
        }
        else if(image->nchannels == 4){
            color_type = PNG_COLOR_TYPE_RGB_ALPHA;
        }else{
            printf("[writeImagePNG] invalide color space: %d\n",image->colorSpace);
            return;
        }
    }

    png_set_IHDR(png_ptr, info_ptr, width, height,
                 bit_depth, color_type, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    png_write_info(png_ptr, info_ptr);


    /* write bytes */
    if (setjmp(png_jmpbuf(png_ptr))){
        printf("[write_png_file] Error during writing bytes\n");
    }

    png_bytep * row_pointers;
    row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);
    for (size_t y=0; y<height; y++) {
        row_pointers[y] = (png_byte *) malloc(png_get_rowbytes(png_ptr, info_ptr));
    }

    for (size_t y=0; y<height; y++) {
        png_byte* row = row_pointers[y];
        for (size_t x=0; x<width; x++) {
            png_byte* ptr = &(row[x*image->nchannels]);
            for (int c = 0; c < image->nchannels; ++c) {
                ptr[c] = (png_byte)(imageValCh(image,x,y,c));
            }
        }
    }

    png_write_image(png_ptr, row_pointers);

    /* end write */
    if (setjmp(png_jmpbuf(png_ptr))){
        printf("[write_png_file] Error during end of write");
    }

    png_write_end(png_ptr, NULL);

    /* cleanup heap allocation */
    for (size_t y=0; y<height; y++){
        free(row_pointers[y]);
    }
    free(row_pointers);

    fclose(fp);
    png_destroy_write_struct(&png_ptr,&info_ptr);
}

void writeImagePGM(Image *image,char *filename, char *magicNumber){

    if (strcmp(magicNumber,"P2")==0){
        writeImageP2(image,filename);
    }
    else if(strcmp(magicNumber,"P5")==0){
        writeImageP5(image,filename);
    }else{
        printf("unknown magic number");
    }
}

void writeImagePPM(Image *image,char *filename, char *magicNumber){

    if (strcmp(magicNumber,"P3")==0){
        writeImageP3(image,filename);
    }
    else if(strcmp(magicNumber,"P6")==0){
        writeImageP6(image,filename);
    }else{
        printf("unknown magic number");
    }
}


void writeImageP2(Image *image,char *filename){
    int i, j;
    FILE *fp = fopen(filename, "w");
    (void) fprintf(fp, "P2\n%d %d\n%d\n", image->nx,image->ny, image->scalingFactor);
    int k = 0;


    float* channel0 = image->channel[0];
    for (i = 0; i < image->ny; i++){
        for (j = 0; j < image->nx; j++){
            (void) fprintf(fp,"%d ",(int)channel0[k]);
            k++;
        }
        (void) fprintf(fp,"\n");
    }


    (void) fclose(fp);
}

void writeImageP3(Image *image,char *filename){
    int i, j;
    FILE *fp = fopen(filename, "w");
    (void) fprintf(fp, "P3\n%d %d\n%d\n", image->nx,image->ny, image->scalingFactor);
    int k=0;

    float* channel0 = image->channel[0];
    float* channel1 = image->channel[1];
    float* channel2 = image->channel[2];
    for (i = 0; i < image->ny; i++){
        for (j = 0; j < image->nx; j++){
            (void) fprintf(fp,"%d %d %d ", (int)(channel0[k]),(int)channel1[k],(int)channel2[k]);
            k++;
        }
        (void) fprintf(fp,"\n");
    }


    (void) fclose(fp);

}

void writeImageP5(Image *image,char *filename){
    FILE *fp = fopen(filename, "wb");
    (void) fprintf(fp, "P5\n%d %d\n%d\n", image->nx,image->ny, image->scalingFactor);
    int k = 0;

    float* channel0 = image->channel[0];
    for (int i = 0; i < image->ny; i++){
        for (int j = 0; j < image->nx; j++){
            static unsigned char gray[1];
            gray[0] = (int)channel0[k];
            (void) fwrite(gray, 1, 1, fp);
            k++;
        }
    }

    (void) fclose(fp);
}


void writeImageP6(Image *image,char *filename){
    int i, j;
    FILE *fp = fopen(filename, "wb");
    (void) fprintf(fp, "P6\n%d %d\n%d\n", image->nx,image->ny, image->scalingFactor);
    int k=0;

    float* channel0 = image->channel[0];
    float* channel1 = image->channel[1];
    float* channel2 = image->channel[2];
    for (i = 0; i < image->ny; i++){
        for (j = 0; j < image->nx; j++){
            static unsigned char color[3];
            color[0] = channel0[k];
            color[1] = channel1[k];
            color[2] = channel2[k];
            (void) fwrite(color, 1, 3, fp);
            k++;
        }
    }

    (void) fclose(fp);
}

void writeImage(Image* image,char *filename){
    char ext[10],*pos;
    pos = strrchr(filename,'.') + 1;
    sscanf(pos,"%s",ext);

    if (strcmp(ext,"pgm")==0){
        writeImagePGM(image,filename,"P5");
    }
    else if (strcmp(ext,"ppm")==0){
        writeImagePPM(image,filename,"P6");
    }else if(strcmp(ext,"png")==0){
        writeImagePNG(image,filename);
    }else if(strcmp(ext,"jpg")==0){
        writeImageJPEG(image,filename);
    } else if(strcmp(ext,"jpeg")==0){
        writeImageJPEG(image,filename);
    }
    else{
        printf("[writeImage] unsuported format\n");
    }
}



Image *convertRGBtoYCbCr(Image *rgbImage)
{
    Image *ycbcrImage;
    float a = 16.0;
    float b = 128.0;

    ycbcrImage = createImage(rgbImage->nx,rgbImage->ny,3);
    ycbcrImage->scalingFactor = rgbImage->scalingFactor;
    for (int k=0; k < rgbImage->numberPixels; k++){
        float lum  = (0.257*rgbImage->channel[0][k])+
                     (0.504*rgbImage->channel[1][k])+
                     (0.098*rgbImage->channel[2][k])+a;
        float Cb = (-0.148*rgbImage->channel[0][k]) +
                   (-0.291*rgbImage->channel[1][k]) +
                   (0.439*rgbImage->channel[2][k])+b;
        float Cr = (0.439*rgbImage->channel[0][k]) +
                   (-0.368*rgbImage->channel[1][k]) +
                   (-0.071*rgbImage->channel[2][k]) +b;

        if (lum < 0)   lum = 0.0;
        if (lum > ycbcrImage->scalingFactor) lum = ycbcrImage->scalingFactor;
        if (Cb < 0)   Cb = 0.0;
        if (Cb > ycbcrImage->scalingFactor) Cb = ycbcrImage->scalingFactor;
        if (Cr < 0)   Cr = 0.0;
        if (Cr > ycbcrImage->scalingFactor) Cr = ycbcrImage->scalingFactor;

        ycbcrImage->channel[0][k] = lum;
        ycbcrImage->channel[1][k] = Cb;
        ycbcrImage->channel[2][k] = Cr;
    }

    return(ycbcrImage);
}


float sumUpAllPixelsValues(Image *image, bool normalize){
    double sum = 0;
    for (int k = 0; k < image->numberPixels; ++k) {
        for (int c = 0; c < image->nchannels; ++c) {
            sum += image->channel[c][k];
        }
    }
    if(normalize){
        sum /= ( image->numberPixels * image->nchannels);
    }
    return ((float)sum);
}

Image* extractImageChannel(Image* image, int channel){
    Image* channelImage = createImage(image->nx,image->ny,1);
    channelImage->scalingFactor = image->scalingFactor;
    for (int k = 0; k < image->numberPixels; ++k) {
        channelImage->channel[0][k] = image->channel[channel][k];
    }
    return channelImage;
}

void printImage(Image* image){
    int index=0;
    for (int y = 0; y < image->ny; ++y) {
        for (int x = 0; x < image->nx; ++x) {
            for (int c = 0; c < image->nchannels; ++c) {
                index = (y*image->nx)+x;
                printf("%.2f ",image->channel[c][index]);
            }
        }
        printf("\n");
    }
}

Image* createAlphaChannel(Image* image,float alpha){
    Image *output = createImage(image->nx,image->ny,image->nchannels+1);
    output->scalingFactor = image->scalingFactor;
    for (int k = 0; k < image->numberPixels; ++k) {
        for (int c = 0; c < image->nchannels; ++c) {
            output->channel[c][k] = image->channel[c][k];
        }
        output->channel[image->nchannels][k] = alpha;
    }
    return output;
}

void addUniformNoise(Image* image, float uniformValue, double probability){
    for (int k = 0; k < image->numberPixels; ++k) {
        if(randomNormalized() <=  probability){
            for (int c = 0; c < image->nchannels; ++c) {
                image->channel[c][k] += uniformValue;
                if(image->channel[c][k] > image->scalingFactor){
                    image->channel[c][k] = image->scalingFactor;
                }
            }
        }
    }
}

void addSaltAndPepperNoise(Image* image, double probability){
    for (int k = 0; k < image->numberPixels; ++k) {
        if(randomNormalized() <=  probability){
            for (int c = 0; c < image->nchannels; ++c) {
                //salt
                if(randomNormalized() >= 0.5){
                    image->channel[c][k] = image->scalingFactor;
                }else{
                    image->channel[c][k] = 0;
                }
            }
        }

    }
}

void addAdditiveGaussianNoise(Image* image, double mean, double variance){


    double standardDeviation = sqrt(variance);
    for (int y = 0; y < image->ny; ++y) {
        for (int x = 0; x < image->nx; ++x) {
            double value = generateGaussianNoise(0,standardDeviation);
            value = value*image->scalingFactor;
            imageVal(image,x,y) += value;
            if(imageVal(image,x,y) < 0){
                imageVal(image,x,y) = 0;
            }
            if(imageVal(image,x,y) > image->scalingFactor){
                imageVal(image,x,y) = image->scalingFactor;
            }
        }
    }
}

uint8_t* convertImage2IntergerArray8bits(Image* imageRGBA){
    uint8_t* output = (uint8_t*)calloc(imageRGBA->nchannels*imageRGBA->nx*imageRGBA->ny ,sizeof(uint8_t));
    long index = 0;
    for (int k = 0; k < imageRGBA->numberPixels; ++k) {
        for (int c = 0; c < imageRGBA->nchannels; ++c) {
            output[index] = (uint8_t)imageRGBA->channel[c][k];
            index++;
        }
    }
    return output;
}

Image* copyImage(Image* image, bool copyChannels){
    Image* imageCopy = createImage(image->nx,image->ny,image->nz,image->nchannels);
    imageCopy->colorSpace = image->colorSpace;
    imageCopy->dx = image->dx;
    imageCopy->dy = image->dy;
    imageCopy->dz = image->dz;
    imageCopy->dataTypeId = image->dataTypeId;
    imageCopy->scalingFactor = image->scalingFactor;
    for (int i = 0; i < 10; ++i) {
        imageCopy->unid[i] = image->unid[i];
    }
    if(copyChannels) {
        for (int k = 0; k < image->numberPixels; ++k) {
            for (int c = 0; c < image->nchannels; ++c) {
                imageCopy->channel[c][k] = image->channel[c][k];
            }
        }
    }
    return imageCopy;
}

Image *imageSubtraction(Image *image1, Image *image2, bool saturation){
    Image *outputImage = NULL;
    if(isImagesSameDomain(image1,image2)){
        outputImage = createImage(image1->nx,image1->ny,image1->nchannels);
        int n = image1->nx*image1->ny;
        int result;
        for (int i = 0; i < n; ++i) {
            for(int c=0; c<image1->nchannels; c++){
                result = image1->channel[c][i] - image2->channel[c][i];
                if(saturation){
                    result = (result<0)?0:result;
                }
                outputImage->channel[c][i] = result;
            }

        }
    }
    return outputImage;
}

bool isImagesSameDomain(Image* image1, Image* image2){
    if(image1->nx == image2->nx){
        if(image1->ny == image2->ny){
            if(image1->nz == image2->nz){
                return true;
            }
        }
    }
    return false;
}

Image* packImagesFromDirectory(DirectoryManager* directoryManager){
    Image* imagePack = NULL;
    Image* currentImage = NULL;
    Image* merged = NULL;
    for (int i = 0; i < (int)directoryManager->nfiles; ++i) {
        currentImage = readImage(directoryManager->files[i]->path);
        merged = mergeImages(imagePack,currentImage);
        destroyImage(&currentImage);
        destroyImage(&imagePack);
        imagePack = copyImage(merged,true);
        destroyImage(&merged);
    }
    return imagePack;
}

Image* convertGrayImage2RGBImage(Image* image_ppm){
    Image* image_rgb = createImage(image_ppm->nx,image_ppm->ny,3);
    for (int k = 0; k < image_rgb->numberPixels; ++k) {
        image_rgb->channel[0][k] = image_ppm->channel[0][k];
        image_rgb->channel[1][k] = image_ppm->channel[0][k];
        image_rgb->channel[2][k] = image_ppm->channel[0][k];
    }
    return image_rgb;
}

Image* mergeImages(Image* image1, Image* image2){
    Image* imageMerged = NULL;
    if(image1 == NULL && image2 == NULL){
        printf("NULL images\n");
        return NULL;
    }
    if(image1 == NULL){
        imageMerged = copyImage(image2,true);
        return imageMerged;
    }
    if(image2 == NULL){
        imageMerged = copyImage(image1,true);
        return imageMerged;
    }

    if(image1->nx != image2->nx || image1->ny != image2->ny){
        printf("images dimensions mismatch\n");
        return NULL;
    }
    if(image1->nchannels != image2->nchannels){
        printf("number of channels mismatch\n");
        return NULL;
    }
    int totalSlices = image1->nz+image2->nz;
    imageMerged = createImage(image1->nx,image1->ny,totalSlices,image1->nchannels);

#pragma omp parallel for
    for (int k = 0; k < image1->numberPixels; ++k) {
        for (int c = 0; c < image1->nchannels; ++c) {
            imageMerged->channel[c][k] = image1->channel[c][k];
        }
    }

    for (int k = 0,kShift = image1->numberPixels; k < image2->numberPixels; ++k,++kShift) {
        for (int c = 0; c < image2->nchannels; ++c) {
            imageMerged->channel[c][kShift] = image2->channel[c][k];
        }
    }
    return imageMerged;
}

Image* getSlice(Image* image, int sliceIndex){
    if(sliceIndex >= image->nz || sliceIndex < 0){
        printf("invalid sliceIndex (z-coordinate) value\n");
        return NULL;
    }
    Image* slice = createImage(image->nx,image->ny,image->nchannels);
    for (int y = 0; y < image->ny; ++y) {
        for (int x = 0; x < image->nx; ++x) {
            for (int c = 0; c < image->nchannels; ++c) {
                imageValCh(slice,x,y,c) = imageVolumeCh(image,x,y,sliceIndex,c);
            }
        }
    }
    return slice;
}

Image* extractSubImage(Image*image, int xCoord,int yCoord, int xsize, int ysize, bool zeroPadding){
    Image *subImage = NULL;
    int imageX = xCoord+xsize;
    int imageY = yCoord+ysize;
    int subImageSizeX = 0;
    int subImageSizeY = 0;

    subImageSizeX = xsize;
    subImageSizeY = ysize;

    if(image->nx <= xCoord && image->ny <= yCoord){
        printf("Warning: xCoord and yCoord are outbound.\n");
    }

    if(image->nx < imageX){
        if(!zeroPadding){
            subImageSizeX = image->nx-xCoord;
            if(subImageSizeX <= 0){
                printf("Impossible create a patch without zero padding\n");
                return NULL;
            }
        }
    }
    if(image->ny < imageY){
        if(!zeroPadding){
            subImageSizeY = image->ny-yCoord;
            if(subImageSizeY <= 0){
                printf("Impossible create a patch without zero padding\n");
                return NULL;
            }
        }
    }

    subImage = createImage(subImageSizeX,subImageSizeY,image->nchannels);
    for (int y = 0; y < subImage->ny; ++y) {
        for (int x = 0; x < subImage->nx; ++x) {
            for (int c = 0; c < subImage->nchannels; ++c) {
                if(isValidPixelCoordinate(image,xCoord+x,yCoord+y)){
                    imageValCh(subImage,x,y,c) = imageValCh(image,(xCoord+x),(yCoord+y),c);
                }
            }
        }
    }
    return subImage;
}

bool isValidPixelCoordinate(Image *image,int pixelCoordinateX,int pixelCoordinateY){
    if(pixelCoordinateX < 0 || pixelCoordinateY < 0){
        return false;
    }
    if(pixelCoordinateX >= image->nx || pixelCoordinateY >= image->ny){
        return false;
    }
    return true;
}

void putSlice(Image* volume, Image* image, int sliceIndex){
    if( (volume->nx != image->nx) || (volume->ny != image->ny)){
        printf("[putSlice] images dimension mismatch\n");
    }

    if( (volume->nchannels != image->nchannels)){
        printf("[putSlice] number of channels mismatch\n");
    }

    if(sliceIndex < 0 || sliceIndex >= volume->nz){
        printf("[putSlice] invalid sliceIndex\n");
    }
    for (int y = 0; y < image->ny; ++y) {
        for (int x = 0; x < image->nx; ++x) {
            for (int c = 0; c < volume->nchannels; ++c) {
                imageVolumeCh(volume,x,y,sliceIndex,c) = imageValCh(image,x,y,c);
            }
        }
    }
}
