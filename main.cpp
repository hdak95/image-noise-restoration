#include<iostream>
#include<fstream>
#include<cstdlib>
#include<ctime>
#include<algorithm>
#include<cmath>
#include<string>
#include<random>

using namespace std;

#define height 512
#define width 512
#define inf 9999
#define pi 3.14

float PSNR(unsigned char* src, unsigned char* dst);//PSNR calculator
void gaussianNoise(unsigned char* img,int mean, int std);//gaussian noise generator
void saltPepperNoise(unsigned char* img,int prob);//salt & pepper noise generator
void gaussianFilter(unsigned char* src, unsigned char* dst, int size, int std);//gaussian filter
void medianFilter(unsigned char* src, unsigned char* dst,int size);//median filter
void alphaTrimmedMeanFilter(unsigned char* src, unsigned char* dst, int size, double a);//alpha trimmed mean filter
void padding(unsigned char* img,unsigned char* pad, int size);//copy padding

int main(){
	ifstream fInput;//file reader
	fInput.open("./lena(512x512).raw", ios::binary);
	unsigned char* grayBuffer = new unsigned char[height * width];
	fInput.read((char*)grayBuffer, height * width);
	fInput.close();

	//noise parameter
	int gNoise[3][2] = { {0,5},{0,10},{0,15} };
	int spNoise[3] = { 20,30,40 };

	//filter parameter
	int gFilter[3][2] = { {3,1},{9,3}, {15,5} };
	int mFilter[2] = { 3,5 };
	double aFilter[4][2] = { {3,0.2},{3,0.4},{5,0.1},{5,0.4} };

	cout << fixed;
	cout.precision(4);//psnr 4 decimal point

	ofstream fOutput;//file writer

	//gaussian noise + gaussian filter
	for (int i = 0; i < 3; i++) {
		unsigned char* gaussian = new unsigned char[height * width];
		copy(grayBuffer, grayBuffer + height * width, gaussian);
		gaussianNoise(gaussian, gNoise[i][0], gNoise[i][1]);	//noise generator

		//image restoration
		for (int j = 0; j < 3; j++) {
			unsigned char* gPadding = new unsigned char[(height + gFilter[j][0] - 1) * (width + gFilter[j][0] - 1)];
			unsigned char* dst = new unsigned char[height * width];
			padding(gaussian, gPadding, gFilter[j][0]);
			gaussianFilter(gPadding, dst, gFilter[j][0], gFilter[j][1]);
			delete[] gPadding;
			cout << "gaussian noise mean " << gNoise[i][0]<< " std " << gNoise[i][1] << " gaussian filter " << gFilter[j][0] << "x" << gFilter[j][0] << " std = " << gFilter[j][1] << " PSNR " << PSNR(grayBuffer, dst) << " dB" << endl;
			fOutput.open("./gaussian_mean" + to_string(gNoise[i][0]) + "std" + to_string(gNoise[i][1]) + "_gaussian_filter_size" + to_string(gFilter[j][0]) + "_std" + to_string(gFilter[j][1]) + "(512x512).raw", ios::binary);
			fOutput.write((char*)dst, height * width);
			fOutput.close();
			delete[] dst;
		}
	}

	//gaussian noise + median filter
	for (int i = 0; i < 3; i++) {
		unsigned char* gaussian = new unsigned char[height * width];
		copy(grayBuffer, grayBuffer + height * width, gaussian);
		gaussianNoise(gaussian, gNoise[i][0], gNoise[i][1]);	//noise generator

		//image restoration
		for (int j = 0; j < 2; j++) {
			unsigned char* gPadding = new unsigned char[(height + mFilter[j] - 1) * (width + mFilter[j] - 1)];
			unsigned char* dst = new unsigned char[height * width];
			padding(gaussian, gPadding, mFilter[j]);
			medianFilter(gPadding, dst, mFilter[j]);
			delete[] gPadding;
			cout << "gaussian noise mean " << gNoise[i][0] << " std " << gNoise[i][1] << " median filter " << mFilter[j] << "x" << mFilter[j] << " PSNR " << PSNR(grayBuffer, dst) << "dB" << endl;
			fOutput.open("./gaussian_mean" + to_string(gNoise[i][0]) + "std" + to_string(gNoise[i][1]) + "_median_filter_size" + to_string(mFilter[j]) + "(512x512).raw", ios::binary);
			fOutput.write((char*)dst, height * width);
			fOutput.close();
			delete[] dst;
		}
	}

	//gaussian noise + alpha-trimmed mean filter
	for (int i = 0; i < 3; i++) {
		unsigned char* gaussian = new unsigned char[height * width];
		copy(grayBuffer, grayBuffer + height * width, gaussian);
		gaussianNoise(gaussian, gNoise[i][0], gNoise[i][1]);	//noise generator

		//image restoration
		for (int j = 0; j < 4; j++) {
			unsigned char* gPadding = new unsigned char[(height + int(aFilter[j][0]) - 1) * (width + int(aFilter[j][0]) - 1)];
			unsigned char* dst = new unsigned char[height * width];
			padding(gaussian, gPadding, int(aFilter[j][0]));
			alphaTrimmedMeanFilter(gPadding, dst, int(aFilter[j][0]), aFilter[j][1]);
			delete[] gPadding;
			cout << "gaussian noise mean " << gNoise[i][0] << " std " << gNoise[i][1] << " alpha-trimmed mean filter " << int(aFilter[j][0]) << "x" << int(aFilter[j][0]) << " a = 0." << int(10 * aFilter[j][1]) << " PSNR " << PSNR(grayBuffer, dst) << "dB" << endl;
			fOutput.open("./gaussian_mean" + to_string(gNoise[i][0]) + "std" + to_string(gNoise[i][1]) + "_alpha-trimmed_mean_filter_size" + to_string(int(aFilter[j][0])) + "_alpha0." + to_string(int(aFilter[j][1] * 10)) + "(512x512).raw", ios::binary);
			fOutput.write((char*)dst, height * width);
			fOutput.close();
			delete[] dst;
		}
	}

	//salt & pepper noise + gaussian filter
	for (int i = 0; i < 3; i++) {
		unsigned char* saltpepper = new unsigned char[height * width];
		copy(grayBuffer, grayBuffer + height * width, saltpepper);
		saltPepperNoise(saltpepper, spNoise[i]);	//noise generator

		//image restoration
		for (int j = 0; j < 3; j++) {
			unsigned char* spPadding = new unsigned char[(height + gFilter[j][0] - 1) * (width + gFilter[j][0] - 1)];
			unsigned char* dst = new unsigned char[height * width];
			padding(saltpepper, spPadding, gFilter[j][0]);
			gaussianFilter(spPadding, dst, gFilter[j][0], gFilter[j][1]);
			delete[] spPadding;
			cout << "salt & pepper noise " << spNoise[i] << "% gaussian filter " << gFilter[j][0] << "x" << gFilter[j][0] << " std = " << gFilter[j][1] << " PSNR " << PSNR(grayBuffer, dst) << " dB" << endl;
			fOutput.open("./salt&pepper_" + to_string(spNoise[i]) + "_gaussian_filter_size" + to_string(gFilter[j][0]) + "_std" + to_string(gFilter[j][1]) + "(512x512).raw", ios::binary);
			fOutput.write((char*)dst, height * width);
			fOutput.close();
			delete[] dst;
		}
	}

	//salt & pepper noise + median filter
	for (int i = 0; i < 3; i++) {
		unsigned char* saltpepper = new unsigned char[height * width];
		copy(grayBuffer, grayBuffer + height * width, saltpepper);
		saltPepperNoise(saltpepper, spNoise[i]);	//noise generator
		
		//image restoration
		for (int j = 0; j < 2; j++) {
			unsigned char* spPadding = new unsigned char[(height + mFilter[j] - 1) * (width + mFilter[j] - 1)];
			unsigned char* dst = new unsigned char[height * width];
			padding(saltpepper, spPadding, mFilter[j]);
			medianFilter(spPadding, dst, mFilter[j]);
			delete[] spPadding;
			cout << "salt & pepper noise " << spNoise[i] << "% median filter " << mFilter[j] << "x" << mFilter[j] << " PSNR " << PSNR(grayBuffer, dst) << "dB" << endl;
			fOutput.open("./salt&pepper_" + to_string(spNoise[i]) + "_median_filter_size" + to_string(mFilter[j]) + "(512x512).raw", ios::binary);
			fOutput.write((char*)dst, height * width);
			fOutput.close();
			delete[] dst;
		}
	}
	
	//salt & pepper noise + alpha-trimmed mean filter
	for (int i = 0; i < 3; i++) {
		unsigned char* saltpepper = new unsigned char[height * width];
		copy(grayBuffer, grayBuffer + height * width, saltpepper);
		saltPepperNoise(saltpepper, spNoise[i]);	//noise generator
		
		//image restoration
		for (int j = 0; j < 4; j++) {
			unsigned char* spPadding = new unsigned char[(height + int(aFilter[j][0]) - 1) * (width + int(aFilter[j][0]) - 1)];
			unsigned char* dst = new unsigned char[height * width];
			padding(saltpepper, spPadding, int(aFilter[j][0]));
			alphaTrimmedMeanFilter(spPadding, dst, int(aFilter[j][0]), aFilter[j][1]);
			delete[] spPadding;
			cout << "salt & pepper noise " << spNoise[i] << "% alpha-trimmed mean filter " << int(aFilter[j][0]) << "x" << int(aFilter[j][0]) << " a = 0." << int(10 * aFilter[j][1]) << " PSNR " << PSNR(grayBuffer, dst) << "dB" << endl;
			fOutput.open("./salt&pepper_" + to_string(spNoise[i]) + "_alpha-trimmed_mean_filter_size" + to_string(int(aFilter[j][0])) + "_alpha0." + to_string(int(aFilter[j][1] * 10)) + "(512x512).raw", ios::binary);
			fOutput.write((char*)dst, height * width);
			fOutput.close();
			delete[] dst;
		}
	}
	
	delete[] grayBuffer;

	return 0;
}


float PSNR(unsigned char* src,unsigned char* dst){
	int sum = 0;
	double MSE,psnr;
	for (int i = 0; i < width * height; i++)
		sum += pow(src[i] - dst[i],2);
	MSE = sum / (width * height);
	if (MSE == 0)	//original image
		return inf;
	else
		return 10 * log10(65025 / MSE); //dB value
}


void gaussianNoise(unsigned char* img, int mean, int std){
	unsigned int seed = static_cast<unsigned int>(time(NULL));//normal distribution time seed
	default_random_engine generator(seed);
	normal_distribution<double> distribution(mean, std);
	
	for (int i = 0; i < height * width; i++) {
		double temp = img[i] + distribution(generator);//original image + noise
		if (temp > 255)
			img[i] = 255;
		else if (temp < 0)
			img[i] = 0;
		else
			img[i] = temp;
	}
}


void saltPepperNoise(unsigned char* img,int prob){
	srand(time(NULL));//random seed
	
	for (int i = 0; i < width * height; i++) {
		if (rand() % 100 < prob) {//salt and peeper addition condition
			img[i] = 255 * (rand() % 2); //remainder 0->0 / 1->255
		}
	}
}


void gaussianFilter(unsigned char* src,unsigned char* dst,int size, int std){
	double* filter = new double[size * size]; //gaussian mask
	double sum = 0;
	
	//2d gaussian mask
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			filter[i * size + j] = exp(-(pow((j - int(size / 2)), 2) + pow((i - int(size / 2)), 2)) / (2 * std * std)) / (2 * pi * std * std);
			sum += filter[i * size + j];
		}
	}
	
	//convolution
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			sum = 0;
			for (int k = 0; k < size; k++) {
				for (int l = 0; l < size; l++) {
					sum += src[(i + k) * (width + size - 1) + j + l] * filter[k * size + l];
				}
			}
			dst[i * width + j] = sum;
		}
	}
	delete[] filter;
}


void medianFilter(unsigned char* src, unsigned char* dst, int size){
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			int* arr = new int[size * size];
			for (int k = 0; k < size; k++) {
				for (int l = 0; l < size; l++) {
					arr[k * size + l] = src[(i + k) * (width + size - 1) + j + l];//copy
				}
			}
			sort(arr, arr + size * size);//sort ascending order
			dst[i * width + j] = arr[size * size / 2];//median
		}
	}
}


void alphaTrimmedMeanFilter(unsigned char* src, unsigned char* dst, int size, double a){
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			int* arr = new int[size * size];
			for (int k = 0; k < size; k++) {
				for (int l = 0; l < size; l++) {
					arr[k * size + l] = src[(i + k) * (width + size - 1) + j + l];//copy
				}
			}
			sort(arr, arr + size * size);//sort ascending order
			int sum = 0;
			for (int k = size * size * a; k < size * size - size * size * a; k++)//trim array
				sum += arr[k];
			dst[i * width + j] = sum/(size*size-int(2*floor(size*size*a)));//mean value of trimmed elements
		}
	}
}


void padding(unsigned char* img, unsigned char* pad, int size) {
	//upper copy padding
	for (int i = 0; i < size / 2; i++) {
		int j = 0;
		for (j; j < size / 2; j++) {
			pad[i * (width + size - 1) + j] = img[0];
		}
		for (j; j < size / 2 + width; j++) {
			pad[i * (width + size - 1) + j] = img[j - size / 2];
		}
		for (j; j < size + width - 1; j++) {
			pad[i * (width + size - 1) + j] = img[width - 1];
		}
	}
	//left right copy padding
	for (int i = size / 2; i < size / 2 + height; i++) {
		//left
		for (int j = 0; j < size / 2; j++) {
			pad[i * (width + size - 1) + j] = img[(i - size / 2) * width];
		}
		for (int j = size / 2; j < size / 2 + width; j++) {
			pad[i * (width + size - 1) + j] = img[(i - size / 2) * width + j - size / 2];
		}
		//right
		for (int j = size / 2 + width; j < size + width - 1; j++) {
			pad[i * (width + size - 1) + j] = img[(i - size / 2 + 1) * width - 1];
		}
	}
	//down copy padding
	for (int i = size / 2 + height; i < size + height - 1; i++) {
		int j = 0;
		for (j; j < size / 2; j++) {
			pad[i * (width + size - 1) + j] = img[(height - 1) * width];
		}
		for (j; j < size / 2 + width; j++) {
			pad[i * (width + size - 1) + j] = img[(height - 1) * width + j - size / 2];
		}
		for (j; j < size + width - 1; j++) {
			pad[i * (width + size - 1) + j] = img[height * width - 1];
		}
	}
}