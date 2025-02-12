#include <stdio.h>
#include <stdint.h>
#include <math.h>

#define PI  3.14159265359
typedef struct {
	double real;
	double imag;
}Complex_t;

typedef struct {
	double gain;
	double phase;
}Polar_t;


#define SPEED_LIGHT  300000000.0 //Speed of light
#define FREQUENCY 1000000000.0 //frequency in HZ
#define WAVE_LENGTH  (SPEED_LIGHT / FREQUENCY)//wave length
#define	WAVE_NUMBER  (2 * PI / WAVE_LENGTH)  // wave number
#define X_ELEMNTS 8
#define Y_ELEMNTS 8
#define d_X_ELEMNT 0.5
#define d_Y_ELEMNT 0.5

#define SAMPLING_RATE 1000000.0   //Hz
#define F_TONE 20000 // Hz
#define N_SAMPLES 1000

#define SIGNAL_DURATION N_SAMPLES / SAMPLING_RATE;
#define NUMBER_OF_ELEMENTS (X_ELEMNTS*Y_ELEMNTS)

Complex_t sv[64];

Complex_t COMPLEX_TONE[N_SAMPLES];
Complex_t transmitted_signal[NUMBER_OF_ELEMENTS][N_SAMPLES];
Complex_t X_weighted[N_SAMPLES];

Complex_t expC(double phase)
{
	Complex_t ret_value;
	ret_value.real = cos(phase);
	ret_value.imag = sin(phase);

	return ret_value;
}

Complex_t  CMul(Complex_t  in1, Complex_t in2)
{
	Complex_t ret_val;
	ret_val.real = in1.real * in2.real - in1.imag * in2.imag;
	ret_val.imag = in1.real * in2.imag + in1.imag * in2.real;
	return ret_val;
}


void RectToPolar(Complex_t in[], int size, Polar_t out[])
{
	double gain, phase;
	double ming = 21543826;

	for (int i = 0; i < size; i++)
	{
		double gainbias, gfainscale;
		gain = in[i].real * in[i].real + in[i].imag * in[i].imag;
		phase = atan2(in[i].imag, in[i].real) * 180 / PI / 5.6;
		out[i].gain = gain;
		out[i].phase = phase;
	}

}

Complex_t  ConjCMul(Complex_t  in1, Complex_t in2)
{
	Complex_t ret_val;
	ret_val.real = in1.real * in2.real + in1.imag * in2.imag;
	ret_val.imag = in1.real * in2.imag - in1.imag * in2.real;
	return ret_val;
}
void PrintComplex(Complex_t inp[], int size)
{
	int i;
	for (i = 0; i < size; i++)
		printf("%lf %lfj\n", inp[i].real, inp[i].imag);
}

void PrintPolar(Polar_t inp[], int size)
{
	int i;
	for (i = 0; i < size; i++)
		printf("Gain :%lf, Phase: %lfj\n", inp[i].gain, inp[i].phase);
}

Complex_t  ComplexAdd(Complex_t in1, Complex_t in2)
{
	Complex_t ret_val;
	ret_val.real = in1.real + in2.real;
	ret_val.imag = in1.imag + in2.imag;
	return ret_val;

}

//----------------------------------------------------------------------

Complex_t SteeringVector[X_ELEMNTS * Y_ELEMNTS];
Complex_t PhaseNGain[X_ELEMNTS * Y_ELEMNTS];
void  steering_vector(const float theta, const float phi, Complex_t  sv[])
{
	int i, j;
	Complex_t isv;
	double theta_rad;
	double phi_rad;
	double xposition[X_ELEMNTS];
	double yposition[Y_ELEMNTS];
	double phase;

	for (i = 0; i < X_ELEMNTS; i++)
	{
		xposition[i] = i * d_X_ELEMNT * WAVE_LENGTH;
		yposition[i] = i * d_Y_ELEMNT * WAVE_LENGTH;
	}
	theta_rad = theta * (PI / 180.0);
	phi_rad = phi * (PI / 180.0);

	for (i = 0; i < Y_ELEMNTS; i++)
	{
		for (j = 0; j < X_ELEMNTS; j++)
		{
			phase = WAVE_NUMBER * (sin(theta_rad) * cos(phi_rad) * xposition[j] + sin(theta_rad) * sin(phi_rad) * yposition[i]);
			sv[i * X_ELEMNTS + j] = expC(phase);
			if ((i == 1) && (j == 7))
			{
				sv[i * X_ELEMNTS + j].real = sv[i * X_ELEMNTS + j].real * 0.8;
				sv[i * X_ELEMNTS + j].imag = sv[i * X_ELEMNTS + j].imag * 0.8;
			}
		}
	}
}
//----------------------------------------------------------------------


int main(void)
{
	int i, j;
	double azimuth, elevation;
	Complex_t sum;
	double maxVal = 0;
	double index_azi, index_ele;

	for (i = 0; i < N_SAMPLES; i++)
	{
		double frequency = 2.0 * PI * F_TONE / SAMPLING_RATE;
		static double phase = 0;

		COMPLEX_TONE[i] = expC(phase);

		phase += frequency;

		if (phase > PI)
		{
			phase -= 2.0 * PI;
		}
		else if (phase < -PI)
		{
			phase += 2.0 * PI;
		}
	}

	steering_vector(40.0, 50.0, SteeringVector);


	for (i = 0; i < NUMBER_OF_ELEMENTS; i++)
	{
		for (j = 0; j < N_SAMPLES; j++)
		{
			transmitted_signal[i][j] = CMul(COMPLEX_TONE[j], SteeringVector[i]);
		}
	}

	for (azimuth = 0; azimuth < 360; azimuth += 3.6)
	{
		for (elevation = 0; elevation < 90; elevation += 0.9)
		{
			double magnitude_square;
			printf("azimuth - %lf, elevation - %lf\n", azimuth, elevation);
			steering_vector(elevation, azimuth, SteeringVector);		

			/*
			* transmitted_signal	:	64x10000	input
			* SteeringVector		:	64x1		steering vector
			* X_weighted			:	1x10000		output
			*/	
			magnitude_square = 0;
			for (j = 0; j < N_SAMPLES; j++)
			{
				sum.real = 0;
				sum.imag = 0;
				for (i = 0; i < NUMBER_OF_ELEMENTS; i++)
				{
					sum = ComplexAdd(sum, ConjCMul(SteeringVector[i], transmitted_signal[i][j]));
				}
				X_weighted[j] = sum;
				magnitude_square += sum.real * sum.real + sum.imag * sum.imag;
			}
			if (maxVal < magnitude_square)
			{
				maxVal = magnitude_square;
				index_azi = azimuth;
				index_ele = elevation;
			}	
		}
	}
	printf("azimuth : %lf, elevation: %lf\n", index_azi, index_ele);
	steering_vector(index_ele, index_azi, SteeringVector);
	RectToPolar(SteeringVector, 64, PhaseNGain);
	PrintPolar(PhaseNGain, 64);
	return 0;
}