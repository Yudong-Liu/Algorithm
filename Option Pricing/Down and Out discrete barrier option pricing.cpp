
//Down and Out discrete barrier option pricing
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <random>
#include <chrono>
#include <tuple>
#include "normdist.h"
using namespace std;
#define pi 3.141592653589793

double risk_free_rate, strike_price, barrier_price, R;
double initial_stock_price, expiration_time, volatility;
int no_of_barriers;
int no_of_trials;

double max(double a, double b) {
	return (b < a) ? a : b;
}

// u.i.i.d. generator
unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
double get_uniform()
{
	std::uniform_real_distribution <double> distribution(0.0, 1.0);
	double number = distribution(generator);
	return (number);
}

double N(const double& z) {
	if (z > 6.0) { return 1.0; }; // this guards against overflow 
	if (z < -6.0) { return 0.0; };
	double b1 = 0.31938153;
	double b2 = -0.356563782;
	double b3 = 1.781477937;
	double b4 = -1.821255978;
	double b5 = 1.330274429;
	double p = 0.2316419;
	double c2 = 0.3989423;
	double a = fabs(z);
	double t = 1.0 / (1.0 + a * p);
	double b = c2 * exp((-z)*(z / 2.0));
	double n = ((((b5*t + b4)*t + b3)*t + b2)*t + b1)*t;
	n = 1.0 - b * n;
	if (z < 0.0) n = 1.0 - n;
	return n;
};

double prob(double stock_price) {	
	double probability= 1.0;
	double mean, variance;
	for (int j = 1; j <= no_of_barriers; j++) {
		mean = initial_stock_price + (((double)j) / ((double)no_of_barriers)*(stock_price - initial_stock_price));

		variance = (((double)j) / ((double)no_of_barriers))*expiration_time*(1.0 - ((double)j) / ((double)no_of_barriers));

		probability *= (1.0 - N((barrier_price - mean) / sqrt(variance)));
	}
	return probability;
}


tuple<double, double, double, double> calculate() {

	double delta_t = expiration_time / ((double)no_of_barriers);
	double delta_R = (risk_free_rate - 0.5*pow(volatility, 2))*delta_t;
	double delta_SD = volatility * sqrt(delta_t);
	double S1 = initial_stock_price;
	double S2 = initial_stock_price;
	double S3 = initial_stock_price;
	double S4 = initial_stock_price;
	double f1 = 1.0, f2 = 1.0, f3 = 1.0, f4 = 1.0;
	double call_option_price = 0.0;
	double put_option_price = 0.0;
	double barriers_call_option_price = 0.0;
	double barriers_put_option_price = 0.0;
	for (int i = 1; i <= no_of_barriers; i++) {
		double x = get_uniform();
		double y = get_uniform();
		double a = sqrt(-2.0*log(x)) * cos(2 * pi*y);
		double b = sqrt(-2.0*log(x)) * sin(2 * pi*y);

		S1 = S1 * exp(delta_R + delta_SD * a);
		S2 = S2 * exp(delta_R - delta_SD * a);
		S3 = S3 * exp(delta_R + delta_SD * b);
		S4 = S4 * exp(delta_R - delta_SD * b);
		
		if (S1 <= barrier_price)	f1 = 0.0;
		if (S2 <= barrier_price)	f2 = 0.0;
		if (S3 <= barrier_price)	f3 = 0.0;
		if (S4 <= barrier_price)	f4 = 0.0;
	}

	call_option_price += f1 * max(0.0, S1 - strike_price);
	call_option_price += f2 * max(0.0, S2 - strike_price);
	call_option_price += f3 * max(0.0, S3 - strike_price);
	call_option_price += f4 * max(0.0, S4 - strike_price);

	put_option_price += f1 * max(0.0, strike_price - S1);
	put_option_price += f2 * max(0.0, strike_price - S2);
	put_option_price += f3 * max(0.0, strike_price - S3);
	put_option_price += f4 * max(0.0, strike_price - S4);

	call_option_price = exp(-risk_free_rate * expiration_time)*call_option_price / 4.0;
	put_option_price = exp(-risk_free_rate * expiration_time)*put_option_price / 4.0;

	barriers_call_option_price += exp(-risk_free_rate * expiration_time)*max(0.0, (S1 - strike_price)) *prob(S1)*f1;
	barriers_call_option_price += exp(-risk_free_rate * expiration_time)*max(0.0, (S2 - strike_price)) *prob(S2)*f2;
	barriers_call_option_price += exp(-risk_free_rate * expiration_time)*max(0.0, (S3 - strike_price)) *prob(S3)*f3;
	barriers_call_option_price += exp(-risk_free_rate * expiration_time)*max(0.0, (S4 - strike_price)) *prob(S4)*f4;
	barriers_call_option_price = barriers_call_option_price / 4.0;

	barriers_put_option_price += exp(-risk_free_rate * expiration_time)*max(0.0, (strike_price - S1)) *prob(S1)*f1;
	barriers_put_option_price += exp(-risk_free_rate * expiration_time)*max(0.0, (strike_price - S2)) *prob(S2)*f2;
	barriers_put_option_price += exp(-risk_free_rate * expiration_time)*max(0.0, (strike_price - S3)) *prob(S3)*f3;
	barriers_put_option_price += exp(-risk_free_rate * expiration_time)*max(0.0, (strike_price - S4)) *prob(S4)*f4;
	barriers_put_option_price = barriers_put_option_price / 4.0;

	return make_tuple(call_option_price, put_option_price, barriers_call_option_price, barriers_put_option_price);
}


int main(int argc, char* argv[])
{
	sscanf_s(argv[1], "%lf", &expiration_time);
	sscanf_s(argv[2], "%lf", &risk_free_rate);
	sscanf_s(argv[3], "%lf", &volatility);
	sscanf_s(argv[4], "%lf", &initial_stock_price);
	sscanf_s(argv[5], "%lf", &strike_price);
	sscanf_s(argv[6], "%d", &no_of_trials);
	sscanf_s(argv[7], "%d", &no_of_barriers);
	sscanf_s(argv[8], "%lf", &barrier_price);

	cout << "European Down and Out Discrete Barrier Option Pricing" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "Barrier Price = " << barrier_price << endl;
	cout << "Number of Trials = " << no_of_trials << endl;
	cout << "Number of Discrete Barriers = " << no_of_barriers << endl;

	double call_option_price = 0.0;
	double put_option_price = 0.0;
	double barriers_call_option_price = 0.0;
	double barriers_put_option_price = 0.0;
	double p1, p2, p3, p4;
	for (int i = 0; i < no_of_trials; i++)
	{
		tie(p1, p2, p3, p4) = calculate(); //return the results of prices
		call_option_price += p1;
		put_option_price += p2;
		barriers_call_option_price += p3;
		barriers_put_option_price += p4;
	}
	call_option_price = call_option_price / (double)no_of_trials;
	put_option_price = put_option_price / (double)no_of_trials;
	barriers_call_option_price = barriers_call_option_price / (double)no_of_trials;
	barriers_put_option_price = barriers_put_option_price / (double)no_of_trials;

	cout << "--------------------------------------" << endl;
	cout << "Price of an European Down and Out Call Option via explicit simulation = " << call_option_price << endl;
	cout << "Price of an European Down and Out Discrete Barrier Call Option = " << barriers_call_option_price << endl;
	cout << "--------------------------------------" << endl;
	cout << "Price of an European Down and Out Put Option via explicit simulation = " << put_option_price << endl;
	cout << "Price of an European Down and Out Discrete Barrier Put Option = " << barriers_put_option_price << endl;
	cout << "--------------------------------------" << endl;

	system("pause");
}
