//Down and Out Barrier option pricing using Monte Carlo 
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
int no_of_divisions;
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

double option_price_put_black_scholes(const double& S,      // spot price
	const double& K,      // Strike (exercise) price,
	const double& r,      // interest rate
	const double& sigma,  // volatility
	const double& time)
{
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r * time) / (sigma*time_sqrt) + 0.5*sigma*time_sqrt;
	double d2 = d1 - (sigma*time_sqrt);
	return K * exp(-r * time)*N(-d2) - S * N(-d1);
};

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
	const double& K,       // strike (exercise) price,
	const double& r,       // interest rate
	const double& sigma,   // volatility 
	const double& time)	  // time to maturity 
{
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r * time) / (sigma*time_sqrt) + 0.5*sigma*time_sqrt;
	double d2 = d1 - (sigma*time_sqrt);
	return S * N(d1) - K * exp(-r * time)*N(d2);
};

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

double closed_form_down_and_out_european_call_option()
{
	double K = (2 * risk_free_rate) / (volatility*volatility);
	double A = option_price_call_black_scholes(initial_stock_price, strike_price,risk_free_rate, volatility, expiration_time);
	double B = (barrier_price*barrier_price) / initial_stock_price;
	double C = pow(initial_stock_price / barrier_price, -(K - 1));
	double D = option_price_call_black_scholes(B, strike_price, risk_free_rate, volatility, expiration_time);
	return (A - D * C);
}

double closed_form_down_and_in_european_put_option()
{
	double S = initial_stock_price;
	double r = risk_free_rate;
	double T = expiration_time;
	double sigma = volatility;
	double H = barrier_price;
	double X = strike_price;

	double lambda = (r + ((sigma*sigma) / 2)) / (sigma*sigma);
	double temp = 2 * lambda - 2.0;
	double x1 = (log(S / H) / (sigma*sqrt(T))) + (lambda*sigma*sqrt(T));
	double y = (log(H*H / (S*X)) / (sigma*sqrt(T))) + (lambda*sigma*sqrt(T));
	double y1 = (log(H / S) / (sigma*sqrt(T))) + (lambda*sigma*sqrt(T));
	return (-S * N(-x1) + X * exp(-r * T)*N(-x1 + sigma * sqrt(T)) +
			S * pow(H / S, 2 * lambda)*(N(y) - N(y1)) -
			X * exp(-r * T)*pow(H / S, temp)*(N(y - sigma * sqrt(T)) - N(y1 - sigma * sqrt(T))));
}

double closed_form_down_and_out_european_put_option()
{
	double vanilla_put = option_price_put_black_scholes(initial_stock_price, strike_price,risk_free_rate, volatility, expiration_time);
	double put_down_in = closed_form_down_and_in_european_put_option();
	return (vanilla_put - put_down_in);
}

double adjusted_prob(double stock_price) {
	double probability_stock_path_has_hit_barrier_by_the_time_it_got_here = 
		exp((-2.0 / (volatility*volatility*expiration_time))*log(initial_stock_price / barrier_price)*log(stock_price / barrier_price));
	if (initial_stock_price > barrier_price && stock_price > barrier_price)
		return (1.0 - probability_stock_path_has_hit_barrier_by_the_time_it_got_here);
	else
		return 0.0;
}

tuple<double, double, double, double> Monte_Carlo() {
	double delta_t = expiration_time / ((double)no_of_divisions);
	double delta_R = (risk_free_rate - 0.5*pow(volatility, 2))*delta_t;
	double delta_SD = volatility * sqrt(delta_t);
	double S1 = initial_stock_price;
	double S2 = initial_stock_price;
	double S3 = initial_stock_price;
	double S4 = initial_stock_price;
	double f1 = 1, f2 = 1, f3 = 1, f4 = 1;
	for (int i = 0; i < no_of_divisions; i++) {
		double x = get_uniform();
		double y = get_uniform();
		double a = sqrt(-2.0*log(x)) * cos(2 * pi*y);
		double b = sqrt(-2.0*log(x)) * sin(2 * pi*y);

		S1 = S1 * exp(delta_R + delta_SD * a);
		S2 = S2 * exp(delta_R - delta_SD * a);
		S3 = S3 * exp(delta_R + delta_SD * b);
		S4 = S4 * exp(delta_R - delta_SD * b);

		if (S1 <= barrier_price)	f1 = 0;
		if (S2 <= barrier_price)	f2 = 0;
		if (S3 <= barrier_price)	f3 = 0;
		if (S4 <= barrier_price)	f4 = 0;
					
	}
	double call_option_price = 0.0;
	double put_option_price = 0.0;
	double adjusted_call_option_price = 0.0;
	double adjusted_put_option_price = 0.0;
	
	call_option_price += f1 * max(0.0, S1 - strike_price);
	call_option_price += f2 * max(0.0, S2 - strike_price);
	call_option_price += f3 * max(0.0, S3 - strike_price);
	call_option_price += f4 * max(0.0, S4 - strike_price);

	put_option_price += f1 * max(0.0, strike_price - S1);
	put_option_price += f2 * max(0.0, strike_price - S2);
	put_option_price += f3 * max(0.0, strike_price - S3);
	put_option_price += f4 * max(0.0, strike_price - S4);

	adjusted_call_option_price +=  max(0.0, S1 - strike_price)*adjusted_prob(S1);
	adjusted_call_option_price +=  max(0.0, S2 - strike_price)*adjusted_prob(S2);
	adjusted_call_option_price +=  max(0.0, S3 - strike_price)*adjusted_prob(S3);
	adjusted_call_option_price +=  max(0.0, S4 - strike_price)*adjusted_prob(S4);

	adjusted_put_option_price +=  max(0.0, strike_price - S1)*adjusted_prob(S1);
	adjusted_put_option_price +=  max(0.0, strike_price - S2)*adjusted_prob(S2);
	adjusted_put_option_price +=  max(0.0, strike_price - S3)*adjusted_prob(S3);
	adjusted_put_option_price +=  max(0.0, strike_price - S4)*adjusted_prob(S4);

	call_option_price = exp(-risk_free_rate * expiration_time)*call_option_price / 4.0;
	put_option_price = exp(-risk_free_rate * expiration_time)*put_option_price / 4.0;
	adjusted_call_option_price= exp(-risk_free_rate * expiration_time)*adjusted_call_option_price / 4.0;
	adjusted_put_option_price = exp(-risk_free_rate * expiration_time)*adjusted_put_option_price / 4.0;
	return make_tuple(call_option_price, put_option_price, adjusted_call_option_price, adjusted_put_option_price) ;
}


int main(int argc, char* argv[])
{
	sscanf_s(argv[1], "%lf", &expiration_time);	
	sscanf_s(argv[2], "%lf", &risk_free_rate);
	sscanf_s(argv[3], "%lf", &volatility);
	sscanf_s(argv[4], "%lf", &initial_stock_price);
	sscanf_s(argv[5], "%lf", &strike_price);
	sscanf_s(argv[6], "%d", &no_of_trials);
	sscanf_s(argv[7], "%d", &no_of_divisions);
	sscanf_s(argv[8], "%lf", &barrier_price);
	
	cout << "European Down-and-Out Barriers Options Pricing" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "Barrier Price = " << barrier_price << endl;
	cout << "Number of Trials = " << no_of_trials << endl;
	cout << "Number of Divisions = " << no_of_divisions << endl;
	cout << "--------------------------------------" << endl;

	double call_option_price = 0.0;
	double put_option_price = 0.0;
	double adjusted_call_option_price = 0.0;
	double adjusted_put_option_price = 0.0;
	double p1, p2, p3, p4;
	for (int i = 0; i < no_of_trials; i++)
	{
		tie(p1, p2, p3, p4) = Monte_Carlo(); //return the results of prices
		call_option_price += p1;
		put_option_price += p2;
		adjusted_call_option_price += p3;
		adjusted_put_option_price += p4;
	}
	double average_call_option_price = call_option_price / (double)no_of_trials;
	double average_put_option_price = put_option_price / (double)no_of_trials;
	double ad_call_option_price = adjusted_call_option_price / (double)no_of_trials;
	double ad_put_option_price = adjusted_put_option_price / (double)no_of_trials;

	cout << "--------------------------------------" << endl;
	cout << "Price of an European Down and Out Call Option = " << average_call_option_price << endl;
	cout << "The Call Price using the (1-p)-adjustment term = " << ad_call_option_price << endl;
	cout << "Price of an European Down and Out Call Option from Theory = " <<closed_form_down_and_out_european_call_option() << endl;
	cout << "--------------------------------------" << endl;
	cout << "Price of an European Down and Out Put Option = " << average_put_option_price << endl;
	cout << "The Put Price using the (1-p)-adjustment term = " << ad_put_option_price << endl;
	cout << "Price of an European Down and Out Put Option from Theory = " <<closed_form_down_and_out_european_put_option() << endl;
	cout << "--------------------------------------" << endl;

}
