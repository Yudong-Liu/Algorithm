//American Option Pricing with memoization
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <map>
#include <vector>
using namespace std;

double up_factor, down_factor, uptick_prob, downtick_prob, risk_free_rate, strike_price;
double initial_stock_price, expiration_time, volatility, R;
int no_of_divisions;
vector <vector <double> > american_call_memo;
vector <vector <double> > american_put_memo;
vector <vector <double> > european_call_memo;
vector <vector <double> > european_put_memo;


double max(double a, double b) {
	return (b < a) ? a : b;
}

double american_call_option(int k, int i, double current_stock_price) {	
	double up_price, down_price, stable_price;
	if (american_call_memo[k][no_of_divisions + i] > -1)
		return american_call_memo[k][no_of_divisions + i];

	if (k == no_of_divisions) {
		american_call_memo[k][no_of_divisions + i] = max(0.0, (current_stock_price - strike_price));
		return american_call_memo[k][no_of_divisions + i];
	}
	else {		
		up_price = uptick_prob*american_call_option(k + 1, i + 1, current_stock_price*up_factor);
		stable_price = (1 - uptick_prob - downtick_prob)*american_call_option(k + 1, i, current_stock_price);
		down_price = downtick_prob*american_call_option(k + 1, i - 1, current_stock_price * down_factor);
		american_call_memo[k][no_of_divisions + i] = max((current_stock_price - strike_price), (up_price + down_price + stable_price) / R);
		return 	american_call_memo[k][no_of_divisions + i];
	}
}

double american_put_option(int k, int i, double current_stock_price) {

	double up_price, down_price, stable_price;
	if (american_put_memo[k][no_of_divisions + i] > -1)
		return american_put_memo[k][no_of_divisions + i];
	if (k == no_of_divisions)
		return max(0.0, (strike_price - current_stock_price));
	else {
		up_price = uptick_prob*american_put_option(k + 1, i + 1, current_stock_price*up_factor);
		stable_price = (1 - uptick_prob - downtick_prob)*american_put_option(k + 1, i, current_stock_price);
		down_price = downtick_prob*american_put_option(k + 1, i - 1, current_stock_price * down_factor);
		american_put_memo[k][no_of_divisions + i]=max((strike_price - current_stock_price), (up_price + down_price + stable_price) / R);
		return american_put_memo[k][no_of_divisions + i];
	}

}

double european_call_option(int k, int i) {
	double up_price, down_price, stable_price;
	if (european_call_memo[k][no_of_divisions + i] > -1)
		return european_call_memo[k][no_of_divisions + i];

	if (k == no_of_divisions) {
		european_call_memo[k][no_of_divisions + i] = max(0.0, (initial_stock_price*pow(up_factor, ((double)i))) - strike_price);
		return european_call_memo[k][no_of_divisions + i];
	}
	else {
		up_price = uptick_prob*european_call_option(k + 1, i + 1);
		stable_price = (1 - uptick_prob - downtick_prob)*european_call_option(k + 1, i);
		down_price = downtick_prob*european_call_option(k + 1, i - 1);
		european_call_memo[k][no_of_divisions + i] = (up_price + down_price + stable_price) / R;
		return 	european_call_memo[k][no_of_divisions + i];
	}
}

double european_put_option(int k, int i) {
	double up_price, down_price, stable_price;
	if (european_put_memo[k][no_of_divisions + i] > -1)
		return european_put_memo[k][no_of_divisions + i];

	if (k == no_of_divisions) {
		european_put_memo[k][no_of_divisions + i] = max(0.0, strike_price - initial_stock_price*pow(up_factor, ((double)i)));
		return european_put_memo[k][no_of_divisions + i];
	}
	else {
		up_price = uptick_prob*european_put_option(k + 1, i + 1);
		stable_price = (1 - uptick_prob - downtick_prob)*european_put_option(k + 1, i);
		down_price = downtick_prob*european_put_option(k + 1, i - 1);
		european_put_memo[k][no_of_divisions + i] = (up_price + down_price + stable_price) / R;
		return 	european_put_memo[k][no_of_divisions + i];
	}
}

double N(const double& z) {
	if (z > 6.0) { return 1.0; }; 
	if (z < -6.0) { return 0.0; };
	double b1 = 0.31938153;
	double b2 = -0.356563782;
	double b3 = 1.781477937;
	double b4 = -1.821255978;
	double b5 = 1.330274429;
	double p = 0.2316419;
	double c2 = 0.3989423;
	double a = fabs(z);
	double t = 1.0 / (1.0 + a*p);
	double b = c2*exp((-z)*(z / 2.0));
	double n = ((((b5*t + b4)*t + b3)*t + b2)*t + b1)*t;
	n = 1.0 - b*n;
	if (z < 0.0) n = 1.0 - n;
	return n;
};

double option_price_put_black_scholes(const double& S,      // spot price
	const double& K,      // Strike (exercise) price,
	const double& r,      // interest rate
	const double& sigma,  // volatility
	const double& time) {
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r*time) / (sigma*time_sqrt) + 0.5*sigma*time_sqrt;
	double d2 = d1 - (sigma*time_sqrt);
	return K*exp(-r*time)*N(-d2) - S*N(-d1);
};

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
	const double& K,       // strike (exercise) price,
	const double& r,       // interest rate
	const double& sigma,   // volatility 
	const double& time) 
{  // time to maturity 
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r*time) / (sigma*time_sqrt) + 0.5*sigma*time_sqrt;
	double d2 = d1 - (sigma*time_sqrt);
	return S*N(d1) - K*exp(-r*time)*N(d2);
};




int main(int argc, char* argv[])
{
	sscanf_s(argv[1], "%lf", &expiration_time);
	sscanf_s(argv[2], "%d", &no_of_divisions);
	sscanf_s(argv[3], "%lf", &risk_free_rate);
	sscanf_s(argv[4], "%lf", &volatility);
	sscanf_s(argv[5], "%lf", &initial_stock_price);
	sscanf_s(argv[6], "%lf", &strike_price);

	up_factor = exp(volatility*sqrt(2 * (expiration_time / ((double)no_of_divisions))));;
	down_factor = (double)1.0 / up_factor;
	R = exp(risk_free_rate*expiration_time / ((double)no_of_divisions));
	uptick_prob = pow((sqrt(R) - sqrt(down_factor)) / (sqrt(up_factor) - sqrt(down_factor)), 2);
	downtick_prob = pow((sqrt(up_factor) - sqrt(R)) / (sqrt(up_factor) - sqrt(down_factor)), 2);
	double american_call_price, american_put_price, european_call_price, european_put_price;

	vector <double> temp(2 * no_of_divisions + 1, -1);
	for (int i = 0; i <= no_of_divisions; i++)
	{
		american_call_memo.push_back(temp);
		american_put_memo.push_back(temp);
		european_call_memo.push_back(temp);
		european_put_memo.push_back(temp);
	}
	american_call_price = american_call_option(0, 0, initial_stock_price);
	american_put_price = american_put_option(0, 0, initial_stock_price);
	european_call_price = european_call_option(0, 0);
	european_put_price = european_put_option(0, 0);


	cout << "Recursive Trinomial American-Asian Option Pricing with Memoization" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Number of Divisions = " << no_of_divisions << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "--------------------------------------" << endl;
	cout << "R = " << R << endl;
	cout << "Up Factor = " << up_factor << endl;
	cout << "Uptick Probability = " << uptick_prob << endl;
	cout << "Downtick Probability = " << downtick_prob << endl;
	cout << "Notick Probability = " << 1 - uptick_prob - downtick_prob << endl;
	cout << "--------------------------------------" << endl;
	cout << "Trinomial Price of an American Call Option = " << american_call_price << endl;
	cout << "Trinomial Price of an American Put Option = " << american_put_price << endl;
	cout << "--------------------------------------" << endl;
	cout << "Trinomial Price of an Europeam Call Option = " << european_call_price << endl;
	cout << "Call Price according to Black-Scholes = " <<
		option_price_call_black_scholes(initial_stock_price, strike_price, risk_free_rate,
			volatility, expiration_time) << endl;
	cout << "--------------------------------------" << endl;
	cout << "Trinomial Price of an Europeam Put Option = " << european_put_price << endl;
	cout << "Put Price according to Black-Scholes = " <<
		option_price_put_black_scholes(initial_stock_price, strike_price, risk_free_rate,
			volatility, expiration_time) << endl;
	cout << "--------------------------------------" << endl;
	

}

