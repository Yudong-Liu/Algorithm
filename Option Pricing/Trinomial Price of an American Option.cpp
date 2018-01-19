//Trinomial Price of an American Option
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
using namespace std;

float up_factor,down_factor, uptick_prob, downtick_prob, risk_free_rate, strike_price;
float initial_stock_price, expiration_time, volatility, R;
int no_of_divisions;

float max(float a, float b) {
	return (b < a) ? a : b;
}

float american_call_option(int k, int i, float current_stock_price) {
	float up_price, down_price, stable_price;
	if (k == no_of_divisions)
		return max(0.0, (current_stock_price - strike_price));
	else{
		up_price = uptick_prob*american_call_option(k + 1, i + 1, current_stock_price*up_factor);
		stable_price = (1 - uptick_prob - downtick_prob)*american_call_option(k + 1, i, current_stock_price);
		down_price = downtick_prob*american_call_option(k + 1, i - 1, current_stock_price * down_factor);

		return max((current_stock_price - strike_price), (up_price + down_price + stable_price) / R);
	}		
}

float american_put_option(int k, int i, float current_stock_price) {
	float up_price, down_price, stable_price;
	if (k == no_of_divisions)
		return max(0.0, (strike_price - current_stock_price));
	else {
		up_price = uptick_prob*american_put_option(k + 1, i + 1, current_stock_price*up_factor);
		stable_price = (1 - uptick_prob - downtick_prob)*american_put_option(k + 1, i, current_stock_price);
		down_price = downtick_prob*american_put_option(k + 1, i - 1, current_stock_price * down_factor);
		return max((strike_price - current_stock_price), (up_price + down_price + stable_price) / R);
	}
		
}

int main(int argc, char* argv[])
{

	sscanf_s(argv[1], "%f", &expiration_time);
	sscanf_s(argv[2], "%d", &no_of_divisions);
	sscanf_s(argv[3], "%f", &risk_free_rate);
	sscanf_s(argv[4], "%f", &volatility);
	sscanf_s(argv[5], "%f", &initial_stock_price);
	sscanf_s(argv[6], "%f", &strike_price);

	up_factor = exp(volatility*sqrt(2 * (expiration_time / ((float)no_of_divisions))));;
	down_factor = (float)1.0 / up_factor;
	R = exp(risk_free_rate*expiration_time / ((float)no_of_divisions));
	uptick_prob = pow((sqrt(R) - sqrt(down_factor)) / (sqrt(up_factor) - sqrt(down_factor)),2);
	downtick_prob = pow((sqrt(up_factor) - sqrt(R)) / (sqrt(up_factor) - sqrt(down_factor)), 2);
	
	cout << "Recursive Trinomial American-Asian Option Pricing" << endl;
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
	cout << "Trinomial Price of an American Call Option = " << american_call_option(0, 0, initial_stock_price) << endl;
	cout << "Trinomial Price of an American Put Option = " << american_put_option(0, 0, initial_stock_price) << endl;
	
	system("pause");
}

