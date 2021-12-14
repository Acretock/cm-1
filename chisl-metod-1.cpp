#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>

using namespace std;

const double PI = 3.141592653589793238;

double function(double x) {
    return x * sin(x) * cos(2.0 * x);
}

double Lagrange(std::vector <double>& X, std::vector <double>& Y, double xp, double n) {
	double res = 0;

	for (int i = 0; i <= n; i++) {
		double tempProd = 1;
		for (int j = 0; j <= n; j++) {
			if (i != j) {
				tempProd *= (xp - X[j]) / (X[i] - X[j]);
			}
		}
		tempProd *= Y[i];
		res += tempProd;
	}
	return res;
}

double chebyshevNode(double i, double n0) {
	auto a = 2.0;
	auto b = 0.0;
	return ((a + b) / 2.0 + (a - b) / 2.0 * cos((PI * (2.0 * i + 1.0) / (2.0 * n0))));
}

double newtonianInterpolation(std::vector<double>& coeffs, std::vector<double>& xValues, double x) {
	double result = 0;
	int n = coeffs.size();
	result += coeffs[0];
	for (int i = 1; i < n; i++) {
		double product = 1;
		for (int k = 0; k < i; k++) {
			product *= (x - xValues[k]);
		}
		result += product * coeffs[i];
	}
	return result;
}

void task_1() {
	std::ofstream out("task_1_output.txt");
	std::ofstream out10("task_1_y.txt");
	double minimumDelta = 1e-7;
double power = 15;

	auto n0 = -1;
	auto start = 2.0;
	auto end = 0.0;
	for (double power = 1; power <= 15; power++) {
	std::vector<double> xValues;
	std::vector<double> yValues;

	double Delta = 0;

	//Setting up nodes grid
	double h = (start - end) / power;
	xValues.push_back(0.0);
	yValues.push_back(function(0));
	for (int i = 1; i <= power; i++) {
		xValues.push_back(xValues[i - 1] + h);
		yValues.push_back(function(xValues[i]));
	}

	//Lagrange interpolation function
	std::vector<double> result;
	for (int i = 0; i <= power; i++) {
		result.push_back(Lagrange(xValues, yValues, xValues[i], power));
	}

	/*if (power == 15) {
		for (int i = 0; i <= power; i++) {
			std::cout << std::fixed << std::setprecision(8) << xValues[i] << "\t\t" << yValues[i] << "\t\t" << result[i] << std::endl;
		}
	}*/

	//Error margin 
	double N = 100'000;

	double step = (start - end) / N;
	double x0 = 0;
	for (int i = 0; i < N; i++) {
		double a = function(x0);
		double b = Lagrange(xValues, yValues, x0, power);
		out10 << b << std::endl;
		Delta = max(abs(a - b), Delta);
		if (power == 15)
			out << std::fixed << std::setprecision(20) << b << std::endl;
		x0 += step;
	}

	std::cout << Delta << std::endl;
	if (Delta < minimumDelta) {
		minimumDelta = Delta;
		n0 = power;
	}
	}
	std::cout << "Minimal delta when n0 = " << n0 << std::endl;
	std::cin.clear();
	out.close();
	out10.close();
}

void task_2() {
	int n0 = 15;
	std::ofstream out2("task_2_cheb.txt");
	
	//Chebyshev nodes
	std::vector<double> xValues;
	std::vector<double> yValues;

	//Uniform nodes 
	std::vector<double> xValuesR;
	std::vector<double> yValuesR;


	std::vector<double> result;
	double s = 1.0 / n0;
	double xs = 0;
	for (double i = 0; i < n0; i++) {
		xValuesR.push_back(xs);
		yValuesR.push_back(function(xValuesR[i]));
		xValues.push_back(chebyshevNode(i, n0));
		yValues.push_back(function(xValues[i]));
		xs += s;
	}
	// chebyshevNode has reversed order
	std::reverse(xValues.begin(), xValues.end());
	std::reverse(yValues.begin(), yValues.end());
	for (int i = 0; i < n0; i++) {
		result.push_back(Lagrange(xValues, yValues, xValues[i], n0 - 1));
	}
	for (int i = 0; i < n0; i++) {
		std::cout << xValues[i] << "\t" << yValues[i]  << "\t" << result[i] << std::endl;
	}

	double N = 100'000;
	double step = (2.0 - 0.0) / N;
	double x = 0, Delta = 0, Delta1 = 0;
	for (int i = 0; i < N; i++) {
		double a = function(x);
		double b = Lagrange(xValues, yValues, x, n0 - 1);
		double c = Lagrange(xValuesR, yValuesR, x, n0 - 1);
		Delta = max(abs(a - b), Delta);
		Delta1 = max(abs(a - c), Delta1);
		out2 << c << std::endl;
		//out2 << a << "\t\t" << b << "\t\t" << abs(a - b) << std::endl;
		x += step;
	}

	std::cout << "Cheibyshev error margin: " << Delta1 << std::endl;
	std::cout << "Uniform error margin: " << Delta << std::endl;
	out2.close();
	std::cin.clear();
}
void task_3() {
	int n0 = 15;
	double task2error = 3.74638e-10;
	std::ofstream out3("task_3_Newtonian.txt");

	std::vector<double> separatedDifferences;
	std::vector<double> xValues;
	std::vector<double> yValues;
	std::vector<double> result;
	for (double i = 0; i < n0; i++) {
		xValues.push_back(chebyshevNode(i, n0));
		yValues.push_back(function(xValues[i]));
	}

	std::reverse(xValues.begin(), xValues.end());
	std::reverse(yValues.begin(), yValues.end());
	separatedDifferences.push_back(yValues[0]);
	for (int n = 1; n < n0; n++) {
		double summ = 0;
		for (int j = 0; j <= n; j++) {
			double prodact = 1;
			for (int i = 0; i <= n; i++) {
				if (i == j) continue;
				prodact *= (xValues[j] - xValues[i]);
			}
			summ += yValues[j] / prodact;
		}
		separatedDifferences.push_back(summ);
	}
	double N = 100000;
	double step = (2.0 - 0.0) / N;
	double x = 0, Delta = 0;
	for (int i = 0; i < N; i++) {
		double x = chebyshevNode(i, N);
		double b = function(x);
		double a = newtonianInterpolation(separatedDifferences, xValues, x);
		Delta = max(abs(b - a), Delta);
		out3 << b << std::endl;
		//out2 << a << "\t\t" << b << "\t\t" << abs(a - b) << std::endl;
		x += step;
	}

	std::cout << "Error estimate for Newton's polynomial : \t\t" << Delta << std::endl;
	std::cout << "Error estimate for the Lagrange polynomial : \t\t" << task2error << std::endl;
	std::cin.clear();
	out3.close();
}

#pragma endregion
#pragma region TASK6

double secondDerrivative(double x) {
	return (-4 * sin(x) - 4 * x * cos(x)) * sin(2 * x) + (2 * cos(x) - 5 * x * sin(x)) * cos(2 * x);
}


std::vector<double> TridiagonalSolve(int n, std::vector<double>& a, std::vector<double>& c, std::vector<double>& b, std::vector<double>& f)
{
	std::vector<double> x(n);
	double m;
	for (int i = 1; i < n; i++)
	{
		m = a[i] / c[i - 1];
		c[i] = c[i] - m * b[i - 1];
		f[i] = f[i] - m * f[i - 1];
	}

	x[n - 1] = f[n - 1] / c[n - 1];

	for (int i = n - 2; i >= 0; i--)
		x[i] = (f[i] - b[i] * x[i + 1]) / c[i];
	return x;

}
double Spline(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d, double x, double xi, double i) {
	return(a[i] + b[i] * (x - xi) + c[i] * pow(x - xi, 2) + d[i] * pow(x - xi, 3));
}
void task_6() {

	std::ofstream out6("task_6_out.txt");
	for (int N = 10; N < 1000; N++) {
		std::vector<double> xValues;
		std::vector<double> yValues;
		std::vector<double> result;
		for (double i = 0; i < N; i++) {
			xValues.push_back(chebyshevNode(i, N));
			yValues.push_back(function(xValues[i]));
		}

		std::reverse(xValues.begin(), xValues.end());
		std::reverse(yValues.begin(), yValues.end());

		std::vector<std::vector<double>> matrix(N);
		std::vector<double> diag(N);
		std::vector<double> lower(N - 1);
		std::vector<double> upper(N - 1);
		//свеободные члены
		std::vector<double> r(N);
		//вектор ответа
		std::vector<double> c(N);


		std::vector<double> b(N);
		for (int i = 0; i < N; i++) {
			matrix[i].resize(N);
			for (int k = 0; k < N; k++) {
				matrix[i][k] = 0;
			}
		}
		diag[0] = 1; diag[N - 1] = 1;
		lower[0] = 0; lower[N - 2] = 0;
		upper[0] = 0; upper[N - 2] = 0;
		for (int i = 1; i < N - 1; i++) {
			diag[i] = 2 * (xValues[i + 1] - xValues[i - 1]);
			if (i < N - 2) {
				lower[i] = xValues[i] - xValues[i - 1];
				upper[i] = xValues[i + 1] - xValues[i];
				r[i] = 3 * ((yValues[i + 1] - yValues[i]) / xValues[i + 1] - (yValues[i] - yValues[i - 1]) / xValues[i]);
			}
		}

		c = TridiagonalSolve(N - 1, lower, diag, upper, r);

		/*for (auto v : c) {
			std::cout << std::fixed << std::setprecision(10) << v << std::endl;
		}*/
		std::vector<double> a(N, 0);
		std::vector<double> d(N, 0);
		for (int i = 1; i < N - 1; i++) {
			a[i] = yValues[i];
			d[i] = (c[i] - c[i - 1]) / (3 * (xValues[i] - xValues[i - 1]));
			b[i] = (a[i] - a[i - 1]) / (xValues[i] - xValues[i - 1]) + (2 * c[i] + c[i - 1]) / 3 * (xValues[i] - xValues[i - 1]);
		}
		double Count = 10000;
		result.resize(0);

		double delta = -10000;

		double step = abs(xValues[0] - xValues[N - 1]) / Count;
		for (int i = 1; i < N - 1; i++) {

			for (double start = xValues[i - 1]; start < xValues[i]; start += step) {
				double temp1 = function(start);
				double temp2 = Spline(a, b, c, d, start, xValues[i], i);
				double error = abs(temp1 - temp2);
				delta = max(delta, error);
				//std::cout << error << std::endl;
				result.push_back(error);

			}

		}
		if (delta < 1e-5) {
			std::cout << "Required number of grid nodes : " << N << std::endl;
			std::reverse(result.begin(), result.end());
			for (auto n : result)
				out6 << std::fixed << std::setprecision(10) << n << std::endl;
			break;
		}
	}

	out6.close();



}
#pragma endregion

int main()
{
		int c = -1;
		while (c != 0) {
			cout << "enter number for task, 0 for exit" << endl;
			cin >> c;
			switch (c) {
			case 1:
				task_1();
				break;
			case 2:
				task_2();
				break;
			case 3:
				task_3();
				break;
			case 6:
				task_6();
				break;
			default:
				break;
			}
		}
}
