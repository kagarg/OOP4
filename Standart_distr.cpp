#include "Standart_distr.h"

void Std_distr::set_shift(double shift) {
	Std_distr::u = shift;
}
void Std_distr::set_scale(double scale) {
	Std_distr::l = scale;
}

void Std_distr::set_form(double form) {
	Std_distr::v = form;
}

double Std_distr::get_shift() const {
	return Std_distr::u;
}
double Std_distr::get_scale() const {
	return Std_distr::l;
}
double Std_distr::get_form() const {
	return Std_distr::v;
}

double Std_distr::Standart(double x) const {
	double xmod = ((x - Std_distr::u) / Std_distr::l);
	return xmod;
}

double Std_distr::density(double x) const {
	double xmod = Standart(x);
	double mod_density = (Std_distr::v / (2 * tgamma(1 / Std_distr::v)) * exp(-pow(abs(xmod), Std_distr::v))) / Std_distr::l;
	return mod_density;
}

double Std_distr::expected_value() const {
	double mod_expected_value = Std_distr::u;
	return mod_expected_value;
}

double Std_distr::dispersion() const {
	double mod_dispersion = (tgamma(3 / Std_distr::v) / tgamma(1 / Std_distr::v)) * pow(Std_distr::l, 2);
	return mod_dispersion;
}

double Std_distr::excess() const {
	double excess = ((tgamma(5 / Std_distr::v) * tgamma(1 / Std_distr::v) / pow(tgamma(3 / Std_distr::v), 2)) - 3);
	return excess;
}

double Std_distr::asymmetry() const { return 0; }

double Std_distr::Randomizer() const {
	double r;
	do r = (double)rand() / RAND_MAX; while (r == 0 || r == 1);
	return r;
}

// part of an algorithm for v in range [1;2)
double Std_distr::Random_item12() const {
	double r = Randomizer();
	double a = (1 / Std_distr::v) - 1;
	double b = 1 / (pow(Std_distr::v, 1 / Std_distr::v));
	double x = 0;
	if (r <= 0.5) {
		x = b * log(2 * r);
	}
	else if (r > 0.5) {
		x = -b * log(2 * (1 - r));
	}
	double r2 = Randomizer();
	if (log(r2) <= (exp(Std_distr::v * (-log(abs(x)))) + (abs(x) / b) + a)) { return x; }
	else { return  Random_item12(); }
}

// part of an algorithm for v in range [2; inf)
double Std_distr::Random_item2()const {
	double a = (1 / Std_distr::v) - 0.5;
	double b = 1 / (pow(Std_distr::v, 1 / Std_distr::v));
	double c = 2 * pow(b, 2);
	double r = Randomizer();
	double r2 = Randomizer();
	double x = b * sqrt(-2 * log(r)) * cos(2 * M_PI * r2);
	double r3 = Randomizer();
	if (log(r3) <= (exp(Std_distr::v * (-log(abs(x)))) + (pow(x, 2) / c) + a)) { return x; }
	else { return Random_item2(); }
}
// function for the whole algorithm use
double Std_distr::rand_var() const {
	double result = 0;/////////////////////////////////
	if ((1 <= Std_distr::v) && (Std_distr::v < 2)) { result = Random_item12(); }
	else if (Std_distr::v >= 2) { result = Random_item2(); }
	return ((result * Std_distr::l + Std_distr::u));
}

std::vector<double> Std_distr::generate_selection(const int n) const {
	std::vector<double> result;
	for (int i = 0; i < n; i++) {
		result.push_back(rand_var());
	}
	sort(result.begin(), result.end());
	return result;
}

std::vector<std::pair<double, double>> Std_distr::generate_graph_selection(const std::vector<double>& vec, const int n) const {
	std::vector<std::pair<double, double>> result;
	for (int i = 0; i < n; ++i) {
		result.push_back(std::make_pair(vec[i], density(vec[i])));
	}
	return result;
}

void Std_distr::save_in_file(std::ofstream& file) {
	file << v << std::endl << u << std::endl << l << std::endl;
}

void Std_distr::load_from_file(std::ifstream& file) {
	double v, u, l;
	if (!file.is_open()) {
		throw 0;
	}
	file >> v >> u >> l;
	if (v <= 0 || l <= 0) {
		throw 1;
	}
	this->v = v;
	this->u = u;
	this->l = l;
}
