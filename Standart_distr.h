#pragma once
#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include <time.h>
#include <string>
#include <fstream>
#include "Distributions.h"

class Std_distr : public IDistribution, public IPresistend {
public:

	Std_distr() :u(0), l(1), v(1) {}; // конструктор по умолчанию

	Std_distr(double form0, double shift0, double scale0) :v(form0), u(shift0), l(scale0) {}; // Конструктор лего

	Std_distr(std::string File_name) { // конструктор из файла 
		double symbol;
		double mass[3];
		std::ifstream in(File_name); // открываем файл для чтения
		if (in.is_open())
		{
			for (int i = 0; i < 3; i++)
			{
				in >> symbol;
				while (in.fail())
				{
					in.clear();
					in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
					std::cout << "Ошибка в данных из входного файла. Производится выход из программы. " << std::endl;
					exit(0);
				}
				mass[i] = symbol;
			}
		}
		else { std::cout << "Неверное имя файла"; exit(0); }
		in.close();     // закрываем файл
		Std_distr::set_form(mass[2]);
		if (mass[2] >= 1) { Std_distr::set_form(mass[2]); }
		else {
			std::cout << "Ошибка в значении параметра формы. Производится выход из программы. ";
			exit(0);
		}
		Std_distr::set_shift(mass[0]);
		Std_distr::set_scale(mass[1]);
	}

	~Std_distr() // деструктор
	{

	}

	double density(const double x) const override;
	double expected_value() const override;
	double dispersion() const override;
	double excess() const override;
	double asymmetry() const override;
	double rand_var() const override;
	std::vector<double> generate_selection(const int n) const override;
	std::vector<std::pair<double, double>> generate_graph_selection(const std::vector<double>& selection, const int n) const override;

	void set_form(const double n);
	void set_shift(const double mu);
	void set_scale(const double lambda);

	double get_form() const;
	double get_shift() const;
	double get_scale() const;

	void load_from_file(std::ifstream& file) override;
	void save_in_file(std::ofstream& file) override;
private:
	double v, u, l;
	double Standart(double x)const;

	double Randomizer()const;
	double Random_item12()const;
	double Random_item2()const;
};
