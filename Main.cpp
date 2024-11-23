#define CATCH_CONFIG_RUNNER
#include <iostream>
#include <fstream>
#include "Standart_distr.h"
#include "Empiric.h"
#include "Mixed_distr.cpp"
#include "catch.hpp"
using namespace std;

int run_unit_tests(int argc, char** argv) {
	int result = Catch::Session().run(argc, argv);
	return result;
}

int main(int argc, char** argv) {
	setlocale(LC_ALL, "ru");
	srand((unsigned)time(0));
	ofstream file1("sdt_distr.txt");
	ofstream file2("emp_distr.txt");

	try {

		Std_distr ld1(1, -6, 2);
		Std_distr ld2(1, -2, 1);
		Std_distr ld3(1, 2, 1);
		Std_distr ld4(1, 6, 2);
		MixtureDistribution<Std_distr, Std_distr> md1(ld1, ld2, 0.5);
		MixtureDistribution<Std_distr, Std_distr> md2(ld3, ld4, 0.5);
		MixtureDistribution<MixtureDistribution<Std_distr, Std_distr>, MixtureDistribution<Std_distr, Std_distr>> md(md2, md1, 0.5);
		auto selection = md.generate_selection(1000);
		auto graph1 = md.generate_graph_selection(selection, 1000);
		EmpiricalDistribution ed(selection);
		auto graph2 = ed.generate_graph_selection(selection, 1000);
		EmpiricalDistribution ed2(ed, 1000, 1);
		auto graph3 = ed2.generate_graph_selection(selection, 1000);

		cout << "параметры стандартных распределений:" << std::endl << "n1 = " << md.component1().component1().get_form() << ", mu1 = " << md.component1().component1().get_shift() << ", lambda1 = " << md.component1().component1().get_scale() <<
			"\nn2 = " << md.component1().component2().get_form() << ", mu2 = " << md.component1().component2().get_shift() << ", lambda2 = " << md.component1().component2().get_scale() <<
			"\nn3 = " << md.component2().component1().get_form() << ", mu3 = " << md.component2().component1().get_shift() << ", lambda3 = " << md.component2().component1().get_scale() <<
			"\nn4 = " << md.component2().component2().get_form() << ", mu4 = " << md.component2().component2().get_shift() << ", lambda4 = " << md.component2().component2().get_scale() << ", p = " << md.get_p() << endl << endl;;
		cout << "смесь двух смесий:\nмат. ожидание:" << md.expected_value() << "\nдисперсия: " << md.dispersion() << "\nассиметрия: " << md.asymmetry() << "\nэксцесс:" << md.excess() << endl << endl;
		cout << "эмпирическое распределение для смеси1:\nмат. ожидание:" << ed.expected_value() << "\nдисперсия: " << ed.dispersion() << "\nассиметрия: " << ed.asymmetry() << "\nэксцесс:" << ed.excess() << endl << endl;
		cout << "2:\nмат. ожидание:" << ed2.expected_value() << "\nдисперсия: " << ed2.dispersion() << "\nассиметрия: " << ed2.asymmetry() << "\nэксцесс:" << ed2.excess() << endl;



		for (int i = 0; i < 1000; ++i) {
			file1 << graph1[i].first << "\t" << graph1[i].second << endl;
			file2 << graph2[i].first << "\t" << graph2[i].second << endl;
		}
		///
	}
	catch (const int error) {
		if (error == 0) {
			cout << "ошибка открытия params.txt!!!" << endl;
		}
		else if (error == 1) {
			cout << "ошибка в данных в файле params.txt" << endl;
		}
	}
	file1.close();
	file2.close();
	run_unit_tests(argc, argv);
}
