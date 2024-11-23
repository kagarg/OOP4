#include "catch.hpp"
#include "Standart_distr.h"
#include "Mixed_distr.cpp"
#include "Empiric.h"

bool equal(const double& x, const double& y) {
    if (abs(x - y) <= 0.1) {
        return true;
    }
    else return false;
}

TEST_CASE("basic methods") {
    Std_distr distr;
    CHECK(distr.get_form() == 1);
    CHECK(distr.get_shift() == 0);
    CHECK(distr.get_scale() == 1);
}

TEST_CASE("standart distribution") {
    Std_distr distr;
    CHECK(distr.density(0) == 0.5);
    CHECK(distr.expected_value() == 0);
    CHECK(distr.dispersion() == 2);
    CHECK(distr.asymmetry() == 0);
    CHECK(distr.excess() == 3);
}

TEST_CASE("shift scale transformation") {
    Std_distr distr;
    distr.set_scale(2);
    distr.set_shift(2);
    CHECK(equal(distr.density(0), 0.091) == true);
    CHECK(distr.expected_value() == 2);
    CHECK(distr.dispersion() == 8);
    CHECK(distr.asymmetry() == 0);
    CHECK(distr.excess() == 3);
}

TEST_CASE("mixture distribution") {
    Std_distr distr1;
    Std_distr distr2;
    distr1.set_scale(2);
    distr2.set_scale(2);
    distr1.set_shift(2);
    distr2.set_shift(2);
    MixtureDistribution< Std_distr, Std_distr> mdistr(distr1, distr2, 0.5);

    CHECK(equal(mdistr.density(0), 0.091) == true);
    CHECK(mdistr.expected_value() == 2);
    CHECK(mdistr.dispersion() == 8);
    CHECK(mdistr.asymmetry() == 0);
    CHECK(equal(mdistr.excess(), 2.953) == true);
}

TEST_CASE("mixture distribution expected") {
    Std_distr distr1;
    Std_distr distr2;
    distr1.set_scale(2);
    distr2.set_scale(2);
    distr1.set_shift(1);
    distr2.set_shift(2);
    MixtureDistribution< Std_distr, Std_distr> mdistr(distr1, distr2, 0.5);
    CHECK(mdistr.expected_value() == 1.5);
}

TEST_CASE("mixture distribution dispersion") {
    Std_distr distr1;
    Std_distr distr2;
    distr1.set_scale(1);
    distr2.set_scale(3);
    MixtureDistribution< Std_distr, Std_distr> mdistr(distr1, distr2, 0.5);
    CHECK(mdistr.dispersion() == 10);
}

TEST_CASE("late binding mechanism") {
    Std_distr distr1;
    Std_distr distr2;
    distr1.set_scale(1);
    distr2.set_scale(3);

    MixtureDistribution< Std_distr, Std_distr> mdistr(distr1, distr2, 0.5);
    MixtureDistribution< Std_distr, MixtureDistribution< Std_distr, Std_distr>> mdistr2(distr1, mdistr, 0.5);
    CHECK(mdistr2.component1().get_form() == 1);
    CHECK(mdistr2.component1().get_shift() == 0);
    CHECK(mdistr2.component1().get_scale() == 1);
    CHECK(mdistr2.component2().component1().get_form() == 1);
    CHECK(mdistr2.component2().component1().get_shift() == 0);
    CHECK(mdistr2.component2().component1().get_scale() == 1);
    CHECK(mdistr2.component2().component2().get_form() == 1);
    CHECK(mdistr2.component2().component2().get_shift() == 0);
    CHECK(mdistr2.component2().component2().get_scale() == 3);
}

TEST_CASE("empirical distribution") {
    Std_distr distr1;
    EmpiricalDistribution ed(distr1, 200, 1);
    CHECK(ed.get_size() == 200);
    CHECK(ed.get_intrevals_number() == 8);
    /*CHECK(equal(ed.expected_value(), 0.0) == true);*/
}
