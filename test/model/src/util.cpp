//
// Created by juli on 18.05.22.
//

#include "util.hpp"

using namespace epi::options;

std::ostream & operator<<(std::ostream & os, const epi::ModelTestComparedModel & model) {
    os << model.model << "[" << model.options << "],";
    os << model.mean_runtime << ",";
    os << -std::log(model.p_value_from_score) << ",";
    os << model.mean_correlation << ",";
    if (model.predict) {
        os << model.mean_prediction_time << ",";
        os << model.prediction_score << "\n";
    }
    else {
        os << "NA,NA\n";
    }
    return os;
}
