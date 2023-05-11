//
// Created by juli on 18.05.22.
//
#include "types.hpp"
#include <sstream>
std::ostream &operator<<(std::ostream &os, const epi::options::EpistasisModel &model) {
    switch (model) {
        case epi::options::EpistasisModel::BAYESIAN_MODEL:
            os << "BAYESIAN";
            break;
        case epi::options::EpistasisModel::PENETRANCE_MODEL:
            os << "PENETRANCE";
            break;
        case epi::options::EpistasisModel::REGRESSION_MODEL:
            os << "REGRESSION";
            break;
        case epi::options::EpistasisModel::VARIANCE_MODEL:
            os << "VARIANCE";
            break;
    }
    return os;
}

std::ostream &operator<<(std::ostream &os, const epi::options::EpistasisScore &model) {
    switch (model) {
        case epi::options::EpistasisScore::BAYESIAN:
            os << "BAYESIAN";
            break;
        case epi::options::EpistasisScore::PENETRANCE_NLL:
            os << "PENETRANCE_NLL";
            break;
        case epi::options::EpistasisScore::PENETRANCE_LLH:
            os << "PENETRANCE_LLH";
            break;
        case epi::options::EpistasisScore::PENETRANCE_AIC:
            os << "PENETRANCE_AIC";
            break;
        case epi::options::EpistasisScore::PENETRANCE_BIC:
            os << "PENETRANCE_BIC";
            break;
        case epi::options::EpistasisScore::REGRESSION_NLL:
            os << "REGRESSION_NLL";
            break;
        case epi::options::EpistasisScore::REGRESSION_LLH:
            os << "REGRESSION_LLH";
            break;
        case epi::options::EpistasisScore::REGRESSION_AIC:
            os << "REGRESSION_AIC";
            break;
        case epi::options::EpistasisScore::REGRESSION_BIC:
            os << "REGRESSION_BIC";
            break;
        case epi::options::EpistasisScore::REGRESSION_NLL_GAIN:
            os << "REGRESSION_NLL-GAIN";
            break;
        case epi::options::EpistasisScore::REGRESSION_LLH_GAIN:
            os << "REGRESSION_LLH-GAIN";
            break;
        case epi::options::EpistasisScore::REGRESSION_AIC_GAIN:
            os << "REGRESSION_AIC-GAIN";
            break;
        case epi::options::EpistasisScore::REGRESSION_BIC_GAIN:
            os << "REGRESSION_BIC-GAIN";
            break;
        case epi::options::EpistasisScore::VARIANCE:
            os << "VARIANCE";
            break;
    }
    return os;
}

namespace epi {
    namespace options {
        EpistasisModel epistasis_model_from_epistasis_score(EpistasisScore score) {
            switch (score) {
                case EpistasisScore::BAYESIAN:
                    return EpistasisModel::BAYESIAN_MODEL;
                case EpistasisScore::PENETRANCE_NLL:
                case EpistasisScore::PENETRANCE_LLH:
                case EpistasisScore::PENETRANCE_AIC:
                case EpistasisScore::PENETRANCE_BIC:
                    return EpistasisModel::PENETRANCE_MODEL;
                case EpistasisScore::REGRESSION_NLL:
                case EpistasisScore::REGRESSION_LLH:
                case EpistasisScore::REGRESSION_AIC:
                case EpistasisScore::REGRESSION_BIC:
                case EpistasisScore::REGRESSION_NLL_GAIN:
                case EpistasisScore::REGRESSION_LLH_GAIN:
                case EpistasisScore::REGRESSION_AIC_GAIN:
                case EpistasisScore::REGRESSION_BIC_GAIN:
                    return EpistasisModel::REGRESSION_MODEL;
                case EpistasisScore::VARIANCE:
                    return EpistasisModel::VARIANCE_MODEL;
                default:
                    return EpistasisModel::BAYESIAN_MODEL;
            }
        }

        EpistasisModel epistasis_model_from_string(const std::string &model_string) {
            if (model_string == "BAYESIAN") {
                return EpistasisModel::BAYESIAN_MODEL;
            } else if (model_string == "PENETRANCE") {
                return EpistasisModel::PENETRANCE_MODEL;
            } else if (model_string.substr(0, 11) == "PENETRANCE_") {
                std::vector<std::string> scores = {"NLL", "LLH", "AIC", "BIC"};
                if (std::find(scores.begin(), scores.end(), model_string.substr(12)) == scores.end()) {
                    throw epi::Error("Unknown PENETRANCE score " + model_string);
                }
                return EpistasisModel::PENETRANCE_MODEL;
            } else if (model_string == "REGRESSION") {
                return EpistasisModel::REGRESSION_MODEL;
            } else if (model_string.substr(0, 11) == "REGRESSION_") {
                std::vector<std::string> scores = {"NLL", "LLH", "AIC", "BIC", "NLL-GAIN", "LLH-GAIN", "AIC-GAIN",
                                                   "BIC-GAIN"};
                if (std::find(scores.begin(), scores.end(), model_string.substr(12)) == scores.end()) {
                    throw epi::Error("Unknown REGRESSION score " + model_string);
                }
                return EpistasisModel::REGRESSION_MODEL;
            } else if (model_string == "VARIANCE") {
                return EpistasisModel::VARIANCE_MODEL;
            } else {
                throw epi::Error("Unknown epistasis model " + model_string);
            }
        }

        EpistasisScore epistasis_score_from_string(const std::string &model_string) {
            if (model_string == "BAYESIAN") {
                return EpistasisScore::BAYESIAN;
            } else if (model_string == "PENETRANCE" || model_string == "PENETRANCE_NLL") {
                return EpistasisScore::PENETRANCE_NLL;
            } else if (model_string == "PENETRANCE_LLH") {
                return EpistasisScore::PENETRANCE_LLH;
            } else if (model_string == "PENETRANCE_AIC") {
                return EpistasisScore::PENETRANCE_AIC;
            } else if (model_string == "PENETRANCE_BIC") {
                return EpistasisScore::PENETRANCE_BIC;
            } else if (model_string == "REGRESSION" || model_string == "REGRESSION_NLL") {
                return EpistasisScore::REGRESSION_NLL;
            } else if (model_string == "REGRESSION_LLH") {
                return EpistasisScore::REGRESSION_LLH;
            } else if (model_string == "REGRESSION_AIC") {
                return EpistasisScore::REGRESSION_AIC;
            } else if (model_string == "REGRESSION_BIC") {
                return EpistasisScore::REGRESSION_BIC;
            } else if (model_string == "REGRESSION_NLL-GAIN") {
                return EpistasisScore::REGRESSION_NLL_GAIN;
            } else if (model_string == "REGRESSION_LLH-GAIN") {
                return EpistasisScore::REGRESSION_LLH_GAIN;
            } else if (model_string == "REGRESSION_AIC-GAIN") {
                return EpistasisScore::REGRESSION_AIC_GAIN;
            } else if (model_string == "REGRESSION_BIC-GAIN") {
                return EpistasisScore::REGRESSION_BIC_GAIN;
            } else if (model_string == "VARIANCE") {
                return EpistasisScore::VARIANCE;
            } else {
                throw epi::Error("Unknown epistasis score " + model_string);
            }
        }

        std::vector<std::string> get_all_epistasis_scores() {
            // return { "PENETRANCE" };
            return std::vector<std::string>{"BAYESIAN",
                   "PENETRANCE_NLL", "PENETRANCE_LLH", "PENETRANCE_AIC", "PENETRANCE_BIC",
                   "REGRESSION_NLL", "REGRESSION_LLH", "REGRESSION_AIC", "REGRESSION_BIC",
                   "REGRESSION_NLL-GAIN", "REGRESSION_LLH-GAIN", "REGRESSION_AIC-GAIN", "REGRESSION_BIC-GAIN",
                   "VARIANCE"};
        }

        std::string epistasis_score_to_string(EpistasisScore model) {
            std::stringstream ss;
            ss << model;
            return ss.str();
        }

        std::string epistasis_model_to_string(EpistasisModel model) {
            std::stringstream ss;
            ss << model;
            return ss.str();
        }

        bool is_epistasis_model_string(const std::string &model_string) {
            if (model_string == "BAYESIAN" || model_string == "VARIANCE") {
                return true;
            } else if (model_string.substr(0, 11) == "PENETRANCE_") {
                std::vector<std::string> scores = {"NLL", "LLH", "AIC", "BIC"};
                return std::find(scores.begin(), scores.end(), model_string.substr(11)) != scores.end();
            } else if (model_string.substr(0, 11) == "REGRESSION_") {
                std::vector<std::string> scores = {"NLL", "LLH", "AIC", "BIC", "NLL-GAIN", "LLH-GAIN", "AIC-GAIN",
                                                   "BIC-GAIN"};
                return std::find(scores.begin(), scores.end(), model_string.substr(11)) != scores.end();
            } else {
                return false;
            }
        }
    }
}
