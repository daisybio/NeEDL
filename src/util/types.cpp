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
        case epi::options::EpistasisScore::BAYESIAN_COV:
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
        case epi::options::EpistasisScore::PENETRANCE_COV_NLL:
            os << "PENETRANCE_COV_NLL";
            break;
        case epi::options::EpistasisScore::PENETRANCE_COV_LLH:
            os << "PENETRANCE_COV_LLH";
            break;
        case epi::options::EpistasisScore::PENETRANCE_COV_AIC:
            os << "PENETRANCE_COV_AIC";
            break;
        case epi::options::EpistasisScore::PENETRANCE_COV_BIC:
            os << "PENETRANCE_COV_BIC";
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
        case epi::options::EpistasisScore::REGRESSION_COV_NLL:
            os << "REGRESSION_COV_NLL";
            break;
        case epi::options::EpistasisScore::REGRESSION_COV_LLH:
            os << "REGRESSION_COV_LLH";
            break;
        case epi::options::EpistasisScore::REGRESSION_COV_AIC:
            os << "REGRESSION_COV_AIC";
            break;
        case epi::options::EpistasisScore::REGRESSION_COV_BIC:
            os << "REGRESSION_COV_BIC";
            break;
        case epi::options::EpistasisScore::REGRESSION_CLG_Q_LC:
            os << "REGRESSION_CLG_Q_LC";
            break;
        case epi::options::EpistasisScore::REGRESSION_CLG_Q_QC:
            os << "REGRESSION_CLG_Q_QC";
            break;
        case epi::options::EpistasisScore::REGRESSION_CLG_L_LC:
            os << "REGRESSION_CLG_L_LC";
            break;
        case epi::options::EpistasisScore::VARIANCE:
            os << "VARIANCE";
            break;
        case epi::options::EpistasisScore::VARIANCE_COV:
            os << "VARIANCE_COV";
            break;
    }
    return os;
}

namespace epi {
    namespace options {
        EpistasisModel epistasis_model_from_epistasis_score(EpistasisScore score) {
            switch (score) {
                case EpistasisScore::BAYESIAN:
                case EpistasisScore::BAYESIAN_COV:
                    return EpistasisModel::BAYESIAN_MODEL;
                case EpistasisScore::PENETRANCE_NLL:
                case EpistasisScore::PENETRANCE_LLH:
                case EpistasisScore::PENETRANCE_AIC:
                case EpistasisScore::PENETRANCE_BIC:
                case EpistasisScore::PENETRANCE_COV_NLL:
                case EpistasisScore::PENETRANCE_COV_LLH:
                case EpistasisScore::PENETRANCE_COV_AIC:
                case EpistasisScore::PENETRANCE_COV_BIC:
                    return EpistasisModel::PENETRANCE_MODEL;
                case EpistasisScore::REGRESSION_NLL:
                case EpistasisScore::REGRESSION_LLH:
                case EpistasisScore::REGRESSION_AIC:
                case EpistasisScore::REGRESSION_BIC:
                case EpistasisScore::REGRESSION_NLL_GAIN:
                case EpistasisScore::REGRESSION_LLH_GAIN:
                case EpistasisScore::REGRESSION_AIC_GAIN:
                case EpistasisScore::REGRESSION_BIC_GAIN:
                case EpistasisScore::REGRESSION_COV_NLL:
                case EpistasisScore::REGRESSION_COV_LLH:
                case EpistasisScore::REGRESSION_COV_AIC:
                case EpistasisScore::REGRESSION_COV_BIC:
                case EpistasisScore::REGRESSION_CLG_Q_LC:
                case EpistasisScore::REGRESSION_CLG_Q_QC:
                case EpistasisScore::REGRESSION_CLG_L_LC:
                    return EpistasisModel::REGRESSION_MODEL;
                case EpistasisScore::VARIANCE:
                case EpistasisScore::VARIANCE_COV:
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
                std::vector<std::string> scores = {"NLL", "LLH", "AIC", "BIC", "COV_NLL", "COV_LLH", "COV_AIC", "COV_BIC"};
                if (std::find(scores.begin(), scores.end(), model_string.substr(12)) == scores.end()) {
                    throw epi::Error("Unknown PENETRANCE score " + model_string);
                }
                return EpistasisModel::PENETRANCE_MODEL;
            } else if (model_string == "REGRESSION") {
                return EpistasisModel::REGRESSION_MODEL;
            } else if (model_string.substr(0, 11) == "REGRESSION_") {
                std::vector<std::string> scores = {"NLL", "LLH", "AIC", "BIC", "COV_NLL", "COV_LLH", "COV_AIC", "COV_BIC", "NLL-GAIN", "LLH-GAIN", "AIC-GAIN",
                                                   "BIC-GAIN", "CLG-L-LC", "CLG-Q-QC", "CLG-Q-LC"};
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
            } else if (model_string == "BAYESIAN_COV") {
                return EpistasisScore::BAYESIAN_COV;
            } else if (model_string == "PENETRANCE" || model_string == "PENETRANCE_NLL") {
                return EpistasisScore::PENETRANCE_NLL;
            } else if (model_string == "PENETRANCE_LLH") {
                return EpistasisScore::PENETRANCE_LLH;
            } else if (model_string == "PENETRANCE_AIC") {
                return EpistasisScore::PENETRANCE_AIC;
            } else if (model_string == "PENETRANCE_BIC") {
                return EpistasisScore::PENETRANCE_BIC;
            } else if (model_string == "PENETRANCE_COV" || model_string == "PENETRANCE_COV_NLL") {
                return EpistasisScore::PENETRANCE_COV_NLL;
            } else if (model_string == "PENETRANCE_COV_LLH") {
                return EpistasisScore::PENETRANCE_COV_LLH;
            } else if (model_string == "PENETRANCE_COV_AIC") {
                return EpistasisScore::PENETRANCE_COV_AIC;
            } else if (model_string == "PENETRANCE_COV_BIC") {
                return EpistasisScore::PENETRANCE_COV_BIC;
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
            } else if (model_string == "REGRESSION_COV" || model_string == "REGRESSION_COV_NLL") {
                return EpistasisScore::REGRESSION_COV_NLL;
            } else if (model_string == "REGRESSION_COV_LLH") {
                return EpistasisScore::REGRESSION_COV_LLH;
            } else if (model_string == "REGRESSION_COV_AIC") {
                return EpistasisScore::REGRESSION_COV_AIC;
            } else if (model_string == "REGRESSION_COV_BIC") {
                return EpistasisScore::REGRESSION_COV_BIC;
            } else if (model_string == "REGRESSION_CLG_Q_LC") {
                return EpistasisScore::REGRESSION_CLG_Q_LC;
            } else if (model_string == "REGRESSION_CLG_Q_QC") {
                return EpistasisScore::REGRESSION_CLG_Q_QC;
            } else if (model_string == "REGRESSION_CLG_L_LC") {
                return EpistasisScore::REGRESSION_CLG_L_LC;
            } else if (model_string == "VARIANCE") {
                return EpistasisScore::VARIANCE;
            } else if (model_string == "VARIANCE_COV") {
                return EpistasisScore::VARIANCE_COV;
            } else {
                throw epi::Error("Unknown epistasis score " + model_string);
            }
        }

        std::vector<std::string> get_all_epistasis_scores(bool include_covariates) {
            // return { "PENETRANCE" };
            auto without_cov = std::vector<std::string>{
                "BAYESIAN",
                "PENETRANCE_NLL", "PENETRANCE_LLH", "PENETRANCE_AIC", "PENETRANCE_BIC",
                "REGRESSION_NLL", "REGRESSION_LLH", "REGRESSION_AIC", "REGRESSION_BIC",
                "REGRESSION_NLL-GAIN", "REGRESSION_LLH-GAIN", "REGRESSION_AIC-GAIN", "REGRESSION_BIC-GAIN",
                "VARIANCE"
            };

            if (include_covariates) {
                auto with_cov = std::vector<std::string>{
                    "BAYESIAN_COV",
                    "PENETRANCE_COV_NLL", "PENETRANCE_COV_LLH", "PENETRANCE_COV_AIC", "PENETRANCE_COV_BIC",
                    "REGRESSION_COV_NLL", "REGRESSION_COV_LLH", "REGRESSION_COV_AIC", "REGRESSION_COV_BIC",
                    "REGRESSION_CLG_Q_LC", "REGRESSION_CLG_Q_QC","REGRESSION_CLG_L_LC",
                    "VARIANCE_COV"
                };

                without_cov.insert(without_cov.begin(), with_cov.begin(), with_cov.end());
            }

            return without_cov;
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
            if (model_string == "BAYESIAN" || model_string == "VARIANCE" || model_string == "BAYESIAN_COV" || model_string == "VARIANCE_COV") {
                return true;
            } else if (model_string.substr(0, 11) == "PENETRANCE_") {
                std::vector<std::string> scores = {"NLL", "LLH", "AIC", "BIC", "COV_NLL", "COV_LLH", "COV_AIC", "COV_BIC"};
                return std::find(scores.begin(), scores.end(), model_string.substr(12)) != scores.end();
            } else if (model_string.substr(0, 11) == "REGRESSION_") {
                std::vector<std::string> scores = {"NLL", "LLH", "AIC", "BIC", "COV_NLL", "COV_LLH", "COV_AIC", "COV_BIC", "NLL-GAIN", "LLH-GAIN", "AIC-GAIN",
                                                   "BIC-GAIN", "CLG-L-LC", "CLG-Q-QC", "CLG-Q-LC"};
                return std::find(scores.begin(), scores.end(), model_string.substr(12)) != scores.end();
            } else {
                return false;
            }
        }
    }
}
