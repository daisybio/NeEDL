//
// Created by juli on 19.05.22.
//

#include "regression_model.hpp"

namespace epi {
    template<>
    QuantitativePhenoType
    RegressionModel<QuantitativePhenoType>::
    predict_(Ind ind) const {
        Eigen::MatrixXd interaction_features;
        construct_features_(snp_set_, ind, true, interaction_features);
        return (interaction_features * interaction_model_)(0, 0);
    }

    template<>
    CategoricalPhenoType
    RegressionModel<CategoricalPhenoType>::
    predict_(Ind ind) const {
        Eigen::MatrixXd interaction_features;
        construct_features_(snp_set_, ind, true, interaction_features);
        if (this->instance_->num_categories() == 2) {
            if (1.0 / (1.0 + std::exp(-(interaction_features * interaction_model_)(0, 0))) < 0.5) {
                return 0;
            }
            return 1;
        }
        Eigen::MatrixXd interaction_hypothesis(1, this->instance_->num_categories());
        interaction_hypothesis = (interaction_features * interaction_model_).array().exp().matrix();
        interaction_hypothesis /= interaction_hypothesis.sum();
        CategoricalPhenoType prediction{0};
        interaction_hypothesis.row(0).maxCoeff(&prediction);
        return prediction;
    }

    template<>
    void
    RegressionModel<QuantitativePhenoType>::
    compute_linear_model_sigma_(const Eigen::MatrixXd & hypothesis, double & sigma) const {
        double error{0.0};
        sigma = 0;
        for (Ind ind{0}; ind < this->instance_->num_inds(); ind++) {
            error = hypothesis(ind, 0) - this->instance_->phenotype(ind);
            sigma += error * error;
        }
        sigma /= static_cast<double>(this->instance_->num_inds());
    }

    template<>
    void
    RegressionModel<QuantitativePhenoType>::
    fit_linear_model_(const Eigen::MatrixXd & feature_matrix, Eigen::MatrixXd & model, Eigen::MatrixXd & hypothesis) const {
        std::size_t num_features{static_cast<std::size_t>(feature_matrix.cols())};
        std::size_t num_inds{this->instance_->num_inds()};
        model.resize(num_features, 1);
        hypothesis.resize(num_inds, 1);
        double loss{std::numeric_limits<double>::max()};
        double loss_old{loss};
        Eigen::MatrixXd gradient(num_features, 1);
        model.setOnes();
        double error{0};
        for (std::size_t itr{0}; itr < max_itrs_; itr++) {
            loss_old = loss;
            hypothesis = feature_matrix * model;
            loss = 0;
            for (Ind ind{0}; ind < num_inds; ind++) {
                error = hypothesis(ind, 0) - this->instance_->phenotype(ind);
                loss += error * error;
            }
            if (loss_old - loss < epsilon_) {
                break;
            }
            gradient = Eigen::VectorXd::Zero(num_features);
            for (Ind ind{0}; ind < num_inds; ind++) {
                gradient += (hypothesis(ind, 0) - this->instance_->phenotype(ind)) * feature_matrix.row(ind).transpose();
            }
            gradient /= static_cast<double>(num_inds);
            model -= learning_rate_ * gradient;
        }
    }

    template<>
    void
    RegressionModel<CategoricalPhenoType>::
    fit_logistic_model_(const Eigen::MatrixXd & feature_matrix, Eigen::MatrixXd & model, Eigen::MatrixXd & hypothesis) const {
        std::size_t num_features{static_cast<std::size_t>(feature_matrix.cols())};
        std::size_t num_inds{this->instance_->num_inds()};
        model.resize(num_features, 1);
        hypothesis.resize(num_inds, 1);
        double loss{std::numeric_limits<double>::max()};
        double loss_old{loss};
        Eigen::MatrixXd gradient(num_features, 1);
        model.setOnes();
        for (std::size_t itr{0}; itr < max_itrs_ and loss > epsilon_; itr++) {
            loss_old = loss;
            hypothesis = ((-(feature_matrix * model)).array().exp() + 1.0).inverse().matrix();
            loss = 0;
            for (Ind ind{0}; ind < num_inds; ind++) {
                if (this->instance_->phenotype(ind) == 0) {
                    loss -= std::log(1.0 - hypothesis(ind, 0));
                }
                else {
                    loss -= std::log(hypothesis(ind, 0));
                }
            }
            if (loss_old - loss < epsilon_) {
                break;
            }
            gradient.setZero();
            for (Ind ind{0}; ind < num_inds; ind++) {
                if (this->instance_->phenotype(ind) == 0) {
                    gradient += hypothesis(ind, 0) * feature_matrix.row(ind).transpose();
                }
                else {
                    gradient += (hypothesis(ind, 0) - 1.0) * feature_matrix.row(ind).transpose();
                }
            }
            gradient /= static_cast<double>(num_inds);
            model -= learning_rate_ * gradient;
        }
        Eigen::MatrixXd temp(hypothesis);
        hypothesis.resize(num_inds, 2);
        hypothesis.col(0) = Eigen::MatrixXd::Ones(num_inds, 1) - temp.col(0);
        hypothesis.col(1) = temp.col(0);
    }

    template<>
    void
    RegressionModel<QuantitativePhenoType>::
    fit_model_(const Eigen::MatrixXd & feature_matrix, Eigen::MatrixXd & model, Eigen::MatrixXd & hypothesis, double & sigma) const {
        fit_linear_model_(feature_matrix, model, hypothesis);
        compute_linear_model_sigma_(hypothesis, sigma);
    }

    template<>
    void
    RegressionModel<CategoricalPhenoType>::
    fit_multinomial_model_(const Eigen::MatrixXd & feature_matrix, Eigen::MatrixXd & model, Eigen::MatrixXd & hypothesis) const {
        std::size_t num_features{static_cast<std::size_t>(feature_matrix.cols())};
        std::size_t num_inds{this->instance_->num_inds()};
        std::size_t num_categories{this->instance_->num_categories()};
        model.resize(num_features, num_categories);
        hypothesis.resize(num_inds, num_categories);
        double loss{std::numeric_limits<double>::max()};
        double loss_old{loss};
        Eigen::MatrixXd gradient(num_features, num_categories);
        model.setOnes();
        for (std::size_t itr{0}; itr < max_itrs_ and loss > epsilon_; itr++) {
            loss_old = loss;
            hypothesis = (feature_matrix * model).array().exp().matrix();
            for (Ind ind{0}; ind < num_inds; ind++) {
                hypothesis.row(ind) /= hypothesis.row(ind).sum();
            }
            loss = 0;
            for (Ind ind{0}; ind < num_inds; ind++) {
                loss -= std::log(hypothesis(ind, this->instance_->phenotype(ind)));
            }
            if (loss_old - loss < epsilon_) {
                break;
            }
            gradient.setZero();
            for (Ind ind{0}; ind < num_inds; ind++) {
                for (std::size_t k{0}; k < num_categories; k++) {
                    if (k == static_cast<std::size_t>(this->instance_->phenotype(ind))) {
                        gradient.col(k) += (hypothesis(ind, k) - 1.0) * feature_matrix.row(ind).transpose();
                    }
                    else {
                        gradient.col(k) += hypothesis(ind, k) * feature_matrix.row(ind).transpose();
                    }
                }
            }
            gradient /= static_cast<double>(num_inds);
            model -= learning_rate_ * gradient;
        }
    }

    template<>
    void
    RegressionModel<CategoricalPhenoType>::
    fit_model_(const Eigen::MatrixXd & feature_matrix, Eigen::MatrixXd & model, Eigen::MatrixXd & hypothesis, double & sigma) const {
        if (this->instance_->num_categories() == 2) {
            fit_logistic_model_(feature_matrix, model, hypothesis);
        }
        else {
            fit_multinomial_model_(feature_matrix, model, hypothesis);
        }
    }

    template<>
    double
    RegressionModel<QuantitativePhenoType>::
    likelihood_(Ind ind, const Eigen::MatrixXd & hypothesis, double sigma) const {
        return misc::normal_pdf(this->instance_->phenotype(ind), hypothesis(ind, 0), sigma);
    }

    template<>
    double
    RegressionModel<CategoricalPhenoType>::
    likelihood_(Ind ind, const Eigen::MatrixXd & hypothesis, double sigma) const {
        return misc::ensure_valid_continuous_probability(hypothesis(ind, this->instance_->phenotype(ind)));
    }

    template<>
    std::size_t
    RegressionModel<QuantitativePhenoType>::
    degrees_of_freedom_per_feature_() const {
        return 1;
    }

    template<>
    std::size_t
    RegressionModel<CategoricalPhenoType>::
    degrees_of_freedom_per_feature_() const {
        return this->instance_->num_categories() - 1;
    }

    template<>
    std::string
    RegressionModel<QuantitativePhenoType>::
    meta_info_() const {
        return "phenotype = quantitative\n";
    }

    template<>
    std::string
    RegressionModel<CategoricalPhenoType>::
    meta_info_() const {
        return std::string("phenotype = categorical\ncategories = ") + std::to_string(this->instance_->num_categories()) + "\n";
    }

    template<>
    void
    RegressionModel<QuantitativePhenoType>::
    check_loaded_model_() const {
        std::size_t num_features{snp_set_.size() + 1 + (snp_set_.size() * (snp_set_.size() - 1)) / 2};
        if (static_cast<std::size_t>(interaction_model_.rows()) != snp_set_.size() + 1 + (snp_set_.size() * (snp_set_.size() - 1)) / 2) {
            throw Error(std::string("Wrong number of features in loaded model for SNP set of size ") + std::to_string(snp_set_.size()) + ". Expected: " + std::to_string(num_features) + ". Actual: " + std::to_string(interaction_model_.rows()) + ".");
        }
        if (interaction_model_.cols() != 1) {
            throw Error(std::string("Wrong number of parameters per feature. Expected: 1. Actual: ") + std::to_string(interaction_model_.cols()) + ".");
        }
    }


    template<>
    void
    RegressionModel<CategoricalPhenoType>::
    check_loaded_model_() const {
        std::size_t num_features{snp_set_.size() + 1 + (snp_set_.size() * (snp_set_.size() - 1)) / 2};
        if (static_cast<std::size_t>(interaction_model_.rows()) != snp_set_.size() + 1 + (snp_set_.size() * (snp_set_.size() - 1)) / 2) {
            throw Error(std::string("Wrong number of features in loaded model for SNP set of size ") + std::to_string(snp_set_.size()) + ". Expected: " + std::to_string(num_features) + ". Actual: " + std::to_string(interaction_model_.rows()) + ".");
        }
        if (this->instance_->num_categories() == 2) {
            if (interaction_model_.cols() != 1) {
                throw Error(std::string("Wrong number of parameters per feature. Expected: 1. Actual: ") + std::to_string(interaction_model_.cols()) + ".");
            }
        }
        else if (static_cast<std::size_t>(interaction_model_.cols()) != this->instance_->num_categories()) {
            throw Error(std::string("Wrong number of parameters per feature. Expected: ") + std::to_string(this->instance_->num_categories()) + ". Actual: " + std::to_string(interaction_model_.cols()) + ".");
        }
    }

}