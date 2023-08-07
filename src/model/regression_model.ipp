/*******************************************************************************
 *                                                                             *
 *   Copyright (C) 2020 by David B. Blumenthal                                 *
 *                                                                             *
 *   This file is part of GenEpiSeeker.                                        *
 *                                                                             *
 *   GenEpiSeeker is free software: you can redistribute it and/or modify it   *
 *   under the terms of the GNU General Public License as published by         *
 *   the Free Software Foundation, either version 3 of the License, or         *
 *   (at your option) any later version.                                       *
 *                                                                             *
 *   GenEpiSeeker is distributed in the hope that it will be useful,           *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the              *
 *   GNU General Public License for more details.                              *
 *                                                                             *
 *   You should have received a copy of the GNU General Public License         *
 *   along with GenEpiSeeker. If not, see <http://www.gnu.org/licenses/>.      *
 *                                                                             *
 ******************************************************************************/


/*!
 * @file regression_model.ipp
 * @brief Definition of epi::RegressionModel.
 */

#ifndef SRC_MODEL_REGRESSION_MODEL_IPP_
#define SRC_MODEL_REGRESSION_MODEL_IPP_

#include "regression_model.hpp"

namespace epi {

    template<class PhenoType>
    RegressionModel<PhenoType>::
    RegressionModel(Instance<PhenoType> * instance) :
            EpistasisModel<PhenoType>(instance),
            score_("NLL"),
            cov_model_(false),
            incl_cov_(false),
            learning_rate_{0.1},
            epsilon_{10e-5},
            max_itrs_{50000},
            snp_set_(),
            interaction_model_() {}

    template<class PhenoType>
    RegressionModel<PhenoType>::
    ~RegressionModel() {}

    template<class PhenoType>
    options::ModelSense
    RegressionModel<PhenoType>::
    model_sense_() const {
        if (gain_score_() or score_ == "LLH") {
            return options::ModelSense::MAXIMIZE;
        }
        return options::ModelSense::MINIMIZE;
    }

    template<class PhenoType>
    bool
    RegressionModel<PhenoType>::
    is_predictive_() const {
        return true;
    }



    template<class PhenoType>
    bool
    RegressionModel<PhenoType>::
    parse_option_(const std::string & option, const std::string & arg) {
        if (option == "score") {
            score_ = arg;
            if (score_ != "LLH" and score_ != "NLL" and score_ != "AIC" and score_ != "BIC" and not gain_score_() and not cov_score_()) {
                throw Error(std::string("Invalid argument \"") + arg  + "\" for option \"--" + option + "\". Usage: options = \"[--" + option + " LLH|LLH-GAIN|NLL|NLL-GAIN|AIC|AIC-GAIN|BIC|BIC-GAIN|CLG-L-LC|CLG-Q-QC|CLG-Q-LC] [...]\"");
            }
            return true;
        }
        if (option == "max-itrs") {
            try {
                max_itrs_ = std::stoul(arg);
            }
            catch (...) {
                throw Error(std::string("Invalid argument \"") + arg  + "\" for option \"--" + option + "\". Usage: options = \"[--" + option + " <convertible to int greater equal 0>] [...]\"");
            }
            if (max_itrs_ == 0) {
                max_itrs_ = std::numeric_limits<std::size_t>::max();
            }
            return true;
        }
        if (option == "learning-rate") {
            try {
                learning_rate_ = std::stod(arg);
            }
            catch (...) {
                throw Error(std::string("Invalid argument \"") + arg  + "\" for option \"--" + option + "\". Usage: options = \"[--" + option + " <convertible to double greater 0>] [...]\"");
            }
            if (learning_rate_ <= 0) {
                throw Error(std::string("Invalid argument \"") + arg  + "\" for option \"--" + option + "\". Usage: options = \"[--" + option + " <convertible to double greater 0>] [...]\"");
            }
            return true;
        }
        if (option == "epsilon") {
            try {
                epsilon_ = std::stod(arg);
            }
            catch (...) {
                throw Error(std::string("Invalid argument \"") + arg  + "\" for option \"--" + option + "\". Usage: options = \"[--" + option + " <convertible to double greater 0>] [...]\"");
            }
            if (epsilon_ <= 0) {
                throw Error(std::string("Invalid argument \"") + arg  + "\" for option \"--" + option + "\". Usage: options = \"[--" + option + " <convertible to double greater 0>] [...]\"");
            }
            return true;
        }
        return false;
    }

    template<class PhenoType>
    void
    RegressionModel<PhenoType>::
    set_default_options_() {
        score_ = "NLL";
        incl_cov_ = false;
        learning_rate_ = 0.1;
        epsilon_ = 10e-5;
        max_itrs_ = 50000;
    }


    template<class PhenoType>
    std::string
    RegressionModel<PhenoType>::
    valid_options_() const {
        return "[--score LLH|LLH-GAIN|NLL|NLL-GAIN|AIC|AIC-GAIN|BIC|BIC-GAIN|CLG-L-LC|CLG-Q-QC|CLG-Q-LC] [--max-itrs <convertible to int greater equal 0>] [--learning-rate <convertible to double greater 0>] [--epsilon <convertible to double greater 0>] [--cov]";
    }

    template<class PhenoType>
    void
    RegressionModel<PhenoType>::
    save_(const std::string & filename) const {
        std::ofstream file(filename.c_str());
        file << "[meta]\n";
        file << "model-type = RegressionModel\n";
        file << "model-cov = " << this->get_cov_status() << "\n";
        file << meta_info_();
        file << "snp-set = ";
        std::size_t pos_next_item{1};
        for (SNP snp : snp_set_) {
            file << snp;
            if (pos_next_item++ < snp_set_.size()) {
                file << ",";
            }
        }
        file << "\n";
        file << "[parameters]\n";
        file << "num-features = " << interaction_model_.rows() << "\n";
        file << "num-parameters-per-feature = " << interaction_model_.cols() << "\n";
        for (std::size_t feature_id{0}; feature_id < static_cast<std::size_t>(interaction_model_.rows()); feature_id++) {
            file << feature_id << " = ";
            for (std::size_t param_id{0}; param_id < static_cast<std::size_t>(interaction_model_.cols()); param_id++) {
                file << interaction_model_(feature_id, param_id);
                if (param_id < static_cast<std::size_t>(interaction_model_.cols()) - 1) {
                    file << ",";
                }
            }
            file << "\n";
        }
        file.close();
    }

    template<class PhenoType>
    void
    RegressionModel<PhenoType>::
    load_(const std::string & filename) {
        boost::property_tree::ptree pt;
        boost::property_tree::ini_parser::read_ini(filename, pt);
        if (pt.get<std::string>("meta.model-type") != "RegressionModel") {
            throw Error(std::string("Unexpected model type ") + pt.get<std::string>("meta.model-type") + ".Expected: RegressionModel.");
        }
        std::vector<std::string> tokens;
        misc::tokenize(pt.get<std::string>("meta.snp-set"), ',', '\'', true, tokens);
        snp_set_.clear();
        for (const auto & snp_as_string : tokens) {
            snp_set_.emplace_back(std::stoul(snp_as_string));
        }
        interaction_model_.resize(pt.get<std::size_t>("parameters.num-features"), pt.get<std::size_t>("parameters.num-parameters-per-feature"));
        for (std::size_t feature_id{0}; feature_id < static_cast<std::size_t>(interaction_model_.rows()); feature_id++) {
            misc::tokenize(pt.get<std::string>(std::string("parameters.") + std::to_string(feature_id)), ',', '\'', true, tokens);
            for (std::size_t param_id{0}; param_id < static_cast<std::size_t>(interaction_model_.cols()); param_id++) {
                interaction_model_(feature_id, param_id) = std::stod(tokens.at(param_id));
            }
        }
        check_loaded_model_(pt.get<std::string>("meta.model-cov") == "true");
    }

    template<class PhenoType>
    double
    RegressionModel<PhenoType>::
    evaluate_(const std::vector<SNP> & snp_set, bool prepare_prediction) {
        if (this->get_cov_status() and not this->instance_->has_cov()){
            throw epi::Error("No Covariates specified.");
        }
        // Compute interaction model.
        // for CLG scores: here the second (bigger) model is created, e.g., CLG-L-LC here, LC (linear with covariates)
        Eigen::MatrixXd interaction_features;
        Eigen::MatrixXd interaction_model;
        Eigen::MatrixXd interaction_hypothesis;
        double interaction_sigma{0};

        construct_feature_matrix_(snp_set, true, this->get_cov_status(), interaction_features);
        fit_model_(interaction_features, interaction_model, interaction_hypothesis, interaction_sigma);

        // Compute additive model if required.
        Eigen::MatrixXd alt_features;
        Eigen::MatrixXd alt_model;
        Eigen::MatrixXd alt_hypothesis;
        double alt_sigma{0};
        Eigen::MatrixXd prop_features;
        Eigen::MatrixXd prop_model;
        Eigen::MatrixXd prop_hypothesis;
        double prop_sigma{0};
        if (gain_score_()) {
            construct_feature_matrix_(snp_set, false, this->get_cov_status(), alt_features);
            fit_model_(alt_features, alt_model, alt_hypothesis, alt_sigma);
        } else if (score_ == "CLG-Q-QC" or score_ == "CLG-Q-LC" or score_ == "CLG-L-LC") {
            if (this->get_cov_status()){
                if (score_ == "CLG-Q-QC") {
                    construct_feature_matrix_(snp_set, true, false, alt_features);
                    fit_model_(alt_features, alt_model, alt_hypothesis, alt_sigma);
                } else if (score_ == "CLG-Q-LC" or score_ == "CLG-L-LC") {
                    // Q or L model
                    construct_feature_matrix_(snp_set, score_ == "CLG-Q-LC", false, alt_features);
                    fit_model_(alt_features, alt_model, alt_hypothesis, alt_sigma);

                    //LC model
                    construct_feature_matrix_(snp_set, false, true, prop_features);
                    fit_model_(prop_features, prop_model, prop_hypothesis, prop_sigma);
                }
            } else {
                if (score_ == "CLG-Q-QC" or score_ == "CLG-Q-LC"){
                    construct_feature_matrix_(snp_set, score_ == "CLG-Q-QC", true, alt_features);
                    fit_model_(alt_features, alt_model, alt_hypothesis, alt_sigma);
                } else if (score_ == "CLG-L-LC") {
                    // L model
                    construct_feature_matrix_(snp_set, false, false, alt_features);
                    fit_model_(alt_features, alt_model, alt_hypothesis, alt_sigma);

                    //LC model
                    construct_feature_matrix_(snp_set, false, true, prop_features);
                    fit_model_(prop_features, prop_model, prop_hypothesis, prop_sigma);
                }

            }
            construct_feature_matrix_(snp_set, false, this->get_cov_status(), alt_features);
            fit_model_(alt_features, alt_model, alt_hypothesis, alt_sigma);
        }

        // Prepare prediction if required.
        if (prepare_prediction) {
            snp_set_ = snp_set;
            interaction_model_ = interaction_model;
        }

        // If the likelihood (LLH) is the selected score, compute and return it.
        if (score_ == "LLH" or score_ == "LLH-GAIN" or (this->get_cov_status() and score_ == "CLG-Q-QC")) {
            return comp_llh_or_cov(interaction_hypothesis, interaction_sigma, alt_hypothesis, alt_sigma);
        } else if ((this->get_cov_status() and (score_ == "CLG-Q-LC" or score_ == "CLG-L-LC")) or (!this->get_cov_status() and score_ == "CLG-L-LC")){
            return comp_llh_or_cov(prop_hypothesis, prop_sigma, alt_hypothesis, alt_sigma);
        } else if ((!this->get_cov_status() and (score_ == "CLG-Q-QC" or score_ == "CLG-Q-LC"))) {
            return comp_llh_or_cov(alt_hypothesis, alt_sigma, interaction_hypothesis, interaction_sigma);
        }

        // Compute the negative log-likelihood.
        double interaction_score{0.0};
        double additive_score{0.0};
        for (Ind ind{0}; ind < this->instance_->num_inds(); ind++) {
            interaction_score -= std::log(likelihood_(ind, interaction_hypothesis, interaction_sigma));
            if (gain_score_()) {
                additive_score -= std::log(likelihood_(ind, alt_hypothesis, alt_sigma));
            }
        }

        // If the negative log-likelihood (NLL) or the NLL gain is the selected score, return it.
        if (score_ == "NLL") {
            return interaction_score;
        }
        if (score_ == "NLL-GAIN") {
            return additive_score - interaction_score;
        }

        // Compute the degrees of freedom of the model, i.e., the number of estimated parameters.
        double interaction_degrees_of_freedom{static_cast<double>(interaction_features.cols() * degrees_of_freedom_per_feature_())};
        double additive_degrees_of_freedom{static_cast<double>(alt_features.cols() * degrees_of_freedom_per_feature_())};
        if (this->instance_->quantitative_phenotypes()) {
            interaction_degrees_of_freedom += 1.0;
            additive_degrees_of_freedom += 1.0;
        }

        // If the Aikake information criterion (AIC) or the AIC gain is the selected score, return it.
        if (score_ == "AIC") {
            return 2.0 * interaction_score + 2.0 * interaction_degrees_of_freedom;
        }
        if (score_ == "AIC-GAIN") {
            return 2.0 * additive_score + 2.0 * additive_degrees_of_freedom - 2.0 * interaction_score + 2.0 * interaction_degrees_of_freedom;
        }

        // If the Bayesian information criterion (BIC) or the BIC gain is the selected score, return it.
        double log_num_inds{std::log(static_cast<double>(this->instance_->num_inds()))};
        if (score_ == "BIC") {
            return 2.0 * interaction_score + log_num_inds * interaction_degrees_of_freedom;
        }
        return 2.0 * additive_score + log_num_inds * additive_degrees_of_freedom - 2.0 * interaction_score + log_num_inds * interaction_degrees_of_freedom;
    }

    template<class PhenoType>
    void
    RegressionModel<PhenoType>::
    construct_feature_matrix_(const std::vector<SNP> & snp_set, bool construct_interaction_features, bool construct_cov_features, Eigen::MatrixXd & feature_matrix) const {
        std::size_t num_features{snp_set.size() + 1};
        if (construct_cov_features) {
            num_features += this->instance_->num_covs() + 10;
        }

        if (construct_interaction_features) {
            num_features += (snp_set.size() * (snp_set.size() - 1)) / 2;
        }

        feature_matrix.resize(this->instance_->num_inds(), num_features);
        feature_matrix.col(0).setOnes();
        for (std::size_t snp_pos{0}; snp_pos < snp_set.size(); snp_pos++) {
            for (Ind ind{0}; ind < this->instance_->num_inds(); ind++) {
                feature_matrix(ind, snp_pos + 1) = this->instance_->genotype_at_snp(snp_set.at(snp_pos), ind);
            }
        }

        // we start here, where we stopped above at snp_set.size() as snp_pos
        std::size_t index{snp_set.size() + 1};

        if (construct_cov_features){
            feature_matrix.block(0, snp_set.size(), this->instance_->get_covariates().rows(), this->instance_->num_covs()) = this->instance_->get_covariates();
            index += this->instance_->num_covs();
        }

        if (construct_interaction_features) {
            Eigen::ArrayXd col_1(this->instance_->num_inds());
            Eigen::ArrayXd col_2(this->instance_->num_inds());
            for (std::size_t feature_pos_1{0}; feature_pos_1 < snp_set.size() - 1; feature_pos_1++) {
                col_1 = feature_matrix.col(feature_pos_1 + 1).array();
                for (std::size_t feature_pos_2{feature_pos_1 + 1}; feature_pos_2 < snp_set.size(); feature_pos_2++) {
                    col_2 = feature_matrix.col(feature_pos_2 + 1).array();
                    feature_matrix.col(index++) = (col_1 * col_2).matrix();
                }
            }
        }

        // Mixed Ancestry Optionalisieren -> durch parameter setzen.
        // Im EpiJSON filterbar machen???
        if (construct_cov_features) {
            // center genotype data + compute covariance matrix
            Eigen::MatrixXd centered = feature_matrix.rowwise() - feature_matrix.colwise().mean();
            Eigen::MatrixXd cov = (centered.adjoint() * centered) / double(this->instance_->num_inds() - 1);

            // compute PCs
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);
            Eigen::VectorXd eigenvalues = eig.eigenvalues();
            Eigen::MatrixXd eigenvectors = eig.eigenvectors();

            // sort PCs by eigenvalue
            std::vector<std::pair<double, int>> pc_order;
            for (int i = 0; i < eigenvalues.size(); i++) {
                pc_order.push_back(std::make_pair(eigenvalues(i), i));
            }
            std::sort(pc_order.rbegin(), pc_order.rend());

            // Select top 10 principal components to use as covariates as suggested by Hartebrodt et al. (2022).
            for (int i = 0; i < 10; i++) {
                int pc_idx = pc_order[i].second;
                Eigen::VectorXd pc = eigenvectors.col(pc_idx);
                feature_matrix.col(index++) = centered * pc;
            }
        }
    }

    template<class PhenoType>
    void
    RegressionModel<PhenoType>::
    construct_features_(const std::vector<SNP> & snp_set, Ind ind, bool construct_interaction_features, bool construct_cov_features, Eigen::MatrixXd & features) const {
        std::size_t num_features{snp_set.size() + 1};
        if (construct_interaction_features) {
            num_features += (snp_set.size() * (snp_set.size() - 1)) / 2;
        }
        if (construct_cov_features) {
            num_features += this->instance_->num_covs() + 10;
        }

        features = Eigen::VectorXd::Ones(num_features).transpose();
        features(0,0) = 1;
        for (std::size_t snp_pos{0}; snp_pos < snp_set.size(); snp_pos++) {
            features(0, snp_pos + 1) = this->instance_->genotype_at_snp(snp_set.at(snp_pos), ind);
        }

        std::size_t index{snp_set.size() + 1};
        if (construct_interaction_features) {
            double feature_1{0.0};
            double feature_2{0.0};
            for (std::size_t snp_pos_1{0}; snp_pos_1 < snp_set.size() - 1; snp_pos_1++) {
                feature_1 = features(0, snp_pos_1 + 1);
                for (std::size_t snp_pos_2{snp_pos_1 + 1}; snp_pos_2 < snp_set.size(); snp_pos_2++) {
                    feature_2 = features(0, snp_pos_2 + 1);
                    features(0, index++) = feature_1 * feature_2;
                }
            }
        }

        //Insert the covariate featues if construct_cov_features == true
        for (std::size_t cov_pos{0}; cov_pos < this->instance_->num_covs() && construct_cov_features; cov_pos++) {
            features(0, index++) = this->instance_->get_covariates_at_ind(ind)[cov_pos];
        }

    }

    template<class PhenoType>
    double
    RegressionModel<PhenoType>::
    comp_llh_or_cov(Eigen::MatrixXd hypothesis1, double sigma1, Eigen::MatrixXd hypothesis2, double sigma2) const {
        double interaction_score{1.0};
        double additive_score{1.0};
        for (Ind ind{0}; ind < this->instance_->num_inds(); ind++) {
            interaction_score *= likelihood_(ind, hypothesis1, sigma1);
            if (gain_score_() or cov_score_()) {
                additive_score *= likelihood_(ind, hypothesis2, sigma2);
            }
        }
        if (gain_score_() or cov_score_()) {
            return interaction_score - additive_score;
        }
        return interaction_score;
    }

    template<class PhenoType>
    bool
    RegressionModel<PhenoType>::
    cov_score_() const {
        return score_ == "CLG-L-LC" or score_ == "CLG-Q-QC" or score_ == "CLG-Q-LC";
    }

    template<class PhenoType>
    bool
    RegressionModel<PhenoType>::
    get_cov_status() const {
        return incl_cov_;
    }

    template<class PhenoType>
    void
    RegressionModel<PhenoType>::
    cov_activate() {
        if (this->instance_->has_cov() == false){
            incl_cov_ = false;
            throw epi::Error("No covariates specified.");
        } else {
            incl_cov_ = true;
        }
    }

    template<class PhenoType>
    void
    RegressionModel<PhenoType>::
    cov_deactivate() {
        incl_cov_ = false;
    }

    template<class PhenoType>
    bool
    RegressionModel<PhenoType>::
    get_mixed_ancestry_status() const {
        return mixed_ancestry_;
    }

    template<class PhenoType>
    void
    RegressionModel<PhenoType>::
    mixed_ancestry_activate() {
        incl_cov_ = true;
    }

    template<class PhenoType>
    void
    RegressionModel<PhenoType>::
    mixed_ancestry_deactivate() {
        incl_cov_ = false;
    }

    template<class PhenoType>
    bool
    RegressionModel<PhenoType>::
    gain_score_() const {
        return score_ == "LLH-GAIN" or score_ == "NLL-GAIN" or score_ == "AIC-GAIN" or score_ == "BIC-GAIN";
    }


    template<class PhenoType>
    const std::string &RegressionModel<PhenoType>::get_score() {
        return score_;
    }

    template<class PhenoType>
    const Eigen::MatrixXd &RegressionModel<PhenoType>::get_interaction_model_() {
        return interaction_model_;
    }

}


#endif /* SRC_MODEL_REGRESSION_MODEL_IPP_ */
