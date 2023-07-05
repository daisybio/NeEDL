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
 * @file types.hpp
 * @brief Declarations and inclusions of types used by various classes.
 */

#ifndef SRC_UTIL_TYPES_HPP_
#define SRC_UTIL_TYPES_HPP_

// Include standard libraries.

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cctype>
#include <cassert>
#include <map>
#include <set>
#include <stdexcept>
#include <algorithm>
#include <type_traits>
#include <cmath>
#include <random>
#include <functional>
#include <dirent.h>
#include <fstream>
#include <cstdlib>
#include <strings.h>
#include <fstream>

// Include Boost Graph Libraries
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/named_function_params.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/config.hpp>


// Include Boost libraries.
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>
// Include OpenMP.
#ifdef _OPENMP
#include <omp.h>
#endif

// Include Eigen.
#include <Eigen/Dense>
#include "Logger.hpp"

/*!
 * @brief Global namespace for GenEpiSeeker.
 */
namespace epi {


/*!
 * @brief Type of SNPs.
 * This currently MUST be a 32 bit integer to use the fast edge unique check implemented in class SNPNetwork
 */
    typedef uint_fast32_t SNP;


/*!
 * @brief Type of SNPs.
 * This currently MUST be a 32 bit integer to use the fast edge unique check implemented in class SNPNetwork
 */
    typedef uint_fast32_t Cov;

/*!
 * @brief Type of individuals.
 */
    typedef std::size_t Ind;

/*!
 * @brief Type of entries in the genotype matrix.
 */
    typedef std::uint_fast8_t GenoType;

/*!
 * @brief Use this type if you want to use GenEpiSeeker for quantitative phenotypes.
 */
    typedef double QuantitativePhenoType;

/*!
 * @brief Use this type if you want to use GenEpiSeeker for categorical phenotypes.
 */
    typedef std::uint_fast8_t CategoricalPhenoType;

/*!
 * @brief Type for reporting runtimes in seconds.
 */
    typedef std::chrono::duration<double> Seconds;

// Constant expressions.
/*!
 * @brief Used to denote undefined unsigned integers.
 * @return An unsigned integer used to denote undefined integers.
 */
    constexpr std::size_t undefined_uint() {return std::numeric_limits<std::size_t>::max();}

/*!
 * @brief Used to denote undefined doubles.
 * @return A double used to denote undefined doubles.
 */
    constexpr double undefined_double() {return std::numeric_limits<double>::infinity();}

/*!
 * @brief The constant @f$\pi@f$.
 * @return The constant @f$\pi@f$.
 */
    constexpr double pi() {return 3.14159265358979323846;}

// Error class.
/*!
 * @brief Exceptions thrown by GenEpiSeeker are of this type.
 */
    struct Error : public std::runtime_error {
        /*!
         * @brief Constructor.
         * @param[in] message Error message
         */
        explicit Error(const std::string & message) : std::runtime_error(message) {}
    };

    struct SNPNotFoundError : public Error {
        SNPNotFoundError(const std::string & message, const std::string & name_) : Error(message) {
            name = name_;
        }

    public:
        std::string get_name() const {
            return name;
        }
    private:
        std::string name;
    };

/*!
 * @brief Contains enum classes used to specify options.
 */
    namespace options {

/*!
 * @brief Specifies the SNP type.
 */
        enum class SNPType {
            NON_CODING,       //!< Non-coding SNPs.
            CODING_SYNONYMOUS,//!< Coding, synonymous SNPs.
            CODING_MISSENSE,  //!< Coding, missense SNPs.
            CODING_NONSENSE   //!< Coding, nonsense SNPs.
        };

/*!
 * @brief Specifies the input format.
 */
        enum class InputFormat {
            JSON_EPIGEN,              //!< Input is JSON file in the format produced by EpiGEN.
            NEEDL_BIN,                //!< A binary file that can be created with NeEDL.
            CSV_SNPS_AS_ROWS_FIRST,   //!< Input is CSV file where the SNPs are represented by the rows and the first column contains SNP information.
            CSV_SNPS_AS_ROWS_LAST,    //!< Input is CSV file where the SNPs are represented by the rows and the last column contains SNP information.
            CSV_SNPS_AS_COLUMNS_FIRST,//!< Input is CSV file where the SNPs are represented by the columns and the first row contains SNP information.
            CSV_SNPS_AS_COLUMNS_LAST,  //!< Input is CSV file where the SNPs are represented by the columns and the last row contains SNP information.
            CSV_COV
        };

/*!
 * @brief Selects the epistasis model.
 */
        enum class EpistasisModel {
            BAYESIAN_MODEL,  //!< Selects epi::BayesianModel.
            PENETRANCE_MODEL,//!< Selects epi::PenetranceModel.
            REGRESSION_MODEL,//!< Selects epi::RegressionModel.
            VARIANCE_MODEL   //!< Selects epi::VarianceModel.
        };

        enum class EpistasisScore {
            BAYESIAN,
            BAYESIAN_COV,
            VARIANCE,
            VARIANCE_COV,
            PENETRANCE_NLL,
            PENETRANCE_LLH,
            PENETRANCE_AIC,
            PENETRANCE_BIC,
            PENETRANCE_COV_NLL,
            PENETRANCE_COV_LLH,
            PENETRANCE_COV_AIC,
            PENETRANCE_COV_BIC,
            REGRESSION_NLL,
            REGRESSION_LLH,
            REGRESSION_AIC,
            REGRESSION_BIC,
            REGRESSION_NLL_GAIN,
            REGRESSION_LLH_GAIN,
            REGRESSION_AIC_GAIN,
            REGRESSION_BIC_GAIN,
            REGRESSION_COV_NLL,
            REGRESSION_COV_LLH,
            REGRESSION_COV_AIC,
            REGRESSION_COV_BIC,
            REGRESSION_CLG_Q_LC,
            REGRESSION_CLG_Q_QC,
            REGRESSION_CLG_L_LC

        };
#define NUM_EPISTASIS_SCORES 27
        EpistasisModel epistasis_model_from_string(const std::string & model_string);
        EpistasisModel epistasis_model_from_epistasis_score(EpistasisScore score);
        EpistasisScore epistasis_score_from_string(const std::string & score_string);
        std::vector<std::string> get_all_epistasis_scores();
        std::string epistasis_model_to_string(EpistasisModel model);
        std::string epistasis_score_to_string(EpistasisScore model);
        bool is_epistasis_model_string(const std::string & model_string);

/*!
 * @brief Specifies the sense of epistasis models a.k.a. objective functions.
 */
        enum class ModelSense {
            MINIMIZE,//!< The objective should be minimized.
            MAXIMIZE //!< The objective should be maximized.
        };

/*!
 * @brief Specified whether the phenotypes are quantitative or categorical.
 */
        enum class PhenoType {
            QUANTITATIVE,//!< Quantitative phenotypes.
            CATEGORICAL  //!< Categorical phenotypes.
        };

/*!
 * @brief Specifies whether the data loaded into the instance should be used for training or validation.
 */
        enum class DataPurpose {
            TRAINING, //!< The data should be used for training.
            VALIDATION//!< The data should be used for validation.
        };

/*!
 * @brief Specifies direction of statistical test.
 */
        enum class TestDirection {
            TWO_TAILED,  //!< Two tailed test, i.e. p-value = 2 * Pr(D > |test_statistic|), with D distributed according null hypothesis.
            LOWER_TAILED,//!< Lower tailed test, i.e. p-value = Pr(D < test_statistic), with D distributed according null hypothesis
            UPPER_TAILED //!< Upper tailed test, i.e. p-value = Pr(D > test_statistic), with D distributed according null hypothesis
        };

    }



}

/*!
 * @brief Streams epi::Options::EpistasisModel.
 * @param[in] os Output stream.
 * @param[in] model Epistasis model selector.
 * @return Output stream.
 */
std::ostream & operator<<(std::ostream & os, const epi::options::EpistasisModel & model);
std::ostream & operator<<(std::ostream & os, const epi::options::EpistasisScore & model);


#ifdef HEADER_ONLY
#include "types.cpp"
#endif

#endif /* SRC_UTIL_TYPES_HPP_ */
