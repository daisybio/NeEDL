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
 * @file  misc.hpp
 * @brief Declarations of miscellaneous helper functions.
 */

#ifndef SRC_UTIL_MISC_HPP_
#define SRC_UTIL_MISC_HPP_

#include "types.hpp"

namespace epi {

/*!
 * @brief Prints a vector to output stream.
 * @param[in,out] os Output stream.
 * @param vector The vector that should be printed.
 * @return Output stream.
 */
template<class T>
std::ostream& operator<<(std::ostream & os, const std::vector<T> & vector);

/*!
 * @brief Contains miscellaneous helper functions.
 */
namespace misc {

/*!
 * @brief Separates a sentence into tokens separated by @p sep (unless contained in quotes marked by @p quote).
 * @param[in] sentence The sentence that should be tokenized.
 * @param[in] sep The separator. Must be different from @p quote.
 * @param[in] quote The quotation mark. Must be different from @p sep.
 * @param[in] trim If @p true, leading and trailing white-spaces are removed from the tokens.
 * @param[out] tokens The obtained tokens are appended to this vector.
 * @param[in] clear_tokens If @p true, the @p tokens is cleared before reading the sentence.
 * @return The number of tokens in the sentence.
 */
std::size_t tokenize(const std::string & sentence, char sep, char quote, bool trim, std::vector<std::string> & tokens, bool clear_tokens = true, bool skip_empty_tokens = true);

/*!
 * @brief Remove leading characters from a string.
 * @param[in,out] token The string from which should be trimmed.
 * @param[in] rem The character that should be removed.
 */
void trim_left(std::string & token, char rem = ' ');

/*!
 * @brief Remove trailing characters from a string.
 * @param[in,out] token The string from which should be trimmed.
 * @param[in] rem The character that should be removed.
 */
void trim_right(std::string & token, char rem = ' ');

/*!
 * @brief Remove leading and trailing characters from a string.
 * @param[in,out] token The string from which should be trimmed.
 * @param[in] rem The character that should be removed.
 */
void trim_both(std::string & token, char rem = ' ');

/*!
 * @brief Computes integer ID of a genotype.
 * @details The integer ID of genotype @f$(g_s)^{n-1}_{s=0}@f$ is defined as @f$\sum^{n-1}_{s=0}g_s\cdot 3^{n-1-s}@f$.
 * For instance, the genotype @f$(1,0,0,2)@f$ yields the ID @f$1\cdot3^3 + 0 \cdot 3^2 + 0 \cdot 3^1 + 2 \cdot 3^0 = 29@f$.
 * @param[in] genotype The genotype whose ID should be computed.
 * @return The integer ID of @p genotype from the range @f$\{0,\ldots, 3^{n} - 1\}@f$, where @f$n@f$ is the size of the genotype.
 */
std::size_t genotype_to_id(const std::vector<GenoType> & genotype);

/*!
 * @brief Power function for unsigned integers.
 * @param[in] base The base.
 * @param[in] exponent The exponent.
 * @return @f$\mathtt{base}^\mathtt{exponent}@f$
 */
std::size_t uint_pow(std::size_t base, std::size_t exponent);

/*!
 * @brief Retrieves the genotype from a genotype ID.
 * @param[in] genotype_id The genotype ID for which the genotype should be retrieved.
 * @param[in] genotype_size The target size of the retrieved genotype.
 * @param[out] genotype The retrieved genotype.
 */
void id_to_genotype(std::size_t genotype_id, std::size_t genotype_size, std::vector<GenoType> & genotype);

/*!
 * @brief Returns size of penetrance table for given SNP set size.
 * @param[in] size_snp_set SNP set size.
 * @return The size of the penetrance table for a SNP set of size @p snp_set_size.
 */
std::size_t size_penetrance_table(std::size_t size_snp_set);

/*!
 * @brief Transforms an options string into an options map.
 * @param[in] options_string Options string of the form <tt>[\--@<option@> @<arg@>] [...]</tt>.
 * @param[out] options_map Map with one key-value pair (<tt>@<option@></tt>, <tt>@<arg@></tt>) for each option contained in the string.
 */
void options_string_to_options_map(const std::string & options_string, std::map<std::string, std::string> & options_map);

/*!
 * @brief Checks whether a word is an option name and, if so, removes the leading dashes.
 * @param[in,out] word
 * @return True if @p word is of the form "--<option>".
 */
bool is_option_name(std::string & word);

/*!
 * @brief Probability density function (PDF) of normal distribution.
 * @param[in] x Value at which the PDF should be evaluated.
 * @param[in] mu Expectation of normal distribution.
 * @param[in] sigma Standard deviation of normal distribution.
 * @return PDF of N(@p mu, @p sigma) at @p x.
 */
double normal_pdf(double x, double mu, double sigma);

/*!
 * @brief Cumulative density function (CDF) of normal distribution.
 * @param[in] x Value at which the CDF should be evaluated.
 * @param[in] mu Expectation of normal distribution.
 * @param[in] sigma Standard deviation of normal distribution.
 * @return CDF of N(@p mu, @p sigma) at @p x.
 */
double normal_cdf(double x, double mu, double sigma);

/*!
 * @brief Probability density function (PDF) of F distribution.
 * @param[in] x Value at which the PDF should be evaluated.
 * @param[in] numerator_degrees_of_freedom Degrees of freedom of the numerator.
 * @param[in] denominator_degrees_of_freedom Degrees of freedom of the denominator.
 * @return PDF of F(@p numerator_degrees_of_freedom, @p denominator_degrees_of_freedom) at @p x.
 */
double f_pdf(double x, std::size_t numerator_degrees_of_freedom, std::size_t denominator_degrees_of_freedom);

/*!
 * @brief Cumulative density function (CDF) of F distribution.
 * @param[in] x Value at which the CDF should be evaluated.
 * @param[in] numerator_degrees_of_freedom Degrees of freedom of the numerator.
 * @param[in] denominator_degrees_of_freedom Degrees of freedom of the denominator.
 * @return CDF of F(@p numerator_degrees_of_freedom, @p denominator_degrees_of_freedom) at @p x.
 */
double f_cdf(double x, std::size_t numerator_degrees_of_freedom, std::size_t denominator_degrees_of_freedom);

/*!
 * @brief Probability density function (PDF) of @f$\chi^2@f$-distribution.
 * @param[in] x Value at which the PDF should be evaluated.
 * @param[in] degrees_of_freedom Degrees of freedom.
 * @return PDF of @f$\chi^2@f$(@p degrees_of_freedom) at @p x.
 */
double chi_squared_pdf(double x, std::size_t degrees_of_freedom);

/*!
 * @brief Cumulative density function (CDF) of @f$\chi^2@f$-distribution.
 * @param[in] x Value at which the CDF should be evaluated.
 * @param[in] degrees_of_freedom Degrees of freedom.
 * @return CDF of @f$\chi^2@f$(@p degrees_of_freedom) at @p x.
 */
double chi_squared_cdf(double x, std::size_t degrees_of_freedom);

/*!
 * @brief Probability density function (PDF) of Student's t-distribution.
 * @param[in] x Value at which the PDF should be evaluated.
 * @param[in] degrees_of_freedom Degrees of freedom.
 * @return PDF of Student's t-distribution at @p x.
 */
double students_t_pdf(double x, std::size_t degrees_of_freedom);

/*!
 * @brief Cumulative density function (CDF) of Student's t-distribution
 * @param[in] x Value at which the CDF should be evaluated.
 * @param[in] degrees_of_freedom Degrees of freedom.
 * @return CDF of Student's t-distribution at @p x.
 */
double students_t_cdf(double x, std::size_t degrees_of_freedom);

/*!
 * @brief Computes p-value of one sample t-test.
 * @details
 * - Alternative hypothesis if @p direction == epi::options::TestDirection::TWO_TAILED: the true mean of
 *   of the distribution from which @p sample is drawn is different from @p expected_mean.
 * - Alternative hypothesis if @p direction == epi::options::TestDirection::LOWER_TAILED: the true mean of
 *   of the distribution from which @p sample is drawn  is smaller than @p expected_mean.
 * - Alternative hypothesis if @p direction == epi::options::TestDirection::UPPER_TAILED: the true mean of
 *   of the distribution from which @p sample is drawn  is greater than @p expected_mean.
 * @param[in] sample The test sample.
 * @param[in] expected_mean The expected mean.
 * @param[in] direction The direction of the test.
 * @return The p-value.
 */
double one_sample_t_test(const std::vector<double> & sample, double expected_mean, options::TestDirection direction);

/*!
 * @brief Ensures that the results returned by *_pdf() and *_cdf() are indeed values which are strictly greater than 0 and smaller than 1.
 * @param[in] possibly_invalid_continuous_probability The possibly invalid probability computed by boost::math::pdf() or boost::math::cdf().
 * @return If @p possibly_invalid_continuous_probability is <= 0 or >= 1, the function returns, respectively, std::numeric_limits<double>::min() or 1 - std::numeric_limits<double>::min().
 */
double ensure_valid_continuous_probability(double possibly_invalid_continuous_probability);

}

}

/*!
 * @brief returns the string at a certrain position, in a string seperated by seperator
 * @param[in] txt_information String to be splitted
 * @param[in] seperator delimiter to seperate string
 * @param[in] place which position should be returned, 0 based, must be smaller than number of seperators
 * @return returns item of txt_information on place
 */
std::string inf_split(std::string txt_information, char seperator, int place);


#ifdef HEADER_ONLY
#include "misc.cpp"
#endif

#endif /* SRC_UTIL_MISC_HPP_ */
