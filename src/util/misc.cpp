/*******************************************************************************
 *                                                                             *
 *   Copyright (C) 2020 by David B. Blumenthal                                 *
 *                                                                             *
 *   This file is part of NeEDL.                                        *
 *                                                                             *
 *   NeEDL is free software: you can redistribute it and/or modify it   *
 *   under the terms of the GNU General Public License as published by         *
 *   the Free Software Foundation, either version 3 of the License, or         *
 *   (at your option) any later version.                                       *
 *                                                                             *
 *   NeEDL is distributed in the hope that it will be useful,           *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the              *
 *   GNU General Public License for more details.                              *
 *                                                                             *
 *   You should have received a copy of the GNU General Public License         *
 *   along with NeEDL. If not, see <http://www.gnu.org/licenses/>.      *
 *                                                                             *
 ******************************************************************************/


/*!
 * @file  misc.ipp
 * @brief Definitions of miscellaneous helper functions.
 */

#ifndef SRC_UTIL_MISC_IPP_
#define SRC_UTIL_MISC_IPP_

#include "misc.hpp"

namespace epi {

template<class T>
std::ostream &
operator<<(std::ostream & os, const std::vector<T> & vector) {
	os << "[";
	for (std::size_t i{0}; i < vector.size(); ++i) {
		os << i << ":" << vector[i];
		if (i != vector.size() - 1) {
			os << ", ";
		}
	}
	os << "]\n";
	return os;
}

namespace misc {

std::size_t
tokenize(const std::string & sentence, char sep, char quote, bool trim, std::vector<std::string> & tokens, bool clear_tokens, bool skip_empty_tokens) {
	if (clear_tokens) {
		tokens.clear();
	}
	std::size_t num_tokens{0};
	if (sep == quote) {
		throw Error("Separator and quotation mark are identical.");
	}
	bool outside_quotes{true};
	std::size_t word_length{0};
	std::size_t pos_word_start{0};
	std::string token;
	for (std::size_t pos{0}; pos < sentence.size(); pos++) {
		if (sentence.at(pos) == quote) {
			if (not outside_quotes and pos < sentence.size() - 1) {
				if (sentence.at(pos + 1) != sep) {
					throw Error("Sentence contains closing quote which is followed by a char different from " + std::to_string(sep) + ".");
				}
			}
			word_length++;
			outside_quotes = not outside_quotes;
		}
		else if (outside_quotes and sentence.at(pos) == sep) {
			if (word_length > 0) {
				token = sentence.substr(pos_word_start, word_length);
				if (trim) {
					trim_both(token);
				}
				tokens.push_back(token);
				num_tokens++;
			} else if (!skip_empty_tokens) {
                tokens.emplace_back("");
                num_tokens++;
            }
			pos_word_start = pos + 1;
			word_length = 0;
		}
		else {
			word_length++;
		}
	}
	if (not outside_quotes) {
		throw Error("Sentence contains unbalanced quotes.");
	}
	if (word_length > 0) {
		token = sentence.substr(pos_word_start, word_length);
		if (trim) {
			trim_both(token);
		}
		tokens.push_back(token);
		num_tokens++;
	} else if (!skip_empty_tokens) {
        tokens.emplace_back("");
        num_tokens++;
    }
	return num_tokens;
}

void
trim_left(std::string & token, char rem) {
	token.erase(token.begin(), std::find_if(token.begin(), token.end(), [rem] (char c) {return c != rem;}));
}

void
trim_right(std::string & token, char rem) {
	token.erase(std::find_if(token.rbegin(), token.rend(), [rem] (char c) {return c != rem;}).base(), token.end());
}

void
trim_both(std::string & token, char rem) {
	trim_left(token);
	trim_right(token);
}

std::size_t
genotype_to_id(const std::vector<GenoType> & genotype) {
	std::size_t genotype_id{0};
	std::size_t exponent{genotype.size() - 1};
	for (GenoType geno_at_snp : genotype) {
		genotype_id += geno_at_snp * uint_pow(3, exponent--);
	}
	return genotype_id;
}

std::size_t
uint_pow(std::size_t base, std::size_t exponent) {
	if (exponent == 0) {
		return 1;
	}
	if (exponent == 1) {
		return base;
	}
	std::size_t temp{uint_pow(base, exponent/2)};
	if (exponent % 2 == 0) {
		return temp * temp;
	}
	return base * temp * temp;
}

std::size_t
size_penetrance_table(std::size_t size_snp_set) {
	return uint_pow(3, size_snp_set);
}

void
options_string_to_options_map(const std::string & options_string, std::map<std::string, std::string> & options_map) {
	if (options_string == "") return;
	options_map.clear();
	std::vector<std::string> words;
	tokenize(options_string, ' ', '\'', true, words);
	std::string option_name;
	bool expect_option_name{true};
	for (auto word : words) {
		if (expect_option_name) {
			if (is_option_name(word)) {
				option_name = word;
				if (options_map.find(option_name) != options_map.end()) {
					throw Error("Multiple specification of option \"" + option_name + "\".");
				}
				options_map[option_name] = "";
			}
			else {
				throw Error("Invalid options \"" + options_string + "\". Usage: options = \"[--<option> <arg>] [...]\"");
			}
		}
		else {
			if (is_option_name(word)) {
				throw Error("Invalid options \"" + options_string + "\". Usage: options = \"[--<option> <arg>] [...]\"");
			}
			else {
				options_map[option_name] = word;
			}
		}
		expect_option_name = not expect_option_name;
	}
}

bool
is_option_name(std::string & word) {
	if (word.at(0) == '\'') {
		word = word.substr(1, word.size() - 2);
		return false;
	}
	if (word.size() < 3) {
		return false;
	}
	if ((word.at(0) == '-') and (word.at(1) == '-') and (word.at(2) != '-')) {
		word = word.substr(2);
		return true;
	}
	return false;
}

void
id_to_genotype(std::size_t genotype_id, std::size_t genotype_size, std::vector<GenoType> & genotype) {
	if (genotype_id >= static_cast<std::size_t>(std::pow(3, genotype_size))) {
		throw Error("Invalid genotype ID for genotype size " + std::to_string(genotype_size) + ". Expected: value between 0 and " + std::to_string(std::pow(3, genotype_size) - 1) + ".");
	}
	genotype.clear();
	while (genotype_id > 0) {
		genotype.emplace_back(static_cast<GenoType>(genotype_id % 3));
		genotype_id /= 3;
	}
	while (genotype.size() < genotype_size) {
		genotype.emplace_back(0);
	}
	std::reverse(genotype.begin(), genotype.end());
}

double normal_pdf(double x, double mu, double sigma) {
	boost::math::normal_distribution<> dist(mu, sigma);
	return ensure_valid_continuous_probability(boost::math::pdf(dist, x));
}

double normal_cdf(double x, double mu, double sigma) {
	boost::math::normal_distribution<> dist(mu, sigma);
	return ensure_valid_continuous_probability(boost::math::cdf(dist, x));
}

double f_pdf(double x, std::size_t numerator_degrees_of_freedom, std::size_t denominator_degrees_of_freedom) {
	boost::math::fisher_f_distribution<> dist(static_cast<double>(numerator_degrees_of_freedom), static_cast<double>(denominator_degrees_of_freedom));
	return ensure_valid_continuous_probability(boost::math::pdf(dist, x));
}

double f_cdf(double x, std::size_t numerator_degrees_of_freedom, std::size_t denominator_degrees_of_freedom) {
	boost::math::fisher_f_distribution<> dist(static_cast<double>(numerator_degrees_of_freedom), static_cast<double>(denominator_degrees_of_freedom));
	return ensure_valid_continuous_probability(boost::math::cdf(dist, x));
}

double chi_squared_pdf(double x, std::size_t degrees_of_freedom) {
	boost::math::chi_squared_distribution<> dist(static_cast<double>(degrees_of_freedom));
	return ensure_valid_continuous_probability(boost::math::pdf(dist, x));
}

double chi_squared_cdf(double x, std::size_t degrees_of_freedom) {
	boost::math::chi_squared_distribution<> dist(static_cast<double>(degrees_of_freedom));
	return ensure_valid_continuous_probability(boost::math::cdf(dist, x));
}

double students_t_pdf(double x, std::size_t degrees_of_freedom) {
	boost::math::students_t dist(static_cast<double>(degrees_of_freedom));
	return ensure_valid_continuous_probability(boost::math::pdf(dist, x));
}

double students_t_cdf(double x, std::size_t degrees_of_freedom) {
	boost::math::students_t dist(static_cast<double>(degrees_of_freedom));
	return ensure_valid_continuous_probability(boost::math::cdf(dist, x));
}

double one_sample_t_test(const std::vector<double> & sample, double expected_mean, options::TestDirection direction) {

	// Throw error if the samples contains less than two elements.
	if (sample.size() < 2) {
		throw Error(std::string("Cannot carry out one sample t-test for a sample of size ") + std::to_string(sample.size()) + ".");
	}

	// Compute sample mean and standard deviation.
	double mu{0};
	for (double value : sample) {
		mu += value;
	}
	mu /= static_cast<double>(sample.size());
	double sigma{0.0};
	for (double value : sample) {
		sigma += (value - mu) * (value - mu);
	}
	sigma /= static_cast<double>(sample.size() - 1);
	sigma = std::sqrt(sigma);

	// Compute test statistic and return appropriate p-values if the standard deviation is 0.
	double t_stat{(mu - expected_mean) * std::sqrt(static_cast<double>(sample.size()))};
	if (sigma > 0.0) {
		t_stat /= sigma;
	}
	else if (t_stat < 0.0) {
		if (direction == options::TestDirection::UPPER_TAILED) {
			return 1.0;
		}
		else {
			return 0.0;
		}
	}
	else if (t_stat > 0.0) {
		if (direction == options::TestDirection::LOWER_TAILED) {
			return 1.0;
		}
		else {
			return 0.0;
		}
	}

	// Compute and return p-values.
	double return_value{0};

	switch (direction) {
	case options::TestDirection::TWO_TAILED:
		return_value = 2.0 * (1.0 - students_t_cdf(std::fabs(t_stat), sample.size() - 1));
		break;
	case options::TestDirection::LOWER_TAILED:
		return_value = students_t_cdf(t_stat, sample.size() - 1);
		break;
	case options::TestDirection::UPPER_TAILED:
		return_value = 1.0 - students_t_cdf(t_stat, sample.size() - 1);
		break;
	}
	return ensure_valid_continuous_probability(return_value);
}

double ensure_valid_continuous_probability(double possibly_invalid_continuous_probability) {
	if (possibly_invalid_continuous_probability <= 0) {
		return std::numeric_limits<double>::min();
	}
	if (possibly_invalid_continuous_probability >= 1) {
		return 1 - std::numeric_limits<double>::min();
	}
	return possibly_invalid_continuous_probability;
}

    std::string inf_split(std::string txt_information, char seperator, int place)
    {
        int count = std::count(txt_information.begin(), txt_information.end(), seperator);
        if (place > count)
        {
            std::string output = "False location - place must be < ";
            output.append(std::to_string(count));
            output.append("! Place is 0-based!");
            return output;
        }

        int start_pos = 0;
        int end_pos = 0;

        for (int i = 0; i < count; i++)
        {
            end_pos = txt_information.find(seperator);
            if (i == place)
            {
                int pos_diff = end_pos - start_pos;
                return txt_information.substr(start_pos, pos_diff);
            }
            else
            {
                txt_information = txt_information.substr(end_pos + 1);

            }

        }

        return txt_information;
    }

}

}




#endif /* SRC_UTIL_MISC_IPP_ */
