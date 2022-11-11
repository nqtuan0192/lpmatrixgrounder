#include <algorithm>
#include <array>
#include <bitset>
#include <boost/bimap.hpp>
#include <boost/dynamic_bitset.hpp>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <unordered_set>
#include <vector>

// #include <clingo.hh>

#include "lpmatrixgrounder.h"

template <class MapType>
void print_map(const MapType& map, const std::string& separator,
               std::ostream& os) {
  typedef typename MapType::const_iterator const_iterator;

  for (const_iterator i = map.begin(), iend = map.end(); i != iend; ++i) {
    for (auto it : i->first) {
      os << it << " ";
    }
    os << separator << i->second << std::endl;
  }
}

void standardizeRules(LPI_Format& rules, uint32_t& atoms_count,
                      LPI_Format& lpi_format_standardized) {
  // sort input rules by head atom
  std::sort(rules.begin(), rules.end(),
            [](const std::vector<int>& a, const std::vector<int>& b) {
              return a[LPI_IDX::RULE_HEAD] < b[LPI_IDX::RULE_HEAD];
            });
  // sort body of each rule
  std::for_each(rules.begin(), rules.end(), [](LPI_Rule& r) {
    std::sort(
        r.begin() + LPI_IDX::RULE_BODY, r.end(),
        [](const int& a, const int& b) { return std::abs(a) < std::abs(b); });
  });

  typedef boost::bimap<LPI_Rule, int> body_rule_bimap;
  typedef body_rule_bimap::value_type body_rule;
  body_rule_bimap duplicate_body_mapping;

  // standardize steps
  // - identify rules with the same head
  // - introduce new atoms for same head rules that the body length is larger
  // than 1
  // - merge rules with the same head and
  // std::map<int, int> atoms_mapping;
  auto idx = rules.begin();
  while (idx != rules.end()) {
    // count number of rules with the same head atom
    int count = 0;
    for (auto it = idx + 1; it != rules.end(); ++it) {
      if ((*idx)[LPI_IDX::RULE_HEAD] != (*it)[LPI_IDX::RULE_HEAD]) {
        break;
      }
      ++count;
    }
    if (count > 0) {  // if there are rules with the same head atom
      LPI_Rule new_or_rule = {LPI_RULE_TYPE::OR_RULE,
                              (*idx)[LPI_IDX::RULE_HEAD]};
      for (auto it = idx; it != idx + count + 1; ++it) {
        if (it->size() >
            LPI_IDX::RULE_BODY + 1) {  // body length is larger than 1
          LPI_Rule new_and_rule = {LPI_RULE_TYPE::AND_RULE,
                                   (*idx)[LPI_IDX::RULE_HEAD]};
          new_and_rule.insert(new_and_rule.end(),
                              it->begin() + LPI_IDX::RULE_BODY, it->end());

          int head_atom;
          auto body_it = (it + LPI_IDX::RULE_BODY);
          body_rule_bimap::left_iterator found_it =
              duplicate_body_mapping.left.find(*body_it);
          if (found_it != duplicate_body_mapping.left.end()) {
            head_atom = duplicate_body_mapping.left.at(*body_it);
          } else {
            head_atom = atoms_count + 1;
            duplicate_body_mapping.insert(body_rule(*body_it, head_atom));
            ++atoms_count;
          }
          new_and_rule[LPI_IDX::RULE_HEAD] = head_atom;
          lpi_format_standardized.emplace_back(new_and_rule);
          new_or_rule.emplace_back(head_atom);
        } else {
          new_or_rule.emplace_back((*it)[LPI_IDX::RULE_BODY]);
        }
      }
      lpi_format_standardized.emplace_back(new_or_rule);
    } else {
      lpi_format_standardized.emplace_back(*idx);
    }
    idx = idx + count + 1;
  }
}

#define N 32
// typedef std::bitset<N> BitVector;
typedef boost::dynamic_bitset<> BitVector;
typedef std::vector<BitVector> BitMatrix;

void printBitMatrix(BitMatrix mat) {
  for (auto v : mat) {
    std::cout << v << std::endl;
  }
}

void printLPBitMatrix(BitMatrix mat, BitVector ruletype, uint32_t atoms_count) {
  for (auto it = mat.begin(); it != mat.end(); ++it) {
    std::cout << *it
              << "\t Type: " << (ruletype[it - mat.begin()] ? " OR" : "AND")
              << "\t";
    std::cout << it - mat.begin() + 1 << " <= ";
    for (uint32_t bit_idx = 0; bit_idx < atoms_count; ++bit_idx) {
      if ((*it)[bit_idx]) {
        std::cout << bit_idx + 1 << ", ";
      }
    }
    for (uint32_t bit_idx = atoms_count; bit_idx < it->size(); ++bit_idx) {
      if ((*it)[bit_idx]) {
        std::cout << (int)-(bit_idx - atoms_count + 1) << ", ";
      }
    }
    std::cout << std::endl;
  }
}

int main(int argc, char const** argv) {
  ClingoGrounder cgrounder = ClingoGrounder();
  cgrounder.input_files(
      {"../samples/model_graphcoloring.lp", "../samples/graph_6nodes.lp"});
  // cgrounder.input_files({"../samples/example.lp"});
  cgrounder.ground();
  for (auto statement : cgrounder.obs_.rule_statements_) {
    std::cout << statement << std::endl;
  }
  for (auto statement : cgrounder.obs_.output_statements_) {
    std::cout << statement << std::endl;
  }

  LPI_Format rules;
  uint32_t atoms_count = 0;
  std::unordered_set<int> atoms_neg;
  for (ASPIF_RuleStatement<int> statement : cgrounder.obs_.rule_statements_) {
    if (!statement.head.rule_type) {
      std::cout << "Disjunction rule: " << statement << std::endl;

      std::cout << "------------head: ";
      if (statement.head.elements.size() == 0) {
        std::cout << "None"
                  << "===>>> constraint" << std::endl;
      } else {
        for (uint32_t i = 0; i < statement.head.elements.size(); ++i) {
          std::cout << statement.head.elements[i] << " ";
        }
        std::cout << std::endl;
      }

      std::cout << "------------body: ";
      if (statement.body.elements.size() == 0) {
        std::cout << "None"
                  << "===>>> fact" << std::endl;
      } else {
        for (uint32_t i = 0; i < statement.body.elements.size(); ++i) {
          std::cout << statement.body.elements[i].first << " ";
        }
        std::cout << std::endl;
      }

      if (statement.head.elements.size() == 0) {
        std::cout << "skip constraint" << std::endl;
      } else {
        for (uint32_t h = 0; h < statement.head.elements.size(); ++h) {
          std::vector<int> rule = {LPI_RULE_TYPE::AND_RULE};
          int atom = statement.head.elements[h];
          rule.emplace_back(atom);
          if (std::abs(atom) > atoms_count) {
            atoms_count = std::abs(atom);
          }
          for (uint32_t i = 0; i < statement.body.elements.size(); ++i) {
            int atom = statement.body.elements[i].first;
            rule.emplace_back(atom);
            if (std::abs(atom) > atoms_count) {
              atoms_count = std::abs(atom);
            }
            if (atom < 0) {
              atoms_neg.insert(-atom);
            }
          }
          std::cout << "Added rule: ";
          for (uint32_t i = 0; i < rule.size(); ++i) {
            std::cout << rule[i] << " ";
          }
          std::cout << std::endl;
          rules.emplace_back(rule);
        }
      }
    } else {
      std::cout << "Choice rule: " << statement << std::endl;
      std::cout << statement << std::endl;
    }
  }

  std::cout << "Max atoms count = " << atoms_count << std::endl;
  std::cout << "Input rules" << std::endl;
  printRules(rules);

  LPI_Format rules_standardized;
  standardizeRules(rules, atoms_count, rules_standardized);
  std::cout << "Max atoms count = " << atoms_count << std::endl;
  std::cout << "Standardized rules" << std::endl;
  printRules(rules_standardized);

  BitMatrix lp_matrix, lp_matrix_reduct;
  BitVector empty_bit_vector(2 * atoms_count);
  for (uint32_t i = 0; i < atoms_count; ++i) {
    lp_matrix.emplace_back(empty_bit_vector);
  }
  // printBitMatrix(lp_matrix);
  empty_bit_vector = BitVector(atoms_count);
  for (uint32_t i = 0; i < atoms_count; ++i) {
    lp_matrix_reduct.emplace_back(empty_bit_vector);
  }
  // printBitMatrix(lp_matrix_reduct);

  BitVector rule_type_vector(atoms_count);
  // std::cout << "Rule type vector" << std::endl;
  // std::cout << rule_type_vector << std::endl;

  for (auto rule : rules_standardized) {
    // printRule(rule); std::cout << std::endl;
    if (rule[LPI_IDX::RULE_TYPE]) {
      rule_type_vector[rule[LPI_IDX::RULE_HEAD - 1]] = 1;
    }
    for (uint32_t b = LPI_IDX::RULE_BODY; b < rule.size(); ++b) {
      uint32_t i = rule[LPI_IDX::RULE_HEAD] - 1, j;
      if (rule[b] > 0) {
        j = rule[b] - 1;
        lp_matrix_reduct[i][j] = 1;
      } else {
        j = atoms_count - rule[b] - 1;
      }
      // std::cout << "assign pos " << i << " " << j << " to 1" << std::endl;
      lp_matrix[i][j] = 1;
    }
    // std::cout << "bit " << lp_matrix[rule[LPI_IDX::RULE_HEAD] - 1] <<
    // std::endl;
  }

  // printBitMatrix(lp_matrix);
  // printBitMatrix(lp_matrix_reduct);
  // std::cout << "Rule type vector" << std::endl;
  // std::cout << rule_type_vector << std::endl;

  printLPBitMatrix(lp_matrix, rule_type_vector, atoms_count);
  printLPBitMatrix(lp_matrix_reduct, rule_type_vector, atoms_count);

  return 0;
}
