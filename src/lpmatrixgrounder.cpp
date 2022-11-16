#include "lpmatrixgrounder.h"

#include <chrono>
#include <fstream>
#include <unordered_set>

std::ostream& operator<<(std::ostream& os, ASPIF_Statement<int>& data) {
  if (data.statement_type == ASPIF_Statement<int>::STATEMENT_TYPE_RULE) {
    ASPIF_RuleStatement<int>* statement = (ASPIF_RuleStatement<int>*)&data;
    os << statement->to_string();
  } else if (data.statement_type ==
             ASPIF_Statement<int>::STATEMENT_TYPE_OUTPUT) {
    ASPIF_OutputStatement<int>* statement = (ASPIF_OutputStatement<int>*)&data;
    os << statement->to_string();
  }
  return os;
}

void ClingoObserver::rule(bool choice, Clingo::AtomSpan head,
                          Clingo::LiteralSpan body) {
  ASPIF_RuleStatement<int> statement;
  statement.statement_type = ASPIF_Statement<int>::STATEMENT_TYPE_RULE;
  statement.head.rule_type = choice;
  for (auto i : head) {
    statement.head.elements.emplace_back(i);
  }
  statement.body.rulebody_type =
      ASPIF_RuleBody<int>::STATEMENT_RULE_BODY_TYPE_NORMAL;
  for (auto i : body) {
    statement.body.elements.emplace_back(std::make_pair<int, int>((int)i, 0));
  }
  rule_statements_.emplace_back(statement);
}

void ClingoObserver::weight_rule(bool choice, Clingo::AtomSpan head,
                                 Clingo::weight_t lower_bound,
                                 Clingo::WeightedLiteralSpan body) {
  ASPIF_RuleStatement<int> statement;
  statement.statement_type = ASPIF_Statement<int>::STATEMENT_TYPE_RULE;
  statement.head.rule_type = choice;
  for (auto i : head) {
    statement.head.elements.emplace_back(i);
  }
  statement.body.rulebody_type =
      ASPIF_RuleBody<int>::STATEMENT_RULE_BODY_TYPE_WEIGHT;
  statement.body.lbound = lower_bound;
  for (auto i : body) {
    statement.body.elements.emplace_back(
        std::make_pair<int, int>(i.literal(), i.weight()));
  }
  rule_statements_.emplace_back(statement);
}

void ClingoObserver::output_atom(Clingo::Symbol symbol, Clingo::atom_t atom) {
  ASPIF_OutputStatement<int> statement;
  statement.statement_type = ASPIF_Statement<int>::STATEMENT_TYPE_OUTPUT;
  statement.str = symbol.to_string();
  if (atom > 0) {
    statement.literals.emplace_back(atom);
  }
  output_statements_.emplace_back(statement);
}

void ClingoObserver::output_term(Clingo::Symbol symbol,
                                 Clingo::LiteralSpan condition) {
  std::cout << "output_term: "
            << "symbol = " << symbol << "\t condition = " << condition
            << std::endl;
}

void ClingoObserver::output_csp(Clingo::Symbol symbol, int value,
                                Clingo::LiteralSpan condition) {
  std::cout << "output_csp: "
            << "symbol = " << symbol << "\t value = " << value
            << "\t condition = " << condition << std::endl;
}

ClingoGrounder::ClingoGrounder() {
  this->logger_ = [](Clingo::WarningCode, char const* message) {
    std::cerr << message << std::endl;
  };
  this->ctl_ = Clingo::Control({nullptr, size_t(nullptr)}, this->logger_, 20);
  this->ctl_.register_observer(this->obs_);
}

void ClingoGrounder::input_files(
    std::initializer_list<std::string> listof_files) {
  std::string line;
  for (auto filename : listof_files) {
    std::ifstream ifile(filename);
    while (std::getline(ifile, line)) {
      this->ctl_.add("base", {}, line.c_str());
    }
    ifile.close();
  }
}

void ClingoGrounder::input_files(const char** files, size_t no_files) {
  std::string line;
  for (size_t i = 0; i < no_files; ++i) {
    std::ifstream ifile(files[i]);
    while (std::getline(ifile, line)) {
      this->ctl_.add("base", {}, line.c_str());
    }
    ifile.close();
  }
}

void ClingoGrounder::ground() { this->ctl_.ground({{"base", {}}}); }

std::ostream& operator<<(std::ostream& os, std::vector<int>& rule) {
  os << rule[LPI_IDX::RULE_HEAD];
  if (rule.size() > LPI_IDX::RULE_BODY) {
    os << " <= ";
    for (uint32_t i = LPI_IDX::RULE_BODY; i < rule.size() - 1; ++i) {
      os << rule[i] << ", ";
    }
    os << rule[rule.size() - 1];
  }
  os << ".";
  return os;
}

std::ostream& operator<<(std::ostream& os,
                         std::vector<std::vector<int>>& rules) {
  for (std::vector<int> r : rules) {
    os << r << std::endl;
  }
  return os;
}

void printRule(std::vector<int> rule) {
  std::cout << (rule[LPI_IDX::RULE_TYPE] ? "OR-rule:  " : "AND-rule: ");
  std::cout << rule;
}

void printRules(std::vector<std::vector<int>> rules) { std::cout << rules; }

// preallocating the buffer to be used in grounding
char* cstr = new char[100'000'000];
struct_ground_ret ground_return(ClingoGrounder& cgrounder) {
  auto start = std::chrono::high_resolution_clock::now();
  cgrounder.ground();
  auto duration_clingo = std::chrono::high_resolution_clock::now() - start;

  start = std::chrono::high_resolution_clock::now();
  LPI_Format rules;
  uint32_t atoms_count = 0;
  std::unordered_set<int> atoms_neg;
  for (ASPIF_RuleStatement<int> statement : cgrounder.obs_.rule_statements_) {
    if (!statement.head.rule_type) {
      if (statement.head.elements.size() == 0) {
        std::vector<int> rule = {LPI_RULE_TYPE::AND_RULE};
        rule.emplace_back(0);
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
        rules.emplace_back(rule);
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
          rules.emplace_back(rule);
        }
      }
    } else {
      for (uint32_t h = 0; h < statement.head.elements.size(); ++h) {
        std::vector<int> rule = {LPI_RULE_TYPE::AND_RULE};
        int atom = statement.head.elements[h];
        rule.emplace_back(atom);
        for (uint32_t hh = 0; hh < statement.head.elements.size(); ++hh) {
          if (h != hh) {
            atom = statement.head.elements[hh];
            rule.emplace_back(-atom);
          }
        }
        rules.emplace_back(rule);
      }
    }
  }

  std::stringstream ss;
  ss << rules;
  const std::string my_str = ss.str();
  auto len = my_str.length();
  strcpy(cstr, my_str.c_str());
  cstr[len] = 0;
  auto duration = std::chrono::high_resolution_clock::now() - start;

  struct_ground_ret ret;
  ret.rawdata = cstr;
  ret.duration_clingo = std::chrono::duration<float>(duration_clingo).count();
  ret.duration_internal = std::chrono::duration<float>(duration).count();

  return ret;
}

struct_ground_ret ground_single(const char* filename) {
  ClingoGrounder cgrounder = ClingoGrounder();
  cgrounder.input_files({filename});
  return ground_return(cgrounder);
}

struct_ground_ret ground_list(const char** files, size_t no_files) {
  ClingoGrounder cgrounder = ClingoGrounder();
  cgrounder.input_files(files, no_files);
  return ground_return(cgrounder);
}