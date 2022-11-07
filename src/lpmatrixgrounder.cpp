#include "lpmatrixgrounder.h"

#include <fstream>

std::ostream& operator<<(std::ostream& os, ASPIF_Statement<int>& data) {
    if (data.statement_type == ASPIF_Statement<int>::STATEMENT_TYPE_RULE) {
        ASPIF_RuleStatement<int>* statement = (ASPIF_RuleStatement<int>*)&data;
        os << statement->to_string();
    } else if (data.statement_type == ASPIF_Statement<int>::STATEMENT_TYPE_OUTPUT) {
        ASPIF_OutputStatement<int>* statement = (ASPIF_OutputStatement<int>*)&data;
        os << statement->to_string();
    }
    return os;
}

void ClingoObserver::rule(bool choice, Clingo::AtomSpan head, Clingo::LiteralSpan body) {
    // std::cout << "rule: " << "choice = " << choice << "\t head = " << head << "\t body = " << body << std::endl;

    ASPIF_RuleStatement<int> statement;
    statement.statement_type = ASPIF_Statement<int>::STATEMENT_TYPE_RULE;
    statement.head.rule_type = choice;
    for (auto i : head) {
        statement.head.elements.emplace_back(i);
    }
    statement.body.rulebody_type = ASPIF_RuleBody<int>::STATEMENT_RULE_BODY_TYPE_NORMAL;
    for (auto i : body) {
        statement.body.elements.emplace_back(std::make_pair<int, int>((int)i, 0));
    }
    rule_statements_.emplace_back(statement);
}

void ClingoObserver::weight_rule(bool choice, Clingo::AtomSpan head, Clingo::weight_t lower_bound, Clingo::WeightedLiteralSpan body) {
    // std::cout << "weight_rule: " << "choice = " << choice << "\t head = " << head << "\t weight = " << lower_bound << "\t body = ";
    // for (auto i : body) {
    //     std::cout << "<literal = " << i.literal()  << ", weight = " << i.weight() << ">; ";
    // }
    // std::cout << std::endl;

    ASPIF_RuleStatement<int> statement;
    statement.statement_type = ASPIF_Statement<int>::STATEMENT_TYPE_RULE;
    statement.head.rule_type = choice;
    for (auto i : head) {
        statement.head.elements.emplace_back(i);
    }
    statement.body.rulebody_type = ASPIF_RuleBody<int>::STATEMENT_RULE_BODY_TYPE_WEIGHT;
    statement.body.lbound = lower_bound;
    for (auto i : body) {
        statement.body.elements.emplace_back(std::make_pair<int, int>(i.literal(), i.weight()));
    }
    rule_statements_.emplace_back(statement);
}

void ClingoObserver::output_atom(Clingo::Symbol symbol, Clingo::atom_t atom) {
    // std::cout << "output_atom: " << "symbol = " << symbol << "\t atom = " << atom << std::endl;

    ASPIF_OutputStatement<int> statement;
    statement.statement_type = ASPIF_Statement<int>::STATEMENT_TYPE_OUTPUT;
    statement.str = symbol.to_string();
    if (atom > 0) {
        statement.literals.emplace_back(atom);
    }
    output_statements_.emplace_back(statement);
}

void ClingoObserver::output_term(Clingo::Symbol symbol, Clingo::LiteralSpan condition) {
    std::cout << "output_term: "
              << "symbol = " << symbol << "\t condition = " << condition << std::endl;
}

void ClingoObserver::output_csp(Clingo::Symbol symbol, int value, Clingo::LiteralSpan condition) {
    std::cout << "output_csp: "
              << "symbol = " << symbol << "\t value = " << value << "\t condition = " << condition << std::endl;
}

ClingoGrounder::ClingoGrounder() {
    this->logger_ = [](Clingo::WarningCode, char const* message) {
        std::cerr << message << std::endl;
    };
    // this->ctl_ = Clingo::Control({argv + 1, size_t(argc - 1)}, this->logger_, 20);
    this->ctl_ = Clingo::Control({nullptr, size_t(nullptr)}, this->logger_, 20);
    this->ctl_.register_observer(this->obs_);
}

void ClingoGrounder::input_files(std::initializer_list<std::string> listof_files) {
    std::string line;
    for (auto filename : listof_files) {
        std::cout << "Reading file " << filename << std::endl;
        std::ifstream ifile(filename);
        while (std::getline(ifile, line)) {
            this->ctl_.add("base", {}, line.c_str());
        }
        ifile.close();
    }
}

void ClingoGrounder::ground() {
    this->ctl_.ground({{"base", {}}});
}