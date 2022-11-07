#ifndef LPMATRIXGROUNDER_LIBRARY_H
#define LPMATRIXGROUNDER_LIBRARY_H

#include <clingo.hh>
#include <sstream>
#include <vector>

template <typename T> struct ASPIF_Statement {
    char statement_type;

    enum STATEMENT_TYPE {
        STATEMENT_TYPE_RULE = 1,
        STATEMENT_TYPE_OUTPUT = 4
    };

    virtual std::string to_string() = 0;
};

std::ostream& operator<<(std::ostream& os, ASPIF_Statement<int>& data);

template <typename T> struct ASPIF_RuleHead {
    bool rule_type;
    std::vector<T> elements;
};

template <typename T> struct ASPIF_RuleBody {
    char rulebody_type;
    T lbound;
    std::vector<std::pair<T, T>> elements;

    enum STATEMENT_RULE_BODY_TYPE {
        STATEMENT_RULE_BODY_TYPE_NORMAL = 0,
        STATEMENT_RULE_BODY_TYPE_WEIGHT = 1
    };
};

template <typename T> struct ASPIF_RuleStatement : ASPIF_Statement<T> {
    ASPIF_RuleHead<T> head;
    ASPIF_RuleBody<T> body;

    std::string to_string() {
        std::stringstream ss;
        ss << (int)this->statement_type;
        ss << " " << (this->head.rule_type ? "1" : "0") << " " << this->head.elements.size();
        for (auto i : this->head.elements) {
            ss << " " << i;
        }
        ss << " " << (int)this->body.rulebody_type;
        if (this->body.rulebody_type == ASPIF_RuleBody<int>::STATEMENT_RULE_BODY_TYPE_NORMAL) {
            ss << " " << this->body.elements.size();
            for (auto i : this->body.elements) {
                ss << " " << i.first;
            }
        } else if (this->body.rulebody_type == ASPIF_RuleBody<int>::STATEMENT_RULE_BODY_TYPE_WEIGHT) {
            ss << " " << this->body.lbound << " " << this->body.elements.size();
            for (auto i : this->body.elements) {
                ss << " " << i.first << " " << i.second;
            }
        }
        return ss.str();
    }
};

template <typename T> struct ASPIF_OutputStatement : ASPIF_Statement<T> {
    std::string str;
    std::vector<T> literals;

    std::string to_string() {
        std::stringstream ss;
        ss << (int)this->statement_type;
        ss << " " << this->str.length() << " " << this->str << " " << this->literals.size();
        for (auto i : this->literals) {
            ss << " " << i;
        }
        return ss.str();
    }
};

class ClingoObserver : public Clingo::GroundProgramObserver {
   public:
    // rule statement
    void rule(bool choice, Clingo::AtomSpan head, Clingo::LiteralSpan body);
    void weight_rule(bool choice, Clingo::AtomSpan head, Clingo::weight_t lower_bound, Clingo::WeightedLiteralSpan body);

    // output statement
    void output_atom(Clingo::Symbol symbol, Clingo::atom_t atom);
    void output_term(Clingo::Symbol symbol, Clingo::LiteralSpan condition);
    void output_csp(Clingo::Symbol symbol, int value, Clingo::LiteralSpan condition);

    // private:
    std::vector<ASPIF_RuleStatement<int>> rule_statements_;
    std::vector<ASPIF_OutputStatement<int>> output_statements_;
};

class ClingoGrounder {
   public:
    ClingoGrounder();
    void input_files(std::initializer_list<std::string> listof_files);
    void ground();
    // private:
    Clingo::Logger logger_;
    Clingo::Control ctl_;
    ClingoObserver obs_;
};

#endif  //LPMATRIXGROUNDER_LIBRARY_H
