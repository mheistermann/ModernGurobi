/*
 * Copyright 2016 Martin Heistermann <code()mheistermann.de>
 *
 * This file is part of ModernGurobi.
 *
 * ModernGurobi s free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation, in version 3 of the License.
 *
 * ModernGurobi is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ModernGurobi.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MODERNGUROBI_HH
#define MODERNGUROBI_HH


#include <vector>
#include <string>
#include <stdexcept>
#include <limits>
#include <memory>
#include <algorithm>
#include <map>
#include <cassert>
#include <iostream>
#include <cmath>
#include <utility>

namespace ModernGurobi {
extern "C" {
#include <gurobi_c.h> // NOLINT(build/include_order)
}

class Model;
class AffineConstraint;

class GurobiException : public std::runtime_error
{
public:
  explicit GurobiException(const std::string &msg, int error = 0)
      : std::runtime_error(msg), error_(error) {}
  explicit GurobiException(const char *msg, int error = 0)
      : std::runtime_error(msg), error_(error) {}
  const std::string getMessage() const { return std::string(what()) + ", error = " + std::to_string(error_); }
  int getErrorCode() const {return error_; }
private:
  int error_;
};

class GurobiInfeasibleException : public GurobiException {
public:
    GurobiInfeasibleException() : GurobiException("Infeasible") {}
};

void throw_if_err(
        int error,
        std::string msg,
        std::string extra,
        std::string filename,
        std::string function,
        unsigned int lineno);

#define EXCEPTWRAP_EXTRA(CALL, EXTRA) do { \
        int error = (CALL); \
        if (error) { \
            throw_if_err(error, #CALL, EXTRA, __FILE__, __PRETTY_FUNCTION__, __LINE__); \
        } \
    } while(0)

#define EXCEPTWRAP(CALL) EXCEPTWRAP_EXTRA(CALL, "")

class Env
{
public:
    explicit Env(const std::string &logfilename = "") {
        EXCEPTWRAP(GRBloadenv(&env_, logfilename.c_str()));
        assert(env_ != nullptr);
    }
    ~Env() {
        GRBfreeenv(env_);
    }
    Env(const Env &other) = delete;
    Env& operator=(const Env &other) = delete;
private:
    friend class Model;
    GRBenv *env() {return env_;}
    GRBmodel *newModel() const {
        GRBmodel *model;
        EXCEPTWRAP(GRBnewmodel(env_, &model, "modelname_dummy", 0, nullptr, nullptr, nullptr, nullptr, nullptr));
        return model;
    }
    GRBenv *env_;
};

class Var
{
public:
    ~Var() = default;
    Var(const Var &other) = delete;
    Var &operator=(const Var &other) = delete;

    std::string to_string() const {
        return name_;
    }

    bool getBinary() const {
         // TODO: is there a better way to get binary vars out?
        return get() > 0.5;
    }

    double get() const {
        if (!computed) {
            throw GurobiException("GRBVar::get(): no value computed yet, maybe call GRBModel::optimize()?");
        }
        return solution_;
    }

protected:
    explicit Var(size_t idx, const std::string &name) : idx_(idx), name_(name) {}

private:
    friend class Model;
    friend class AffineConstraint;
    void set(double solution) {
        computed = true;
        solution_ = solution;
    }
    size_t idx() const {return idx_;}

    size_t idx_ = -1;
    std::string name_; // this is a bit wasteful, also stored in gurobi, but we need the model to access it
    bool computed = false;
    double solution_ = std::numeric_limits<double>::quiet_NaN();
};

/**
 * @brief The VarPtr class
 * Using std::shared_ptr<Var> directly as VarPtr and handing these out to clients
 * to use in expressions is problematic, as shared_ptr implements == and <, returning
 * bools as opposed to errors or the intended AffineConstraint.
 *
 * Wrap it to remove all non-needed access.
 */
class VarPtr
{
private:
    using T = Var;
    using Ptr = std::shared_ptr<T>;
    Ptr ptr_;
public:
    VarPtr() {}
    explicit VarPtr(Ptr &&ptr)
        : ptr_(std::move(ptr))
    {}

    inline T* operator->() const { return ptr_.operator->(); }

    /** Useful for std::map. */
    struct ptr_lessthan {
        bool operator()(const VarPtr &first, const VarPtr &second) {
            return first.ptr_ < second.ptr_;
        }
    };
};

class LinearExpr {
public:
    LinearExpr(VarPtr v) // NOLINT(runtime/explicit)
        : coeffmap_{{v, 1}}
    {}
    /**
     * @brief LinearExpr default constructor: value is 0
     */
    LinearExpr() {}
    ~LinearExpr() = default;
    LinearExpr(const LinearExpr &other) = default;
    LinearExpr &operator=(const LinearExpr &other) = default;

    static LinearExpr zero() { return LinearExpr(); }


    std::string to_string() const {
        if (coeffmap_.empty()) {
            return "0";
        }
        std::string ret = "";
        for(auto &entry: coeffmap_) {
            ret += " + ";
            if (entry.second != 1.0) {
                ret += std::to_string(entry.second) + " * ";
            }
            ret += entry.first->to_string();
        }
        return ret;
    }

    LinearExpr& operator*=(double d) {
        for(auto entry: coeffmap_) { entry.second *=d; }
        return *this;
    }

    LinearExpr operator-() const {
        LinearExpr ret;
        auto &cmap = ret.coeffmap();
        for(auto &entry: coeffmap_) { cmap[entry.first] = -entry.second; }
        return ret;
    }

    LinearExpr& operator+=(const ModernGurobi::LinearExpr &other) {
        for(auto &entry: other.coeffmap()) { coeffmap_[entry.first] += entry.second; }
        return *this;
    }

    LinearExpr& operator-=(const ModernGurobi::LinearExpr &other) {
        for(auto &entry: other.coeffmap()) { coeffmap_[entry.first] -= entry.second; }
        return *this;
    }

private:
    friend class AffineConstraint;
    LinearExpr(VarPtr v, double d): coeffmap_{{v,d}} {}

    friend LinearExpr operator+(LinearExpr x, const LinearExpr& y);
    friend LinearExpr operator-(LinearExpr x, const LinearExpr& y);
    friend LinearExpr operator*(const LinearExpr &x, double d);
    friend LinearExpr operator*(double d, const LinearExpr& x);

    using CoeffMap = std::map<VarPtr, double, VarPtr::ptr_lessthan>;
    const CoeffMap &coeffmap() const {return coeffmap_;}
    CoeffMap &coeffmap() {return coeffmap_;}

    CoeffMap coeffmap_;
};

LinearExpr operator+(LinearExpr x, const LinearExpr& y);
LinearExpr operator-(LinearExpr x, const LinearExpr& y);
LinearExpr operator*(const LinearExpr &x, double d);
LinearExpr operator*(double d, const LinearExpr& x);

class AffineExpr {
public:
    AffineExpr(double d)       // NOLINT(runtime/explicit)
        : linPart_(LinearExpr::zero()),
          constant_(d)
    {}
    AffineExpr(int i)          // NOLINT(runtime/explicit)
        : AffineExpr(static_cast<double>(i))
    {}
    AffineExpr(LinearExpr lin) // NOLINT(runtime/explicit)
        : linPart_(lin),
          constant_(0)
    {}
    ~AffineExpr() = default;
    AffineExpr(const AffineExpr &other) = default;
    AffineExpr &operator=(const AffineExpr &other) = default;
private:
    friend class AffineConstraint;
    friend AffineExpr operator+(AffineExpr x, const AffineExpr& y);
    friend AffineExpr operator-(AffineExpr x, const AffineExpr& y);
    friend AffineExpr operator*(AffineExpr expr, double d);
    friend AffineExpr operator*(double d, AffineExpr expr);

    friend AffineConstraint operator<=(const AffineExpr &x, const AffineExpr &y);
    friend AffineConstraint operator>=(const AffineExpr &x, const AffineExpr &y);
    friend AffineConstraint operator==(const AffineExpr &x, const AffineExpr &y);

    const double &getConstant() const {return constant_; }
    double &getConstant() {return constant_; }
    const LinearExpr &getLinearPart() const { return linPart_; }
    LinearExpr &getLinearPart() { return linPart_; }

    LinearExpr linPart_;
    double constant_;
};

AffineExpr operator+(AffineExpr x, const AffineExpr& y);
AffineExpr operator-(AffineExpr x, const AffineExpr& y);
AffineExpr operator*(AffineExpr expr, double d);
AffineExpr operator*(double d, AffineExpr expr);




class Constraint {
public:
    Constraint() {}
};

class AffineConstraint: public Constraint {
public:
    ~AffineConstraint() = default;
    AffineConstraint(const AffineConstraint &other) = default;
    AffineConstraint &operator=(const AffineConstraint &other) = default;

    std::string to_string() const {
        return name_ + ": " + expr_.to_string() + " " + sense_ + " " + std::to_string(rhs_);
    }
    const std::string &get_name() const {return name_;}
    void set_name(const std::string &name) {name_ = name;}

private:
    friend class Model;
    friend AffineConstraint operator<=(const AffineExpr &x, const AffineExpr &y);
    friend AffineConstraint operator>=(const AffineExpr &x, const AffineExpr &y);
    friend AffineConstraint operator==(const AffineExpr &x, const AffineExpr &y);

    AffineConstraint(const LinearExpr &expr, char sense, double rhs)
        : expr_(expr),
          sense_(sense),
          rhs_(rhs)
    {}

    void add_to_model(GRBmodel *model) const;

    LinearExpr expr_;
    char sense_;
    double rhs_;
    std::string name_;
};

AffineConstraint operator<=(const AffineExpr &x, const AffineExpr &y);
AffineConstraint operator>=(const AffineExpr &x, const AffineExpr &y);
AffineConstraint operator==(const AffineExpr &x, const AffineExpr &y);

inline AffineConstraint operator<=(const VarPtr &x, const VarPtr &y) {
    return LinearExpr(x) <= LinearExpr(y);
}

inline AffineConstraint operator==(const VarPtr &x, const VarPtr &y) {
    return LinearExpr(x) == LinearExpr(y);
}

inline AffineConstraint operator>=(const VarPtr &x, const VarPtr &y) {
    return LinearExpr(x) >= LinearExpr(y);
}

class Model
{
public:
    explicit Model(const Env &env)
        : model_(env.newModel())
    {
        assert(model_ != nullptr);
    }
    ~Model() {
        GRBfreemodel(model_);
    }
    Model(const Model &other) = delete;
    Model& operator=(const Model &other) = delete;

    /**
     * @brief addVar
     * @param lb        lower bound
     * @param ub        upper bound
     * @param obj       objective coefficient
     * @param vtype     variable type
     * @param vname     variable name (optional)
     * @return a VarPtr to the newly added var
     */
    VarPtr addVar(double lb, double ub, double obj, char vtype, const std::string &vname="")
    {
        struct make_shared_enabler : public Var { make_shared_enabler(size_t idx, const std::string &name) : Var(idx, name) {}};
        size_t idx = vars_.size();
        std::string name = vname;
        if (name.empty()) {
            name = std::string("v_") + std::to_string(idx);
        }
        vars_.emplace_back(std::make_shared<make_shared_enabler>(idx, name));
        assert(std::isfinite(obj));
        EXCEPTWRAP_EXTRA(GRBaddvar(model_,
                  0,
                  nullptr,
                  nullptr,
                  obj, lb, ub,
                  vtype,
                  name.c_str()),
                  (std::string("lb = ") + std::to_string(lb)
                   + ", ub = " + std::to_string(ub)
                   + ", obj = " + std::to_string(obj)
                   + ", name = " + vname
                   ));
        return vars_.back();
    }

    VarPtr addBinaryVar(double obj, std::string vname="") {
        return addVar(0, 1, obj, GRB_BINARY, vname);
    }

    void addConstraint(const AffineConstraint &constr, const std::string &name = "")
    {
        constr.add_to_model(model_);
        affine_constrs_.emplace_back(constr);
        affine_constrs_.back().set_name(name);
    }

    void update() {
         EXCEPTWRAP(GRBupdatemodel(model_));
    }

    void optimize() {
        EXCEPTWRAP(GRBoptimize(model_));
        int status = getOptimStatus();
        if (status == GRB_INFEASIBLE) {
            throw GurobiInfeasibleException();
        }
        const std::vector<double> solution = getSolution();
        for (size_t i=0; i < vars_.size(); i++) {
            vars_[i]->set(solution[i]);
        }
    }

    int getOptimStatus() const {
        int optimstatus;
        EXCEPTWRAP(GRBgetintattr(model_, GRB_INT_ATTR_STATUS, &optimstatus));
        // TODO: translate optimization status to proper enum?
        return optimstatus;
    }

    double getObjective() const {
        double objval;
        EXCEPTWRAP(GRBgetdblattr(model_, GRB_DBL_ATTR_OBJVAL, &objval));
        return objval;
    }

    void debugWhyInfeasible() const {
        EXCEPTWRAP(GRBcomputeIIS(model_));
        for(size_t idx = 0; idx < affine_constrs_.size(); ++idx) {
            int in_iis;
            EXCEPTWRAP(GRBgetintattrelement(model_, GRB_INT_ATTR_IIS_CONSTR, static_cast<int>(idx), &in_iis));
            if (in_iis) {
                std::cerr << "Constr in IIS: " << affine_constrs_[idx].to_string() << std::endl;
            }
        }
    }

    void setMIPGap(double gap) {
        GRBenv *env = GRBgetenv(model_);
        EXCEPTWRAP(GRBsetdblparam(env, GRB_DBL_PAR_MIPGAP, gap));
    }

    void setTimeLimitSeconds(const double limit) {
        GRBenv *env = GRBgetenv(model_);
        EXCEPTWRAP(GRBsetdblparam(env, GRB_DBL_PAR_TIMELIMIT, limit));
    }


private:
    const std::vector<double> &getSolution() {
        solution_.resize(vars_.size());
        EXCEPTWRAP(GRBgetdblattrarray(model_, GRB_DBL_ATTR_X,
                                       0, static_cast<int>(vars_.size()),
                                       &solution_[0]));
        return solution_;
    }
    mutable GRBmodel *model_ = 0;
    std::vector<VarPtr> vars_;
    std::vector<AffineConstraint> affine_constrs_; // For debug, e.g. infeasible constraints
    std::vector<double> solution_;
};


} // namespace ModernGurobi

#endif // MODERNGUROBI_HH
