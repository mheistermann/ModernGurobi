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

namespace ModernGurobi {

extern "C" {
#include "gurobi_c.h"
}


class Model;
class AffineConstraint;

class GurobiException : public std::runtime_error
{
public:
  GurobiException(const std::string &msg, int error = 0)
      : std::runtime_error(msg), error_(error) {}
  GurobiException(const char *msg, int error = 0)
      : std::runtime_error(msg), error_(error) {}
  const std::string getMessage() const { return std::string(what()) + ", error = " + std::to_string(error_); }
  int getErrorCode() const {return error_; }
private:
  int error_;
};

class Env
{
public:
    explicit Env(const std::string &logfilename = "gurobi.log") {
        int error = GRBloadenv(&env_, logfilename.c_str());
        if (error || env_ == nullptr) {
          throw GurobiException("GRBEnv::GRBEnv():  GRBloadenv failed", error);
        }
    }
    ~Env() {
        GRBfreeenv(env_);
    }
    Env(const Env &other) = delete;
    Env& operator=(const Env &other) = delete;

private:
    friend class Model;
    GRBenv *env_;
};

class Var
{
public:
    Var(const Var &other) = delete;
    Var &operator=(const Var &other) = delete;

    double get() {
        if (!computed) {
            throw GurobiException("GRBVar::get(): no value computed yet, maybe call GRBModel::optimize()?");
        }
        return solution_;
    }

protected:
    Var(size_t idx) : idx_(idx) {}

private:
    friend class Model;
    friend class AffineConstraint;
    void set(double solution) {
        computed = true;
        solution_ = solution;
    }
    size_t idx() const {return idx_;}

    size_t idx_ = -1;
    bool computed = false;
    double solution_ = std::numeric_limits<double>::quiet_NaN();
};

using VarPtr = std::shared_ptr<Var>;


class LinearExpr {
public:
    LinearExpr(VarPtr v): coeffmap_{{v, 1}} {}
    static LinearExpr zero() { return LinearExpr(); }
    LinearExpr& operator*=(double d) {
        for(auto entry: coeffmap_) { entry.second *=d; }
        return *this;
    }

    LinearExpr& operator-(const ModernGurobi::LinearExpr &y) {
        for(auto entry: coeffmap_) { entry.second = -entry.second; }
        return *this;
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
    LinearExpr() {}
    LinearExpr(VarPtr v, double d): coeffmap_{{v,d}} {}

    friend LinearExpr operator+(LinearExpr x, const LinearExpr& y);
    friend LinearExpr operator-(LinearExpr x, const LinearExpr& y);
    friend LinearExpr operator*(const LinearExpr &x, double d);
    friend LinearExpr operator*(double d, const LinearExpr& x);

    using CoeffMap = std::map<VarPtr, double>;
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
    AffineExpr(double d)
        : linPart_(LinearExpr::zero()),
          constant_(d)
    {}
    AffineExpr(int i) // just for convenient auto-conversion
        : AffineExpr(static_cast<double>(i))
    {}
    AffineExpr(LinearExpr lin)
        : linPart_(lin),
          constant_(0)
    {}
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
protected:
    friend class Model;
    virtual void add_to_model(GRBmodel *model) const = 0;
};

class AffineConstraint: public Constraint {
private:
    friend AffineConstraint operator<=(const AffineExpr &x, const AffineExpr &y);
    friend AffineConstraint operator>=(const AffineExpr &x, const AffineExpr &y);
    friend AffineConstraint operator==(const AffineExpr &x, const AffineExpr &y);

    AffineConstraint(const LinearExpr &expr, char sense, double rhs)
        : expr_(expr),
          sense_(sense),
          rhs_(rhs)
    {}

    virtual void add_to_model(GRBmodel *model) const;

    LinearExpr expr_;
    char sense_;
    double rhs_;
};

AffineConstraint operator<=(const AffineExpr &x, const AffineExpr &y);
AffineConstraint operator>=(const AffineExpr &x, const AffineExpr &y);
AffineConstraint operator==(const AffineExpr &x, const AffineExpr &y);

class Model
{
public:
    explicit Model(const Env &env) {
        int error = GRBnewmodel(env.env_, &model_, "modelname_dummy", 0, nullptr, nullptr, nullptr, nullptr, nullptr);
        if (error || model_ == nullptr) {
          throw GurobiException("GRBModel::GRBModel(): GRBnewmodel() failed", error);
        }
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
     * @return std::shared_ptr<GRBVar>
     */
    VarPtr addVar(double lb, double ub, double obj, char vtype, std::string vname="")
    {
        struct make_shared_enabler : public Var { make_shared_enabler(size_t idx) : Var(idx) {}};
        size_t idx = vars_.size();
        vars_.emplace_back(std::make_shared<make_shared_enabler>(idx));
        int  error = GRBaddvar(model_,
                  0,
                  nullptr,
                  nullptr,
                  obj, lb, ub,
                  vtype,
                  vname.c_str());
        if (error) {
            throw GurobiException("GRBModel::addVar(): GRBaddvar() failed", error);
        }
        return vars_.back();
    }

    VarPtr addBinaryVar(double obj, std::string vname="") {
        return addVar(0, 1, obj, GRB_BINARY, vname);
    }

    void addConstraint(const Constraint &constr)
    {
        constr.add_to_model(model_);
    }

    void update() {
         int error = GRBupdatemodel(model_);
         if (error) {
             throw GurobiException("GRBModel::update(): GRBupdatemodel() failed", error);
         }
    }

    void optimize() {
        int error = GRBoptimize(model_);
        if (error) {
            throw GurobiException("GRBModel::optimize(): GRBoptimize() failed", error);
        }
        const std::vector<double> solution = getSolution();
        for (size_t i=0; i < vars_.size(); i++) {
            vars_[i]->set(solution[i]);
        }
    }

    int getOptimStatus()  {
        int optimstatus;
        int error = GRBgetintattr(model_, GRB_INT_ATTR_STATUS, &optimstatus);
        if (error) {
            throw GurobiException("GRBModel::getOptimStatus: GRBgetintattr(model, GRB_INT_ATTR_STATUS, &...) failed", error);
        }
        // TODO: translate error status to proper enum?
        return optimstatus;
    }

    double getObjective() {
        double objval;
        int error = GRBgetdblattr(model_, GRB_DBL_ATTR_OBJVAL, &objval);
        if (error) {
            throw GurobiException("GRBModel::getObjective: GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &...) failed", error);
        }
        return objval;
    }

private:
    const std::vector<double> &getSolution() {
        solution_.resize(vars_.size());
        int error = GRBgetdblattrarray(model_, GRB_DBL_ATTR_X,
                                       0, static_cast<int>(vars_.size()),
                                       &solution_[0]);
        if (error) {
            throw GurobiException("GRBModel::getSolution: GRBgetdblattrarray(model, GRB_DBL_ATTR_X, &...) failed", error);
        }
        return solution_;
    }
    GRBmodel *model_ = 0;
    std::vector<std::shared_ptr<Var>> vars_;
    std::vector<Constraint> constrs_;
    std::vector<double> solution_;
};


} // namespace

#endif // MODERNGUROBI_HH
