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

#include "ModernGurobi.hh"

namespace ModernGurobi {

LinearExpr operator+(LinearExpr x, const LinearExpr &y)
{ return x += y; }

LinearExpr operator-(LinearExpr x, const LinearExpr &y)
{ return x -= y; }

LinearExpr operator*(const LinearExpr &x, double d)
{
    LinearExpr ret;
    LinearExpr::CoeffMap &cmap = ret.coeffmap();
    for(auto& entry: x.coeffmap()) {
        cmap[entry.first] = d * entry.second;
    }
    return ret;
}

LinearExpr operator*(double d, const LinearExpr &x)
{ return x*d; }

AffineExpr operator+(AffineExpr x, const AffineExpr& y) {
    x.getLinearPart() += y.getLinearPart();
    x.getConstant() += y.getConstant();
    return x;
}

AffineExpr operator-(AffineExpr x, const AffineExpr& y) {
    x.getLinearPart() -= y.getLinearPart();
    x.getConstant() -= y.getConstant();
    return x;
}

AffineExpr operator*(AffineExpr expr, double d) {
    expr.getLinearPart() *= d;
    expr.getConstant() *= d;
    return expr;
}

AffineExpr operator*(double d, AffineExpr expr) {
    return expr * d;
}

AffineConstraint operator<=(const AffineExpr &x, const AffineExpr &y) {
    //     xl + xc <= yl + yc
    // <=> xl - yl <= yc - xc
    return AffineConstraint(x.getLinearPart()-y.getLinearPart(),
                            GRB_LESS_EQUAL,
                            y.getConstant() - x.getConstant());
}

AffineConstraint operator>=(const AffineExpr &x, const AffineExpr &y) {
    return y <= x;
}

AffineConstraint operator==(const AffineExpr &x, const AffineExpr &y) {
    return AffineConstraint(x.getLinearPart()-y.getLinearPart(),
                            GRB_EQUAL,
                            y.getConstant() - x.getConstant());
}

void AffineConstraint::add_to_model(GRBmodel *model) const {
    std::vector<int> indices;
    std::vector<double> coeffs;
    const LinearExpr::CoeffMap &parts = expr_.coeffmap();
    if (parts.empty()) {
        // Avoid an "attempt to access an element in an empty container." if STL bounds checking is enabled.
        // As we pass a size of 0, these values will never be read!
        // however we can't just skip empty constraints, as they might be a contradiction,
        // or even if not, it would mess up the equivalence of our and gurobi's constraint indices
        indices.push_back(0);
        coeffs.push_back(0);
    }
    for(auto &entry: parts) {
        indices.push_back(static_cast<int>(entry.first->idx()));
        coeffs.push_back(entry.second);
    }
    int error = GRBaddconstr(model,
                             static_cast<int>(parts.size()),
                             &indices.front(),
                             &coeffs.front(),
                             sense_, rhs_, nullptr);
    if (error) {
        throw GurobiException("LinearConstraint::add_to_model(): GRBaddconstr() failed", error);
    }
}

void throw_if_err(
        int error,
        std::string msg,
        std::string extra,
        std::string filename,
        std::string function,
        unsigned int lineno) {
    switch(error) {
    case 0:
        return; // no error, no exception!
        // TODO: maybe return/raise specialized exceptions?
#define _GRB_ERROR_CASE(ERR) case ERR: msg = msg + " (" + # ERR + ")"; break;
        _GRB_ERROR_CASE(GRB_ERROR_OUT_OF_MEMORY)
                _GRB_ERROR_CASE(GRB_ERROR_NULL_ARGUMENT)
                _GRB_ERROR_CASE(GRB_ERROR_INVALID_ARGUMENT)
                _GRB_ERROR_CASE(GRB_ERROR_UNKNOWN_ATTRIBUTE)
                _GRB_ERROR_CASE(GRB_ERROR_DATA_NOT_AVAILABLE)
                _GRB_ERROR_CASE(GRB_ERROR_INDEX_OUT_OF_RANGE)
                _GRB_ERROR_CASE(GRB_ERROR_UNKNOWN_PARAMETER)
                _GRB_ERROR_CASE(GRB_ERROR_VALUE_OUT_OF_RANGE)
                _GRB_ERROR_CASE(GRB_ERROR_NO_LICENSE)
                _GRB_ERROR_CASE(GRB_ERROR_SIZE_LIMIT_EXCEEDED)
                _GRB_ERROR_CASE(GRB_ERROR_CALLBACK)
                _GRB_ERROR_CASE(GRB_ERROR_FILE_READ)
                _GRB_ERROR_CASE(GRB_ERROR_FILE_WRITE)
                _GRB_ERROR_CASE(GRB_ERROR_NUMERIC)
                _GRB_ERROR_CASE(GRB_ERROR_IIS_NOT_INFEASIBLE)
                _GRB_ERROR_CASE(GRB_ERROR_NOT_FOR_MIP)
                _GRB_ERROR_CASE(GRB_ERROR_OPTIMIZATION_IN_PROGRESS)
                _GRB_ERROR_CASE(GRB_ERROR_DUPLICATES)
                _GRB_ERROR_CASE(GRB_ERROR_NODEFILE)
                _GRB_ERROR_CASE(GRB_ERROR_Q_NOT_PSD)
                _GRB_ERROR_CASE(GRB_ERROR_QCP_EQUALITY_CONSTRAINT)
                _GRB_ERROR_CASE(GRB_ERROR_NETWORK)
                _GRB_ERROR_CASE(GRB_ERROR_JOB_REJECTED)
                _GRB_ERROR_CASE(GRB_ERROR_NOT_SUPPORTED)
                _GRB_ERROR_CASE(GRB_ERROR_EXCEED_2B_NONZEROS)
                _GRB_ERROR_CASE(GRB_ERROR_INVALID_PIECEWISE_OBJ)
                _GRB_ERROR_CASE(GRB_ERROR_UPDATEMODE_CHANGE)
        #undef _GRB_ERROR_CASE
            default: msg += "(unknown!)";
    }

    msg += " in function " + function + " in file " + filename + ":" + std::to_string(lineno);
    if (!extra.empty()) {
        msg += std::string(" (") +  extra + ")";
    }
    throw GurobiException(msg, error);
}

} // namespace ModernGurobi
