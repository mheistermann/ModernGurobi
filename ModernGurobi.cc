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
    for(auto entry: x.coeffmap()) {
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
    for(auto &entry: parts) {
        indices.push_back(static_cast<int>(entry.first->idx()));
        coeffs.push_back(entry.second);
    }
    int error = GRBaddconstr(model,
                             static_cast<int>(parts.size()),
                             &indices[0],
                             &coeffs[0],
                             sense_, rhs_, nullptr);
    if (error) {
        throw GurobiException("LinearConstraint::add_to_model(): GRBaddconstr() failed", error);
    }
}

} // namespace
