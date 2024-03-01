import sympy
from sympy import symbols, factor, expand
import gurobipy as gp
from gurobipy import GRB

try:
    # Create a new model
    m = gp.Model("4NodeGraph")

    # Create variables
    x_1 = m.addVar(vtype=GRB.BINARY, name="x_1")
    x_2 = m.addVar(vtype=GRB.BINARY, name="x_2")
    x_3 = m.addVar(vtype=GRB.BINARY, name="x_3")
    x_4 = m.addVar(vtype=GRB.BINARY, name="x_4")
    x_5 = m.addVar(vtype=GRB.BINARY, name="x_5")
    x_6 = m.addVar(vtype=GRB.BINARY, name="x_6")
    x_7 = m.addVar(vtype=GRB.BINARY, name="x_7")
    x_8 = m.addVar(vtype=GRB.BINARY, name="x_8")

    x_9 = m.addVar(vtype=GRB.BINARY, name="x_9")
    x_10 = m.addVar(vtype=GRB.BINARY, name="x_10")
    x_11 = m.addVar(vtype=GRB.BINARY, name="x_11")
    x_12 = m.addVar(vtype=GRB.BINARY, name="x_12")
    x_13 = m.addVar(vtype=GRB.BINARY, name="x_13")

    # x_14 = m.addVar(vtype=GRB.BINARY, name="x_14")
    # x_15 = m.addVar(vtype=GRB.BINARY, name="x_15")
    y_r1_1 = m.addVar(vtype=GRB.BINARY, name="y_r1_1")
    y_r1 = m.addVar(vtype=GRB.BINARY, name="y_r1")
    # y_r2 = m.addVar(vtype=GRB.BINARY, name="y_r2")

    y_r1_2 = m.addVar(vtype=GRB.BINARY, name="y_r1_2")
    y_r1_3 = m.addVar(vtype=GRB.BINARY, name="y_r1_3")
    y_r1_4 = m.addVar(vtype=GRB.BINARY, name="y_r1_4")
    y_r1_5 = m.addVar(vtype=GRB.BINARY, name="y_r1_5")
    y_r1_6 = m.addVar(vtype=GRB.BINARY, name="y_r1_6")
    y_r1_7 = m.addVar(vtype=GRB.BINARY, name="y_r1_7")
    y_r1_8 = m.addVar(vtype=GRB.BINARY, name="y_r1_8")
    y_r1_9 = m.addVar(vtype=GRB.BINARY, name="y_r1_9")

    # y_c1 = m.addVar(vtype=GRB.BINARY, name="y_c1")
    objective = x_1 * x_2 + 2 * x_1 * y_r1_1 - x_10 * y_r1_8 + 2 * x_10 * y_r1_9 + 2 * x_11 * y_r1 - x_11 * y_r1_9 - 20 * x_12 * x_13 + 20 * x_13 * y_r1 + 2 * x_2 * y_r1_1 - x_3 * y_r1_1 + 2 * x_3 * y_r1_2 + x_3 - x_4 * y_r1_2 + 2 * x_4 * y_r1_3 + x_4 - x_5 * y_r1_3 + 2 * x_5 * y_r1_4 + x_5 - x_6 * y_r1_4 + 2 * x_6 * y_r1_5 + x_6 - x_7 * y_r1_5 + 2 * x_7 * y_r1_6 + x_7 - x_8 * y_r1_6 + 2 * x_8 * y_r1_7 - x_9 * y_r1_7 + 2 * x_9 * y_r1_8 - 2 * y_r1 * y_r1_9 - 9 * y_r1 - 2 * y_r1_1 * y_r1_2 - 2 * y_r1_2 * y_r1_3 + 2 * y_r1_2 - 2 * y_r1_3 * y_r1_4 + 2 * y_r1_3 - 2 * y_r1_4 * y_r1_5 + 2 * y_r1_4 - 2 * y_r1_5 * y_r1_6 + 2 * y_r1_5 - 2 * y_r1_6 * y_r1_7 + 2 * y_r1_6 - 2 * y_r1_7 * y_r1_8 + 2 * y_r1_7 - 2 * y_r1_8 * y_r1_9 + 2 * y_r1_8 + 2 * y_r1_9 + 21

    # Set objective
    m.setObjective(objective, GRB.MINIMIZE)
    # k1*(x_2 + x_3 + x_5) + k2*((1 - x_8)) + k3*((1 - x_6 - y_r1 + 2 * x_6 * y_r1) + (1 - x_7 - y_r2 + 2 * x_7 * y_r2))
    # + k4*((x_8 + y_c1 - 2 * x_8 * y_c1)) + k5*(((1 - x_1) * (1 - x_2) - 2 * (1 - x_1) * y_r1_1 - 2 * (1 - x_2) * y_r1_1 + 3 *
    # y_r1_1) + ((1 - y_r1_1) * (1 - x_3) - 2 * (1 - y_r1_1) * y_r1 - 2 * (1 - x_3) * y_r1 + 3 * y_r1) + ((1 - x_4) * (1 - x_5)
    # - 2 * (1 - x_4) * y_r2 - 2 * (1 - x_5) * y_r2 + 3 * y_r2)) + k6*((x_6 * x_7 - 2 * x_6 * y_c1 - 2 * x_7 * y_c1 + 3 * y_c1))

    # is the same as
    # obj = gp.LinExpr()
    # obj += x
    # obj += y
    # obj += 2 * z
    # m.setObjective(obj, GRB.MAXIMIZE)

    # Add constraint: x + 2 y + 3 z <= 4
    ##m.addConstr(x + 2 * y + 3 * z <= 4, "c0")

    # Add constraint: x + y >= 1
    ## m.addConstr(x + y >= 1, "c1")

    # Optimize model
    m.optimize()

    for v in m.getVars():
        print('%s %g' % (v.VarName, v.X))

    print('Obj: %g' % m.ObjVal)

except gp.GurobiError as e:
    print('Error code ' + str(e.errno) + ': ' + str(e))

except AttributeError:
    print('Encountered an attribute error')