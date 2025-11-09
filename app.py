import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import combinations
from matplotlib.patches import Polygon
from typing import List, Tuple, Optional
import re

st.set_page_config(
    page_title="Solveur Big M",
    page_icon="image.ico",
    layout="wide",
    initial_sidebar_state="expanded"
)

class BigM_Solver:
    def __init__(self, c: np.ndarray, A: np.ndarray, b: np.ndarray, 
                 constraint_types: List[str], M: float = 1e6, is_maximize: bool = False):
        self.c_original = np.array(c, dtype=float)
        if is_maximize:
            self.c_original = -self.c_original
        self.is_maximize = is_maximize
        
        self.A_original = np.array(A, dtype=float)
        self.b = np.array(b, dtype=float)
        self.constraint_types = constraint_types
        self.M = M
        
        self.n_vars = len(c)
        self.n_constraints = len(b)
        
        self.var_names = []
        self.base_vars = []
        self.tableau = None
        self.iteration = 0
        self.tableaux_history = []
        
    def _initialize_tableau(self):
        n_slack = sum(1 for ct in self.constraint_types if ct == '<=')
        n_surplus = sum(1 for ct in self.constraint_types if ct == '>=')
        n_artificial = sum(1 for ct in self.constraint_types if ct in ['>=', '='])
        
        total_vars = self.n_vars + n_slack + n_surplus + n_artificial
        self.tableau = np.zeros((self.n_constraints + 1, total_vars + 1))
        
        self.tableau[:-1, :self.n_vars] = self.A_original
        self.tableau[:-1, -1] = self.b
        
        self.var_names = [f'x{i+1}' for i in range(self.n_vars)]
        self.base_vars = []
        
        slack_idx = surplus_idx = artificial_idx = 0
        col_idx = self.n_vars
        
        for i, ct in enumerate(self.constraint_types):
            row_idx = i
            
            if ct == '<=':
                var_name = f's{slack_idx+1}'
                self.var_names.append(var_name)
                self.tableau[row_idx, col_idx] = 1
                self.base_vars.append(var_name)
                slack_idx += 1
                col_idx += 1
                
            elif ct == '>=':
                surplus_name = f'e{surplus_idx+1}'
                self.var_names.append(surplus_name)
                self.tableau[row_idx, col_idx] = -1
                surplus_idx += 1
                col_idx += 1
                
                artificial_name = f'A{artificial_idx+1}'
                self.var_names.append(artificial_name)
                self.tableau[row_idx, col_idx] = 1
                self.base_vars.append(artificial_name)
                artificial_idx += 1
                col_idx += 1
                
            elif ct == '=':
                artificial_name = f'A{artificial_idx+1}'
                self.var_names.append(artificial_name)
                self.tableau[row_idx, col_idx] = 1
                self.base_vars.append(artificial_name)
                artificial_idx += 1
                col_idx += 1
        
        z_row = -1
        self.tableau[z_row, :self.n_vars] = self.c_original
        
        for i, var_name in enumerate(self.var_names):
            if var_name.startswith('A'):
                self.tableau[z_row, i] = self.M
        
        for base_var in self.base_vars:
            if base_var.startswith('A'):
                base_idx = self.var_names.index(base_var)
                for row in range(self.n_constraints):
                    if self.tableau[row, base_idx] == 1:
                        self.tableau[z_row, :] -= self.M * self.tableau[row, :]
                        break
    
    def _get_tableau_df(self):
        col_names = self.var_names + ['RHS']
        row_names = self.base_vars + ['Z']
        df = pd.DataFrame(self.tableau, columns=col_names, index=row_names)
        return df
    
    def _find_entering_variable(self) -> Optional[int]:
        z_row = self.tableau[-1, :-1]
        min_idx = np.argmin(z_row)
        if z_row[min_idx] >= -1e-10:
            return None
        return min_idx
    
    def _find_leaving_variable(self, entering_idx: int) -> Optional[int]:
        pivot_column = self.tableau[:-1, entering_idx]
        rhs = self.tableau[:-1, -1]
        ratios = np.full(len(pivot_column), np.inf)
        for i, val in enumerate(pivot_column):
            if val > 1e-10:
                ratios[i] = rhs[i] / val
        if np.all(ratios == np.inf):
            return None
        return np.argmin(ratios)
    
    def _pivot(self, pivot_row: int, pivot_col: int):
        pivot_element = self.tableau[pivot_row, pivot_col]
        self.tableau[pivot_row, :] /= pivot_element
        for i in range(len(self.tableau)):
            if i != pivot_row:
                factor = self.tableau[i, pivot_col]
                self.tableau[i, :] -= factor * self.tableau[pivot_row, :]
        entering_var = self.var_names[pivot_col]
        self.base_vars[pivot_row] = entering_var
    
    def solve(self, store_history: bool = True) -> dict:
        self._initialize_tableau()
        self.iteration = 0
        if store_history:
            self.tableaux_history.append({
                'iteration': 0,
                'tableau': self._get_tableau_df().copy(),
                'base_vars': self.base_vars.copy()
            })
        max_iterations = 1000
        while self.iteration < max_iterations:
            self.iteration += 1
            entering_idx = self._find_entering_variable()
            if entering_idx is None:
                return self._extract_solution('optimal')
            leaving_row = self._find_leaving_variable(entering_idx)
            if leaving_row is None:
                return {'status': 'unbounded', 'message': 'Le problème est non borné'}
            self._pivot(leaving_row, entering_idx)
            if store_history:
                self.tableaux_history.append({
                    'iteration': self.iteration,
                    'tableau': self._get_tableau_df().copy(),
                    'base_vars': self.base_vars.copy()
                })
        return {'status': 'max_iterations', 'message': 'Nombre maximum d\'itérations atteint'}
    
    def _extract_solution(self, status: str) -> dict:
        for base_var in self.base_vars:
            if base_var.startswith('A'):
                row_idx = self.base_vars.index(base_var)
                if self.tableau[row_idx, -1] > 1e-6:
                    return {'status': 'infeasible', 'message': 'Le problème est infaisable'}
        x_solution = np.zeros(self.n_vars)
        for i, var_name in enumerate(self.var_names[:self.n_vars]):
            if var_name in self.base_vars:
                row_idx = self.base_vars.index(var_name)
                x_solution[i] = self.tableau[row_idx, -1]
        objective_value = -self.tableau[-1, -1]
        if self.is_maximize:
            objective_value = -objective_value
        return {
            'status': status,
            'x': x_solution,
            'objective_value': objective_value,
            'iterations': self.iteration,
            'base_vars': self.base_vars.copy(),
            'final_tableau': self._get_tableau_df()
        }

def parse_coefficients(text: str) -> List[float]:
    try:
        return [float(x.strip()) for x in text.split(',') if x.strip()]
    except ValueError:
        return None

def validate_inputs(obj_coeffs, constraints_data):
    if not obj_coeffs or len(obj_coeffs) == 0:
        return False, "Les coefficients de la fonction objectif sont invalides"
    if not constraints_data or len(constraints_data) == 0:
        return False, "Veuillez ajouter au moins une contrainte"
    n_vars = len(obj_coeffs)
    for i, c in enumerate(constraints_data):
        if len(c['coeffs']) != n_vars:
            return False, f"La contrainte {i+1} doit avoir {n_vars} coefficients"
        if c['rhs'] is None:
            return False, f"La contrainte {i+1} a un RHS invalide"
    return True, ""

def plot_graphical_solution(c, constraints_data):
    inequalities = []
    equality = None
    for cons in constraints_data:
        a, b = cons['coeffs']
        rhs = cons['rhs']
        t = cons['type']
        if t == '<=':
            inequalities.append((a, b, rhs))
        elif t == '>=':
            inequalities.append((-a, -b, -rhs))
        elif t == '=':
            equality = (a, b, rhs)

    def intersection(c1, c2):
        a1, b1, c1v = c1
        a2, b2, c2v = c2
        det = a1*b2 - a2*b1
        if abs(det) < 1e-10:
            return None
        x = (c1v*b2 - c2v*b1) / det
        y = (a1*c2v - a2*c1v) / det
        return (x, y)

    points = []
    if equality:
        for cineq in inequalities:
            p = intersection(cineq, equality)
            if p:
                points.append(p)

    feasible_points = []
    for (x, y) in points:
        if all(a*x + b*y <= c + 1e-6 for (a,b,c) in inequalities):
            feasible_points.append((x, y))
    if len(feasible_points) == 0:
        st.warning("Aucune zone réalisable trouvée pour la représentation graphique.")
        return
    feasible_points = np.array(feasible_points)
    Z_values = feasible_points @ np.array(c)
    opt_index = np.argmin(Z_values)
    x_opt, y_opt = feasible_points[opt_index]
    z_opt = Z_values[opt_index]
    st.info(f"Point optimal (graphique) : x1 = {x_opt:.2f}, x2 = {y_opt:.2f}, Z = {z_opt:.2f}")

    fig, ax = plt.subplots(figsize=(6, 6))
    x = np.linspace(0, max(feasible_points[:,0])*1.5, 400)
    for (a, b, c_val) in inequalities:
        if b != 0:
            ax.plot(x, (c_val - a*x)/b, 'k-', alpha=0.7)
        else:
            ax.axvline(c_val/a, color='k', alpha=0.7)
    if equality:
        ax.plot(x, (equality[2] - equality[0]*x)/equality[1], 'b-', linewidth=2)
    polygon = Polygon(feasible_points, color="#90EE90", alpha=0.5)
    ax.add_patch(polygon)
    ax.plot(feasible_points[:,0], feasible_points[:,1], 'ro')
    ax.plot(x_opt, y_opt, 'bo', markersize=10)
    ax.set_xlabel("x1")
    ax.set_ylabel("x2")
    ax.set_title("Méthode Graphique — Zone réalisable")
    ax.grid(True)
    st.pyplot(fig)

def main():
    st.title("Master Optimization Big M Way")
    st.markdown("---")
    with st.sidebar:
        st.header("Configuration")
        optimization_type = st.radio("Type d'optimisation", ["Minimiser", "Maximiser"], index=0)
        is_maximize = (optimization_type == "Maximiser")
        M_value = st.number_input("Valeur de M (Big M)", min_value=1000.0, max_value=1e10, value=1e6, format="%.0e")
        show_iterations = st.checkbox("Afficher les itérations", value=True)
        st.markdown("---")

    obj_input = st.text_input(f"{optimization_type} Z =", value="60, 80")
    obj_coeffs = parse_coefficients(obj_input)
    st.markdown("---")
    st.header("Contraintes")
    if 'constraints' not in st.session_state:
        st.session_state.constraints = [
            {'coeffs_text': '1, 1', 'type': '=', 'rhs': 8.0},
            {'coeffs_text': '1, 0', 'type': '>=', 'rhs': 2.0},
            {'coeffs_text': '0, 1', 'type': '>=', 'rhs': 3.0},
            {'coeffs_text': '1, 0', 'type': '<=', 'rhs': 5.0},
            {'coeffs_text': '0, 1', 'type': '<=', 'rhs': 6.0},
        ]
    constraints_to_remove = []
    for i, constraint in enumerate(st.session_state.constraints):
        col1, col2, col3, col4, col5 = st.columns([3, 1, 1, 1, 1])
        with col1:
            constraint['coeffs_text'] = st.text_input(f"Coefficients C{i+1}", value=constraint['coeffs_text'], key=f"coeffs_{i}")
        with col2:
            constraint['type'] = st.selectbox(f"Type C{i+1}", ['<=', '>=', '='], index=['<=', '>=', '='].index(constraint['type']), key=f"type_{i}")
        with col3:
            constraint['rhs'] = st.number_input(f"RHS C{i+1}", value=float(constraint['rhs']), key=f"rhs_{i}")
        with col5:
            if st.button("🗑️", key=f"remove_{i}"):
                constraints_to_remove.append(i)
    for i in reversed(constraints_to_remove):
        st.session_state.constraints.pop(i)
        st.rerun()
    if st.button("➕ Ajouter une contrainte"):
        if obj_coeffs:
            default_coeffs = ", ".join(['0'] * len(obj_coeffs))
            st.session_state.constraints.append({'coeffs_text': default_coeffs, 'type': '<=', 'rhs': 0.0})
            st.rerun()

    if st.button("Résoudre le problème", type="secondary", use_container_width=True):
        constraints_data = []
        for constraint in st.session_state.constraints:
            coeffs = parse_coefficients(constraint['coeffs_text'])
            if coeffs:
                constraints_data.append({'coeffs': coeffs, 'type': constraint['type'], 'rhs': constraint['rhs']})
        valid, error_msg = validate_inputs(obj_coeffs, constraints_data)
        if not valid:
            st.error(error_msg)
        else:
            c = np.array(obj_coeffs)
            A = np.array([c['coeffs'] for c in constraints_data])
            b = np.array([c['rhs'] for c in constraints_data])
            types = [c['type'] for c in constraints_data]
            solver = BigM_Solver(c, A, b, types, M=M_value, is_maximize=is_maximize)
            result = solver.solve(store_history=show_iterations)
            if result['status'] == 'optimal':
                st.subheader("Solution optimale")
                st.write(pd.DataFrame({'Variable': [f'x{i+1}' for i in range(len(result['x']))], 'Valeur': result['x']}))
                st.write("Valeur optimale :", result['objective_value'])
                if len(obj_coeffs) == 2:
                    st.markdown("---")
                    st.header("Représentation Graphique (2 variables)")
                    plot_graphical_solution(obj_coeffs, constraints_data)
            elif result['status'] == 'infeasible':
                st.error("Problème infaisable")
            elif result['status'] == 'unbounded':
                st.error("Problème non borné")
            else:
                st.error(result.get('message', 'Erreur inconnue'))

if __name__ == "__main__":
    main()
