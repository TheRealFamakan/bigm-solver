import streamlit as st
import numpy as np
import pandas as pd
from typing import List, Tuple, Optional
import re

# Configuration de la page
st.set_page_config(
    page_title="Solveur Big M",
    page_icon="image.ico",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ============================================================================
# PARTIE 3 : SOLVEUR BIG M (Backend)
# ============================================================================

class BigM_Solver:
    """
    Implémentation robuste de la méthode Big M pour la programmation linéaire.
    """
    
    def __init__(self, c: np.ndarray, A: np.ndarray, b: np.ndarray, 
                 constraint_types: List[str], M: float = 1e6, is_maximize: bool = False):
        """
        Initialise le solveur Big M.
        
        Parameters:
        -----------
        c : np.ndarray
            Coefficients de la fonction objectif
        A : np.ndarray
            Matrice des coefficients des contraintes
        b : np.ndarray
            Vecteur du côté droit des contraintes
        constraint_types : List[str]
            Types de contraintes: '<=', '>=' ou '='
        M : float
            Valeur de pénalité Big M
        is_maximize : bool
            True pour maximisation, False pour minimisation
        """
        self.c_original = np.array(c, dtype=float)
        if is_maximize:
            self.c_original = -self.c_original  # Convertir max en min
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
        """Initialise le tableau simplex."""
        n_slack = sum(1 for ct in self.constraint_types if ct == '<=')
        n_surplus = sum(1 for ct in self.constraint_types if ct == '>=')
        n_artificial = sum(1 for ct in self.constraint_types if ct in ['>=', '='])
        
        total_vars = self.n_vars + n_slack + n_surplus + n_artificial
        self.tableau = np.zeros((self.n_constraints + 1, total_vars + 1))
        
        self.tableau[:-1, :self.n_vars] = self.A_original
        self.tableau[:-1, -1] = self.b
        
        self.var_names = [f'x{i+1}' for i in range(self.n_vars)]
        self.base_vars = []
        
        slack_idx = 0
        surplus_idx = 0
        artificial_idx = 0
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
        """Retourne le tableau sous forme de DataFrame."""
        col_names = self.var_names + ['RHS']
        row_names = self.base_vars + ['Z']
        
        df = pd.DataFrame(
            self.tableau,
            columns=col_names,
            index=row_names
        )
        return df
    
    def _find_entering_variable(self) -> Optional[int]:
        """Trouve la variable entrante."""
        z_row = self.tableau[-1, :-1]
        min_idx = np.argmin(z_row)
        
        if z_row[min_idx] >= -1e-10:
            return None
        
        return min_idx
    
    def _find_leaving_variable(self, entering_idx: int) -> Optional[int]:
        """Trouve la variable sortante."""
        pivot_column = self.tableau[:-1, entering_idx]
        rhs = self.tableau[:-1, -1]
        
        ratios = np.full(len(pivot_column), np.inf)
        for i, val in enumerate(pivot_column):
            if val > 1e-10:
                ratios[i] = rhs[i] / val
        
        if np.all(ratios == np.inf):
            return None
        
        leaving_row = np.argmin(ratios)
        return leaving_row
    
    def _pivot(self, pivot_row: int, pivot_col: int):
        """Effectue l'opération de pivot."""
        pivot_element = self.tableau[pivot_row, pivot_col]
        self.tableau[pivot_row, :] /= pivot_element
        
        for i in range(len(self.tableau)):
            if i != pivot_row:
                factor = self.tableau[i, pivot_col]
                self.tableau[i, :] -= factor * self.tableau[pivot_row, :]
        
        leaving_var = self.base_vars[pivot_row]
        entering_var = self.var_names[pivot_col]
        self.base_vars[pivot_row] = entering_var
    
    def solve(self, store_history: bool = True) -> dict:
        """Résout le problème."""
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
                    'base_vars': self.base_vars.copy(),
                    'entering': self.var_names[entering_idx],
                    'leaving': self.tableaux_history[-1]['base_vars'][leaving_row]
                })
        
        return {'status': 'max_iterations', 'message': 'Nombre maximum d\'itérations atteint'}
    
    def _extract_solution(self, status: str) -> dict:
        """Extrait la solution finale."""
        for base_var in self.base_vars:
            if base_var.startswith('A'):
                base_idx = self.var_names.index(base_var)
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

# ============================================================================
# PARTIE 4 : INTÉGRATION FRONTEND-BACKEND
# ============================================================================

def parse_coefficients(text: str) -> List[float]:
    """Parse une chaîne de coefficients séparés par des virgules."""
    try:
        coeffs = [float(x.strip()) for x in text.split(',') if x.strip()]
        return coeffs
    except ValueError:
        return None

def validate_inputs(obj_coeffs, constraints_data):
    """Valide les entrées utilisateur."""
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

# ============================================================================
# APPLICATION STREAMLIT PRINCIPALE
# ============================================================================

def main():
    st.title("Master Optimization Big M Way")
    st.markdown("---")
    
    # Sidebar pour les paramètres
    with st.sidebar:
        st.header("Configuration")
        
        optimization_type = st.radio(
            "Type d'optimisation",
            ["Minimiser", "Maximiser"],
            index=0
        )
        is_maximize = (optimization_type == "Maximiser")
        
        M_value = st.number_input(
            "Valeur de M (Big M)",
            min_value=1000.0,
            max_value=1e10,
            value=1e6,
            format="%.0e"
        )
        
        show_iterations = st.checkbox("Afficher les itérations", value=True)
        
        st.markdown("---")
        st.markdown("### Guide rapide")
        st.markdown("""
        1. Entrez les coefficients de la fonction objectif
        2. Ajoutez vos contraintes une par une
        3. Cliquez sur 'Résoudre'
        4. Explorez les résultats !
        """)
    
    # Section principale
    col1, col2 = st.columns([2, 1])
    
    with col1:
        st.header("Fonction Objectif")
        obj_input = st.text_input(
            f"{optimization_type} Z =",
            value="60, 80",
            help="Entrez les coefficients séparés par des virgules (ex: 3, 5, 2)"
        )
        
        obj_coeffs = parse_coefficients(obj_input)
        if obj_coeffs:
            st.success(f"✓ {len(obj_coeffs)} variables détectées: " + 
                      ", ".join([f"{c}x{i+1}" for i, c in enumerate(obj_coeffs)]))
        else:
            st.error("❌ Format invalide pour la fonction objectif")
    
    with col2:
        st.header("Nombre de variables")
        if obj_coeffs:
            st.metric("Variables", len(obj_coeffs))
        else:
            st.metric("Variables", "—")
    
    st.markdown("---")
    
    # Section des contraintes
    st.header("Contraintes")
    
    # Initialiser le state pour les contraintes
    if 'constraints' not in st.session_state:
        st.session_state.constraints = [
            {'coeffs_text': '1, 1', 'type': '=', 'rhs': 8.0},
            {'coeffs_text': '1, 0', 'type': '>=', 'rhs': 2.0},
            {'coeffs_text': '0, 1', 'type': '>=', 'rhs': 3.0},
            {'coeffs_text': '1, 0', 'type': '<=', 'rhs': 5.0},
            {'coeffs_text': '0, 1', 'type': '<=', 'rhs': 6.0},
        ]
    
    # Afficher les contraintes existantes
    constraints_to_remove = []
    for i, constraint in enumerate(st.session_state.constraints):
        col1, col2, col3, col4, col5 = st.columns([3, 1, 1, 1, 1])
        
        with col1:
            constraint['coeffs_text'] = st.text_input(
                f"Coefficients C{i+1}",
                value=constraint['coeffs_text'],
                key=f"coeffs_{i}",
                help="Coefficients séparés par des virgules"
            )
        
        with col2:
            constraint['type'] = st.selectbox(
                f"Type C{i+1}",
                ['<=', '>=', '='],
                index=['<=', '>=', '='].index(constraint['type']),
                key=f"type_{i}"
            )
        
        with col3:
            constraint['rhs'] = st.number_input(
                f"RHS C{i+1}",
                value=float(constraint['rhs']),
                key=f"rhs_{i}"
            )
        
        with col4:
            st.write("")
            st.write("")
            coeffs = parse_coefficients(constraint['coeffs_text'])
            if coeffs:
                st.success("✓")
            else:
                st.error("✗")
        
        with col5:
            st.write("")
            st.write("")
            if st.button("🗑️", key=f"remove_{i}"):
                constraints_to_remove.append(i)
    
    # Supprimer les contraintes marquées
    for i in reversed(constraints_to_remove):
        st.session_state.constraints.pop(i)
        st.rerun()
    
    # Bouton pour ajouter une nouvelle contrainte
    col1, col2 = st.columns([1, 4])
    with col1:
        if st.button("➕ Ajouter une contrainte"):
            if obj_coeffs:
                default_coeffs = ", ".join(['0'] * len(obj_coeffs))
                st.session_state.constraints.append({
                    'coeffs_text': default_coeffs,
                    'type': '<=',
                    'rhs': 0.0
                })
                st.rerun()
    
    st.markdown("---")
    
    # Bouton de résolution
    if st.button("Résoudre le problème", type="secondary", use_container_width=True):
        # Parser toutes les contraintes
        constraints_data = []
        for constraint in st.session_state.constraints:
            coeffs = parse_coefficients(constraint['coeffs_text'])
            if coeffs:
                constraints_data.append({
                    'coeffs': coeffs,
                    'type': constraint['type'],
                    'rhs': constraint['rhs']
                })
        
        # Valider les entrées
        valid, error_msg = validate_inputs(obj_coeffs, constraints_data)
        
        if not valid:
            st.error(f"❌ {error_msg}")
        else:
            # Préparer les données pour le solveur
            c = np.array(obj_coeffs)
            A = np.array([c['coeffs'] for c in constraints_data])
            b = np.array([c['rhs'] for c in constraints_data])
            types = [c['type'] for c in constraints_data]
            
            # Résoudre
            with st.spinner('🔄 Résolution en cours...'):
                try:
                    solver = BigM_Solver(c, A, b, types, M=M_value, is_maximize=is_maximize)
                    result = solver.solve(store_history=show_iterations)
                    
                    # Afficher les résultats
                    st.success("Résolution terminée !")
                    st.markdown("---")
                    
                    if result['status'] == 'optimal':
                        # Solution optimale
                        st.header("Solution Optimale")
                        
                        col1, col2, col3 = st.columns(3)
                        with col1:
                            st.metric(
                                "Valeur Optimale",
                                f"{result['objective_value']:.6f}"
                            )
                        with col2:
                            st.metric(
                                "Itérations",
                                result['iterations']
                            )
                        with col3:
                            st.metric(
                                "Statut",
                                "Optimal"
                            )
                        
                        # Variables de décision
                        st.subheader("Variables de décision")
                        solution_data = {
                            'Variable': [f'x{i+1}' for i in range(len(result['x']))],
                            'Valeur': result['x']
                        }
                        st.dataframe(
                            pd.DataFrame(solution_data),
                            use_container_width=True,
                            hide_index=True
                        )
                        
                        # Tableau final
                        st.subheader("Tableau Final")
                        st.dataframe(
                            result['final_tableau'].style.format("{:.4f}"),
                            use_container_width=True
                        )
                        
                        # Historique des itérations
                        if show_iterations and hasattr(solver, 'tableaux_history'):
                            st.markdown("---")
                            st.header("Historique des Itérations")
                            
                            for record in solver.tableaux_history:
                                with st.expander(f"Itération {record['iteration']}", 
                                               expanded=(record['iteration'] == 0)):
                                    if record['iteration'] > 0:
                                        st.info(f"🔄 {record['entering']} entre, {record['leaving']} sort")
                                    st.dataframe(
                                        record['tableau'].style.format("{:.4f}"),
                                        use_container_width=True
                                    )
                    
                    elif result['status'] == 'infeasible':
                        st.error("❌ Problème infaisable")
                        st.warning(result['message'])
                    
                    elif result['status'] == 'unbounded':
                        st.error("❌ Problème non borné")
                        st.warning(result['message'])
                    
                    else:
                        st.error(f"❌ Erreur: {result.get('message', 'Erreur inconnue')}")
                
                except Exception as e:
                    st.error(f"❌ Erreur lors de la résolution: {str(e)}")
                    st.exception(e)

if __name__ == "__main__":

    main()
