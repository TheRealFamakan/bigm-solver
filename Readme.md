# Solveur de Programmation Linéaire - Méthode Big M

Application web interactive pour résoudre des problèmes de programmation linéaire en utilisant la méthode Big M avec l'algorithme du simplexe.

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![Streamlit](https://img.shields.io/badge/Streamlit-1.28+-red.svg)
![NumPy](https://img.shields.io/badge/NumPy-1.24+-green.svg)
![Pandas](https://img.shields.io/badge/Pandas-2.0+-yellow.svg)

## Fonctionnalités

- **Minimisation et Maximisation** - Résout les deux types de problèmes
- **Contraintes mixtes** - Support complet de `<=`, `>=` et `=`
- **Interface interactive** - Ajout/suppression dynamique de contraintes
- **Visualisation détaillée** - Affichage de toutes les itérations du simplexe
- **Validation en temps réel** - Vérification instantanée des entrées
- **Résultats complets** - Solution optimale, valeurs des variables, tableaux

## Installation rapide

### Prérequis

- Python 3.8 ou supérieur
- pip (gestionnaire de paquets Python)

### Étapes d'installation

1. **Cloner ou télécharger le projet**
   ```bash
   git clone https://github.com/therealfamakan/bigm-solver.git
   cd bigm-solver
   ```

2. **Installer les dépendances**
   ```bash
   pip install -r requirements.txt
   ```

3. **Lancer l'application**
   ```bash
   streamlit run app.py
   ```

4. **Ouvrir dans le navigateur**
   - L'application s'ouvre automatiquement sur `http://localhost:8501`
   - Si ce n'est pas le cas, ouvrez manuellement ce lien dans votre navigateur

## Guide d'utilisation

### Exemple simple

**Problème :** Minimiser `60x₁ + 80x₂`

**Sous contraintes :**
- x₁ + x₂ = 8
- x₁ ≥ 2
- x₂ ≥ 3
- x₁ ≤ 5
- x₂ ≤ 6

**Étapes dans l'application :**

1. **Configuration** (sidebar)
   - Sélectionnez "Minimiser"
   - Gardez M = 1e6 (valeur par défaut)

2. **Fonction objectif**
   - Entrez : `60, 80`

3. **Contraintes**
   - Contrainte 1 : `1, 1` | `=` | `8`
   - Contrainte 2 : `1, 0` | `>=` | `2`
   - Contrainte 3 : `0, 1` | `>=` | `3`
   - Contrainte 4 : `1, 0` | `<=` | `5`
   - Contrainte 5 : `0, 1` | `<=` | `6`

4. **Résoudre**
   - Cliquez sur le bouton "Résoudre le problème"

5. **Résultats**
   - Solution optimale : x₁ = 2, x₂ = 3
   - Valeur optimale : Z = 360

## Structure du projet

```
bigm-solver/
│
├── app.py                  # Application Streamlit complète
├── requirements.txt        # Dépendances Python
├── README.md              # Ce fichier
└── .gitignore             # Fichiers à ignorer par Git (optionnel)
```

## 🔧 Configuration avancée

### Paramètres ajustables

Dans la sidebar de l'application :

- **Type d'optimisation** : Minimiser ou Maximiser
- **Valeur de M** : Pénalité Big M (défaut: 1e6)
- **Afficher les itérations** : Voir tous les tableaux intermédiaires

### Format des entrées

**Fonction objectif :**
```
3, 5, 2         # Pour 3x₁ + 5x₂ + 2x₃
```

**Contraintes :**
- Coefficients : `2, 3, 1` (séparés par des virgules)
- Type : `<=`, `>=` ou `=`
- RHS : Nombre décimal (ex: 10.5)

## Tests

### Test 1 : Problème standard
```
Minimiser: 3x₁ + 2x₂
Contraintes:
  2x₁ + x₂ ≤ 10
  x₁ + 2x₂ ≤ 8
Solution attendue: x₁ = 4, x₂ = 2, Z = 16
```

### Test 2 : Problème avec contraintes >=
```
Minimiser: x₁ + 2x₂
Contraintes:
  x₁ + x₂ ≥ 3
  2x₁ + x₂ ≥ 4
Solution attendue: x₁ = 1, x₂ = 2, Z = 5
```

## Technologies utilisées

- **[Streamlit](https://streamlit.io/)** - Framework web interactif
- **[NumPy](https://numpy.org/)** - Calculs matriciels et algèbre linéaire
- **[Pandas](https://pandas.pydata.org/)** - Manipulation et affichage des données

## Méthode Big M

La méthode Big M est une technique pour résoudre des problèmes de programmation linéaire avec des contraintes mixtes (`<=`, `>=`, `=`). Elle fonctionne en :

1. **Ajoutant des variables artificielles** pour les contraintes `>=` et `=`
2. **Pénalisant ces variables** avec un grand coefficient M dans la fonction objectif
3. **Utilisant l'algorithme du simplexe** pour trouver la solution optimale
4. **Garantissant** que les variables artificielles sortent de la base

## Limitations

- L'application est conçue pour des problèmes de **taille raisonnable** (< 20 variables, < 30 contraintes)
- Tous les coefficients du membre de droite (RHS) doivent être **non-négatifs**
- Les variables sont implicitement **non-négatives** (xᵢ ≥ 0)

## Dépannage

### Erreur : "Module not found"
```bash
pip install --upgrade streamlit numpy pandas
```

### L'application ne démarre pas
- Vérifiez que Python 3.8+ est installé : `python --version`
- Vérifiez que les dépendances sont installées : `pip list`

### Erreur de calcul
- Vérifiez que tous les coefficients sont numériques
- Vérifiez que le nombre de coefficients correspond au nombre de variables
- Assurez-vous que les valeurs RHS sont positives

### Port déjà utilisé
```bash
streamlit run app.py --server.port 8502
```

## Contact & Support

Pour toute question ou suggestion :
- Email : camarafamakan2@gmail.com
- Issues : [GitHub Issues](https://github.com/TheRealFamakan/bigm-solver)

## Licence

Ce projet est sous licence MIT. Voir le fichier [LICENSE](LICENSE) pour plus de détails.

## Remerciements

- Méthode Big M développée par George Dantzig
- Interface construite avec [Streamlit](https://streamlit.io/)
- Algorithmes implémentés en Python avec NumPy et Pandas
- Mes collegues de travail: Nawal Ait-Tami, Hiba El Hamdani et Chaimae EL Mounjali.

## Mises à jour

### Version 1.0.0 (2024)
- Version initiale
- Support complet de la méthode Big M
- Interface Streamlit interactive
- Validation des entrées en temps réel
- Affichage des itérations

---


**Made with ❤️ and Python by The best Team IID2 Ensa Khouribga**
