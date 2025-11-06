# Modélisation Géométrique - UE-INF2315M

Ce dépôt contient le code source et le rapport du projet réalisé dans le cadre de l'UE-INF2315M (Modélisation Géométrique) à l'Université Claude Bernard Lyon 1.

**Auteur :** Alexandre COTTIER

## 📜 À propos du projet

Ce projet explore différentes techniques de modélisation géométrique et de génération procédurale basées sur les surfaces implicites. Il est divisé en trois parties principales, plus une section bonus :

1.  **Modélisation à l'aide de Blobs** (Metaballs)
2.  **Modélisation à l'aide de Distances Signées** (SDF)
3.  **Génération Procédurale** (Érosion et Sphere Tracing)
4.  **Bonus :** Génération de terrain procédural (Heightmap avec bruit de Perlin)

---

## 1. 💧 Modélisation à l'aide de Blobs

Cette section implémente un champ scalaire implicite défini comme la somme des contributions de plusieurs primitives. La surface est définie par une valeur iso.

### Fonctionnalités

* **Primitives Points :** `AddPoint(P, R, w, c)` ajoute un "blob" centré en P, avec un rayon d'influence R et un poids W.
* **Primitives Segments :** `AddSegment(A, B, R, w, c)` crée un tube lisse entre deux points A et B.
* **Fonction d'influence :** Utilisation du **noyau borné de Wyvill** (cubique et lisse) pour une fusion douce entre les primitives.
* **Fusion :** Les blobs fusionnent naturellement lorsque leurs rayons d'influence se chevauchent, ou lorsqu'un point a un poids plus élevé, "tirant" la surface vers lui.

### Galerie (Blobs)

| 3 Points + 2 Segments | Fusion progressive (Rayon croissant) |
| :---: | :---: |
|  |  |
| Figure 1 du rapport | Figures 3, 4, 5, 6 du rapport |

---

## 2. 🧊 Modélisation à l'aide de Distances Signées (SDF)

Cette partie implémente un système de modélisation basé sur les champs de distances signés (SDF). Les formes sont définies par une fonction retournant la distance la plus courte à la surface (négative si à l'intérieur, positive à l'extérieur).

### Primitives (SDFNode)

* `SphereNode`
* `BoxNode`
* `CapsuleNode` (Segment gonflé)
* `TorusNode`

### Opérateurs

* **Union :** `UnionNode` - `min(A, B)`
* **Intersection :** `IntersectionNode` - `max(A, B)`
* **Différence :** `DifferenceNode` - `max(A, -B)`
* **Fusion Lisse :** `BlendNode` - `smooth-min(A, B, k)`

---

## 3. ⚙️ Génération Procédurale (Surfaces Implicites)

Cette section utilise les SDF pour des techniques de génération procédurale.

* **Sphere Tracing :** Implémentation d'un algorithme de *sphere tracing* (ou *ray marching*) pour calculer l'intersection d'un rayon avec la surface implicite (SDF). Le rayon avance par pas garantis de ne toucher aucune surface.
* **Érosion Procédurale :** Simulation d'usure ou d'érosion en combinant une forme de base (ex: une sphère) avec un ensemble de sphères d'érosion via des opérateurs (smoothmax, smoothmin, max). Cela permet de créer des cratères, des bosses et des formes complexes.

---

## 4. 🏔️ Bonus : Génération de Terrain Procédural

Une section supplémentaire, hors-sujet du cours principal, explore la génération de terrain procédural à l'aide de **bruit de Perlin** (ou ValueNoise).

* **Technique :** Le terrain est un maillage de type `heightmap`, où l'altitude $h(x, y)$ est calculée via une fonction de bruit.
* **Biomes :** La coloration du maillage dépend de l'altitude pour simuler différents biomes (eau, herbe, roche, neige).
* **Paramètres :** Le générateur est contrôlé par des paramètres tels que l'amplitude, la fréquence, les octaves, `warpStrength` (pour casser la symétrie) et `heightBias` (niveau de la mer).

---

## 📊 Analyse des Performances

Le projet inclut une analyse des performances des différentes techniques.

### Performances (Blobs)

* **Complexité :** Le temps de calcul montre une complexité en $O(n^3 \times m)$, où `n` est la résolution de la grille et `m` le nombre de primitives (blobs).
* **Scalabilité (Résolution) :** À nombre de primitives fixe, le temps de calcul est multiplié par **environ 8** à chaque doublement de `n` (ex: 314ms en n=32 vs 2200ms en n=64 ; 12127ms en n=128), ce qui est cohérent avec le $O(n^3)$.
* **Scalabilité (Primitives) :** Le temps de calcul croît de manière **linéaire** avec le nombre de primitives `m`.
* **Triangles :** Le nombre de triangles générés décroît à très haute densité de primitives, car les blobs fusionnent en un volume englobant plus simple.

### Performances (SDF)

Un benchmark sur 107 appels a été réalisé pour évaluer les opérations SDF.

| Opération | Ns / appel | M calls / s |
| :--- | :--- | :--- |
| **Sphere** | 106.5 ns | 9.39 M |
| **Torus** | 115.7 ns | 8.64 M |
| **Capsule** | 134.7 ns | 7.42 M |
| **Box** | 141.7 ns | 7.06 M |
| **Blend** (lissé) | 141.1 ns | 7.08 M |
| **Union** (brute) | 157.9 ns | 6.33 M |
| **Différence** | 158.6 ns | 6.30 M |
| **Intersection** | 161.7 ns | 6.18 M |

**Conclusion notable :** L'opération `Blend` (fusion lisse) est plus rapide que l'union simple (`min`/`max`). Cela s'explique par le fait que les CPU/GPU préfèrent les calculs arithmétiques continus (add, mul) aux branchements conditionnels (if, min, max) qui cassent le pipeline d'exécution et la vectorisation.

### Performances (Érosion)

* L'évaluation **incrémentale** (ajout d'impacts successifs) est environ **10-15% plus rapide** que l'évaluation **batch** (application de tous les impacts en une seule fois).

---
