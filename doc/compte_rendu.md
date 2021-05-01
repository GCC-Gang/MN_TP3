# MN_TP3
### FODOR Gergely, PELISSE VERDOUX Cyprien , VIALLET Camille

## Parallélisation

Les fonctions BLAS sont souvent utilisées par de divers programmes et avec des jeux de donnés très grands. Il est important 
de réduire un maximum le coup de ces fonctions. Dans le TP2, nous faisons attention à cela en réduisant la complexité 
des algorithmes implémentés ou encore en utilisant des variables de type `register` afin de minimiser les lectures en mémoire.

Dans ce TP, nous continuerons à augmenter les performances des fonctions BLAS mais cette fois en utilisant les fonctionnalités
offertes par le matériel. Nous allons tout d'abord répartir les calculs sur plusieurs threads, utilisant ainsi le plein
potentiel du processeur. Nous allons ensuite vectoriser nos donnés, nous permettant ainsi de traité plusieurs variables
à la fois. *Cette dernière partie étant optionnelle, seul une fonction n'a été vectorisé afin de voir les différences de 
performance.*

## OpenMP

La librairie OpenMP nous permet de faire du multi-threading sans avoir à utiliser les mutex et autres appels systèmes.
Grâce à des primitives de type `#pragma [OPTION]`, nous pouvons facilement répartir certaine partie de l'exécution 
sur plusieurs threads. Dans notre cas, ce sont les boucles sur des vecteurs et matrices de grande taille que nous allons 
paralléliser.

## Vectorisation

Comme spécifié dans l'introduction, cette partie du TP est optionnel. Néanmoins, nous avons voulu voir le potentiel
de la vectorisation et nous l'avons implémenté sur la fonction `scopy()`. Cette implémentation étant en conflit avec le 
reste du TP, je vous invite à regarder le code dans **test_copy_vect.c**.

La fonction `scopy()` copie un vecteur de float dans un autre vecteurs de float. Ici, notre but est de modifier le format
de donné afin que chaque case du tableau contiennent 4 float. On utilise le type `__m128`. L'initialisation du vecteur 
devient :

```c
#define VECSIZE 65536

typedef __m128 vfloat[VECSIZE];

void vector_init(vfloat V, float x) {
    unsigned int i;

    float tab[4] __attribute__ ((aligned(16))) = {x, x, x, x};

    for (i = 0; i < VECSIZE; i++)
        V[i] = _mm_load_ps(tab);

    return;
}
```

La fonction `scopy()` change à peine car la copie reste la même. Seul le prototype de la fonction devient :

```c
void mncblas_scopy(const int N, const __m128 *X, const int incX, __m128 *Y, const int incY);
```

Cette simple modification permet d'avoir une copie à **30 GB/s**, contre le **11 GB/s** avec seulement la parallélisation 
OpenMP. 

## Conclusion

Lors du semestre 5, nous avons vu comment rendre un algorithme plus performant en en réduisant la complexité. Nous avons 
notamment vu avec les tries la différence entre du O(n²) et du O(log n). L'algo avancé de ce semestre nous a montré qu'il 
est possible d'obtenir des gains de performances considérables lors de l'implementation de ces algorithmes. Le multi-threading 
et la vectorisation prennent avantage du matériel afin d'obtenir des programmes toujours plus rapides.