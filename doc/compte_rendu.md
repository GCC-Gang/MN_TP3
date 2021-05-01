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

