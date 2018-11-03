/* résoudre g^x = h mod q */

/* décomposition en facteurs de la base B */
decomposition(x, B) = {
     X = lift(x);
     e = [];
     for(i = 1, #B,
          e = concat(e,0);
     );
     for(i = 1, #B,
          while(X%B[i] == 0,
               X \= B[i];
               e[i]++;
          );
     );
     if(X == 1,
          return(e);
     ,
          return([]);
     );
};

pivot(M) = {
     my(modulo, pivot);
     modulo = M[1,1].mod;
     for(i = 1, #M - 1,
          /* on recherche le premier pivot inversible */
          cpt = i;
          while(cpt <= #M~ && gcd(M[cpt,i], modulo) != 1,
               cpt++;
          );
          /* si l'élément trouvé n'est pas sur la ie ligne, on permute */
          if(cpt > #M~,
               error("Matrice non inversible");
          ,
          );
          if(i != cpt,
               tmp = M[cpt,];
               M[cpt,] = M[i,];
               M[i,] = tmp;
          ,
          );
          /* on réduit toutes les lignes en dessous */
          pivot = M[i,i]^-1;
          M[i,] *= pivot;
          for(j = i + 1, #M~,
               M[j,] -= M[j,i]*M[i,];
          );
     );
     /* ici, la matrice est triangulaire inférieure             *
      * on recopie uniquement les éléments qui nous intéressent */
     N = M[1,];
     for(i = 2, #M - 1,
          N = matconcat([N;M[i,]]);
     );
     /* on met des 0 au dessus de la diagonale */
     forstep(i = #M - 1, 1, -1,
          for(j = 1, i - 1,
               N[j,] -= N[j,i]*N[i,];
          );
     );
     return(N);
};

LogDiscret(h, g /* mod q */, a) = {
     if(type(h) != "t_INT", return([]));
     if(type(g) != "t_INTMOD", return([]));
     if(type(a) != "t_INT", return([]));
     if('a == a,
          a = 2;
     ,
          if(a > 2,
               a = floor(a);
          ,
               a = 2;
          );
     );
     \\my(B, relations, puissance, size_relations, v, k, ord, f, s, x, logBase, r);
     ord = g.mod - 1;
     if(verification(g) == 0,
          return([]);
     ,
     );
     /* si on cherche la puissance de g qui envoie sur 1, la seule possible est *
      * celle correspondant à l'ordre du groupe                                 */
     if(h == 1,
          return(ord);
     ,
     );

     /* calcul du plus grand nombre de la base B */
     B_max = exp(sqrt(log(ord)*log(log(ord))/2));

     r = primepi(B_max);

     /* primes(i) renvoie un tableau composé des i premiers nombres premiers  *
      * primepi(i) renvoie le nombre de nombres premiers inférieurs/égaux à i */
     B = primes(r);

     /* création d'une matrice vide dans laquelle *
      * les relations indépendantes sont stockées */
     relations = [];

     nbr_relations = 0;

     /* remplissage des matrices */
     while(nbr_relations < a*r,
          k = random(ord) + 1;
          v = decomposition(g^k, B);
          if(v != 0 && gcd(gcd(v),ord) == 1,
               nbr_relations++;
               v = concat(v, k);
               if(nbr_relations == 1,
                    relations = v;
               ,
                    relations = matconcat([relations; v]);
               );
          ,
          );
     );
     relations = Mod(relations, ord);
     A = pivot(relations);
     logBase = lift(A[,#A]);
     s = 1;

     /* logBase correspond au logarithme discret de chaque élément de B (en base g)*/
     while(decomposition(g^s*h, B) == [],
          s++;
     );
     f = Mod(decomposition(g^s*h, B), ord);
     x = lift(f*logBase) - s;
     return(x);
};

/* g est-il générateur du groupe d'ordre ord et g.mod est-il premier ? */
verification(g) = {
     if(isprime(g.mod),
          ord = g.mod - 1;
          facteur = factor(ord);
          for(i = 1, #mattranspose(facteur),
               if(g^(ord/facteur[i,1]) == 1,
                    print("g n'est pas un générateur !");
                    return(0);
               ,
               );
          );
     ,
          print("Le générateur doit être construit modulo un nombre premier !");
          return(0);
     );
     return(1);
};
