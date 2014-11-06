# Demo of the VineCPP toolbox

This demo should illustrate how the VineCPP toolbox can be used to work with vine copulas. This demo starts with the simulation of data from vine copula models. The focus in this part is on the simulation from vine copulas, where some of the conditional copulas are copulas for which the parameter is a function of the conditioning variable. The second part shows how one can use the functions of the VineCPP toolbox to select simplified vine copula models. Estimation is illustrated in the third part of this demo. The fourth part of the demo shows how one can use the toolbox to test for simplified vine copulas. It is shown how the simplifying assumption can be tested with an sequential approach, based on a statistical tests on vectorial independencies.

## Simulating possibly non-simplified vine copulas

In the following three different data generating processes will be used to present the functioning of the toolbox. The first example is the five-dimensional Clayton copula. From this copula, one can simulate in two different ways, either  by using the laplace transform to simulate directly from the five-dimensional Archimedean copula or by simulating from the corresponding C-Vine representation (like it will be done in the following).

<pre class="codeinput"><span class="comment">% Choose the number of dimensions for the C-Vine.</span>
dimension = 5;

<span class="comment">% Specify the pair-copula families.</span>
families = ones(1,dimension*(dimension-1)/2).*7;

<span class="comment">% Select the structure of the C-Vine (i.e., the order of the nodes).</span>
structure = 1:dimension;

<span class="comment">% Define all copulas being bart of the C-Vine as pair-copulas (i.e.,</span>
<span class="comment">% unconditional bivariate copulas)</span>
simplified = ones(1,(dimension-1)*(dimension-2)/2);

<span class="comment">% Specify the vine copula type (yet C-Vine is the only possible choice)</span>
type = <span class="string">'C-Vine'</span>;

<span class="comment">% Give the parameters for all unconditional and partial copulas being</span>
<span class="comment">% building blocks of the C-Vine copula.</span>
theta = 3;
parameters = repmat(theta,1,dimension-1);
<span class="keyword">for</span> i=2:dimension-1
    parameters = [parameters , repmat(theta/((i-1)*theta+1),1,dimension-i)];
<span class="keyword">end</span>
parameters;

<span class="comment">% Use all the pre-defined properties to construct an object of the</span>
<span class="comment">% VineCopula class</span>
VineCopulaObject1 = VineCopula(dimension,type,simplified,structure,families,parameters)

U = Sim(VineCopulaObject1,500);
figure(<span class="string">'Units'</span>,<span class="string">'normalized'</span>,<span class="string">'Position'</span>,[0.2 0.2 0.6 0.8],<span class="string">'PaperPositionMode'</span>,<span class="string">'auto'</span>);
plotmatrix(U)
</pre><pre class="codeoutput">
VineCopulaObject1 = 

VineCopula

   Vine copula properties:
          Type: 'C-Vine'
     Dimension: 5
     Structure: [1 2 3 4 5]
    simplified: [1 1 1 1 1 1]
        MaxLLs: []

   Building blocks tree by tree:

   1. Tree
    Pair_Copula: 'C_1,2: Clayton copula (theta = 3)'

    Pair_Copula: 'C_1,3: Clayton copula (theta = 3)'

    Pair_Copula: 'C_1,4: Clayton copula (theta = 3)'

    Pair_Copula: 'C_1,5: Clayton copula (theta = 3)'

   2. Tree
    Pair_Copula: 'C_2,3;1: Clayton copula (theta = 0.75)'

    Pair_Copula: 'C_2,4;1: Clayton copula (theta = 0.75)'

    Pair_Copula: 'C_2,5;1: Clayton copula (theta = 0.75)'

   3. Tree
    Pair_Copula: 'C_3,4;1,2: Clayton copula (theta = 0.4286)'

    Pair_Copula: 'C_3,5;1,2: Clayton copula (theta = 0.4286)'

   4. Tree
    Pair_Copula: 'C_4,5;1,2,3: Clayton copula (theta = 0.3)'

</pre><img vspace="5" hspace="5" src="demoVineCPP_01.png" alt="">

The second example is the same vine copula (the C-Vine representation of the five-dimensional Clayton copula), where the last partial copula (i.e., the copula C_45|123 is substituted by a Frank copula with functional parameter theta(x_1) = (4x_1-2)^3. The described C-Vine is constructed as a member of the VineCopula class, which form the central class of the whole toolbox.

<pre class="codeinput"><span class="comment">% Choose the number of dimensions for the C-Vine.</span>
dimension = 5;

<span class="comment">% Specify the pair-copula families</span>
families(end) = 9;

<span class="comment">% Select the structure of the C-Vine (i.e., the order of the nodes).</span>
structure = 1:dimension;

<span class="comment">% The last copula is a conditional copula</span>
simplified(end) = 0;

<span class="comment">% Specify the vine copula type (yet C-Vine is the only possible choice)</span>
type = <span class="string">'C-Vine'</span>;

<span class="comment">% Give the parameters for all pair-copulas copulas being</span>
<span class="comment">% building blocks of the C-Vine copula.</span>
parameters(end) = [];

<span class="comment">% Define the functional parameter as a cell of function-handles.</span>
ParamFunctional = {@(x) (4*x(:,1)-2).^3};

<span class="comment">% Use all the pre-defined properties to construct an object of the</span>
<span class="comment">% VineCopula class</span>
VineCopulaObject2 = VineCopula(dimension,type,simplified,structure,families,parameters,ParamFunctional)
</pre><pre class="codeoutput">
VineCopulaObject2 = 

VineCopula

   Vine copula properties:
          Type: 'C-Vine'
     Dimension: 5
     Structure: [1 2 3 4 5]
    simplified: [1 1 1 1 1 0]
        MaxLLs: []

   Building blocks tree by tree:

   1. Tree
    Pair_Copula: 'C_1,2: Clayton copula (theta = 3)'

    Pair_Copula: 'C_1,3: Clayton copula (theta = 3)'

    Pair_Copula: 'C_1,4: Clayton copula (theta = 3)'

    Pair_Copula: 'C_1,5: Clayton copula (theta = 3)'

   2. Tree
    Pair_Copula: 'C_2,3;1: Clayton copula (theta = 0.75)'

    Pair_Copula: 'C_2,4;1: Clayton copula (theta = 0.75)'

    Pair_Copula: 'C_2,5;1: Clayton copula (theta = 0.75)'

   3. Tree
    Pair_Copula: 'C_3,4;1,2: Clayton copula (theta = 0.4286)'

    Pair_Copula: 'C_3,5;1,2: Clayton copula (theta = 0.4286)'

   4. Tree
    Cond_Copula: 'C_4,5|1,2,3: Frank copula (@(x)(4*x(:,1)-2).^3)'

</pre><p>Having defined the object of the VineCopula class, one can simulate from it by using the method Sim.</p><pre class="codeinput"><span class="comment">% Use the method Sim for objects from the VineCopula class, to simulate</span>
<span class="comment">% from the specified non-simplified C-Vine copula.</span>
V = Sim(VineCopulaObject2,500);
figure(<span class="string">'Units'</span>,<span class="string">'normalized'</span>,<span class="string">'Position'</span>,[0.2 0.2 0.6 0.8],<span class="string">'PaperPositionMode'</span>,<span class="string">'auto'</span>);
plotmatrix(V)
</pre><img vspace="5" hspace="5" src="demoVineCPP_02.png" alt=""> <p>The third example is a arbitrary combination of partial copulas and conditional copulas with functional parameter as building block of a five-dimensional non-simplified C-Vine copula.</p><pre class="codeinput"><span class="comment">% Construct the object of the VineCopula class and apply the method Sim to</span>
<span class="comment">% it.</span>
VineCopulaObject3 = VineCopula(5,<span class="string">'C-Vine'</span>,<span class="keyword">...</span>
    [0,0,1,1,1,1],<span class="keyword">...</span>
    [1 2 3 4 5],<span class="keyword">...</span>
    {<span class="string">'Clayton'</span>,<span class="string">'Gumbel'</span>,<span class="string">'t'</span>,<span class="string">'Indep'</span>,<span class="string">'Frank'</span>,<span class="string">'BB1'</span>,<span class="string">'Gaussian'</span>,<span class="string">'Indep'</span>,<span class="string">'Tawn'</span>,<span class="string">'Indep'</span>; [],[],[],[],[],<span class="string">'r90'</span>,[],[],<span class="string">'r270'</span>,[]},<span class="keyword">...</span>
    [1.2, 3,-0.9,1.6,0.7,5,0.1,0.9],<span class="keyword">...</span>
    {@(u) (4.*u-2).^3, @(u) (u+1), @(u) (u+1)})

W = Sim(VineCopulaObject3,500);
figure(<span class="string">'Units'</span>,<span class="string">'normalized'</span>,<span class="string">'Position'</span>,[0.2 0.2 0.6 0.8],<span class="string">'PaperPositionMode'</span>,<span class="string">'auto'</span>);
plotmatrix(W)
</pre><pre class="codeoutput">
VineCopulaObject3 = 

VineCopula

   Vine copula properties:
          Type: 'C-Vine'
     Dimension: 5
     Structure: [1 2 3 4 5]
    simplified: [0 0 1 1 1 1]
        MaxLLs: []

   Building blocks tree by tree:

   1. Tree
    Pair_Copula: 'C_1,2: Clayton copula (theta = 1.2)'

    Pair_Copula: 'C_1,3: Gumbel copula (theta = 3)'

    Pair_Copula: 'C_1,4: t copula (theta = (-0.9,1.6))'

    Pair_Copula: 'C_1,5: Indep copula'

   2. Tree
    Cond_Copula: 'C_2,3|1: Frank copula (@(u)(4.*u-2).^3)'

    Cond_Copula: 'C_2,4|1: BB1 copula (90&deg;, @(u)(u+1) and @(u)(u+1))'

    Pair_Copula: 'C_2,5;1: Gaussian copula (theta = 0.7)'

   3. Tree
    Pair_Copula: 'C_3,4;1,2: Indep copula'

    Pair_Copula: 'C_3,5;1,2: Tawn copula (270&deg;, theta = (5,0.1,0.9))'

   4. Tree
    Pair_Copula: 'C_4,5;1,2,3: Indep copula'

</pre><img vspace="5" hspace="5" src="demoVineCPP_03.png" alt="">

## Selecting simplified vine copula models

In a next step, the method StructureSelect for objects from the class VineCopula should be used for selecting simplified vine copulas, which are then used as an approximation to the overall distributions. The nodes of the C-vine are chosen in a way that in each tree, the root (i.e. the node, which is connected by a copula to all other nodes) is the variable which has maximal dependence with all other variables. The maximal dependence is found by choosing the variable which has the maximal column sum in the matrix of absolute empirical Kendall&#8217;s &#964; (cf. Schepsmeier, St&ouml;ber, and Brechmann (2013) for an R-function (RVineStructureSelect) of exactly the same procedure and Czado, Schepsmeier, and Min (2012, p. 240) for the theoretical background of the approach). Furthermore, the copula families are chosen according to the AIC criterion and for each pair- copula an independence test is performed (cf. Schepsmeier, St&ouml;ber, and Brechmann (2013) and Brechmann and Schepsmeier (2013) for R-functions (RVineStructureSelect / RVineCopSelect / CDVineCopSelect) and Genest and Favre (2007, p. 351) for the independence test).

<pre class="codeinput">VineCopulaSel3 = StructureSelect(VineCopula(5,<span class="string">'C-Vine'</span>,1),W);
</pre><p>For the first and second example, the structure is selected manually. For the second example with the Frank copula with varying parameter a Independence copula is selected as partial copula to obtain a approximation of the overall distribution, which is a simplified C-Vine copula.</p><pre class="codeinput">VineCopulaSel1 = VineCopula(5,<span class="string">'C-Vine'</span>,<span class="keyword">...</span>
    true,<span class="keyword">...</span>
    [1 2 3 4 5],<span class="keyword">...</span>
    repmat({<span class="string">'Clayton'</span>},1,5*(5-1)/2));

VineCopulaSel2 = VineCopula(5,<span class="string">'C-Vine'</span>,<span class="keyword">...</span>
    true,<span class="keyword">...</span>
    [1 2 3 4 5],<span class="keyword">...</span>
    [repmat({<span class="string">'Clayton'</span>},1,5*(5-1)/2-1) {<span class="string">'Indep'</span>}]);</pre>

## Estimating simplified vine copula models

The next part of this demo demonstrates how the Fit method of the VineCopula class can be used to estimate simplified C-Vine copulas. In the following we want to estimate the selected models from the last section to the simulated data. As estimation method one can choose between a joint estimation and a sequential estimation approach. In the case of a joint estimation, the sequential approach is applied first to obtain starting values for the numerically optimization of the joint / overall log-likelihood. The default method is the joint estimation.

<pre class="codeinput">VineCopulaHat1 = Fit(VineCopulaSel1,U);

VineCopulaHat2 = Fit(VineCopulaSel2,V);

VineCopulaHat3 = Fit(VineCopulaSel3,W);</pre>

## Testing the simplified assumption for vine copulas

Now the simplified assumption should be tested in a sequential manner. Note that the different tests, used in the following are not tests on the simplifying assumption in general but every conditional copula being part of the C-Vine is tested against a partial copula. This can be seen as a sequential testing procedure on the simplifying assumption as the simplifying assumtion is equivalent to the assumption that every conditional copula of a pair-copula construction (PCC) is a partial copula. For the interpretation of the test results, one has to keep in mind that everytime a test on the partial copula is performed, the assumption that one is able to obtain observations from the conditional copula that should be tested is needed. This assumption includes that the pair-copulas (including the assumption of having partial copulas) in the lower trees of the vine copula are correctly specified.

<pre class="codeinput"><span class="comment">% Apply the tests to the data examples.</span>
[pVals1,TestStats1] = SeqTestOnSimplified(VineCopulaHat1,U)
</pre><pre class="codeoutput">
pVals1 =

    0.5750
    0.2100
    0.0610
    0.3970
    0.1650
    0.0100


TestStats1 =

    0.0206
    0.0325
    0.0458
    0.0252
    0.0325
    0.0528

</pre><pre class="codeinput">[pVals2,TestStats2] = SeqTestOnSimplified(VineCopulaHat2,V)
</pre><pre class="codeoutput">
pVals2 =

    0.9680
    0.8340
    0.5640
    0.2480
    0.5440
    0.1460


TestStats2 =

    0.0124
    0.0155
    0.0212
    0.0293
    0.0209
    0.0286

</pre><pre class="codeinput">[pVals3,TestStats3] = SeqTestOnSimplified(VineCopulaHat3,W)
</pre><pre class="codeoutput">
pVals3 =

    0.0520
    0.0070
         0
    0.0410
    0.3390
    0.1640


TestStats3 =

    0.0353
    0.0473
    0.1153
    0.0246
    0.0132
    0.0094</pre>

## References

[1]  Brechmann, E. C. and U. Schepsmeier (2013), "Modeling Dependence      with C- and D-Vine Copulas: The R-Package CDVine", Journal of      Statistical Software 52(3), R package version 1.1-13, pp. 1-27, url: http://CRAN.R-project.org/package=CDVine.

[2]  Czado, C., U. Schepsmeier, and A. Min (2012), "Maximum likelihood      estimation of mixed C-vines with application to exchange rates",      Statistical Modelling 12(3), pp. 229-255.

[3]  Genest, C. and A. Favre (2007), "Everything You Always Wanted to      Know about Copula Modeling but Were Afraid to Ask", Journal of      Hydrologic Engineering 12(4), pp. 347-368.

[4]  Schepsmeier, U., J. St&ouml;ber, and E. C. Brechmann (2013), VineCopula:      Statistical inference of vine copulas, R package version 1.2, url:     http://CRAN.R-project.org/package=VineCopula.

---

Author: Malte Kurz

---


<p class="footer"><br><a href="http://www.mathworks.com/products/matlab/">Published with MATLAB&reg; R2013b</a><br></p>