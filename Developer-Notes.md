#  Developer Notes

Random developer notes. 

## Known materials database implementation
Possible approaches -<br />
    1) Json file - clean and efficient, however commercial modifications could be challeneging. <br />
    2) Hash table - Efficient, common, and scaleable.<br />
    3) Vector - Inefficent but very common and scaleable.<br />
Possible further development - <br />
    1) Live server maintnance - Rapid growth, development, and colabrative improvment. <br />
        -Would be best for further reseach and development with smart mode. <br />
        -Would not be free(options). <br />
        -Would be multi-fauceted for use in other related projects not just tld reader. <br />
    2) ...<br />
    
## First order kinetics implmentation

I(T)=Im*b^(b/(b-1))*expv*((b-1)*(1-xa)*(T/Tm)^2*expv+Zm)^(-b/(b-1)) <br />
xa=2*k*T/E<br />
xb=2*k*Tm/E<br />
expv=exp(E/(k*T)*(T-Tm)/Tm)<br />
Zm=1+(b-1)*xb<br />

where:<br />
b is the kinetic parameter [lies between 1 and 2)]<br />
I is the glow peak intensity<br />
E the activation energy [in ev]<br />
k the Boltzmann constant [in eV/k]<br />
T the temperature in K with constant heating rate [in K/s]<br />
Tm the temperature at maximum thermoluminescence intensity [in K]<br />
Im the maximum intensity<br />

## OTOR 

... Oh boy. 

## Efficiency Reports 

To date: N/A 

Optimization Ideas -<br />
    1) Adapt material recognition in all three models, not just smart mode, for faster recognition and analysis. <br />
    2) ...<br />


## Bin Commands 
./bin/gitupdate

## Resources
http://www.iram.fr/IRAMFR/GILDAS/doc/html/map-html/node37.html
