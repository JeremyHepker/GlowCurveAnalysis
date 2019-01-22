#  Developer Notes

Random developer notes. 

## Known materials database implementation
Possible approaches -
    1) Json file - clean and efficient, however commercial modifications could be challeneging. 
    2) Hash table - Efficient, common, and scaleable.
    3) Vector - Inefficent but very common and scaleable.
Possible further development - 
    1) Live server maintnance - Rapid growth, development, and colabrative improvment. 
        -Would be best for further reseach and development with smart mode. 
        -Would not be free(options). 
        -Would be multi-fauceted for use in other related projects not just tld reader. 
    2) ...
    
## First order kinetics implmentation

I(T)=Im*b^(b/(b-1))*expv*((b-1)*(1-xa)*(T/Tm)^2*expv+Zm)^(-b/(b-1)) <br />
xa=2*k*T/E
xb=2*k*Tm/E
expv=exp(E/(k*T)*(T-Tm)/Tm)
Zm=1+(b-1)*xb

where:
b is the kinetic parameter [lies between 1 and 2)]
I is the glow peak intensity
E the activation energy [in ev]
k the Boltzmann constant [in eV/k]
T the temperature in K with constant heating rate [in K/s]
Tm the temperature at maximum thermoluminescence intensity [in K]
Im the maximum intensity

## OTOR 

... Oh boy. 

## Efficiency Reports 

To date: N/A 

Optimization Ideas -
    1) Adapt material recognition in all three models, not just smart mode, for faster recognition and analysis. 
    2) ...


## Bin Commands 
./bin/gitupdate
