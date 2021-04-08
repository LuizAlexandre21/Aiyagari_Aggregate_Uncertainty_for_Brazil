
<h1> Aiyagari_Aggregate_Uncertainty_for_Brazil </h1>


# **Clonando o repositório**
~~~bash 
$git clone https://github.com/LuizAlexandre21/Aiyagari_Aggregate_Uncertainty_for_Brazil.git
~~~
 
# **Requisitos** 

## **Primeiros Passos no Julia** 

### **Pre-requisitos**
~~~
# Bibliotecas para julia

Distributions 
Plots 
DataFrames 
Random 
ForwardDiff
LinearAlgebra
Interpolations 
DataFrames
Optim 
NBInclude
~~~

### **Instalando Bibliotecas**
~~~julia 
using Pkg 
Pkg.add("Biblioteca")
Pkg.add("Plots")
~~~

# **Modelo de Aiyagari**
Para a execução do modelo de aiyagari no terminal utilizar o seguinte código:
~~~bash
$julia Aiyagari.jl
~~~
# **Importando Bibliotecas e Arquivos Externos**

## **Importando as Bibliotecas** 
~~~julia
using Distributions
using Plots
using DataFrames
using Random
using ForwardDiff
using LinearAlgebra
using Interpolations
using DataFrames
using Optim
using NBInclude
~~~

## **Importando os notebooks**
~~~julia
@nbinclude("utils.ipynb") 
@nbinclude("RBC.ipynb")  
@nbinclude("EGM.ipynb")   
@nbinclude("SteadyState.ipynb") 
@nbinclude("GenBKM.ipynb") 
~~~

##### Notebooks do professor Julian Pascal, com alteração nos parametros e nas funçoes dos juros e salários, no arquivo utils.ipynb


# **Resultados do modelo de aiyagari**

## **Estado Estácionario**

![](https://imgur.com/pmfPC6P.png)

## **Choques no Capital Humano**

![Imgur](https://i.imgur.com/bYhUjR4.png)

## Checando a Linearidade 
~~~julia
~~~
![Imgur](https://i.imgur.com/fVNSbSC.png)


# Contribuidores 

<table> 
    <tr>
    <td align="center"><a http://lattes.cnpq.br/9458204748985902><img src= https://avatars.githubusercontent.com/u/23129808?s=400&u=739133cee3ef3858b8a277dba0b7388145a503fa&v=4 width="100px" ><sub><b>Luiz Alexandre Moreira Barros</b></td>
    <td align="center"><a http://lattes.cnpq.br/6695247122647476 ><img src= http://servicosweb.cnpq.br/wspessoa/servletrecuperafoto?tipo=1&id=K8211460U6 width="100px" ><sub><b> Raphael Douglas de Freitas Lucena</b></td>
    </tr>
</table>
