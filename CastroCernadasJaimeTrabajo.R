###Trabajo R Jaime Castro Cernadas

###Cargamos los datos
datos <- read.table("trabajoR.txt",header=TRUE)

###Pregunta 1: Examinámos los datos. ¿Cuántas variables hay? ¿Cuántos tratamientos?
head(datos)
tail(datos)
summary(datos)
dim(datos)
str(datos)

### Hay 2 variables:  Variable 1 y Variable 2. Hay 5 tratamientos
# diferentes

#PREGUNTA 2: Haz un boxplot para nuestros datos. Uno para cada variable.
#Colorea a Variable 1 y a Variable 2 de forma diferente (guarda esos
#colores para las siguientes gráficas)

boxplot(Variable1~Tratamiento, main= "Boxplot Variable 1", col=("blue"), data = datos)
boxplot(Variable2~Tratamiento, main= "Boxplot Variable 2", col=("green"), data = datos)

#PREGUNTA3: Haz un gráfico de dispersión con las dos variables. Cada
#tratamiento debe de ir de un color distinto.
with(datos, plot(Variable1, Variable2, col=Tratamiento))

#PREGUNTA4: Ponle leyenda al gráfico del apartado anterior, en el margen 
#inferior derecho.

legend(x="bottomright", legend=c("Tto1","Tto2","Tto3","Tto4","Tto5"), fill= c("black","red","green","royalblue","blue"))

#PREGUNTA5:Haz un histograma para cada variable

hist(x=datos$Variable1, col = ("blue"), main = "Histograma Variable 1", xlab = "Pacientes",)
hist(x=datos$Variable2, col = ("green"), main = "Histograma Variable 2", xlab = "Pacientes",)

#PREGUNTA6: Haz un factor en la columna de tratamiento y guárdalo en una variable.

factor(datos$Tratamiento)
Factor <- factor(datos$Tratamiento)

#PREGUNTA7: Cálcula la media y desviación estandar para cada tratamiento.

aggregate(Variable1~Factor,datos,mean)
Factor Variable1
1      1      4.00
2      2      4.90
3      3      8.77
4      4     50.80
5      5     34.90

aggregate(Variable2~Factor,datos,mean)
Factor Variable2
1      1     0.510
2      2     1.300
3      3     5.310
4      4     8.730
5      5     9.018

aggregate(Variable1~Factor,datos,sd)
Factor Variable1
1      1  1.290133
2      2  3.754997
3      3  3.857475
4      4 11.113555
5      5  8.633912

aggregate(Variable2~Factor,datos,sd)
Factor Variable2
1      1 0.2884826
2      2 0.4189935
3      3 1.3568346
4      4 1.3333750
5      5 1.2146769

#PREGUNTA8: Averigua cuantos elementos tiene cada tratamiento
table(datos$Tratamiento, Factor)
Factor
1  2  3  4  5
1 10  0  0  0  0
2  0 10  0  0  0
3  0  0 10  0  0
4  0  0  0 10  0
5  0  0  0  0 10

table(Factor)
Factor
1  2  3  4  5 
10 10 10 10 10 
#PREGUNTA9: Extrae los datos para el tto1 y el tto4 Y guárdalos cada uno en una variable

> tto1 <- datos[1:10,]
> tto1
Tratamiento Variable1 Variable2
1            1       2.2       0.2
2            1       2.3       0.3
3            1       2.4       0.4
4            1       4.5       0.5
5            1       5.4       0.2
6            1       5.4       0.4
7            1       5.4       0.4
8            1       4.0       1.0
9            1       4.0       0.9
10           1       4.4       0.8

tto4 <- datos[31:40,]
> tto4
Tratamiento Variable1 Variable2
31           4        45       9.0
32           4        54       9.8
33           4        46       9.8
34           4        56       9.8
35           4        76       9.9
36           4        54       9.0
37           4        54       7.0
38           4        34       8.0
39           4        45       9.0
40           4        44       6.0

#PREGUNTA10: Nuestra hipótesis nula es que las medias de tratamiento 1 y tratamiento 4 para la
#Variable 1 son iguales. ¿Puedes comprobarlo? Para ello, necesitarás comprobar
#primero si los datos se distribuyen de forma normal. En función del resultado de la
#prueba de normalidad, ¿qué test usarías? ** En general, asumimos que las muestras
#son independientes, pero ¿son sus varianzas iguales? Actúa de acuerdo a tus
#resultados.

#Comprobar distribución normal de Tratamiento
dnormaltto <- shapiro.test(datos$Tratamiento)
pvaluetto <-dnormaltto$p.value
pvaluetto
[1] 0.000204409
#Datos de Tratamiento no se distribuyen de forma normal porque pvalue<0.05

#Distribución normal Variable 1
dnormalVariable1 <- shapiro.test(datos$Variable1)
 pvalueVariable1 <- dnormalVariable1$p.value
 pvalueVariable1
[1] 9.522027e-06
 #Datos de Variable 1 no se distribuyen de forma normal porque pvalue<0.05
 
 #Distribución normal Variable 2
 dnormalVariable2 <- shapiro.test(datos$Variable2)
 pvalueVariable2 <- dnormalVariable2$p.value
 pvalueVariable2
 [1] 1.875741e-05
#Datos de Variable 2 no se distribuyen de forma normal porque pvalue<0.05

#Para comparar 2 grupos independientes que no siguen una distribución normal empleamos 
#el Mann-Withney-Wilcoxon test. 
 
 #Tenemos que extraer las medias de tratamiento y tratamiento 4 para la variable 1, 
 # y después realizar el test sobre estas medias.
 
 #Media tto1 para variable 1
 mediatto1 <- aggregate(Variable1~Tratamiento,tto1,mean)
  mediatto1
 Tratamiento Variable1
 1           1         4
 
 #Media tto4 para variable 1
 mediatto4 <- aggregate(Variable1~Tratamiento,tto4,mean)
 mediatto4
 Tratamiento Variable1
 1           4      50.8
 
 #Mann-Whitney-Wilcoxon test
 MMTO1vsMMTO4 <- wilcox.test(mediatto4$Variable1,mediatto1$Variable1)
  MMTO1vsMMTO4
 
 Wilcoxon rank sum exact test
 
 data:  mediatto4$Variable1 and mediatto1$Variable1
 W = 1, p-value = 1
 alternative hypothesis: true location shift is not equal to 0

 #El p-value=1 nos indica que las medias, de acuerdo con este resultado, son iguales. 
 #Es decir, confirma nuestra H0.
 
 #Al ser muestras independientes, hay que calcular la varianza para el 
 #tto 1 y el tto 4, ya que sabemos que no seguimos una distribución normal
 
 
 #Varianza Tratamiento 1
 Vartto1 <- aggregate(Variable1~Tratamiento,tto1,var)
 Vartto1
 Tratamiento Variable1
 1           1  1.664444
 
 #Varianza Tratamiento 4
 Vartto4 <- aggregate(Variable1~Tratamiento,tto4,var)
 Vartto4
 Tratamiento Variable1
 1           4  123.5111
 
 #Comparamos ambas variables por Mann-Whitney-Wilcoxon
 VTT01vsVTT04 <- wilcox.test(Vartto1$Variable1,Vartto4$Variable1)
 VTT01vsVTT04
 
 Wilcoxon rank sum exact test
 
 data:  Vartto1$Variable1 and Vartto4$Variable1
 W = 0, p-value = 1
 alternative hypothesis: true location shift is not equal to 0

 #El p-value=1 nos indica que estadísticamente las varianzas no son diferentes.
 #Es decir, que los resultados son iguales.