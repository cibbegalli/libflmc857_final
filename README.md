# README #

*************************************************
<h1>Portuguese</h1>

Para compilar a biblioteca acesse, pelo o terminal, o diretorio onde esta a pasta "libfl" e acesse a pasta fazendo
cd libfl

Uma vez dentro da pasta libfl, escreva no terminal o comando "make" (sem aspas) para compilar a biblioteca.

Algumas breves descricões sobre os diretorios da libFL

*pasta include: contém as estruturas e os cabeçalhos das funções da biblioteca.

*pasta src: contém as implementações das funções da pasta include

*pasta obj: quando os arquivos são compilados são gerado os binários que possuem extensões .o ("o" de object).
Então essa pasta contém os binários de todos os arquivos da biblioteca.

*pasta lib: contém o binário da biblioteca (um arquivo com extensão .a).

*pasta pythonUtils: As vezes o professor Falcão que plotar um gráfico, então essa pasta contem alguns
scripts em python que *para fazer a plotagem de gráficos.

*pasta cmake: a IDE que eu uso para programar em C é o CLION (https://www.jetbrains.com/clion/download/), que é uma
IDE paga, porém estudantes ganha uma licença gratuita e podem usar a IDE de graça (basta se cadastrar no site deles).
Para eu poder usar as ferramentas de debug foi necessário eu fazer um cmake.

*pasta demo: a ideia da pasta demo é você criar códigos (demos) que importem a biblioteca libFL para testar suas
funções que estão na biblioteca. Para compilar uma demo basta você digitar o comando "make" (sem aspas) seguido
do nome da demo que você que compilar sem a extensão. Por exemplo na pasta demo tem uma demo chamada rwImage.cpp,
se você quer compilar ela você acessa o diretório (via terminal) das demos e da o comando.

make rwImage

então ele ira compilar a demo.

****************************************************************
<h1>English</h1>
<br>

****************************************************************

<h2>Compiling the lib</h2>
<br>

Open terminal then, type

"cd PATH/TO/LIBFL/FOLDER"

after that, type

"make"

to compile the lib






****************************************************************
<h2>Compiling a demo</h2>
<br>
Open terminal then, type

"cd PATH/TO/DEMO/FOLDER"

after that, type

"make demoFileNameWithNoExtension"

to compile a demo

EXAMPLE:
"make rwImage"




****************************************************************
<h2>LIB FOLDERS</h2>
<br>

 <ul>
  <li>include: contains the functions prototypes files (.h)</li>
  <li>src: contains the functions implementation (.cpp) (most of them are written in C)</li>
  <li>obj: contains the object files (.o)</li>
  <li>lib: contains the static library file (.a)</li>
  <li>pythonUtils: simples python scripts (.py) for general purposes (such as plotting)</li>
  <li>demo: contains the demo files (.cpp)</li>
</ul>







# libflmc857_final
