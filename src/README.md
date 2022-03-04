## Modelo Hash Table

### Parametros

* **Tamanho da tabela**: Calculado baseado na quantidade esperada de k-mers (`min(4^k, sequence_length)`) e fator de carga `alpha`.
* **Função de hash**: Usamos o algoritmo *Fibonacci hashing*.
* **Função de Fingerprint**: Ao invés de guardar o k-mer, guardamos sua *fingerprint*. Usamos os 3 bits mais significativos da multiplicação do *Fibonacci hashing*.

### Estrutura

Cada espaço na tabela é um byte organizado da seguinte forma:

* **Bit 0**: *Dirty bit*. Indica se o espaço está ocupado ou não. (1 = ocupado, 0 = livre)
* **Bits 1-3**: *Fingerprint*.
* **Bits 4-7**: *out edges* encontradas.
  * **Bit 4**: Indica a presença de uma aresta de saída 'A'.
  * **Bit 5**: Indica a presença de uma aresta de saída 'C'.
  * **Bit 6**: Indica a presença de uma aresta de saída 'G'.
  * **Bit 7**: Indica a presença de uma aresta de saída 'T'.

### Operações

#### Inserção

Caso a posição dada pelo resultado da função de *hash* esteja vazia (*dirty bit* = 0), então a *fingerprint* é guardada e o *dirty bit* é marcado como 1. Do contrário, é realizado *linear probing* até que um espaço livre seja encontrado. A *fingerprint* é usada para identificar se o valor já foi guardado, de forma que caso o valor na posição dada pelo *hash*, ou algum outro valor encontrado durante o *linear probing* tenha a *fingerprint* do valor sendo inserido, ele é considerado como já tendo sido inserido antes e nada é feito.

#### Atualização de Arestas

Quando um k-mer é inserido em sequência a um anterior, é possível atualizar as arestas de saída do k-mer anterior. Para isso, sua entrada na *hash table* é encontrada e o bit correspondente à aresta relevante é setado como 1.

#### Consulta

A consulta é feita como na inserção ou na atualização, primeiro com a posição dada pela função de *hash*, sendo realizado *linear probing* até que a *fingerprint* do valor seja encontrada. É retornado, então um byte com a mascara de bits correspondente às arestas de saída. Caso o valor não seja encontrado, um byte com todos os bits setados (ou seja, o valor 255) é retornado.