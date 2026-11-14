# Fótons isolados em jatos — produção fatiada em pThat

Mede quanto a **técnica de isolamento** consegue limpar uma amostra de fótons
energéticos, e separa o **sinal** (fótons diretos, do processo duro) do que
o imita: fragmentação e decaimento de hádrons.

Pythia 8 + FastJet dentro do ROOT (`TPythia8`, `TClonesArray`, `TParticle`).
pp a 14 TeV, anti-kT com R = 0,4.

São **duas amostras separadas**, geradas por dois códigos irmãos com os
mesmos cortes:

| | processo | o que mede |
|---|---|---|
| `jatos_foton_v2.cc` | `HardQCD:all` | o **ruído** — fótons de π⁰ + η, de outro hádron e de fragmentação |
| `foton_direto.cc` | `PromptPhoton:all` | o **sinal** — fótons que saem do processo duro |

**As três origens do ruído:**

| categoria | mãe do fóton |
|---|---|
| **π⁰ + η** | PDG 111 ou 221 |
| outro hádron | qualquer outro, PDG > 100 — **já sem o η** |
| fragmentação | quark, PDG 1 a 6 |

O η está junto com o π⁰ porque faz a mesma coisa: decai em dois fótons
(η → γγ, ~39 %) e, quando um dos dois fica mole, o outro se disfarça de fóton
direto. É o mesmo mecanismo de contaminação, e a marcação por massa invariante
pega os dois — basta olhar também a janela de massa do η.

Os dois usam a mesma `soma_cone`, o mesmo corte de isolamento (ΣpT < 2 GeV/c
em R < 0,4), a mesma aceitação (|η| < 0,5) e o **mesmo binning** no espectro
de pT (35 bins de 5 a 40 GeV/c), de propósito: os histogramas podem ser
sobrepostos bin a bin.

### Cuidado com a palavra "sinal"

Há **duas convenções** na literatura, e elas dão números diferentes:

| convenção | sinal | fundo |
|---|---|---|
| *prompt* | direto **+ fragmentação** | só decaimento de hádron |
| *direto* | só direto | fragmentação + decaimento |

**Este repositório usa a segunda**: sinal = fóton direto, e a fragmentação
conta como fundo. É a convenção pedida pelo orientador.

A primeira é a mais comum em medida experimental — ATLAS, CMS e ALICE medem
*isolated prompt photon* — porque a separação entre direto e fragmentação
**não é observável**: depende do esquema e da escala de fatorização em NLO.
O fóton de fragmentação também nasce do espalhamento duro, não de decaimento.

Ao citar número daqui num texto, diga qual convenção está usando.

---

## Resultado

Ruído: **13 milhões** de eventos `HardQCD` — 1 milhão em cada uma das faixas
4–7, 7–10 e 10–15, e **10 milhões** na faixa > 15, que é a que manda na cauda
de pT alto. Sinal: 400 mil de `PromptPhoton`, 100 mil por faixa. Tudo pesado
por σ/N, então as estatísticas desiguais se compensam sozinhas.

**Fração isolada, por origem** — a ordenação é o resultado central:

| origem | fração isolada | ⟨E_iso⟩ |
|---|---|---|
| **fóton direto** (sinal) | **89,2 %** | 0,74 GeV/c |
| fragmentação | 73,1 % | 1,70 GeV/c |
| outro hádron | 11,0 % | 7,32 GeV/c |
| π⁰ + η | 11,0 % | 7,60 GeV/c |

O isolamento guarda quase nove em cada dez fótons diretos e descarta nove em
cada dez de decaimento. É exatamente o que deveria fazer: o fóton direto
nasce sozinho no processo duro, o de decaimento vem dentro de um jato.

(Note que aqui a fragmentação é o que está sendo enriquecido, porque esta
amostra é só `HardQCD:all` e não contém fóton direto nenhum. O sinal de
verdade está na seção seguinte.)

**De onde vem esse número.** A fração isolada é medida do histograma de
`E_iso`, que conta **todo** fóton energético. Até esta versão ela também era
impressa a partir dos histogramas de `z`, e ali saía errada para a
fragmentação — **67 % em vez de 73 %**. A causa: o histograma de `z` vai de 0
a 1 e o `Integral()` do ROOT, chamado sem argumentos, descarta o overflow. Um
fóton de fragmentação isolado é praticamente o jato inteiro, ou seja z ≈ 1,
e caía justamente lá. Eram 26 % dos fótons de fragmentação isolados fora da
conta. Corrigido com `Integral(0, nbins+1)`; as duas tabelas agora fecham
dígito por dígito, o que vale como teste de consistência.

O viés era só da fragmentação porque só ela tem z ≈ 1 — o fóton de decaimento
divide o jato com os outros hádrons e fica longe da borda.

---

## O sinal: fótons diretos

`foton_direto.cc`, `PromptPhoton:all`, pThat > 5 GeV, 100 mil eventos:

| | |
|---|---|
| σ | 0,00747141 ± 0,0000115 mb |
| fótons energéticos | 7541 |
| **diretos** | 7278 (96,5 %) |
| **diretos e isolados** | 6349 (**87,2 %** dos diretos) |

Somando as 4 faixas de pThat com peso σ/N, a fração isolada do sinal é
**89,2 %**. O número acima, 87,2 %, é da rodada única de pThat > 5.
| σ no gráfico, 5–40 GeV/c | 4,718 × 10⁻⁴ mb |

A fração isolada põe as quatro origens em ordem:

| origem | fração isolada |
|---|---|
| **direto** (sinal) | **89,2 %** |
| fragmentação | 73,1 % |
| outro hádron | 11,0 % |
| π⁰ + η | 11,0 % |

É exatamente o que o isolamento deveria fazer: o fóton direto nasce sozinho no
processo duro, então quase sempre passa; o de decaimento vem dentro de um jato
de hádrons e quase nunca passa.

Os 3,5 % de fótons energéticos que **não** são diretos nessa amostra são
decaimento no jato de recuo — o teste não aprova tudo que vê.

### Como o fóton direto é reconhecido

No Pythia 8 as duas partículas que **saem** do processo duro estão nas linhas
5 e 6. O `TPythia8::ImportParticles` descarta a linha 0 (a pseudo-partícula
`system`, PDG 90) e desloca todos os índices em −1, mães inclusive — então no
`TClonesArray` elas ficam nos índices **4 e 5**. Isso está conferido pelo
diagnóstico que a macro imprime nos 3 primeiros eventos, não suposto.

O fóton final quase nunca **é** essa linha: o chuveiro o recopia. Nos eventos
de teste as cadeias tinham 3, 5 e 6 cópias:

```
[320] <- [264] <- [72] <- [47] <- [17] <- [8] <- [4]
```

Por isso o teste sobe pelas mães **enquanto a mãe for PDG 22** e pergunta se o
fóton mais antigo da cadeia é o índice 4 ou o 5.

**A ordem importa.** Logo acima do topo vem um quark ou um glúon — o párton
que entrou na colisão. A `origem_foton` do `jatos_foton_v2.cc` chamaria isso de
*fragmentação*. Por isso o teste do processo duro vem **antes** de qualquer
teste de "a mãe é quark".

### Comparar sinal e ruído

Na janela de 5 a 40 GeV/c, o sinal vale 4,72 × 10⁻⁴ mb e o ruído isolado
somado (π⁰ + fragmentação + outro hádron) vale 1,28 × 10⁻² mb.

**Cuidado ao usar essa razão:** o ruído foi gerado com pThat > 4 e o sinal com
pThat > 5, então os dois não estão na mesma janela de pThat. Para uma razão
honesta, rode o `foton_direto.cc` nas mesmas 4 faixas do ruído — a macro já
aceita `pt_hat_max` para isso.

---

## Massa invariante: quanto do fundo dá para marcar

Um fóton energético isolado pode ser, na verdade, filho de um π⁰ (ou de um η)
cujo irmão ficou mole demais para ser notado. Para pegar esses, cada fóton
energético é combinado com **todo outro fóton** do evento que um calorímetro
veria (|η| < 0,9 e pT > 0,5 GeV/c) e a massa do par é calculada.

**Duas janelas**, e o fóton é marcado se o par cair em qualquer uma delas:

| janela | massa | meia-largura |
|---|---|---|
| π⁰ | 0,1349768 GeV/c² | ± 20 MeV |
| η | 0,547862 GeV/c² | ± 20 MeV |

As duas não se encostam (0,115–0,155 e 0,528–0,568), então um mesmo par nunca
cai nas duas. Mas um mesmo **fóton** pode ter um parceiro em cada janela — por
isso o total de marcados não é a soma das duas colunas, e a contagem que vale
é a dos espectros `hist_pt_pos_*`, que guardam quem de fato sobreviveu.

O código imprime as duas janelas separadas, para dar para ver qual delas está
marcando o quê — em especial quantos fótons de **fragmentação** cada uma marca
por engano, que é o preço da marcação.

O espectro reproduz o PDG na terceira casa:

| pico medido | partícula | PDG |
|---|---|---|
| 0,1350 | π⁰ | 0,13498 |
| 0,5470 | η | 0,54786 |

**O resultado da marcação** (somado nas 4 faixas, pesado por σ):

| origem | marcado | **sem parceiro visível** |
|---|---|---|
| π⁰ + η | 35,9 % | **70,9 %** |
| outro hádron | 1,9 % | 99,9 % |
| fragmentação | **0,3 %** | 99,9 % |

A marcação tira **28,0 %** da amostra isolada e custa **0,3 %** em
fragmentação. A fração de fragmentação no que sobra sobe de 20,0 % para
27,8 % — esperado, já que a marcação remove decaimento e quase não toca
fragmentação, mas vale notar: o ruído que resta é proporcionalmente mais
"tipo-sinal" do que o de antes.

Separando as duas janelas, quem cobra o preço é o **η**: dos fótons de
fragmentação marcados por engano, a janela do η responde por 8,9 × 10⁻⁶ mb
contra 1,1 × 10⁻⁶ mb da do π⁰ — oito vezes mais. O motivo é geométrico: a
janela do η fica em massa mais alta, e um par aleatório de fótons chega mais
fácil a 0,548 GeV/c² do que a 0,135.

Mas o número que importa é o outro: **70,9 % dos fótons de π⁰ e η não têm
parceiro visível nenhum**. Esses não há técnica de massa que resolva — é o
**fundo irredutível**, e existe exatamente pelo motivo que motiva o trabalho:
o decaimento do π⁰ tem assimetria plana, então uma fração constante dos
decaimentos joga quase toda a energia num fóton só.

O limiar de 0,5 GeV/c no parceiro é o parâmetro que manda nesse número. Ele
imita o limiar de cluster de um calorímetro real. Baixá-lo faria a marcação
parecer melhor do que é.

### Os marcados são removidos

Cada rodada grava **dois** espectros pareados:

| histograma | conteúdo |
|---|---|
| `hist_pt_iso_<origem>` | isolado — antes da marcação |
| `hist_pt_pos_<origem>` | isolado **e não marcado** — depois |

A ordem é isolamento primeiro, massa invariante depois.

### A marcação agora é aplicada nas duas amostras

Numa versão anterior a marcação era aplicada **só ao ruído**, e o sinal
entrava idêntico nas duas etapas da pureza. Aquilo era uma **suposição**: que a
marcação não remove fóton direto nenhum.

Agora o `foton_direto.cc` também marca, com os mesmos critérios, e guarda o
`hist_pt_direto_pos`. A perda de sinal deixou de ser suposta e passou a ser
**medida**. É importante que seja: a pureza depois da marcação é

```
P = S_depois / (S_depois + B_depois)
```

e um corte que tira 20 % do fundo mas também 20 % do sinal não melhora nada.
Com as duas amostras marcadas, a conta diz a verdade sozinha — o `pureza.cc`
imprime a razão `ruído tirado / sinal perdido`, que só vale a pena acima de 1.

O que a marcação pega no sinal é sempre **acaso**: ali não existe fóton de
decaimento, então todo fóton direto marcado formou massa de π⁰ ou de η com um
fóton qualquer do evento, sem ter nada a ver com ele.

**E o acaso é pequeno.** Medido direto na amostra de sinal, pesado por σ:

> **0,39 %** dos fótons diretos isolados são marcados por acaso

Em contagem bruta são 98 fótons em 27 485, e a separação por janela repete o
padrão do ruído: 73 pela janela do η contra 25 pela do π⁰, quase 3× mais.

Isso fecha com folga a estimativa anterior, que era um limite de Poisson
inferido por analogia com os fótons de fragmentação (zero marcados em 386,
logo < 0,8 % a 95 % de confiança). O limite valia; agora há a medida.

A conta que decide se a marcação vale a pena é a razão entre as duas:

> ruído tirado **28,0 %** ÷ sinal perdido **0,39 %** = **72,6**

Um corte que tirasse 28 % do fundo e 28 % do sinal não melhoraria nada. Este
tira 72 vezes mais fundo que sinal, e por isso a pureza sobe.

---

## Isolamento por origem: a distribuição de E_iso

`E_iso` é a soma do pT das partículas finais no cone de R < 0,4 em volta do
fóton — a mesma soma que o corte de isolamento compara com 2 GeV/c. Guardada
num histograma por origem (0 a 20 GeV/c, 100 bins de 0,2 GeV), para **todo**
fóton energético, isolado ou não.

A largura do bin foi escolhida para que o corte de 2,0 GeV/c caia numa
**borda**: assim a integral abaixo dele é exatamente a eficiência do
isolamento, sem partir bin nenhum ao meio.

O gráfico `pureza_eiso.png` põe as **quatro** curvas juntas — o sinal e as três
origens do ruído — cada uma normalizada por fóton (área 1), com uma linha
vertical no corte. É o gráfico que explica o isolamento: mostra, para cada tipo
de fóton, como a atividade em volta dele se distribui. O fóton direto nasce
sozinho e empilha perto de zero; o de decaimento vem dentro de um jato e
espalha para a direita. A fração de cada curva à esquerda da linha é a
eficiência do isolamento para aquele tipo — os mesmos números da tabela de
fração isolada, agora com a distribuição inteira visível por trás deles.

A versão só com o ruído, por faixa de pThat, sai em `v2_<rodada>_eiso.png`; a
somada, em `soma_eiso.png`.

---

## Três escalas do mesmo espectro

O mesmo espectro de pT aparece em três eixos verticais diferentes. **São só
mudanças de escala**, não medidas diferentes:

| eixo | o que é | para quê |
|---|---|---|
| `σ [mb] por bin` | o padrão, o que a soma σ/N produz | é o único que serve para a **pureza** |
| `(1/N) dN/dX` | cada curva dividida pela própria área (área 1) | comparar **formas** quando as estatísticas são muito diferentes |
| `(1/2π) d²σ/(dpT dy)` | σ do bin ÷ (2π · largura do bin · dy) | é a forma que se compara com **dado publicado**, porque não depende de quão largo o bin foi escolhido |

Na invariante, `dy = 1,0`: os fótons são selecionados com |η| < 0,5 e, para
fóton (que não tem massa), y = η, então dy = 0,5 − (−0,5) = 1.

Os dois parâmetros estão no começo do `soma_pthat.cc` e do `pureza.cc`, e
**têm que ser iguais nos dois**:

```cpp
const Double_t SP_DY = 1.0;              // intervalo de rapidez
const Bool_t   SP_USA_1_SOBRE_PT = kFALSE;   // liga o fator 1/pT
```

Ligando a chave, o eixo passa a ser `E d³σ/dp³ = (1/2π pT) d²σ/(dpT dy)`, em
mb/(GeV/c)², que é a outra convenção da literatura. Está em `kFALSE` **até
confirmar com o orientador** qual das duas usar.

**A pureza continua sendo calculada em seção de choque, sempre.** Razão de
curvas normalizadas não é pureza: a normalização joga fora justamente a
informação de tamanho que distingue sinal de fundo.

---

## Pureza: quanto do sinal é ruído?

```
pureza(pT) = S(pT) / [ S(pT) + B(pT) ]
```

Com as **duas amostras somadas nas mesmas 4 faixas de pThat**, cada uma pesada
pela sua própria σ. Isso não é detalhe: uma rodada de `pTHatMin = X` é
incompleta perto de X, e comparar faixas diferentes estragaria a pureza
justamente em pT baixo, onde há mais estatística.

Integrado de 5 a 40 GeV/c:

| | só isolamento | + massa invariante |
|---|---|---|
| sinal (γ direto) | 6,4368 × 10⁻⁴ mb | 6,4120 × 10⁻⁴ mb |
| ruído total | 1,7366 × 10⁻² mb | 1,2500 × 10⁻² mb |
| — π⁰ + η | 1,3517 × 10⁻² mb | 8,6674 × 10⁻³ mb |
| — outro hádron | 3,6829 × 10⁻⁴ mb | 3,6141 × 10⁻⁴ mb |
| — fragmentação | 3,4806 × 10⁻³ mb | 3,4711 × 10⁻³ mb |
| **pureza** | **3,57 %** | **4,88 %** |

A marcação por massa tira 28,0 % do ruído e dá **+1,31 ponto** de pureza. O
sinal **também** cai de uma coluna para a outra, só muito menos (0,39 %) —
essa perda é medida, não suposta.

Na convenção *prompt* (fragmentação conta como sinal), só com isolamento:
**22,90 %**.

Diferencial em pT — a pureza **sobe** com o pT, porque o espectro do sinal
cai mais devagar que o do fundo:

| pT [GeV/c] | só isolamento | + massa | N_B |
|---|---|---|---|
| 5 – 10 | 2,9 % ± 0,2 % | 4,0 % ± 0,3 % | 257 |
| 10 – 15 | 18,3 % ± 2,1 % | 20,7 % ± 2,5 % | 50 |
| 15 – 20 | 34,1 % ± 4,5 % | 35,8 % ± 5,0 % | 25 |
| 20 – 25 | 30,0 % ± 6,6 % | 39,4 % ± 3,0 % | 10 |
| 25 – 30 | 48,7 % ± 5,5 % | 48,6 % ± 5,5 % | 21 |
| 30 – 35 | 39,1 % ± 6,2 % | 39,1 % ± 6,2 % | 15 |
| 35 – 40 | 46,4 % ± 9,5 % | 46,4 % ± 9,5 % | 7 |

**Estatística.** O fundo tem **5651 eventos** de Monte Carlo no espectro
inteiro — eram 116 na primeira produção e 1212 com 1 milhão por faixa —
contra 27 485 do sinal. Os 10 milhões na faixa > 15 foram para cá: o bin de
30–35 GeV/c passou de 2 eventos efetivos para 15, e o de 35–40 de 1 para 7.
Os erros em pT alto caíram de ±15,5 % e ±23,5 % para ±6,2 % e ±9,5 %.

`N_B` é o número **efetivo** de eventos por bin, (conteúdo/erro)². Acima de
30 GeV/c a curva ainda se apoia em poucas unidades, mas já não em unidades
isoladas.

**Um aviso sobre o bin de 20–25 GeV/c.** É o único em que o erro **diminui**
depois da marcação (±6,6 % → ±3,0 %), o que parece errado e não é. O `N_B`
efetivo daquele bin sobe de 10,2 para 64,0 quando a marcação entra. A causa:
cerca de 30 % do conteúdo vinha de **um único evento** de Monte Carlo da faixa
10–15, cujo peso σ/N é 6,4 × 10⁻⁶ mb, umas 30 vezes o peso de um evento da
faixa > 15. Esse evento sozinho dominava o conteúdo e a variância; a marcação o
removeu, e a variância caiu muito mais do que o conteúdo.

A aritmética está certa — é como variância de histograma pesado se comporta —
mas o ponto "só isolamento" nesse bin é o mais frágil da curva inteira, porque
depende de um evento só. Ao citar a pureza diferencial num texto, vale dizer
isso, ou juntar 20–30 num bin só. É também o argumento mais direto para a
próxima rodada de estatística ser nas faixas **de baixo**, não na de cima: a
faixa > 15 já está firme, e são os eventos de peso alto que ainda mandam na
incerteza em pT intermediário.

Bins sem nenhum evento de fundo apareceriam **vazios**, não como pureza de
100 %: zero eventos gerados não quer dizer seção de choque zero.

---

## Como abrir os resultados sem rodar nada

Todos os `.root` estão versionados. Entre na pasta e escolha o que quer ver:

| quero ver | comando |
|---|---|
| a soma do ruído, 18 canvases | `root -l soma_pthat.cc` |
| a soma com eixo em ocorrências | `root -l 'soma_pthat.cc(100000)'` |
| as 4 faixas + a soma | `root -l ve_faixas.cc` |
| as 4 faixas sem a soma | `root -l 've_faixas.cc(1)'` |
| as 4 faixas em ocorrências | `root -l 've_faixas.cc(0, 100000)'` |
| o sinal (fótons diretos) | `root -l 'abre_jatos_foton_v2.cc("foton_direto_pthat04_07.root")'` |
| **a pureza** | `root -l pureza.cc` |
| a pureza sem reagrupar | `root -l 'pureza.cc(1)'` |
| as 4 sobrepostas (com teto) | `root -l abre_todos.cc` |
| as 4 sobrepostas (sem teto) | `root -l 'abre_todos.cc(kTRUE, kFALSE)'` |
| o teste de consistência | `root -l compara_pthat.cc` |
| uma faixa por inteiro | `root -l 'abre_jatos_foton_v2.cc("jatos_foton_pthat15_inf.root")'` |

Nenhum deles simula; todos leem os `.root` em segundos.

**Unidade do eixo vertical.** O padrão agora é **seção de choque, `σ [mb] por
bin`** — é o que a soma pesada produz naturalmente e é como o resultado se
compara com dado experimental e com o sinal. O eixo em número de ocorrências
continua existindo, mas é preciso pedir: passe o número de eventos
equivalente. As tabelas impressas saem sempre em mb, e também na unidade do
eixo quando ela for outra.

---

## Por que o pThat é fatiado

Fóton energético é raro em evento mole. Rodando só `pTHatMin = 4`, quase todo
evento é de pThat baixo e a cauda de pT alto fica vazia.

A tabela abaixo é a justificativa inteira do método:

| faixa | σ [mb] | % da seção de choque | % dos fótons energéticos |
|---|---|---|---|
| 4 – 7 | 114,28 | 80,9 % | 14,3 % |
| 7 – 10 | 18,39 | 13,0 % | 19,2 % |
| 10 – 15 | 6,44 | 4,6 % | 22,3 % |
| > 15 | 2,20 | 1,6 % | **44,2 %** |

Quase metade dos fótons energéticos vem de 1,6 % dos eventos. Sobre-amostrar o pThat alto e
corrigir com o peso σ/N mede a mesma física com **30× mais precisão** pelo mesmo
custo de CPU.

### Com teto e sem teto

Há **dois** conjuntos de 4 rodadas, e a diferença importa:

- **Com teto** (`pthat04_07`, `pthat07_10`, `pthat10_15`, `pthat15_inf`) — usam
  `PhaseSpace:pTHatMax`, são faixas **exclusivas** e **podem ser somadas**.
- **Sem teto** (`pthat04`, `pthat07`, `pthat10`, `pthat15`) — só `pTHatMin`, são
  faixas **encaixadas**. Servem para comparar entre si, mas somá-las conta a
  região dura quatro vezes. O `compara_pthat.cc` mede esse fator: dá 4,00.

A última faixa não tem teto nos dois conjuntos — ela precisa ficar aberta por
cima, senão o espectro acabaria num corte artificial. São 3 tetos, não 4.

### A soma

```
h_total = Σ_i (σ_i / N_i) · h_i
```

O peso σ/N é lido do **bin 5** do `hist_sigma`, que cada rodada grava no seu
próprio `.root`. Cada faixa entra com a **sua própria** σ — não há média de
seções de choque em lugar nenhum.

Conferência: as 4 σ exclusivas somam 141,314 mb contra 140,818 mb da rodada
inclusiva de pThat > 4 — 0,35 % de diferença. A faixa > 15, com 10 milhões de
eventos, deu σ = 2,19964 ± 0,00039 mb contra 2,20348 ± 0,00389 mb da rodada de
1 milhão: mesmo valor dentro do erro, com a incerteza 10× menor, que é o √10
esperado. A soma em ocorrências reproduz uma rodada inclusiva de 100 mil
eventos em 0,1 % a 0,8 %.

---

## Arquivos

### Simulação

| arquivo | o que faz |
|---|---|
| `jatos_foton_v2.cc` | o RUÍDO — `HardQCD:all`, classifica origem, isola |
| `foton_direto.cc` | o SINAL — `PromptPhoton:all`, fótons do processo duro |
| `producao_direto.sh` | roda o sinal nas 4 faixas de pThat |
| `producao.sh` | roda as 4 faixas em paralelo (`./producao.sh` para os dois conjuntos) |
| `jatos_foton.cc` | versão v1, parte do jato (mantida para referência) |

### Análise e visualização

| arquivo | o que faz |
|---|---|
| `soma_pthat.cc` | soma as 4 faixas com teto, pesadas por σ/N (conta e imprime quantos histogramas somou) |
| `pureza.cc` | S/(S+B) em pT, sinal e ruído nas mesmas faixas |
| `ve_faixas.cc` | as 4 faixas sobrepostas, com a soma em magenta |
| `abre_todos.cc` | as 4 sobrepostas, 42 histogramas, sem somar |
| `compara_pthat.cc` | teste de consistência das faixas encaixadas |
| `abre_jatos_foton_v2.cc` | abre uma rodada qualquer, 12 canvases |
| `abre_jatos_foton.cc` | o mesmo, para o formato do v1 |

### Resultados

| arquivo | conteúdo |
|---|---|
| `jatos_foton_pthat04_07/07_10/10_15.root` | as 3 faixas de baixo, 1 milhão cada |
| `jatos_foton_pthat15_inf.root` | a faixa de cima, **10 milhões** |
| `jatos_foton_pthat04/07/10/15.root` | as 4 sem teto, 100 mil cada |
| `soma_pthat.root` | a soma das 4 faixas com teto |
| `foton_direto_pthat04_07/07_10/10_15/15_inf.root` | o sinal nas 4 faixas |
| `pureza.root` | espectros de S e B e a curva de pureza |
| `foton_direto_pthat05.root` | o sinal, 100 mil eventos de pThat > 5 |
| `jatos_foton_isolado.root` | 1 milhão de eventos, pThat > 4 |
| `jatos_foton.root` | 1 milhão, código v1 |
| `figuras/` | os PNGs |

### As figuras

| arquivo | eixo vertical |
|---|---|
| `soma_gerais.png` | σ [mb] por bin |
| `soma_gerais_norm.png` | **(1/N) dN/dX**, área 1 |
| `soma_pt_iso_por_origem.png` | σ [mb] por bin |
| `soma_pt_iso_norm.png` | **(1/N) dN/dpT**, área 1 |
| `soma_pt_iso_invariante.png` | **(1/2π) d²σ/(dpT dy)** |
| `soma_eiso.png` | E_iso das 3 origens do ruído, área 1 |
| `soma_massa_invariante.png` | 4 painéis: espectro, zoom π⁰, melhor par na janela do π⁰ e na do η |
| `pureza_espectros.png` | σ [mb] por bin, sinal e ruído |
| `pureza_espectros_norm.png` | área 1, sinal e ruído |
| `pureza_invariante.png` | invariante, sinal e ruído, antes e depois da marcação |
| `pureza_eiso.png` | **as 4 curvas de E_iso**, área 1, com o corte |
| `pureza_curva.png` | a pureza em função do pT |
| `foton_direto_<rodada>_massa.png` | massa invariante na amostra de sinal |
| `foton_direto_<rodada>_eiso.png` | E_iso do fóton direto |

---

## Dois detalhes do código que custaram caro

**Fótons recopiados no chuveiro.** Quando um quark irradia um fóton, o Pythia
registra o quark como mãe. Mas o chuveiro continua, e sempre que esse fóton
serve de recuo numa emissão seguinte ele é **recopiado**, com o fóton antigo
(PDG 22) como mãe. Olhando só a mãe direta, todo fóton de fragmentação caía em
"não classificado" e a categoria ficava **vazia por construção**. A
`origem_foton` sobe pelas cópias antes de classificar.

**`TFile::Close()` apaga os histogramas.** Histogramas criados depois de abrir um
`TFile` pertencem a ele, e fechá-lo os destrói — os canvases ficam em branco.
Por isso as macros chamam `Write()` mas **não** `Close()`.

**Fóton direto parece fragmentação.** A mãe de um fóton direto também é um
quark ou um glúon. Qualquer teste de origem que pergunte "a mãe é quark?"
antes de perguntar "isto veio do processo duro?" classifica o sinal inteiro
como fragmentação. No `foton_direto.cc` o teste do processo duro vem primeiro.

---

## Reproduzir do zero

Na ordem:

```bash
# 1) teste rápido, para ver se nada quebrou  (~3 min)
./producao.sh 5000 teto 5000
./producao_direto.sh 5000

# 2) produção de verdade
#    ruído: 1 milhão nas faixas 4-7, 7-10 e 10-15, e 10 MILHÕES em > 15
./producao.sh                      # ~20 h, dominado pela faixa > 15
./producao_direto.sh 100000        # o sinal nas 4 faixas (~13 min)

# 3) soma e pureza  (segundos, só leem os .root)
root -l -b -q soma_pthat.cc
root -l -b -q pureza.cc
```

**Por que 10 milhões só na faixa de cima.** É ela que alimenta a região de pT
alto do espectro, justamente onde a pureza mais sobe e onde a estatística era
de 5, 2 e 1 evento efetivo por bin. As três faixas de baixo já têm estatística
de sobra com 1 milhão.

Isso **não** desequilibra a soma: o peso de cada faixa é σ/N, com o N aceito
lido do `hist_sigma` do próprio arquivo — tanto o `soma_pthat.cc` como o
`pureza.cc` fazem isso, no bin 5. Gerar dez vezes mais eventos numa faixa
divide o peso dela por dez; a seção de choque somada continua a mesma, e só a
barra de erro diminui.

**Tempo.** Na produção anterior, 1 milhão na faixa > 15 levou ~2 h 10 min
dividindo a máquina com as outras três. Dez milhões, com as outras três
terminando em 1–2 h e deixando a máquina livre, devem levar **~18 a 22 horas**.
Rode com `nohup`:

```bash
nohup ./producao.sh > nohup_producao.log 2>&1 &
```

O sinal **precisa** ser rodado de novo: os histogramas de `E_iso`
(`hist_eiso_direto`), de massa invariante e o espectro pós-marcação
(`hist_pt_direto_pos`) são novos, e o `pureza.cc` para com erro sem eles.
São 13 min, não é problema.

Precisa de ROOT com `TPythia8` e do FastJet no `ROOT_INCLUDE_PATH`.
