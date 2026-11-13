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
| `jatos_foton_v2.cc` | `HardQCD:all` | o **ruído** — fótons de π⁰, de outro hádron e de fragmentação |
| `foton_direto.cc` | `PromptPhoton:all` | o **sinal** — fótons que saem do processo duro |

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

Ruído: **4 milhões** de eventos `HardQCD`, 1 milhão por faixa de pThat.
Sinal: 400 mil de `PromptPhoton`, 100 mil por faixa. Tudo pesado por σ/N.

**Fração isolada, por origem** — a ordenação é o resultado central:

| origem | fração isolada |
|---|---|
| **fóton direto** (sinal) | **89,2 %** |
| fragmentação | 67,9 % |
| outro hádron | 12,4 % |
| π⁰ | 10,5 % |

O isolamento guarda quase nove em cada dez fótons diretos e descarta nove em
cada dez de decaimento. É exatamente o que deveria fazer: o fóton direto
nasce sozinho no processo duro, o de decaimento vem dentro de um jato.

(Note que aqui a fragmentação é o que está sendo enriquecido, porque esta
amostra é só `HardQCD:all` e não contém fóton direto nenhum. O sinal de
verdade está na seção seguinte.)

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
| fragmentação | 67,9 % |
| outro hádron | 12,4 % |
| π⁰ | 10,5 % |

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

Um fóton energético isolado pode ser, na verdade, filho de um π⁰ cujo irmão
ficou mole demais para ser notado. Para pegar esses, cada fóton energético é
combinado com **todo outro fóton** do evento que um calorímetro veria
(|η| < 0,9 e pT > 0,5 GeV/c) e a massa do par é calculada. Se algum par cair
na massa do π⁰, o fóton é marcado.

O espectro reproduz o PDG na terceira casa:

| pico medido | partícula | PDG |
|---|---|---|
| 0,1350 | π⁰ | 0,13498 |
| 0,5470 | η | 0,54786 |

**O resultado da marcação** (somado nas 4 faixas, pesado por σ):

| origem | marcado | **sem parceiro visível** |
|---|---|---|
| π⁰ | 33,8 % | **62,3 %** |
| outro hádron | 0,7 % | 95,1 % |
| fragmentação | **0,0 %** | 98,5 % |

A marcação tira **22,1 %** da amostra isolada e **não custa nada** em
fragmentação — nenhum fóton de fragmentação foi marcado por engano. A fração
de fragmentação sobe de 20,4 % para 26,2 %.

Mas o número que importa é o outro: **62,3 % dos fótons de π⁰ não têm
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

### Suposição declarada: a marcação não remove fóton direto

A marcação é aplicada **só ao ruído**. Na curva de pureza o sinal é o mesmo
histograma nas duas etapas.

O motivo é que a massa invariante existe para retirar fóton de **decaimento**
de π⁰, e na amostra de sinal a verdade da simulação já garante que todo fóton
é direto — não há decaimento para retirar. Aplicá-la ali mediria outra coisa:
quantos fótons diretos cairiam na janela do π⁰ **por acaso**, formando par com
um fóton qualquer do evento.

Esse efeito deve ser desprezível. Para um par dar a massa do π⁰ os dois fótons
precisam estar muito próximos em ângulo, e o fóton direto isolado, por
definição, tem pouca coisa em volta.

**Como conferir o tamanho da suposição:** o fóton de **fragmentação** do ruído
é o análogo mais próximo do direto — também não vem de decaimento e também
sobrevive ao isolamento. A coluna `marcado %` da tabela de marcação mostra
quantos deles são marcados por engano.

**O número, com 4 milhões de eventos:** nenhum dos **386** fótons de
fragmentação isolados foi marcado. Pela regra de três de Poisson, zero em
386 dá um limite de

> marcação acidental em fótons tipo-sinal **< 0,8 %** a 95 % de confiança

A suposição está sustentada pelos dados, não só pelo argumento. Com 400 mil
eventos o limite era de ~7 %, frouxo demais para citar.

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
| sinal (γ direto) | 6,437 × 10⁻⁴ mb | 6,437 × 10⁻⁴ mb |
| ruído total | 1,732 × 10⁻² mb | 1,349 × 10⁻² mb |
| — π⁰ | 1,128 × 10⁻² mb | 7,463 × 10⁻³ mb |
| — outro hádron | 2,516 × 10⁻³ mb | 2,497 × 10⁻³ mb |
| — fragmentação | 3,526 × 10⁻³ mb | 3,526 × 10⁻³ mb |
| **pureza** | **3,58 %** | **4,56 %** |

A marcação por massa tira 22,1 % do ruído e dá **+0,97 ponto** de pureza. O
sinal é o mesmo nas duas colunas, pela suposição declarada acima.

Na convenção *prompt* (fragmentação conta como sinal), só com isolamento:
**23,22 %**.

Diferencial em pT — a pureza **sobe** com o pT, porque o espectro do sinal
cai mais devagar que o do fundo:

| pT [GeV/c] | só isolamento | + massa | N_B |
|---|---|---|---|
| 5 – 10 | 2,9 % ± 0,2 % | 3,8 % ± 0,3 % | 255 |
| 10 – 15 | 18,7 % ± 2,3 % | 20,8 % ± 2,6 % | 44 |
| 15 – 20 | 28,6 % ± 4,6 % | 30,2 % ± 5,0 % | 20 |
| 20 – 25 | 29,5 % ± 8,3 % | 37,3 % ± 8,9 % | 6 |
| 25 – 30 | 28,5 % ± 9,1 % | 28,5 % ± 9,1 % | 5 |
| 30 – 35 | 32,5 % ± 15,5 % | 32,5 % ± 15,5 % | 2 |
| 35 – 40 | 37,7 % ± 23,5 % | 37,7 % ± 23,5 % | 1 |

**Estatística.** O fundo tem **1212 eventos** de Monte Carlo no espectro
inteiro (eram 116 antes do 1 milhão por faixa), contra 27 485 do sinal. Todos
os 7 bins agora têm fundo — antes três apareciam vazios. `N_B` é o número
**efetivo** de eventos por bin, (conteúdo/erro)², e mostra que acima de
25 GeV/c a curva ainda se apoia em poucas unidades. Ali o erro citado é
grande de propósito; não é para ser lido como medida precisa.

Bins sem nenhum evento de fundo apareceriam **vazios**, não como pureza de
100 %: zero eventos gerados não quer dizer seção de choque zero.

---

## Como abrir os resultados sem rodar nada

Todos os `.root` estão versionados. Entre na pasta e escolha o que quer ver:

| quero ver | comando |
|---|---|
| a soma do ruído, 13 canvases | `root -l soma_pthat.cc` |
| a soma com eixo em ocorrências | `root -l 'soma_pthat.cc(100000)'` |
| as 4 faixas + a soma | `root -l ve_faixas.cc` |
| as 4 faixas sem a soma | `root -l 've_faixas.cc(1)'` |
| as 4 faixas em ocorrências | `root -l 've_faixas.cc(0, 100000)'` |
| o sinal (fótons diretos) | `root -l 'abre_jatos_foton_v2.cc("foton_direto_pthat05.root")'` |
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
| 4 – 7 | 114,56 | 80,9 % | 9,4 % |
| 7 – 10 | 18,38 | 13,0 % | 19,4 % |
| 10 – 15 | 6,44 | 4,5 % | 23,6 % |
| > 15 | 2,20 | 1,6 % | **47,6 %** |

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

Conferência: as 4 σ exclusivas somam 141,577 mb contra 140,818 mb da rodada
inclusiva de pThat > 4, e a soma em ocorrências reproduz uma rodada inclusiva de
100 mil eventos em 0,1 % a 0,8 %.

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
| `soma_pthat.cc` | soma as 4 faixas com teto, pesadas por σ/N (68 histogramas) |
| `pureza.cc` | S/(S+B) em pT, sinal e ruído nas mesmas faixas |
| `ve_faixas.cc` | as 4 faixas sobrepostas, com a soma em magenta |
| `abre_todos.cc` | as 4 sobrepostas, 42 histogramas, sem somar |
| `compara_pthat.cc` | teste de consistência das faixas encaixadas |
| `abre_jatos_foton_v2.cc` | abre uma rodada qualquer, 12 canvases |
| `abre_jatos_foton.cc` | o mesmo, para o formato do v1 |

### Resultados

| arquivo | conteúdo |
|---|---|
| `jatos_foton_pthat04_07/07_10/10_15/15_inf.root` | as 4 faixas com teto, 100 mil cada |
| `jatos_foton_pthat04/07/10/15.root` | as 4 sem teto, 100 mil cada |
| `soma_pthat.root` | a soma, 68 histogramas |
| `foton_direto_pthat04_07/07_10/10_15/15_inf.root` | o sinal nas 4 faixas |
| `pureza.root` | espectros de S e B e a curva de pureza |
| `foton_direto_pthat05.root` | o sinal, 100 mil eventos de pThat > 5 |
| `jatos_foton_isolado.root` | 1 milhão de eventos, pThat > 4 |
| `jatos_foton.root` | 1 milhão, código v1 |
| `figuras/` | os PNGs |

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

```bash
# o ruido
./producao.sh              # os dois conjuntos, 100 mil cada (~45 min)
./producao.sh 5000         # teste rápido
./producao.sh 100000 teto  # só as faixas com teto (~23 min)
root -l -b -q soma_pthat.cc

# o sinal
root -l -b -q 'foton_direto.cc(5., -1., 43, 5000, "teste")'   # teste, ~50 s
root -l -b -q 'foton_direto.cc(5., -1., 43, 100000, "")'      # ~16 min
./producao_direto.sh 100000   # o sinal nas 4 faixas (~13 min)

# a pureza, depois dos dois
root -l -b -q pureza.cc
```

Precisa de ROOT com `TPythia8` e do FastJet no `ROOT_INCLUDE_PATH`.
