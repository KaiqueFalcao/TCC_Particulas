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

Soma das 4 faixas de pThat, 400 mil eventos, em ocorrências equivalentes a
100 mil eventos de pThat > 4 GeV:

| origem | isolados | total | fração isolada |
|---|---|---|---|
| π⁰ | 6,00 ± 1,70 | 70,5 ± 3,0 | **8,5 %** |
| outro hádron | 1,08 ± 0,30 | 12,8 ± 1,1 | **8,5 %** |
| fragmentação | 1,49 ± 0,37 | 2,08 ± 0,41 | **71,6 %** |

O isolamento guarda **71,6 %** dos fótons de fragmentação e só **8,5 %** dos de
decaimento. A fração de **fragmentação** na amostra sobe de **2,44 %** para
**17,41 %**, um ganho de **7,1×**.

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
| σ no gráfico, 5–40 GeV/c | 4,718 × 10⁻⁴ mb |

A fração isolada põe as quatro origens em ordem:

| origem | fração isolada |
|---|---|
| **direto** (sinal) | **87,2 %** |
| fragmentação | 71,6 % |
| π⁰ | 8,5 % |
| outro hádron | 8,5 % |

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
| `producao.sh` | roda as 4 faixas em paralelo (`./producao.sh` para os dois conjuntos) |
| `jatos_foton.cc` | versão v1, parte do jato (mantida para referência) |

### Análise e visualização

| arquivo | o que faz |
|---|---|
| `soma_pthat.cc` | soma as 4 faixas com teto, pesadas por σ/N (64 histogramas) |
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
| `soma_pthat.root` | a soma, 64 histogramas |
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
```

Precisa de ROOT com `TPythia8` e do FastJet no `ROOT_INCLUDE_PATH`.
