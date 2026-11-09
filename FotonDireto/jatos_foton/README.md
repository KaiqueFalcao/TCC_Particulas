# Fótons isolados em jatos — produção fatiada em pThat

Mede quanto a **técnica de isolamento** consegue limpar uma amostra de fótons
energéticos, separando os de **fragmentação** (sinal) dos de **decaimento de
hádrons** (fundo).

Pythia 8 + FastJet dentro do ROOT (`TPythia8`, `TClonesArray`, `TParticle`).
pp a 14 TeV, `HardQCD:all`, anti-kT com R = 0,4.

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
decaimento. A pureza do sinal sobe de **2,44 %** para **17,41 %**, um ganho de
**7,1×**.

---

## Como abrir os resultados sem rodar nada

Todos os `.root` estão versionados. Entre na pasta e escolha o que quer ver:

| quero ver | comando |
|---|---|
| a soma, 12 canvases | `root -l soma_pthat.cc` |
| a soma com eixo em mb | `root -l 'soma_pthat.cc(0)'` |
| as 4 faixas + a soma | `root -l ve_faixas.cc` |
| as 4 faixas sem a soma | `root -l 've_faixas.cc(1)'` |
| as 4 sobrepostas (com teto) | `root -l abre_todos.cc` |
| as 4 sobrepostas (sem teto) | `root -l 'abre_todos.cc(kTRUE, kFALSE)'` |
| o teste de consistência | `root -l compara_pthat.cc` |
| uma faixa por inteiro | `root -l 'abre_jatos_foton_v2.cc("jatos_foton_pthat15_inf.root")'` |

Nenhum deles simula; todos leem os `.root` em segundos.

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

Quase metade do sinal vem de 1,6 % dos eventos. Sobre-amostrar o pThat alto e
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
| `jatos_foton_v2.cc` | a simulação — parte do fóton, classifica origem, isola |
| `producao.sh` | roda as 4 faixas em paralelo (`./producao.sh` para os dois conjuntos) |
| `jatos_foton.cc` | versão v1, parte do jato (mantida para referência) |

### Análise e visualização

| arquivo | o que faz |
|---|---|
| `soma_pthat.cc` | soma as 4 faixas com teto, pesadas por σ/N |
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
| `soma_pthat.root` | a soma, 61 histogramas |
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

---

## Reproduzir do zero

```bash
./producao.sh              # os dois conjuntos, 100 mil cada (~25 min)
./producao.sh 5000         # teste rápido
./producao.sh 100000 teto  # só as faixas com teto
root -l -b -q soma_pthat.cc
```

Precisa de ROOT com `TPythia8` e do FastJet no `ROOT_INCLUDE_PATH`.
