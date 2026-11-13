# TCC — Fótons diretos e fragmentação do píon neutro

Estudo de Monte Carlo em colisões próton-próton a 14 TeV que mede a
**contaminação na identificação de fótons isolados**: de tudo que passa num
corte de isolamento, que fração é fóton direto e que fração é fundo?

O fundo dominante vem do decaimento π⁰ → γγ. O trabalho gera as duas amostras
separadamente — `PromptPhoton:all` para o sinal, `HardQCD:all` para o fundo —
e as compara em seção de choque.

**Resposta:** a pureza é **4,76 %**, ou **24,22 %** se a fragmentação contar
como sinal. A contaminação é de 95,24 %.

---

## Etapa atual: [`jatos_foton/`](jatos_foton/)

A análise parte do **fóton**, não do jato. São dois códigos irmãos com os
mesmos cortes: `jatos_foton_v2.cc` (`HardQCD:all`) gera o fundo e
`foton_direto.cc` (`PromptPhoton:all`) gera o sinal.

### Três resultados

**1. O isolamento funciona, e ordena as origens como deveria.**

| origem | fração isolada |
|---|---|
| fóton direto (sinal) | **87,2 %** |
| fragmentação | 71,6 % |
| π⁰ | 8,5 % |
| outro hádron | 8,5 % |

**2. A massa invariante marca parte do fundo — mas só parte.**

Pareando cada fóton energético com os outros fótons do evento, o espectro
reproduz o PDG na terceira casa: π⁰ em 0,1350 (PDG 0,13498) e η em 0,5470
(PDG 0,54786). A marcação remove 41,2 % dos fótons de π⁰ sem custar **nenhum**
fóton de fragmentação.

Mas **57,8 % dos fótons de π⁰ não têm parceiro visível nenhum**. Esses são o
fundo irredutível, e existem exatamente pelo motivo descrito mais abaixo: a
assimetria do decaimento é plana.

**3. A pureza é baixa.** 4,76 % na convenção em que só o fóton direto é sinal;
24,22 % contando a fragmentação junto.

A produção é fatiada em faixas de pThat e somada com peso σ/N, o que mede a
mesma física com 30× mais precisão pelo mesmo custo de CPU. Sinal e fundo usam
as **mesmas** faixas — comparar faixas diferentes estragaria a pureza em pT
baixo. Todos os `.root` estão versionados: dá para abrir os resultados sem
rodar nada.

**Limite conhecido:** o fundo tem só 116 eventos de Monte Carlo no espectro
inteiro. O valor integrado é confiável; a curva bin a bin acima de ~15 GeV/c
não é.

Documentação completa em [`jatos_foton/README.md`](jatos_foton/README.md).

---

## Contexto físico

A função de fragmentação mede a fração do momento do jato carregada por um
hádron:

```
z = pT(π⁰) / pT(jato)
```

Um π⁰ decai em dois fótons (BR = 98,8%). A **assimetria** do decaimento

```
A = |E₁ − E₂| / (E₁ + E₂)
```

é **uniforme em [0,1]**, porque o π⁰ tem spin 0 e decai isotropicamente. Isso
significa que uma fração constante e não desprezível dos decaimentos é muito
assimétrica: um fóton leva quase toda a energia e o outro fica soft demais para
ser detectado. Esse fóton solitário é indistinguível de um fóton direto — é o
**fundo irredutível** que motiva o trabalho.

Para um π⁰ de 10 GeV, o fóton mole fica com ~0,5 MeV.

---

## Requisitos

| Pacote | Versão usada |
|---|---|
| Pythia | 8.317 |
| FastJet | 3.5.1 |
| ROOT | 6.40.02 |

Variáveis de ambiente esperadas (no `.bashrc`):

```bash
export PYTHIA8=/caminho/para/pythia8317
export PYTHIA8DATA=$PYTHIA8/share/Pythia8/xmldoc
export FASTJET=/caminho/para/fastjet-install
export PATH=$FASTJET/bin:$PATH
export LD_LIBRARY_PATH=$FASTJET/lib:$PYTHIA8/lib:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=$PYTHIA8/include:$FASTJET/include
source /caminho/para/root/bin/thisroot.sh
```

O `ROOT_INCLUDE_PATH` e o `LD_LIBRARY_PATH` são o que permitem usar FastJet
dentro de macros ROOT sem `gSystem->Load()`.

---

## Programas

### `jatos_pi0.cc` — versão compilada, Pythia nativo

Usa a API nativa do Pythia (`Pythia8::Pythia`). O π⁰ decai, e é reencontrado
subindo a cadeia de mães. Salva histogramas e uma `TTree` num arquivo `.root`.

```bash
g++ jatos_pi0.cc -o jatos_pi0 -std=c++17 \
    -I$PYTHIA8/include -L$PYTHIA8/lib -Wl,-rpath,$PYTHIA8/lib -lpythia8 \
    $(fastjet-config --cxxflags --libs) $(root-config --cflags --libs)

./jatos_pi0 [nEventos] [pTHatMin] [eCM] [R] [pTJatoMin]
./jatos_pi0 5000 20 13000 0.4 20
```

### `jatos_pi0_root.cc` — macro ROOT, π⁰ estável

Usa `TPythia8` + `TClonesArray`/`TParticle`, sem Pythia nativo. Com
`111:mayDecay = off`, o π⁰ é partícula final e entra direto como constituinte
do jato — o caminho mais simples para medir a função de fragmentação.

```bash
root -l jatos_pi0_root.cc
```

Gera `jatos_pi0.png` (pT dos jatos, multiplicidade, pT dos π⁰, função de
fragmentação em log) e `mapa_eta_phi.png`.

### `jatos_pi0_decai.cc` — macro ROOT, π⁰ decaindo

Igual ao anterior, mas com o π⁰ decaindo (comportamento físico real). O π⁰ não
é mais constituinte: quem entra no jato são os dois fótons, e o π⁰ é
reencontrado pela função `mae_pi0()`, que sobe a cadeia de mães.

Além dos gráficos do anterior, calcula a **assimetria** do decaimento e a
**função de fragmentação do fóton líder**.

```bash
root -l jatos_pi0_decai.cc
```

Gera também `fotons_decai.png`.

---

## Parâmetros padrão

```
√s              = 14000 GeV
processo        = HardQCD:all
pTHatMin        = 20 GeV
algoritmo       = anti-kT, R = 0,4
pT mínimo jato  = 20 GeV
aceitação       = |η| < 0,9 (partículas), |y| < 0,5 (jatos)
semente         = 43
```

---

## Resultados (5000 eventos)

| | π⁰ estável | π⁰ decaindo |
|---|---|---|
| jatos reconstruídos | 881 | 885 |
| π⁰ dentro de jatos | 3122 | 3806 |
| π⁰ por jato | 3,54 | 4,30 |
| ⟨pT⟩ dos jatos | 30,6 GeV/c | 30,3 GeV/c |
| ⟨z⟩ dos π⁰ | 0,0682 | 0,0620 |
| ⟨A⟩ assimetria | — | 0,4919 |
| ⟨z⟩ do fóton líder | — | 0,0469 |

A diferença no número de π⁰ é metodológica, não um erro: com o π⁰ estável ele
precisa estar ele próprio dentro do jato; com o π⁰ decaindo ele é contado se
**pelo menos um** dos dois fótons for constituinte — critério mais frouxo.

### Validações físicas

Três testes independentes confirmam que a navegação do registro de eventos está
correta:

1. **98,7%** dos π⁰ produziram par γγ — o PDG diz BR(γγ) = 98,8%. O resto é
   Dalitz (π⁰ → γe⁺e⁻), corretamente descartado.
2. **Assimetria plana**: ⟨A⟩ = 0,4919 e RMS = 0,2814 ≈ 1/√12 = 0,2887.
3. **⟨z_γ⟩ / ⟨z_π⁰⟩ = 0,756**, contra a previsão ⟨(1+A)/2⟩ = **0,746**.

---

## Figuras

Em [`figuras/`](figuras/). As imagens geradas na raiz ao rodar são ignoradas
pelo `.gitignore`; copie para `figuras/` as versões que quiser versionar.

---

## Etapas concluídas

- [x] Seção de choque normalizada por `sigmaGen()`/`nAccepted()` — o peso σ/N
      de cada rodada fica no bin 5 do `hist_sigma`, dentro do próprio `.root`
- [x] Massa invariante combinatória de pares de fótons — picos do π⁰ e do η
      reproduzidos na terceira casa
- [x] π⁰-tagging — remove 41,2 % dos fótons de π⁰ sem custar fragmentação;
      57,8 % ficam sem parceiro visível e são irredutíveis
- [x] Corte de isolamento (ΣpT < 2 GeV/c num cone ΔR < 0,4)
- [x] Classificação por verdade do MC: direto / fragmentação / decaimento /
      outro hádron, subindo pelas cópias do chuveiro
- [x] Pureza vs. pT, sinal e fundo nas mesmas faixas de pThat

## O que falta

- [ ] **Mais estatística de `HardQCD`.** O fundo tem 116 eventos de MC contra
      27 485 do sinal. Para a curva de pureza valer acima de 15 GeV/c seriam
      precisos uns 10× mais eventos, cerca de 4 h de CPU.
- [ ] Eficiência do corte de isolamento **em função do pT** (hoje só o valor
      integrado por origem)
- [ ] Teste de fechamento formal da soma fatiada (hoje há a conferência das σ:
      141,577 mb somados contra 140,818 mb da rodada inclusiva, 0,54 %)
- [ ] Incerteza sistemática do limiar de 0,5 GeV/c no parceiro, que é o
      parâmetro que manda na eficiência da marcação por massa
