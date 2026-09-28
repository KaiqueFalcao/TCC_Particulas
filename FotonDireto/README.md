# TCC — Fótons diretos e fragmentação do píon neutro

Estudo de Monte Carlo em colisões próton-próton para a reconstrução de jatos e
a medida da função de fragmentação do π⁰, como etapa preparatória para uma
análise de **fótons diretos**.

O fundo dominante para fótons diretos vem do decaimento π⁰ → γγ. Este
repositório constrói as ferramentas para quantificar esse fundo.

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

## Próximas etapas

- [ ] Seção de choque do fóton líder (normalizar por `sigmaGen()`/`nAccepted()`)
- [ ] Massa invariante combinatória de pares de fótons (picos do π⁰ e do η)
- [ ] π⁰-tagging: remover fótons que formam massa de π⁰ com algum parceiro
- [ ] Corte de isolamento (ΣpT num cone ΔR < 0,4)
- [ ] Classificação por verdade do MC: direto / fragmentação / decaimento
- [ ] Pureza e eficiência vs. pT; teste de fechamento
