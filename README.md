# TCC — Física de partículas com Pythia 8, FastJet e ROOT

Simulações de Monte Carlo de colisões próton-próton, usando Pythia 8 como
gerador de eventos, FastJet para reconstrução de jatos e ROOT para análise e
histogramas.

O repositório reúne três etapas do trabalho, da mais recente para a mais
antiga.

---

## [`FotonDireto/`](FotonDireto/) — projeto atual

Reconstrução de jatos e medida da função de fragmentação do π⁰, como etapa
preparatória para uma análise de **fótons diretos**. Inclui o cálculo da
assimetria do decaimento π⁰ → γγ, que é a origem do fundo irredutível dessa
análise.

Documentação detalhada em [`FotonDireto/README.md`](FotonDireto/README.md).

| Arquivo | O que faz |
|---|---|
| `jatos_pi0.cc` | Pythia nativo, compilado com `g++`, salva `TTree` |
| `jatos_pi0_root.cc` | macro ROOT (`TPythia8`/`TParticle`), π⁰ estável |
| `jatos_pi0_decai.cc` | macro ROOT, π⁰ decaindo; assimetria e fóton líder |

---

## [`DeltaE/`](DeltaE/) — conservação de energia em colisões pp

Mede `ΔE = (energia inicial dos feixes) − (energia final dos feixes)` e
relaciona com multiplicidade e energia cinética.

O problema central: em `HardQCD` os prótons **não sobrevivem intactos**, então
"energia final do feixe" precisa de uma definição. Quatro métodos foram
testados.

| Arquivo | O que faz |
|---|---|
| `teste_investiga_protons.cpp` | comparação dos quatro métodos |
| `investiga_protons.cpp` | duas definições de energia final do feixe |
| `investiga_protons_status.cpp` | tentativa via status code |
| `teste_status.cc` | Pythia nativo, status codes reais e genealogia |
| `origem_particulas.cpp` | cada partícula final descende do feixe ou do processo duro? |

**Resultado principal (negativo, e por isso relevante):** ~59% das partículas
finais descendem do remanescente do feixe *e* do processo duro ao mesmo tempo.
Na hadronização, a string de cor liga os dois, e os hádrons nascem ao longo
dela. Atribuir um hádron final a um párton inicial não é uma pergunta bem
posta — os métodos cinemáticos são a via operacional.

Essa conclusão é exatamente a motivação para o conceito de **jato**, usado no
projeto seguinte.

---

## [`Distribuicoes/`](Distribuicoes/) — exercícios iniciais

Distribuições de energia, momento, rapidez e conservação por evento.

| Arquivo | O que faz |
|---|---|
| `teste.cpp` | distribuições de energia e momento |
| `teste_Ajustado.cpp` | momentos por partícula e conservação por evento |

---

## Requisitos

| Pacote | Versão usada |
|---|---|
| Pythia | 8.317 |
| FastJet | 3.5.1 |
| ROOT | 6.40.02 |

Variáveis de ambiente esperadas no `.bashrc`:

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

## Convenções

- Cada pasta tem suas figuras em `figuras/`.
- Macros ROOT rodam com `root -l arquivo.cc`; programas com `main()` são
  compilados com `g++`.
- Binários, arquivos `.root` e figuras geradas fora de `figuras/` não são
  versionados (ver `.gitignore`).
