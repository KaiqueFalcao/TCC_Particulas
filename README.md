# TCC — Física de partículas com Pythia 8, FastJet e ROOT

Simulações de Monte Carlo de colisões próton-próton, usando Pythia 8 como
gerador de eventos, FastJet para reconstrução de jatos e ROOT para análise e
histogramas.

O repositório reúne três etapas do trabalho, da mais recente para a mais
antiga.

---

## [`FotonDireto/`](FotonDireto/) — projeto atual

Mede a **contaminação na identificação de fótons isolados**: de tudo que passa
num corte de isolamento, que fração é fóton direto e que fração é fundo de
decaimento de hádrons?

**Resposta: a pureza é 4,76 %** (24,22 % se a fragmentação contar como sinal).

Sinal e fundo são gerados como duas amostras separadas e comparados em
seção de choque:

| | processo | papel |
|---|---|---|
| `foton_direto.cc` | `PromptPhoton:all` | sinal — fótons do processo duro |
| `jatos_foton_v2.cc` | `HardQCD:all` | fundo — π⁰, outro hádron, fragmentação |

O trabalho inclui o corte de isolamento (ΣpT < 2 GeV/c em ΔR < 0,4), a
classificação de origem por verdade do MC, a marcação por massa invariante de
pares de fótons (picos do π⁰ e do η reproduzidos na terceira casa) e a curva
de pureza em pT.

A produção é fatiada em faixas de pThat e somada com peso σ/N. Os `.root`
estão versionados, então os resultados abrem sem precisar simular de novo.

Documentação detalhada em [`FotonDireto/README.md`](FotonDireto/README.md) e,
para a etapa atual, em
[`FotonDireto/jatos_foton/README.md`](FotonDireto/jatos_foton/README.md).

### Etapas anteriores, na raiz de `FotonDireto/`

Função de fragmentação do π⁰ e assimetria do decaimento π⁰ → γγ, que é a
origem do fundo irredutível medido acima.

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
- Binários e figuras geradas fora de `figuras/` não são versionados.
- Arquivos `.root` são ignorados por padrão, **com exceção** dos resultados de
  produção em `FotonDireto/jatos_foton/`, que são versionados de propósito
  para abrir em outra máquina sem repetir rodadas de 15 min a 1 h cada. As
  exceções estão listadas uma por linha no `.gitignore` daquela pasta — o git
  não aceita comentário no fim de uma linha de padrão.
