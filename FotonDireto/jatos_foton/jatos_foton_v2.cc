// Contaminacao de fotons energeticos -- copia do jatos_foton.cc.
//
// MUDANCA DE PONTO DE VISTA: o jatos_foton.cc parte do JATO e classifica os
// fotons que sao constituintes dele. Aqui a analise parte do FOTON: varremos
// TODAS as particulas finais do evento atras de fotons energeticos, e so'
// depois procuramos em que jato cada um caiu.
//
// A pergunta: entre os fotons energeticos que carregam quase todo o pT do
// proprio jato, quantos sao contaminacao de decaimento (pi0 e outros hadrons)
// e quantos sao de fragmentacao? E o quanto o isolamento melhora essa razao?
//
// CONJUNTOS DE HISTOGRAMAS:
//   C) lider x nao-lider dos pi0, sem corte            (6 faixas)  -- como antes
//   D) lider x nao-lider dos pi0, z > corte_z          (6 faixas)  -- como antes
//   E) lider x nao-lider dos pi0, z > corte_z e A      (6 faixas)  -- como antes
//   G) z dos fotons energeticos ISOLADOS,     por origem (3 x 6 faixas)  -- NOVO
//   H) z dos fotons energeticos NAO ISOLADOS, por origem (3 x 6 faixas)  -- NOVO
//
// AS TRES ORIGENS:
//   0 = pi0 + eta       decaimento de pi0 (PDG 111) OU de eta (PDG 221)
//   1 = outro hadron    qualquer outro hadron, JA' SEM o eta
//   2 = fragmentacao    foton irradiado por um quark no chuveiro
//
// O eta foi tirado de "outro hadron" e juntado ao pi0 porque ele faz a MESMA
// coisa que o pi0: decai em dois fotons (eta -> gama gama, cerca de 39% das
// vezes) e, quando um dos dois fica mole, o outro se disfarca de foton direto.
// E' o mesmo mecanismo de contaminacao, e a marcacao por massa invariante
// pega os dois do mesmo jeito -- basta olhar tambem a janela de massa do eta.
//
// Os conjuntos A, B e F do jatos_foton.cc foram removidos: G e H os substituem,
// agora fatiados por pT do FOTON e separados pelas tres origens.
//
// As faixas sao as mesmas de antes: 5-8, 8-10, 10-12, 12-16, 16-20 e >20 GeV.
// A ultima e' implementada como 20 a 1000, ou seja, "acima de 20".
//
// ============================================================================
// RODADAS FATIADAS EM pThat
// ============================================================================
// Foton com pT alto e' raro: numa unica rodada com pTHatMin = 4 quase todo
// evento e' de pThat baixo, e as faixas de cima ficam vazias. A solucao e'
// rodar varias vezes em faixas de pThat e depois SOMAR com peso.
//
// ATENCAO -- as faixas tem que ser EXCLUSIVAS. pTHatMin = 4 nao gera eventos
// "perto de 4": gera TODOS os eventos com pThat > 4, inclusive os de pThat
// > 15. Se usarmos so' pTHatMin = 4, 7, 10 e 15, as amostras ficam encaixadas
// uma na outra, e somar com peso sigma_i/N_i conta a regiao de pThat alto
// quatro vezes. Por isso existe o parametro pt_hat_max: ele fecha a faixa por
// cima com PhaseSpace:pTHatMax, e as quatro amostras passam a ser disjuntas.
//
// Producao padrao (4 processos em paralelo, sementes diferentes):
//   root -l -b -q 'jatos_foton_v2.cc( 4.,  7., 43, 100000, "pthat04_07")'
//   root -l -b -q 'jatos_foton_v2.cc( 7., 10., 44, 100000, "pthat07_10")'
//   root -l -b -q 'jatos_foton_v2.cc(10., 15., 45, 100000, "pthat10_15")'
//   root -l -b -q 'jatos_foton_v2.cc(15., -1., 46, 100000, "pthat15_inf")'
// pt_hat_max = -1 quer dizer "sem teto", usado so' na ultima faixa.
//
// COMO SOMAR DEPOIS: cada rodada salva o seu proprio .root e, dentro dele, o
// histograma hist_sigma com a secao de choque e o peso da rodada. O espectro
// combinado e' a soma ponderada
//      h_total = soma_i ( sigma_i / N_i ) * h_i
// com sigma_i em mb lido de hist_sigma. Como as faixas sao exclusivas, essa
// soma reproduz o espectro de pThat > 4 com a estatistica reforcada no alto.
// Os histogramas tem Sumw2 ligado, entao as incertezas acompanham o Scale().

#include "TSystem.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TPythia8.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include <vector>
#include <set>

const Int_t NFAIXAS = 6;

// Massas do pi0 e do eta, PDG 2024. Sao os dois alvos da marcacao por massa
// invariante: os dois decaem em dois fotons e os dois contaminam a amostra.
const Double_t M_PI0 = 0.1349768;   // GeV
const Double_t M_ETA = 0.547862;    // GeV
const Int_t NORIGENS = 3;   // 0 = pi0 + eta, 1 = outro hadron, 2 = fragmentacao

// Ponteiros em escopo de arquivo: continuam acessiveis no prompt do ROOT
// depois que a macro termina, entao da' para mexer nos histogramas.
TH1F   *hist_pt_jato     = 0;
TH1F   *hist_n_part      = 0;
TH1F   *hist_pt_foton    = 0;
TH2F   *hist_eta_phi     = 0;
TH1F   *hist_assimetria  = 0;
TH1F   *hist_z_gama      = 0;   // foton LIDER do decaimento do pi0 (como antes)
TH1F   *hist_soma_cone   = 0;   // soma de pT no cone, dos fotons energeticos

// Espectro de pT dos fotons ISOLADOS, um histograma por origem. Mesmo
// binning do hist_pt_direto_iso da foton_direto.cc, para o RUIDO de cada
// origem poder ser sobreposto ao SINAL dos fotons diretos.
TH1F   *hist_pt_iso_pi0eta = 0;
TH1F   *hist_pt_iso_frag   = 0;
TH1F   *hist_pt_iso_outro  = 0;

// E_iso = soma de pT no cone de isolamento, um histograma por origem, para
// TODO foton energetico (isolado ou nao). E' a distribuicao cujo pedaco
// abaixo do corte de 2 GeV/c e' a eficiencia do isolamento daquela origem.
// A foton_direto.cc tem o mesmo histograma para o foton direto, com o mesmo
// binning, para as quatro curvas poderem ir ao mesmo grafico.
TH1F   *hist_eiso[NORIGENS];

// MASSA INVARIANTE dos pares de fotons. Serve para marcar e descartar o
// foton energetico que, na verdade, e' filho de um pi0 cujo irmao ficou
// mole demais para ser notado.
TH1F   *hist_massa_gg        = 0;   // todos os pares, todos os energeticos
TH1F   *hist_massa_melhor    [NORIGENS];  // par mais proximo do pi0, isolados
TH1F   *hist_massa_melhor_eta[NORIGENS];  // par mais proximo do eta, isolados

// Espectro de pT DEPOIS da marcacao: os fotons marcados como filhos de pi0
// sao jogados fora. E' a amostra que sobra para a analise de verdade.
TH1F   *hist_pt_pos_pi0eta = 0;
TH1F   *hist_pt_pos_frag   = 0;
TH1F   *hist_pt_pos_outro  = 0;

// Informacao da rodada: secao de choque, numero de eventos e peso.
// TH1D (nao TH1F) porque sigma em mb precisa de precisao dupla.
TH1D   *hist_sigma       = 0;

// (C) sem corte ; (D) z > corte_z ; (E) z e A
TH2F   *hist_td[NFAIXAS] = {0, 0, 0, 0, 0, 0};
TH2F   *hist_Z [NFAIXAS] = {0, 0, 0, 0, 0, 0};
TH2F   *hist_AZ[NFAIXAS] = {0, 0, 0, 0, 0, 0};

// (G) e (H): z dos fotons energeticos, por origem e por faixa de pT do foton
TH1F   *hist_z_iso   [NORIGENS][NFAIXAS];
TH1F   *hist_z_naoiso[NORIGENS][NFAIXAS];

TCanvas *c1  = 0;
TCanvas *c2  = 0;
TCanvas *c3  = 0;
TCanvas *c4  = 0;
TCanvas *c5  = 0;
TCanvas *c6  = 0;
TCanvas *c7  = 0;
TCanvas *c8  = 0;
TCanvas *c9  = 0;
TCanvas *c10 = 0;
TCanvas *c11 = 0;
TCanvas *c12 = 0;
TCanvas *c13 = 0;   // massa invariante
TCanvas *c14 = 0;   // E_iso por origem

Bool_t eh_neutrino(Int_t pdg) // devolve true se o PDG for de neutrino
{
  Int_t a = TMath::Abs(pdg);
  return (a == 12 || a == 14 || a == 16 || a == 18);
}

// Soma do pT das particulas finais num cone de raio R_iso em volta do foton
// de indice i_foton. Percorre o evento INTEIRO, nao so' o jato: o isolamento
// olha tudo que esta' em volta do foton, esteja ou nao no mesmo jato.
Double_t soma_cone(TClonesArray *lista, Int_t i_foton,
                   Double_t R_iso, Double_t eta_part_max)
{
  Int_t n = lista->GetEntriesFast();

  TParticle *foton  = (TParticle*) lista->At(i_foton);
  Double_t eta_foton = foton->Eta();
  Double_t phi_foton = foton->Phi();

  Double_t soma_pt_cone = 0;

  for (Int_t i = 0; i < n; i++)
  {
    if (i == i_foton) continue;                       // o proprio foton

    TParticle *p = (TParticle*) lista->At(i);

    if (p->GetStatusCode() <= 0) continue;            // intermediaria: contaria a mesma energia duas vezes
    if (eh_neutrino(p->GetPdgCode())) continue;
    if (TMath::Abs(p->Eta()) > eta_part_max) continue;

    Double_t d_eta = p->Eta() - eta_foton;

    // Phi() devolve [0, 2pi); a diferenca e' trazida para [-pi, pi]
    Double_t d_phi = p->Phi() - phi_foton;
    if (d_phi >  TMath::Pi()) d_phi -= 2 * TMath::Pi();
    if (d_phi < -TMath::Pi()) d_phi += 2 * TMath::Pi();

    Double_t R_i = TMath::Sqrt(d_eta * d_eta + d_phi * d_phi);
    if (R_i > R_iso) continue;

    soma_pt_cone += p->Pt();
  }
  return soma_pt_cone;
}

// Massa invariante de um par de fotons.
//
//   m^2 = (E1 + E2)^2 - |p1 + p2|^2
//
// Para dois fotons, que nao tem massa, isso se reduz a
//   m^2 = 2 E1 E2 (1 - cos(theta))
// com theta o angulo entre eles. Ou seja: a massa do par e' grande quando
// os dois sao energeticos E abrem bastante. Um pi0 muito energetico decai
// em dois fotons quase colineares, e esses pares de theta pequeno sao
// justamente os dificeis de separar num calorimetro real.
//
// A conta e' feita somando os quadrivetores, que e' menos sujeita a erro
// do que a formula do angulo.
//
// A mesma funcao serve para o eta: o que muda e' so' o valor de massa que a
// gente procura (0,548 em vez de 0,135).
Double_t massa_par(TParticle *a, TParticle *b)
{
  Double_t E  = a->Energy() + b->Energy();
  Double_t px = a->Px()     + b->Px();
  Double_t py = a->Py()     + b->Py();
  Double_t pz = a->Pz()     + b->Pz();

  Double_t m2 = E*E - px*px - py*py - pz*pz;

  // m2 pode dar levemente negativo por arredondamento quando os dois
  // fotons sao quase colineares. Nesse caso a massa e' zero, nao NaN.
  return (m2 > 0) ? TMath::Sqrt(m2) : 0.;
}

// Devolve em qual faixa um valor de pT cai (-1 se estiver fora de todas)
Int_t acha_faixa(Double_t pt, Double_t *lo, Double_t *hi)
{
  for (Int_t f = 0; f < NFAIXAS; f++)
    if (pt >= lo[f] && pt < hi[f]) return f;
  return -1;
}

// Classifica a origem do foton pelo PDG da MAE DIRETA.
// GetFirstMother() devolve um INDICE na lista, nao um PDG.
//    0 = decaimento de pi0 OU de eta (|PDG| == 111 ou == 221)
//    1 = decaimento de outro hadron  (|PDG| > 100, sem o pi0 e sem o eta)
//    2 = fragmentacao                (|PDG| de 1 a 6: irradiado por um quark)
//   -1 = nao classificado: sem mae valida, ou mae que nao e' nenhum dos tres
//        (gluon, lepton, boson). Esses sao contados a parte, nao descartados
//        em silencio.
//
// POR QUE E' PRECISO SUBIR PELAS COPIAS DO PROPRIO FOTON:
// quando um quark irradia um foton no chuveiro, o Pythia registra o foton
// tendo o quark como mae. Mas o chuveiro CONTINUA, e sempre que esse foton
// participa de uma emissao seguinte como particula de recuo ele e' RECOPIADO:
// a copia nova tem o foton antigo (PDG 22) como mae. Olhando so' a mae direta,
// esses fotons nao sao pi0, nem hadron, nem quark, e caem todos em "nao
// classificado" -- a categoria fragmentacao fica vazia por construcao.
// Medido com 50 mil eventos: os 2 fotons energeticos sem origem tinham mae
// PDG 22 e, subida a cadeia, ancestral PDG +-4 (quark charm). Nenhum foton
// energetico tinha quark como mae DIRETA, ou seja, todos passaram por copia.
// Por isso aqui subimos enquanto a mae for ela mesma um foton, e so' entao
// classificamos pelo PDG do primeiro ancestral que nao e' foton.
//
// A ordem dos testes importa: 111 e 221 sao maiores que 100, entao o pi0 e o
// eta tem que ser testados ANTES do caso geral de hadron -- senao os dois
// cairiam em "outro hadron".
Int_t origem_foton(TClonesArray *lista, Int_t i_foton)
{
  Int_t n = lista->GetEntriesFast();

  TParticle *foton = (TParticle*) lista->At(i_foton);
  Int_t mae = foton->GetFirstMother();

  // sobe pelas copias: enquanto a mae for um foton, passa para a mae dela.
  // A trava de 50 passos protege contra qualquer laco no registro do evento.
  Int_t passos = 0;
  while (mae >= 0 && mae < n && passos < 50)
  {
    Int_t pdg_mae = TMath::Abs(((TParticle*) lista->At(mae))->GetPdgCode());
    if (pdg_mae != 22) break;           // chegou num ancestral que nao e' foton
    mae = ((TParticle*) lista->At(mae))->GetFirstMother();
    passos++;
  }

  if (mae < 0 || mae >= n) return -1;

  Int_t pdg = TMath::Abs(((TParticle*) lista->At(mae))->GetPdgCode());

  if (pdg == 111 || pdg == 221) return 0;   // pi0 ou eta
  if (pdg >= 1 && pdg <= 6)     return 2;   // quark -> fragmentacao
  if (pdg > 100)                return 1;   // outro hadron (ja' sem o eta)
  return -1;                            // inclui o caso de estourar a trava
}

// pt_hat_min / pt_hat_max : faixa de pThat desta rodada (pt_hat_max = -1
//                           quer dizer sem teto, so' para a ultima faixa)
// semente                 : Random:seed, diferente em cada rodada
// nev                     : eventos a gerar
// sufixo                  : entra no nome do .root e dos PNGs, para as
//                           rodadas em paralelo nao se sobreporem. Deixando
//                           vazio, e' montado a partir do pTHatMin.
void jatos_foton_v2(Double_t pt_hat_min = 4.,
                    Double_t pt_hat_max = -1.,
                    Int_t    semente    = 43,
                    Int_t    nev        = 100000,
                    const char *sufixo  = "")
{
  // ---- parametros
  Double_t R            = 0.4;
  Double_t pt_jato_min  = 5.;      // casa com a primeira faixa
  Double_t eta_part_max = 0.9;
  Double_t y_jato_max   = 0.9 - R;
  Double_t energia      = 14000.;
  Double_t corte_A      = 0.9;     // assimetria minima
  Double_t corte_z      = 0.8;     // fracao minima do pT do jato
  Double_t R_iso        = 0.4;     // raio do cone de isolamento
  Double_t pt_iso_max   = 2.0;     // GeV/c: isolado se a soma no cone for menor
  Double_t pt_foton_min = 5.0;     // GeV/c: define "foton energetico"

  // ---- marcacao por massa invariante (parametros NOVOS, nenhum corte
  // ---- antigo foi tocado)
  //
  // pt_parceiro_min e' o que mais importa aqui, e merece atencao.
  // O Pythia entrega TODOS os fotons, inclusive os de 1 MeV. Um
  // calorimetro real nao: abaixo de alguns centesimos de GeV o deposito
  // some no ruido. Se deixassemos o parceiro ser arbitrariamente mole,
  // acharíamos o irmao de QUALQUER pi0 e concluiriamos que a contaminacao
  // e' toda removivel -- o que e' falso, e seria falso justamente pelo
  // motivo que da' origem ao problema: o fundo irredutivel existe PORQUE
  // o irmao e' mole demais para ser visto.
  // 0,5 GeV/c e' a ordem de grandeza de um limiar de cluster tipico.
  Double_t pt_parceiro_min = 0.5;  // GeV/c: limiar para "ver" o parceiro
  Double_t janela_pi0      = 0.02; // GeV: meia-largura da janela do pi0
  Double_t janela_eta      = 0.02; // GeV: meia-largura da janela do eta

  // O foton e' MARCADO se existir um parceiro visivel cujo par caia em
  // QUALQUER UMA das duas janelas -- a do pi0 ou a do eta. As duas janelas
  // nao se encostam (0,115-0,155 e 0,528-0,568), entao nao ha' ambiguidade:
  // um mesmo par nunca cai nas duas.

  // ---- E_iso: binning do histograma da soma de pT no cone, por origem.
  // 0 a 20 GeV/c em 100 bins = 0,2 GeV por bin. A largura foi escolhida para
  // que o corte de isolamento (2,0 GeV/c) caia exatamente numa BORDA de bin:
  // assim a integral abaixo do corte e' a eficiencia do isolamento, sem
  // partir bin nenhum ao meio.
  Double_t eiso_max   = 20.;
  Int_t    eiso_nbins = 100;

  // seis faixas; a ultima vai ate 1000, ou seja, "acima de 20 GeV"
  Double_t faixa_lo[NFAIXAS] = { 5.,  8., 10., 12., 16.,   20.};
  Double_t faixa_hi[NFAIXAS] = { 8., 10., 12., 16., 20., 1000.};

  // nomes das tres origens, usados nos titulos e no resumo
  const char *nome_origem[NORIGENS] = {"#pi^{0} + #eta", "outro hadron",
                                       "fragmentacao"};
  const char *curto_origem[NORIGENS] = {"pi0eta", "outro", "frag"};

  // ---- buffers para montar textos com a familia printf
  char comando[200], nome_arquivo[200], png[200], rodada[60];

  // ---- nome desta rodada. Com sufixo vazio ele e' montado a partir do
  // pTHatMin -- "pthat04", "pthat07", "pthat10", "pthat15" -- entao duas
  // rodadas de pThat diferente nunca escrevem no mesmo arquivo, mesmo
  // esquecendo de batizar cada uma. Para pTHatMin nao inteiro, passe o
  // sufixo na mao: o %02.0f arredonda (7,5 viraria "08").
  if (sufixo[0] != '\0') sprintf(rodada, "%s", sufixo);
  else                    sprintf(rodada, "pthat%02.0f", pt_hat_min);

  // ---- arquivo de saida: criado ANTES dos histogramas, para que todos
  // ---- sejam associados a ele e salvos automaticamente no Write().
  // O sufixo no nome e' o que permite rodar as 4 faixas ao mesmo tempo.
  sprintf(nome_arquivo, "jatos_foton_%s.root", rodada);
  TFile *arquivo = new TFile(nome_arquivo, "RECREATE");

  // ---- Pythia: pTHatMin, pTHatMax e semente vem dos parametros
  TPythia8 *pythia8 = new TPythia8();
  pythia8->ReadString("HardQCD:all = on");
  pythia8->ReadString("Random:setSeed = on");

  sprintf(comando, "PhaseSpace:pTHatMin = %.1f", pt_hat_min);
  pythia8->ReadString(comando);

  // pt_hat_max < 0 quer dizer "sem teto": nao manda pTHatMax nenhum,
  // e o Pythia usa o limite cinematico da colisao.
  if (pt_hat_max > 0)
  {
    sprintf(comando, "PhaseSpace:pTHatMax = %.1f", pt_hat_max);
    pythia8->ReadString(comando);
  }

  sprintf(comando, "Random:seed = %d", semente);
  pythia8->ReadString(comando);

  pythia8->Initialize(2212, 2212, energia);

  TClonesArray *particulas = new TClonesArray("TParticle", 10000);

  // ---- FastJet
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
  std::vector<fastjet::PseudoJet> entradas;

  // ---- histogramas gerais
  char nome[100], titulo[250], rotulo[60];

  // Sumw2 em TODOS os histogramas criados daqui para baixo. E' indispensavel
  // para a soma ponderada: sem ele, o ROOT recalcula a incerteza como
  // sqrt(conteudo) e o Scale(sigma/N) estraga as barras de erro. Com ele, a
  // soma dos quadrados dos pesos e' guardada e escalada junto.
  TH1::SetDefaultSumw2(kTRUE);

  // Informacao desta rodada, preenchida DEPOIS do laco (sigmaGen converge ao
  // longo da geracao, entao so' o valor final vale). Bins com rotulo para
  // dar para ler no prompt com hist_sigma->GetXaxis()->GetBinLabel(i).
  hist_sigma = new TH1D("hist_sigma",
      "informacao da rodada;;valor", 7, 0, 7);
  hist_sigma->GetXaxis()->SetBinLabel(1, "sigma_mb");
  hist_sigma->GetXaxis()->SetBinLabel(2, "sigma_err_mb");
  hist_sigma->GetXaxis()->SetBinLabel(3, "n_aceitos");
  hist_sigma->GetXaxis()->SetBinLabel(4, "n_pedidos");
  hist_sigma->GetXaxis()->SetBinLabel(5, "peso_mb_por_ev");
  hist_sigma->GetXaxis()->SetBinLabel(6, "pt_hat_min");
  hist_sigma->GetXaxis()->SetBinLabel(7, "pt_hat_max");

  hist_pt_jato  = new TH1F("hist_pt_jato",
      "pT dos jatos;pT^{jato} [GeV/c];jatos", 100, 0, 100);
  hist_n_part   = new TH1F("hist_n_part",
      "particulas por jato;N_{const};jatos", 60, 0, 60);
  hist_pt_foton = new TH1F("hist_pt_foton",
      "pT dos #gamma;pT^{#gamma} [GeV/c];#gamma", 100, 0, 40);
  hist_eta_phi  = new TH2F("hist_eta_phi",
      "mapa #eta-#phi dos #gamma;#eta;#phi [rad]",
      60, -eta_part_max, eta_part_max, 60, -TMath::Pi(), TMath::Pi());

  hist_assimetria = new TH1F("hist_assimetria",
      "assimetria do decaimento;A = |E_{1}-E_{2}|/(E_{1}+E_{2});#pi^{0}",
      50, 0, 1);
  hist_z_gama = new TH1F("hist_z_gama",
      "fragmentacao do foton lider;z_{#gamma} = pT^{#gamma}/pT^{jato};#gamma",
      50, 0, 1);

  // isolamento: 0,5 GeV por bin, entao o corte de 2 GeV cai numa borda de bin
  sprintf(titulo, "soma de pT no cone, #gamma energeticos (R < %.1f);"
      "#Sigma pT^{cone} [GeV/c];#gamma", R_iso);
  hist_soma_cone = new TH1F("hist_soma_cone", titulo, 80, 0, 40);

  // ---- espectro de pT dos fotons ISOLADOS, um histograma por origem.
  // 35 bins de 5 a 40 GeV/c, ou seja 1 GeV por bin. E' EXATAMENTE o binning
  // do hist_pt_direto_iso da foton_direto.cc, de proposito: assim o ruido
  // de cada origem e o sinal dos fotons diretos podem ser desenhados no
  // mesmo grafico, bin por bin, sem reamostrar nada.
  hist_pt_iso_pi0eta = new TH1F("hist_pt_iso_pi0eta",
      "#gamma isolado de #pi^{0} + #eta;p_{T}^{#gamma} [GeV/c];#gamma",
      35, 5., 40.);
  hist_pt_iso_frag   = new TH1F("hist_pt_iso_frag",
      "#gamma isolado de fragmentacao;p_{T}^{#gamma} [GeV/c];#gamma",
      35, 5., 40.);
  hist_pt_iso_outro  = new TH1F("hist_pt_iso_outro",
      "#gamma isolado de outro hadron;p_{T}^{#gamma} [GeV/c];#gamma",
      35, 5., 40.);

  // ---- E_iso: soma de pT no cone, um histograma por origem.
  // Preenchido para TODO foton energetico, isolado ou nao: e' a forma da
  // distribuicao que interessa, e e' ela que explica por que o isolamento
  // separa fragmentacao de decaimento.
  for (Int_t o = 0; o < NORIGENS; o++)
  {
    sprintf(nome, "hist_eiso_%s", curto_origem[o]);
    sprintf(titulo, "E_{iso} de %s;#Sigma pT^{cone} (R < %.1f) [GeV/c];#gamma",
            nome_origem[o], R_iso);
    hist_eiso[o] = new TH1F(nome, titulo, eiso_nbins, 0., eiso_max);
  }

  // ---- massa invariante dos pares de fotons.
  // 500 bins de 0 a 1 GeV, 2 MeV por bin. A largura foi escolhida para que
  // a massa do pi0 (0,1349768) caia no CENTRO de um bin, [0,134 ; 0,136]:
  // com 5 MeV ela caia bem na borda e o pico aparecia deslocado para
  // 0,1325, o que confunde na hora de ler o grafico. O alcance ate' 1 GeV
  // mostra tambem o eta em 0,548, que agora e' alvo da marcacao junto com
  // o pi0. A massa do eta (0,547862) tambem cai no centro de um bin.
  hist_massa_gg = new TH1F("hist_massa_gg",
      "massa invariante de pares de #gamma;m_{#gamma#gamma} [GeV/c^{2}];pares", 500, 0., 1.);

  // Dois conjuntos: um com o par mais proximo da massa do PI0 e outro com o
  // par mais proximo da massa do ETA. Sao coisas diferentes -- o parceiro que
  // melhor imita um pi0 nao e' o mesmo que melhor imita um eta -- e e' por
  // isso que cada janela precisa do seu proprio histograma.
  for (Int_t o = 0; o < NORIGENS; o++)
  {
    sprintf(nome, "hist_massa_melhor_%s", curto_origem[o]);
    sprintf(titulo, "melhor par (#pi^{0}) de %s, #gamma ISOLADO;m_{#gamma#gamma} [GeV/c^{2}];#gamma", nome_origem[o]);
    hist_massa_melhor[o] = new TH1F(nome, titulo, 500, 0., 1.);

    sprintf(nome, "hist_massa_melhor_eta_%s", curto_origem[o]);
    sprintf(titulo, "melhor par (#eta) de %s, #gamma ISOLADO;m_{#gamma#gamma} [GeV/c^{2}];#gamma", nome_origem[o]);
    hist_massa_melhor_eta[o] = new TH1F(nome, titulo, 500, 0., 1.);
  }

  // ---- espectro de pT DEPOIS da marcacao, mesmo binning do de antes.
  // Par a par com os hist_pt_iso_*: um e' a amostra isolada, o outro e' a
  // mesma amostra com os marcados removidos. A razao entre os dois e' a
  // eficiencia do corte de massa em funcao do pT.
  hist_pt_pos_pi0eta = new TH1F("hist_pt_pos_pi0eta",
      "#gamma de #pi^{0} + #eta, isolado e NAO marcado;p_{T}^{#gamma} [GeV/c];#gamma", 35, 5., 40.);
  hist_pt_pos_frag   = new TH1F("hist_pt_pos_frag",
      "#gamma de fragmentacao, isolado e NAO marcado;p_{T}^{#gamma} [GeV/c];#gamma", 35, 5., 40.);
  hist_pt_pos_outro  = new TH1F("hist_pt_pos_outro",
      "#gamma de outro hadron, isolado e NAO marcado;p_{T}^{#gamma} [GeV/c];#gamma", 35, 5., 40.);

  // ---- histogramas fatiados por faixa.
  // Nome e titulo sao montados com sprintf porque mudam de faixa para faixa.
  // A ultima faixa recebe o rotulo "> 20" em vez de "20 - 1000".
  for (Int_t f = 0; f < NFAIXAS; f++)
  {
    // --- conjuntos C, D e E: rotulo pelo pT do JATO (como antes)
    if (f == NFAIXAS - 1) sprintf(rotulo, "pT^{jato} > %.0f GeV", faixa_lo[f]);
    else sprintf(rotulo, "%.0f < pT^{jato} < %.0f GeV", faixa_lo[f], faixa_hi[f]);

    // (C) lider x nao-lider, sem corte
    sprintf(nome, "hist_td_%d", f);
    sprintf(titulo, "%s, todos os #pi^{0};"
        "pT^{#gamma lider} [GeV/c];pT^{#gamma nao-lider} [GeV/c]", rotulo);
    hist_td[f] = new TH2F(nome, titulo, 50, 0, 20, 50, 0, 20);

    // (D) lider x nao-lider, so z
    sprintf(nome, "hist_Z_%d", f);
    sprintf(titulo, "%s, z > %.1f;"
        "pT^{#gamma lider} [GeV/c];pT^{#gamma nao-lider} [GeV/c]",
        rotulo, corte_z);
    hist_Z[f] = new TH2F(nome, titulo, 50, 0, 20, 50, 0, 20);

    // (E) lider x nao-lider, z e assimetria
    sprintf(nome, "hist_AZ_%d", f);
    sprintf(titulo, "%s, z > %.1f e A > %.1f;"
        "pT^{#gamma lider} [GeV/c];pT^{#gamma nao-lider} [GeV/c]",
        rotulo, corte_z, corte_A);
    hist_AZ[f] = new TH2F(nome, titulo, 50, 0, 20, 50, 0, 20);

    // --- conjuntos G e H: rotulo pelo pT do FOTON
    if (f == NFAIXAS - 1) sprintf(rotulo, "pT^{#gamma} > %.0f GeV", faixa_lo[f]);
    else sprintf(rotulo, "%.0f < pT^{#gamma} < %.0f GeV",
                 faixa_lo[f], faixa_hi[f]);

    for (Int_t o = 0; o < NORIGENS; o++)
    {
      // (G) isolados
      sprintf(nome, "hist_z_iso_%s_%d", curto_origem[o], f);
      sprintf(titulo, "#gamma ISOLADO de %s, %s;"
          "z_{#gamma} = pT^{#gamma}/pT^{jato};#gamma", nome_origem[o], rotulo);
      hist_z_iso[o][f] = new TH1F(nome, titulo, 50, 0, 1);

      // (H) nao isolados
      sprintf(nome, "hist_z_naoiso_%s_%d", curto_origem[o], f);
      sprintf(titulo, "#gamma NAO isolado de %s, %s;"
          "z_{#gamma} = pT^{#gamma}/pT^{jato};#gamma", nome_origem[o], rotulo);
      hist_z_naoiso[o][f] = new TH1F(nome, titulo, 50, 0, 1);
    }
  }

  // ---- contadores
  Long64_t n_jatos = 0, n_fotons = 0, n_gama = 0;
  Long64_t n_td[NFAIXAS] = {0, 0, 0, 0, 0, 0};
  Long64_t n_Z [NFAIXAS] = {0, 0, 0, 0, 0, 0};
  Long64_t n_AZ[NFAIXAS] = {0, 0, 0, 0, 0, 0};

  // fotons energeticos: por origem e por faixa de pT do foton
  Long64_t n_en    [NORIGENS][NFAIXAS];
  Long64_t n_en_iso[NORIGENS][NFAIXAS];
  Long64_t n_iso_zalto[NORIGENS];      // isolados com z > corte_z
  for (Int_t o = 0; o < NORIGENS; o++)
  {
    n_iso_zalto[o] = 0;
    for (Int_t f = 0; f < NFAIXAS; f++) { n_en[o][f] = 0; n_en_iso[o][f] = 0; }
  }

  // marcacao por massa: isolados que tem um parceiro em alguma das janelas.
  // Separados por janela, para dar para ver qual das duas esta' marcando o
  // que -- em especial quantos fotons de FRAGMENTACAO cada uma marca por
  // engano, que e' a medida do preco da marcacao.
  Long64_t n_marcado    [NORIGENS] = {0, 0, 0};   // pi0 OU eta
  Long64_t n_marcado_pi0[NORIGENS] = {0, 0, 0};   // so' a janela do pi0
  Long64_t n_marcado_eta[NORIGENS] = {0, 0, 0};   // so' a janela do eta
  Long64_t n_iso_tot    [NORIGENS] = {0, 0, 0};   // isolados, para a fracao

  Long64_t n_sem_jato   = 0;   // foton energetico que nao caiu em jato nenhum
  Long64_t n_sem_origem = 0;   // nem pi0, nem hadron, nem quark na ancestralidade
  Long64_t n_diag       = 0;   // quantas linhas de diagnostico ja' sairam

  // ---- laco de eventos
  for (Int_t iev = 0; iev < nev; iev++)
  {
    if (iev % 2000 == 0) printf("  evento %d / %d\n", iev, nev);

    pythia8->GenerateEvent();
    pythia8->ImportParticles(particulas, "All");
    Int_t np = particulas->GetEntriesFast();

    // monta a lista de entrada do FastJet
    entradas.clear();
    for (Int_t ip = 0; ip < np; ip++)
    {
      TParticle *part = (TParticle*) particulas->At(ip);

      if (part->GetStatusCode() <= 0) continue;
      if (eh_neutrino(part->GetPdgCode())) continue;
      if (part->Pt() < 1e-9) continue;
      if (TMath::Abs(part->Eta()) > eta_part_max) continue;

      fastjet::PseudoJet pj(part->Px(), part->Py(), part->Pz(), part->Energy());
      pj.set_user_index(ip);     // ponte de volta para o TClonesArray
      entradas.push_back(pj);
    }
    if (entradas.size() < 2) continue;

    // agrupa: sempre com pT minimo, senao cada particula mole vira um "jato"
    fastjet::ClusterSequence cluster_seq(entradas, jet_def);
    std::vector<fastjet::PseudoJet> jatos =
        fastjet::sorted_by_pt(cluster_seq.inclusive_jets(pt_jato_min));

    Int_t n_jatos_ev = jatos.size();

    // ================================================================
    //  PARTE 1 -- parte do JATO: histogramas gerais e conjuntos C, D, E
    // ================================================================
    for (Int_t ij = 0; ij < n_jatos_ev; ij++)
    {
      if (TMath::Abs(jatos[ij].rap()) > y_jato_max) continue;

      Double_t pt_jato = jatos[ij].pt();
      n_jatos++;

      std::vector<fastjet::PseudoJet> constituintes = jatos[ij].constituents();
      Int_t n_const = constituintes.size();

      hist_pt_jato->Fill(pt_jato);
      hist_n_part ->Fill(n_const);

      // faixa pelo pT do JATO (vale para os conjuntos C, D e E)
      Int_t faixa_jato = acha_faixa(pt_jato, faixa_lo, faixa_hi);

      // pi0 de onde vieram os fotons deste jato. O set so' serve para nao
      // repetir o pi0 na assimetria quando os DOIS fotons dele estao no jato.
      std::set<Int_t> indices_pi0;

      for (Int_t ic = 0; ic < n_const; ic++)
      {
        Int_t indice = constituintes[ic].user_index();
        TParticle *foton = (TParticle*) particulas->At(indice);

        if (foton->GetPdgCode() != 22) continue;   // so' fotons
        n_fotons++;

        // TParticle::Phi() devolve [0, 2pi); o eixo do mapa vai de -pi a +pi
        Double_t phi = foton->Phi();
        if (phi > TMath::Pi()) phi -= 2 * TMath::Pi();

        hist_pt_foton->Fill(foton->Pt());
        hist_eta_phi ->Fill(foton->Eta(), phi);

        // guarda o pi0 de origem, se for o caso
        Int_t mae = foton->GetFirstMother();
        if (mae >= 0 && mae < np)
          if (((TParticle*) particulas->At(mae))->GetPdgCode() == 111)
            indices_pi0.insert(mae);
      }

      // ---- assimetria e lider x nao-lider, um por pi0 (sem repetir)
      std::set<Int_t>::iterator it;
      for (it = indices_pi0.begin(); it != indices_pi0.end(); it++)
      {
        TParticle *part = (TParticle*) particulas->At(*it);

        Double_t z_pi0 = part->Pt() / pt_jato;   // usado nos cortes D e E

        // A mae direta do foton e' a copia do pi0 que decaiu,
        // entao as filhas dela sao os dois fotons.
        Int_t d1 = part->GetFirstDaughter();
        Int_t d2 = part->GetLastDaughter();

        if (d1 < 0 || d2 <= d1 || d2 >= np) continue;

        TParticle *gama1 = (TParticle*) particulas->At(d1);
        TParticle *gama2 = (TParticle*) particulas->At(d2);

        if (gama1->GetPdgCode() != 22) continue;
        if (gama2->GetPdgCode() != 22) continue;

        Double_t e1 = gama1->Energy();
        Double_t e2 = gama2->Energy();

        Double_t assimetria = TMath::Abs(e1 - e2) / (e1 + e2);
        hist_assimetria->Fill(assimetria);

        // separa o foton lider (mais energetico) do nao-lider
        TParticle *gama_lider     = gama1;
        TParticle *gama_nao_lider = gama2;
        if (e2 > e1) { gama_lider = gama2; gama_nao_lider = gama1; }

        hist_z_gama->Fill(gama_lider->Pt() / pt_jato);
        n_gama++;

        if (faixa_jato >= 0)
        {
          Double_t pt_lider     = gama_lider->Pt();
          Double_t pt_nao_lider = gama_nao_lider->Pt();

          // (C) sem corte
          hist_td[faixa_jato]->Fill(pt_lider, pt_nao_lider);
          n_td[faixa_jato]++;

          // (D) so z
          if (z_pi0 > corte_z) {
            hist_Z[faixa_jato]->Fill(pt_lider, pt_nao_lider);
            n_Z[faixa_jato]++;
          }

          // (E) z e assimetria
          if (z_pi0 > corte_z && assimetria > corte_A) {
            hist_AZ[faixa_jato]->Fill(pt_lider, pt_nao_lider);
            n_AZ[faixa_jato]++;
          }
        }
      }
    }

    // ================================================================
    //  PARTE 2 -- parte do FOTON: a analise nova
    //  Varre TODAS as particulas do evento, nao so' as que viraram
    //  constituintes de jato.
    // ================================================================
    for (Int_t ip = 0; ip < np; ip++)
    {
      TParticle *foton = (TParticle*) particulas->At(ip);

      if (foton->GetStatusCode() <= 0) continue;   // so' estado final
      if (foton->GetPdgCode() != 22)    continue;  // so' fotons

      Double_t pt_foton = foton->Pt();
      if (pt_foton < pt_foton_min) continue;       // so' os energeticos

      // O cone de isolamento precisa caber INTEIRO dentro da aceitacao.
      // Se o foton estiver na borda, a parte do cone que cai fora nao tem
      // particulas para somar e ele parece mais isolado do que e'.
      if (TMath::Abs(foton->Eta()) > eta_part_max - R_iso) continue;

      // faixa pelo pT do FOTON
      Int_t faixa = acha_faixa(pt_foton, faixa_lo, faixa_hi);
      if (faixa < 0) continue;

      // origem pelo PDG da mae direta
      Int_t origem = origem_foton(particulas, ip);

      // ---- DIAGNOSTICO: para os primeiros 20 fotons que AINDA ficam sem
      // origem, imprime o PDG da mae direta e, se a mae for um foton, tambem
      // o PDG da avo. Serve para conferir que a subida pelas copias resolveu.
      // PDG 0 na impressao quer dizer "sem mae valida na lista".
      // Para desligar o diagnostico, troque o 20 por 0.
      if (origem < 0 && n_diag < 20)
      {
        Int_t i_mae = foton->GetFirstMother();
        Int_t pdg_mae = (i_mae >= 0 && i_mae < np)
                      ? ((TParticle*) particulas->At(i_mae))->GetPdgCode() : 0;

        printf("  [diag %2lld] pT = %6.2f   mae PDG = %6d",
               n_diag + 1, pt_foton, pdg_mae);

        if (TMath::Abs(pdg_mae) == 22)
        {
          Int_t i_avo = ((TParticle*) particulas->At(i_mae))->GetFirstMother();
          Int_t pdg_avo = (i_avo >= 0 && i_avo < np)
                        ? ((TParticle*) particulas->At(i_avo))->GetPdgCode() : 0;
          printf("   avo PDG = %6d", pdg_avo);
        }
        printf("\n");
        n_diag++;
      }

      if (origem < 0) { n_sem_origem++; continue; }

      // isolamento
      Double_t soma_pt_cone = soma_cone(particulas, ip, R_iso, eta_part_max);
      Bool_t isolado = (soma_pt_cone < pt_iso_max);

      hist_soma_cone->Fill(soma_pt_cone);

      // ---- E_iso desta origem. Entra SEMPRE, isolado ou nao: o grafico
      // mostra a distribuicao inteira, e o corte de isolamento aparece nele
      // como uma linha vertical em 2 GeV/c.
      hist_eiso[origem]->Fill(soma_pt_cone);

      // ---- espectro de pT do foton ISOLADO, por origem.
      // Preenchido AQUI: logo depois de decidir o isolamento e ANTES de
      // procurar o jato associado. Nao depende de jato nenhum, entao conta
      // tambem os fotons isolados que nao cairam em nenhum jato
      // reconstruido (os que mais abaixo entram em n_sem_jato).
      if (isolado)
      {
        if      (origem == 0) hist_pt_iso_pi0eta->Fill(pt_foton);
        else if (origem == 1) hist_pt_iso_outro ->Fill(pt_foton);
        else if (origem == 2) hist_pt_iso_frag  ->Fill(pt_foton);
      }

      // ---- MASSA INVARIANTE: procura o parceiro.
      // Para cada foton energetico, combina com TODO outro foton final do
      // evento que um calorimetro conseguiria ver (dentro da aceitacao e
      // acima do limiar) e calcula a massa do par. Se algum par cair na
      // massa do pi0 OU na do eta, este foton e' filho de um desses dois
      // e deveria ser descartado.
      //
      // Guarda o par mais proximo de CADA uma das duas massas, que e' o que
      // um algoritmo de marcacao de verdade usaria para decidir.
      Double_t massa_pi0_melhor = -1,  dist_pi0_melhor = 1e9;
      Double_t massa_eta_melhor = -1,  dist_eta_melhor = 1e9;

      for (Int_t jp = 0; jp < np; jp++)
      {
        if (jp == ip) continue;                        // o proprio foton

        TParticle *parceiro = (TParticle*) particulas->At(jp);

        if (parceiro->GetStatusCode() <= 0) continue;  // so' estado final
        if (parceiro->GetPdgCode() != 22)   continue;  // so' fotons
        if (parceiro->Pt() < pt_parceiro_min) continue;
        if (TMath::Abs(parceiro->Eta()) > eta_part_max) continue;

        Double_t m = massa_par(foton, parceiro);

        hist_massa_gg->Fill(m);

        Double_t d_pi0 = TMath::Abs(m - M_PI0);
        if (d_pi0 < dist_pi0_melhor)
        { dist_pi0_melhor = d_pi0; massa_pi0_melhor = m; }

        Double_t d_eta = TMath::Abs(m - M_ETA);
        if (d_eta < dist_eta_melhor)
        { dist_eta_melhor = d_eta; massa_eta_melhor = m; }
      }

      // So' os ISOLADOS entram na estatistica de marcacao: sao eles que
      // sobrevivem ao isolamento e contaminam a amostra de sinal.
      if (isolado)
      {
        n_iso_tot[origem]++;

        // MARCADO: existe um parceiro visivel E o par cai na janela do pi0
        // ou na do eta. Sem parceiro visivel o foton NAO e' marcado -- e e'
        // justamente esse o caso irredutivel: o irmao ficou abaixo do
        // limiar e nao ha' como saber que ele existiu.
        Bool_t marcado_pi0 = (massa_pi0_melhor >= 0 &&
                              dist_pi0_melhor < janela_pi0);
        Bool_t marcado_eta = (massa_eta_melhor >= 0 &&
                              dist_eta_melhor < janela_eta);
        Bool_t marcado     = (marcado_pi0 || marcado_eta);

        if (massa_pi0_melhor >= 0) hist_massa_melhor    [origem]->Fill(massa_pi0_melhor);
        if (massa_eta_melhor >= 0) hist_massa_melhor_eta[origem]->Fill(massa_eta_melhor);

        // Os tres contadores: um por janela e um do total. O total NAO e' a
        // soma dos dois, porque um mesmo foton pode ter um parceiro na
        // janela do pi0 e outro, diferente, na janela do eta.
        if (marcado_pi0) n_marcado_pi0[origem]++;
        if (marcado_eta) n_marcado_eta[origem]++;
        if (marcado)     n_marcado    [origem]++;

        // ---- AQUI o foton marcado e' de fato REMOVIDO.
        // Ate' esta versao o codigo so' CONTAVA os marcados; o espectro
        // seguia com eles dentro. Agora ha' dois espectros: hist_pt_iso_*
        // (antes) e hist_pt_pos_* (depois), e e' o segundo que entra na
        // pureza pos-marcacao.
        if (!marcado)
        {
          if      (origem == 0) hist_pt_pos_pi0eta->Fill(pt_foton);
          else if (origem == 1) hist_pt_pos_outro ->Fill(pt_foton);
          else if (origem == 2) hist_pt_pos_frag  ->Fill(pt_foton);
        }
      }

      // ---- jato associado: aquele que CONTEM este foton.
      // Sem o corte de rapidez do jato: o corte de eta ja' foi aplicado
      // no proprio foton, logo acima.
      Double_t pt_jato_assoc = -1;
      for (Int_t ij = 0; ij < n_jatos_ev; ij++)
      {
        std::vector<fastjet::PseudoJet> cons = jatos[ij].constituents();
        Int_t n_cons = cons.size();

        for (Int_t ic = 0; ic < n_cons; ic++)
          if (cons[ic].user_index() == ip) { pt_jato_assoc = jatos[ij].pt(); break; }

        if (pt_jato_assoc > 0) break;
      }
      if (pt_jato_assoc < 0) { n_sem_jato++; continue; }

      Double_t z = pt_foton / pt_jato_assoc;

      n_en[origem][faixa]++;

      if (isolado)
      {
        hist_z_iso[origem][faixa]->Fill(z);
        n_en_iso[origem][faixa]++;
        if (z > corte_z) n_iso_zalto[origem]++;
      }
      else
      {
        hist_z_naoiso[origem][faixa]->Fill(z);
      }
    }
  }

  // ================================================================
  //  SECAO DE CHOQUE DESTA RODADA
  //  Lida só agora: sigmaGen() e' a estimativa Monte Carlo da integral e
  //  vai convergindo evento a evento, entao o valor que vale e' o final.
  //  O peso da rodada usa n_aceitos (e nao nev), porque sigmaGen e' definida
  //  em relacao aos eventos que o Pythia de fato aceitou.
  // ================================================================
  Double_t sigma     = pythia8->Pythia8()->info.sigmaGen();   // mb
  Double_t sigma_err = pythia8->Pythia8()->info.sigmaErr();   // mb
  Long64_t n_aceitos = (Long64_t) pythia8->Pythia8()->info.nAccepted();
  Double_t peso      = n_aceitos ? sigma / n_aceitos : 0.;

  hist_sigma->SetBinContent(1, sigma);
  hist_sigma->SetBinError  (1, sigma_err);
  hist_sigma->SetBinContent(2, sigma_err);
  hist_sigma->SetBinContent(3, (Double_t) n_aceitos);
  hist_sigma->SetBinContent(4, (Double_t) nev);
  hist_sigma->SetBinContent(5, peso);
  hist_sigma->SetBinContent(6, pt_hat_min);
  hist_sigma->SetBinContent(7, pt_hat_max);

  // ================================================================
  //  CANVASES
  // ================================================================

  // ---- c1: histogramas gerais
  c1 = new TCanvas("c1", "Jatos e fotons -- gerais", 1200, 1000);
  c1->Divide(2, 2);
  c1->cd(1); hist_pt_jato ->Draw();
  c1->cd(2); hist_n_part  ->Draw();
  c1->cd(3); hist_pt_foton->Draw();
  c1->cd(4); gPad->SetLogy(); hist_soma_cone->Draw();
  c1->Update();
  sprintf(png, "v2_%s_gerais.png", rodada);  c1->SaveAs(png);

  // ---- c2: mapa eta-phi
  c2 = new TCanvas("c2", "Mapa eta-phi dos fotons", 900, 700);
  hist_eta_phi->Draw("COLZ");
  c2->Update();
  sprintf(png, "v2_%s_mapa_eta_phi.png", rodada);  c2->SaveAs(png);

  // ---- c3: assimetria e foton lider
  c3 = new TCanvas("c3", "Fotons do decaimento do pi0", 1200, 500);
  c3->Divide(2, 1);
  c3->cd(1); hist_assimetria->Draw();
  c3->cd(2); gPad->SetLogy(); hist_z_gama->Draw();
  c3->Update();
  sprintf(png, "v2_%s_assimetria.png", rodada);  c3->SaveAs(png);

  // ---- c4, c5, c6: conjuntos C, D e E
  c4 = new TCanvas("c4", "Lider x nao-lider (todos)", 1500, 1000);
  c4->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c4->cd(f + 1); hist_td[f]->Draw("COLZ"); }
  c4->Update();
  sprintf(png, "v2_%s_lider_todos.png", rodada);  c4->SaveAs(png);

  c5 = new TCanvas("c5", "Lider x nao-lider (z alto)", 1500, 1000);
  c5->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c5->cd(f + 1); hist_Z[f]->Draw("COLZ"); }
  c5->Update();
  sprintf(png, "v2_%s_lider_z.png", rodada);  c5->SaveAs(png);

  c6 = new TCanvas("c6", "Lider x nao-lider (z e A altos)", 1500, 1000);
  c6->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c6->cd(f + 1); hist_AZ[f]->Draw("COLZ"); }
  c6->Update();
  sprintf(png, "v2_%s_lider_z_e_assimetria.png", rodada);  c6->SaveAs(png);

  // ---- c7, c8, c9: z dos ISOLADOS, uma origem por canvas
  c7 = new TCanvas("c7", "z: isolado, de pi0 + eta", 1500, 1000);
  c7->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c7->cd(f + 1); hist_z_iso[0][f]->Draw(); }
  c7->Update();
  sprintf(png, "v2_%s_z_iso_pi0eta.png", rodada);  c7->SaveAs(png);

  c8 = new TCanvas("c8", "z: isolado, de outro hadron", 1500, 1000);
  c8->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c8->cd(f + 1); hist_z_iso[1][f]->Draw(); }
  c8->Update();
  sprintf(png, "v2_%s_z_iso_outro.png", rodada);  c8->SaveAs(png);

  c9 = new TCanvas("c9", "z: isolado, de fragmentacao", 1500, 1000);
  c9->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c9->cd(f + 1); hist_z_iso[2][f]->Draw(); }
  c9->Update();
  sprintf(png, "v2_%s_z_iso_frag.png", rodada);  c9->SaveAs(png);

  // ---- c10, c11, c12: z dos NAO ISOLADOS, para comparar
  c10 = new TCanvas("c10", "z: nao isolado, de pi0 + eta", 1500, 1000);
  c10->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c10->cd(f + 1); hist_z_naoiso[0][f]->Draw(); }
  c10->Update();
  sprintf(png, "v2_%s_z_naoiso_pi0eta.png", rodada);  c10->SaveAs(png);

  c11 = new TCanvas("c11", "z: nao isolado, de outro hadron", 1500, 1000);
  c11->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c11->cd(f + 1); hist_z_naoiso[1][f]->Draw(); }
  c11->Update();
  sprintf(png, "v2_%s_z_naoiso_outro.png", rodada);  c11->SaveAs(png);

  c12 = new TCanvas("c12", "z: nao isolado, de fragmentacao", 1500, 1000);
  c12->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c12->cd(f + 1); hist_z_naoiso[2][f]->Draw(); }
  c12->Update();
  sprintf(png, "v2_%s_z_naoiso_frag.png", rodada);  c12->SaveAs(png);

  // ---- c13: MASSA INVARIANTE
  // Painel 1: todos os pares, alcance inteiro -- aparecem o pi0 e o eta.
  // Painel 2: zoom na janela do pi0.
  // Painel 3: zoom na janela do eta, que agora tambem marca.
  // Painel 4: o melhor par de cada foton ISOLADO, separado por origem.
  //           E' o painel que mostra se a marcacao funciona: os de
  //           pi0 + eta tem que empilhar em 0,135 e os de fragmentacao nao.
  c13 = new TCanvas("c13", "Massa invariante dos pares de fotons", 1400, 1000);
  c13->Divide(2, 2);

  c13->cd(1);
  gPad->SetLogy();
  gPad->SetLeftMargin(0.16);
  hist_massa_gg->SetLineColor(kBlack);
  hist_massa_gg->SetLineWidth(2);
  hist_massa_gg->Draw("HIST E");

  // Os zooms sao COPIAS: mexer no eixo do original mudaria o painel 1
  // e o que vai para o .root.
  c13->cd(2);
  gPad->SetLeftMargin(0.16);
  TH1F *zoom = (TH1F*) hist_massa_gg->Clone("hist_massa_gg_zoom");
  zoom->SetDirectory(0);
  zoom->SetTitle("zoom na janela do #pi^{0}");
  zoom->GetXaxis()->SetRangeUser(0., 0.3);
  zoom->SetLineColor(kBlue);
  zoom->SetLineWidth(2);
  zoom->Draw("HIST E");

  c13->cd(3);
  gPad->SetLeftMargin(0.16);
  TH1F *zoom_eta = (TH1F*) hist_massa_gg->Clone("hist_massa_gg_zoom_eta");
  zoom_eta->SetDirectory(0);
  zoom_eta->SetTitle("zoom na janela do #eta");
  zoom_eta->GetXaxis()->SetRangeUser(0.35, 0.75);
  zoom_eta->SetLineColor(kGreen + 2);
  zoom_eta->SetLineWidth(2);
  zoom_eta->Draw("HIST E");

  c13->cd(4);
  gPad->SetLogy();
  gPad->SetLeftMargin(0.16);
  {
    Int_t cor_m[NORIGENS] = {kBlue, kGreen + 2, kRed};
    Double_t topo_m = 0;
    for (Int_t o = 0; o < NORIGENS; o++)
      if (hist_massa_melhor[o]->GetMaximum() > topo_m)
        topo_m = hist_massa_melhor[o]->GetMaximum();

    TLegend *leg_m = new TLegend(0.45, 0.70, 0.89, 0.89);
    leg_m->SetTextSize(0.035);

    for (Int_t o = 0; o < NORIGENS; o++)
    {
      hist_massa_melhor[o]->SetTitle("melhor par, #gamma ISOLADO");
      hist_massa_melhor[o]->SetStats(0);
      hist_massa_melhor[o]->SetLineColor(cor_m[o]);
      hist_massa_melhor[o]->SetLineWidth(2);
      hist_massa_melhor[o]->SetMaximum(topo_m * 2.);
      hist_massa_melhor[o]->GetXaxis()->SetRangeUser(0., 0.7);
      hist_massa_melhor[o]->Draw(o == 0 ? "HIST E" : "HIST E SAME");

      sprintf(titulo, "%s  (%.0f)", curto_origem[o],
              hist_massa_melhor[o]->GetEntries());
      leg_m->AddEntry(hist_massa_melhor[o], titulo, "l");
    }
    leg_m->Draw();
  }
  c13->Update();
  sprintf(png, "v2_%s_massa_invariante.png", rodada);  c13->SaveAs(png);

  // ---- c14: E_iso (soma de pT no cone) por origem.
  // Cada curva e' dividida pela SUA PROPRIA area, entao o eixo nao e' mais
  // contagem: e' (1/N_gama) dN/dE_iso, e a area debaixo de cada uma vale 1.
  // Assim as tres formas podem ser comparadas mesmo tendo estatisticas
  // muito diferentes -- o que interessa aqui e' a FORMA, nao o tamanho.
  // A linha vertical e' o corte de isolamento: o que esta' a esquerda dela
  // e' aceito, o que esta' a direita e' jogado fora.
  c14 = new TCanvas("c14", "E_iso por origem", 900, 700);
  gPad->SetLogy();
  gPad->SetLeftMargin(0.16);
  {
    Int_t cor_e[NORIGENS] = {kBlue, kGreen + 2, kRed};
    TH1F *norm_e[NORIGENS];

    Double_t topo_e = 0, piso_e = 1e30;

    for (Int_t o = 0; o < NORIGENS; o++)
    {
      sprintf(nome, "hist_eiso_%s_norm", curto_origem[o]);
      norm_e[o] = (TH1F*) hist_eiso[o]->Clone(nome);
      norm_e[o]->SetDirectory(0);     // copia de grafico: nao vai para o .root

      // Scale(c, "width") divide cada bin pela largura dele e multiplica por
      // c. Com c = 1/integral, a area debaixo da curva passa a valer 1.
      Double_t area = norm_e[o]->Integral();
      if (area > 0) norm_e[o]->Scale(1. / area, "width");

      if (norm_e[o]->GetMaximum() > topo_e) topo_e = norm_e[o]->GetMaximum();

      for (Int_t b = 1; b <= norm_e[o]->GetNbinsX(); b++)
      {
        Double_t v = norm_e[o]->GetBinContent(b);
        if (v > 0 && v < piso_e) piso_e = v;
      }
    }
    if (piso_e > 1e29) piso_e = 1e-6;

    TLegend *leg_e = new TLegend(0.50, 0.70, 0.89, 0.89);
    leg_e->SetTextSize(0.035);

    for (Int_t o = 0; o < NORIGENS; o++)
    {
      norm_e[o]->SetTitle("isolamento por origem, normalizado por #gamma");
      norm_e[o]->SetYTitle("(1/N_{#gamma}) dN/dE_{iso}  [c/GeV]");
      norm_e[o]->SetStats(0);
      norm_e[o]->SetLineColor(cor_e[o]);
      norm_e[o]->SetLineWidth(2);
      norm_e[o]->SetMinimum(piso_e * 0.5);
      norm_e[o]->SetMaximum(topo_e * 3.);
      norm_e[o]->Draw(o == 0 ? "HIST E" : "HIST E SAME");
      leg_e->AddEntry(norm_e[o], nome_origem[o], "l");
    }

    TLine *linha_corte = new TLine(pt_iso_max, piso_e * 0.5,
                                   pt_iso_max, topo_e * 3.);
    linha_corte->SetLineColor(kBlack);
    linha_corte->SetLineStyle(2);
    linha_corte->SetLineWidth(2);
    linha_corte->Draw();

    leg_e->AddEntry(linha_corte, "corte de isolamento", "l");
    leg_e->Draw();
  }
  c14->Update();
  sprintf(png, "v2_%s_eiso.png", rodada);  c14->SaveAs(png);

  // ================================================================
  //  RESUMO -- vem ANTES do Write(), porque le os histogramas
  // ================================================================
  printf("\n===================== RESUMO =====================\n");
  printf("  rodada                    : %s\n", rodada);
  if (pt_hat_max > 0)
    printf("  faixa de pThat            : %.1f a %.1f GeV\n",
           pt_hat_min, pt_hat_max);
  else
    printf("  faixa de pThat            : acima de %.1f GeV (sem teto)\n",
           pt_hat_min);
  printf("  semente                   : %d\n", semente);
  printf("  eventos gerados           : %d\n", nev);
  printf("  eventos aceitos           : %lld\n", n_aceitos);
  printf("  secao de choque           : %.6g +- %.3g mb\n", sigma, sigma_err);
  printf("  peso = sigma / n_aceitos  : %.6g mb por evento\n", peso);
  printf("  jatos reconstruidos       : %lld\n", n_jatos);
  printf("  fotons em jatos           : %lld\n", n_fotons);
  printf("  <pT> dos jatos            : %.2f GeV/c\n", hist_pt_jato->GetMean());
  printf("  pares gama gama           : %lld\n", n_gama);
  printf("  <A> assimetria            : %.4f   (esperado ~0,5 se plana)\n",
         hist_assimetria->GetMean());
  printf("  RMS da assimetria         : %.4f   (esperado ~0,289 = 1/sqrt12)\n",
         hist_assimetria->GetRMS());

  printf("--------------------------------------------------\n");
  printf("  POR FAIXA DE pT DO JATO (lider x nao-lider dos pi0)\n");
  printf("  faixa [GeV]     todos  z_pi0>%.1f   z e A\n", corte_z);
  for (Int_t f = 0; f < NFAIXAS; f++)
  {
    if (f == NFAIXAS - 1) sprintf(rotulo, "> %.0f", faixa_lo[f]);
    else sprintf(rotulo, "%.0f - %.0f", faixa_lo[f], faixa_hi[f]);
    printf("  %-11s %9lld %9lld %8lld\n", rotulo, n_td[f], n_Z[f], n_AZ[f]);
  }

  // ---- (a) tabela por faixa de pT do foton
  printf("==================================================\n");
  printf("  FOTONS ENERGETICOS (pT > %.1f GeV/c)\n", pt_foton_min);
  printf("  (isolado: soma de pT no cone R < %.1f menor que %.1f GeV/c)\n",
         R_iso, pt_iso_max);
  printf("--------------------------------------------------\n");
  printf("  (a) CONTAGENS POR FAIXA DE pT DO FOTON\n");
  printf("  %-11s %22s %22s\n", "", "--- energeticos ---", "---- isolados ----");
  printf("  %-11s %7s %7s %7s %7s %7s %7s\n", "faixa [GeV]",
         "pi0", "outro", "frag", "pi0", "outro", "frag");

  Long64_t t_en[NORIGENS], t_iso[NORIGENS];
  for (Int_t o = 0; o < NORIGENS; o++) { t_en[o] = 0; t_iso[o] = 0; }

  for (Int_t f = 0; f < NFAIXAS; f++)
  {
    if (f == NFAIXAS - 1) sprintf(rotulo, "> %.0f", faixa_lo[f]);
    else sprintf(rotulo, "%.0f - %.0f", faixa_lo[f], faixa_hi[f]);

    printf("  %-11s %7lld %7lld %7lld %7lld %7lld %7lld\n", rotulo,
           n_en[0][f], n_en[1][f], n_en[2][f],
           n_en_iso[0][f], n_en_iso[1][f], n_en_iso[2][f]);

    for (Int_t o = 0; o < NORIGENS; o++)
    {
      t_en[o]  += n_en[o][f];
      t_iso[o] += n_en_iso[o][f];
    }
  }
  printf("  %-11s %7lld %7lld %7lld %7lld %7lld %7lld\n", "total",
         t_en[0], t_en[1], t_en[2], t_iso[0], t_iso[1], t_iso[2]);

  // ---- (b) fracao isolada de cada origem
  Long64_t t_en_tot = t_en[0] + t_en[1] + t_en[2];
  Long64_t t_iso_tot = t_iso[0] + t_iso[1] + t_iso[2];

  printf("--------------------------------------------------\n");
  printf("  (b) FRACAO ISOLADA = isolados / energeticos\n");
  for (Int_t o = 0; o < NORIGENS; o++)
    printf("    %-14s %7lld / %-7lld = %6.2f%%\n", curto_origem[o],
           t_iso[o], t_en[o], t_en[o] ? 100. * t_iso[o] / t_en[o] : 0.);
  printf("    %-14s %7lld / %-7lld = %6.2f%%\n", "total",
         t_iso_tot, t_en_tot, t_en_tot ? 100. * t_iso_tot / t_en_tot : 0.);

  // ---- (c) contaminacao: quanto da amostra vem de decaimento
  printf("--------------------------------------------------\n");
  printf("  (c) CONTAMINACAO = (pi0 + eta + outro hadron) / total\n");
  printf("      Compara a amostra sem corte com a amostra isolada.\n");
  printf("  %-11s %14s %14s %12s\n",
         "faixa [GeV]", "sem isolamento", "com isolamento", "reducao");
  for (Int_t f = 0; f < NFAIXAS; f++)
  {
    if (f == NFAIXAS - 1) sprintf(rotulo, "> %.0f", faixa_lo[f]);
    else sprintf(rotulo, "%.0f - %.0f", faixa_lo[f], faixa_hi[f]);

    Long64_t tot_f = n_en[0][f] + n_en[1][f] + n_en[2][f];
    Long64_t tot_i = n_en_iso[0][f] + n_en_iso[1][f] + n_en_iso[2][f];

    Double_t c_sem = tot_f ? 100. * (n_en[0][f] + n_en[1][f]) / tot_f : 0.;
    Double_t c_com = tot_i ? 100. * (n_en_iso[0][f] + n_en_iso[1][f]) / tot_i : 0.;

    printf("  %-11s %13.2f%% %13.2f%% %11.2f%%\n",
           rotulo, c_sem, c_com, c_sem - c_com);
  }

  Double_t c_sem_tot = t_en_tot  ? 100. * (t_en[0]  + t_en[1])  / t_en_tot  : 0.;
  Double_t c_com_tot = t_iso_tot ? 100. * (t_iso[0] + t_iso[1]) / t_iso_tot : 0.;
  printf("  %-11s %13.2f%% %13.2f%% %11.2f%%\n",
         "total", c_sem_tot, c_com_tot, c_sem_tot - c_com_tot);

  // ---- (d) isolados com z alto
  printf("--------------------------------------------------\n");
  printf("  (d) ISOLADOS COM z > %.1f, por origem\n", corte_z);
  for (Int_t o = 0; o < NORIGENS; o++)
    printf("    %-14s %7lld\n", curto_origem[o], n_iso_zalto[o]);
  printf("    %-14s %7lld\n", "total",
         n_iso_zalto[0] + n_iso_zalto[1] + n_iso_zalto[2]);

  // ---- (e) MARCACAO POR MASSA INVARIANTE
  printf("--------------------------------------------------\n");
  printf("  (e) MARCACAO POR MASSA INVARIANTE (pi0 e eta)\n");
  printf("      Um gama ISOLADO e' marcado se existir outro gama no\n");
  printf("      evento, com pT > %.2f GeV/c e |eta| < %.1f, tal que o par\n",
         pt_parceiro_min, eta_part_max);
  printf("      caia em UMA DAS DUAS janelas de massa:\n");
  printf("        pi0 : %.4f +- %.0f MeV\n", M_PI0, janela_pi0 * 1000.);
  printf("        eta : %.4f +- %.0f MeV\n", M_ETA, janela_eta * 1000.);
  printf("      A coluna frag e' o PRECO da marcacao: foton de fragmentacao\n");
  printf("      nao vem de decaimento nenhum, entao tudo que e' marcado la'\n");
  printf("      foi marcado por engano.\n");
  printf("  %-8s %9s %9s %9s %9s %9s\n", "origem", "isolados",
         "m. pi0", "m. eta", "m. total", "total %");

  for (Int_t o = 0; o < NORIGENS; o++)
    printf("  %-8s %9lld %9lld %9lld %9lld %8.2f%%\n", curto_origem[o],
           n_iso_tot[o], n_marcado_pi0[o], n_marcado_eta[o], n_marcado[o],
           n_iso_tot[o] ? 100. * n_marcado[o] / n_iso_tot[o] : 0.);
  printf("      (m. total nao e' a soma das duas janelas: um mesmo foton\n");
  printf("       pode ter um parceiro em cada janela.)\n");

  // O que sobra depois de jogar fora os marcados. E' a conta que interessa:
  // a marcacao so' vale a pena se tirar mais fundo do que sinal.
  {
    Long64_t sobra[NORIGENS];
    Long64_t antes = 0, depois = 0;
    for (Int_t o = 0; o < NORIGENS; o++)
    {
      sobra[o] = n_iso_tot[o] - n_marcado[o];
      antes  += n_iso_tot[o];
      depois += sobra[o];
    }

    printf("      ------------------------------------------\n");
    printf("      amostra isolada ANTES da marcacao : %lld\n", antes);
    printf("      amostra isolada DEPOIS            : %lld\n", depois);

    // "fragmentacao" e' indice 2; pi0 e outro hadron sao decaimento
    Double_t f_antes  = antes  ? 100. * n_iso_tot[2] / antes  : 0.;
    Double_t f_depois = depois ? 100. * sobra[2]     / depois : 0.;

    printf("      fracao de fragmentacao ANTES      : %6.2f%%\n", f_antes);
    printf("      fracao de fragmentacao DEPOIS     : %6.2f%%\n", f_depois);
    printf("      (nesta amostra nao ha' foton direto: HardQCD:all so'\n");
    printf("       produz decaimento e fragmentacao.)\n");
  }

  // ---- (f) E_ISO: quanto cada origem se concentra abaixo do corte.
  // A fracao abaixo do corte e' a EFICIENCIA do isolamento para aquela
  // origem. Como o corte de 2 GeV/c cai numa borda de bin, a integral de
  // 1 ate' o bin do corte conta exatamente o que e' aceito.
  printf("--------------------------------------------------\n");
  printf("  (f) E_iso = soma de pT no cone (R < %.1f), por origem\n", R_iso);
  printf("  %-8s %10s %12s %12s %10s\n", "origem", "energet.",
         "<E_iso>", "abaixo corte", "eficiencia");
  {
    Int_t b_corte = hist_eiso[0]->FindBin(pt_iso_max) - 1;  // ultimo aceito

    for (Int_t o = 0; o < NORIGENS; o++)
    {
      Int_t nbe = hist_eiso[o]->GetNbinsX();
      Double_t tudo  = hist_eiso[o]->Integral(0, nbe + 1);
      Double_t aceito = hist_eiso[o]->Integral(1, b_corte);

      printf("  %-8s %10.0f %12.3f %12.0f %9.2f%%\n", curto_origem[o],
             tudo, hist_eiso[o]->GetMean(), aceito,
             tudo ? 100. * aceito / tudo : 0.);
    }
  }

  // ---- (g) fotons que nao cairam em jato nenhum
  printf("--------------------------------------------------\n");
  printf("  (g) fotons energeticos sem jato : %lld\n", n_sem_jato);
  printf("      origem nao classificada     : %lld\n", n_sem_origem);
  printf("==================================================\n");

  // ---- salva tudo no .root (POR ULTIMO, depois de ler os histogramas).
  // Sem Close(): fechar o TFile destroi os histogramas que pertencem a ele,
  // e os canvases ficariam em branco. O Write() sozinho ja' grava em disco.
  arquivo->Write();
  printf("  Arquivo salvo: %s\n", nome_arquivo);
}
