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
#include "TMath.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include <vector>
#include <set>

const Int_t NFAIXAS = 6;
const Int_t NORIGENS = 3;   // 0 = pi0, 1 = outro hadron, 2 = fragmentacao

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
TH1F   *hist_pt_iso_pi0   = 0;
TH1F   *hist_pt_iso_frag  = 0;
TH1F   *hist_pt_iso_outro = 0;

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

// Devolve em qual faixa um valor de pT cai (-1 se estiver fora de todas)
Int_t acha_faixa(Double_t pt, Double_t *lo, Double_t *hi)
{
  for (Int_t f = 0; f < NFAIXAS; f++)
    if (pt >= lo[f] && pt < hi[f]) return f;
  return -1;
}

// Classifica a origem do foton pelo PDG da MAE DIRETA.
// GetFirstMother() devolve um INDICE na lista, nao um PDG.
//    0 = decaimento de pi0           (|PDG| == 111)
//    1 = decaimento de outro hadron  (|PDG| > 100)
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
// A ordem dos testes importa: 111 e' maior que 100, entao o pi0 tem que ser
// testado ANTES do caso geral de hadron.
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

  if (pdg == 111)           return 0;   // pi0
  if (pdg >= 1 && pdg <= 6) return 2;   // quark -> fragmentacao
  if (pdg > 100)            return 1;   // outro hadron
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

  // seis faixas; a ultima vai ate 1000, ou seja, "acima de 20 GeV"
  Double_t faixa_lo[NFAIXAS] = { 5.,  8., 10., 12., 16.,   20.};
  Double_t faixa_hi[NFAIXAS] = { 8., 10., 12., 16., 20., 1000.};

  // nomes das tres origens, usados nos titulos e no resumo
  const char *nome_origem[NORIGENS] = {"#pi^{0}", "outro hadron", "fragmentacao"};
  const char *curto_origem[NORIGENS] = {"pi0", "outro", "frag"};

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
  hist_pt_iso_pi0   = new TH1F("hist_pt_iso_pi0",
      "#gamma isolado de #pi^{0};p_{T}^{#gamma} [GeV/c];#gamma",
      35, 5., 40.);
  hist_pt_iso_frag  = new TH1F("hist_pt_iso_frag",
      "#gamma isolado de fragmentacao;p_{T}^{#gamma} [GeV/c];#gamma",
      35, 5., 40.);
  hist_pt_iso_outro = new TH1F("hist_pt_iso_outro",
      "#gamma isolado de outro hadron;p_{T}^{#gamma} [GeV/c];#gamma",
      35, 5., 40.);

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

      // ---- espectro de pT do foton ISOLADO, por origem.
      // Preenchido AQUI: logo depois de decidir o isolamento e ANTES de
      // procurar o jato associado. Nao depende de jato nenhum, entao conta
      // tambem os fotons isolados que nao cairam em nenhum jato
      // reconstruido (os que mais abaixo entram em n_sem_jato).
      if (isolado)
      {
        if      (origem == 0) hist_pt_iso_pi0  ->Fill(pt_foton);
        else if (origem == 1) hist_pt_iso_outro->Fill(pt_foton);
        else if (origem == 2) hist_pt_iso_frag ->Fill(pt_foton);
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
  c7 = new TCanvas("c7", "z: isolado, de pi0", 1500, 1000);
  c7->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c7->cd(f + 1); hist_z_iso[0][f]->Draw(); }
  c7->Update();
  sprintf(png, "v2_%s_z_iso_pi0.png", rodada);  c7->SaveAs(png);

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
  c10 = new TCanvas("c10", "z: nao isolado, de pi0", 1500, 1000);
  c10->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c10->cd(f + 1); hist_z_naoiso[0][f]->Draw(); }
  c10->Update();
  sprintf(png, "v2_%s_z_naoiso_pi0.png", rodada);  c10->SaveAs(png);

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
  printf("  (c) CONTAMINACAO = (pi0 + outro hadron) / total\n");
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

  // ---- (e) fotons que nao cairam em jato nenhum
  printf("--------------------------------------------------\n");
  printf("  (e) fotons energeticos sem jato : %lld\n", n_sem_jato);
  printf("      origem nao classificada     : %lld\n", n_sem_origem);
  printf("==================================================\n");

  // ---- salva tudo no .root (POR ULTIMO, depois de ler os histogramas).
  // Sem Close(): fechar o TFile destroi os histogramas que pertencem a ele,
  // e os canvases ficariam em branco. O Write() sozinho ja' grava em disco.
  arquivo->Write();
  printf("  Arquivo salvo: %s\n", nome_arquivo);
}
