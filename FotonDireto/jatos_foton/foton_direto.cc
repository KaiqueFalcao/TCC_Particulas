// ============================================================================
//  foton_direto.cc  --  o SINAL: fotons DIRETOS (prompt photons)
// ============================================================================
//
// Codigo irmao do jatos_foton_v2.cc. Os dois medem o espectro de pT de fotons
// isolados, com os MESMOS cortes, mas de amostras diferentes:
//
//   jatos_foton_v2.cc   HardQCD:all        -> o RUIDO
//                       fotons de decaimento (pi0, outro hadron) e de
//                       fragmentacao, que imitam um foton direto
//
//   foton_direto.cc     PromptPhoton:all   -> o SINAL
//                       fotons que SAEM DO PROCESSO DURO, e' o que a gente
//                       quer de verdade medir
//
// Os dois graficos tem o MESMO BINNING (35 bins de 5 a 40 GeV/c) justamente
// para poderem ser sobrepostos depois, e os dois guardam hist_sigma com o
// mesmo layout de bins, para a normalizacao por secao de choque ser a mesma.
//
// Nao usa FastJet: o isolamento nao precisa de jato nenhum, ele soma o pT das
// particulas num cone em volta do foton. O jato so' aparecia no v2 para
// calcular z = pT_foton / pT_jato, e aqui nao medimos z.
//
// ----------------------------------------------------------------------------
//  A MARCACAO POR MASSA INVARIANTE TAMBEM E' APLICADA AQUI
// ----------------------------------------------------------------------------
// A marcacao existe para retirar foton de DECAIMENTO: se o foton energetico
// formar, com algum parceiro visivel, um par de massa igual a do pi0 ou a do
// eta, ele e' filho desse hadron e vai fora.
//
// Na amostra de sinal nao ha' decaimento para retirar -- a verdade da
// simulacao ja' garante que todo foton selecionado aqui saiu do processo
// duro. Entao tudo que a marcacao marcar neste codigo foi marcado POR ACASO:
// um foton direto que, junto com um foton qualquer do evento, caiu numa das
// janelas de massa sem ter nada a ver com ela.
//
// Rodar a marcacao nas duas amostras com os MESMOS criterios e' o que permite
// MEDIR esse acaso em vez de supor que ele e' pequeno. A fracao marcada aqui
// e' a perda de sinal do corte; a fracao marcada no ruido e' o ganho. A
// pureza pos-marcacao usa as duas.
//
// ----------------------------------------------------------------------------
//  COMO SE RECONHECE UM FOTON DIRETO
// ----------------------------------------------------------------------------
// No Pythia 6 as linhas 7 e 8 do registro eram as duas particulas que SAEM do
// processo duro. No Pythia 8 elas estao nas linhas 5 e 6, e o
// TPythia8::ImportParticles joga fora a linha 0 (a pseudo-particula "system",
// PDG 90) e desloca todos os indices em -1, mae inclusive. Entao no
// TClonesArray as duas saidas do processo duro cairiam nos indices 4 e 5.
//
// Isso esta' CONFERIDO pelo diagnostico la' embaixo (os 3 primeiros eventos
// imprimem o processo duro e a cadeia de maes de cada foton), e nao chutado.
//
// O PROBLEMA DA COPIA: o foton que chega ao estado final quase nunca e' a
// mesma linha do processo duro. O chuveiro RECOPIA a particula a cada emissao
// em que ela participa como recuo, e a copia nova tem o foton antigo (PDG 22)
// como mae. Por isso, para saber se um foton e' direto, subimos pelas maes
// ENQUANTO a mae for PDG 22: o foton mais antigo da cadeia e' o "topo". Se o
// topo for o indice 4 ou o 5, o foton nasceu no processo duro -- e' DIRETO.
//
// A ORDEM DO TESTE IMPORTA: a mae do foton direto tambem e' um quark ou um
// gluon (o parton que entrou na colisao). Se testassemos "a mae e' quark"
// primeiro, como faz a origem_foton do jatos_foton_v2.cc, o foton direto
// seria classificado como fragmentacao. Por isso o teste de "topo e' o
// processo duro" vem ANTES de qualquer outro.
//
// ----------------------------------------------------------------------------
//  RODAR
// ----------------------------------------------------------------------------
//   teste (5 mil eventos, so' para ver o diagnostico):
//     root -l -b -q 'foton_direto.cc(5., -1., 43, 5000, "teste")'
//
//   rodada de verdade (100 mil eventos):
//     root -l -b -q 'foton_direto.cc(5., -1., 43, 100000, "")'
//
//   nas mesmas faixas de pThat do ruido, se precisar:
//     root -l -b -q 'foton_direto.cc( 4.,  7., 43, 100000, "pthat04_07")'
//     root -l -b -q 'foton_direto.cc( 7., 10., 44, 100000, "pthat07_10")'
//     root -l -b -q 'foton_direto.cc(10., 15., 45, 100000, "pthat10_15")'
//     root -l -b -q 'foton_direto.cc(15., -1., 46, 100000, "pthat15_inf")'
//   pt_hat_max = -1 quer dizer "sem teto".
//
// Saida: foton_direto_<sufixo>.root  e  tres PNGs:
//   foton_direto_<sufixo>.png            espectro de pT, antes e depois
//   foton_direto_<sufixo>_massa.png      massa invariante dos pares
//   foton_direto_<sufixo>_eiso.png       soma de pT no cone (E_iso)

#include "TSystem.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TPythia8.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TROOT.h"

// Os dois indices do TClonesArray que correspondem as saidas do processo
// duro. Vem das linhas 5 e 6 do registro do Pythia 8, menos 1 pela linha
// "system" que o ImportParticles descarta. Conferido pelo diagnostico.
const Int_t I_DURO_A = 4;
const Int_t I_DURO_B = 5;

// Massas do pi0 e do eta, PDG 2024. As mesmas do jatos_foton_v2.cc: os dois
// codigos tem que marcar com o mesmo criterio, senao sinal e ruido estariam
// medidos com reguas diferentes.
const Double_t M_PI0 = 0.1349768;   // GeV
const Double_t M_ETA = 0.547862;    // GeV

// Ponteiros em escopo de arquivo: continuam acessiveis no prompt do ROOT
// depois que a macro termina.
TH1F *hist_pt_direto_iso = 0;   // bruto, em CONTAGEM
TH1F *hist_pt_direto_pos = 0;   // o mesmo, DEPOIS da marcacao por massa
TH1F *hist_pt_direto_mb  = 0;   // a copia multiplicada pelo peso, em mb
TH1F *hist_eiso_direto   = 0;   // soma de pT no cone, todo #gamma direto
TH1F *hist_massa_gg      = 0;   // todos os pares de fotons
TH1F *hist_massa_melhor     = 0;   // par mais proximo do pi0, #gamma isolado
TH1F *hist_massa_melhor_eta = 0;   // par mais proximo do eta, #gamma isolado
TH1D *hist_sigma         = 0;   // informacao da rodada
TCanvas *c1 = 0;
TCanvas *c2 = 0;   // massa invariante
TCanvas *c3 = 0;   // E_iso

Bool_t eh_neutrino(Int_t pdg) // devolve true se o PDG for de neutrino
{
  Int_t a = TMath::Abs(pdg);
  return (a == 12 || a == 14 || a == 16 || a == 18);
}

// Soma do pT das particulas finais num cone de raio R_iso em volta do foton
// de indice i_foton. IDENTICA a do jatos_foton_v2.cc -- tem que ser, senao
// sinal e ruido estariam medidos com criterios diferentes.
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

    // Phi() devolve [0, 2pi); a diferenca e' trazida para [-pi, pi]
    Double_t d_phi = p->Phi() - phi_foton;
    if (d_phi >  TMath::Pi()) d_phi -= 2 * TMath::Pi();
    if (d_phi < -TMath::Pi()) d_phi += 2 * TMath::Pi();

    Double_t d_eta = p->Eta() - eta_foton;

    Double_t R_i = TMath::Sqrt(d_eta * d_eta + d_phi * d_phi);
    if (R_i > R_iso) continue;

    soma_pt_cone += p->Pt();
  }
  return soma_pt_cone;
}

// Massa invariante de um par de fotons. IDENTICA a do jatos_foton_v2.cc.
//
//   m^2 = (E1 + E2)^2 - |p1 + p2|^2
//
// Para dois fotons, que nao tem massa, isso se reduz a
//   m^2 = 2 E1 E2 (1 - cos(theta))
// com theta o angulo entre eles. A conta e' feita somando os quadrivetores,
// que e' menos sujeita a erro do que a formula do angulo.
Double_t massa_par(TParticle *a, TParticle *b)
{
  Double_t E  = a->Energy() + b->Energy();
  Double_t px = a->Px()     + b->Px();
  Double_t py = a->Py()     + b->Py();
  Double_t pz = a->Pz()     + b->Pz();

  Double_t m2 = E*E - px*px - py*py - pz*pz;

  // m2 pode dar levemente negativo por arredondamento quando os dois fotons
  // sao quase colineares. Nesse caso a massa e' zero, nao NaN.
  return (m2 > 0) ? TMath::Sqrt(m2) : 0.;
}

// Sobe pelas COPIAS do proprio foton e devolve o indice do foton mais antigo
// da cadeia -- o "topo". Se o foton nao tiver sido recopiado, o topo e' ele
// mesmo. A trava de 50 passos protege contra qualquer laco no registro.
//
// Repare na diferenca em relacao a origem_foton do jatos_foton_v2.cc: la' nos
// queriamos o primeiro ancestral que NAO e' foton (para saber de que particula
// o foton saiu); aqui queremos o ultimo que AINDA e' foton (para saber se essa
// linha e' uma das saidas do processo duro).
Int_t topo_cadeia_foton(TClonesArray *lista, Int_t i_foton)
{
  Int_t n = lista->GetEntriesFast();
  Int_t atual = i_foton;
  Int_t passos = 0;

  while (passos < 50)
  {
    Int_t mae = ((TParticle*) lista->At(atual))->GetFirstMother();

    if (mae < 0 || mae >= n) break;    // nao tem mae valida: para aqui

    Int_t pdg_mae = TMath::Abs(((TParticle*) lista->At(mae))->GetPdgCode());
    if (pdg_mae != 22) break;          // a mae nao e' foton: o topo e' o atual

    atual = mae;                       // a mae e' foton: era uma copia, sobe
    passos++;
  }
  return atual;
}

// true se o foton for DIRETO, ou seja, se o topo da cadeia de copias for uma
// das duas linhas de saida do processo duro.
Bool_t eh_direto(TClonesArray *lista, Int_t i_foton)
{
  Int_t topo = topo_cadeia_foton(lista, i_foton);
  return (topo == I_DURO_A || topo == I_DURO_B);
}

// DIAGNOSTICO: imprime a cadeia de indices e PDGs de um foton ate' o topo do
// registro, para conferir a mao de onde ele veio.
void imprime_cadeia(TClonesArray *lista, Int_t i_foton)
{
  Int_t n = lista->GetEntriesFast();
  Int_t atual = i_foton;
  Int_t passos = 0;

  printf("        cadeia de maes:");

  while (passos < 50)
  {
    TParticle *p = (TParticle*) lista->At(atual);
    printf("  [%d]PDG=%d", atual, p->GetPdgCode());

    Int_t mae = p->GetFirstMother();
    if (mae < 0 || mae >= n) { printf("  (sem mae)"); break; }

    printf("  <-");
    atual = mae;
    passos++;
  }

  Int_t topo = topo_cadeia_foton(lista, i_foton);
  printf("\n        topo da cadeia de fotons = indice %d   =>  %s\n",
         topo, (topo == I_DURO_A || topo == I_DURO_B) ? "DIRETO" : "nao direto");
}

// pt_hat_min / pt_hat_max : faixa de pThat desta rodada (pt_hat_max = -1
//                           quer dizer sem teto)
// semente                 : Random:seed
// nev                     : eventos a gerar
// sufixo                  : entra no nome do .root e do PNG. Deixando vazio,
//                           e' montado a partir do pTHatMin.
void foton_direto(Double_t pt_hat_min = 5.,
                  Double_t pt_hat_max = -1.,
                  Int_t    semente    = 43,
                  Int_t    nev        = 100000,
                  const char *sufixo  = "")
{
  // ---- parametros: OS MESMOS do jatos_foton_v2.cc, nenhum mudou
  Double_t eta_part_max = 0.9;
  Double_t energia      = 14000.;
  Double_t R_iso        = 0.4;     // raio do cone de isolamento
  Double_t pt_iso_max   = 2.0;     // GeV/c: isolado se a soma no cone for menor
  Double_t pt_foton_min = 5.0;     // GeV/c: define "foton energetico"

  // ---- marcacao por massa invariante: OS MESMOS parametros do
  // ---- jatos_foton_v2.cc, para a regua ser a mesma nas duas amostras.
  Double_t pt_parceiro_min = 0.5;  // GeV/c: limiar para "ver" o parceiro
  Double_t janela_pi0      = 0.02; // GeV: meia-largura da janela do pi0
  Double_t janela_eta      = 0.02; // GeV: meia-largura da janela do eta

  // ---- E_iso: binning igual ao do jatos_foton_v2.cc. 0 a 20 GeV/c em 100
  // bins (0,2 GeV por bin), de modo que o corte de 2,0 GeV/c caia numa
  // BORDA de bin e a integral abaixo dele seja exatamente o que e' aceito.
  Double_t eiso_max   = 20.;
  Int_t    eiso_nbins = 100;

  // quantos eventos imprimem o diagnostico. Troque por 0 para desligar.
  Int_t n_ev_diag = 3;

  // ---- buffers para montar textos
  char comando[200], nome_arquivo[200], png[200], rodada[60], titulo[250];

  if (sufixo[0] != '\0') sprintf(rodada, "%s", sufixo);
  else                    sprintf(rodada, "pthat%02.0f", pt_hat_min);

  // ---- arquivo de saida: criado ANTES dos histogramas, para que eles
  // ---- sejam associados a ele e salvos no Write().
  sprintf(nome_arquivo, "foton_direto_%s.root", rodada);
  TFile *arquivo = new TFile(nome_arquivo, "RECREATE");

  // ---- Pythia: PromptPhoton no lugar de HardQCD. E' a unica diferenca de
  // ---- configuracao em relacao ao jatos_foton_v2.cc.
  TPythia8 *pythia8 = new TPythia8();
  pythia8->ReadString("PromptPhoton:all = on");
  pythia8->ReadString("Random:setSeed = on");

  sprintf(comando, "PhaseSpace:pTHatMin = %.1f", pt_hat_min);
  pythia8->ReadString(comando);

  if (pt_hat_max > 0)
  {
    sprintf(comando, "PhaseSpace:pTHatMax = %.1f", pt_hat_max);
    pythia8->ReadString(comando);
  }

  sprintf(comando, "Random:seed = %d", semente);
  pythia8->ReadString(comando);

  pythia8->Initialize(2212, 2212, energia);

  TClonesArray *particulas = new TClonesArray("TParticle", 10000);

  // ---- Sumw2 em todos os histogramas criados daqui para baixo, para que o
  // ---- Scale(peso) nao estrague as barras de erro.
  TH1::SetDefaultSumw2(kTRUE);

  // ---- informacao da rodada, com o MESMO layout de bins do jatos_foton_v2.cc
  hist_sigma = new TH1D("hist_sigma", "informacao da rodada;;valor", 7, 0, 7);
  hist_sigma->GetXaxis()->SetBinLabel(1, "sigma_mb");
  hist_sigma->GetXaxis()->SetBinLabel(2, "sigma_err_mb");
  hist_sigma->GetXaxis()->SetBinLabel(3, "n_aceitos");
  hist_sigma->GetXaxis()->SetBinLabel(4, "n_pedidos");
  hist_sigma->GetXaxis()->SetBinLabel(5, "peso_mb_por_ev");
  hist_sigma->GetXaxis()->SetBinLabel(6, "pt_hat_min");
  hist_sigma->GetXaxis()->SetBinLabel(7, "pt_hat_max");

  // ---- o histograma do sinal. 35 bins de 5 a 40 GeV/c, ou seja, 1 GeV por
  // ---- bin. EXATAMENTE o mesmo binning dos hist_pt_iso_* do
  // ---- jatos_foton_v2.cc, para os dois poderem ser sobrepostos.
  hist_pt_direto_iso = new TH1F("hist_pt_direto_iso",
      "#gamma DIRETO isolado;p_{T}^{#gamma} [GeV/c];#gamma",
      35, 5., 40.);

  // ---- o mesmo espectro DEPOIS da marcacao por massa invariante: so' os
  // que NAO foram marcados. A razao entre os dois e' a eficiencia da
  // marcacao para o sinal, que antes era suposta igual a 1 e agora e' medida.
  hist_pt_direto_pos = new TH1F("hist_pt_direto_pos",
      "#gamma DIRETO isolado e NAO marcado;p_{T}^{#gamma} [GeV/c];#gamma",
      35, 5., 40.);

  // ---- E_iso: soma de pT no cone, para TODO foton direto energetico
  // (isolado ou nao). Mesmo binning dos hist_eiso_* do jatos_foton_v2.cc,
  // para as quatro curvas irem ao mesmo grafico na pureza.cc.
  sprintf(titulo, "E_{iso} do #gamma direto;#Sigma pT^{cone} (R < %.1f) [GeV/c];#gamma", R_iso);
  hist_eiso_direto = new TH1F("hist_eiso_direto", titulo, eiso_nbins, 0., eiso_max);

  // ---- massa invariante dos pares. 500 bins de 0 a 1 GeV, 2 MeV por bin:
  // a largura foi escolhida para que 0,1349768 (pi0) e 0,547862 (eta) caiam
  // no CENTRO de um bin, e os picos aparecam no lugar certo.
  hist_massa_gg = new TH1F("hist_massa_gg",
      "massa invariante de pares de #gamma;m_{#gamma#gamma} [GeV/c^{2}];pares",
      500, 0., 1.);
  hist_massa_melhor = new TH1F("hist_massa_melhor",
      "melhor par (#pi^{0}), #gamma direto ISOLADO;m_{#gamma#gamma} [GeV/c^{2}];#gamma",
      500, 0., 1.);
  hist_massa_melhor_eta = new TH1F("hist_massa_melhor_eta",
      "melhor par (#eta), #gamma direto ISOLADO;m_{#gamma#gamma} [GeV/c^{2}];#gamma",
      500, 0., 1.);

  // ---- contadores
  Long64_t n_energ      = 0;   // fotons finais energeticos na aceitacao
  Long64_t n_direto     = 0;   // desses, quantos sao diretos
  Long64_t n_direto_iso = 0;   // desses, quantos sao isolados

  // marcacao por massa, nas duas janelas. Aqui tudo que e' marcado e' engano:
  // na amostra de sinal nao existe foton de decaimento.
  Long64_t n_marcado     = 0;  // pi0 OU eta
  Long64_t n_marcado_pi0 = 0;
  Long64_t n_marcado_eta = 0;
  Long64_t n_com_parceiro = 0; // isolados que tem algum parceiro visivel

  // ---- laco de eventos
  for (Int_t iev = 0; iev < nev; iev++)
  {
    if (iev % 2000 == 0) printf("  evento %d / %d\n", iev, nev);

    pythia8->GenerateEvent();
    pythia8->ImportParticles(particulas, "All");
    Int_t np = particulas->GetEntriesFast();

    // ================================================================
    //  DIAGNOSTICO dos primeiros eventos: confere na mao quais indices
    //  do TClonesArray sao o processo duro.
    // ================================================================
    if (iev < n_ev_diag)
    {
      printf("\n########## DIAGNOSTICO, evento %d ##########\n", iev);

      printf("\n-- process.list() do Pythia (o processo duro) --\n");
      pythia8->Pythia8()->process.list();

      // As 8 primeiras linhas do TClonesArray, lado a lado com a linha
      // correspondente do registro do Pythia. Se o deslocamento for de -1,
      // o indice i do TClonesArray tem que casar com event[i+1].
      printf("\n-- TClonesArray x event[] do Pythia --\n");
      printf("  %5s %8s %8s | %8s %8s %8s\n",
             "idx", "PDG", "final?", "ev[i+1]", "PDG", "status");

      Int_t n_mostra = (np < 8) ? np : 8;
      for (Int_t i = 0; i < n_mostra; i++)
      {
        TParticle *p = (TParticle*) particulas->At(i);
        printf("  %5d %8d %8d | %8d %8d %8d\n",
               i, p->GetPdgCode(), p->GetStatusCode(),
               i + 1,
               pythia8->Pythia8()->event[i + 1].id(),
               pythia8->Pythia8()->event[i + 1].status());
      }
      printf("  (status 23 = saida do subprocesso duro mais energetico)\n");

      // E agora a cadeia de cada foton final energetico
      printf("\n-- cadeia de maes dos fotons finais energeticos --\n");
      Int_t achou = 0;
      for (Int_t ip = 0; ip < np; ip++)
      {
        TParticle *f = (TParticle*) particulas->At(ip);
        if (f->GetStatusCode() <= 0)  continue;
        if (f->GetPdgCode()   != 22)  continue;
        if (f->Pt() < pt_foton_min)   continue;

        printf("      foton final [%d], pT = %.2f, eta = %+.3f\n",
               ip, f->Pt(), f->Eta());
        imprime_cadeia(particulas, ip);
        achou++;
      }
      if (!achou) printf("      (nenhum neste evento)\n");
      printf("############################################\n\n");
    }

    // ================================================================
    //  ANALISE
    // ================================================================
    for (Int_t ip = 0; ip < np; ip++)
    {
      TParticle *foton = (TParticle*) particulas->At(ip);

      if (foton->GetStatusCode() <= 0) continue;   // so' estado final
      if (foton->GetPdgCode() != 22)    continue;  // so' fotons

      Double_t pt_foton = foton->Pt();
      if (pt_foton < pt_foton_min) continue;       // so' os energeticos

      // O cone de isolamento precisa caber INTEIRO dentro da aceitacao, senao
      // a parte do cone que cai fora nao tem particulas para somar e o foton
      // parece mais isolado do que e'. Com eta_part_max = 0,9 e R_iso = 0,4
      // isso da' |eta| < 0,5.
      if (TMath::Abs(foton->Eta()) > eta_part_max - R_iso) continue;

      n_energ++;

      // ---- O TESTE DO FOTON DIRETO, antes de qualquer outro
      if (!eh_direto(particulas, ip)) continue;
      n_direto++;

      // ---- isolamento, com a mesma funcao e o mesmo corte do v2
      Double_t soma_pt_cone = soma_cone(particulas, ip, R_iso, eta_part_max);

      // E_iso entra SEMPRE, isolado ou nao: o grafico mostra a distribuicao
      // inteira, e o corte aparece nele como uma linha vertical.
      hist_eiso_direto->Fill(soma_pt_cone);

      if (soma_pt_cone >= pt_iso_max) continue;    // nao isolado

      n_direto_iso++;
      hist_pt_direto_iso->Fill(pt_foton);

      // ---- MARCACAO POR MASSA INVARIANTE, com os mesmos criterios do
      // jatos_foton_v2.cc: combina este foton com todo outro foton final
      // que um calorimetro veria, e guarda o par mais proximo de cada uma
      // das duas massas procuradas.
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

      Bool_t marcado_pi0 = (massa_pi0_melhor >= 0 &&
                            dist_pi0_melhor < janela_pi0);
      Bool_t marcado_eta = (massa_eta_melhor >= 0 &&
                            dist_eta_melhor < janela_eta);
      Bool_t marcado     = (marcado_pi0 || marcado_eta);

      if (massa_pi0_melhor >= 0)
      {
        n_com_parceiro++;
        hist_massa_melhor    ->Fill(massa_pi0_melhor);
        hist_massa_melhor_eta->Fill(massa_eta_melhor);
      }

      if (marcado_pi0) n_marcado_pi0++;
      if (marcado_eta) n_marcado_eta++;
      if (marcado)     n_marcado++;

      // ---- o marcado e' REMOVIDO, igual ao que o ruido faz
      if (!marcado) hist_pt_direto_pos->Fill(pt_foton);
    }
  }

  // ================================================================
  //  SECAO DE CHOQUE DESTA RODADA
  //  Lida so' agora: sigmaGen() vai convergindo evento a evento, entao o
  //  valor que vale e' o final. O peso usa n_aceitos (e nao nev) porque e'
  //  em relacao aos eventos aceitos que a sigmaGen esta' definida.
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
  //  GRAFICO: a copia multiplicada pelo peso, em secao de choque
  //  O SetDirectory(0) mantem essa copia FORA do .root -- o que vai para o
  //  arquivo e' o histograma bruto, em contagem, mais o hist_sigma. Quem
  //  quiser o grafico em mb refaz o Scale com o peso que esta' no bin 5.
  // ================================================================
  hist_pt_direto_mb = (TH1F*) hist_pt_direto_iso->Clone("hist_pt_direto_mb");
  hist_pt_direto_mb->SetDirectory(0);
  hist_pt_direto_mb->Scale(peso);
  hist_pt_direto_mb->SetTitle("#gamma DIRETO isolado -- secao de choque");
  hist_pt_direto_mb->SetXTitle("p_{T}^{#gamma} [GeV/c]");
  hist_pt_direto_mb->SetYTitle("#sigma [mb] por bin");
  hist_pt_direto_mb->SetLineColor(kBlue);
  hist_pt_direto_mb->SetLineWidth(2);
  hist_pt_direto_mb->SetStats(0);

  // A mesma copia escalada para o espectro DEPOIS da marcacao, para os dois
  // aparecerem no mesmo grafico e dar para ver quanto sinal o corte custou.
  TH1F *hist_pt_pos_mb = (TH1F*) hist_pt_direto_pos->Clone("hist_pt_pos_mb");
  hist_pt_pos_mb->SetDirectory(0);
  hist_pt_pos_mb->Scale(peso);
  hist_pt_pos_mb->SetLineColor(kRed);
  hist_pt_pos_mb->SetLineWidth(2);
  hist_pt_pos_mb->SetLineStyle(2);
  hist_pt_pos_mb->SetStats(0);

  c1 = new TCanvas("c1", "Fotons diretos isolados", 900, 700);
  gPad->SetLogy();
  gPad->SetLeftMargin(0.16);           // senao os expoentes saem cortados
  hist_pt_direto_mb->Draw("HIST E");
  hist_pt_pos_mb   ->Draw("HIST E SAME");

  TLegend *leg1 = new TLegend(0.45, 0.74, 0.89, 0.89);
  leg1->SetTextSize(0.032);
  leg1->AddEntry(hist_pt_direto_mb, "isolado", "l");
  leg1->AddEntry(hist_pt_pos_mb,    "isolado e nao marcado", "l");
  leg1->Draw();

  c1->Update();
  sprintf(png, "foton_direto_%s.png", rodada);
  c1->SaveAs(png);

  // ---- c2: MASSA INVARIANTE dos pares de fotons desta amostra.
  // Painel 1: todos os pares. Diferente do ruido, aqui NAO deve haver pico
  //           em 0,135 nem em 0,548 vindo do foton energetico -- se houver,
  //           e' de pi0 e eta do resto do evento, nao do foton selecionado.
  // Painel 2: o melhor par de cada foton direto isolado, nas duas janelas.
  //           E' o grafico que mede o acaso: as janelas estao marcadas.
  c2 = new TCanvas("c2", "Massa invariante -- sinal", 1200, 500);
  c2->Divide(2, 1);

  c2->cd(1);
  gPad->SetLogy();  gPad->SetLeftMargin(0.16);
  hist_massa_gg->SetLineColor(kBlack);
  hist_massa_gg->SetLineWidth(2);
  hist_massa_gg->Draw("HIST E");

  c2->cd(2);
  gPad->SetLogy();  gPad->SetLeftMargin(0.16);
  hist_massa_melhor->SetTitle("melhor par do #gamma direto isolado");
  hist_massa_melhor->SetStats(0);
  hist_massa_melhor->SetLineColor(kBlue);
  hist_massa_melhor->SetLineWidth(2);
  hist_massa_melhor->GetXaxis()->SetRangeUser(0., 0.7);
  hist_massa_melhor->Draw("HIST E");

  hist_massa_melhor_eta->SetStats(0);
  hist_massa_melhor_eta->SetLineColor(kGreen + 2);
  hist_massa_melhor_eta->SetLineWidth(2);
  hist_massa_melhor_eta->Draw("HIST E SAME");

  TLegend *leg2 = new TLegend(0.40, 0.74, 0.89, 0.89);
  leg2->SetTextSize(0.032);
  leg2->AddEntry(hist_massa_melhor,     "par mais proximo do #pi^{0}", "l");
  leg2->AddEntry(hist_massa_melhor_eta, "par mais proximo do #eta",    "l");
  leg2->Draw();

  c2->Update();
  sprintf(png, "foton_direto_%s_massa.png", rodada);
  c2->SaveAs(png);

  // ---- c3: E_iso do foton direto, normalizado por foton (area 1).
  // A linha vertical e' o corte de isolamento. Comparada com as curvas das
  // tres origens do ruido (ver pureza.cc), esta e' a mais empilhada perto de
  // zero: o foton direto nasce sozinho no evento.
  c3 = new TCanvas("c3", "E_iso -- sinal", 900, 700);
  gPad->SetLogy();  gPad->SetLeftMargin(0.16);
  {
    TH1F *norm = (TH1F*) hist_eiso_direto->Clone("hist_eiso_direto_norm");
    norm->SetDirectory(0);            // copia de grafico: fora do .root

    // Scale(c, "width") divide cada bin pela largura e multiplica por c.
    // Com c = 1/integral, a area debaixo da curva passa a valer 1.
    Double_t area = norm->Integral();
    if (area > 0) norm->Scale(1. / area, "width");

    Double_t piso = 1e30;
    for (Int_t b = 1; b <= norm->GetNbinsX(); b++)
    {
      Double_t v = norm->GetBinContent(b);
      if (v > 0 && v < piso) piso = v;
    }
    if (piso > 1e29) piso = 1e-6;

    norm->SetTitle("E_{iso} do #gamma direto, normalizado por #gamma");
    norm->SetYTitle("(1/N_{#gamma}) dN/dE_{iso}  [c/GeV]");
    norm->SetStats(0);
    norm->SetLineColor(kBlue);
    norm->SetLineWidth(2);
    norm->SetMinimum(piso * 0.5);
    norm->Draw("HIST E");

    TLine *linha = new TLine(pt_iso_max, piso * 0.5,
                             pt_iso_max, norm->GetMaximum());
    linha->SetLineColor(kBlack);
    linha->SetLineStyle(2);
    linha->SetLineWidth(2);
    linha->Draw();
  }
  c3->Update();
  sprintf(png, "foton_direto_%s_eiso.png", rodada);
  c3->SaveAs(png);

  // ================================================================
  //  RESUMO -- vem ANTES do Write(), porque le os histogramas
  // ================================================================
  printf("\n===================== RESUMO =====================\n");
  printf("  processo                  : PromptPhoton:all (SINAL)\n");
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
  printf("--------------------------------------------------\n");
  printf("  CORTES (os mesmos do jatos_foton_v2.cc)\n");
  printf("    pT do foton             : > %.1f GeV/c\n", pt_foton_min);
  printf("    aceitacao do foton      : |eta| < %.1f\n",
         eta_part_max - R_iso);
  printf("    isolado                 : soma de pT em R < %.1f abaixo de %.1f GeV/c\n",
         R_iso, pt_iso_max);
  printf("--------------------------------------------------\n");
  printf("  fotons energeticos        : %lld\n", n_energ);
  printf("  desses, DIRETOS           : %lld", n_direto);
  if (n_energ) printf("   (%.2f%%)", 100. * n_direto / n_energ);
  printf("\n");
  printf("  desses, ISOLADOS          : %lld", n_direto_iso);
  if (n_direto) printf("   (%.2f%%)", 100. * n_direto_iso / n_direto);
  printf("\n");
  printf("  fracao isolada            : %.2f%%\n",
         n_direto ? 100. * n_direto_iso / n_direto : 0.);
  printf("--------------------------------------------------\n");
  printf("  MARCACAO POR MASSA INVARIANTE (pi0 e eta)\n");
  printf("    parceiro                : pT > %.2f GeV/c, |eta| < %.1f\n",
         pt_parceiro_min, eta_part_max);
  printf("    janela do pi0           : %.4f +- %.0f MeV\n",
         M_PI0, janela_pi0 * 1000.);
  printf("    janela do eta           : %.4f +- %.0f MeV\n",
         M_ETA, janela_eta * 1000.);
  printf("    isolados com parceiro   : %lld", n_com_parceiro);
  if (n_direto_iso) printf("   (%.2f%%)", 100. * n_com_parceiro / n_direto_iso);
  printf("\n");
  printf("    marcados pela janela pi0: %lld\n", n_marcado_pi0);
  printf("    marcados pela janela eta: %lld\n", n_marcado_eta);
  printf("    marcados no total       : %lld", n_marcado);
  if (n_direto_iso) printf("   (%.2f%%)", 100. * n_marcado / n_direto_iso);
  printf("\n");
  printf("    sobram depois           : %lld\n", n_direto_iso - n_marcado);
  printf("    TUDO que foi marcado aqui foi marcado POR ACASO: nesta\n");
  printf("    amostra nao existe foton de decaimento. Essa fracao e' a\n");
  printf("    PERDA DE SINAL do corte de massa.\n");
  printf("--------------------------------------------------\n");
  printf("  integral do grafico       : %.6g mb\n",
         hist_pt_direto_mb->Integral());
  printf("  (so' conta o que caiu entre 5 e 40 GeV/c)\n");
  printf("==================================================\n");

  // ---- salva (POR ULTIMO). Sem Close(): fechar o TFile destroi os
  // ---- histogramas que pertencem a ele e o canvas ficaria em branco.
  arquivo->Write();
  printf("  Arquivo salvo: %s  (histogramas brutos + hist_sigma)\n", nome_arquivo);
  printf("  Figuras salvas: foton_direto_%s.png, _massa.png e _eiso.png\n",
         rodada);
}
