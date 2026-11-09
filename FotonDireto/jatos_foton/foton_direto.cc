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
// Saida: foton_direto_<sufixo>.root  e  foton_direto_<sufixo>.png

#include "TSystem.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TPythia8.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TMath.h"
#include "TROOT.h"

// Os dois indices do TClonesArray que correspondem as saidas do processo
// duro. Vem das linhas 5 e 6 do registro do Pythia 8, menos 1 pela linha
// "system" que o ImportParticles descarta. Conferido pelo diagnostico.
const Int_t I_DURO_A = 4;
const Int_t I_DURO_B = 5;

// Ponteiros em escopo de arquivo: continuam acessiveis no prompt do ROOT
// depois que a macro termina.
TH1F *hist_pt_direto_iso = 0;   // bruto, em CONTAGEM
TH1F *hist_pt_direto_mb  = 0;   // a copia multiplicada pelo peso, em mb
TH1D *hist_sigma         = 0;   // informacao da rodada
TCanvas *c1 = 0;

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

  // quantos eventos imprimem o diagnostico. Troque por 0 para desligar.
  Int_t n_ev_diag = 3;

  // ---- buffers para montar textos
  char comando[200], nome_arquivo[200], png[200], rodada[60];

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

  // ---- contadores
  Long64_t n_energ      = 0;   // fotons finais energeticos na aceitacao
  Long64_t n_direto     = 0;   // desses, quantos sao diretos
  Long64_t n_direto_iso = 0;   // desses, quantos sao isolados

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
      if (soma_pt_cone >= pt_iso_max) continue;    // nao isolado

      n_direto_iso++;
      hist_pt_direto_iso->Fill(pt_foton);
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

  c1 = new TCanvas("c1", "Fotons diretos isolados", 900, 700);
  gPad->SetLogy();
  gPad->SetLeftMargin(0.16);           // senao os expoentes saem cortados
  hist_pt_direto_mb->Draw("HIST E");
  c1->Update();
  sprintf(png, "foton_direto_%s.png", rodada);
  c1->SaveAs(png);

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
  printf("  integral do grafico       : %.6g mb\n",
         hist_pt_direto_mb->Integral());
  printf("  (so' conta o que caiu entre 5 e 40 GeV/c)\n");
  printf("==================================================\n");

  // ---- salva (POR ULTIMO). Sem Close(): fechar o TFile destroi os
  // ---- histogramas que pertencem a ele e o canvas ficaria em branco.
  arquivo->Write();
  printf("  Arquivo salvo: %s  (histograma bruto + hist_sigma)\n", nome_arquivo);
  printf("  Figura salva : %s\n", png);
}
