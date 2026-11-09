// Abre o jatos_foton_isolado.root e redesenha TODOS os canvases do
// jatos_foton_v2.cc, sem rodar a simulacao de novo.
//
// E' o par do abre_jatos_foton.cc, mas para o arquivo do v2: os conjuntos
// A, B e F sairam e entraram G e H (z por origem, isolados e nao isolados),
// entao os nomes dos histogramas mudaram e o visualizador antigo nao serve.
//
// RODAR
//   root -l abre_jatos_foton_v2.cc
//
// Depois os ponteiros ficam disponiveis no prompt:
//   hist_z_iso[0][2]->Draw()        // pi0, faixa 10-12 GeV
//   hist_td[0]->Draw("COLZ")
//   c7->SaveAs("outro_nome.png")
//
// O arquivo NAO e' fechado de proposito: fechar um TFile destroi todos os
// histogramas que pertencem a ele, e os canvases ficariam em branco.

#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"

const Int_t NFAIXAS  = 6;
const Int_t NORIGENS = 3;   // 0 = pi0, 1 = outro hadron, 2 = fragmentacao

TFile *arquivo = 0;

TH1F *hist_pt_jato = 0, *hist_n_part = 0, *hist_pt_foton = 0;
TH1F *hist_assimetria = 0, *hist_z_gama = 0, *hist_soma_cone = 0;
TH2F *hist_eta_phi = 0;

TH2F *hist_td[NFAIXAS], *hist_Z[NFAIXAS], *hist_AZ[NFAIXAS];
TH1F *hist_z_iso[NORIGENS][NFAIXAS], *hist_z_naoiso[NORIGENS][NFAIXAS];

TCanvas *c1=0,*c2=0,*c3=0,*c4=0,*c5=0,*c6=0;
TCanvas *c7=0,*c8=0,*c9=0,*c10=0,*c11=0,*c12=0;

// desenha os seis paineis de um conjunto por faixa
void seis(TCanvas *c, TH1F **h)
{
  c->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c->cd(f + 1); if (h[f]) h[f]->Draw(); }
  c->Update();
}
void seis2(TCanvas *c, TH2F **h)
{
  c->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c->cd(f + 1); if (h[f]) h[f]->Draw("COLZ"); }
  c->Update();
}

void abre_jatos_foton_v2(const char *nome_arq = "jatos_foton_isolado.root")
{
  arquivo = TFile::Open(nome_arq);
  if (!arquivo || arquivo->IsZombie()) {
    printf("\n  ERRO: nao consegui abrir %s\n", nome_arq);
    printf("  Rode de dentro da pasta jatos_foton, ou passe o caminho:\n");
    printf("    root -l 'abre_jatos_foton_v2.cc(\"/caminho/arquivo.root\")'\n\n");
    return;
  }

  // ---- histogramas simples
  hist_pt_jato    = (TH1F*) arquivo->Get("hist_pt_jato");
  hist_n_part     = (TH1F*) arquivo->Get("hist_n_part");
  hist_pt_foton   = (TH1F*) arquivo->Get("hist_pt_foton");
  hist_assimetria = (TH1F*) arquivo->Get("hist_assimetria");
  hist_z_gama     = (TH1F*) arquivo->Get("hist_z_gama");
  hist_soma_cone  = (TH1F*) arquivo->Get("hist_soma_cone");
  hist_eta_phi    = (TH2F*) arquivo->Get("hist_eta_phi");

  // ---- os que vem em seis faixas
  const char *curto_origem[NORIGENS] = {"pi0", "outro", "frag"};
  char n[100];

  for (Int_t f = 0; f < NFAIXAS; f++)
  {
    sprintf(n, "hist_td_%d", f);  hist_td[f] = (TH2F*) arquivo->Get(n);
    sprintf(n, "hist_Z_%d",  f);  hist_Z [f] = (TH2F*) arquivo->Get(n);
    sprintf(n, "hist_AZ_%d", f);  hist_AZ[f] = (TH2F*) arquivo->Get(n);

    for (Int_t o = 0; o < NORIGENS; o++)
    {
      sprintf(n, "hist_z_iso_%s_%d", curto_origem[o], f);
      hist_z_iso[o][f] = (TH1F*) arquivo->Get(n);

      sprintf(n, "hist_z_naoiso_%s_%d", curto_origem[o], f);
      hist_z_naoiso[o][f] = (TH1F*) arquivo->Get(n);
    }
  }

  // ---- os 12 canvases, na mesma ordem e layout do jatos_foton_v2.cc
  c1 = new TCanvas("c1", "Jatos e fotons -- gerais", 1200, 1000);
  c1->Divide(2, 2);
  c1->cd(1); hist_pt_jato ->Draw();
  c1->cd(2); hist_n_part  ->Draw();
  c1->cd(3); hist_pt_foton->Draw();
  c1->cd(4); gPad->SetLogy(); hist_soma_cone->Draw();
  c1->Update();

  c2 = new TCanvas("c2", "Mapa eta-phi dos fotons", 900, 700);
  hist_eta_phi->Draw("COLZ");
  c2->Update();

  c3 = new TCanvas("c3", "Fotons do decaimento do pi0", 1200, 500);
  c3->Divide(2, 1);
  c3->cd(1); hist_assimetria->Draw();
  c3->cd(2); gPad->SetLogy(); hist_z_gama->Draw();
  c3->Update();

  c4 = new TCanvas("c4", "Lider x nao-lider (todos)",       1500, 1000);
  seis2(c4, hist_td);
  c5 = new TCanvas("c5", "Lider x nao-lider (z alto)",      1500, 1000);
  seis2(c5, hist_Z);
  c6 = new TCanvas("c6", "Lider x nao-lider (z e A altos)", 1500, 1000);
  seis2(c6, hist_AZ);

  c7 = new TCanvas("c7", "z: isolado, de pi0",           1500, 1000);
  seis(c7, hist_z_iso[0]);
  c8 = new TCanvas("c8", "z: isolado, de outro hadron",  1500, 1000);
  seis(c8, hist_z_iso[1]);
  c9 = new TCanvas("c9", "z: isolado, de fragmentacao",  1500, 1000);
  seis(c9, hist_z_iso[2]);

  c10 = new TCanvas("c10", "z: nao isolado, de pi0",          1500, 1000);
  seis(c10, hist_z_naoiso[0]);
  c11 = new TCanvas("c11", "z: nao isolado, de outro hadron", 1500, 1000);
  seis(c11, hist_z_naoiso[1]);
  c12 = new TCanvas("c12", "z: nao isolado, de fragmentacao", 1500, 1000);
  seis(c12, hist_z_naoiso[2]);

  // ---- recontagem a partir dos histogramas, para conferir a produção
  Long64_t t_en[NORIGENS], t_iso[NORIGENS];
  for (Int_t o = 0; o < NORIGENS; o++)
  {
    t_iso[o] = 0; t_en[o] = 0;
    for (Int_t f = 0; f < NFAIXAS; f++)
    {
      if (hist_z_iso[o][f])    t_iso[o] += (Long64_t) hist_z_iso[o][f]->GetEntries();
      if (hist_z_naoiso[o][f]) t_en [o] += (Long64_t) hist_z_naoiso[o][f]->GetEntries();
    }
    t_en[o] += t_iso[o];   // energeticos = isolados + nao isolados
  }
  Long64_t s_en  = t_en [0] + t_en [1] + t_en [2];
  Long64_t s_iso = t_iso[0] + t_iso[1] + t_iso[2];

  printf("\n============== ARQUIVO ABERTO ==============\n");
  printf("  %s\n", nome_arq);
  printf("  jatos no hist_pt_jato : %.0f\n", hist_pt_jato->GetEntries());
  printf("  <pT> dos jatos        : %.2f GeV/c\n", hist_pt_jato->GetMean());
  printf("  <A> assimetria        : %.4f\n", hist_assimetria->GetMean());
  printf("--------------------------------------------\n");
  printf("  FOTONS ENERGETICOS   energeticos   isolados\n");
  for (Int_t o = 0; o < NORIGENS; o++)
    printf("    %-14s %10lld %10lld\n", curto_origem[o], t_en[o], t_iso[o]);
  printf("    %-14s %10lld %10lld\n", "total", s_en, s_iso);
  if (s_en > 0 && s_iso > 0)
  {
    printf("  contaminacao sem isolamento : %6.2f%%\n",
           100. * (t_en[0]  + t_en[1])  / s_en);
    printf("  contaminacao com isolamento : %6.2f%%\n",
           100. * (t_iso[0] + t_iso[1]) / s_iso);
  }
  printf("--------------------------------------------\n");
  printf("  12 canvases redesenhados (c1 a c12)\n");
  printf("  No prompt voce pode usar:\n");
  printf("    hist_z_iso[0][2]->Draw()   pi0, faixa 10-12 GeV\n");
  printf("    hist_td[0]->Draw(\"COLZ\")\n");
  printf("    c7->SaveAs(\"nome.png\")     .q para sair\n");
  printf("============================================\n\n");
}
