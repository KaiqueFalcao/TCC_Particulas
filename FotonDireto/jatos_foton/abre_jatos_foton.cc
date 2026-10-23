// Abre o jatos_foton.root e redesenha TODOS os canvases, sem rodar a simulacao.
// Use depois de uma rodada longa para inspecionar os resultados quantas vezes
// quiser -- leva segundos em vez de uma hora.
//
// RODAR
//   root -l abre_jatos_foton.cc
//
// Depois os ponteiros ficam disponiveis no prompt:
//   hist_z->Draw()
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

const Int_t NFAIXAS = 6;

TFile  *arquivo = 0;

TH1F *hist_pt_jato = 0, *hist_n_part = 0, *hist_pt_foton = 0, *hist_z = 0;
TH1F *hist_assimetria = 0, *hist_z_gama = 0, *hist_soma_cone = 0;
TH2F *hist_eta_phi = 0, *hist_z_vs_ptjato = 0;

TH1F *hist_z_jato[NFAIXAS], *hist_z_foton[NFAIXAS];
TH2F *hist_td[NFAIXAS], *hist_Z[NFAIXAS], *hist_AZ[NFAIXAS];
TH1F *hist_z_iso_pi0[NFAIXAS], *hist_z_iso_nao_pi0[NFAIXAS];
TH1F *hist_z_naoiso_pi0[NFAIXAS], *hist_z_naoiso_nao_pi0[NFAIXAS];

TCanvas *c1=0,*c2=0,*c3=0,*c4=0,*c5=0,*c6=0,*c7=0;
TCanvas *c8=0,*c9=0,*c10=0,*c11=0,*c12=0,*c13=0,*c14=0;

// desenha os seis paineis de um conjunto de histogramas por faixa
void seis(TCanvas *c, TH1F **h, const char *opt)
{
  c->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c->cd(f + 1); if (h[f]) h[f]->Draw(opt); }
  c->Update();
}
void seis2(TCanvas *c, TH2F **h, const char *opt)
{
  c->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c->cd(f + 1); if (h[f]) h[f]->Draw(opt); }
  c->Update();
}

void abre_jatos_foton(const char *nome = "jatos_foton.root")
{
  arquivo = TFile::Open(nome);
  if (!arquivo || arquivo->IsZombie()) {
    printf("\n  ERRO: nao consegui abrir %s\n", nome);
    printf("  Rode de dentro da pasta jatos_foton, ou passe o caminho:\n");
    printf("    root -l 'abre_jatos_foton.cc(\"/caminho/jatos_foton.root\")'\n\n");
    return;
  }

  // ---- histogramas simples
  hist_pt_jato     = (TH1F*) arquivo->Get("hist_pt_jato");
  hist_n_part      = (TH1F*) arquivo->Get("hist_n_part");
  hist_pt_foton    = (TH1F*) arquivo->Get("hist_pt_foton");
  hist_z           = (TH1F*) arquivo->Get("hist_z");
  hist_assimetria  = (TH1F*) arquivo->Get("hist_assimetria");
  hist_z_gama      = (TH1F*) arquivo->Get("hist_z_gama");
  hist_soma_cone   = (TH1F*) arquivo->Get("hist_soma_cone");
  hist_eta_phi     = (TH2F*) arquivo->Get("hist_eta_phi");
  hist_z_vs_ptjato = (TH2F*) arquivo->Get("hist_z_vs_ptjato");

  // ---- os que vem em seis faixas
  char n[100];
  for (Int_t f = 0; f < NFAIXAS; f++)
  {
    sprintf(n, "hist_z_jato_%d", f);          hist_z_jato[f]  = (TH1F*) arquivo->Get(n);
    sprintf(n, "hist_z_foton_%d", f);         hist_z_foton[f] = (TH1F*) arquivo->Get(n);
    sprintf(n, "hist_td_%d", f);              hist_td[f]      = (TH2F*) arquivo->Get(n);
    sprintf(n, "hist_Z_%d", f);               hist_Z[f]       = (TH2F*) arquivo->Get(n);
    sprintf(n, "hist_AZ_%d", f);              hist_AZ[f]      = (TH2F*) arquivo->Get(n);
    sprintf(n, "hist_z_iso_pi0_%d", f);       hist_z_iso_pi0[f]       = (TH1F*) arquivo->Get(n);
    sprintf(n, "hist_z_iso_nao_pi0_%d", f);   hist_z_iso_nao_pi0[f]   = (TH1F*) arquivo->Get(n);
    sprintf(n, "hist_z_naoiso_pi0_%d", f);    hist_z_naoiso_pi0[f]    = (TH1F*) arquivo->Get(n);
    sprintf(n, "hist_z_naoiso_nao_pi0_%d", f);hist_z_naoiso_nao_pi0[f]= (TH1F*) arquivo->Get(n);
  }

  // ---- os 14 canvases, na mesma ordem e layout do jatos_foton.cc
  c1 = new TCanvas("c1", "Jatos e fragmentacao do foton", 1200, 1000);
  c1->Divide(2, 2);
  c1->cd(1); hist_pt_jato ->Draw();
  c1->cd(2); hist_n_part  ->Draw();
  c1->cd(3); hist_pt_foton->Draw();
  c1->cd(4); hist_z->Draw();
  c1->Update();

  c2 = new TCanvas("c2", "Mapa eta-phi dos fotons", 900, 700);
  hist_eta_phi->Draw("COLZ");
  c2->Update();

  c3 = new TCanvas("c3", "Fotons do decaimento do pi0", 1200, 500);
  c3->Divide(2, 1);
  c3->cd(1); hist_assimetria->Draw();
  c3->cd(2); gPad->SetLogy(); hist_z_gama->Draw();
  c3->Update();

  c4 = new TCanvas("c4", "z x pT do jato", 900, 700);
  gPad->SetLogz();
  hist_z_vs_ptjato->Draw("COLZ");
  c4->Update();

  c5  = new TCanvas("c5",  "z do foton por faixa de pT do jato",  1500, 1000);
  seis(c5, hist_z_jato, "");
  c6  = new TCanvas("c6",  "z do foton por faixa de pT do foton", 1500, 1000);
  seis(c6, hist_z_foton, "");

  c7  = new TCanvas("c7",  "Lider x nao-lider (todos)",       1500, 1000);
  seis2(c7, hist_td, "COLZ");
  c8  = new TCanvas("c8",  "Lider x nao-lider (z alto)",      1500, 1000);
  seis2(c8, hist_Z, "COLZ");
  c9  = new TCanvas("c9",  "Lider x nao-lider (z e A altos)", 1500, 1000);
  seis2(c9, hist_AZ, "COLZ");

  c10 = new TCanvas("c10", "Isolamento dos fotons", 900, 700);
  gPad->SetLogy();
  hist_soma_cone->Draw();
  c10->Update();

  c11 = new TCanvas("c11", "z: isolado e de pi0",         1500, 1000);
  seis(c11, hist_z_iso_pi0, "");
  c12 = new TCanvas("c12", "z: isolado e nao de pi0",     1500, 1000);
  seis(c12, hist_z_iso_nao_pi0, "");
  c13 = new TCanvas("c13", "z: nao isolado e de pi0",     1500, 1000);
  seis(c13, hist_z_naoiso_pi0, "");
  c14 = new TCanvas("c14", "z: nao isolado e nao de pi0", 1500, 1000);
  seis(c14, hist_z_naoiso_nao_pi0, "");

  // ---- resumo do que foi carregado
  printf("\n============== ARQUIVO ABERTO ==============\n");
  printf("  %s\n", nome);
  printf("  jatos no hist_pt_jato : %.0f\n", hist_pt_jato->GetEntries());
  printf("  fotons no hist_z      : %.0f\n", hist_z->GetEntries());
  printf("  <pT> dos jatos        : %.2f GeV/c\n", hist_pt_jato->GetMean());
  printf("  <A> assimetria        : %.4f\n", hist_assimetria->GetMean());
  printf("  14 canvases redesenhados (c1 a c14)\n");
  printf("============================================\n");
  printf("  No prompt voce pode usar:\n");
  printf("    hist_z->Draw()           hist_td[0]->Draw(\"COLZ\")\n");
  printf("    c7->SaveAs(\"nome.png\")   .q para sair\n\n");
}
