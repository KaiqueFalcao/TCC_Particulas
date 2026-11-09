// Soma as 4 faixas EXCLUSIVAS de pThat, cada uma pesada pela SUA PROPRIA
// secao de choque, e redesenha todos os graficos da jatos_foton_v2.cc.
//
// RODAR
//   root -l -b -q soma_pthat.cc
//
// DE ONDE VEM O PESO: cada rodada da jatos_foton_v2.cc guarda, no bin 5 do
// histograma hist_sigma, o valor sigma/N -- a secao de choque daquela faixa
// de pThat dividida pelo numero de eventos aceitos. E' quantos milibarn vale
// um evento daquela rodada. Multiplicando o histograma por esse peso, o eixo
// vertical deixa de ser contagem e passa a ser secao de choque.
//
// POR QUE A SOMA ESTA' CERTA AQUI: as 4 faixas sao exclusivas (4-7, 7-10,
// 10-15 e acima de 15), ou seja, nenhum evento pertence a duas delas. Entao
//      h_total = soma_i ( sigma_i / N_i ) * h_i
// reconstroi o espectro de pThat > 4 inteiro, sem contar nada duas vezes, e
// com a estatistica reforcada na parte de pThat alto.
// NAO existe media de secoes de choque nesta conta: cada faixa entra com a
// sua propria sigma.
//
// Os histogramas da jatos_foton_v2.cc tem Sumw2 ligado, entao o Scale e o Add
// propagam as incertezas corretamente.

#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TMath.h"
#include "TROOT.h"

const Int_t NRODADAS = 4;
const Int_t NFAIXAS  = 6;
const Int_t NORIGENS = 3;   // 0 = pi0, 1 = outro hadron, 2 = fragmentacao

// As 4 rodadas, na ordem das faixas de pThat
const char *sp_sufixo[NRODADAS] = {"pthat04_07", "pthat07_10",
                                   "pthat10_15", "pthat15_inf"};

TFile    *sp_arquivo[NRODADAS];
Double_t  sp_peso   [NRODADAS];

// Fator final aplicado a todos os somados, e o rotulo do eixo vertical.
// Em modo "ocorrencias" o fator converte mb de volta para contagem; em
// modo "mb" ele vale 1 e nada muda.
Double_t    sp_fator = 1.;
const char *sp_eixo  = "mb por bin";

// ---- histogramas SOMADOS. Mesmos nomes da jatos_foton_v2.cc, para a
// ---- correspondencia ficar obvia. Em escopo de arquivo, entao continuam
// ---- acessiveis no prompt do ROOT depois que a macro termina.
TH1F *hist_pt_jato    = 0;
TH1F *hist_n_part     = 0;
TH1F *hist_pt_foton   = 0;
TH2F *hist_eta_phi    = 0;
TH1F *hist_assimetria = 0;
TH1F *hist_z_gama     = 0;
TH1F *hist_soma_cone  = 0;

TH2F *hist_td[NFAIXAS] = {0, 0, 0, 0, 0, 0};
TH2F *hist_Z [NFAIXAS] = {0, 0, 0, 0, 0, 0};
TH2F *hist_AZ[NFAIXAS] = {0, 0, 0, 0, 0, 0};

TH1F *hist_z_iso   [NORIGENS][NFAIXAS];
TH1F *hist_z_naoiso[NORIGENS][NFAIXAS];

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

// Soma um histograma de UMA dimensao das 4 rodadas, cada uma multiplicada
// pelo peso da sua faixa. Devolve 0 se o histograma nao existir em algum
// dos arquivos.
//
// O clone recebe um nome temporario e so' depois e' rebatizado: se clonasse
// ja' com o nome final, ele nasceria dentro do TFile aberto, com nome
// repetido. O SetDirectory(0) solta o clone de qualquer arquivo, para que
// ele sobreviva ao Close la' no fim.
TH1F *soma_1d(const char *nome_hist)
{
  TH1F *total = 0;

  for (Int_t i = 0; i < NRODADAS; i++)
  {
    TH1F *h = (TH1F*) sp_arquivo[i]->Get(nome_hist);

    if (!h)
    {
      printf("  ERRO: %s nao existe em %s\n", nome_hist, sp_sufixo[i]);
      return 0;
    }

    if (i == 0)
    {
      total = (TH1F*) h->Clone("sp_tmp_total");
      total->SetDirectory(0);
      total->SetName(nome_hist);
      total->Scale(sp_peso[0]);
    }
    else
    {
      TH1F *parcela = (TH1F*) h->Clone("sp_tmp_parcela");
      parcela->SetDirectory(0);
      parcela->Scale(sp_peso[i]);
      total->Add(parcela);
      delete parcela;
    }
  }

  // Depois do peso o conteudo e' secao de choque, nao contagem. O rotulo
  // do eixo tem que dizer isso, senao o grafico mente.
  total->Scale(sp_fator);
  total->SetYTitle(sp_eixo);
  return total;
}

// Mesma coisa para os de DUAS dimensoes
TH2F *soma_2d(const char *nome_hist)
{
  TH2F *total = 0;

  for (Int_t i = 0; i < NRODADAS; i++)
  {
    TH2F *h = (TH2F*) sp_arquivo[i]->Get(nome_hist);

    if (!h)
    {
      printf("  ERRO: %s nao existe em %s\n", nome_hist, sp_sufixo[i]);
      return 0;
    }

    if (i == 0)
    {
      total = (TH2F*) h->Clone("sp_tmp_total2");
      total->SetDirectory(0);
      total->SetName(nome_hist);
      total->Scale(sp_peso[0]);
    }
    else
    {
      TH2F *parcela = (TH2F*) h->Clone("sp_tmp_parcela2");
      parcela->SetDirectory(0);
      parcela->Scale(sp_peso[i]);
      total->Add(parcela);
      delete parcela;
    }
  }

  total->Scale(sp_fator);
  total->SetZTitle(sp_eixo);
  return total;
}

// eventos_equiv > 0  : eixo em NUMERO DE OCORRENCIAS -- quantas haveria numa
//                      rodada inclusiva de pThat > 4 com esse tanto de
//                      eventos. E' o default, 100 mil, igual ao que cada
//                      faixa gerou.
// eventos_equiv = 0  : eixo em mb (secao de choque).
void soma_pthat(Double_t eventos_equiv = 100000.)
{
  const char *curto_origem[NORIGENS] = {"pi0", "outro", "frag"};

  char nome[200], png[200];

  // ================================================================
  //  1) ABRE OS 4 ARQUIVOS E LE OS PESOS
  // ================================================================
  printf("\n============== RODADAS E PESOS ==============\n");
  printf("  %-13s %12s %12s %16s\n",
         "rodada", "sigma[mb]", "N", "peso[mb/ev]");

  Double_t sigma_total = 0;

  for (Int_t i = 0; i < NRODADAS; i++)
  {
    sprintf(nome, "jatos_foton_%s.root", sp_sufixo[i]);
    sp_arquivo[i] = TFile::Open(nome);

    if (!sp_arquivo[i] || sp_arquivo[i]->IsZombie())
    {
      printf("  ERRO: nao consegui abrir %s\n", nome);
      printf("  Rode as 4 producoes antes de somar:\n");
      printf("    root -l -b -q 'jatos_foton_v2.cc(4., 7.,43,100000,\"pthat04_07\")'\n");
      printf("    root -l -b -q 'jatos_foton_v2.cc(7.,10.,44,100000,\"pthat07_10\")'\n");
      printf("    root -l -b -q 'jatos_foton_v2.cc(10.,15.,45,100000,\"pthat10_15\")'\n");
      printf("    root -l -b -q 'jatos_foton_v2.cc(15.,-1.,46,100000,\"pthat15_inf\")'\n");
      return;
    }

    TH1D *info = (TH1D*) sp_arquivo[i]->Get("hist_sigma");

    if (!info)
    {
      printf("  ERRO: %s nao tem hist_sigma.\n", nome);
      printf("  Foi gerado por uma versao antiga da jatos_foton_v2.cc?\n");
      return;
    }

    sp_peso[i] = info->GetBinContent(5);   // bin 5 = sigma / n_aceitos
    sigma_total += info->GetBinContent(1);

    printf("  %-13s %12.4f %12.0f %16.6e\n", sp_sufixo[i],
           info->GetBinContent(1), info->GetBinContent(3), sp_peso[i]);
  }

  printf("---------------------------------------------\n");
  printf("  sigma total das 4 faixas : %.4f mb\n", sigma_total);
  printf("  Cada faixa entra com a SUA sigma. Nao ha' media de sigma.\n");

  // ---- escolhe a unidade do eixo vertical
  // A soma pesada nasce em mb. Para virar NUMERO DE OCORRENCIAS basta
  // dividir pela sigma total (dando ocorrencias por mb) e multiplicar pelo
  // tamanho da amostra equivalente. Com eventos_equiv = 100000, o resultado
  // e' "quantas ocorrencias haveria em 100 mil eventos de pThat > 4" -- o
  // mesmo eixo de uma rodada unica, mas com a estatistica das 4 faixas.
  if (eventos_equiv > 0 && sigma_total > 0)
  {
    sp_fator = eventos_equiv / sigma_total;
    sp_eixo  = "ocorrencias";
    printf("  EIXO EM OCORRENCIAS, equivalente a %.0f eventos de pThat > 4\n",
           eventos_equiv);
    printf("  (fator aplicado: %.4f ocorrencias por mb)\n", sp_fator);
  }
  else
  {
    sp_fator = 1.;
    sp_eixo  = "mb por bin";
    printf("  EIXO EM mb (milibarn), ou seja, secao de choque.\n");
  }

  // ================================================================
  //  2) SOMA TODOS OS HISTOGRAMAS
  //  O hist_sigma fica de fora: ele guarda informacao da rodada,
  //  nao e' um espectro para somar.
  // ================================================================
  hist_pt_jato    = soma_1d("hist_pt_jato");
  hist_n_part     = soma_1d("hist_n_part");
  hist_pt_foton   = soma_1d("hist_pt_foton");
  hist_assimetria = soma_1d("hist_assimetria");
  hist_z_gama     = soma_1d("hist_z_gama");
  hist_soma_cone  = soma_1d("hist_soma_cone");

  hist_eta_phi    = soma_2d("hist_eta_phi");

  for (Int_t f = 0; f < NFAIXAS; f++)
  {
    sprintf(nome, "hist_td_%d", f);  hist_td[f] = soma_2d(nome);
    sprintf(nome, "hist_Z_%d",  f);  hist_Z [f] = soma_2d(nome);
    sprintf(nome, "hist_AZ_%d", f);  hist_AZ[f] = soma_2d(nome);

    for (Int_t o = 0; o < NORIGENS; o++)
    {
      sprintf(nome, "hist_z_iso_%s_%d", curto_origem[o], f);
      hist_z_iso[o][f] = soma_1d(nome);

      sprintf(nome, "hist_z_naoiso_%s_%d", curto_origem[o], f);
      hist_z_naoiso[o][f] = soma_1d(nome);
    }
  }

  // Confere que nada falhou antes de desenhar
  if (!hist_pt_jato || !hist_eta_phi || !hist_z_iso[2][5])
  {
    printf("  Algum histograma nao foi somado. Parando aqui.\n");
    return;
  }

  // ================================================================
  //  3) FECHA OS ARQUIVOS DE ENTRADA
  //  Pode fechar: todos os somados levaram SetDirectory(0) e estao
  //  na memoria, soltos de qualquer arquivo.
  // ================================================================
  for (Int_t i = 0; i < NRODADAS; i++) sp_arquivo[i]->Close();

  printf("  4 arquivos lidos e fechados. 61 histogramas somados.\n");

  // ================================================================
  //  4) ARQUIVO DE SAIDA
  //  Criado agora, DEPOIS da soma. Os somados sao entregues a ele
  //  com SetDirectory, para que o Write() os grave.
  // ================================================================
  TFile *arquivo = new TFile("soma_pthat.root", "RECREATE");

  hist_pt_jato   ->SetDirectory(arquivo);
  hist_n_part    ->SetDirectory(arquivo);
  hist_pt_foton  ->SetDirectory(arquivo);
  hist_assimetria->SetDirectory(arquivo);
  hist_z_gama    ->SetDirectory(arquivo);
  hist_soma_cone ->SetDirectory(arquivo);
  hist_eta_phi   ->SetDirectory(arquivo);

  for (Int_t f = 0; f < NFAIXAS; f++)
  {
    hist_td[f]->SetDirectory(arquivo);
    hist_Z [f]->SetDirectory(arquivo);
    hist_AZ[f]->SetDirectory(arquivo);

    for (Int_t o = 0; o < NORIGENS; o++)
    {
      hist_z_iso   [o][f]->SetDirectory(arquivo);
      hist_z_naoiso[o][f]->SetDirectory(arquivo);
    }
  }

  // ================================================================
  //  5) CANVASES -- os mesmos da jatos_foton_v2.cc, com prefixo soma_
  // ================================================================

  // ---- c1: histogramas gerais
  c1 = new TCanvas("c1", "SOMA -- jatos e fotons, gerais", 1200, 1000);
  c1->Divide(2, 2);
  c1->cd(1); hist_pt_jato ->Draw("HIST E");
  c1->cd(2); hist_n_part  ->Draw("HIST E");
  c1->cd(3); hist_pt_foton->Draw("HIST E");
  c1->cd(4); gPad->SetLogy(); hist_soma_cone->Draw("HIST E");
  c1->Update();
  c1->SaveAs("soma_gerais.png");

  // ---- c2: mapa eta-phi
  c2 = new TCanvas("c2", "SOMA -- mapa eta-phi dos fotons", 900, 700);
  hist_eta_phi->Draw("COLZ");
  c2->Update();
  c2->SaveAs("soma_mapa_eta_phi.png");

  // ---- c3: assimetria e foton lider
  c3 = new TCanvas("c3", "SOMA -- fotons do decaimento do pi0", 1200, 500);
  c3->Divide(2, 1);
  c3->cd(1); hist_assimetria->Draw("HIST E");
  c3->cd(2); gPad->SetLogy(); hist_z_gama->Draw("HIST E");
  c3->Update();
  c3->SaveAs("soma_assimetria.png");

  // ---- c4, c5, c6: lider x nao-lider, os 3 niveis de corte
  c4 = new TCanvas("c4", "SOMA -- lider x nao-lider (todos)", 1500, 1000);
  c4->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c4->cd(f + 1); hist_td[f]->Draw("COLZ"); }
  c4->Update();
  c4->SaveAs("soma_lider_todos.png");

  c5 = new TCanvas("c5", "SOMA -- lider x nao-lider (z alto)", 1500, 1000);
  c5->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c5->cd(f + 1); hist_Z[f]->Draw("COLZ"); }
  c5->Update();
  c5->SaveAs("soma_lider_z.png");

  c6 = new TCanvas("c6", "SOMA -- lider x nao-lider (z e A altos)", 1500, 1000);
  c6->Divide(3, 2);
  for (Int_t f = 0; f < NFAIXAS; f++) { c6->cd(f + 1); hist_AZ[f]->Draw("COLZ"); }
  c6->Update();
  c6->SaveAs("soma_lider_z_e_assimetria.png");

  // ---- c7, c8, c9: z dos ISOLADOS, uma origem por canvas
  TCanvas *canvas_iso[NORIGENS];
  const char *arq_iso[NORIGENS] = {"soma_z_iso_pi0.png",
                                   "soma_z_iso_outro.png",
                                   "soma_z_iso_frag.png"};
  const char *tit_iso[NORIGENS] = {"SOMA -- z: isolado, de pi0",
                                   "SOMA -- z: isolado, de outro hadron",
                                   "SOMA -- z: isolado, de fragmentacao"};

  for (Int_t o = 0; o < NORIGENS; o++)
  {
    sprintf(nome, "c%d", 7 + o);
    canvas_iso[o] = new TCanvas(nome, tit_iso[o], 1500, 1000);
    canvas_iso[o]->Divide(3, 2);
    for (Int_t f = 0; f < NFAIXAS; f++)
    {
      canvas_iso[o]->cd(f + 1);
      gPad->SetLeftMargin(0.16);          // senao os numeros saem cortados
      hist_z_iso[o][f]->Draw("HIST E");
    }
    canvas_iso[o]->Update();
    canvas_iso[o]->SaveAs(arq_iso[o]);
  }
  c7 = canvas_iso[0];  c8 = canvas_iso[1];  c9 = canvas_iso[2];

  // ---- c10, c11, c12: z dos NAO ISOLADOS, para comparar
  TCanvas *canvas_nao[NORIGENS];
  const char *arq_nao[NORIGENS] = {"soma_z_naoiso_pi0.png",
                                   "soma_z_naoiso_outro.png",
                                   "soma_z_naoiso_frag.png"};
  const char *tit_nao[NORIGENS] = {"SOMA -- z: nao isolado, de pi0",
                                   "SOMA -- z: nao isolado, de outro hadron",
                                   "SOMA -- z: nao isolado, de fragmentacao"};

  for (Int_t o = 0; o < NORIGENS; o++)
  {
    sprintf(nome, "c%d", 10 + o);
    canvas_nao[o] = new TCanvas(nome, tit_nao[o], 1500, 1000);
    canvas_nao[o]->Divide(3, 2);
    for (Int_t f = 0; f < NFAIXAS; f++)
    {
      canvas_nao[o]->cd(f + 1);
      gPad->SetLeftMargin(0.16);
      hist_z_naoiso[o][f]->Draw("HIST E");
    }
    canvas_nao[o]->Update();
    canvas_nao[o]->SaveAs(arq_nao[o]);
  }
  c10 = canvas_nao[0];  c11 = canvas_nao[1];  c12 = canvas_nao[2];

  // ================================================================
  //  RESUMO
  // ================================================================
  printf("=============================================\n");
  printf("  ESPECTRO SOMADO (integral em %s)\n", sp_eixo);
  printf("    pT dos jatos    : %.6g\n", hist_pt_jato->Integral());
  printf("    pT dos fotons   : %.6g\n", hist_pt_foton->Integral());
  printf("---------------------------------------------\n");
  printf("  FOTONS ENERGETICOS SOMADOS, por origem (%s)\n", sp_eixo);
  printf("  %-8s %14s %14s\n", "origem", "isolado", "nao isolado");

  for (Int_t o = 0; o < NORIGENS; o++)
  {
    Double_t s_iso = 0, s_nao = 0;
    for (Int_t f = 0; f < NFAIXAS; f++)
    {
      s_iso += hist_z_iso   [o][f]->Integral();
      s_nao += hist_z_naoiso[o][f]->Integral();
    }
    printf("  %-8s %14.4f %14.4f\n", curto_origem[o], s_iso, s_nao);
  }

  printf("=============================================\n");

  // Sem Close(): fechar o TFile de saida destroi os histogramas que agora
  // pertencem a ele, e os canvases ficariam em branco. O Write() grava.
  arquivo->Write();
  printf("  Arquivo salvo: soma_pthat.root\n");
  printf("  12 PNGs salvos com prefixo soma_\n");
}
