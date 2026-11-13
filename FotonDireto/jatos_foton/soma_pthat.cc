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
//
// UNIDADE DO EIXO VERTICAL: o padrao e' SECAO DE CHOQUE, em mb por bin. E' o
// que a soma pesada produz naturalmente, e e' a unidade em que o resultado se
// compara com dado experimental e com a foton_direto.cc. Quem quiser o eixo em
// numero de ocorrencias passa o numero de eventos equivalente:
//   root -l -b -q soma_pthat.cc                 // mb  (padrao)
//   root -l -b -q 'soma_pthat.cc(100000)'       // ocorrencias em 100 mil ev

#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TMath.h"
#include "TLegend.h"
#include "TROOT.h"

const Int_t NRODADAS = 4;
const Int_t NFAIXAS  = 6;
const Int_t NORIGENS = 3;   // 0 = pi0, 1 = outro hadron, 2 = fragmentacao

// Massa do pi0 e a janela usada para marcar, as mesmas da jatos_foton_v2.cc
const Double_t M_PI0      = 0.1349768;   // GeV
const Double_t JANELA_PI0 = 0.02;        // GeV

// As 4 rodadas, na ordem das faixas de pThat
const char *sp_sufixo[NRODADAS] = {"pthat04_07", "pthat07_10",
                                   "pthat10_15", "pthat15_inf"};

TFile    *sp_arquivo[NRODADAS];
Double_t  sp_peso   [NRODADAS];

// Fator final aplicado a todos os somados, e o rotulo do eixo vertical.
// Em modo "ocorrencias" o fator converte mb de volta para contagem; em
// modo "mb" ele vale 1 e nada muda.
Double_t    sp_fator = 1.;
const char *sp_eixo  = "#sigma [mb] por bin";

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

// Espectro de pT dos fotons ISOLADOS, um por origem (35 bins de 5 a 40
// GeV/c). Mesmo binning do hist_pt_direto_iso da foton_direto.cc.
TH1F *hist_pt_iso_pi0   = 0;
TH1F *hist_pt_iso_frag  = 0;
TH1F *hist_pt_iso_outro = 0;

// Massa invariante dos pares de fotons
TH1F *hist_massa_gg = 0;
TH1F *hist_massa_melhor[NORIGENS];

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
TCanvas *c13 = 0;   // espectro de pT dos isolados, por origem
TCanvas *c14 = 0;   // massa invariante

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

// eventos_equiv = 0  : eixo em SECAO DE CHOQUE, mb por bin. E' o PADRAO.
// eventos_equiv > 0  : eixo em NUMERO DE OCORRENCIAS -- quantas haveria numa
//                      rodada inclusiva de pThat > 4 com esse tanto de
//                      eventos. Passe 100000 para o equivalente ao que cada
//                      faixa gerou.
void soma_pthat(Double_t eventos_equiv = 0.)
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
    sp_eixo  = "#sigma [mb] por bin";
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

  // espectro de pT dos isolados, por origem
  hist_pt_iso_pi0   = soma_1d("hist_pt_iso_pi0");
  hist_pt_iso_frag  = soma_1d("hist_pt_iso_frag");
  hist_pt_iso_outro = soma_1d("hist_pt_iso_outro");

  // massa invariante dos pares de fotons
  hist_massa_gg = soma_1d("hist_massa_gg");
  for (Int_t o = 0; o < NORIGENS; o++)
  {
    sprintf(nome, "hist_massa_melhor_%s", curto_origem[o]);
    hist_massa_melhor[o] = soma_1d(nome);
  }

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

  // Os tres espectros por origem sao novos: se faltarem, o .root foi gerado
  // por uma versao da jatos_foton_v2.cc anterior a eles.
  if (!hist_pt_iso_pi0 || !hist_pt_iso_frag || !hist_pt_iso_outro ||
      !hist_massa_gg   || !hist_massa_melhor[2])
  {
    printf("  Faltam histogramas novos (hist_pt_iso_* ou hist_massa_*).\n");
    printf("  Rode a producao de novo:\n");
    printf("    ./producao.sh 100000 teto\n");
    return;
  }

  // ================================================================
  //  3) FECHA OS ARQUIVOS DE ENTRADA
  //  Pode fechar: todos os somados levaram SetDirectory(0) e estao
  //  na memoria, soltos de qualquer arquivo.
  // ================================================================
  for (Int_t i = 0; i < NRODADAS; i++) sp_arquivo[i]->Close();

  printf("  4 arquivos lidos e fechados. 68 histogramas somados.\n");

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
  hist_pt_iso_pi0  ->SetDirectory(arquivo);
  hist_pt_iso_frag ->SetDirectory(arquivo);
  hist_pt_iso_outro->SetDirectory(arquivo);
  hist_massa_gg    ->SetDirectory(arquivo);
  for (Int_t o = 0; o < NORIGENS; o++)
    hist_massa_melhor[o]->SetDirectory(arquivo);
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

  // ---- c13: espectro de pT dos fotons ISOLADOS, um painel por origem.
  // Este e' o RUIDO separado por origem, no mesmo binning do sinal da
  // foton_direto.cc -- da' para por os quatro no mesmo grafico depois.
  TH1F *espectro[NORIGENS] = {hist_pt_iso_pi0, hist_pt_iso_frag,
                              hist_pt_iso_outro};
  const char *tit_esp[NORIGENS] = {"#gamma isolado de #pi^{0}",
                                   "#gamma isolado de fragmentacao",
                                   "#gamma isolado de outro hadron"};
  Int_t cor_esp[NORIGENS] = {kBlue, kRed, kGreen + 2};

  c13 = new TCanvas("c13", "SOMA -- pT dos isolados por origem", 1500, 500);
  c13->Divide(3, 1);

  for (Int_t o = 0; o < NORIGENS; o++)
  {
    c13->cd(o + 1);
    gPad->SetLogy();
    gPad->SetLeftMargin(0.16);        // senao os expoentes saem cortados

    espectro[o]->SetTitle(tit_esp[o]);
    espectro[o]->SetXTitle("p_{T}^{#gamma} [GeV/c]");
    espectro[o]->SetYTitle(sp_eixo);
    espectro[o]->SetLineColor(cor_esp[o]);
    espectro[o]->SetLineWidth(2);
    espectro[o]->SetStats(0);
    espectro[o]->Draw("HIST E");
  }
  c13->Update();
  c13->SaveAs("soma_pt_iso_por_origem.png");

  // ---- c14: MASSA INVARIANTE, ja' somada e pesada por secao de choque
  c14 = new TCanvas("c14", "SOMA -- massa invariante", 1500, 500);
  c14->Divide(3, 1);

  c14->cd(1);
  gPad->SetLogy();  gPad->SetLeftMargin(0.16);
  hist_massa_gg->SetStats(0);
  hist_massa_gg->SetLineColor(kBlack);
  hist_massa_gg->SetLineWidth(2);
  hist_massa_gg->SetYTitle(sp_eixo);
  hist_massa_gg->Draw("HIST E");

  // zoom numa COPIA, para nao mexer no eixo do que vai para o .root
  c14->cd(2);
  gPad->SetLeftMargin(0.16);
  TH1F *zoom = (TH1F*) hist_massa_gg->Clone("sp_massa_zoom");
  zoom->SetDirectory(0);
  zoom->SetTitle("zoom no #pi^{0}");
  zoom->GetXaxis()->SetRangeUser(0., 0.3);
  zoom->SetLineColor(kBlue);
  zoom->Draw("HIST E");

  c14->cd(3);
  gPad->SetLogy();  gPad->SetLeftMargin(0.16);
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
      hist_massa_melhor[o]->SetYTitle(sp_eixo);
      hist_massa_melhor[o]->GetXaxis()->SetRangeUser(0., 0.5);
      hist_massa_melhor[o]->Draw(o == 0 ? "HIST E" : "HIST E SAME");
      leg_m->AddEntry(hist_massa_melhor[o], curto_origem[o], "l");
    }
    leg_m->Draw();
  }
  c14->Update();
  c14->SaveAs("soma_massa_invariante.png");

  // ================================================================
  //  RESUMO
  // ================================================================
  // Os somados ja' foram multiplicados por sp_fator. Para voltar a mb basta
  // desfazer essa multiplicacao -- assim a tabela sai nas DUAS unidades e
  // nao fica dependendo de em que modo a macro foi chamada.
  Bool_t outra_unidade = (sp_fator != 1.);

  printf("=============================================\n");
  printf("  ESPECTRO SOMADO (integral)\n");
  printf("    pT dos jatos    : %.6e mb", hist_pt_jato->Integral() / sp_fator);
  if (outra_unidade) printf("   = %.6g %s", hist_pt_jato->Integral(), sp_eixo);
  printf("\n");
  printf("    pT dos fotons   : %.6e mb", hist_pt_foton->Integral() / sp_fator);
  if (outra_unidade) printf("   = %.6g %s", hist_pt_foton->Integral(), sp_eixo);
  printf("\n");

  printf("---------------------------------------------\n");
  printf("  FOTONS ENERGETICOS SOMADOS, por origem\n");

  if (outra_unidade)
    printf("  %-8s %14s %14s %13s %13s\n", "origem",
           "iso [mb]", "nao-iso [mb]", "iso [eixo]", "nao-iso[eixo]");
  else
    printf("  %-8s %14s %14s\n", "origem", "iso [mb]", "nao-iso [mb]");

  for (Int_t o = 0; o < NORIGENS; o++)
  {
    Double_t s_iso = 0, s_nao = 0;
    for (Int_t f = 0; f < NFAIXAS; f++)
    {
      s_iso += hist_z_iso   [o][f]->Integral();
      s_nao += hist_z_naoiso[o][f]->Integral();
    }

    if (outra_unidade)
      printf("  %-8s %14.6e %14.6e %13.4f %13.4f\n", curto_origem[o],
             s_iso / sp_fator, s_nao / sp_fator, s_iso, s_nao);
    else
      printf("  %-8s %14.6e %14.6e\n", curto_origem[o],
             s_iso / sp_fator, s_nao / sp_fator);
  }

  // ---- MARCACAO POR MASSA INVARIANTE, pesada por secao de choque.
  // Esta e' a versao que vale: os contadores crus que cada rodada
  // imprime sao por faixa de pThat e nao podem ser somados direto.
  printf("---------------------------------------------\n");
  printf("  MARCACAO POR MASSA INVARIANTE\n");
  printf("  Marcado = existe parceiro com |m - m_pi0| < %.0f MeV\n",
         JANELA_PI0 * 1000.);
  // ATENCAO ao denominador. O hist_massa_melhor so' e' preenchido quando
  // existe ALGUM parceiro visivel. O foton cujo irmao ficou abaixo do
  // limiar nao entra nele -- e e' justamente esse o fundo IRREDUTIVEL.
  // Usar o hist_massa_melhor como denominador daria a eficiencia "entre
  // os que tem parceiro", que e' otimista e nao responde nada. O
  // denominador certo e' o total de isolados, do hist_pt_iso_*.
  printf("  %-8s %12s %12s %12s %8s %9s\n", "origem",
         "isolado[mb]", "c/ parc[mb]", "marcado[mb]",
         "marcado%", "s/ parc%");

  {
    // na ordem de curto_origem: pi0, outro hadron, fragmentacao
    TH1F *iso_org[NORIGENS] = {hist_pt_iso_pi0, hist_pt_iso_outro,
                               hist_pt_iso_frag};

    Int_t b_lo = hist_massa_melhor[0]->FindBin(M_PI0 - JANELA_PI0);
    Int_t b_hi = hist_massa_melhor[0]->FindBin(M_PI0 + JANELA_PI0);

    Double_t antes_tot = 0, sobra_tot = 0, sobra_frag = 0, frag_tot = 0;

    for (Int_t o = 0; o < NORIGENS; o++)
    {
      // Integral(0, n+1) inclui o overflow: o foton de pT acima de 40
      // GeV/c existe e tem que entrar na conta.
      Int_t nb = iso_org[o]->GetNbinsX();
      Double_t tot  = iso_org[o]->Integral(0, nb + 1) / sp_fator;

      Double_t comp = hist_massa_melhor[o]->Integral() / sp_fator;
      Double_t marc = hist_massa_melhor[o]->Integral(b_lo, b_hi) / sp_fator;

      printf("  %-8s %12.4e %12.4e %12.4e %7.2f%% %8.2f%%\n",
             curto_origem[o], tot, comp, marc,
             tot ? 100. * marc / tot : 0.,
             tot ? 100. * (tot - comp) / tot : 0.);

      antes_tot += tot;
      sobra_tot += tot - marc;
      if (o == 2) { sobra_frag = tot - marc; frag_tot = tot; }
    }

    printf("  ------------------------------------------\n");
    printf("  amostra isolada ANTES : %.6e mb\n", antes_tot);
    printf("  amostra isolada DEPOIS: %.6e mb\n", sobra_tot);
    printf("  rejeitado pela marcacao: %6.2f%%\n",
           antes_tot ? 100. * (antes_tot - sobra_tot) / antes_tot : 0.);
    printf("  fragmentacao ANTES    : %6.2f%%\n",
           antes_tot ? 100. * frag_tot   / antes_tot : 0.);
    printf("  fragmentacao DEPOIS   : %6.2f%%\n",
           sobra_tot ? 100. * sobra_frag / sobra_tot : 0.);
  }

  // ---- espectro de pT dos isolados, por origem (os histogramas novos).
  // Esta integral conta TODO foton isolado, inclusive os sem jato, entao e'
  // naturalmente maior que a coluna "iso" da tabela de cima, que vem dos
  // histogramas de z e por isso exige um jato associado.
  printf("---------------------------------------------\n");
  printf("  ESPECTRO DE pT DOS ISOLADOS, 5 a 40 GeV/c\n");
  printf("  (conta tambem os que nao cairam em jato)\n");

  TH1F *esp_res[NORIGENS] = {hist_pt_iso_pi0, hist_pt_iso_frag,
                             hist_pt_iso_outro};
  const char *nome_res[NORIGENS] = {"pi0", "frag", "outro"};

  for (Int_t o = 0; o < NORIGENS; o++)
  {
    printf("  %-8s %14.6e mb", nome_res[o],
           esp_res[o]->Integral() / sp_fator);
    if (outra_unidade) printf("   = %.4f %s", esp_res[o]->Integral(), sp_eixo);
    printf("   (%.0f entradas brutas)\n", esp_res[o]->GetEntries());
  }

  printf("=============================================\n");

  // Sem Close(): fechar o TFile de saida destroi os histogramas que agora
  // pertencem a ele, e os canvases ficariam em branco. O Write() grava.
  arquivo->Write();
  printf("  Arquivo salvo: soma_pthat.root\n");
  printf("  14 PNGs salvos com prefixo soma_\n");
}
