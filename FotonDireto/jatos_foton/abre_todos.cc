// Abre as 4 rodadas de pThat AO MESMO TEMPO e desenha TODOS os histogramas
// de uma dimensao com as quatro versoes SOBREPOSTAS, para comparar de olho.
//
// RODAR
//   root -l abre_todos.cc            // normalizado pela secao de choque
//   root -l 'abre_todos.cc(kFALSE)'  // contagem bruta, sem normalizar
//
// NORMALIZADO (o padrao) multiplica cada rodada pelo seu peso w = sigma/N,
// lido do hist_sigma de cada arquivo. E' o unico jeito de comparar: sem isso
// as quatro curvas so' mostram quantos eventos cada uma gerou, nao fisica.
//
// SAO 42 HISTOGRAMAS SOBREPOSTOS, em 9 canvases:
//   tod_c1  espectros gerais       pt_jato, pt_foton, n_part, soma_cone
//   tod_c2  decaimento do pi0      assimetria, z do foton lider
//   tod_c3  as 4 secoes de choque
//   tod_c4  z dos ISOLADOS de pi0           (6 faixas)
//   tod_c5  z dos ISOLADOS de outro hadron  (6 faixas)
//   tod_c6  z dos ISOLADOS de fragmentacao  (6 faixas)
//   tod_c7  z dos NAO isolados de pi0           (6 faixas)
//   tod_c8  z dos NAO isolados de outro hadron  (6 faixas)
//   tod_c9  z dos NAO isolados de fragmentacao  (6 faixas)
//
// OS DE DUAS DIMENSOES NAO DA' PARA SOBREPOR: quatro mapas de cor um por
// cima do outro nao se leem. Para esses existe a funcao tod_mapa, que
// desenha as 4 rodadas LADO A LADO. Chame no prompt:
//   tod_mapa("hist_eta_phi")
//   tod_mapa("hist_td_0")     // lider x nao-lider, faixa 5-8 GeV
//   tod_mapa("hist_AZ_2")     // z e assimetria, faixa 10-12 GeV
//
// Os ponteiros ficam disponiveis no prompt depois que a macro termina:
//   tod_z_iso[2][0]           // fragmentacao, faixa 5-8: vetor das 4 rodadas
//   tod_pt_foton[3]->Draw()   // rodada pthat15
//   tod_c4->SaveAs("x.png")
//
// Os arquivos NAO sao fechados de proposito: fechar um TFile destroi os
// histogramas que pertencem a ele e os canvases ficariam em branco.
//
// Nomes com prefixo tod_ para nao colidir com o abre_jatos_foton_v2.cc,
// caso os dois sejam carregados na mesma sessao do ROOT.

#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TROOT.h"

const Int_t NRODADAS = 4;
const Int_t NFAIXAS  = 6;
const Int_t NORIGENS = 3;   // 0 = pi0, 1 = outro hadron, 2 = fragmentacao

TFile *tod_arquivo[NRODADAS];

// Os vetores tem a rodada como ULTIMO indice de proposito: assim
// tod_z_iso[o][f] ja' e' o vetor das 4 rodadas, pronto para sobrepor.
TH1F *tod_pt_jato   [NRODADAS];
TH1F *tod_pt_foton  [NRODADAS];
TH1F *tod_n_part    [NRODADAS];
TH1F *tod_soma_cone [NRODADAS];
TH1F *tod_assimetria[NRODADAS];
TH1F *tod_z_gama    [NRODADAS];

TH1F *tod_z_iso   [NORIGENS][NFAIXAS][NRODADAS];
TH1F *tod_z_naoiso[NORIGENS][NFAIXAS][NRODADAS];

TH1D *tod_sigmas = 0;

TCanvas *tod_c1 = 0;
TCanvas *tod_c2 = 0;
TCanvas *tod_c3 = 0;
TCanvas *tod_c4 = 0;
TCanvas *tod_c5 = 0;
TCanvas *tod_c6 = 0;
TCanvas *tod_c7 = 0;
TCanvas *tod_c8 = 0;
TCanvas *tod_c9 = 0;

// Preenchido pela abre_todos conforme o conjunto de arquivos escolhido
const char *tod_rotulo[NRODADAS] = {"", "", "", ""};
Int_t tod_cor[NRODADAS] = {kBlack, kBlue, kGreen + 2, kRed};

// Guarda se os histogramas foram multiplicados pelo peso sigma/N. O eixo
// vertical muda de significado conforme isso, e o rotulo tem que dizer qual.
Bool_t tod_normalizado = kTRUE;

// Desenha o mesmo histograma das 4 rodadas sobrepostos no painel atual.
// logy: 1 liga a escala logaritmica no eixo vertical.
void tod_sobrepoe(TH1F **h, Int_t logy)
{
  char txt[120];

  // Os paineis de um canvas dividido 3x2 sao estreitos e os numeros do eixo
  // vertical saem cortados na borda. Mais margem resolve.
  gPad->SetLeftMargin(0.16);

  // O rotulo do eixo vertical TEM que dizer o que esta' plotado. Depois do
  // Scale(sigma/N) o conteudo nao e' mais contagem de particulas, e sim
  // secao de choque em mb. Deixar escrito "#gamma" seria mentira.
  const char *eixo_y = tod_normalizado ? "mb por bin" : "contagem";

  Double_t topo = 0, total = 0;
  for (Int_t i = 0; i < NRODADAS; i++)
  {
    if (!h[i]) return;
    if (h[i]->GetMaximum() > topo) topo = h[i]->GetMaximum();
    total += h[i]->GetEntries();
  }

  // Painel sem nenhuma entrada nas 4 rodadas: desenha o quadro e escreve
  // isso na cara, em vez de deixar um eixo de 0 a 1 que nao quer dizer nada.
  if (total <= 0)
  {
    h[0]->SetYTitle(eixo_y);
    h[0]->SetMaximum(1.);
    h[0]->Draw("HIST");

    TLatex *aviso = new TLatex();
    aviso->SetNDC();
    aviso->SetTextAlign(22);
    aviso->SetTextSize(0.07);
    aviso->SetTextColor(kGray + 2);
    aviso->DrawLatex(0.5, 0.5, "sem entradas nas 4 rodadas");
    return;
  }

  if (logy) gPad->SetLogy();

  TLegend *leg = new TLegend(0.45, 0.65, 0.88, 0.88);
  leg->SetTextSize(0.033);

  for (Int_t i = 0; i < NRODADAS; i++)
  {
    h[i]->SetYTitle(eixo_y);

    // O maximo tem que ser o maior das quatro, senao a primeira curva
    // desenhada define a escala e as outras saem cortadas.
    h[i]->SetMaximum(topo * 1.3);
    if (i == 0) h[i]->Draw("HIST E");
    else        h[i]->Draw("HIST E SAME");

    // O numero de entradas BRUTAS entra na legenda. Sem isso nao da' para
    // saber se uma curva e' uma distribuicao ou tres particulas soltas.
    sprintf(txt, "%s  (%.0f entradas)", tod_rotulo[i], h[i]->GetEntries());
    leg->AddEntry(h[i], txt, "l");
  }
  leg->Draw();
}

// Desenha um histograma de DUAS dimensoes das 4 rodadas LADO A LADO,
// porque sobrepor mapa de cor nao se le. Chame no prompt do ROOT.
void tod_mapa(const char *nome_hist)
{
  char nome[200];

  if (!tod_arquivo[0])
  {
    printf("  Rode a abre_todos primeiro.\n");
    return;
  }

  sprintf(nome, "tod_mapa_%s", nome_hist);
  TCanvas *c = new TCanvas(nome, nome_hist, 1400, 1000);
  c->Divide(2, 2);

  for (Int_t i = 0; i < NRODADAS; i++)
  {
    TH2F *h = (TH2F*) tod_arquivo[i]->Get(nome_hist);
    if (!h) { printf("  %s nao existe na rodada %d\n", nome_hist, i); continue; }

    c->cd(i + 1);
    h->SetStats(0);
    sprintf(nome, "%s -- %s", nome_hist, tod_rotulo[i]);
    h->SetTitle(nome);
    h->Draw("COLZ");
  }
  c->Update();
  printf("  %s: 4 rodadas lado a lado\n", nome_hist);
}

// exclusivas = kTRUE  le as faixas fechadas  (pthat04_07, 07_10, 10_15, 15_inf)
// exclusivas = kFALSE le as rodadas encaixadas (pthat04, 07, 10, 15), que sao
//                     as antigas, sem pTHatMax. Elas nao podem ser somadas.
void abre_todos(Bool_t normalizar = kTRUE, Bool_t exclusivas = kTRUE)
{
  const char *sufixo[NRODADAS];
  if (exclusivas)
  {
    sufixo[0] = "pthat04_07";  sufixo[1] = "pthat07_10";
    sufixo[2] = "pthat10_15";  sufixo[3] = "pthat15_inf";
    tod_rotulo[0] = "pThat 4-7 GeV";   tod_rotulo[1] = "pThat 7-10 GeV";
    tod_rotulo[2] = "pThat 10-15 GeV"; tod_rotulo[3] = "pThat > 15 GeV";
  }
  else
  {
    sufixo[0] = "pthat04";  sufixo[1] = "pthat07";
    sufixo[2] = "pthat10";  sufixo[3] = "pthat15";
    tod_rotulo[0] = "pThat > 4 GeV";  tod_rotulo[1] = "pThat > 7 GeV";
    tod_rotulo[2] = "pThat > 10 GeV"; tod_rotulo[3] = "pThat > 15 GeV";
  }
  const char *curto[NORIGENS]  = {"pi0", "outro", "frag"};
  const char *nome_origem[NORIGENS] = {"#pi^{0}", "outro hadron", "fragmentacao"};

  char nome[200], titulo[300];

  tod_normalizado = normalizar;   // o rotulo do eixo y depende disso

  // Com peso em mb os valores sao da ordem de 1e-3. Sem isso o ROOT escreve
  // "0.0018" e os digitos nao cabem na lateral de um painel estreito.
  TGaxis::SetMaxDigits(3);

  tod_sigmas = new TH1D("tod_sigmas",
      "secao de choque de cada rodada;;#sigma [mb]", NRODADAS, 0, NRODADAS);
  tod_sigmas->SetDirectory(0);

  printf("\n============== ARQUIVOS ABERTOS ==============\n");
  printf("  %-9s %12s %10s %14s\n", "rodada", "sigma[mb]", "N", "peso[mb/ev]");

  for (Int_t i = 0; i < NRODADAS; i++)
  {
    sprintf(nome, "jatos_foton_%s.root", sufixo[i]);
    tod_arquivo[i] = TFile::Open(nome);

    if (!tod_arquivo[i] || tod_arquivo[i]->IsZombie())
    {
      printf("  ERRO: nao consegui abrir %s\n", nome);
      printf("  Rode as 4 producoes antes de abrir.\n");
      return;
    }

    TH1D *info = (TH1D*) tod_arquivo[i]->Get("hist_sigma");
    if (!info)
    {
      printf("  ERRO: %s nao tem hist_sigma (versao antiga da macro?)\n", nome);
      return;
    }

    Double_t sigma = info->GetBinContent(1);
    Double_t n_ev  = info->GetBinContent(3);
    Double_t peso  = info->GetBinContent(5);

    tod_sigmas->SetBinContent(i + 1, sigma);
    tod_sigmas->SetBinError  (i + 1, info->GetBinContent(2));
    tod_sigmas->GetXaxis()->SetBinLabel(i + 1, sufixo[i]);

    printf("  %-9s %12.4f %10.0f %14.4e\n", sufixo[i], sigma, n_ev, peso);

    Double_t escala = normalizar ? peso : 1.;

    // ---- os 6 gerais
    const char *lista[6] = {"hist_pt_jato", "hist_pt_foton", "hist_n_part",
                            "hist_soma_cone", "hist_assimetria", "hist_z_gama"};
    TH1F *destino[6];

    for (Int_t k = 0; k < 6; k++)
    {
      TH1F *orig = (TH1F*) tod_arquivo[i]->Get(lista[k]);
      sprintf(nome, "tod_%s_%s", lista[k], sufixo[i]);

      // Clonar evita mexer na copia que pertence ao arquivo.
      // SetDirectory(0) deixa o clone na memoria, solto de qualquer TFile.
      destino[k] = (TH1F*) orig->Clone(nome);
      destino[k]->SetDirectory(0);
      destino[k]->SetLineColor(tod_cor[i]);
      destino[k]->SetLineWidth(2);
      destino[k]->SetStats(0);      // 4 caixas de estatistica tapariam o painel
      destino[k]->Scale(escala);
    }

    tod_pt_jato   [i] = destino[0];
    tod_pt_foton  [i] = destino[1];
    tod_n_part    [i] = destino[2];
    tod_soma_cone [i] = destino[3];
    tod_assimetria[i] = destino[4];
    tod_z_gama    [i] = destino[5];

    // ---- os 36 de z: 2 conjuntos x 3 origens x 6 faixas
    for (Int_t o = 0; o < NORIGENS; o++)
    {
      for (Int_t f = 0; f < NFAIXAS; f++)
      {
        sprintf(nome, "hist_z_iso_%s_%d", curto[o], f);
        TH1F *a = (TH1F*) tod_arquivo[i]->Get(nome);
        sprintf(nome, "tod_z_iso_%s_%d_%s", curto[o], f, sufixo[i]);
        tod_z_iso[o][f][i] = (TH1F*) a->Clone(nome);
        tod_z_iso[o][f][i]->SetDirectory(0);
        tod_z_iso[o][f][i]->SetLineColor(tod_cor[i]);
        tod_z_iso[o][f][i]->SetLineWidth(2);
        tod_z_iso[o][f][i]->SetStats(0);
        tod_z_iso[o][f][i]->Scale(escala);

        sprintf(nome, "hist_z_naoiso_%s_%d", curto[o], f);
        TH1F *b = (TH1F*) tod_arquivo[i]->Get(nome);
        sprintf(nome, "tod_z_naoiso_%s_%d_%s", curto[o], f, sufixo[i]);
        tod_z_naoiso[o][f][i] = (TH1F*) b->Clone(nome);
        tod_z_naoiso[o][f][i]->SetDirectory(0);
        tod_z_naoiso[o][f][i]->SetLineColor(tod_cor[i]);
        tod_z_naoiso[o][f][i]->SetLineWidth(2);
        tod_z_naoiso[o][f][i]->SetStats(0);
        tod_z_naoiso[o][f][i]->Scale(escala);
      }
    }
  }

  printf("---------------------------------------------\n");
  if (normalizar)
    printf("  NORMALIZADO: cada rodada multiplicada pelo seu peso sigma/N.\n"
           "  O eixo vertical e' secao de choque, nao contagem.\n");
  else
    printf("  BRUTO: contagem crua, sem peso. Serve para ver a estatistica\n"
           "  de cada rodada, NAO para comparar fisica.\n");

  // ---- c1: os espectros gerais
  tod_c1 = new TCanvas("tod_c1", "4 rodadas -- espectros", 1300, 1000);
  tod_c1->Divide(2, 2);
  tod_c1->cd(1); tod_sobrepoe(tod_pt_jato,   1);
  tod_c1->cd(2); tod_sobrepoe(tod_pt_foton,  1);
  tod_c1->cd(3); tod_sobrepoe(tod_n_part,    0);
  tod_c1->cd(4); tod_sobrepoe(tod_soma_cone, 1);
  tod_c1->Update();

  // ---- c2: decaimento do pi0
  tod_c2 = new TCanvas("tod_c2", "4 rodadas -- decaimento do pi0", 1300, 500);
  tod_c2->Divide(2, 1);
  tod_c2->cd(1); tod_sobrepoe(tod_assimetria, 0);
  tod_c2->cd(2); tod_sobrepoe(tod_z_gama,     1);
  tod_c2->Update();

  // ---- c3: as secoes de choque
  tod_c3 = new TCanvas("tod_c3", "Secoes de choque", 800, 600);
  gPad->SetLogy();
  tod_sigmas->SetMarkerStyle(20);
  tod_sigmas->SetStats(0);
  tod_sigmas->Draw("P E");
  tod_c3->Update();

  // ---- c4 a c9: os 36 de z, uma origem por canvas, 6 faixas por canvas
  TCanvas *canvas_z[6];

  for (Int_t o = 0; o < NORIGENS; o++)
  {
    sprintf(nome, "tod_c%d", 4 + o);
    sprintf(titulo, "z: ISOLADO de %s -- 4 rodadas", nome_origem[o]);
    canvas_z[o] = new TCanvas(nome, titulo, 1500, 1000);
    canvas_z[o]->Divide(3, 2);
    for (Int_t f = 0; f < NFAIXAS; f++)
    {
      canvas_z[o]->cd(f + 1);
      tod_sobrepoe(tod_z_iso[o][f], 0);
    }
    canvas_z[o]->Update();

    sprintf(nome, "tod_c%d", 7 + o);
    sprintf(titulo, "z: NAO isolado de %s -- 4 rodadas", nome_origem[o]);
    canvas_z[3 + o] = new TCanvas(nome, titulo, 1500, 1000);
    canvas_z[3 + o]->Divide(3, 2);
    for (Int_t f = 0; f < NFAIXAS; f++)
    {
      canvas_z[3 + o]->cd(f + 1);
      tod_sobrepoe(tod_z_naoiso[o][f], 0);
    }
    canvas_z[3 + o]->Update();
  }

  tod_c4 = canvas_z[0];  tod_c5 = canvas_z[1];  tod_c6 = canvas_z[2];
  tod_c7 = canvas_z[3];  tod_c8 = canvas_z[4];  tod_c9 = canvas_z[5];

  printf("---------------------------------------------\n");
  printf("  42 histogramas sobrepostos em 9 canvases:\n");
  printf("    tod_c1 espectros   tod_c2 pi0   tod_c3 sigmas\n");
  printf("    tod_c4/c5/c6  z ISOLADO     de pi0 / outro / frag\n");
  printf("    tod_c7/c8/c9  z NAO isolado de pi0 / outro / frag\n");
  printf("  Os de 2 dimensoes nao da' para sobrepor; use lado a lado:\n");
  printf("    tod_mapa(\"hist_eta_phi\")\n");
  printf("    tod_mapa(\"hist_td_0\")    tod_mapa(\"hist_AZ_2\")\n");
  printf("  .q para sair\n");
  printf("=============================================\n");
}
