// As 4 faixas de pThat sobrepostas, com a SOMA por cima.
//
// RODAR
//   root -l ve_faixas.cc        // as 4 faixas + a soma (magenta grossa)
//   root -l 've_faixas.cc(1)'   // so' as 4 faixas, sem a soma
//
// Cores: preto, azul, verde e vermelho para as faixas, como nos outros
// visualizadores. A soma e' magenta e mais grossa, para destacar.
//
// SAO 42 HISTOGRAMAS em 8 canvases:
//   v1  espectros gerais   pt_jato, pt_foton, n_part, soma_cone
//   v2  decaimento do pi0  assimetria, z do foton lider
//   v3  v4  v5   z dos ISOLADOS     de pi0 / outro hadron / fragmentacao
//   v6  v7  v8   z dos NAO isolados de pi0 / outro hadron / fragmentacao
// Cada painel de v3 a v8 e' uma faixa de pT do foton (6 por canvas).
//
// Todas as faixas ja' vem multiplicadas pelo seu peso sigma/N, lido do bin 5
// do hist_sigma, entao o eixo vertical e' SECAO DE CHOQUE em mb por bin --
// nao contagem. Esse e' o padrao. Para ver o eixo em numero de ocorrencias,
// passe o numero de eventos equivalente no segundo argumento:
//   root -l ve_faixas.cc              // mb  (padrao)
//   root -l 've_faixas.cc(0, 100000)' // ocorrencias em 100 mil eventos
//
// PARA OS DEMAIS GRAFICOS:
//   so' a soma, nos 12 canvases originais    root -l soma_pthat.cc
//   uma faixa inteira, com os 2D             root -l 'abre_jatos_foton_v2.cc("jatos_foton_pthat15_inf.root")'

#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TGaxis.h"
#include "TLatex.h"

const Int_t NRODADAS = 4;
const Int_t NFAIXAS  = 6;
const Int_t NORIGENS = 3;

const char *vf_sufixo[NRODADAS] = {"pthat04_07", "pthat07_10",
                                   "pthat10_15", "pthat15_inf"};
const char *vf_rotulo[NRODADAS] = {"pThat 4-7 GeV",   "pThat 7-10 GeV",
                                   "pThat 10-15 GeV", "pThat > 15 GeV"};

// As cores de sempre: preto, azul, verde, vermelho
Int_t vf_cor[NRODADAS] = {kBlack, kBlue, kGreen + 2, kRed};
Int_t vf_cor_soma = kMagenta + 2;

TFile    *vf_arquivo[NRODADAS];
Double_t  vf_peso   [NRODADAS];
Bool_t    vf_com_soma = kTRUE;

// Fator final e rotulo do eixo vertical. Em "ocorrencias" o fator converte
// mb de volta para contagem; em "mb" vale 1 e nada muda.
Double_t    vf_fator = 1.;
const char *vf_eixo  = "#sigma [mb] por bin";

TCanvas *v1 = 0, *v2 = 0, *v3 = 0, *v4 = 0;
TCanvas *v5 = 0, *v6 = 0, *v7 = 0, *v8 = 0;

// Le um histograma das 4 rodadas, ja' multiplicado pelo peso de cada uma,
// e monta a soma em h[NRODADAS] (uma posicao extra no fim do vetor).
// O SetDirectory(0) solta os clones do arquivo, para sobreviverem.
void vf_le(const char *nome_hist, TH1F **h)
{
  char nome[200];

  for (Int_t i = 0; i < NRODADAS; i++)
  {
    TH1F *orig = (TH1F*) vf_arquivo[i]->Get(nome_hist);
    sprintf(nome, "vf_%s_%d", nome_hist, i);
    h[i] = (TH1F*) orig->Clone(nome);
    h[i]->SetDirectory(0);
    h[i]->Scale(vf_peso[i] * vf_fator);
  }

  sprintf(nome, "vf_%s_soma", nome_hist);
  h[NRODADAS] = (TH1F*) h[0]->Clone(nome);
  h[NRODADAS]->SetDirectory(0);
  for (Int_t i = 1; i < NRODADAS; i++) h[NRODADAS]->Add(h[i]);
}

// Mesma coisa, mas somando as 6 faixas de pT do foton de um conjunto z.
void vf_le_z(const char *conjunto, const char *origem, Int_t faixa, TH1F **h)
{
  char nome[200], chave[200];

  for (Int_t i = 0; i < NRODADAS; i++)
  {
    sprintf(chave, "hist_z_%s_%s_%d", conjunto, origem, faixa);
    TH1F *orig = (TH1F*) vf_arquivo[i]->Get(chave);
    sprintf(nome, "vf_z_%s_%s_%d_%d", conjunto, origem, faixa, i);
    h[i] = (TH1F*) orig->Clone(nome);
    h[i]->SetDirectory(0);
    h[i]->Scale(vf_peso[i] * vf_fator);
  }

  sprintf(nome, "vf_z_%s_%s_%d_soma", conjunto, origem, faixa);
  h[NRODADAS] = (TH1F*) h[0]->Clone(nome);
  h[NRODADAS]->SetDirectory(0);
  for (Int_t i = 1; i < NRODADAS; i++) h[NRODADAS]->Add(h[i]);
}

// Desenha as 4 faixas sobrepostas no painel atual, e a soma por cima
// se vf_com_soma estiver ligado.
void vf_desenha(TH1F **h, Int_t logy)
{
  char txt[150];

  gPad->SetLeftMargin(0.16);

  // Quantos desenhar: 4 faixas, ou 5 contando a soma
  Int_t n = vf_com_soma ? NRODADAS + 1 : NRODADAS;

  Double_t topo = 0, total = 0;
  for (Int_t i = 0; i < n; i++)
  {
    if (h[i]->GetMaximum() > topo) topo = h[i]->GetMaximum();
    if (i < NRODADAS) total += h[i]->GetEntries();
  }

  // Painel sem nenhuma entrada: escreve isso em vez de deixar um eixo
  // de 0 a 1 que nao quer dizer nada.
  if (total <= 0)
  {
    h[0]->SetYTitle(vf_eixo);
    h[0]->SetMaximum(1.);
    h[0]->SetStats(0);
    h[0]->Draw("HIST");

    TLatex *aviso = new TLatex();
    aviso->SetNDC();
    aviso->SetTextAlign(22);
    aviso->SetTextSize(0.07);
    aviso->SetTextColor(kGray + 2);
    aviso->DrawLatex(0.5, 0.5, "sem entradas nas 4 faixas");
    return;
  }

  if (logy) gPad->SetLogy();

  TLegend *leg = new TLegend(0.45, 0.62, 0.89, 0.89);
  leg->SetTextSize(0.032);

  for (Int_t i = 0; i < n; i++)
  {
    Bool_t eh_soma = (i == NRODADAS);

    h[i]->SetYTitle(vf_eixo);
    h[i]->SetStats(0);
    h[i]->SetFillStyle(0);
    h[i]->SetLineColor(eh_soma ? vf_cor_soma : vf_cor[i]);
    h[i]->SetLineWidth(eh_soma ? 3 : 2);

    // O maximo tem que ser o maior de todos, senao o primeiro desenhado
    // define a escala e os outros saem cortados.
    h[i]->SetMaximum(topo * 1.5);
    if (i == 0) h[i]->Draw("HIST E");
    else        h[i]->Draw("HIST E SAME");

    if (eh_soma) sprintf(txt, "SOMA das 4");
    else sprintf(txt, "%s  (%.0f ent.)", vf_rotulo[i], h[i]->GetEntries());
    leg->AddEntry(h[i], txt, "l");
  }
  leg->Draw();
}

// sem_soma      : 1 esconde a curva da soma, deixando so' as 4 faixas
// eventos_equiv : 0 (PADRAO) deixa o eixo em SECAO DE CHOQUE, mb por bin.
//                 Um valor > 0 poe o eixo em NUMERO DE OCORRENCIAS,
//                 equivalente a uma rodada inclusiva de pThat > 4 com esse
//                 tanto de eventos.
void ve_faixas(Int_t sem_soma = 0, Double_t eventos_equiv = 0.)
{
  const char *curto[NORIGENS]  = {"pi0", "outro", "frag"};
  const char *bonito[NORIGENS] = {"#pi^{0}", "outro hadron", "fragmentacao"};

  char nome[200];
  vf_com_soma = !sem_soma;
  TGaxis::SetMaxDigits(3);

  // ---- abre os 4 arquivos e le os pesos
  printf("\n============== FAIXAS DE pThat ==============\n");
  printf("  %-16s %12s %16s\n", "faixa", "sigma[mb]", "peso[mb/ev]");

  Double_t sigma_total = 0;

  for (Int_t i = 0; i < NRODADAS; i++)
  {
    sprintf(nome, "jatos_foton_%s.root", vf_sufixo[i]);
    vf_arquivo[i] = TFile::Open(nome);

    if (!vf_arquivo[i] || vf_arquivo[i]->IsZombie())
    {
      printf("  ERRO: nao consegui abrir %s\n", nome);
      printf("  Rode  ./producao.sh  antes.\n");
      return;
    }

    TH1D *info = (TH1D*) vf_arquivo[i]->Get("hist_sigma");
    if (!info) { printf("  ERRO: %s sem hist_sigma\n", nome); return; }

    vf_peso[i] = info->GetBinContent(5);
    sigma_total += info->GetBinContent(1);
    printf("  %-16s %12.4f %16.6e\n",
           vf_rotulo[i], info->GetBinContent(1), vf_peso[i]);
  }

  printf("---------------------------------------------\n");

  if (eventos_equiv > 0 && sigma_total > 0)
  {
    vf_fator = eventos_equiv / sigma_total;
    vf_eixo  = "ocorrencias";
    printf("  EIXO EM OCORRENCIAS, equivalente a %.0f eventos de pThat > 4\n",
           eventos_equiv);
  }
  else
  {
    vf_fator = 1.;
    vf_eixo  = "#sigma [mb] por bin";
    printf("  EIXO EM mb (milibarn), ou seja, secao de choque.\n");
  }

  if (vf_com_soma)
    printf("  4 faixas sobrepostas + a SOMA em magenta, mais grossa.\n");
  else
    printf("  So' as 4 faixas. Para a soma, rode sem argumento.\n");

  // ---- v1: espectros gerais
  TH1F *a[NRODADAS + 1], *b[NRODADAS + 1], *c[NRODADAS + 1], *d[NRODADAS + 1];
  vf_le("hist_pt_jato",   a);
  vf_le("hist_pt_foton",  b);
  vf_le("hist_n_part",    c);
  vf_le("hist_soma_cone", d);

  v1 = new TCanvas("v1", "Faixas de pThat -- espectros", 1300, 1000);
  v1->Divide(2, 2);
  v1->cd(1); vf_desenha(a, 1);
  v1->cd(2); vf_desenha(b, 1);
  v1->cd(3); vf_desenha(c, 0);
  v1->cd(4); vf_desenha(d, 1);
  v1->Update();
  v1->SaveAs("faixas_gerais.png");

  // ---- v2: decaimento do pi0
  TH1F *e[NRODADAS + 1], *g[NRODADAS + 1];
  vf_le("hist_assimetria", e);
  vf_le("hist_z_gama",     g);

  v2 = new TCanvas("v2", "Faixas de pThat -- decaimento do pi0", 1300, 500);
  v2->Divide(2, 1);
  v2->cd(1); vf_desenha(e, 0);
  v2->cd(2); vf_desenha(g, 1);
  v2->Update();
  v2->SaveAs("faixas_assimetria.png");

  // ---- v3 a v8: os 36 de z, uma origem por canvas, 6 faixas de pT cada
  TCanvas *canvas_z[6];
  const char *conjunto[2] = {"iso", "naoiso"};
  const char *titulo_c[2] = {"ISOLADO", "NAO isolado"};

  for (Int_t k = 0; k < 2; k++)
  {
    for (Int_t o = 0; o < NORIGENS; o++)
    {
      Int_t idx = k * NORIGENS + o;
      sprintf(nome, "v%d", 3 + idx);

      char titulo[250];
      sprintf(titulo, "z: %s de %s -- 4 faixas", titulo_c[k], bonito[o]);
      canvas_z[idx] = new TCanvas(nome, titulo, 1500, 1000);
      canvas_z[idx]->Divide(3, 2);

      for (Int_t f = 0; f < NFAIXAS; f++)
      {
        TH1F *h[NRODADAS + 1];
        vf_le_z(conjunto[k], curto[o], f, h);
        canvas_z[idx]->cd(f + 1);
        vf_desenha(h, 0);
      }
      canvas_z[idx]->Update();

      sprintf(nome, "faixas_z_%s_%s.png", conjunto[k], curto[o]);
      canvas_z[idx]->SaveAs(nome);
    }
  }

  v3 = canvas_z[0];  v4 = canvas_z[1];  v5 = canvas_z[2];
  v6 = canvas_z[3];  v7 = canvas_z[4];  v8 = canvas_z[5];

  printf("---------------------------------------------\n");
  printf("  42 histogramas em 8 canvases:\n");
  printf("    v1 espectros    v2 decaimento do pi0\n");
  printf("    v3 v4 v5  z ISOLADO     de pi0 / outro / frag\n");
  printf("    v6 v7 v8  z NAO isolado de pi0 / outro / frag\n");
  printf("  Os de 2 dimensoes estao na soma_pthat.cc e no\n");
  printf("  abre_jatos_foton_v2.cc de cada faixa.\n");
  printf("=============================================\n");
}
