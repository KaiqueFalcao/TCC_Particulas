// ============================================================================
//  pureza.cc  --  quanto da amostra isolada e' sinal de verdade?
// ============================================================================
//
// Responde a pergunta central do trabalho:
//
//     de tudo que passa no corte de isolamento, que fracao e' foton DIRETO
//     (o que eu quero medir) e que fracao e' ruido (o que me engana)?
//
//     pureza(pT) = S(pT) / [ S(pT) + B(pT) ]
//
// com S e B em SECAO DE CHOQUE, nao em contagem. Tem que ser em secao de
// choque: sinal e ruido vem de amostras diferentes, geradas com numeros de
// eventos diferentes, e so' a secao de choque as poe na mesma escala fisica.
//
// Calcula a pureza em DOIS momentos:
//
//   ANTES  -- so' com o corte de isolamento
//   DEPOIS -- tambem com a marcacao por massa invariante, que joga fora o
//             foton que forma massa de pi0 com algum parceiro visivel
//
// ----------------------------------------------------------------------------
//  SUPOSICAO DECLARADA: a marcacao nao remove foton direto
// ----------------------------------------------------------------------------
// A marcacao por massa invariante e' aplicada SO' AO RUIDO. Entre as duas
// etapas da pureza, o sinal e' o mesmo histograma.
//
// O motivo e' que a massa invariante existe para retirar foton de DECAIMENTO
// de pi0, e na amostra de sinal a identificacao pela verdade da simulacao ja'
// garante que todo foton e' direto: nao ha' decaimento para retirar. Aplicar
// a marcacao la' mediria outra coisa -- quantos fotons diretos cairiam na
// janela do pi0 POR ACASO, formando par com um foton qualquer do evento.
//
// Esse efeito deve ser desprezivel. Para um par dar a massa do pi0 os dois
// fotons precisam estar muito proximos em angulo, e o foton direto isolado,
// por definicao, tem pouca coisa em volta.
//
// COMO CONFERIR O TAMANHO DA SUPOSICAO: o foton de FRAGMENTACAO do ruido e' o
// analogo mais proximo do direto -- tambem nao vem de decaimento e tambem
// sobrevive ao isolamento. A tabela de marcacao da soma_pthat.cc mostra
// quantos deles sao marcados por engano. Se esse numero for compativel com
// zero, a suposicao esta' sustentada pelos proprios dados.
//
// ----------------------------------------------------------------------------
//  POR QUE AS DUAS AMOSTRAS PRECISAM DAS MESMAS FAIXAS DE pThat
// ----------------------------------------------------------------------------
// Uma rodada de pTHatMin = X e' INCOMPLETA para qualquer observavel perto de
// X: faltam os eventos de pThat abaixo de X que flutuaram para cima. Se o
// ruido vier de pThat > 4 e o sinal de pThat > 5, a razao entre os dois mistura
// essa incompletude com a fisica, e a pureza sai errada justamente na regiao
// de pT baixo, que e' onde ha' mais estatistica.
//
// Por isso este codigo exige as MESMAS 4 faixas exclusivas nas duas amostras,
// cada uma pesada pela sua propria secao de choque:
//
//     h_total = soma_i ( sigma_i / N_i ) * h_i
//
// ----------------------------------------------------------------------------
//  RODAR
// ----------------------------------------------------------------------------
//   root -l pureza.cc          // bins de 5 GeV/c na curva de pureza
//   root -l 'pureza.cc(1)'     // sem reagrupar, 1 GeV/c por bin (ruidoso)
//
// Precisa de 8 arquivos, as 4 faixas de cada amostra:
//   jatos_foton_pthat04_07.root ...   (ruido,  ./producao.sh 1000000 teto)
//   foton_direto_pthat04_07.root ...  (sinal,  ./producao_direto.sh 100000)

#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TROOT.h"
#include "TLatex.h"

const Int_t NRODADAS = 4;
const Int_t NORIGENS = 3;   // 0 = pi0, 1 = outro hadron, 2 = fragmentacao

const char *pz_sufixo[NRODADAS] = {"pthat04_07", "pthat07_10",
                                   "pthat10_15", "pthat15_inf"};

// Histogramas em escopo de arquivo: continuam acessiveis no prompt do ROOT
TH1F *hist_sinal     = 0;          // foton direto isolado -- o MESMO nas
                                   // duas etapas, por suposicao
TH1F *hist_fundo     = 0;          // soma das 3 origens, ANTES
TH1F *hist_fundo_pos = 0;          // o mesmo, DEPOIS
TH1F *hist_ruido    [NORIGENS];    // pi0, outro hadron, fragmentacao, ANTES
TH1F *hist_ruido_pos[NORIGENS];    // o mesmo, DEPOIS

TH1F *hist_pureza        = 0;      // S/(S+B) antes
TH1F *hist_pureza_pos    = 0;      // S/(S+B) depois
TH1F *hist_pureza_prompt = 0;      // antes, contando fragmentacao como sinal

TCanvas *p1 = 0;
TCanvas *p2 = 0;

// Le um histograma das 4 faixas de uma amostra, cada uma multiplicada pelo
// peso sigma/N da sua faixa, e devolve a soma. O prefixo diz de que amostra:
// "jatos_foton" para o ruido, "foton_direto" para o sinal.
//
// Devolve 0 se faltar arquivo ou histograma, depois de dizer qual.
TH1F *soma_faixas(const char *prefixo, const char *nome_hist,
                  const char *apelido)
{
  TH1F *total = 0;
  char arq[200];

  for (Int_t i = 0; i < NRODADAS; i++)
  {
    sprintf(arq, "%s_%s.root", prefixo, pz_sufixo[i]);

    TFile *f = TFile::Open(arq);
    if (!f || f->IsZombie())
    {
      printf("  ERRO: nao consegui abrir %s\n", arq);
      if (total) delete total;
      return 0;
    }

    TH1D *info = (TH1D*) f->Get("hist_sigma");
    TH1F *h    = (TH1F*) f->Get(nome_hist);

    if (!info || !h)
    {
      printf("  ERRO: %s nao tem %s\n", arq, info ? nome_hist : "hist_sigma");
      f->Close();
      if (total) delete total;
      return 0;
    }

    Double_t peso = info->GetBinContent(5);   // bin 5 = sigma / n_aceitos

    // O clone nasce com nome temporario e so' depois e' rebatizado: clonar
    // ja' com o nome final o criaria dentro do TFile, com nome repetido.
    TH1F *copia = (TH1F*) h->Clone("pz_tmp");
    copia->SetDirectory(0);          // solta do arquivo, sobrevive ao Close
    copia->Scale(peso);              // contagem -> mb

    if (!total) { total = copia; total->SetName(apelido); }
    else        { total->Add(copia); delete copia; }

    f->Close();
  }
  return total;
}

// Preenche "saida" com a pureza bin a bin, a partir das quatro componentes.
//
//     P = X / (X + Y),  X = sinal,  Y = fundo
//
// Os erros sao propagados na mao:
//
//     dP/dX =  Y / (X+Y)^2        dP/dY = -X / (X+Y)^2
//     sigma_P^2 = [ Y^2 sigma_X^2 + X^2 sigma_Y^2 ] / (X+Y)^4
//
// NAO da' para usar TH1::Divide com a opcao "B" (binomial): binomial vale
// quando o numerador e' um subconjunto do denominador medido na MESMA
// amostra. Aqui sinal e fundo vem de amostras independentes.
//
// frag_e_sinal muda de convencao: se for verdadeiro, a fragmentacao conta
// como sinal em vez de fundo.
void calcula_pureza(TH1F *S, TH1F *Bp, TH1F *Bo, TH1F *Bf,
                    TH1F *saida, Bool_t frag_e_sinal)
{
  for (Int_t b = 1; b <= S->GetNbinsX(); b++)
  {
    Double_t s  = S ->GetBinContent(b), es  = S ->GetBinError(b);
    Double_t bp = Bp->GetBinContent(b), ebp = Bp->GetBinError(b);
    Double_t bo = Bo->GetBinContent(b), ebo = Bo->GetBinError(b);
    Double_t bf = Bf->GetBinContent(b), ebf = Bf->GetBinError(b);

    Double_t tot = s + bp + bo + bf;
    if (tot <= 0) continue;              // bin vazio nas duas amostras

    Double_t X, Y, eX2, eY2;

    if (frag_e_sinal)
    {
      X = s + bf;   eX2 = es*es + ebf*ebf;
      Y = bp + bo;  eY2 = ebp*ebp + ebo*ebo;
    }
    else
    {
      X = s;             eX2 = es*es;
      Y = bp + bo + bf;  eY2 = ebp*ebp + ebo*ebo + ebf*ebf;
    }

    // SE O FUNDO NAO TEM NENHUM EVENTO DE MONTE CARLO NESTE BIN, a pureza
    // nao e' 100%: e' DESCONHECIDA. Zero eventos gerados nao quer dizer
    // zero secao de choque, quer dizer que a amostra e' pequena demais para
    // dizer alguma coisa aqui. O ponto fica de fora do grafico.
    if (Y <= 0) continue;

    saida->SetBinContent(b, X / tot);
    saida->SetBinError  (b, TMath::Sqrt(Y*Y*eX2 + X*X*eY2) / (tot*tot));
  }
}

// rebin : quantos bins de 1 GeV/c agrupar na curva de pureza. O espectro
//         e' desenhado sempre na granularidade original; so' a pureza e'
//         reagrupada, porque ela e' uma razao e fica ilegivel com pouca
//         estatistica por bin.
void pureza(Int_t rebin = 5)
{
  const char *curto[NORIGENS] = {"pi0", "outro", "frag"};
  char nome_hist[100], apelido[100];

  printf("\n================= PUREZA =================\n");

  // ================================================================
  //  1) SINAL
  // ================================================================
  // Um histograma so'. Ele entra identico nas duas etapas da pureza: a
  // marcacao por massa invariante nao toca no sinal (ver o cabecalho).
  hist_sinal = soma_faixas("foton_direto", "hist_pt_direto_iso",
                           "hist_sinal");

  if (!hist_sinal)
  {
    printf("\n  Falta a producao do SINAL nas 4 faixas de pThat. Rode:\n\n");
    printf("    ./producao_direto.sh 100000\n\n");
    return;
  }

  // ================================================================
  //  2) RUIDO, antes e depois
  // ================================================================
  for (Int_t o = 0; o < NORIGENS; o++)
  {
    sprintf(nome_hist, "hist_pt_iso_%s", curto[o]);
    sprintf(apelido,   "hist_ruido_%s",  curto[o]);
    hist_ruido[o] = soma_faixas("jatos_foton", nome_hist, apelido);

    sprintf(nome_hist, "hist_pt_pos_%s",    curto[o]);
    sprintf(apelido,   "hist_ruido_pos_%s", curto[o]);
    hist_ruido_pos[o] = soma_faixas("jatos_foton", nome_hist, apelido);

    if (!hist_ruido[o] || !hist_ruido_pos[o])
    {
      printf("\n  Falta a producao do RUIDO, ou ela e' anterior a marcacao.\n");
      printf("  Rode:\n    ./producao.sh 1000000 teto\n");
      return;
    }
  }

  hist_fundo = (TH1F*) hist_ruido[0]->Clone("hist_fundo");
  hist_fundo->SetDirectory(0);
  hist_fundo->Add(hist_ruido[1]);
  hist_fundo->Add(hist_ruido[2]);

  hist_fundo_pos = (TH1F*) hist_ruido_pos[0]->Clone("hist_fundo_pos");
  hist_fundo_pos->SetDirectory(0);
  hist_fundo_pos->Add(hist_ruido_pos[1]);
  hist_fundo_pos->Add(hist_ruido_pos[2]);

  // ================================================================
  //  3) PUREZA
  // ================================================================
  // Copias reagrupadas, usadas so' na curva. O espectro do painel 1
  // continua na granularidade original, de 1 GeV/c.
  TH1F *S  = (TH1F*) hist_sinal    ->Clone("pz_S");   S ->SetDirectory(0);
  TH1F *Bp = (TH1F*) hist_ruido[0] ->Clone("pz_Bp");  Bp->SetDirectory(0);
  TH1F *Bo = (TH1F*) hist_ruido[1] ->Clone("pz_Bo");  Bo->SetDirectory(0);
  TH1F *Bf = (TH1F*) hist_ruido[2] ->Clone("pz_Bf");  Bf->SetDirectory(0);

  // O sinal da etapa "depois" e' o MESMO da etapa "antes".
  TH1F *Bpd = (TH1F*) hist_ruido_pos[0]->Clone("pz_Bpd"); Bpd->SetDirectory(0);
  TH1F *Bod = (TH1F*) hist_ruido_pos[1]->Clone("pz_Bod"); Bod->SetDirectory(0);
  TH1F *Bfd = (TH1F*) hist_ruido_pos[2]->Clone("pz_Bfd"); Bfd->SetDirectory(0);

  if (rebin > 1)
  {
    S->Rebin(rebin);  Bp ->Rebin(rebin); Bo ->Rebin(rebin); Bf ->Rebin(rebin);
    Bpd->Rebin(rebin); Bod->Rebin(rebin); Bfd->Rebin(rebin);
  }

  hist_pureza = (TH1F*) S->Clone("hist_pureza");
  hist_pureza->SetDirectory(0);  hist_pureza->Reset();

  hist_pureza_pos = (TH1F*) S->Clone("hist_pureza_pos");
  hist_pureza_pos->SetDirectory(0);  hist_pureza_pos->Reset();

  hist_pureza_prompt = (TH1F*) S->Clone("hist_pureza_prompt");
  hist_pureza_prompt->SetDirectory(0);  hist_pureza_prompt->Reset();

  calcula_pureza(S,  Bp,  Bo,  Bf,  hist_pureza,        kFALSE);
  calcula_pureza(S,  Bpd, Bod, Bfd, hist_pureza_pos,    kFALSE);
  calcula_pureza(S,  Bp,  Bo,  Bf,  hist_pureza_prompt, kTRUE);

  // ================================================================
  //  4) GRAFICOS
  // ================================================================
  // --- p1: os espectros, sinal e ruido no mesmo eixo
  p1 = new TCanvas("p1", "Sinal e ruido", 900, 700);
  gPad->SetLogy();
  gPad->SetLeftMargin(0.16);

  Int_t cor_r[NORIGENS] = {kBlue, kGreen + 2, kRed};

  hist_fundo->SetTitle("#gamma isolados: sinal e ruido");
  hist_fundo->SetXTitle("p_{T}^{#gamma} [GeV/c]");
  hist_fundo->SetYTitle("#sigma [mb] por bin");
  hist_fundo->SetStats(0);
  hist_fundo->SetLineColor(kBlack);
  hist_fundo->SetLineWidth(3);

  // O piso do eixo tem que caber o MENOR bin nao vazio das duas amostras.
  // Sem isso o ROOT escolhe o piso pelo histograma desenhado primeiro (o
  // fundo) e a cauda do sinal, que vai mais fundo, some embaixo do eixo.
  {
    Double_t piso = 1e30;
    for (Int_t b = 1; b <= hist_fundo->GetNbinsX(); b++)
    {
      Double_t cs = hist_sinal->GetBinContent(b);
      Double_t cb = hist_fundo->GetBinContent(b);
      if (cs > 0 && cs < piso) piso = cs;
      if (cb > 0 && cb < piso) piso = cb;
    }
    if (piso < 1e29) hist_fundo->SetMinimum(piso * 0.4);
  }

  hist_fundo->Draw("HIST E");

  for (Int_t o = 0; o < NORIGENS; o++)
  {
    hist_ruido[o]->SetLineColor(cor_r[o]);
    hist_ruido[o]->SetLineWidth(1);
    hist_ruido[o]->SetLineStyle(2);
    hist_ruido[o]->Draw("HIST E SAME");
  }

  hist_sinal->SetLineColor(kMagenta + 2);
  hist_sinal->SetLineWidth(3);
  hist_sinal->Draw("HIST E SAME");

  TLegend *l1 = new TLegend(0.50, 0.63, 0.89, 0.89);
  l1->SetTextSize(0.030);
  l1->AddEntry(hist_sinal,    "SINAL: #gamma direto", "l");
  l1->AddEntry(hist_fundo,    "RUIDO: soma das 3",    "l");
  l1->AddEntry(hist_ruido[0], "ruido: #pi^{0}",       "l");
  l1->AddEntry(hist_ruido[1], "ruido: outro hadron",  "l");
  l1->AddEntry(hist_ruido[2], "ruido: fragmentacao",  "l");
  l1->Draw();

  p1->Update();
  p1->SaveAs("pureza_espectros.png");

  // --- p2: a curva de pureza, antes e depois da marcacao
  p2 = new TCanvas("p2", "Pureza", 900, 700);
  gPad->SetLeftMargin(0.16);
  gPad->SetGridy();

  hist_pureza->SetTitle("pureza da amostra de #gamma isolados");
  hist_pureza->SetXTitle("p_{T}^{#gamma} [GeV/c]");
  hist_pureza->SetYTitle("pureza = S / (S + B)");
  hist_pureza->SetStats(0);
  hist_pureza->SetMinimum(0.);
  hist_pureza->SetMaximum(1.05);
  hist_pureza->SetLineColor(kMagenta + 2);
  hist_pureza->SetLineWidth(3);
  hist_pureza->SetMarkerStyle(20);
  hist_pureza->SetMarkerColor(kMagenta + 2);
  hist_pureza->Draw("E1");

  hist_pureza_pos->SetLineColor(kAzure + 2);
  hist_pureza_pos->SetLineWidth(3);
  hist_pureza_pos->SetMarkerStyle(21);
  hist_pureza_pos->SetMarkerColor(kAzure + 2);
  hist_pureza_pos->Draw("E1 SAME");

  TLegend *l2 = new TLegend(0.40, 0.76, 0.89, 0.89);
  l2->SetTextSize(0.030);
  l2->AddEntry(hist_pureza,     "so' isolamento",          "lp");
  l2->AddEntry(hist_pureza_pos, "isolamento + massa m_{#gamma#gamma}", "lp");
  l2->SetHeader("sinal identico nas duas");
  l2->Draw();

  p2->Update();
  p2->SaveAs("pureza_curva.png");

  // ================================================================
  //  5) TABELA
  // ================================================================
  printf("  Sinal  : foton direto isolado, PromptPhoton:all\n");
  printf("  Ruido  : pi0 + outro hadron + fragmentacao, HardQCD:all\n");
  printf("  As duas amostras somam as MESMAS 4 faixas de pThat, cada uma\n");
  printf("  pesada pela sua propria secao de choque, e passam pelo MESMO\n");
  printf("  corte de massa invariante.\n");
  printf("------------------------------------------\n");
  printf("  %10s %11s %11s %7s %19s %19s\n", "pT [GeV/c]",
         "S [mb]", "B [mb]", "N_B", "pureza so' isol.", "pureza + massa");
  printf("  (N_B = numero EFETIVO de eventos MC de fundo no bin, (c/e)^2)\n");

  for (Int_t b = 1; b <= S->GetNbinsX(); b++)
  {
    Double_t s  = S->GetBinContent(b);
    Double_t bg = Bp->GetBinContent(b) + Bo->GetBinContent(b)
                + Bf->GetBinContent(b);
    if (s + bg <= 0) continue;

    // Numero efetivo de eventos num histograma com pesos:
    //    n_ef = (soma w)^2 / (soma w^2) = (conteudo / erro)^2
    // E' quantos eventos nao pesados dariam a mesma incerteza relativa.
    Double_t eb = TMath::Sqrt(Bp->GetBinError(b)*Bp->GetBinError(b)
                            + Bo->GetBinError(b)*Bo->GetBinError(b)
                            + Bf->GetBinError(b)*Bf->GetBinError(b));
    Double_t nB = (eb > 0) ? (bg / eb) * (bg / eb) : 0.;

    char faixa[40];
    sprintf(faixa, "%.0f - %.0f",
            S->GetBinLowEdge(b), S->GetBinLowEdge(b) + S->GetBinWidth(b));

    printf("  %10s %11.4e %11.4e %7.0f", faixa, s, bg, nB);

    if (bg <= 0) printf(" %19s", "sem fundo MC");
    else         printf(" %12.1f%% +- %.1f%%",
                        100. * hist_pureza->GetBinContent(b),
                        100. * hist_pureza->GetBinError(b));

    if (hist_pureza_pos->GetBinContent(b) <= 0) printf(" %19s", "-");
    else printf(" %12.1f%% +- %.1f%%",
                100. * hist_pureza_pos->GetBinContent(b),
                100. * hist_pureza_pos->GetBinError(b));
    printf("\n");
  }

  // ---- numeros globais, integrando tudo de 5 a 40 GeV/c
  Int_t nb = hist_sinal->GetNbinsX();

  Double_t S_tot = hist_sinal->Integral(0, nb + 1);
  Double_t B_tot = hist_fundo->Integral(0, nb + 1);
  Double_t F_tot = hist_ruido[2]->Integral(0, nb + 1);

  Double_t S_pos = S_tot;                  // por suposicao, o sinal nao muda
  Double_t B_pos = hist_fundo_pos->Integral(0, nb + 1);

  printf("------------------------------------------\n");
  printf("  INTEGRADO de 5 a 40 GeV/c\n");
  printf("  %-22s %13s %13s\n", "", "so' isolamento", "+ massa");
  printf("  %-22s %13.4e %13.4e\n", "sinal (direto) [mb]", S_tot, S_pos);
  printf("  %-22s\n", "  (o sinal e' o mesmo: a marcacao so' atua no ruido)");
  printf("  %-22s %13.4e %13.4e\n", "ruido total [mb]",    B_tot, B_pos);
  printf("  %-22s %13.4e %13.4e\n", "  pi0",
         hist_ruido[0]->Integral(0, nb + 1),
         hist_ruido_pos[0]->Integral(0, nb + 1));
  printf("  %-22s %13.4e %13.4e\n", "  outro hadron",
         hist_ruido[1]->Integral(0, nb + 1),
         hist_ruido_pos[1]->Integral(0, nb + 1));
  printf("  %-22s %13.4e %13.4e\n", "  fragmentacao",
         F_tot, hist_ruido_pos[2]->Integral(0, nb + 1));

  printf("------------------------------------------\n");
  if (S_tot + B_tot > 0)
  {
    printf("  PUREZA so' com isolamento   : %6.2f%%\n",
           100. * S_tot / (S_tot + B_tot));
    printf("  PUREZA com massa invariante : %6.2f%%\n",
           (S_pos + B_pos) > 0 ? 100. * S_pos / (S_pos + B_pos) : 0.);
    printf("  ganho                       : %6.2f pontos\n",
           ((S_pos + B_pos) > 0 ? 100. * S_pos / (S_pos + B_pos) : 0.)
           - 100. * S_tot / (S_tot + B_tot));
    printf("\n");
    printf("  O que o corte de massa fez:\n");
    printf("    ruido tirado  : %6.2f%%\n",
           B_tot ? 100. * (B_tot - B_pos) / B_tot : 0.);
    printf("    sinal perdido : 0%% POR SUPOSICAO, nao por medida.\n");
    printf("      A marcacao nao e' aplicada ao sinal. Ver o cabecalho\n");
    printf("      do pureza.cc, e a coluna de fragmentacao marcada na\n");
    printf("      soma_pthat.cc, que e' a conferencia dessa suposicao.\n");
    printf("\n");
    printf("  PUREZA contando fragmentacao como sinal (so' isolamento):\n");
    printf("    %6.2f%%\n", 100. * (S_tot + F_tot) / (S_tot + B_tot));
  }

  // ---- quanta estatistica de Monte Carlo sustenta tudo isso
  {
    Double_t nB_cru = 0;
    for (Int_t o = 0; o < NORIGENS; o++) nB_cru += hist_ruido[o]->GetEntries();

    printf("------------------------------------------\n");
    printf("  ESTATISTICA DE MONTE CARLO\n");
    printf("    eventos no sinal : %.0f\n", hist_sinal->GetEntries());
    printf("    eventos no fundo : %.0f\n", nB_cru);

    if (nB_cru < 500)
    {
      printf("\n    AVISO: o fundo tem POUCOS eventos de Monte Carlo.\n");
      printf("    O valor integrado e' dominado pela regiao de pT baixo,\n");
      printf("    onde ha' estatistica, e la' ele e' confiavel. Mas a\n");
      printf("    CURVA bin a bin em pT alto nao e': um unico evento de\n");
      printf("    Monte Carlo carrega um peso grande, e os bins sem nenhum\n");
      printf("    evento de fundo ficam vazios no grafico em vez de\n");
      printf("    aparecer como pureza de 100%%, que seria falso.\n");
    }
  }
  printf("==========================================\n");

  // ================================================================
  //  6) SALVA
  // ================================================================
  TFile *arquivo = new TFile("pureza.root", "RECREATE");

  hist_sinal        ->SetDirectory(arquivo);
  hist_fundo        ->SetDirectory(arquivo);
  hist_fundo_pos    ->SetDirectory(arquivo);
  hist_pureza       ->SetDirectory(arquivo);
  hist_pureza_pos   ->SetDirectory(arquivo);
  hist_pureza_prompt->SetDirectory(arquivo);

  for (Int_t o = 0; o < NORIGENS; o++)
  {
    hist_ruido    [o]->SetDirectory(arquivo);
    hist_ruido_pos[o]->SetDirectory(arquivo);
  }

  // Sem Close(): fechar o TFile destroi os histogramas que agora pertencem
  // a ele, e os canvases ficariam em branco. O Write() grava.
  arquivo->Write();
  printf("  Arquivo salvo: pureza.root\n");
  printf("  Figuras: pureza_espectros.png  pureza_curva.png\n");
}
