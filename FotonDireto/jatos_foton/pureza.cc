// ============================================================================
//  pureza.cc  --  quanto da amostra isolada e' sinal de verdade?
// ============================================================================
//
// Responde a pergunta central do trabalho:
//
//     de tudo que passa no corte de isolamento, que fracao e' fóton DIRETO
//     (o que eu quero medir) e que fracao e' ruido (o que me engana)?
//
//     pureza(pT) = S(pT) / [ S(pT) + B(pT) ]
//
// com S e B em SECAO DE CHOQUE, nao em contagem. Tem que ser em secao de
// choque: sinal e ruido vem de amostras diferentes, geradas com numeros de
// eventos diferentes, e so' a secao de choque as poe na mesma escala fisica.
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
//   root -l 'pureza.cc(7)'     // bins de 7 GeV/c
//
// Precisa de 8 arquivos, as 4 faixas de cada amostra:
//   jatos_foton_pthat04_07.root ...   (ruido,  ./producao.sh 100000 teto)
//   foton_direto_pthat04_07.root ...  (sinal,  ver producao_direto.sh)

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
TH1F *hist_sinal   = 0;            // foton direto isolado
TH1F *hist_fundo   = 0;            // soma das 3 origens de ruido
TH1F *hist_ruido[NORIGENS];        // pi0, outro hadron, fragmentacao
TH1F *hist_pureza  = 0;            // S / (S + B), convencao "direto"
TH1F *hist_pureza_prompt = 0;      // (S + frag) / (S + B), convencao "prompt"

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

// rebin : quantos bins de 1 GeV/c agrupar na curva de pureza. O espectro
//         e' desenhado sempre na granularidade original; so' a pureza e'
//         reagrupada, porque ela e' uma razao e fica ilegivel com pouca
//         estatistica por bin.
void pureza(Int_t rebin = 5)
{
  const char *curto[NORIGENS] = {"pi0", "outro", "frag"};

  printf("\n================= PUREZA =================\n");

  // ================================================================
  //  1) SINAL: foton direto isolado, as 4 faixas somadas
  // ================================================================
  hist_sinal = soma_faixas("foton_direto", "hist_pt_direto_iso", "hist_sinal");

  if (!hist_sinal)
  {
    printf("\n  Falta a producao do SINAL nas 4 faixas de pThat.\n");
    printf("  O sinal em pThat > 5 que ja' existe NAO serve aqui: o ruido\n");
    printf("  e' pThat > 4, e comparar faixas diferentes estraga a pureza\n");
    printf("  justamente em pT baixo. Rode:\n\n");
    printf("    ./producao_direto.sh 100000\n\n");
    return;
  }

  // ================================================================
  //  2) RUIDO: as 3 origens, cada uma somada nas 4 faixas
  // ================================================================
  char apelido[100];
  for (Int_t o = 0; o < NORIGENS; o++)
  {
    char nome_hist[100];
    sprintf(nome_hist, "hist_pt_iso_%s", curto[o]);
    sprintf(apelido,   "hist_ruido_%s",  curto[o]);

    hist_ruido[o] = soma_faixas("jatos_foton", nome_hist, apelido);
    if (!hist_ruido[o])
    {
      printf("\n  Falta a producao do RUIDO. Rode:\n");
      printf("    ./producao.sh 100000 teto\n");
      return;
    }
  }

  hist_fundo = (TH1F*) hist_ruido[0]->Clone("hist_fundo");
  hist_fundo->SetDirectory(0);
  hist_fundo->Add(hist_ruido[1]);
  hist_fundo->Add(hist_ruido[2]);

  // ================================================================
  //  3) PUREZA, bin a bin
  //
  //  P = S / (S + B), com S e B vindos de amostras INDEPENDENTES.
  //
  //  Por isso NAO da' para usar TH1::Divide com a opcao "B" (binomial):
  //  binomial vale quando o numerador e' um subconjunto do denominador,
  //  medido na mesma amostra. Aqui nao e' o caso. A propagacao correta e'
  //
  //     dP/dS =  B / (S+B)^2        dP/dB = -S / (S+B)^2
  //     sigma_P^2 = [ B^2 sigma_S^2 + S^2 sigma_B^2 ] / (S+B)^4
  //
  //  Os dois histogramas tem Sumw2 ligado, entao sigma_S e sigma_B ja'
  //  vem corretos depois do Scale e do Add.
  // ================================================================
  // Copias reagrupadas, usadas so' na curva de pureza. O espectro do
  // painel 1 continua na granularidade original, de 1 GeV/c.
  TH1F *S  = (TH1F*) hist_sinal   ->Clone("pz_S");   S ->SetDirectory(0);
  TH1F *B  = (TH1F*) hist_fundo   ->Clone("pz_B");   B ->SetDirectory(0);
  TH1F *Bp = (TH1F*) hist_ruido[0]->Clone("pz_Bp");  Bp->SetDirectory(0);
  TH1F *Bo = (TH1F*) hist_ruido[1]->Clone("pz_Bo");  Bo->SetDirectory(0);
  TH1F *Bf = (TH1F*) hist_ruido[2]->Clone("pz_Bf");  Bf->SetDirectory(0);

  if (rebin > 1)
  {
    S->Rebin(rebin);  B->Rebin(rebin);
    Bp->Rebin(rebin); Bo->Rebin(rebin); Bf->Rebin(rebin);
  }

  hist_pureza = (TH1F*) S->Clone("hist_pureza");
  hist_pureza->SetDirectory(0);
  hist_pureza->Reset();

  hist_pureza_prompt = (TH1F*) S->Clone("hist_pureza_prompt");
  hist_pureza_prompt->SetDirectory(0);
  hist_pureza_prompt->Reset();

  for (Int_t b = 1; b <= S->GetNbinsX(); b++)
  {
    // As quatro componentes e os seus erros. Vem de amostras
    // INDEPENDENTES, entao os erros se somam em quadratura sem correlacao.
    Double_t s  = S ->GetBinContent(b), es  = S ->GetBinError(b);
    Double_t bp = Bp->GetBinContent(b), ebp = Bp->GetBinError(b);
    Double_t bo = Bo->GetBinContent(b), ebo = Bo->GetBinError(b);
    Double_t bf = Bf->GetBinContent(b), ebf = Bf->GetBinError(b);

    Double_t tot = s + bp + bo + bf;
    if (tot <= 0) continue;         // bin vazio nas duas amostras

    // SE O FUNDO NAO TEM NENHUM EVENTO DE MONTE CARLO NESTE BIN, a pureza
    // nao e' 100%: e' DESCONHECIDA. Zero eventos gerados nao quer dizer
    // zero secao de choque, quer dizer que a amostra e' pequena demais
    // para dizer alguma coisa aqui. Deixar o ponto em 1,00 com erro zero
    // seria mentira; o bin fica vazio no grafico.
    if (bp + bo + bf <= 0) continue;

    // Funcao auxiliar escrita na mao, porque e' usada duas vezes com
    // reparticoes diferentes entre sinal e fundo:
    //
    //    P = X / (X + Y)
    //    sigma_P = sqrt( Y^2 sigma_X^2 + X^2 sigma_Y^2 ) / (X + Y)^2

    // --- convencao "direto": so' o foton direto e' sinal
    {
      Double_t X = s;
      Double_t Y = bp + bo + bf;
      Double_t eX2 = es*es;
      Double_t eY2 = ebp*ebp + ebo*ebo + ebf*ebf;

      hist_pureza->SetBinContent(b, X / tot);
      hist_pureza->SetBinError  (b,
          TMath::Sqrt(Y*Y*eX2 + X*X*eY2) / (tot*tot));
    }

    // --- convencao "prompt": direto + fragmentacao sao sinal.
    // O denominador e' o mesmo; muda o que conta como sinal. Como o fundo
    // aqui e' so' pi0 + outro hadron, ele pode estar vazio num bin em que
    // o fundo da outra convencao nao estava. Mesma regra: sem evento de
    // Monte Carlo no fundo, o ponto nao e' 100%, e' desconhecido.
    if (bp + bo > 0)
    {
      Double_t X = s + bf;
      Double_t Y = bp + bo;
      Double_t eX2 = es*es + ebf*ebf;
      Double_t eY2 = ebp*ebp + ebo*ebo;

      hist_pureza_prompt->SetBinContent(b, X / tot);
      hist_pureza_prompt->SetBinError  (b,
          TMath::Sqrt(Y*Y*eX2 + X*X*eY2) / (tot*tot));
    }
  }

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
  l1->AddEntry(hist_sinal,   "SINAL: #gamma direto",  "l");
  l1->AddEntry(hist_fundo,   "RUIDO: soma das 3",     "l");
  l1->AddEntry(hist_ruido[0], "ruido: #pi^{0}",       "l");
  l1->AddEntry(hist_ruido[1], "ruido: outro hadron",  "l");
  l1->AddEntry(hist_ruido[2], "ruido: fragmentacao",  "l");
  l1->Draw();

  p1->Update();
  p1->SaveAs("pureza_espectros.png");

  // --- p2: a curva de pureza
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

  hist_pureza_prompt->SetLineColor(kOrange + 7);
  hist_pureza_prompt->SetLineWidth(2);
  hist_pureza_prompt->SetLineStyle(2);
  hist_pureza_prompt->SetMarkerStyle(24);
  hist_pureza_prompt->SetMarkerColor(kOrange + 7);
  hist_pureza_prompt->Draw("E1 SAME");

  TLegend *l2 = new TLegend(0.45, 0.75, 0.89, 0.89);
  l2->SetTextSize(0.030);
  l2->AddEntry(hist_pureza,        "sinal = direto",                "lp");
  l2->AddEntry(hist_pureza_prompt, "sinal = direto + fragmentacao", "lp");
  l2->Draw();

  p2->Update();
  p2->SaveAs("pureza_curva.png");

  // ================================================================
  //  5) TABELA
  // ================================================================
  printf("  Sinal  : foton direto isolado, PromptPhoton:all\n");
  printf("  Ruido  : pi0 + outro hadron + fragmentacao, HardQCD:all\n");
  printf("  As duas amostras somam as MESMAS 4 faixas de pThat,\n");
  printf("  cada uma pesada pela sua propria secao de choque.\n");
  printf("------------------------------------------\n");
  printf("  %10s %12s %12s %7s %7s %18s\n",
         "pT [GeV/c]", "S [mb]", "B [mb]", "N_S", "N_B", "pureza");
  printf("  (N_S e N_B = numero EFETIVO de eventos MC no bin, (c/e)^2)\n");

  for (Int_t b = 1; b <= S->GetNbinsX(); b++)
  {
    Double_t s  = S->GetBinContent(b), es = S->GetBinError(b);
    Double_t bg = B->GetBinContent(b), eb = B->GetBinError(b);
    if (s + bg <= 0) continue;

    // Numero efetivo de eventos num histograma com pesos:
    //    n_ef = (soma w)^2 / (soma w^2) = (conteudo / erro)^2
    // E' quantos eventos nao pesados dariam a mesma incerteza relativa.
    Double_t nS = (es > 0) ? (s / es) * (s / es) : 0.;
    Double_t nB = (eb > 0) ? (bg / eb) * (bg / eb) : 0.;

    char faixa[40];
    sprintf(faixa, "%.0f - %.0f",
            S->GetBinLowEdge(b), S->GetBinLowEdge(b) + S->GetBinWidth(b));

    if (bg <= 0)
      printf("  %10s %12.4e %12.4e %7.0f %7.0f %18s\n", faixa, s, bg,
             nS, nB, "sem fundo MC");
    else
      printf("  %10s %12.4e %12.4e %7.0f %7.0f %10.1f%% +- %.1f%%\n",
             faixa, s, bg, nS, nB,
             100. * hist_pureza->GetBinContent(b),
             100. * hist_pureza->GetBinError(b));
  }

  // ---- numero global, integrando tudo de 5 a 40 GeV/c
  Double_t S_tot = hist_sinal->Integral();
  Double_t B_tot = hist_fundo->Integral();
  Double_t F_tot = hist_ruido[2]->Integral();

  printf("------------------------------------------\n");
  printf("  INTEGRADO de 5 a 40 GeV/c\n");
  printf("    sinal (direto)      : %.6e mb\n", S_tot);
  printf("    ruido total         : %.6e mb\n", B_tot);
  printf("      do qual pi0       : %.6e mb\n", hist_ruido[0]->Integral());
  printf("      outro hadron      : %.6e mb\n", hist_ruido[1]->Integral());
  printf("      fragmentacao      : %.6e mb\n", F_tot);
  printf("------------------------------------------\n");

  if (S_tot + B_tot > 0)
  {
    printf("    PUREZA (direto)         : %.2f%%\n",
           100. * S_tot / (S_tot + B_tot));
    printf("    contaminacao            : %.2f%%\n",
           100. * B_tot / (S_tot + B_tot));
    printf("    PUREZA (direto + frag)  : %.2f%%\n",
           100. * (S_tot + F_tot) / (S_tot + B_tot));
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
      printf("    Para a curva valer em pT alto, e' preciso mais estatistica\n");
      printf("    de HardQCD -- foton isolado e energetico e' raro.\n");
    }
  }
  printf("==========================================\n");

  // ================================================================
  //  6) SALVA
  // ================================================================
  TFile *arquivo = new TFile("pureza.root", "RECREATE");

  hist_sinal        ->SetDirectory(arquivo);
  hist_fundo        ->SetDirectory(arquivo);
  hist_pureza       ->SetDirectory(arquivo);
  hist_pureza_prompt->SetDirectory(arquivo);
  for (Int_t o = 0; o < NORIGENS; o++) hist_ruido[o]->SetDirectory(arquivo);

  // Sem Close(): fechar o TFile destroi os histogramas que agora pertencem
  // a ele, e os canvases ficariam em branco. O Write() grava.
  arquivo->Write();
  printf("  Arquivo salvo: pureza.root\n");
  printf("  Figuras: pureza_espectros.png  pureza_curva.png\n");
}
